# =============================================================================
# utils_models.R —— k-TSP (multiclassPairs) wrappers, anchor selection, baseline
# methods, open-set scoring and admixture mixtures.
# Note: argument names and return structures of multiclassPairs may change between
# versions; run assert_multiclasspairs_api() before use and adapt if needed.
# =============================================================================

assert_multiclasspairs_api <- function() {
  need <- list(
    ReadData                = c("Data", "Labels"),
    filter_genes_TSP        = c("data_object", "filter", "featureNo"),
    train_one_vs_rest_TSP   = c("data_object", "filtered_genes", "k_range", "include_pivot"),
    predict_one_vs_rest_TSP = c("classifier", "Data", "tolerate_missed_genes", "weighted_votes")
  )
  for (fn in names(need)) {
    f <- getExportedValue("multiclassPairs", fn)
    miss <- setdiff(need[[fn]], names(formals(f)))
    check(length(miss) == 0, "multiclassPairs::", fn, " has no argument ", paste(miss, collapse = ","),
          ": API differs between versions; read the docs and source before editing the wrapper")
  }
  log_msg("multiclassPairs version ", as.character(utils::packageVersion("multiclassPairs")), ", API check passed")
  invisible(TRUE)
}

# ---- Training ------------------------------------------------------------------------
# anchors = NULL: main model, using the package's own within-fold filtering (featureNo genes per class)
# anchors = gene vector: anchors × whole-universe pairing (include_pivot = TRUE)
ktsp_train <- function(tv, anchors = NULL, include_pivot = TRUE, seed,
                       featureNo = CFG$genome_featureNo, k_range = CFG$k_range) {
  obj <- multiclassPairs::ReadData(Data = tv$expr, Labels = tv$meta$source,
                                   Platform = tv$meta$study, verbose = FALSE)
  fg <- multiclassPairs::filter_genes_TSP(data_object = obj, filter = "one_vs_rest",
                                          platform_wise = TRUE, featureNo = featureNo,
                                          UpDown = TRUE, verbose = FALSE)
  n_anchor <- NA_integer_
  if (!is.null(anchors)) {
    a <- intersect(anchors, rownames(tv$expr))
    check(length(a) >= 4, "only ", length(a), " anchors left in the training universe")
    # VERIFY corrected 2026-10-10: nested structure fg$OnevsrestScheme$filtered_genes
  # (original assumption $filtered_genes was wrong; E2-type lesson): structure of filtered_genes; the assumption below is fg$filtered_genes[[class name]] = gene vector
    check(is.list(fg$OnevsrestScheme$filtered_genes) && all(CFG$classes %in% names(fg$OnevsrestScheme$filtered_genes)),
          "filtered_genes structure does not match the assumption: run str(fg) to verify before changing this function")
    for (cl in CFG$classes) fg$OnevsrestScheme$filtered_genes[[cl]] <- a
    n_anchor <- length(a)
  } else include_pivot <- FALSE
  clf <- multiclassPairs::train_one_vs_rest_TSP(
    data_object = obj, filtered_genes = fg, k_range = k_range,
    include_pivot = include_pivot, one_vs_one_scores = TRUE,
    platform_wise_scores = FALSE, seed = seed, verbose = FALSE)
  attr(clf, "n_anchor") <- n_anchor
  clf
}

# ---- Prediction (missing genes skipped at rule level, no filling allowed)----------------------------------
# Fix 2026-10-11 (multiclassPairs 0.4.3 boundary bug found during LOSO runs):
#   at prediction time, if a class has exactly 1 surviving rule, the package's c$TSPs[sum == 2, ] drops to a vector without drop=FALSE,
#   the matrix(nrow = ncol(inputMat), ncol = nrow(TSPs)) call in SWAP.KTSP.Statistics, with
#   nrow(vector) = NULL, reports "non-numeric matrix extent"; with 0 surviving rules the weighted denominator is 0, giving a NaN column.
#   Handling (rule-level skip + counting, E3): classes with n≥2 go to the package as-is (bit-identical to previous results on normal folds; script 05 assertion (6));
#   classes with n=1 call switchBox::SWAP.KTSP.Statistics directly for single-rule scoring (equivalent to the package's weighted formula:
#   for a single rule, score/sum(score) = 1); classes with n=0 are all NA (not evaluable, not imputed; per utils_metrics
#   'NA counts as an error', and counted in the log here).
ktsp_predict <- function(clf, test_expr) {
  rg <- rownames(test_expr)
  n_surv <- vapply(CFG$classes, function(cl) {
    tsp <- as.matrix(clf$classifiers[[cl]]$TSPs)
    if (ncol(tsp) < 2) return(0L)
    sum(tsp[, 1] %in% rg & tsp[, 2] %in% rg)
  }, integer(1))
  eval2  <- CFG$classes[n_surv >= 2]
  single <- CFG$classes[n_surv == 1]
  none   <- CFG$classes[n_surv == 0]
  if (length(none) || length(single)) {
    log_msg("ktsp_predict degenerate classes (E3 count): ", paste0(CFG$classes, "=", n_surv, collapse = " "),
            "; single rule→switchBox: [", paste(single, collapse = ","), "]",
            "; all missing→NA: [", paste(none, collapse = ","), "]")
  }
  score_cols <- list()
  res <- NULL
  if (length(eval2)) {
    res <- multiclassPairs::predict_one_vs_rest_TSP(
      classifier = clf, Data = test_expr, tolerate_missed_genes = TRUE,
      weighted_votes = TRUE, classes = eval2, verbose = FALSE)
    check(all(c("max_score", eval2) %in% colnames(res)),
          "prediction output columns do not match the assumption: run str(res) to verify before changing this function")
    for (cl in eval2) score_cols[[cl]] <- as.numeric(res[[cl]])
  }
  for (cl in single) {
    tsp <- as.matrix(clf$classifiers[[cl]]$TSPs)
    keep <- tsp[, 1] %in% rg & tsp[, 2] %in% rg
    s1 <- clf$classifiers[[cl]]$score[keep][1]
    check(length(s1) == 1 && is.finite(s1), "single-rule class ", cl, " has malformed score field")
    kappa <- switchBox::SWAP.KTSP.Statistics(
      inputMat = as.matrix(test_expr),
      classifier = list(TSPs = tsp[keep, , drop = FALSE], score = s1, name = "TSPs1"))
    score_cols[[cl]] <- as.numeric(rowSums(t(t(kappa$comparisons) * s1)) / s1)
  }
  for (cl in none) score_cols[[cl]] <- rep(NA_real_, ncol(test_expr))
  out <- data.frame(gsm = colnames(test_expr), stringsAsFactors = FALSE)
  if (!length(single) && !length(none)) {
    out$max_score <- as.character(res$max_score)   # standard path: bit-identical to package output
  } else {
    sc <- do.call(cbind, score_cols[CFG$classes])
    out$max_score <- vapply(seq_len(nrow(sc)), function(i) {
      v <- sc[i, ]
      if (all(is.na(v))) NA_character_ else CFG$classes[which.max(replace(v, is.na(v), -Inf))]
    }, character(1))
  }
  out$pred <- out$max_score
  for (cl in CFG$classes) out[[paste0("score_", cl)]] <- score_cols[[cl]]
  # evaluable classes must not contain NaN; degenerate counts are carried with the return value
  sc_mat <- as.matrix(out[, paste0("score_", CFG$classes)])
  nan_ok <- !any(is.nan(sc_mat[, !CFG$classes %in% none, drop = FALSE]))
  check(nan_ok, "ktsp_predict: NaN in scores of evaluable classes")
  attr(out, "n_degenerate_classes") <- length(single) + length(none)
  attr(out, "n_surviving_rules") <- n_surv
  out
}

# ---- Rule table: extract gene pairs from the classifier and set direction empirically from training data ----------------------
# Does not rely on the package's 'which direction votes' semantics (unclear in v1); computed directly from data
ktsp_rules <- function(clf, train_expr, train_labels) {
  # rule-structure assertion (VERIFY marker to be backfilled after script 05 verification, M5 remediation)
  check(!is.null(clf$classifiers) && all(CFG$classes %in% names(clf$classifiers)),
        "classifier structure does not match the assumption: run str(clf, max.level = 3) to verify before changing this function")
  do.call(rbind, lapply(CFG$classes, function(cl) {
    tsp <- clf$classifiers[[cl]]$TSPs
    check(!is.null(tsp) && ncol(as.matrix(tsp)) >= 2, "rule matrix of class ", cl, " is empty or malformed")
    tsp <- as.matrix(tsp)
    # M5 observation: the package keeps NA placeholder rows for ignored rules ('Ignore N rules' under tolerate_missed_genes),
    # drop them before reading rules, or row-name indexing goes out of bounds
    tsp <- tsp[!is.na(tsp[, 1]) & !is.na(tsp[, 2]), , drop = FALSE]
    g1 <- tsp[, 1]; g2 <- tsp[, 2]
    in_cl <- train_labels == cl
    f_cl <- vapply(seq_along(g1), function(i) mean(train_expr[g1[i], in_cl] > train_expr[g2[i], in_cl]), numeric(1))
    f_rest <- vapply(seq_along(g1), function(i) mean(train_expr[g1[i], !in_cl] > train_expr[g2[i], !in_cl]), numeric(1))
    data.frame(class = cl, rule_rank = seq_along(g1), gene1 = g1, gene2 = g2,
               frac_g1_gt_g2_in_class = f_cl, frac_g1_gt_g2_in_rest = f_rest,
               class_direction = ifelse(f_cl >= f_rest, "gene1>gene2", "gene1<gene2"),
               separation = abs(f_cl - f_rest), stringsAsFactors = FALSE)
  }))
}

# ---- Anchors --------------------------------------------------------------------------
row_anova_F <- function(x, g) {
  g <- factor(g); k <- nlevels(g); n <- ncol(x)
  if (k < 2 || n - k < 1) return(rep(NA_real_, nrow(x)))
  mu <- rowMeans(x); ssb <- 0; ssw <- 0
  for (lv in levels(g)) {
    xi <- x[, g == lv, drop = FALSE]; mi <- rowMeans(xi)
    ssb <- ssb + ncol(xi) * (mi - mu)^2
    ssw <- ssw + rowSums((xi - mi)^2)
  }
  (ssb / (k - 1)) / pmax(ssw / (n - k), .Machine$double.eps)
}

# genome anchors: F percentile computed within each training study, then averaged across studies (no cross-study normalisation)
genome_anchor_rank <- function(tv, n = CFG$anchor_size) {
  pr <- lapply(split(seq_len(nrow(tv$meta)), tv$meta$study), function(idx) {
    f <- row_anova_F(tv$expr[, idx, drop = FALSE], tv$meta$source[idx])
    if (all(is.na(f))) return(NULL)
    rank(f, na.last = "keep") / sum(!is.na(f))
  })
  pr <- Filter(Negate(is.null), pr)
  check(length(pr) >= 2, "fewer than 2 training studies allow F computation")
  score <- rowMeans(do.call(cbind, pr), na.rm = TRUE)
  names(score) <- rownames(tv$expr)
  names(sort(score, decreasing = TRUE))[seq_len(n)]
}

# random anchors: match the target set by within-fold expression and detection rate (HOX by default); generate to disk before training
make_random_anchor_sets <- function(feat, target_genes, n_sets, seed_base) {
  tg <- intersect(target_genes, feat$gene)
  check(length(tg) >= 4, "fewer than 4 target anchors in the training universe")
  targets <- data.frame(gene = tg, contrast = "all", stringsAsFactors = FALSE)
  nl <- sample_matched_null(targets, feat, n_sets = n_sets, seed_base = seed_base)
  split(nl$gene, nl$set)
}

# ---- Baseline methods (same fold, same test samples)------------------------------------------
rank_transform <- function(x) {
  apply(x, 2, function(v) (rank(v, ties.method = "average") - 0.5) / length(v))
}

# Uses only genes shared by training ∩ test (only test gene names used, no labels)
baseline_glmnet <- function(tv, test_expr, seed, alpha = 0.5, return_fit = FALSE) {
  genes <- intersect(rownames(tv$expr), rownames(test_expr))
  xtr <- t(rank_transform(tv$expr[genes, , drop = FALSE]))
  xte <- t(rank_transform(test_expr[genes, , drop = FALSE]))
  foldid <- as.integer(factor(tv$meta$cluster_id))
  check(length(unique(foldid)) >= 3, "glmnet grouped cross-validation requires ≥3 training clusters")
  set.seed(seed)
  cv <- glmnet::cv.glmnet(xtr, factor(tv$meta$source, levels = CFG$classes), family = "multinomial",
                          alpha = alpha, foldid = foldid, type.measure = "class")
  pred <- as.character(stats::predict(cv, xte, s = "lambda.1se", type = "class"))
  nz <- length(unique(unlist(lapply(stats::coef(cv, s = "lambda.1se"), function(b) rownames(b)[as.numeric(b) != 0]))))
  out <- list(pred = data.frame(gsm = colnames(test_expr), pred = pred, stringsAsFactors = FALSE),
              n_genes = nz - 1L)   # exclude intercept
  if (return_fit) out$fit <- cv     # used by the final model script (44); default behaviour unchanged
  out
}

baseline_centroid <- function(tv, test_expr) {
  genes <- intersect(rownames(tv$expr), rownames(test_expr))
  rtr <- rank_transform(tv$expr[genes, , drop = FALSE])
  rte <- rank_transform(test_expr[genes, , drop = FALSE])
  cen <- sapply(CFG$classes, function(cl) rowMeans(rtr[, tv$meta$source == cl, drop = FALSE]))
  sim <- stats::cor(rte, cen)
  data.frame(gsm = colnames(test_expr), pred = CFG$classes[max.col(sim, ties.method = "first")],
             stringsAsFactors = FALSE)
}

# ---- Open-set rejection ------------------------------------------------------------------
# after score normalisation, take the 'top1 - top2' difference; the threshold comes from the distribution of correctly predicted LOSO samples
# 2026-10-11: NA tolerance -- score columns of all-missing classes are NA (E3 fix), normalisation uses available classes;
# when fewer than 2 classes are available the sample margin is NA (downstream reject treats it as rejected; NA counts as an error)
score_margin <- function(pred_df) {
  s <- as.matrix(pred_df[, paste0("score_", CFG$classes)])
  s <- s / pmax(rowSums(s, na.rm = TRUE), .Machine$double.eps)
  apply(s, 1, function(v) {
    v <- sort(v[!is.na(v)], decreasing = TRUE)
    if (length(v) < 2) NA_real_ else v[1] - v[2]
  })
}

choose_tau <- function(loso_pred, coverage = 0.95) {
  ok <- !is.na(loso_pred$pred) & loso_pred$pred == loso_pred$true
  unname(stats::quantile(score_margin(loso_pred[ok, ]), 1 - coverage, na.rm = TRUE))
}

apply_reject <- function(pred_df, tau) {
  pred_df$margin <- score_margin(pred_df)
  pred_df$rejected <- pred_df$margin < tau
  pred_df$pred_openset <- ifelse(pred_df$rejected, NA_character_, pred_df$pred)
  pred_df
}

# ---- Mixtures (mix on linear scale, then log; simulation assumption to be stated in the paper)-----------------
mix_log2 <- function(a, b, f, offset = 0) {
  log2(f * (2^a - offset) + (1 - f) * (2^b - offset) + offset)
}
