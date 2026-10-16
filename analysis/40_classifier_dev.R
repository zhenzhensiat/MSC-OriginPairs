# =============================================================================
# 40_classifier_dev.R —— classifier development (development layer only,
# leave-one-group-out cross-validation).
#
# Methods (all evaluated on the same folds and training universe):
#   genome_main   primary model: in-fold gene filtering + all-pair search
#                 (per-fold models saved for the admixture simulation)
#   HOX / published / genome39   ~39 anchors each x all-pair search (include_pivot)
#   random_001..  random anchors matched to HOX on in-fold expression/detection
#                 (null distribution)
#   genome_k3..   primary model at fixed k (accuracy vs number of rules curve)
#   glmnet / centroid   baseline methods (within-sample rank transform)
#
# Smoke test: MSC_SMOKE=1 restricts to four representative folds x 3 anchor sets.
# Resumable: each method writes its own RDS; after code changes delete the
# output directory and re-run.
# =============================================================================
source("scripts/00_setup.R")
need_pkgs(c("multiclassPairs", "glmnet"))

SMOKE <- identical(Sys.getenv("MSC_SMOKE"), "1")
OUT <- if (SMOKE) "data/loso_smoke" else "data/loso_v2r2"   # run1 archived; new directory v2r2

run_step(if (SMOKE) "40_classifier_dev_SMOKE" else "40_classifier_dev", function() {
  assert_multiclasspairs_api()
  studies <- load_studies("D")
  meta <- bind_meta(studies)
  folds <- make_folds(meta)
  folds <- folds[!names(folds) %in% CFG$loso_exclude_clusters]   # D2: pre-screened matrices excluded from fold evaluation
  if (SMOKE) {
    want <- c("G09", "G01", "G12", "G11")
    folds <- folds[names(folds) %in% want]
    check(length(folds) == length(want), "smoke-test specified folds missing: ",
          paste(setdiff(want, names(folds)), collapse = ","))
  }
  n_rand <- if (SMOKE) 3L else CFG$n_random_anchor_sets
  fold_idx <- stats::setNames(seq_along(folds), names(folds))
  dir.create(proj_path(OUT, "models"), recursive = TRUE, showWarnings = FALSE)

  hox <- read_tsv("curation/hox_list.tsv")$gene
  pub <- unique(read_tsv("curation/marker_list.tsv")$gene)

  # training-view cache: all methods share the same training studies and gene universe within a fold
  tv_cache <- new.env(parent = emptyenv())
  get_tv <- function(cl) {
    key <- as.character(cl)
    if (!exists(key, envir = tv_cache, inherits = FALSE)) {
      assign(key, training_view(studies, cl, restrict_genes = test_genes(studies, cl)), envir = tv_cache)
    }
    get(key, envir = tv_cache)
  }

  # ---- 1. in-fold anchors: generate and persist all first, then train------------------------------
  prep_file <- proj_path(OUT, "fold_prep.rds")
  if (!file.exists(prep_file)) {
    prep <- lapply(names(folds), function(f) {
      tv <- get_tv(folds[[f]])
      hu <- length(intersect(hox, tv$genes))
      log_msg("fold ", f, ": universe ", length(tv$genes), " genes (training studies ", length(tv$studies),
              "), HOX retained ", hu)
      list(genome39 = genome_anchor_rank(tv, CFG$anchor_size),
           random = make_random_anchor_sets(gene_features_train(tv), hox, n_rand,
                                            derive_seed(40, fold_idx[[f]] * 1000L)),
           n_train_studies = length(tv$studies), n_universe = length(tv$genes),
           hox_survive = hu)
    })
    names(prep) <- names(folds)
    saveRDS(prep, prep_file)
  }
  prep <- readRDS(prep_file)

  # ---- 2. method definitions ----------------------------------------------------------------
  ktsp_fp <- function(method, anchor_of, k_range = CFG$k_range, save_model = FALSE) {
    force(method); force(anchor_of); force(k_range); force(save_model)
    function(cl) {
      key <- as.character(cl)
      tv <- get_tv(cl)
      anchors <- anchor_of(key)
      clf <- ktsp_train(tv, anchors = anchors, include_pivot = !is.null(anchors),
                        seed = derive_seed(41, fold_idx[[key]]), k_range = k_range)
      if (save_model) saveRDS(clf, proj_path(OUT, "models", paste0(method, "_", key, ".rds")))
      prs <- lapply(test_view(studies, cl), function(tt) ktsp_predict(clf, tt$expr))
      n_deg <- sum(vapply(prs, function(p) attr(p, "n_degenerate_classes"), integer(1)))
      list(pred = do.call(rbind, prs), n_skipped = n_deg)
    }
  }

  methods <- list(
    genome_main = ktsp_fp("genome_main", function(key) NULL, save_model = TRUE),
    HOX         = ktsp_fp("HOX", function(key) hox),
    published   = ktsp_fp("published", function(key) pub),
    genome39    = ktsp_fp("genome39", function(key) prep[[key]]$genome39)
  )
  for (k in CFG$k_range) {
    nm <- paste0("genome_k", k)
    methods[[nm]] <- ktsp_fp(nm, function(key) NULL, k_range = k)
  }
  methods$glmnet <- function(cl) {
    key <- as.character(cl); tv <- get_tv(cl)
    pr <- do.call(rbind, lapply(test_view(studies, cl), function(tt)
      baseline_glmnet(tv, tt$expr, seed = derive_seed(42, fold_idx[[key]]))$pred))
    list(pred = pr, n_skipped = NA)
  }
  methods$centroid <- function(cl) {
    tv <- get_tv(cl)
    pr <- do.call(rbind, lapply(test_view(studies, cl), function(tt) baseline_centroid(tv, tt$expr)))
    list(pred = pr, n_skipped = NA)
  }
  for (s in seq_len(n_rand)) {
    nm <- sprintf("random_%03d", s)
    methods[[nm]] <- ktsp_fp(nm, local({ s0 <- as.character(s); function(key) prep[[key]]$random[[s0]] }))
  }

  # ---- 3. per-method LOSO (resumable)-------------------------------------------------
  res <- list()
  for (nm in names(methods)) {
    f_out <- proj_path(OUT, paste0("loso_", nm, ".rds"))
    if (!file.exists(f_out)) {
      log_msg("LOSO: ", nm)
      saveRDS(run_loso(folds, nm, methods[[nm]], meta), f_out)
    }
    res[[nm]] <- readRDS(f_out)
    mem_guard(nm)
  }
  fold_tab <- do.call(rbind, lapply(res, `[[`, "folds"))
  write_tsv(fold_tab, file.path(OUT, "loso_folds.tsv"))   # includes failed rows

  # ---- 4. compare only on folds where all named methods succeeded (R5)-------------------------------
  rnd <- grep("^random_", names(res), value = TRUE)
  named <- setdiff(names(res), rnd)
  common <- common_ok_folds(lapply(res[named], `[[`, "folds"))
  log_msg("folds with all named methods successful: ", length(common), " / ", length(folds))
  bad <- fold_tab[fold_tab$method %in% named & grepl("NaN", fold_tab$note), ]
  check(nrow(bad) == 0, "named methods raised NaN warnings: ", paste(unique(bad$method), collapse = ","))
  check(all(fold_tab$status[fold_tab$method %in% named] == "ok"), "named methods have failed folds (Q3: 100% required)")
  check(length(common) >= ceiling(0.8 * length(folds)),
        "fewer than 80% common successful folds: fix the failed folds first (see loso_folds.tsv); do not compare directly")

  pred_common <- function(nm) {
    pr <- res[[nm]]$pred
    attach_truth(pr[pr$fold %in% common, , drop = FALSE], meta)
  }
  summ <- do.call(rbind, lapply(named, function(nm) {
    df <- pred_common(nm)
    out <- flat_summary(summarize_predictions(df), nm)
    b <- cluster_bootstrap(df, function(d) balanced_accuracy(d$true, d$pred), seed = derive_seed(43, 1))
    out$ba_lo <- b$lo; out$ba_hi <- b$hi; out$ci_unreliable <- b$unreliable
    out
  }))
  summ$n_common_folds <- length(common)
  write_tsv(summ, file.path(OUT, "loso_summary.tsv"))
  # column-union padding (approved by the requester on 2026-10-12): baseline methods (glmnet/centroid) have no score_* columns; pad with NA before merging
  preds_all <- lapply(named, pred_common)
  all_cols <- Reduce(union, lapply(preds_all, names))
  preds_all <- lapply(preds_all, function(d) {
    for (cc in setdiff(all_cols, names(d))) d[[cc]] <- NA
    d[, all_cols, drop = FALSE]
  })
  write_tsv(do.call(rbind, preds_all), file.path(OUT, "loso_predictions.tsv"))

  # ---- 5. anchor source vs random-anchor null distribution (R6)-----------------------------------------
  rnd_fail <- rnd[vapply(rnd, function(nm) any(res[[nm]]$folds$status != "ok"), logical(1))]
  log_msg("random-anchor sets with failed folds: ", length(rnd_fail), " / ", length(rnd), " (failed sets leave the null distribution; Q3 threshold 10%)")
  check(length(rnd_fail) <= 0.10 * length(rnd), "random-set failure rate exceeds 10%; stop and investigate: ",
        paste(rnd_fail, collapse = ","))
  ok_fold <- function(nm) res[[nm]]$folds$fold[res[[nm]]$folds$status == "ok"]
  rnd_ok <- rnd[vapply(rnd, function(nm) all(common %in% ok_fold(nm)), logical(1))]
  rnd_ba <- vapply(rnd_ok, function(nm) { df <- pred_common(nm); balanced_accuracy(df$true, df$pred) }, numeric(1))
  write_tsv(data.frame(set = rnd_ok, ba = rnd_ba, stringsAsFactors = FALSE), file.path(OUT, "random_anchor_null.tsv"))
  cmp <- summ[summ$model %in% c("HOX", "published", "genome39"), c("model", "ba", "ba_lo", "ba_hi")]
  cmp$random_mean <- mean(rnd_ba)
  cmp$random_q95 <- unname(stats::quantile(rnd_ba, 0.95))
  cmp$percentile_in_random <- vapply(cmp$ba, function(b) mean(rnd_ba <= b), numeric(1))
  cmp$empirical_p <- vapply(cmp$ba, function(b) empirical_p(b, rnd_ba), numeric(1))
  cmp$n_random_valid <- length(rnd_ok); cmp$n_random_total <- length(rnd)
  write_tsv(cmp, file.path(OUT, "anchor_compare.tsv"))

  # ---- 6. rule-count curve (genome_k*)-------------------------------------------------
  pk <- summ[grepl("^genome_k", summ$model), c("model", "ba", "ba_lo", "ba_hi")]
  pk$k_per_class <- as.integer(sub("genome_k", "", pk$model))
  write_tsv(pk[order(pk$k_per_class), ], file.path(OUT, "panel_size_curve.tsv"))

  # ---- 7. open-set thresholds (one copy per named k-TSP method)--------------------------------
  for (nm in c("genome_main", "HOX", "published", "genome39")) {
    saveRDS(choose_tau(pred_common(nm)), proj_path(OUT, paste0("tau_", nm, ".rds")))
  }

  # ---- 8. advisory only; no automatic finalization (final choice confirmed separately)-------------------------------
  cand <- summ[summ$model %in% c("genome_main", "HOX", "published", "genome39", "glmnet"), ]
  best <- max(cand$ba)
  tied <- cand$model[cand$ba >= best - CFG$tie_tolerance]
  log_msg("best primary metric ", round(best, 3), "; within tolerance ", CFG$tie_tolerance, ", methods treated as tied: ",
          paste(tied, collapse = ", "), ". Per the tie rule (smaller panel -> higher detection rate), confirm the final choice")
  invisible(summ)
})
