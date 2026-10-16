# =============================================================================
# 92_fig_extract.R - read-only extraction of figure input tables.
# Inputs: result files plus four historic caches cross-anchored and then read.
# Outputs: fig5a_vote_matrix_v2, fig5a_vote_recompute_v2,
#   fig4b_confusion_v2, fig4d_rule_genes_v2, fig5d_platform_accuracy_v2,
#   fig2f_direction_consistency_v2, fig2c_null_rates_published_v2,
#   fig3a_null_rates_hox_v2, fig4c_random_anchor_meta_v2.
# Hard check: weighted-vote recomputation must equal confirm_predictions_v2 for
#   all 176 samples (NA==NA included); stop otherwise.
# Discipline: read-only extraction; no statistical tests recomputed.
# =============================================================================
source("scripts/00_setup.R")
need_pkgs(c("multiclassPairs", "switchBox"))

run_step("92_fig_extract", function() {

  # ---- 0. Historical caches: cross-anchor -> reuse --------------------------------------
  caches <- c("data/within_study_de_v2.rds", "data/pc_precomputed_v2.rds",
              "data/null_draws_v2.rds", "data/loso_v2r2/fold_prep.rds")
  for (f in caches) log_msg("cache ", f, ": mtime=", format(file.mtime(proj_path(f)), "%Y-%m-%d %H:%M"))
  de <- readRDS(proj_path("data", "within_study_de_v2.rds"))
  q05 <- read_tsv("data/within_study_de_v2_q05.tsv")   # derived artifact (anchor)
  deq <- de[de$q_within < 0.05, , drop = FALSE]
  check(nrow(deq) == nrow(q05), "q05 row-count reconciliation failed: ", nrow(deq), " vs ", nrow(q05))
  numc <- c("logFC", "t", "df", "p", "q_within", "expr_pct", "detect")
  okq <- identical(as.character(deq$study), as.character(q05$study)) &&
    identical(as.character(deq$contrast), as.character(q05$contrast)) &&
    identical(as.character(deq$gene), as.character(q05$gene)) &&
    all(vapply(numc, function(cc) isTRUE(all.equal(deq[[cc]], q05[[cc]], tolerance = 1e-9)),
               logical(1)))
  check(okq, "cached within_study_de does not match the q05 anchor")

  pre <- readRDS(proj_path("data", "pc_precomputed_v2.rds"))
  mt_ref <- read_tsv("data/marker_transfer_v2.tsv")
  j <- match(paste(mt_ref$gene, mt_ref$contrast, sep = "|"),
             paste(pre$gene, pre$contrast, sep = "|"))
  check(all(mt_ref$K[is.na(j)] == 0), "units missing from pre should have K=0 (never measured); actual: ",
        paste(mt_ref$gene[is.na(j)], mt_ref$K[is.na(j)], collapse = " / "))

  draws <- readRDS(proj_path("data", "null_draws_v2.rds"))
  check(nrow(draws$published) == 99000L && nrow(draws$hox) == 117000L, "null_draws row count anomaly")

  # Deep cache reconciliation (strong anchor): recompute both transfer tables with the same code as-is (de+pre joint computation, stricter than the K column alone)
  hoxg <- read_tsv("curation/hox_list.tsv")$gene
  hox_t_full <- expand.grid(gene = hoxg, contrast = CFG$contrasts, stringsAsFactors = FALSE)
  hox_t_full$direction <- 0L; hox_t_full$exclude <- ""
  mk_t_full <- expand_claims(read_tsv("curation/marker_list.tsv"))
  hk <- read_tsv("data/hox_transfer_v2.tsv")
  comp_cols <- c("K", "pc_p", "dir_used", "testable", "pc_q", "transfer", "direction")
  chk2 <- function(chk, ref, name) {
    # Order-independent: align on the quadruple key (claim_id|gene|contrast|direction) and reconcile column by column
    key <- function(d) paste(if ("claim_id" %in% names(d)) as.character(d$claim_id) else "",
                             d$gene, d$contrast, d$direction, sep = "|")
    kc <- key(chk); kf <- key(ref)
    check(!anyDuplicated(kc) && !anyDuplicated(kf), name, ": duplicate quadruple key")
    m <- match(kf, kc)
    check(!anyNA(m), name, ": reference rows missing from the recomputation (extra ",
          paste(kc[!(kc %in% kf)], collapse = ","), ")")
    for (cc in comp_cols) {
      check(isTRUE(all.equal(chk[[cc]][m], ref[[cc]], tolerance = 1e-9)),
            name, ": column ", cc, " does not match the reference table")
    }
    log_msg("cache joint-computation reconciliation passed: ", name, " (", nrow(froz), " rows, quadruple keys aligned; ",
            "row-order differences stem from safe_merge leaving the order of duplicate keys undefined; content is identical)")
  }
  chk2(family_calls(pc_rows(mk_t_full, de, pre)), mt_ref, "marker_transfer")
  chk2(family_calls(pc_rows(hox_t_full, de, pre)), hk, "hox_transfer")

  for (f in caches) check(file.exists(proj_path(f)), "missing input: ", f)

  # ---- 1. Fig 5a: vote matrix (30 rules x 176 samples) + 176/176 hard check ------------------
  clf <- readRDS(proj_path("data", "model_v2.rds"))
  rules <- read_tsv("data/model_v2_rules.tsv")
  check(nrow(rules) == 30L, "Rule count should be 30")
  studies <- load_studies("C")
  vm_rows <- list(); rc_rows <- list()
  for (s in studies) {
    m <- s$meta[s$meta$arm == "include" & s$meta$source %in% CFG$classes, , drop = FALSE]
    E <- s$expr[, m$gsm, drop = FALSE]; rg <- rownames(E)
    sc_mat <- matrix(NA_real_, nrow = ncol(E), ncol = length(CFG$classes),
                     dimnames = list(colnames(E), CFG$classes))
    for (cl in CFG$classes) {
      rr <- rules[rules$class == cl, , drop = FALSE]
      tsp_full <- as.matrix(clf$classifiers[[cl]]$TSPs)
      score_full <- clf$classifiers[[cl]]$score
      check(length(score_full) == nrow(tsp_full), "score and TSPs row counts differ (class ", cl, ")")
      okrow <- !is.na(tsp_full[, 1]) & !is.na(tsp_full[, 2])   # M5: drop NA placeholder rows
      g1m <- unname(tsp_full[okrow, 1]); g2m <- unname(tsp_full[okrow, 2]); scm <- score_full[okrow]
      check(identical(g1m, rr$gene1) && identical(g2m, rr$gene2),
            "rule table and model object mismatch (class ", cl, ")")
      present <- (g1m %in% rg) & (g2m %in% rg)
      sup <- matrix(NA, length(g1m), ncol(E), dimnames = list(NULL, colnames(E)))
      if (any(present)) {
        gt <- E[g1m[present], , drop = FALSE] > E[g2m[present], , drop = FALSE]
        flip <- rr$class_direction[present] == "gene1<gene2"
        gt[flip, ] <- !gt[flip, , drop = FALSE]
        sup[present, ] <- gt
      }
      w <- scm[present]
      if (sum(present) >= 2L) {
        sc_mat[, cl] <- colSums(sup[present, , drop = FALSE] * w) / sum(w)
      } else if (sum(present) == 1L) {
        sc_mat[, cl] <- as.numeric(sup[present, ])
      }   # 0 rules: keep NA
      vm_rows[[paste(s$gse, cl)]] <- data.frame(
        gsm = rep(colnames(E), each = length(g1m)),
        class = cl, rule_rank = seq_along(g1m), gene1 = rep(g1m, times = ncol(E)),
        gene2 = rep(g2m, times = ncol(E)),
        class_direction = rep(rr$class_direction, times = ncol(E)),
        state = ifelse(is.na(as.vector(t(sup))), "missing",
                       ifelse(as.vector(t(sup)), "holds", "not")),
        stringsAsFactors = FALSE)
    }
    pred_rec <- vapply(seq_len(ncol(E)), function(i) {
      v <- sc_mat[i, ]
      if (all(is.na(v))) NA_character_ else CFG$classes[which.max(replace(v, is.na(v), -Inf))]
    }, character(1))
    rc_rows[[s$gse]] <- data.frame(gsm = colnames(E), study = s$gse,
                                   score_AD = sc_mat[, "AD"], score_BM = sc_mat[, "BM"],
                                   score_UC = sc_mat[, "UC"], pred_recomputed = pred_rec,
                                   stringsAsFactors = FALSE)
  }
  vm <- do.call(rbind, vm_rows); rownames(vm) <- NULL
  rc <- do.call(rbind, rc_rows); rownames(rc) <- NULL
  check(nrow(vm) == 30L * 176L && nrow(rc) == 176L, "vote matrix size anomaly")

  cp <- read_tsv("data/confirm_predictions_v2.tsv")
  check(nrow(cp) == 176L, "confirm_predictions should have 176 rows")
  cpc <- cp[, c("gsm", "pred", "score_AD", "score_BM", "score_UC")]
  names(cpc) <- c("gsm", "pred_confirm", "score_AD_confirm", "score_BM_confirm", "score_UC_confirm")
  rc <- safe_merge(rc, cpc, by = "gsm", type = "1:1")
  rc$match <- (rc$pred_recomputed == rc$pred_confirm) |
    (is.na(rc$pred_recomputed) & is.na(rc$pred_confirm))
  n_ok <- sum(rc$match)
  if (n_ok != 176L) {
    bad <- rc[!rc$match, c("gsm", "pred_recomputed", "pred_confirm")]
    print(bad)
    check(FALSE, "vote recomputation disagrees with confirm: ", n_ok, "/176")
  }
  scd <- vapply(CFG$classes, function(cl) {
    v1 <- rc[[paste0("score_", cl)]]; v2 <- rc[[paste0("score_", cl, "_confirm")]]
    same_na <- is.na(v1) & is.na(v2)
    check(all(same_na | (!is.na(v1) & !is.na(v2))), cl, ": NA patterns differ")
    if (all(same_na)) return(0) else max(abs(v1[!same_na] - v2[!same_na]))
  }, numeric(1))
  check(max(scd) < 1e-6, "recomputed scores differ too much from the reference scores: ", signif(max(scd), 3))
  log_msg("Fig5a hard check passed: recomputed vote pred == confirm_predictions 176/176; ",
          "max score difference ", signif(max(scd), 3))
  write_tsv(vm, "data/fig5a_vote_matrix_v2.tsv")
  write_tsv(rc, "data/fig5a_vote_recompute_v2.tsv")

  # ---- 2. Fig 4b: LOSO confusion matrix (with BM<->AD mutual error counts)--------------------------------
  lp <- read_tsv("data/loso_v2r2/loso_predictions.tsv")
  lpm <- lp[lp$method == "genome_main", , drop = FALSE]
  check(nrow(lpm) == 174L, "loso genome_main should have 174 rows")
  check(!anyNA(lpm$pred), "LOSO predictions contain NA (abstain_rate should be 0)")
  lpm$true <- factor(lpm$true, levels = CFG$classes); lpm$pred <- factor(lpm$pred, levels = CFG$classes)
  tb <- as.data.frame(table(true = lpm$true, pred = lpm$pred), stringsAsFactors = FALSE)
  names(tb)[3] <- "n"; tb$n <- as.integer(tb$n)
  tb$true <- as.character(tb$true); tb$pred <- as.character(tb$pred)
  write_tsv(tb, "data/fig4b_confusion_v2.tsv")
  n_bm_ad <- sum(lpm$true == "BM" & lpm$pred == "AD")
  n_ad_bm <- sum(lpm$true == "AD" & lpm$pred == "BM")
  log_msg("Fig4b: BM<->AD mutual errors - BM->AD ", n_bm_ad, "; AD->BM ", n_ad_bm)

  # ---- 3. Fig 4d: annotation of the 56 rule genes -------------------------------------------------
  rg_genes <- sort(unique(c(rules$gene1, rules$gene2)))
  check(length(rg_genes) == 56L, "Rule gene count should be 56, actual ", length(rg_genes))
  mkg <- unique(read_tsv("curation/marker_list.tsv")$hgnc_symbol)
  gi <- data.frame(gene = rg_genes, is_hox = rg_genes %in% hoxg,
                   is_published_marker = rg_genes %in% mkg, stringsAsFactors = FALSE)
  gi$n_rules_in <- vapply(rg_genes, function(g) sum(rules$gene1 == g | rules$gene2 == g), integer(1))
  gi$classes <- vapply(rg_genes, function(g)
    paste(sort(unique(rules$class[rules$gene1 == g | rules$gene2 == g])), collapse = ";"), character(1))
  write_tsv(gi, "data/fig4d_rule_genes_v2.tsv")
  log_msg("Fig4d: HOX among the 56 genes: ", sum(gi$is_hox), " (",
          paste(gi$gene[gi$is_hox], collapse = ","), "); previously reported markers ",
          sum(gi$is_published_marker), " (",
          paste(gi$gene[gi$is_published_marker], collapse = ","), ")")

  # ---- 4. Fig 5d: external validation grouped by platform ------------------------------------------------
  cp$platform <- ifelse(cp$scale == "array", "microarray", "RNA-seq")
  pl <- do.call(rbind, lapply(c("microarray", "RNA-seq"), function(p) {
    d <- cp[cp$platform == p, , drop = FALSE]
    rec <- function(cl) { dd <- d[d$true == cl, , drop = FALSE]
      mean(!is.na(dd$pred) & dd$pred == cl) }
    data.frame(platform = p, n = nrow(d), n_studies = length(unique(d$study)),
               ba = balanced_accuracy(d$true, d$pred),
               recall_AD = rec("AD"), recall_BM = rec("BM"), recall_UC = rec("UC"),
               stringsAsFactors = FALSE)
  }))
  check(sum(pl$n) == 176L, "platform-grouped sample counts should sum to 176")
  write_tsv(pl, "data/fig5d_platform_accuracy_v2.tsv")
  for (i in seq_len(nrow(pl))) log_msg("Fig5d platform ", pl$platform[i], ": n=", pl$n[i],
                                       ", BA=", round(pl$ba[i], 4))

  # ---- 5. Fig 2f: proportion of studies with a consistent direction, per unit --------------------------------------------
  mt <- mt_ref[mt_ref$testable %in% c(TRUE, "TRUE"), , drop = FALSE]
  check(nrow(mt) == 99L, "Testable units should be 99, actual ", nrow(mt))
  mt$direction <- as.integer(mt$direction)
  check(all(mt$direction %in% c(-1L, 1L)), "direction value outside +/-1")
  de_key <- paste(de$gene, de$contrast, sep = "|")
  check(all(mt$dir_used == mt$direction), "dir_used and direction disagree (published units)")
  f2 <- do.call(rbind, lapply(seq_len(nrow(mt)), function(i) {
    sub <- which(de_key == paste(mt$gene[i], mt$contrast[i], sep = "|"))
    ex <- if (nzchar(mt$exclude[i])) trimws(strsplit(mt$exclude[i], "[;,]", perl = TRUE)[[1]]) else character(0)
    n_removed <- sum(de$study[sub] %in% ex)
    check(length(sub) - n_removed == mt$K[i], "K does not match the cached row count (after excluding original discovery studies): ",
          mt$gene[i], "/", mt$contrast[i], " (cache=", length(sub), " removed ", n_removed,
          " K=", mt$K[i], ")")
    sub <- sub[!(de$study[sub] %in% ex)]
    n_same <- sum(sign(de$logFC[sub]) == mt$dir_used[i])
    data.frame(gene = mt$gene[i], contrast = mt$contrast[i], claim_id = mt$claim_id[i],
               direction = as.integer(mt$direction[i]), K = mt$K[i], n_same_direction = n_same,
               frac_same_direction = n_same / mt$K[i], stringsAsFactors = FALSE)
  }))
  rownames(f2) <- NULL
  log_msg("Fig2f: ", sum(nzchar(mt$exclude)), " units counted after excluding original discovery studies")
  write_tsv(f2, "data/fig2f_direction_consistency_v2.tsv")
  log_msg("Fig2f: 99 units; median consistent fraction by contrast ",
          paste(vapply(split(f2$frac_same_direction, f2$contrast),
                       function(v) sprintf("%.3f", stats::median(v)), character(1)),
                collapse = " / "))

  # ---- 6. Fig 2c / 3a: two null-distribution rate sets (null_rates reused as-is)--------------------------
  feat <- pre[, c("gene", "contrast", "K", "expr_pct", "detect")]
  in_feat <- function(t) t[paste(t$gene, t$contrast) %in% paste(feat$gene, feat$contrast), , drop = FALSE]
  hox_in <- in_feat(hox_t_full); mk_in <- in_feat(mk_t_full)
  rp <- null_rates(draws$published, mk_in, de, pre)
  rh <- null_rates(draws$hox, hox_in, de, pre)
  ns <- read_tsv("data/null_summary_v2.tsv")
  chk <- function(g, v) {
    s <- ns[ns$group == g, ]
    check(length(v) == s$n_sets, g, ": set count mismatch (", length(v), " vs ", s$n_sets, ")")
    check(abs(mean(v, na.rm = TRUE) - s$null_mean) < 1e-9 &&
            abs(stats::sd(v, na.rm = TRUE) - s$null_sd) < 1e-9,
          g, ": mean/SD do not match null_summary")
    log_msg("null distribution ", g, ": 1000 rate sets reconcile with null_summary (mean=",
            signif(mean(v, na.rm = TRUE), 4), ")")
  }
  chk("published", rp); chk("hox", rh)
  write_tsv(data.frame(set = as.integer(names(rp)), rate = as.numeric(rp)),
            "data/fig2c_null_rates_published_v2.tsv")
  write_tsv(data.frame(set = as.integer(names(rh)), rate = as.numeric(rh)),
            "data/fig3a_null_rates_hox_v2.tsv")

  # ---- 7. Fig 4c: random-anchor metadata ----------------------------------------------------
  fp <- readRDS(proj_path("data", "loso_v2r2", "fold_prep.rds"))
  ra <- read_tsv("data/loso_v2r2/random_anchor_null.tsv")
  check(nrow(ra) == 100L, "Random anchors should be 100 sets")
  meta <- rbind(
    data.frame(metric = "n_random_sets", key = "all_folds", value = "100", stringsAsFactors = FALSE),
    data.frame(metric = "hox_survive", key = names(fp),
               value = vapply(fp, function(x) as.character(x$hox_survive), character(1)),
               stringsAsFactors = FALSE),
    data.frame(metric = "n_universe", key = names(fp),
               value = vapply(fp, function(x) as.character(x$n_universe), character(1)),
               stringsAsFactors = FALSE))
  write_tsv(meta, "data/fig4c_random_anchor_meta_v2.tsv")
  log_msg("Fig4c: random-anchor set sizes (HOX survivors per fold)=",
          paste(vapply(fp, function(x) as.character(x$hox_survive), character(1)), collapse = ","))

  # ---- 8. Outputs ------------------------------------------------------------------
  outs <- c("data/fig5a_vote_matrix_v2.tsv", "data/fig5a_vote_recompute_v2.tsv",
            "data/fig4b_confusion_v2.tsv", "data/fig4d_rule_genes_v2.tsv",
            "data/fig5d_platform_accuracy_v2.tsv", "data/fig2f_direction_consistency_v2.tsv",
            "data/fig2c_null_rates_published_v2.tsv", "data/fig3a_null_rates_hox_v2.tsv",
            "data/fig4c_random_anchor_meta_v2.tsv")
  for (f in outs) adopt(f)
  log_msg("92 done: 9 tables written + 4 caches cross-anchored; vote hard check 176/176 passed")
  invisible(NULL)
})
