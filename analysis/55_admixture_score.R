# =============================================================================
# 55_admixture_score.R - detecting admixture with a continuous score
# (exploratory; complements the call-flip criterion of script 43, which needs
# ~50% minor component). Uses the normalised score of the minor source as an
# admixture index.
# Inputs: per-group models trained without the held-out group;
#   development-layer multi-source studies (excluding the preselected-matrix
#   group)
# Outputs: data/admixture_scores_v2.tsv (per simulated sample),
#   data/admixture_roc_v2.tsv
# Assumption: linear-scale mixing then log2 (as in 43; computer simulation).
# =============================================================================
source("scripts/00_setup.R")
need_pkgs(c("multiclassPairs", "switchBox"))

auc_rank <- function(pos, neg) {
  r <- rank(c(pos, neg))
  (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) / (length(pos) * length(neg))
}

run_step("55_admixture_score", function() {
  studies <- load_studies("D")
  studies <- Filter(function(s) !s$cluster_id %in% CFG$loso_exclude_clusters, studies)
  fracs <- c(0, 0.05, 0.10, 0.20, 0.30, 0.50)
  max_pairs <- 20L
  rows <- list()
  for (s in studies) {
    m <- s$meta[s$meta$arm == "include" & s$meta$source %in% CFG$classes, , drop = FALSE]
    present <- intersect(CFG$classes, unique(m$source))
    if (length(present) < 2) next
    mf <- proj_path("data", "loso_v2r2", "models", paste0("genome_main_", s$cluster_id, ".rds"))
    check(file.exists(mf), "missing the leave-one-cluster-out model: ", mf)
    clf <- readRDS(mf)
    pairs <- t(utils::combn(present, 2))
    for (pi in seq_len(nrow(pairs))) for (side in 1:2) {
      major <- pairs[pi, side]; minor <- pairs[pi, 3 - side]
      grid <- expand.grid(a = m$gsm[m$source == major], b = m$gsm[m$source == minor], stringsAsFactors = FALSE)
      set.seed(derive_seed(55, length(rows) + 1L))
      if (nrow(grid) > max_pairs) grid <- grid[sample.int(nrow(grid), max_pairs), , drop = FALSE]
      for (f in fracs) {
        mixed <- vapply(seq_len(nrow(grid)), function(i) mix_log2(s$expr[, grid$a[i]], s$expr[, grid$b[i]], 1 - f),
                        numeric(nrow(s$expr)))
        mixed <- matrix(mixed, nrow = nrow(s$expr), dimnames = list(rownames(s$expr), paste0("mix", seq_len(nrow(grid)))))
        p <- ktsp_predict(clf, mixed)
        sc <- as.matrix(p[, paste0("score_", CFG$classes)])
        sc <- sc / pmax(rowSums(sc, na.rm = TRUE), .Machine$double.eps)
        rows[[length(rows) + 1L]] <- data.frame(
          study = s$gse, major = major, minor = minor, minor_frac = f, pair = seq_len(nrow(grid)),
          minor_index = sc[, paste0("score_", minor)], major_index = sc[, paste0("score_", major)],
          pred = p$pred, stringsAsFactors = FALSE)
      }
    }
  }
  sc <- do.call(rbind, rows)
  write_tsv(sc, "data/admixture_scores_v2.tsv")

  # Within each (study, major component, minor component) combination: AUC of pure samples (f=0) vs f>0; threshold = 95th percentile of the minor index among pure samples
  roc <- do.call(rbind, lapply(split(sc, paste(sc$study, sc$major, sc$minor)), function(d) {
    pure <- d$minor_index[d$minor_frac == 0]
    thr <- unname(stats::quantile(pure, 0.95, na.rm = TRUE))
    do.call(rbind, lapply(setdiff(sort(unique(d$minor_frac)), 0), function(f) {
      mix <- d$minor_index[d$minor_frac == f]
      data.frame(study = d$study[1], major = d$major[1], minor = d$minor[1], minor_frac = f,
                 auc = auc_rank(mix, pure), sens_at_spec95 = mean(mix > thr), stringsAsFactors = FALSE)
    }))
  }))
  agg <- stats::aggregate(cbind(auc, sens_at_spec95) ~ minor_frac, data = roc, FUN = stats::median)
  agg$post_hoc <- TRUE
  write_tsv(roc, "data/admixture_roc_by_combo_v2.tsv")
  write_tsv(agg, "data/admixture_roc_v2.tsv")
  log_msg("admixture index: ", paste(sprintf("f=%.2f AUC=%.2f", agg$minor_frac, agg$auc), collapse = "; "))
  invisible(agg)
})
