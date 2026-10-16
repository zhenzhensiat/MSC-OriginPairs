# =============================================================================
# 43_mixture.R —— simulated admixture detection limit (development layer;
# uses models trained without the held-out group to prevent leakage).
# Assumption: two samples mixed proportionally on the linear scale, then log2
# (computer simulation; stated as such in the manuscript).
# Inputs: data/loso_v2r2/models/genome_main_<group>.rds, data/loso_v2r2/tau_genome_main.rds
# Outputs: data/mixture_v2.tsv (per combination), data/mixture_curve_v2.tsv
# =============================================================================
source("scripts/00_setup.R")
need_pkgs(c("multiclassPairs"))

run_step("43_mixture", function() {
  assert_multiclasspairs_api()
  studies <- load_studies("D")
  tau <- readRDS(proj_path("data", "loso_v2r2", "tau_genome_main.rds"))
  fracs <- c(0, 0.05, 0.10, 0.20, 0.30, 0.50)
  max_pairs <- 20L
  rows <- list()

  for (s in studies) {
    if (s$cluster_id %in% CFG$loso_exclude_clusters) next   # excluded clusters have no held-out model (D2: pre-screened matrices)
    m <- s$meta[s$meta$arm == "include" & s$meta$source %in% CFG$classes, , drop = FALSE]
    present <- intersect(CFG$classes, unique(m$source))
    if (length(present) < 2) next
    mf <- proj_path("data", "loso_v2r2", "models", paste0("genome_main_", s$cluster_id, ".rds"))
    check(file.exists(mf), "missing the model trained with this cluster held out: ", mf, " (run 40 first)")
    clf <- readRDS(mf)
    pairs <- t(utils::combn(present, 2))
    for (pi in seq_len(nrow(pairs))) for (side in 1:2) {
      major <- pairs[pi, side]; minor <- pairs[pi, 3 - side]
      grid <- expand.grid(a = m$gsm[m$source == major], b = m$gsm[m$source == minor], stringsAsFactors = FALSE)
      set.seed(derive_seed(43, length(rows) + 1L))
      if (nrow(grid) > max_pairs) grid <- grid[sample.int(nrow(grid), max_pairs), , drop = FALSE]
      for (f in fracs) {
        mixed <- vapply(seq_len(nrow(grid)), function(i) mix_log2(s$expr[, grid$a[i]], s$expr[, grid$b[i]], 1 - f),
                        numeric(nrow(s$expr)))
        mixed <- matrix(mixed, nrow = nrow(s$expr), dimnames = list(rownames(s$expr), paste0("mix", seq_len(nrow(grid)))))
        pr <- apply_reject(ktsp_predict(clf, mixed), tau)
        rows[[length(rows) + 1L]] <- data.frame(
          study = s$gse, major = major, minor = minor, minor_frac = f, n = nrow(pr),
          still_major = mean(pr$pred == major, na.rm = TRUE),
          called_minor = mean(pr$pred == minor, na.rm = TRUE),
          rejected = mean(pr$rejected),
          flagged = mean(pr$rejected | is.na(pr$pred) | pr$pred != major),
          stringsAsFactors = FALSE)
      }
    }
  }
  check(length(rows) > 0, "no multi-source studies in the development layer for admixture simulation")
  mix <- do.call(rbind, rows)
  write_tsv(mix, "data/mixture_v2.tsv")

  curve <- stats::aggregate(cbind(flagged, still_major, rejected) ~ minor_frac, data = mix, FUN = mean)
  hit <- curve$minor_frac[curve$minor_frac > 0 & curve$flagged >= 0.80]
  curve$detection_limit_80 <- if (length(hit)) min(hit) else NA_real_
  write_tsv(curve, "data/mixture_curve_v2.tsv")
  log_msg("admixture 80% detection limit (minor component fraction): ", curve$detection_limit_80[1],
          "; false-positive rate at minor_frac = 0 = ", round(curve$flagged[curve$minor_frac == 0], 3))
  invisible(curve)
})
