# =============================================================================
# 60_stability.R - passage stability of the final rules.
# Inputs: stability = TRUE studies (metadata with donor and numeric passage);
#   data/model_v2_rules.tsv
# Criterion: a rule is "stable" in a source if the same-sign donor fraction is
#   >= 0.90 with >=3 evaluable donors.
# Outputs: pair_stability_passage_v2.tsv (rule x study),
#   pair_stability_pooled_v2.tsv (rule-level summary),
#   pair_stability_lmm_v2.tsv (mixed model slope ~ passage + (1|donor)).
# Single-donor data receive descriptive trajectories only.
# =============================================================================
source("scripts/00_setup.R")
need_pkgs(c("lme4"))

run_step("60_stability", function() {
  rules <- read_tsv("data/model_v2_rules.tsv")
  studies <- load_studies(stability_only = TRUE)

  per <- list(); lmm_rows <- list()
  for (s in studies) {
    m <- s$meta[s$meta$arm == "include" & s$meta$source %in% CFG$classes & !is.na(s$meta$passage), , drop = FALSE]
    if (nrow(m) == 0) { log_msg(s$gse, ": no samples with passage information, skipping"); next }
    donors <- split(m, m$donor)
    donors <- donors[vapply(donors, function(d) length(unique(d$passage)) >= 2, logical(1))]
    if (!length(donors)) { log_msg(s$gse, ": no donor has >=2 passages, skipping"); next }
    src <- unique(m$source)
    check(length(src) == 1, s$gse, ": stability studies should have a single source")

    for (i in seq_len(nrow(rules))) {
      g1 <- rules$gene1[i]; g2 <- rules$gene2[i]
      if (!all(c(g1, g2) %in% rownames(s$expr))) {
        per[[length(per) + 1L]] <- data.frame(study = s$gse, source = src, class = rules$class[i],
          gene1 = g1, gene2 = g2, n_donors = 0L, same_sign_prop = NA_real_,
          expected_dir_early = NA_real_, expected_dir_late = NA_real_, status = "gene_missing",
          stringsAsFactors = FALSE)
        next
      }
      d_all <- s$expr[g1, m$gsm] - s$expr[g2, m$gsm]
      ed <- vapply(donors, function(d) {
        early <- d$gsm[which.min(d$passage)]; late <- d$gsm[which.max(d$passage)]
        c(s$expr[g1, early] - s$expr[g2, early], s$expr[g1, late] - s$expr[g2, late])
      }, numeric(2))
      same <- sign(ed[1, ]) == sign(ed[2, ]) & ed[1, ] != 0
      # whether the rule's 'class direction' holds in this source's samples (meaningful only when the source equals the rule's class)
      want <- if (rules$class_direction[i] == "gene1>gene2") 1 else -1
      per[[length(per) + 1L]] <- data.frame(study = s$gse, source = src, class = rules$class[i],
        gene1 = g1, gene2 = g2, n_donors = ncol(ed), same_sign_prop = mean(same),
        expected_dir_early = mean(sign(ed[1, ]) == want), expected_dir_late = mean(sign(ed[2, ]) == want),
        status = "ok", stringsAsFactors = FALSE)

      dd <- data.frame(d = d_all, passage = m$passage, donor = m$donor)
      if (length(unique(dd$donor)) >= 3) {
        fit <- tryCatch(lme4::lmer(d ~ passage + (1 | donor), data = dd, REML = TRUE),
                        error = function(e) NULL)
        if (!is.null(fit)) {
          est <- lme4::fixef(fit)[["passage"]]
          se <- sqrt(diag(as.matrix(stats::vcov(fit))))[2]
          lmm_rows[[length(lmm_rows) + 1L]] <- data.frame(study = s$gse, class = rules$class[i], gene1 = g1, gene2 = g2,
            slope = est, lo = est - 1.96 * se, hi = est + 1.96 * se,
            drift = (est - 1.96 * se) > 0 | (est + 1.96 * se) < 0, stringsAsFactors = FALSE)
        }
      }
    }
  }
  check(length(per) > 0, "no usable stability data")
  per <- do.call(rbind, per)
  write_tsv(per, "data/pair_stability_passage_v2.tsv")

  ok <- per[per$status == "ok", ]
  key <- paste(ok$class, ok$gene1, ok$gene2, ok$source, sep = "|")
  pooled <- do.call(rbind, lapply(split(ok, key), function(d) {
    n <- sum(d$n_donors)
    prop <- sum(d$same_sign_prop * d$n_donors) / n
    data.frame(class = d$class[1], gene1 = d$gene1[1], gene2 = d$gene2[1], source = d$source[1],
               n_studies = nrow(d), n_donors = n, same_sign_prop = prop,
               stable = n >= CFG$stability$min_donors & prop >= CFG$stability$min_same_sign,
               evaluable = n >= CFG$stability$min_donors, stringsAsFactors = FALSE)
  }))
  write_tsv(pooled, "data/pair_stability_pooled_v2.tsv")
  if (length(lmm_rows)) write_tsv(do.call(rbind, lmm_rows), "data/pair_stability_lmm_v2.tsv")

  ev <- pooled[pooled$evaluable, ]
  log_msg("evaluable rule x source combinations: ", nrow(ev), ", of which stable: ", sum(ev$stable), " (",
          round(100 * mean(ev$stable), 1), "%)")
  invisible(pooled)
})
