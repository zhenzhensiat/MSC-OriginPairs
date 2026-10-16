# =============================================================================
# 21b_sens_phase2.R —— sensitivity analyses (post hoc).
# (1) published markers re-computed with the direction-agnostic convention
#     (matching the HOX convention; all other settings unchanged)
# (2) re-computation excluding GSE38947 (preselected-matrix bias); observed and
#     null sets handled identically.
# Inputs: data/within_study_de_v2.rds, data/pc_precomputed_v2.rds,
#   data/null_draws_v2.rds, curation/marker_list.tsv, hox_list.tsv
# Outputs: data/sensitivity_published_two_sided_v2.tsv,
#   data/sensitivity_no_gse38947_v2.tsv (written). Post hoc (run after the main
#   results were unblinded).
# =============================================================================
source("scripts/00_setup.R")

run_step("21b_sens_phase2", function() {
  de <- readRDS("data/within_study_de_v2.rds")
  pre <- readRDS("data/pc_precomputed_v2.rds")
  draws <- readRDS("data/null_draws_v2.rds")
  feat <- pre[, c("gene", "contrast", "K", "expr_pct", "detect")]
  in_feat <- function(t) t[paste(t$gene, t$contrast) %in% paste(feat$gene, feat$contrast), , drop = FALSE]
  hox <- read_tsv("curation/hox_list.tsv")$gene
  claims <- read_tsv("curation/marker_list.tsv")

  # ---- (i) published markers: two-sided definition (direction changed only; discrete null controls resampled with an independent seed)--------------
  mk2 <- expand_claims(claims)
  mk2$direction <- 0L
  mk2_in <- in_feat(mk2)
  res2 <- family_calls(pc_rows(mk2_in, de, pre))
  d2 <- sample_matched_null(mk2_in, feat, CFG$n_null_sets, derive_seed(21, 200000L))
  v2 <- null_rates(d2, mk2_in, de, pre)
  out2 <- data.frame(
    analysis = "published_two_sided", post_hoc = TRUE, discovery_exclusion = TRUE,
    n_rows = nrow(mk2_in), testable = sum(res2$testable), transfer = sum(res2$transfer),
    rate = transfer_rate(res2),
    null_mean = mean(v2, na.rm = TRUE), null_sd = stats::sd(v2, na.rm = TRUE),
    null_q95 = unname(stats::quantile(v2, 0.95, na.rm = TRUE)),
    empirical_p = empirical_p(transfer_rate(res2), v2),
    n_sets = length(v2), stringsAsFactors = FALSE)
  write_tsv(out2, "data/sensitivity_published_two_sided_v2.tsv")
  log_msg("(i) two-sided published markers: rate ", round(out2$rate, 4), "(testable ", out2$testable,
          "), null ", round(out2$null_mean, 4), "±", round(out2$null_sd, 4),
          ", p=", out2$empirical_p)

  # ---- (ii) exclude GSE38947 (exclusion applied to observed and existing null controls alike)--------------------------
  ex1 <- function(e) ifelse(nzchar(e), paste0(e, ";GSE38947"), "GSE38947")
  hox_t <- expand.grid(gene = hox, contrast = CFG$contrasts, stringsAsFactors = FALSE)
  hox_t$direction <- 0L; hox_t$exclude <- "GSE38947"
  hox_in <- in_feat(hox_t)
  mk1 <- expand_claims(claims); mk1$exclude <- ex1(mk1$exclude)
  mk1_in <- in_feat(mk1)
  r_hox <- family_calls(pc_rows(hox_in, de, pre))
  v_hox <- null_rates(draws$hox, hox_in, de, pre)
  r_mk <- family_calls(pc_rows(mk1_in, de, pre))
  v_mk <- null_rates(draws$published, mk1_in, de, pre)
  out3 <- do.call(rbind, list(
    data.frame(analysis = "hox_no_gse38947", post_hoc = TRUE, excluded = "GSE38947",
               n_rows = nrow(hox_in), testable = sum(r_hox$testable),
               transfer = sum(r_hox$transfer), rate = transfer_rate(r_hox),
               null_mean = mean(v_hox, na.rm = TRUE), null_sd = stats::sd(v_hox, na.rm = TRUE),
               null_q95 = unname(stats::quantile(v_hox, 0.95, na.rm = TRUE)),
               empirical_p = empirical_p(transfer_rate(r_hox), v_hox),
               n_sets = length(v_hox), stringsAsFactors = FALSE),
    data.frame(analysis = "published_directed_no_gse38947", post_hoc = TRUE, excluded = "GSE38947",
               n_rows = nrow(mk1_in), testable = sum(r_mk$testable),
               transfer = sum(r_mk$transfer), rate = transfer_rate(r_mk),
               null_mean = mean(v_mk, na.rm = TRUE), null_sd = stats::sd(v_mk, na.rm = TRUE),
               null_q95 = unname(stats::quantile(v_mk, 0.95, na.rm = TRUE)),
               empirical_p = empirical_p(transfer_rate(r_mk), v_mk),
               n_sets = length(v_mk), stringsAsFactors = FALSE)))
  write_tsv(out3, "data/sensitivity_no_gse38947_v2.tsv")
  log_msg("(ii) exclude GSE38947: HOX rate ", round(out3$rate[1], 4), "; published ", round(out3$rate[2], 4))

  log_msg("21b done: both sensitivity analyses tabulated and written (post-hoc annotated)")
})
