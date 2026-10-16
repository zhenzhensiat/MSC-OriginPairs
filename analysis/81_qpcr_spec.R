# =============================================================================
# 81_qpcr_spec.R - qPCR specification table for the 30 final rules / 56 genes.
# Purpose: a primer-design-ready gene list, flagging low-expression and
#   detection-risk genes and rules.
# Inputs: data/model_v2_rules.tsv; development- and validation-cohort study
#   objects (expression levels only; no re-prediction)
# Outputs: data/qpcr_spec_genes_v2.tsv (per gene), data/qpcr_spec_rules_v2.tsv
#   (per rule). Descriptive criteria: median within-study expression percentile
#   < 0.30 = low-expression risk; measured in < 90% of studies = coverage risk.
# =============================================================================
source("scripts/00_setup.R")

run_step("81_qpcr_spec", function() {
  rules <- read_tsv("data/model_v2_rules.tsv")
  genes <- sort(unique(c(rules$gene1, rules$gene2)))
  studies <- load_studies(c("D", "C"))

  per <- do.call(rbind, lapply(studies, function(s) {
    m <- s$meta[s$meta$arm == "include" & s$meta$source %in% CFG$classes, , drop = FALSE]
    ave <- rowMeans(s$expr[, m$gsm, drop = FALSE])
    pct <- (rank(ave) - 0.5) / length(ave)
    data.frame(gse = s$gse, layer = s$layer, gene = genes, present = genes %in% names(pct),
               expr_pct = unname(pct[genes]), stringsAsFactors = FALSE)
  }))
  gs <- do.call(rbind, lapply(split(per, per$gene), function(d) {
    data.frame(gene = d$gene[1], frac_studies_present = mean(d$present),
               median_expr_pct = stats::median(d$expr_pct, na.rm = TRUE),
               q10_expr_pct = unname(stats::quantile(d$expr_pct, 0.10, na.rm = TRUE)),
               stringsAsFactors = FALSE)
  }))
  gs$low_expression_risk <- gs$median_expr_pct < 0.30
  gs$coverage_risk <- gs$frac_studies_present < 0.90
  write_tsv(gs, "data/qpcr_spec_genes_v2.tsv")

  r <- rules
  r$g1_pct <- gs$median_expr_pct[match(r$gene1, gs$gene)]
  r$g2_pct <- gs$median_expr_pct[match(r$gene2, gs$gene)]
  r$min_pct <- pmin(r$g1_pct, r$g2_pct)
  r$any_risk <- gs$low_expression_risk[match(r$gene1, gs$gene)] | gs$low_expression_risk[match(r$gene2, gs$gene)] |
                gs$coverage_risk[match(r$gene1, gs$gene)] | gs$coverage_risk[match(r$gene2, gs$gene)]
  r <- r[order(r$class, -r$separation), ]
  write_tsv(r, "data/qpcr_spec_rules_v2.tsv")
  log_msg("low-expression risk among the 56 genes: ", sum(gs$low_expression_risk), "; platform coverage risk: ", sum(gs$coverage_risk), "")
  invisible(r)
})
