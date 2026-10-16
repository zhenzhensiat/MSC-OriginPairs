# =============================================================================
# 90_figures.R - assembly of main figures 1-6 (reads result tables only; every
# plotted number traces to a table; figure text is English).
# Written-out examples: the marker-audit and confirmation figure inputs; the
# remaining figure input tables follow the same style.
# =============================================================================
source("scripts/00_setup.R")
source(proj_path("R", "utils_plot.R"), encoding = "UTF-8")
need_pkgs(c("ggplot2", "patchwork"))
library(ggplot2)

# ---- Fig 2: marker reproducibility audit (summarized by evidence type and contrast type; 2026-10-13 instruction: do not name papers)----
fig2_marker_audit <- function() {
  bs <- read_tsv("data/marker_by_stratum_v2.tsv")
  panel_lab <- c(assay_orig = "Evidence type", contrast = "Contrast")
  long <- do.call(rbind, lapply(unique(bs$stratum), function(v) {
    x <- bs[bs$stratum == v, , drop = FALSE]
    do.call(rbind, lapply(c("claimed", "reversed", "not_replicated"), function(k) {
      data.frame(panel = panel_lab[[v]], level = x$level, n_rows = x$n_rows, category = k,
                 frac = x[[k]] / x$n_rows, stringsAsFactors = FALSE)
    }))
  }))
  long$category <- factor(long$category, levels = c("claimed", "reversed", "not_replicated"),
                          labels = c("Replicated (reported direction)", "Discordant direction", "Not replicated"))
  nlab <- unique(long[, c("panel", "level", "n_rows")])
  ggplot(long, aes(x = level, y = frac, fill = category)) +
    geom_col(width = 0.75) +
    geom_text(data = nlab, aes(x = level, y = 1.03, label = paste0("n=", n_rows)),
              inherit.aes = FALSE, size = 2.6) +
    scale_fill_manual(values = c(PAL$claimed, PAL$reversed, PAL$not_rep), name = NULL) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1.08)) +
    facet_wrap(~panel, scales = "free_x") +
    labs(x = NULL, y = "Fraction of testable claims",
         title = "Replication status of published source-marker claims") +
    theme_msc() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

# ---- Fig 5: independent confirmation (per-study accuracy + permutation null)-------------------------------------
fig5_confirmation <- function() {
  cp <- read_tsv("data/confirm_predictions_v2.tsv")
  cr <- read_tsv("data/confirm_results_v2.tsv")
  nl <- read_tsv("data/confirm_null_v2.tsv")
  st <- do.call(rbind, lapply(split(cp, cp$study), function(d) {
    data.frame(study = d$study[1], true = d$true[1], n = nrow(d),
               acc = mean(!is.na(d$pred) & d$pred == d$true),
               abstain = mean(is.na(d$pred)), stringsAsFactors = FALSE)
  }))
  st$study <- factor(st$study, levels = st$study[order(st$true, -st$n)])
  pa <- ggplot(st, aes(x = study, y = acc, size = n, colour = true)) +
    geom_point() +
    geom_point(data = st[st$abstain > 0, ], aes(x = study, y = acc), shape = 4, colour = PAL$accent, size = 3) +
    scale_colour_manual(values = c(AD = PAL$AD, BM = PAL$BM, UC = PAL$UC), name = "Source") +
    scale_size_area(max_size = 5, name = "Samples") +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(x = NULL, y = "Within-study accuracy", title = "Independent confirmation by study (x = all samples abstained)") +
    theme_msc() + theme(axis.text.x = element_text(angle = 60, hjust = 1))
  ba <- as.numeric(cr$value[cr$metric == "balanced_accuracy"])
  pb <- ggplot(nl, aes(x = perm_ba)) +
    geom_histogram(bins = 30, fill = PAL$null) +
    geom_vline(xintercept = ba, colour = PAL$accent, linewidth = 0.8) +
    annotate("text", x = ba, y = Inf, vjust = 1.5, hjust = 1.05, size = 2.8,
             label = sprintf("Observed BA = %.3f", ba)) +
    labs(x = "Balanced accuracy under permuted labels", y = "Count", title = "Chance level") +
    theme_msc()
  patchwork::wrap_plots(pa, pb, widths = c(2.2, 1))
}

# ---- Input table inventory for the remaining main figures (written in the same style as above)--------------------------------
# Fig 1 study design: census pipeline counts, layer_allocation_v2b, c_layer_studies_v2, platform composition
# Fig 3 HOX and families: hox_transfer_v2, null_summary_v2, family_ranking_v2 (post-fix version), HOX subcluster direction matrix
# Fig 4 variance and classification: variance_partition_v2, loso_v2r2/loso_summary, random_anchor_null (state test limitations)
# Fig 6 stability and application: pair_stability_*_v2, arm_robustness_summary_v2, openset_summary_v2, qpcr_spec_rules_v2

run_step("90_figures", function() {
  save_fig(fig2_marker_audit(), "Fig2_marker_audit", 180, 100)
  save_fig(fig5_confirmation(), "Fig5_confirmation", 180, 90)
  invisible(TRUE)
})
