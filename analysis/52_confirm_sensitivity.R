# =============================================================================
# 52_confirm_sensitivity.R - confirmation sensitivity analyses (reads the
# prediction tables only; nothing is re-predicted).
# Inputs: data/confirm_predictions_v2.tsv, data/confirm_results_v2.tsv
# Output: data/confirm_sensitivity_v2.tsv (analysis, metric, value, post_hoc)
# Content: (a) main result recomputed and asserted against the result table;
#   (b) excluding GSE24871; (c) exact binomial bounds for evaluable samples and
#   groups; (d) donor level; (e) leave-one-group-out influence; (f) rejection
#   rate; (g) by data scale.
# =============================================================================
source("scripts/00_setup.R")

run_step("52_confirm_sensitivity", function() {
  df <- read_tsv("data/confirm_predictions_v2.tsv")
  fr <- read_tsv("data/confirm_results_v2.tsv")
  out <- list()
  add <- function(analysis, metric, value, post_hoc = TRUE) {
    out[[length(out) + 1L]] <<- data.frame(analysis = analysis, metric = metric,
                                           value = as.character(value), post_hoc = post_hoc,
                                           stringsAsFactors = FALSE)
  }

  # (a) Recompute cross-check
  s <- summarize_predictions(df)
  ba_ref <- as.numeric(fr$value[fr$metric == "balanced_accuracy"])
  check(abs(s$ba - ba_ref) < 1e-9, "recomputed BA does not match the result table")
  add("main_recheck", "balanced_accuracy", s$ba, FALSE)

  # (b) Exclude GSE24871 (technical full abstention)
  d2 <- df[df$study != "GSE24871", , drop = FALSE]
  s2 <- summarize_predictions(d2)
  add("excl_GSE24871", "balanced_accuracy", s2$ba)
  for (i in seq_len(nrow(s2$per_class))) add("excl_GSE24871", paste0("recall_", s2$per_class$class[i]), s2$per_class$recall[i])

  # (c) Evaluable samples and clusters
  ev <- df[!is.na(df$pred), , drop = FALSE]
  k <- sum(ev$pred == ev$true); n <- nrow(ev)
  bt <- stats::binom.test(k, n)
  add("evaluable_samples", "correct/n", paste0(k, "/", n))
  add("evaluable_samples", "exact_ci_lo", bt$conf.int[1])
  cl_ok <- tapply(ev$pred == ev$true, ev$cluster_id, all)
  bc <- stats::binom.test(sum(cl_ok), length(cl_ok))
  add("evaluable_clusters", "all_correct/n", paste0(sum(cl_ok), "/", length(cl_ok)))
  add("evaluable_clusters", "exact_ci_lo", bc$conf.int[1])

  # (d) Donor level
  dl <- donor_level(df)
  add("donor_level", "balanced_accuracy", balanced_accuracy(dl$true, dl$pred))

  # (e) Leave-one-cluster-out
  loo <- vapply(sort(unique(df$cluster_id)), function(cl) {
    d <- df[df$cluster_id != cl, , drop = FALSE]
    balanced_accuracy(d$true, d$pred)
  }, numeric(1))
  add("leave_one_cluster_out", "ba_min", min(loo))
  add("leave_one_cluster_out", "ba_max", max(loo))
  add("leave_one_cluster_out", "most_influential", names(loo)[which.min(loo)])

  # (f) Rejection (denominator = samples with predictions)
  rj <- df$rejected %in% c(TRUE, "TRUE")
  add("openset", "rejection_rate_among_predicted", sum(rj & !is.na(df$pred)) / sum(!is.na(df$pred)))
  add("openset", "rejected_but_correct", sum(rj & !is.na(df$pred) & df$pred == df$true))

  # (g) By scale
  for (sc in sort(unique(df$scale))) {
    d <- df[df$scale == sc, , drop = FALSE]
    add(paste0("scale_", sc), "balanced_accuracy", balanced_accuracy(d$true, d$pred))
  }

  res <- do.call(rbind, out)
  write_tsv(res, "data/confirm_sensitivity_v2.tsv")
  invisible(res)
})
