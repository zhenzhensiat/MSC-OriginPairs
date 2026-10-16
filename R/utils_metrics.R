# =============================================================================
# utils_metrics.R —— classification metrics and confidence intervals
# Key points:
#   * primary metric: pooled sample-level balanced accuracy (mean of three recalls)
#   * no fold-level macro-F1: absent classes capped it at 0.667
#   * a class recall uses only true samples of that class
#   * NA predictions (no call / no evaluable rules) count as errors; rejection rate
#     is reported separately
# Prediction-table convention: data.frame(gsm, study, cluster_id, donor, true, pred, ...)
# =============================================================================

per_class_metrics <- function(true, pred, classes = CFG$classes) {
  do.call(rbind, lapply(classes, function(cl) {
    is_t <- true == cl
    is_p <- !is.na(pred) & pred == cl
    tp <- sum(is_t & is_p); fn <- sum(is_t & !is_p); fp <- sum(!is_t & is_p)
    rec  <- if (tp + fn > 0) tp / (tp + fn) else NA_real_
    prec <- if (tp + fp > 0) tp / (tp + fp) else NA_real_
    f1 <- if (!is.na(rec) && !is.na(prec) && (rec + prec) > 0) 2 * rec * prec / (rec + prec) else NA_real_
    data.frame(class = cl, n_true = tp + fn, n_pred = tp + fp,
               recall = rec, precision = prec, f1 = f1, stringsAsFactors = FALSE)
  }))
}

# average only over classes present in the sample (absent classes are not counted as 0)
balanced_accuracy <- function(true, pred, classes = CFG$classes) {
  m <- per_class_metrics(true, pred, classes)
  mean(m$recall, na.rm = TRUE)
}

study_accuracy <- function(df) {
  sp <- split(df, df$study)
  data.frame(study = names(sp),
             n = vapply(sp, nrow, integer(1)),
             acc = vapply(sp, function(d) mean(!is.na(d$pred) & d$pred == d$true), numeric(1)),
             stringsAsFactors = FALSE, row.names = NULL)
}

majority_label <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_character_)
  tb <- table(x)
  top <- names(tb)[tb == max(tb)]
  if (length(top) > 1) NA_character_ else top   # ties recorded as NA (counted as errors)
}

# for single-source studies: "per-dataset majority vote"; multi-source studies are grouped by study x source
majority_vote_accuracy <- function(df) {
  key <- paste(df$study, df$true, sep = "|")
  sp <- split(df, key)
  ok <- vapply(sp, function(d) identical(majority_label(d$pred), d$true[1]), logical(1))
  mean(ok)
}

# donor level: take the majority vote first within the samples of one donor (multiple passages), so repeated samples do not inflate the sample size
donor_level <- function(df) {
  key <- paste(df$study, df$donor, df$true, sep = "|")
  sp <- split(df, key)
  data.frame(study = vapply(sp, function(d) d$study[1], character(1)),
             cluster_id = vapply(sp, function(d) d$cluster_id[1], character(1)),
             donor = vapply(sp, function(d) d$donor[1], character(1)),
             true = vapply(sp, function(d) d$true[1], character(1)),
             pred = vapply(sp, function(d) majority_label(d$pred), character(1)),
             stringsAsFactors = FALSE, row.names = NULL)
}

summarize_predictions <- function(df, classes = CFG$classes) {
  check(all(c("study", "cluster_id", "donor", "true", "pred") %in% names(df)), "prediction table is missing required columns")
  sa <- study_accuracy(df)
  dl <- donor_level(df)
  list(
    n_samples   = nrow(df),
    n_studies   = length(unique(df$study)),
    ba          = balanced_accuracy(df$true, df$pred, classes),
    per_class   = per_class_metrics(df$true, df$pred, classes),
    study_acc_median = stats::median(sa$acc),
    study_acc_min    = min(sa$acc),
    majority_acc     = majority_vote_accuracy(df),
    donor_ba         = balanced_accuracy(dl$true, dl$pred, classes),
    abstain_rate     = mean(is.na(df$pred))
  )
}

flat_summary <- function(s, label) {
  pc <- s$per_class
  out <- data.frame(model = label, n_samples = s$n_samples, n_studies = s$n_studies,
                    ba = s$ba, study_acc_median = s$study_acc_median, study_acc_min = s$study_acc_min,
                    majority_acc = s$majority_acc, donor_ba = s$donor_ba, abstain_rate = s$abstain_rate,
                    stringsAsFactors = FALSE)
  for (i in seq_len(nrow(pc))) {
    out[[paste0("recall_", pc$class[i])]] <- pc$recall[i]
    out[[paste0("precision_", pc$class[i])]] <- pc$precision[i]
  }
  out
}

# ---- confidence intervals: bootstrap over independent clusters, stratified by cluster class composition (R2)------------
cluster_bootstrap <- function(df, stat_fun, R = CFG$boot_R, seed, cluster_col = "cluster_id",
                              classes = CFG$classes) {
  cl_comp <- tapply(df$true, df[[cluster_col]], function(v) paste(sort(unique(v)), collapse = "+"))
  strata <- split(names(cl_comp), unname(cl_comp))
  idx_by_cluster <- split(seq_len(nrow(df)), df[[cluster_col]])
  set.seed(seed)
  stats <- vapply(seq_len(R), function(b) {
    picked <- unlist(lapply(strata, function(s) s[sample.int(length(s), length(s), replace = TRUE)]),
                     use.names = FALSE)
    rows <- unlist(idx_by_cluster[picked], use.names = FALSE)
    stat_fun(df[rows, , drop = FALSE])
  }, numeric(1))
  n_cl <- vapply(classes, function(cl) sum(grepl(cl, cl_comp, fixed = TRUE)), integer(1))
  list(est = stat_fun(df),
       lo = unname(stats::quantile(stats, 0.025, na.rm = TRUE)),
       hi = unname(stats::quantile(stats, 0.975, na.rm = TRUE)),
       clusters_per_class = n_cl,
       unreliable = any(n_cl < CFG$min_clusters_reliable),
       boot = stats)
}

# ---- verdict wording---------------------------------------
h2_verdict <- function(s, boot, thr = CFG$h2, valid = TRUE) {
  if (!isTRUE(valid)) return("not testable")
  rec <- s$per_class$recall
  if (anyNA(rec)) return("not testable")          # a class is absent from the confirmation data
  ok <- s$ba >= thr$ba && boot$lo >= thr$ci_lo && all(rec >= thr$recall)
  if (isTRUE(ok)) "supported" else "not supported"
}

# empirical p: upper tail of the observed value in the null distribution
empirical_p <- function(obs, null) (sum(null >= obs, na.rm = TRUE) + 1) / (sum(!is.na(null)) + 1)
