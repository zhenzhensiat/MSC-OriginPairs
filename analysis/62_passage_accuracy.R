# =============================================================================
# 62_passage_accuracy.R - whole-sample accuracy across passages (post hoc,
# exploratory). Grouping rule fixed before running: primary = earliest vs
# latest passage within each study; secondary = absolute bands.
# Question: does the final model still classify whole samples correctly across
# passages (rule-direction stability is checked separately by 60).
# Inputs: stability=TRUE study objects; final model + tau
# Outputs: data/passage_accuracy_pred_v2.tsv (incl. donor/passage/margin for the
#   Fig 6d trajectories); data/passage_accuracy_summary_v2.tsv (written)
# =============================================================================
source("scripts/00_setup.R")
need_pkgs(c("multiclassPairs", "switchBox"))

run_step("62_passage_accuracy", function() {
  clf <- readRDS(proj_path("data", "model_v2.rds"))
  tau <- readRDS(proj_path("data", "model_v2_tau.rds"))
  studies <- load_studies(stability_only = TRUE)
  check(length(studies) == 10L, "stability studies should be 10, actual ", length(studies))

  rows <- list()
  for (s in studies) {
    m <- s$meta[s$meta$arm == "include" & s$meta$source %in% CFG$classes & !is.na(s$meta$passage), , drop = FALSE]
    if (!nrow(m)) { log_msg(s$gse, ": no include samples with passage, skipping"); next }
    p <- apply_reject(ktsp_predict(clf, s$expr[, m$gsm, drop = FALSE]), tau)
    j <- match(p$gsm, m$gsm)
    check(!anyNA(j), s$gse, ": predictions do not match meta")
    pm <- m$passage[j]
    band <- ifelse(pm == min(pm), "early", ifelse(pm == max(pm), "late", "mid"))
    band_abs <- ifelse(pm <= 3, "P_le3", ifelse(pm >= 8, "P_ge8", "P_4_7"))
    rows[[s$gse]] <- data.frame(
      study = s$gse, gsm = p$gsm, donor = m$donor[j], passage = pm,
      source = m$source[j], pred = p$pred, margin = p$margin,
      rejected = p$rejected %in% c(TRUE, "TRUE"),
      band_research = band, band_abs = band_abs, stringsAsFactors = FALSE)
    log_msg(s$gse, ": n=", nrow(m), "; passage range ", min(pm), "-", max(pm),
            " (early ", sum(band == "early"), "/late ", sum(band == "late"), ")")
  }
  pr <- do.call(rbind, rows)
  pr$correct <- !is.na(pr$pred) & pr$pred == pr$source
  write_tsv(pr, "data/passage_accuracy_pred_v2.tsv")

  summ <- function(d, v) {
    do.call(rbind, lapply(split(d, d[[v]]), function(x) {
      data.frame(band_type = v, band = x[[v]][1], n = nrow(x),
                 n_studies = length(unique(x$study)),
                 accuracy = mean(x$correct), abstain = mean(is.na(x$pred)),
                 rejected = mean(x$rejected),
                 median_margin = stats::median(x$margin, na.rm = TRUE),
                 post_hoc = TRUE, stringsAsFactors = FALSE)
    }))
  }
  sm <- rbind(summ(pr, "band_research"), summ(pr, "band_abs"))
  # Per-study primary-definition summary (for Fig 6d)
  smr <- do.call(rbind, lapply(split(pr, pr$study), function(x)
    data.frame(study = x$study[1], n = nrow(x),
               accuracy_early = mean(x$correct[x$band_research == "early"]),
               accuracy_late = mean(x$correct[x$band_research == "late"]),
               med_margin_early = stats::median(x$margin[x$band_research == "early"], na.rm = TRUE),
               med_margin_late = stats::median(x$margin[x$band_research == "late"], na.rm = TRUE),
               post_hoc = TRUE, stringsAsFactors = FALSE)))
  write_tsv(smr, "data/passage_accuracy_by_study_v2.tsv")
  write_tsv(sm, "data/passage_accuracy_summary_v2.tsv")
  log_msg("primary definition: early ", round(100 * mean(pr$correct[pr$band_research == "early"]), 1),
          "% (n=", sum(pr$band_research == "early"), "); late ",
          round(100 * mean(pr$correct[pr$band_research == "late"]), 1),
          "% (n=", sum(pr$band_research == "late"), "); all ",
          round(100 * mean(pr$correct), 1), "% (n=", nrow(pr), ")")
  invisible(sm)
})
