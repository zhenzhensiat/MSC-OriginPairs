# =============================================================================
# 51_confirm_secondary.R - secondary comparator: the
# elastic-net model on the validation cohort (descriptive).
# Precondition: data/confirm_complete exists (primary confirmation done); the
#   elastic-net model was written together with the primary model.
# Inputs: data/model_v2_glmnet.rds; validation-cohort study objects
# Outputs: data/confirm_secondary_glmnet_v2.tsv, data/confirm_secondary_summary_v2.tsv
# Rule: elastic net cannot handle missing genes; a study lacking model genes is
#   recorded as not evaluable; no imputation.
# =============================================================================
source("scripts/00_setup.R")
need_pkgs(c("glmnet"))

run_step("51_confirm_secondary", function() {
  check(file.exists(proj_path("data", "confirm_complete")), "primary confirmation not yet complete")
  gpath <- "data/model_v2_glmnet.rds"
  obj <- readRDS(proj_path(gpath))
  fit <- if (inherits(obj, "cv.glmnet")) obj else obj$fit
  check(inherits(fit, "cv.glmnet"), "glmnet object structure does not match assumptions; inspect with str(obj) first")
  cf <- stats::coef(fit, s = "lambda.1se")
  model_genes <- setdiff(rownames(cf[[1]]), "(Intercept)")   # all features at training time (not only non-zero coefficients)

  studies <- load_studies("C")
  meta <- bind_meta(studies)
  nz_genes <- unique(unlist(lapply(cf, function(b) rownames(b)[as.numeric(b) != 0])))
  nz_genes <- setdiff(nz_genes, "(Intercept)")
  log_msg("glmnet features: ", length(model_genes), ", of which non-zero-coefficient genes: ", length(nz_genes), "")

  rows <- lapply(studies, function(s) {
    m <- s$meta[s$meta$arm == "include" & s$meta$source %in% CFG$classes, , drop = FALSE]
    miss_nz <- setdiff(nz_genes, rownames(s$expr))
    if (length(miss_nz)) {   # missing non-zero-coefficient genes: not evaluable, no imputation performed
      return(data.frame(gsm = m$gsm, pred = NA_character_, status = "not_evaluable",
                        n_missing_nz_genes = length(miss_nz), stringsAsFactors = FALSE))
    }
    # Rank transform is computed on the intersection of model features and this study's genes (coverage usually >96%, close to the training definition, labelled approximate);
    # missing features are all zero-coefficient; 0.5 placeholders exist only to satisfy matrix dimensions and do not affect predictions
    have <- intersect(model_genes, rownames(s$expr))
    r <- rank_transform(s$expr[have, m$gsm, drop = FALSE])
    x <- matrix(0.5, nrow = length(m$gsm), ncol = length(model_genes),
                dimnames = list(m$gsm, model_genes))
    x[, have] <- t(r)
    check(all(setdiff(model_genes, have) %in% setdiff(model_genes, nz_genes)), "non-zero-coefficient gene among placeholder features")
    p <- as.character(stats::predict(fit, x, s = "lambda.1se", type = "class"))
    data.frame(gsm = m$gsm, pred = p, status = "ok", n_missing_nz_genes = 0L,
               rank_coverage = length(have) / length(model_genes), stringsAsFactors = FALSE)
  })
  rows <- lapply(rows, function(r) { if (is.null(r$rank_coverage)) r$rank_coverage <- NA_real_; r })
  pr <- attach_truth(do.call(rbind, rows), meta)
  write_tsv(pr, "data/confirm_secondary_glmnet_v2.tsv")

  ev <- pr[pr$status == "ok", , drop = FALSE]
  s_all <- summarize_predictions(pr)              # not-evaluable scored as error (conservative)
  s_ev <- summarize_predictions(ev)               # only evaluable studies
  sm <- rbind(flat_summary(s_all, "glmnet_all_conservative"), flat_summary(s_ev, "glmnet_evaluable_only"))
  sm$n_studies_not_evaluable <- length(unique(pr$study[pr$status != "ok"]))
  sm$descriptive_only <- TRUE
  write_tsv(sm, "data/confirm_secondary_summary_v2.tsv")
  invisible(sm)
})
