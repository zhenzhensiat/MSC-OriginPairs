# =============================================================================
# 22_marker_reproducibility.R —— three-way classification, per-paper summaries
# and predictors of reproducibility for the published markers.
# Inputs: data/marker_transfer_v2.tsv (directed results),
#   data/within_study_de_v2.rds, data/pc_precomputed_v2.rds,
#   curation/marker_list.tsv,
#   curation/marker_claim_attributes_v2.tsv (claim_id, n_orig_min, assay_orig,
#   medium_orig, passage_orig)
# Outputs: data/marker_three_class_v2r.tsv, data/marker_by_paper_v2.tsv,
#   data/marker_predictors_v2.tsv
# Notes: three-way convention as elsewhere; results are asserted against the
#   reference 36/16/47 counts first. The predictor analysis is exploratory
#   (post_hoc), small (99 rows, 8 papers), descriptive only.
# =============================================================================
source("scripts/00_setup.R")
need_pkgs(c("lme4"))

run_step("22_marker_reproducibility", function() {
  tm  <- read_tsv("data/marker_transfer_v2.tsv")
  de  <- readRDS(proj_path("data", "within_study_de_v2.rds"))
  pre <- readRDS(proj_path("data", "pc_precomputed_v2.rds"))
  ml  <- read_tsv("curation/marker_list.tsv")
  att <- read_tsv("curation/marker_claim_attributes_v2.tsv")

  tm$exclude[is.na(tm$exclude)] <- ""
  tm <- tm[tm$testable %in% c(TRUE, "TRUE"), , drop = FALSE]
  tm$transfer <- tm$transfer %in% c(TRUE, "TRUE")

  # ---- 1. two-sided test on the same rows (within-group BH) to identify 'reproduced opposite' ----------------------
  two <- family_calls(pc_rows(data.frame(gene = tm$gene, contrast = tm$contrast,
                                         direction = 0L, exclude = tm$exclude,
                                         stringsAsFactors = FALSE), de, pre))
  check(nrow(two) == nrow(tm), "two-sided result row counts differ")
  tm$two_q <- two$pc_q
  tm$two_dir <- two$dir_used
  tm$category <- ifelse(tm$transfer, "claimed",
                 ifelse(two$transfer & two$dir_used == -tm$direction, "reversed", "not_replicated"))
  tab <- table(factor(tm$category, levels = c("claimed", "reversed", "not_replicated")))
  log_msg("three-way classification: ", paste(names(tab), tab, sep = "=", collapse = ", "))
  check(identical(as.integer(tab), c(36L, 16L, 47L)),
        "three-way classification disagrees with the reference 7a result (36/16/47): verify definitions before proceeding")

  # ---- 2. detail and summary (2026-10-13 definition adjustment: no paper named in the main text or main figures)------------------
  # itemized detail used as a supplementary table only (data only: claim, source PMID, effect size, class; no commentary)
  tm <- safe_merge(tm, ml[, c("claim_id", "pmid", "evidence")], by = "claim_id", type = "m:1")
  pre_feat <- pre[, c("gene", "contrast", "expr_pct", "detect")]
  tm <- safe_merge(tm, pre_feat, by = c("gene", "contrast"), type = "m:1")
  write_tsv(tm, "data/marker_three_class_v2r.tsv")

  # per-paper summary: internal check only, not for the manuscript (directive 2026-10-13)
  by_paper <- do.call(rbind, lapply(split(tm, tm$pmid), function(d) {
    data.frame(pmid = d$pmid[1], n_rows = nrow(d),
               claimed = sum(d$category == "claimed"),
               reversed = sum(d$category == "reversed"),
               not_replicated = sum(d$category == "not_replicated"),
               median_abs_logFC = stats::median(abs(d$median_logFC)),
               stringsAsFactors = FALSE)
  }))
  by_paper$frac_claimed <- by_paper$claimed / by_paper$n_rows
  by_paper$frac_reversed <- by_paper$reversed / by_paper$n_rows
  by_paper$internal_only <- TRUE   # internal-check flag; must not be used in the manuscript text or main figures
  write_tsv(by_paper[order(-by_paper$frac_reversed), ], "data/marker_by_paper_v2.tsv")

  # ---- 2c. summary by assay and contrast (manuscript presentation definition; data source for Fig 2)----
  d <- safe_merge(tm, att, by = "claim_id", type = "m:1")
  by_stratum <- do.call(rbind, lapply(c("assay_orig", "contrast"), function(v) {
    do.call(rbind, lapply(sort(unique(d[[v]])), function(lv) {
      x <- d[d[[v]] == lv, , drop = FALSE]
      data.frame(stratum = v, level = lv, n_rows = nrow(x),
                 claimed = sum(x$category == "claimed"),
                 reversed = sum(x$category == "reversed"),
                 not_replicated = sum(x$category == "not_replicated"),
                 frac_claimed = mean(x$category == "claimed"),
                 stringsAsFactors = FALSE)
    }))
  }))
  write_tsv(by_stratum, "data/marker_by_stratum_v2.tsv")

  # ---- 3. predictors of reproducibility (exploratory)---------------------------------------------
  d$rep <- as.integer(d$category == "claimed")
  d$log2_n <- log2(pmax(as.numeric(d$n_orig_min), 1))
  uni <- do.call(rbind, lapply(c("assay_orig", "medium_orig", "contrast"), function(v) {
    do.call(rbind, lapply(split(d, d[[v]]), function(x) {
      data.frame(variable = v, level = x[[v]][1], n = nrow(x), rep_rate = mean(x$rep),
                 reversed_rate = mean(x$category == "reversed"), stringsAsFactors = FALSE)
    }))
  }))
  fit_txt <- tryCatch({
    f <- lme4::glmer(rep ~ log2_n + expr_pct + K + (1 | pmid), data = d, family = stats::binomial())
    co <- summary(f)$coefficients
    paste(sprintf("%s: %.3f (SE %.3f)", rownames(co), co[, 1], co[, 2]), collapse = "; ")
  }, error = function(e) paste("glmer did not converge: ", conditionMessage(e)))
  uni$model_summary <- fit_txt
  uni$post_hoc <- TRUE
  write_tsv(uni, "data/marker_predictors_v2.tsv")

  invisible(by_paper)
})
