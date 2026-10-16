# =============================================================================
# 53_arm_robustness.R - source calls on excluded validation-cohort arms
# (diseased donors / treated / differentiated / non-standard culture;
# exploratory). Questions whether the final rules still identify source in
# non-standard states, directly relevant to quality-control scenarios.
# Precondition: confirmation completed; excluded samples are used, the primary
#   confirmation is unaffected.
# Input: curation/c_arm_manifest_v2.tsv (curated before running),
#   one row per GSM: gse, cluster_id, gsm, source, donor, arm_type, arm_detail,
#   passage, medium, include_in_53; arm_type vocabulary: disease_donor /
#   treated / primed / differentiated / culture_variant / other
# Outputs: data/arm_robustness_pred_v2.tsv, data/arm_robustness_summary_v2.tsv
# =============================================================================
source("scripts/00_setup.R")
source(proj_path("R", "utils_build.R"), encoding = "UTF-8")
need_pkgs(c("multiclassPairs", "switchBox"))

run_step("53_arm_robustness", function() {
  man <- read_tsv("curation/c_arm_manifest_v2.tsv")
  man <- man[man$include_in_53 %in% c(TRUE, "TRUE") & man$source %in% CFG$classes, , drop = FALSE]
  check(nrow(man) > 0, "no usable samples in the manifest")
  ok_types <- c("disease_donor", "treated", "primed", "differentiated", "culture_variant", "other")
  check(all(man$arm_type %in% ok_types), "arm_type contains values outside the controlled vocabulary")
  man$arm <- "include"   # only to satisfy the validate_study column convention; the analysis uses arm_type

  clf <- readRDS(proj_path("data", "model_v2.rds"))
  tau <- readRDS(proj_path("data", "model_v2_tau.rds"))
  preds <- lapply(split(man, man$gse), function(mr) {
    s <- tryCatch(build_from_manifest(mr$gse[1], mr, layer = "O", cluster_id = mr$cluster_id[1]),
                  error = function(e) conditionMessage(e))
    if (is.character(s)) {   # write build failures to the table as well
      return(data.frame(gsm = mr$gsm, pred = NA_character_, margin = NA_real_, rejected = NA,
                        build_status = s, stringsAsFactors = FALSE))
    }
    p <- apply_reject(ktsp_predict(clf, s$expr), tau)
    p$build_status <- "ok"
    p
  })
  pr <- do.call(rbind, lapply(preds, function(p) p[, intersect(c("gsm", "pred", "margin", "rejected", "build_status"), names(p))]))
  pr <- safe_merge(pr, man[, c("gsm", "gse", "cluster_id", "source", "donor", "arm_type", "arm_detail")], by = "gsm", type = "1:1")
  names(pr)[names(pr) == "source"] <- "true"
  names(pr)[names(pr) == "gse"] <- "study"
  write_tsv(pr, "data/arm_robustness_pred_v2.tsv")

  sm <- do.call(rbind, lapply(split(pr, pr$arm_type), function(d) {
    data.frame(arm_type = d$arm_type[1], n = nrow(d), n_studies = length(unique(d$study)),
               accuracy = mean(!is.na(d$pred) & d$pred == d$true),
               abstain = mean(is.na(d$pred)),
               rejected = mean(d$rejected %in% c(TRUE, "TRUE")),
               median_margin = stats::median(d$margin, na.rm = TRUE),
               post_hoc = TRUE, stringsAsFactors = FALSE)
  }))
  write_tsv(sm, "data/arm_robustness_summary_v2.tsv")
  invisible(sm)
})
