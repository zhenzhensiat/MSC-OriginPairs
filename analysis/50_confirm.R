# =============================================================================
# 50_confirm.R - external confirmation run (executed once).
# Precondition: the final model files exist (data/model_v2.rds,
#   data/model_v2_tau.rds; details in curation/model_info_v2.json).
# Protection: writes data/confirm_complete on success; afterwards only a
#   deterministic re-run (MSC_CONFIRM_RERUN=1) is allowed, which writes
#   *_rerun.tsv and compares byte-wise with the first results.
# Outputs: confirm_predictions_v2.tsv, confirm_results_v2.tsv,
#   confirm_strata_v2.tsv, confirm_null_v2.tsv
# =============================================================================
source("scripts/00_setup.R")
need_pkgs(c("multiclassPairs", "digest", "jsonlite"))

run_step("50_confirm", function() {
  marker <- proj_path("data", "confirm_complete")
  rerun <- identical(Sys.getenv("MSC_CONFIRM_RERUN"), "1")
  check(!file.exists(marker) || rerun, "Confirmation has already been run once. Deterministic re-run only allowed via MSC_CONFIRM_RERUN=1")
  sfx <- if (rerun) "_rerun" else ""

  assert_multiclasspairs_api()
  fz <- jsonlite::read_json(proj_path("curation", "model_info_v2.json"), simplifyVector = TRUE)

  studies <- load_studies("C")
  overlap <- intersect(vapply(studies, `[[`, "", "cluster_id"), fz$training_clusters)
  check(length(overlap) == 0, "confirmation studies overlap training clusters: ", paste(overlap, collapse = ","))
  meta <- bind_meta(studies)

  clf <- readRDS(proj_path("data", "model_v2.rds"))
  tau <- readRDS(proj_path("data", "model_v2_tau.rds"))
  pr <- do.call(rbind, lapply(studies, function(s) {
    m <- s$meta[s$meta$arm == "include" & s$meta$source %in% CFG$classes, , drop = FALSE]
    ktsp_predict(clf, s$expr[, m$gsm, drop = FALSE])
  }))
  pr <- apply_reject(pr, tau)
  df <- attach_truth(pr, meta)
  sc <- data.frame(study = names(studies), scale = vapply(studies, `[[`, "", "scale"), stringsAsFactors = FALSE)
  df <- safe_merge(df, sc, by = "study", type = "m:1")
  df <- safe_merge(df, meta[, c("gsm", "passage", "medium")], by = "gsm", type = "m:1")
  write_tsv(df, paste0("data/confirm_predictions_v2", sfx, ".tsv"))

  # ---- Primary results and verdict ----------------------------------------------------------------
  s_main <- summarize_predictions(df)
  boot <- cluster_bootstrap(df, function(d) balanced_accuracy(d$true, d$pred), seed = derive_seed(50, 1))
  verdict <- h2_verdict(s_main, boot)
  d_os <- df; d_os$pred <- d_os$pred_openset
  s_os <- summarize_predictions(d_os)

  # ---- Chance level: permute 'cluster -> class' labels among single-source clusters (predictions unchanged)----------------------
  cl_lab <- tapply(df$true, df$cluster_id, function(v) if (length(unique(v)) == 1) v[1] else NA_character_)
  single <- names(cl_lab)[!is.na(cl_lab)]
  set.seed(derive_seed(50, 2))
  perm <- vapply(seq_len(CFG$boot_R), function(b) {
    lab <- cl_lab; lab[single] <- sample(cl_lab[single])
    d <- df; idx <- d$cluster_id %in% single
    d$true[idx] <- lab[d$cluster_id[idx]]
    balanced_accuracy(d$true, d$pred)
  }, numeric(1))
  write_tsv(data.frame(perm_ba = perm), paste0("data/confirm_null_v2", sfx, ".tsv"))

  pc <- s_main$per_class
  res <- data.frame(
    metric = c("balanced_accuracy", "ba_ci_lo", "ba_ci_hi", paste0("recall_", pc$class),
               "study_acc_median", "study_acc_min", "majority_acc", "donor_ba", "abstain_rate",
               "openset_ba", "rejection_rate", "perm_null_mean", "perm_null_q95",
               paste0("clusters_", names(boot$clusters_per_class)), "ci_unreliable", "H2_verdict"),
    value = c(s_main$ba, boot$lo, boot$hi, pc$recall,
              s_main$study_acc_median, s_main$study_acc_min, s_main$majority_acc, s_main$donor_ba, s_main$abstain_rate,
              s_os$ba, mean(df$rejected), mean(perm), unname(stats::quantile(perm, 0.95)),
              boot$clusters_per_class, boot$unreliable, verdict),
    stringsAsFactors = FALSE)
  write_tsv(res, paste0("data/confirm_results_v2", sfx, ".tsv"))

  # ---- Stratification (fixed)--------------------------------------------------------------
  strata <- list(scale = df$scale,
                 passage_bin = ifelse(is.na(df$passage), "unknown", ifelse(df$passage <= 5, "<=P5", ">P5")),
                 medium = ifelse(is.na(df$medium) | !nzchar(df$medium), "unknown", df$medium))
  st <- do.call(rbind, lapply(names(strata), function(v) {
    do.call(rbind, lapply(sort(unique(strata[[v]])), function(lv) {
      d <- df[strata[[v]] == lv, , drop = FALSE]
      data.frame(stratum = v, level = lv, n = nrow(d), n_studies = length(unique(d$study)),
                 ba = balanced_accuracy(d$true, d$pred), stringsAsFactors = FALSE)
    }))
  }))
  write_tsv(st, paste0("data/confirm_strata_v2", sfx, ".tsv"))

  if (rerun) {
    for (f in c("confirm_predictions_v2", "confirm_results_v2", "confirm_strata_v2", "confirm_null_v2")) {
      same <- identical(readLines(proj_path("data", paste0(f, ".tsv")), warn = FALSE),
                        readLines(proj_path("data", paste0(f, "_rerun.tsv")), warn = FALSE))
      log_msg("deterministic re-run ", f, ": ", if (same) "identical output" else "MISMATCH; investigate")
    }
  } else {
    writeLines(c("confirmation executed", format(Sys.time())), marker)
  }
  log_msg("H2 verdict: ", verdict, " (balanced accuracy = ", round(s_main$ba, 3),
          ", 95% CI ", round(boot$lo, 3), "–", round(boot$hi, 3), ")")
  invisible(res)
})
