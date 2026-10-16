# =============================================================================
# 44_final_model.R —— train the final model on the full development layer.
# Precondition: the final method has been confirmed; fill FINAL below accordingly
# and record the choice alongside the model files.
# Outputs: data/model_v2.rds, data/model_v2_rules.tsv, data/model_v2_tau.rds,
#   data/model_v2_glmnet.rds (secondary comparator).
# =============================================================================
source("scripts/00_setup.R")
need_pkgs(c("multiclassPairs", "digest", "jsonlite"))

FINAL <- list(
  method       = "genome_main",   # matches the method names in 40: genome_main / HOX / published / genome39
  anchors_file = NA_character_,   # for the anchor method, provide the anchor gene table (relative path, column 'gene') and add it first
  k_range      = CFG$k_range
)

run_step("44_final_model", function() {
  check(!file.exists(proj_path("data", "model_v2.rds")), "model already written; do not overwrite; save changes as a new version file")
  assert_multiclasspairs_api()
  alloc <- read_tsv("data/layer_allocation_v2.tsv")

  studies <- load_studies("D")
  lay <- alloc$layer[match(names(studies), alloc$gse)]
  check(!anyNA(lay) && all(lay == "D"), "training studies not in layer D: ", paste(names(studies)[is.na(lay) | lay != "D"], collapse = ","))

  # final model trained on the analysis universe: D intersection intersected with every C platform (6,643 genes; rules defined in 44a)
  universe <- read_tsv("data/universe_candidate_v2.tsv")$gene
  tv <- training_view(studies, exclude_cluster = "__none__", restrict_genes = universe)
  check(length(tv$genes) >= 6000, "analysis-universe training set has fewer than 6000 genes: ", length(tv$genes))
  log_msg("analysis-universe training: ", length(tv$genes), " genes / ", ncol(tv$expr), " samples / ",
          length(unique(tv$meta$cluster_id)), " clusters")
  anchors <- NULL
  if (!is.na(FINAL$anchors_file)) {
    anchors <- read_tsv(FINAL$anchors_file)$gene
  }
  clf <- ktsp_train(tv, anchors = anchors, include_pivot = !is.null(anchors),
                    seed = derive_seed(44, 1), k_range = FINAL$k_range)
  saveRDS(clf, proj_path("data", "model_v2.rds"))

  rules <- ktsp_rules(clf, tv$expr, tv$meta$source)
  write_tsv(rules, "data/model_v2_rules.tsv")
  log_msg("rule count: ", nrow(rules), "; genes involved ", length(unique(c(rules$gene1, rules$gene2))), "")

  tau_src <- proj_path("data", "loso_v2r2", paste0("tau_", FINAL$method, ".rds"))
  check(file.exists(tau_src), "missing the LOSO rejection threshold for this method: ", tau_src)
  saveRDS(readRDS(tau_src), proj_path("data", "model_v2_tau.rds"))

  # self-check: accuracy on re-substituted training data is recorded only (not a performance estimate)
  resub <- ktsp_predict(clf, tv$expr)
  log_msg("re-substitution balanced accuracy (self-check only): ",
          round(balanced_accuracy(tv$meta$source, resub$pred), 3))

  # selected k (per class; from one-vs-rest classifier names, e.g. '5TSPS' means 5 rules for that class)
  k_sel <- vapply(clf$classifiers, function(cc) as.integer(sub("TSPS$", "", cc$name)), integer(1))
  nm_k <- names(clf$classifiers)
  if (is.null(nm_k)) nm_k <- paste0("classifier", seq_along(k_sel))
  names(k_sel) <- nm_k
  log_msg("selected k per class: ", paste(names(k_sel), k_sel, sep = "=", collapse = ", "))

  # secondary comparator: glmnet (within-sample rank transform + multinomial lasso, lambda.1se); no rejection threshold, descriptive reporting only
  gl <- baseline_glmnet(tv, tv$expr, seed = derive_seed(44, 2), return_fit = TRUE)
  saveRDS(list(fit = gl$fit, genes = rownames(tv$expr), alpha = 0.5, seed = derive_seed(44, 2),
               n_genes = gl$n_genes,
               note = "glmnet secondary comparator; rank-transformed expression; no reject threshold"),
          proj_path("data", "model_v2_glmnet.rds"))
  log_msg("glmnet comparator re-substitution BA (self-check only): ",
          round(balanced_accuracy(tv$meta$source, gl$pred$pred), 3), "; genes with non-zero coefficients ", gl$n_genes)

  info <- list(
    method = FINAL$method, k_range = FINAL$k_range, selected_k = as.list(k_sel),
    tau = readRDS(proj_path("data", "model_v2_tau.rds")),
    training_studies = names(tv$studies),
    training_clusters = sort(unique(tv$meta$cluster_id)),
    n_training_samples = ncol(tv$expr), n_universe_genes = length(tv$genes),
    multiclassPairs_version = as.character(utils::packageVersion("multiclassPairs")),
    R_version = R.version.string, date = format(Sys.Date()))
  jsonlite::write_json(info, proj_path("curation", "model_info_v2.json"), auto_unbox = TRUE, pretty = TRUE)
  log_msg("model and secondary comparator written (details: curation/model_info_v2.json)")
  invisible(info)
})
