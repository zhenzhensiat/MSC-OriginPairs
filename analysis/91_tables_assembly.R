# =============================================================================
# 91_tables_assembly.R - assembly of supplementary table inputs.
# Inputs (read-only): data/census_adjudication_v2_r3.tsv (authoritative
#   revision), census_geo_v2.tsv, confirm_results_v2.tsv,
#   confirm_sensitivity_v2.tsv; layer sample counts from the study objects and
#   the robustness manifests under curation/.
# Outputs: data/fig1a_screening_counts_v2.tsv (two-part funnel counts),
#   data/table_s1_candidates_v2.tsv (all 1,089 candidates + group/layer),
#   data/platform_composition_v2.tsv (layer x source x platform x samples),
#   data/confirm_sensitivity_summary_v2.tsv (six forest-plot entries).
# The medium-strata table already exists (confirm_strata_curated_v2.tsv).
# =============================================================================
source("scripts/00_setup.R")

run_step("91_tables_assembly", function() {
  a <- read_tsv("data/census_adjudication_v2_r3.tsv")
  g <- read_tsv("data/census_geo_v2.tsv")
  check(!anyDuplicated(a$accession), "the 1089 candidate accessions should be unique")

  # ---- 1. Fig 1a pipeline counts (two-part: verdicts for the v2 new pool + layer assignment inherited from v1)----
  c1 <- data.frame(scope = "v2_census_1089", category = names(table(a$human_verdict)),
                   count = as.integer(table(a$human_verdict)), stringsAsFactors = FALSE)
  c2 <- data.frame(scope = "v1_inherited", category = names(table(g$layer_final)),
                   count = as.integer(table(g$layer_final)), stringsAsFactors = FALSE)
  s1 <- rbind(c1, c2)
  s1$total_check <- ifelse(s1$scope == "v2_census_1089", nrow(a), nrow(g))
  check(sum(c1$count) == 1089L && sum(c2$count) == 49L, "two-pool count totals anomaly")
  write_tsv(s1, "data/fig1a_screening_counts_v2.tsv")

  # ---- 2. Table S1: master table of 1089 adjudications + cluster/layer assignment (join with the v1 legacy table; left keys may be missing, hence m:1)----
  t1 <- safe_merge(a, g[, c("accession", "cluster_id", "layer_final", "stability")],
                   by = "accession", type = "m:1", all.x = TRUE)
  mv <- !is.na(t1$layer_final)
  check(sum(mv) == nrow(g), "hit count should equal the v1 legacy table row count (", nrow(g), "), actual ", sum(mv))
  check(setequal(t1$accession[mv], g$accession), "hit set should be exactly the 49 v1 legacy studies")
  check(all(t1$human_verdict[mv] == "X_EXPLORATORY_PHASE"),
        "hit rows should have human_verdict == X_EXPLORATORY_PHASE (semantic check; machine_hint is a machine pre-screen column, agreeing with manual review in 35/49, not used as the check basis)")
  write_tsv(t1, "data/table_s1_candidates_v2.tsv")

  # ---- 3. Platform composition (study-level; layer x source x platform x sample count; Fig 1d stacked data source)----
  # Sources split into columns (n_AD/n_BM/n_UC/n_other) for stacking; the D layer carries loso_eval/train_used flags
  #  (evaluation 18/174 and training 17/166 definitions); the supp layer includes all 12 studies (10 passage series with passage_series).
  comp <- list()
  scale_plot <- c(array = "microarray", logCPM = "RNA-seq counts", log2TPM = "RNA-seq TPM/FPKM")
  add <- function(s, layer) {
    m <- s$meta[s$meta$arm == "include" & s$meta$source %in% CFG$classes, , drop = FALSE]
    src <- m$source
    nt <- ""
    if (identical(s$gse, "GSE38947")) nt <- "pre-screen matrix: exited LOSO evaluation under D2; not genome-wide (not used in training)"
    if (identical(s$gse, "GSE8954"))  nt <- "focused array; not genome-wide (not used in training; LOSO fold universe 4563 genes)"
    if (identical(layer, "supp")) nt <- if (isTRUE(s$stability)) "passage series" else "non-passage"
    data.frame(layer = layer, gse = s$gse, cluster_id = s$cluster_id,
               n_samples = nrow(m), n_AD = sum(src == "AD"), n_BM = sum(src == "BM"),
               n_UC = sum(src == "UC"), n_other = 0L,
               platform_class = s$scale, platform_plot = unname(scale_plot[s$scale]),
               gpl = g$gpl[match(s$gse, g$gse)], full_genome = isTRUE(s$full_genome),
               loso_eval = if (identical(layer, "D")) !(s$cluster_id %in% CFG$loso_exclude_clusters) else NA,
               train_used = if (identical(layer, "D")) isTRUE(s$full_genome) else NA,
               passage_series = isTRUE(s$stability), note = nt, stringsAsFactors = FALSE)
  }
  d_st <- load_studies("D"); c_st <- load_studies("C"); s_st <- load_studies("supp")
  for (s in d_st) comp[[length(comp) + 1L]] <- add(s, "D")
  for (s in c_st) comp[[length(comp) + 1L]] <- add(s, "C")
  for (s in s_st) comp[[length(comp) + 1L]] <- add(s, "supp")
  pc0 <- do.call(rbind, comp); rownames(pc0) <- NULL
  sum_l <- function(d, l) sum(d$n_samples[d$layer == l])
  check(length(d_st) == 19L && sum_l(pc0, "D") == 180L, "D layer should have 19 studies / 180 samples")
  check(sum(pc0$n_samples[pc0$layer == "D" & pc0$loso_eval %in% TRUE]) == 174L, "D evaluated samples should be 174")
  check(sum(pc0$n_samples[pc0$layer == "D" & pc0$train_used %in% TRUE]) == 166L, "D training samples should be 166")
  check(length(c_st) == 22L && sum_l(pc0, "C") == 176L, "C layer should have 22 studies / 176 samples")
  check(length(s_st) == 12L && sum_l(pc0, "supp") == 214L, "supp layer should have 12 studies / 214 samples")
  check(sum(pc0$passage_series[pc0$layer == "supp"] %in% TRUE) == 10L &&
        sum(pc0$n_samples[pc0$layer == "supp" & pc0$passage_series %in% TRUE]) == 200L,
        "supp-layer passage series should be 10 studies / 200 samples")
  arm <- read_tsv("curation/c_arm_manifest_v2.tsv")
  ooc <- read_tsv("curation/ooc_manifest_v2.tsv")
  arm_i <- arm[arm$include_in_53 %in% c(TRUE, "TRUE"), , drop = FALSE]
  check(nrow(arm_i) == 140L, "samples included in 53 arms should be 140, actual ", nrow(arm_i))
  check(nrow(ooc) == 65L && all(ooc$include %in% c(TRUE, "TRUE")), "OOC should be 65 samples, all included")
  home <- setNames(pc0$layer, pc0$gse); scale_by <- setNames(pc0$platform_class, pc0$gse)
  pca <- do.call(rbind, lapply(split(arm_i, arm_i$gse), function(x) {
    hl <- if (x$gse[1] %in% names(home)) home[[x$gse[1]]] else "not built (arm-only)"
    sc <- unname(scale_by[x$gse[1]])
    data.frame(layer = "O_arm", gse = x$gse[1], cluster_id = x$cluster_id[1],
               n_samples = nrow(x), n_AD = sum(x$source == "AD"), n_BM = sum(x$source == "BM"),
               n_UC = sum(x$source == "UC"), n_other = 0L,
               platform_class = sc, platform_plot = unname(scale_plot[sc]),
               gpl = g$gpl[match(x$gse[1], g$gse)], full_genome = NA,
               loso_eval = NA, train_used = NA, passage_series = NA,
               note = paste0("arm samples; home layer=", hl), stringsAsFactors = FALSE)
  }))
  pco <- do.call(rbind, lapply(split(ooc, ooc$gse), function(x)
    data.frame(layer = "O_ooc", gse = x$gse[1], cluster_id = x$cluster_id[1],
               n_samples = nrow(x), n_AD = 0L, n_BM = 0L, n_UC = 0L, n_other = nrow(x),
               platform_class = NA_character_, platform_plot = NA_character_,
               gpl = g$gpl[match(x$gse[1], g$gse)], full_genome = NA,
               loso_eval = NA, train_used = NA, passage_series = NA,
               note = paste0("Non-target tissue (", paste(sort(unique(x$tissue)), collapse = "/"), ")"),
               stringsAsFactors = FALSE)))
  pc <- rbind(pc0, pca, pco); rownames(pc) <- NULL
  check(sum_l(pc, "O_arm") == 140L && sum_l(pc, "O_ooc") == 65L, "O-layer sample count anomaly")
  check(!anyNA(pc$n_samples) && all(pc$n_samples == pc$n_AD + pc$n_BM + pc$n_UC + pc$n_other),
        "source split counts do not sum to the total")
  write_tsv(pc, "data/platform_composition_v2.tsv")

  # ---- 4. Fig 5d sensitivity summary (six items, ordered as in the figure; display generated from fixed values)----
  cs <- read_tsv("data/confirm_sensitivity_v2.tsv")
  cr <- read_tsv("data/confirm_results_v2.tsv")
  pick <- function(analysis, metric) {
    v <- cs$value[cs$analysis == analysis & cs$metric == metric]
    check(length(v) == 1L, "sensitivity item missing: ", analysis, "/", metric)
    v
  }
  crm <- function(metric) {
    v <- cr$value[cr$metric == metric]
    check(length(v) == 1L, "confirmation result item missing: ", metric)
    v
  }
  f3 <- function(x) sprintf("%.3f", as.numeric(x))
  sm <- data.frame(
    item = c("main", "excl_GSE24871", "leave_one_cluster_min", "donor_level",
             "evaluable_samples_169_169", "evaluable_clusters_20_20"),
    value = c(pick("main_recheck", "balanced_accuracy"),
              pick("excl_GSE24871", "balanced_accuracy"),
              pick("leave_one_cluster_out", "ba_min"),
              pick("donor_level", "balanced_accuracy"),
              pick("evaluable_samples", "exact_ci_lo"),
              pick("evaluable_clusters", "exact_ci_lo")),
    display = c(sprintf("BA %s (95%%CI %s–%s; main analysis)",
                        f3(pick("main_recheck", "balanced_accuracy")),
                        f3(crm("ba_ci_lo")), f3(crm("ba_ci_hi"))),
                sprintf("BA = %s after excluding GSE24871", f3(pick("excl_GSE24871", "balanced_accuracy"))),
                sprintf("lowest BA under leave-one-cluster-out = %s", f3(pick("leave_one_cluster_out", "ba_min"))),
                sprintf("donor-level BA = %s", f3(pick("donor_level", "balanced_accuracy"))),
                sprintf("%s evaluable samples correct (exact binomial 95%% lower bound %s)",
                        pick("evaluable_samples", "correct/n"),
                        f3(pick("evaluable_samples", "exact_ci_lo"))),
                sprintf("%s evaluable clusters correct (exact binomial 95%% lower bound %s)",
                        pick("evaluable_clusters", "all_correct/n"),
                        f3(pick("evaluable_clusters", "exact_ci_lo")))),
    source_table = "data/confirm_sensitivity_v2.tsv / confirm_results_v2.tsv",
    post_hoc = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE), stringsAsFactors = FALSE)
  write_tsv(sm, "data/confirm_sensitivity_summary_v2.tsv")

  log_msg("91: four tables produced; census total ", sum(c1$count), " (expected 1089), legacy pool ", sum(c2$count), " (expected 49); ",
          "platform composition ", nrow(pc), " rows (D ", length(d_st), "/C ", length(c_st), "/supp ", length(s_st),
          "/O_arm ", length(unique(arm_i$gse)), "/O_ooc ", length(unique(ooc$gse)), " studies)")
  invisible(NULL)
})
