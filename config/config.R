# =============================================================================
# Global settings: thresholds, seeds and class definitions used by the pipeline.
# =============================================================================

CFG <- list(
  master_seed   = 20261009L,
  classes       = c("AD", "BM", "UC"),            # fixed alphabetical order; also the class order for ties
  contrasts     = c("BM-AD", "BM-UC", "AD-UC"),   # direction convention: former minus latter

  # ---- resources ----
  mem_limit_gb  = 14,                              # watchdog threshold (hard cap 16 GB)

  # ---- data build ----
  genome_min_genes = 10000,   # below this gene count = not genome-wide: profiling only, never training
  annot_hit_min    = 0.80,    # warn when probe annotation hit rate falls below this
  log_scale_max    = 30,      # matrix max above this counts as not log-transformed: error

  # ---- marker transfer audit (Phase 2) ----
  pc_q        = 0.05,         # partial-conjunction BH q threshold
  min_K       = 2L,           # at least 2 evaluable studies required
  n_null_sets = 1000L,        # number of matched-null sets
  n_family_null_sets = 200L,  # null sets per family in the gene-family ranking
  family_size = c(15L, 80L),  # family-size range entering the ranking
  family_testable_min = 0.80, # lower bound on the testable gene x contrast fraction per family
  expr_bins   = 10L,          # matching: mean-expression deciles
  detect_bins = 5L,           # matching: detection-rate quintiles
  min_pool    = 5L,           # relax stepwise when the matching pool is smaller than this

  # ---- classifier (Phase 4) ----
  anchor_size = 39L,
  n_random_anchor_sets = 100L,
  k_range     = c(3L, 5L, 8L, 10L, 15L),  # checked against the multiclassPairs k_range semantics
  genome_featureNo = 500L,                # per-class gene count in the primary model
  tie_tolerance = 0.02,                   # metric difference below 0.02 counts as a tie
  loso_exclude_clusters = c("cl_kuala_lumpur_18"),  # preselected-matrix series; excluded from LOSO evaluation

  # ---- decision thresholds ----
  h2 = list(ba = 0.85, ci_lo = 0.70, recall = 0.75),
  stability = list(min_same_sign = 0.90, min_donors = 3L),
  boot_R = 1000L,
  min_clusters_reliable = 5L
)
