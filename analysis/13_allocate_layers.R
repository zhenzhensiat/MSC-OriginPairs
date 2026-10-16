# =============================================================================
# 13_allocate_layers.R —— allocation of the development (D) and validation (C)
# layers (metadata only, no expression values).
# Input: data/census_clusters_v2.tsv, one row per GSE with columns:
#   gse, cluster_id, design(multi_source/single_source), data_type(rnaseq/array),
#   full_genome(TRUE/FALSE), in_v1_L1, in_v1_L2, stability
# Output: data/layer_allocation_v2.tsv; written immediately after running
# Rules:
#   1) one study-group code may enter only one layer;
#   2) clusters containing earlier-phase studies -> D; earlier supplementary -> supp;
#   3) new multi-source RNA-seq clusters: >=6 in number, a fixed-seed random third
#      enters C, the rest D; fewer than 6 all enter D;
#   4) new single-source genome-wide clusters -> C.
# =============================================================================
source("scripts/00_setup.R")

run_step("13_allocate_layers", function() {
  out_rel <- "data/layer_allocation_v2.tsv"
  check(!file.exists(proj_path(out_rel)), out_rel, " already exists: allocation runs only once")
  cc <- read_tsv("data/census_clusters_v2.tsv")
  need <- c("gse", "cluster_id", "design", "data_type", "full_genome", "in_v1_L1", "in_v1_L2", "stability")
  check(all(need %in% names(cc)), "census cluster table is missing columns: ", paste(setdiff(need, names(cc)), collapse = ","))
  check(!anyDuplicated(cc$gse), "duplicate GSEs")
  cc$cluster_id <- as.character(cc$cluster_id)
  for (v in c("full_genome", "in_v1_L1", "in_v1_L2", "stability")) {
    check(is.logical(cc[[v]]) && !anyNA(cc[[v]]), "column ", v, " must be TRUE/FALSE")
  }

  by_cl <- split(cc, cc$cluster_id)
  cl <- data.frame(
    cluster_id   = names(by_cl),
    has_v1_L1    = vapply(by_cl, function(d) any(d$in_v1_L1), logical(1)),
    has_v1_L2    = vapply(by_cl, function(d) any(d$in_v1_L2), logical(1)),
    multi_rnaseq = vapply(by_cl, function(d) any(d$design == "multi_source" & d$data_type == "rnaseq" & d$full_genome), logical(1)),
    multi_any    = vapply(by_cl, function(d) any(d$design == "multi_source"), logical(1)),
    single_full  = vapply(by_cl, function(d) all(d$design == "single_source") && all(d$full_genome), logical(1)),
    stringsAsFactors = FALSE, row.names = NULL)

  cl$layer <- NA_character_
  cl$layer[cl$has_v1_L1] <- "D"
  cl$layer[is.na(cl$layer) & cl$has_v1_L2] <- "supp"

  new_multi <- which(is.na(cl$layer) & cl$multi_rnaseq)
  if (length(new_multi) >= 6) {
    set.seed(derive_seed(13, 1))
    to_c <- new_multi[sample.int(length(new_multi), floor(length(new_multi) / 3))]
    cl$layer[to_c] <- "C"
  }
  cl$layer[is.na(cl$layer) & cl$multi_any] <- "D"
  cl$layer[is.na(cl$layer) & cl$single_full] <- "C"
  cl$layer[is.na(cl$layer)] <- "unassigned"   # e.g. non-genome-wide single source: decided manually and recorded

  alloc <- safe_merge(cc, cl[, c("cluster_id", "layer")], by = "cluster_id", type = "m:1")
  write_tsv(alloc, out_rel)

  tab <- table(alloc$layer)
  log_msg("allocation result (by GSE): ", paste(names(tab), tab, sep = "=", collapse = ", "))
  log_msg("Note: validation-layer single-source studies need >=8 independent clusters per class; check in census_gate.md")
  invisible(alloc)
})
