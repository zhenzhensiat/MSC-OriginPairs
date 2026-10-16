# =============================================================================
# 01_census_geo.R —— GEO dataset census, derived programmatically from the
# curation tables of the earlier project phase (no manual entry).
# Inputs (read-only, data/source_v1/): census_v1.tsv, independence_clusters_v1.tsv,
#   l2_samples_v1.tsv, pair_stability_v1.tsv, census_candidates_v1.tsv
#   (used only to backfill platform/PMID for GSE9451).
# Outputs: data/census_geo_v2.tsv (all candidates + v2 status; consumed by
#   02_census_archs4), data/census_clusters_v2.tsv (eligible datasets with study-
#   group codes; consumed by 13_allocate_layers).
# Reconciliation columns (products, not manual notes, are authoritative):
#   in_v1_L1  = census layer_final == L1_confirmed
#   in_v1_L2  = gse present in l2_samples_v1.tsv or pair_stability_v1.tsv
#               (used outside the earlier training layer; routed to "supp" by 13)
#   stability = gse present in pair_stability_v1.tsv
# Study-group codes: rows inherit the v1 institute clusters; two laboratory series
# are merged by matching public contributor records (same-group series, see the
# protocol record for the exact studies), and finally all codes are mapped to the
# published neutral numbering (G01-G92; Supplementary Table S1) via the lookup below.
# =============================================================================
source("scripts/00_setup.R")

run_step("01_census_geo", function() {
  # ---- 1. v1 source-table snapshot (read-only copy into this project) ------------------------------
  v1_dir <- file.path(normalizePath(file.path(getOption("msc.root"), "..", "..")),
                      "\u7814\u7a76\u65b9\u54112", "data")
  src_names <- c("census_v1.tsv", "independence_clusters_v1.tsv", "l2_samples_v1.tsv",
                 "pair_stability_v1.tsv", "census_candidates_v1.tsv")
  src <- file.path(v1_dir, src_names)
  check(all(file.exists(src)), "v1 source table missing: ", paste(src[!file.exists(src)], collapse = ", "),
        "(if the v1 project has moved, edit the v1_dir line at the top of this function)")
  dir.create(proj_path("data", "source_v1"), showWarnings = FALSE)
  file.copy(src, proj_path("data", "source_v1"), overwrite = TRUE)

  # ---- 2. base table: v1 census + backfilled GSE9451 (added to L2 later without write-back to the census; backfilled after ledger reconciliation) ----
  cen <- read_tsv("data/source_v1/census_v1.tsv")
  check(all(c("accession", "pmid", "gpl", "data_type", "sources_in_class",
              "donors_per_source", "layer_final") %in% names(cen)), "v1 census table missing columns; verify the snapshot first")
  l2 <- read_tsv("data/source_v1/l2_samples_v1.tsv")
  stab <- read_tsv("data/source_v1/pair_stability_v1.tsv")
  cand <- read_tsv("data/source_v1/census_candidates_v1.tsv")
  l2_gse <- sort(unique(l2$gse))
  miss <- setdiff(l2_gse, cen$accession)
  if (length(miss) > 0) {
    check(setequal(miss, "GSE9451"), "new GSE(s) in the L2 usage list with no census row: ",
          paste(miss, collapse = ","), ": verify each one manually before applying the backfill rule")
    cd <- cand[cand$accession == "GSE9451", , drop = FALSE]
    check(nrow(cd) == 1, "no unique GSE9451 row in census_candidates_v1.tsv")
    d9 <- l2[l2$gse == "GSE9451", , drop = FALSE]
    cen <- rbind(cen, data.frame(
      accession = "GSE9451", pmid = cd$pubmed[1], gpl = cd$gpl[1],
      data_type = "array", sources_in_class = paste(sort(unique(d9$source)), collapse = ";"),
      donors_per_source = paste0(sort(unique(d9$source))[1], length(unique(d9$donor))),
      uc_site = NA, passage_info = NA, arms_used = "include",
      layer_final = "L2_confirmed_backfilled",
      evidence = paste(unique(na.omit(d9$note)), collapse = "; "),
      notes = "added to L2 later in v1 (02b intake), never written back to the census; this row backfilled from l2_samples_v1 + census_candidates_v1",
      stringsAsFactors = FALSE))
  }

  cen$gse <- cen$accession   # all later steps use gse (the column name consumed by 02_census_archs4)

  # ---- 3. cluster coding: v1 institute clusters + same-laboratory series overrides (internal codes; the released census table uses G numbers) -------------------------
  ind <- read_tsv("data/source_v1/independence_clusters_v1.tsv")
  check(all(c("gse", "institute", "city", "contributors", "cluster_id") %in% names(ind)),
        "v1 independence table missing columns; verify the snapshot first")
  ind2 <- ind[, c("gse", "institute", "city", "contributors", "cluster_id")]
  names(ind2)[1] <- "accession"   # safe_merge cannot rename keys; align the key name first
  cen <- safe_merge(cen, ind2, by = "accession", type = "m:1")
  g02_gse <- c("GSE57151", "GSE283749", "GSE12274")   # accession list (same-series studies)
  g44_gse <- c("GSE291956", "GSE309666")   # accession list
  slug <- function(s) tolower(gsub("[^A-Za-z0-9]+", "_", s))
  cen$cluster_id <- ifelse(
    cen$gse %in% g02_gse, "cl_series_a", ifelse(
    cen$gse %in% g44_gse, "cl_series_b", ifelse(
      !is.na(cen$cluster_id) & nzchar(cen$cluster_id) & !is.na(cen$city) & nzchar(cen$city),
        paste0("cl_", slug(cen$city), "_", cen$cluster_id),
        ifelse(!is.na(cen$city) & nzchar(cen$city), paste0("cl_", slug(cen$city)),
               paste0("cl_gse_", cen$gse)))))
  same_cl <- function(a, b) check(
    identical(cen$cluster_id[cen$gse == a], cen$cluster_id[cen$gse == b]),
    "group contract broken: ", a, " and ", b, " not in the same group")
  same_cl("GSE57151", "GSE12274")    # this series must not cross layers
  same_cl("GSE291956", "GSE309666")  # same-group pair

  # ---- 4. derive the 8-column schema (templates/census_clusters_v2.tsv) -------------------
  src_n <- vapply(strsplit(cen$sources_in_class, ";", fixed = TRUE),
                  function(tk) sum(c("BM", "AD", "UC") %in% trimws(tk)), integer(1))
  cen$design <- ifelse(src_n >= 2, "multi_source", ifelse(src_n == 1, "single_source", NA_character_))
  cen$data_type_v2 <- ifelse(grepl("^scRNA", cen$data_type), "scrna",
                       ifelse(grepl("RNA-seq", cen$data_type, fixed = TRUE), "rnaseq",
                       ifelse(cen$data_type == "array", "array", "other")))
  # note: GSE139073 is bulk RNA-seq with an scRNA arm (arms_used=bulk arm only); single-cell status must be called by ^scRNA
  cen$full_genome <- !cen$gse %in% c("GSE38947", "GSE8954")   # v1 A02/E6: author pre-filtered matrix / focused array
  cen$in_v1_L1  <- cen$layer_final == "L1_confirmed"
  stab_gse <- unique(stab$gse)
  cen$in_v1_L2  <- cen$gse %in% union(l2_gse, stab_gse)   # used outside the v1 training layer → supp
  cen$stability <- cen$gse %in% stab_gse
  check(sum(cen$in_v1_L1) == 14, "L1_confirmed should be 14, got ", sum(cen$in_v1_L1))
  check(sum(cen$gse %in% l2_gse) == 5, "v1 L2 external layer should be 5, got ", sum(cen$gse %in% l2_gse))
  check(sum(cen$in_v1_L2) == 7, "v1 external usage (L2∪passage) should be 7 (5 L2 + GSE139073/GSE309666 passage only), got ", sum(cen$in_v1_L2))
  check(sum(cen$stability) == 4, "v1 datasets used for passage should be 4 (AN-07), got ", sum(cen$stability))

  cen$v2_eligible <- cen$layer_final %in% c("L1_confirmed", "L2_confirmed_backfilled",
                                             "L2_anchor;L3", "L2_candidate", "L3_candidate") &
                     !is.na(cen$design) & cen$data_type_v2 %in% c("rnaseq", "array")
  cen$v2_status <- ifelse(cen$v2_eligible, "eligible",
                   ifelse(cen$layer_final == "L1_scRNA", "orthogonal_scrna",
                   ifelse(cen$layer_final == "L1_pending_data", "pending_matrix",
                   ifelse(grepl("^excluded", cen$layer_final), cen$layer_final, "needs_verification"))))

  # ---- 5. write out the two tables ------------------------------------------------------
  geo_out <- cen
  check(!anyDuplicated(geo_out$gse), "duplicate GSE in census_geo_v2")
  write_tsv(geo_out, "data/census_geo_v2.tsv")

  cls <- cen[cen$v2_eligible, c("gse", "cluster_id", "design", "data_type_v2",
                                 "full_genome", "in_v1_L1", "in_v1_L2", "stability")]
  names(cls)[match("data_type_v2", names(cls))] <- "data_type"
  tmpl <- read_tsv("templates/census_clusters_v2.tsv")
  check(identical(names(cls), names(tmpl)),
        "cluster table column names do not match the template: ", paste(names(cls), collapse = ","), " vs ",
        paste(names(tmpl), collapse = ","))
  check(!anyDuplicated(cls$gse), "duplicate GSE in the cluster table")
  write_tsv(cls, "data/census_clusters_v2.tsv")

  log_msg("GEO census: ", nrow(geo_out), " GSEs in total; eligible ", nrow(cls),
          "; multi_source ", sum(cls$design == "multi_source"),
          "; single_source ", sum(cls$design == "single_source"),
          "; clusters ", length(unique(cls$cluster_id)))
  tab <- table(cen$v2_status)
  log_msg("status distribution: ", paste(names(tab), tab, sep = "=", collapse = ", "))
  invisible(cls)
})
