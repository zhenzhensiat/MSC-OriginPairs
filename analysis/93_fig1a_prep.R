# =============================================================================
# 93_fig1a_prep.R - preparation for the screening funnel: UNCLEAR recoding
# (converged categories) and the English category-name table.
# Inputs: data/census_adjudication_v2_r3.tsv,
#   data/fig1a_screening_counts_v2.tsv
# Outputs: data/fig1a_unclear_recode_v2.tsv,
#   data/fig1a_category_labels_v2.tsv
# Discipline: presentation-layer recoding only; the r3 verdicts are untouched.
# =============================================================================
source("scripts/00_setup.R")

run_step("93_fig1a_prep", function() {
  ad <- read_tsv("data/census_adjudication_v2_r3.tsv")
  un <- ad[ad$human_verdict == "UNCLEAR", , drop = FALSE]
  check(nrow(un) == 50L, "r3 UNCLEAR should be 50, actual ", nrow(un))

  # Recode case by case (unresolved_issue is the secondary note; category uses the condensed class name)
  rec <- data.frame(
    accession = c(
      "GSE219227", "GSE20431", "GSE137173", "GSE61153", "GSE267908", "GSE119096",
      "GSE168055", "GSE2776", "GSE35955", "GSE35956", "GSE35957", "GSE35958",
      "GSE48662", "GSE108915", "GSE111822", "GSE118961", "GSE124859", "GSE124860",
      "GSE142563", "GSE152507", "GSE184225", "GSE194303", "GSE222013", "GSE230534",
      "GSE246107", "GSE261358", "GSE45663", "GSE47869", "GSE52841", "GSE52842",
      "GSE52843", "GSE52844", "GSE60525", "GSE87796", "GSE90131", "GSE92640",
      "GSE92778", "GSE102401", "GSE133346", "GSE193021", "GSE231660", "GSE249452",
      "GSE280881", "GSE288316", "GSE293806", "GSE293949", "GSE294580", "GSE309538",
      "GSE309540", "GSE317531"),
    category = c(
      "Insufficient metadata", "Insufficient metadata", "Insufficient metadata",
      "Insufficient metadata", "Insufficient metadata", "Insufficient metadata",
      "Insufficient metadata", "Insufficient metadata", "Insufficient metadata",
      "Insufficient metadata", "Insufficient metadata", "Insufficient metadata",
      "Insufficient metadata", "Insufficient metadata", "Insufficient metadata",
      "Insufficient metadata", "Insufficient metadata", "Insufficient metadata",
      "Insufficient metadata", "Insufficient metadata", "Insufficient metadata",
      "Duplicate or super-series", "Insufficient metadata", "Insufficient metadata",
      "Insufficient metadata", "Insufficient metadata", "Insufficient metadata",
      "Insufficient metadata", "Insufficient metadata", "Insufficient metadata",
      "Insufficient metadata", "Insufficient metadata", "Insufficient metadata",
      "Insufficient metadata", "Insufficient metadata", "Insufficient metadata",
      "Insufficient metadata", "Insufficient metadata", "Insufficient metadata",
      "Insufficient metadata", "Insufficient metadata", "Insufficient metadata",
      "Insufficient metadata", "Insufficient metadata", "Insufficient metadata",
      "Insufficient metadata", "Insufficient metadata", "Insufficient metadata",
      "Insufficient metadata", "Insufficient metadata"),
    unresolved_issue = c(
      "Arm and tissue source not specified",
      "Untreated control arm unverified",
      "Sequencing target (MSC vs tissue) unresolved",
      "Adult-source tissue and donor count unverified",
      "Control (non-diseased) arm and tissue source unverified",
      "Conventional-culture control arm unverified",
      "Normoxia control arm and tissue source unverified",
      "Untreated control arm unverified",
      "Healthy-donor arm unverified; potential team overlap with discovery cohort",
      "Healthy-donor arm unverified; potential team overlap with discovery cohort",
      "Healthy-donor arm unverified; potential team overlap with discovery cohort",
      "Healthy-donor arm unverified; potential team overlap with discovery cohort",
      "Untreated control arm unverified",
      "Healthy-donor arm unverified",
      "Undifferentiated baseline composition unverified",
      "Untreated control arm unverified",
      "Companion-series composition unverified",
      "Companion-series composition unverified",
      "Bone-arm tissue site and donors unverified",
      "Conventional 2D control arm unverified",
      "Conventional 2D arm and donors unverified",
      "Super-series without independent data (represented by sub-series)",
      "Healthy-donor arm unverified",
      "Healthy-donor arm tissue and donors unverified",
      "Standard-culture (FBS) arm inclusion unverified",
      "Healthy-donor arm unverified",
      "Mono-culture control arm unverified",
      "ADSC-side composition unverified",
      "Normal-weight baseline arm unverified",
      "Normal-weight baseline arm unverified",
      "Normal-weight baseline arm unverified",
      "Normal-weight baseline arm unverified",
      "Standard-confluence control arm unverified",
      "Source tissue (putative BM) unverified",
      "Untreated control arm unverified",
      "Standard-protocol arm and donor composition unverified",
      "Healthy-donor arm unverified",
      "Untreated control arm unverified",
      "Healthy-donor arm and tissue site unverified",
      "Baseline sample inclusion unverified",
      "Donor-batch composition unverified",
      "Product transcriptome composition unverified",
      "Healthy-donor arm unverified",
      "Undifferentiated baseline composition unverified",
      "Healthy-donor arm unverified",
      "Donor and tissue site unverified",
      "FBS baseline composition unverified",
      "Healthy-donor arm unverified",
      "Healthy-donor arm unverified",
      "Undifferentiated baseline composition unverified"),
    stringsAsFactors = FALSE)
  check(!anyDuplicated(rec$accession) && nrow(rec) == 50L, "recode table structure anomaly")
  check(setequal(rec$accession, un$accession), "recode set does not match the r3 UNCLEAR set")
  j <- match(rec$accession, un$accession)
  rec$human_reason <- un$human_reason[j]   # verbatim r3 text (to state reasons in Table S1)
  write_tsv(rec, "data/fig1a_unclear_recode_v2.tsv")
  log_msg("UNCLEAR recode: Insufficient metadata ", sum(rec$category == "Insufficient metadata"),
          "; Duplicate or super-series ", sum(rec$category == "Duplicate or super-series"))

  # Plain-English label table for funnel categories (covering all fig1a categories)
  lab <- rbind(
    data.frame(scope = "v2_census_1089", code = c(
      "C_AD", "C_BM", "C_UC", "MULTI_NEW", "ORTHO_SCRNA", "PASSAGE_SERIES", "UNCLEAR",
      "X_CELLLINE_IPSC", "X_DUP", "X_FETAL_EMBRYO", "X_NO_MATRIX", "X_NOT_CULTURE",
      "X_NOT_EXPR", "X_NOT_HUMAN", "X_NOT_MSC", "X_SMALL_N", "X_SOURCE_OOC",
      "X_TREAT_NO_CTRL", "X_EXPLORATORY_PHASE"),
      label_en = c(
        "Multi-source candidate dataset", "Multi-source candidate dataset",
        "Multi-source candidate dataset", "Multi-source candidate dataset",
        "Single-cell dataset (exploratory)", "Passage-series dataset (stability)",
        "Insufficient metadata (recoded)", "Cell line or iPSC-derived",
        "Duplicate or super-series", "Fetal or embryonic tissue",
        "Expression matrix unavailable", "Not culture-expanded",
        "No expression profiling", "Non-human", "Not MSC", "Too few donors",
        "Tissue outside BM/AD/UC", "Treated, no untreated control",
        "Used in prior analyses (v1)"), stringsAsFactors = FALSE),
    data.frame(scope = "v1_inherited", code = c(
      "excluded_cell_lines", "excluded_coculture_profile", "excluded_commercial",
      "excluded_EV", "excluded_insufficient_donors", "excluded_methylation",
      "excluded_miRNA", "excluded_non_expression", "excluded_not_msc_profile",
      "excluded_single_source", "excluded_SVF_uncultured", "L1_confirmed",
      "L1_pending_data", "L1_scRNA", "L2_anchor;L3", "L2_candidate",
      "L2_confirmed_backfilled", "L3_candidate"),
      label_en = c(
        "Cell lines", "Co-culture (profile only)", "Commercial product (uncultured)",
        "EV-treatment study", "Too few donors", "Methylation platform",
        "miRNA platform", "No expression profiling", "Not MSC profile",
        "Single source (v1)", "Uncultured SVF", "Multi-source, v1-confirmed",
        "Multi-source, data pending", "Single-cell (v1)",
        "Anchor series with passage sub-series", "Single-source candidate (v1)",
        "Single-source, v1-confirmed", "Passage-series candidate (v1)"),
      stringsAsFactors = FALSE))
  f1 <- read_tsv("data/fig1a_screening_counts_v2.tsv")
  want <- unique(f1[, c("scope", "category")])
  have <- lab[, c("scope", "code")]; names(have)[2] <- "category"
  check(setequal(paste(want$scope, want$category), paste(have$scope, have$category)),
        "label table and fig1a category sets differ: missing ",
        paste(setdiff(paste(want$scope, want$category), paste(have$scope, have$category)), collapse = ","))
  write_tsv(lab, "data/fig1a_category_labels_v2.tsv")

  log_msg("93 done: 50 recoded rows (49 IM + 1 DUP) + label table ", nrow(lab), " rows")
  invisible(NULL)
})
