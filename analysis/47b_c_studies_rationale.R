# =============================================================================
# 47b_c_studies_rationale.R - per-study written rationale table for the
# validation cohort (needed for Fig 1 and the Methods).
# Purpose: a written justification per study (retain / exclude / build-blocked),
#   including donor-count definitions and known process notes.
# Inputs: curation/c_layer_studies_v2.tsv, data/c_layer_build_report_v2.tsv,
#   data/census_geo_v2.tsv, data/census_adjudication_v2_r3.tsv
# Output: curation/c_layer_studies_rationale_v2.tsv (per study)
# Note: the rationale texts are written statements based on existing
#   records; no speculation. (Texts are in Chinese as compiled internally.)
# =============================================================================
source("scripts/00_setup.R")

run_step("47b_c_studies_rationale", function() {
  cs <- read_tsv("curation/c_layer_studies_v2.tsv")
  rep0 <- read_tsv("data/c_layer_build_report_v2.tsv")
  cen <- read_tsv("data/census_geo_v2.tsv")

  # Per-study written rationale (drafted from c_layer_studies columns + build audit + census records)
  rat <- list(
    GSE130393 = "AD-class independent cluster; retained the normoxia control arm 4/8 (excluded the hypoxia-treated arm). Blocked by the build safeguard: the matrix is a 188-gene panel, only 10 rows remain after dropping missing entries, and 0 of the 56 final-model rule genes are hit - retention cannot yield valid predictions either way, so it remains excluded together with other non-informative clusters.",
    GSE144605 = "AD-class independent cluster; retained the normal-glucose ASC arm 3/12 (excluded the HRMVPC non-MSC control arm and the high-glucose-treated arm). Build passed (30,044 genes).",
    GSE200086 = "AD-class independent cluster; retained the FBS standard-culture arm 4/8 (excluded the SFD serum-free arm). Build passed.",
    GSE212368 = "AD-class; exited primary confirmation in the same batch - GPL32619 is a NanoString 58-gene platform with no evaluable final-model rules; exit reason recorded in the cohort record (2026-10-12).",
    GSE217709 = "AD-class independent cluster; retained age-matched healthy controls 9/38 (excluded the DKD disease arm). Donor accounting: the final report recorded a dash; here 9 healthy controls were recounted from sample-level metadata (see donor_note). Build passed.",
    GSE222263 = "AD-class independent cluster; retained lean-donor L-ASCs 3/6 (excluded the obese O-ASCs arm). Process note: the r2 review had flagged X_SMALL_N ('2 different cells'); restored after sample-level verification (3 donors) - one of the node-A restoration records. Build passed.",
    GSE265844 = "AD-class independent cluster; retained the adherent standard-culture adhered_MSC arm 4/8 (excluded the CD73-sorted isolated_MSC arm). Build passed (transcript-level from Salmon; ENST-to-symbol aggregation).",
    GSE308101 = "AD-class independent cluster; retained lean donors BMI<25 (SUB25) 5/15 (excluded the overweight/obese SUB35/SUB40 arms). Build passed.",
    GSE90149 = "BM-class independent cluster; retained bmMSC healthy donors 11/25 (excluded the placenta fetal/maternal side arms). Donor accounting: the final report recorded 12; sample-level recount from esummary gives bmMSC = 11 (see donor_note). Build passed.",
    GSE111495 = "BM-class independent cluster; retained the FBS(FS20) standard-culture arm 11/22 (excluded the HPGF(PF10) platelet lysate arm). Build passed.",
    GSE129801 = "BM-class independent cluster; retained healthy donors (HD-3~HD-10) 7/22 (excluded the MM/MGUS disease arm). Build passed.",
    GSE14828 = "BM-class independent cluster; retained the no-FGF control arm 9/18 (3 donors x P1/P4/P7; excluded the FGF-treated arm). Build passed.",
    GSE189524 = "BM-class; retained the non-osteoporosis (CF...NO) arm 8/14 (excluded the osteoporosis FO...O arm). Blocked at build: gene_name/gene_id across all 182,434 rows are UCSC transcript IDs (uc*), with no symbol-level annotation and no mapping possible from existing data sources. Remains excluded.",
    GSE24871 = "BM-class independent cluster; retained 7/16 primary BMSC samples (BMSC-1~7); excluded hESC/CD34+/fibroblast control arms. Build passed but coverage is weak (2,284 genes; 2 of the 56 rule genes hit, 7 samples abstained) - disclosed in the build and confirmation workflow and covered by sensitivity analysis.",
    GSE40033 = "BM-class independent cluster; retained healthy volunteers (healthy subject) 4/14 (excluded the psoriasis/fetal arms). Build passed.",
    GSE63099 = "BM-class independent cluster; retained monoculture controls (Control donor 1-6) 6/12 (excluded the HUVEC co-culture arm). Build passed.",
    GSE107214 = "UC-class independent cluster; retained the non-obese group (N1-N7 x triplicate) 21/42 (excluded the obese O1-O7 arm). Build passed.",
    GSE118808 = "UC-class independent cluster; retained UC-MSC (3-6) 4/12 (excluded the amnion/chorionic plate MSC arm). Build passed.",
    GSE165811 = "UC-class independent cluster; retained the hUC.MSC1-32 untreated arm 32/50 (excluded the blank/2IF companion arms of MSC1-3, keeping the final-report definition). Build passed.",
    GSE169533 = "UC-class; retained the control (no DEP) arm 2/4 (excluded the diesel particulate-treated arm). Blocked at build: only 2 samples marked keep, below the >=3-sample build threshold. Remains excluded.",
    GSE18134 = "UC-class; exited primary confirmation in the same batch - the GPL9193 platform covers only 53% of the D intersection; exit reason recorded in the cohort record (2026-10-12).",
    GSE215041 = "UC-class independent cluster; retained WJ-MSC (1-3) 3/12 (excluded the salivary gland MSC and extracellular vesicle arms). Build passed.",
    GSE224190 = "UC-class independent cluster; retained untreated controls (untreated ctrl 1-3) 3/6 (excluded the cytokine-primed arm). Build passed.",
    GSE269586 = "UC-class independent cluster; retained P0 baseline x RNAseq platform 11/46 (excluded the MACE platform and P2/P3 passage arms). Donor accounting: the final report recorded 12+; recount P0 x RNAseq = 11 (see donor_note). exploratory-phase audit (2026-10-12): the v1 L3 actual usage set does not include this dataset; no leakage; retained in primary confirmation.",
    GSE62665 = "UC-class independent cluster; retained the serum-containing standard-culture arm 3/6 (excluded the chemically defined serum-free arm). Build passed.",
    GSE72597 = "UC-class; retained the WJ+CT arm 5/7 (excluded the CPJ arm). Same series as GSE76295 (PMID 27107345); merged and counted as 1 independent cluster. Build passed.",
    GSE76295 = "UC-class; retained the WJ+CT arm 11/16 (excluded the CPJ arm). Same cluster as GSE72597 (same series). Includes Early/Late passage samples; also listed as passage-stability candidates. Build passed.")

  check(setequal(names(rat), cs$gse), "Rationale table and c_layer_studies study set differ: ",
        paste(setdiff(union(names(rat), cs$gse), intersect(names(rat), cs$gse)), collapse = ","))

  decision <- ifelse(cs$gse %in% c("GSE212368", "GSE18134"), "exit_f4",
              ifelse(cs$gse %in% c("GSE130393", "GSE189524", "GSE169533"), "excluded_at_step9",
              ifelse(cs$evaluable_main == "YES", "keep_main", "check")))
  donor_note <- ifelse(cs$gse %in% c("GSE90149", "GSE217709", "GSE269586"),
                       "final-report donor count differs from sample-level recount (see rationale); sample-level recount prevails",
                       ifelse(is.na(cs$n_donor_final_report), "final report recorded no donor count; recount in n_donor_recount",
                              "final report agrees with recount"))
  out <- data.frame(
    gse = cs$gse, source_class = cs$source_class, cluster_group = cs$cluster_group,
    decision = decision, rationale = unlist(rat[cs$gse]), donor_note = donor_note,
    evidence = "c_layer_studies_v2 + c_layer_build_report_v2 + census_geo_v2 + cohort record (2026-10-12)",
    date = "2026-10-13", stringsAsFactors = FALSE)
  check(nrow(out) == nrow(cs) && !anyNA(out$rationale), "rationale table incomplete")
  write_tsv(out, "curation/c_layer_studies_rationale_v2.tsv")
  log_msg("rationale table: ", nrow(out), " rows (keep_main ", sum(decision == "keep_main"),
          " / exit_f4 ", sum(decision == "exit_f4"),
          " / excluded_at_step9 ", sum(decision == "excluded_at_step9"), ")")
  invisible(out)
})
