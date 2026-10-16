# =============================================================================
# 45_c_layer_manifest.R —— build the validation-cohort sample manifest
# (step 7c of the validation-cohort plan).
# Basis: per validation study, the GSM list, source label, arm, platform,
#   main-evaluability, and the same-cell exclusion rule (same-platform
#   Spearman > 0.98 against the development layer).
# List source: 26 eligible independent groups (AD 8 / BM 8 / UC 10 plus the
#   two-study series 72597/76295); rebuilt 2026-10-12 by verbatim re-check of
#   the original list, without re-adjudication.
# Sample lists: NCBI esummary (db=gds), raw JSON archived in
#   data/c_layer_meta/esummary_{AD,BM,UC}.json; parsed copy
#   data/c_layer_meta/c_layer_gsm_raw_v2.tsv (476 samples).
# Arm rules: mechanical regular expressions on sample titles; re-checked on
#   series-matrix characteristics at build time.
# Outputs: curation/c_layer_manifest_v2.tsv (per-GSM rows;
#   role=keep/exclude/exclude_f4); curation/c_layer_studies_v2.tsv (study level).
# =============================================================================
source("scripts/00_setup.R")

run_step("45_c_layer_manifest", function() {
  raw <- read_tsv("data/c_layer_meta/c_layer_gsm_raw_v2.tsv")
  check(nrow(raw) == 476, "raw GSM table should have 476 rows, actual ", nrow(raw))
  check(all(c("gse", "gsm", "title", "gpl") %in% names(raw)), "raw table is missing columns")

  # ---- study-level rule table (roster = node A final report; reported-number fields = final-report values; recount = counts re-derived here)----
  mk <- function(cls, keep, excl, donor_rep = NA_real_, donor_recount = NA_real_,
                 arm_k, arm_x, cluster = NULL, evaluable = "YES", v1 = "not covered in v1",
                 note = "") {
    list(cls = cls, keep = keep, excl = excl, donor_rep = donor_rep,
         donor_recount = donor_recount, arm_k = arm_k, arm_x = arm_x,
         cluster = cluster, evaluable = evaluable, v1 = v1, note = note)
  }
  rules <- list(
    GSE130393 = mk("AD", "^donor [0-9]+ control$", "hypoxia", 4, 4,
                   "normoxia control (donor 1-4 control)", "hypoxia-treated arm"),
    GSE144605 = mk("AD", "^ASC donor [0-9]+ normal glucose$", "HRMVPC|high glucose", 3, 3,
                   "normal-glucose ASC", "HRMVPC non-MSC control arm + high-glucose treatment arm"),
    GSE200086 = mk("AD", "FBS culture condition$", "SFD", 4, 4,
                   "FBS routine culture", "SFD serum-free arm"),
    GSE212368 = mk("AD", NA_character_, ".*", 5, 5,
                   "(withdrawn from main analysis)", "all (GPL32619 is the NanoString 58-gene panel; no evaluable rules under the final model)",
                   evaluable = "NO", note = "withdrawn from main confirmatory analysis; explained in supplementary material"),
    GSE217709 = mk("AD", "^Control\\.[0-9]+$", "^DKD", NA_real_, 9,
                   "age-matched healthy control (Control.30-45)", "DKD disease arm",
                   note = "final report records donor as '-'; recount found 9 healthy controls"),
    GSE222263 = mk("AD", "^L-ASCs-[0-9]+$", "^O-ASCs-", 3, 3,
                   "lean-donor L-ASCs (1-3)", "obese O-ASCs arm",
                   note = "r2 review flagged X_SMALL_N ('2 different cells'); final report restored after sample-level verification (3 donors)"),
    GSE265844 = mk("AD", "^adhered_MSCs", "^isolated_MSCs", 4, 4,
                   "adherent routine culture adhered_MSC", "CD73-sorted isolated_MSC"),
    GSE308101 = mk("AD", "SUB25", "SUB35|SUB40", 5, 5,
                   "lean donors BMI<25 (SUB25)", "overweight/obese arm (SUB35/SUB40)"),
    GSE90149 = mk("BM", "bone marrow, patient", "placenta", 12, 11,
                  "bmMSC healthy donors (BRD/LQ/WHR/GSK/SDK/BRM/HNC/BEH/ANB/SRF/KKG)",
                  "placenta arm (fetal/maternal side)",
                  note = "final report records 12 donors; sample-level esummary recount gives bmMSC = 11; discrepancy logged for re-check at build time"),
    GSE111495 = mk("BM", "FS20$", "PF10", 5, 5,
                   "FBS(FS20) routine culture", "HPGF(PF10) platelet lysate arm"),
    GSE129801 = mk("BM", "^HD-[0-9]+$", "^MM|^MGUS", 7, 7,
                   "healthy donors (HD-3~HD-10)", "MM/MGUS disease arm"),
    GSE14828 = mk("BM", "^[123]CP[147]$", "FP", 3, 3,
                   "no-FGF control (CP; 3 donors x P1/P4/P7)", "FGF-treated arm (FP)"),
    GSE189524 = mk("BM", "^CF", "^FO", 4, 4,
                   "non-osteoporosis (CF...NO)", "osteoporosis arm (FO...O)"),
    GSE24871 = mk("BM", "^BMSC-[0-9]+", "hESC|CD34|Fibroblast", 7, 7,
                   "primary BMSC (BMSC-1~7)", "hESC/CD34+/fibroblast control"),
    GSE40033 = mk("BM", "healthy subject", "psoriasis|aborted fetus", 4, 4,
                   "healthy volunteers (healthy subject)", "psoriasis/fetal arm"),
    GSE63099 = mk("BM", "^Control donor", "^Treated", 6, 6,
                   "monoculture control (Control donor 1-6)", "HUVEC co-culture arm (Treated)"),
    GSE107214 = mk("UC", "^N[1-7]_[0-9]$", "^O[1-7]_", 7, 7,
                   "non-obese group (N1-N7 x 3 replicates)", "obese arm (O1-O7)"),
    GSE118808 = mk("UC", "^UC-MSC[0-9]+$", "^AM-MSC|^P-MSC", 4, 4,
                   "UC-MSC (3-6)", "amniotic membrane/chorionic plate MSC arm"),
    GSE165811 = mk("UC", "^hUC\\.MSC[0-9]+$", "MSC[0-9]+-(blank|2IF)", 32, 32,
                   "hUC.MSC1-32 untreated", "blank/2IF paired arms of MSC1-3 (kept excluded per the final-report definition)"),
    GSE169533 = mk("UC", "Control", "Diesel", 2, 2,
                   "control (no DEP)", "diesel exhaust particle-treated arm"),
    GSE18134 = mk("UC", NA_character_, ".*", 4, 4,
                   "(withdrawn from main analysis)", "all (GPL9193 covers only 53% of the D intersection)",
                   evaluable = "NO", note = "withdrawn from main confirmatory analysis; explained in supplementary material"),
    GSE215041 = mk("UC", "wharton's jelly mesenchymal stem cells", "salivary|extracellular vesicle", 3, 3,
                   "WJ-MSC (1-3)", "salivary gland MSC and extracellular vesicle arm"),
    GSE224190 = mk("UC", "untreated, ctrl", "primed", 3, 3,
                   "untreated control (untreated ctrl 1-3)", "cytokine-primed arm (primed)"),
    GSE269586 = mk("UC", "passage 0.*RNAseq", "MACE|passage [23]", 12, 11,
                   "P0 baseline x RNAseq platform (donor 190-206)", "MACE platform and P2/P3 passage arms",
                   v1 = "not used in v1 (verified)",
                   note = paste("final report records 12+ donors; recount P0 x RNAseq = 11. ",
                                "exploratory-phase verification (2026-10-12): l2_samples_v1 does not contain this series; ",
                                "v1 L3 actually-used set = {139073,275980,291956,309666}; v1 data_raw did not download this series; ",
                                "X_EXPLORATORY_PHASE was mislabeled by the mechanical rule in 01c as 'listed in the v1 final table'; no leakage; kept in the main confirmatory analysis")),
    GSE62665 = mk("UC", "serum containing", "serum free", 3, 3,
                  "serum-containing routine culture", "chemically defined serum-free arm"),
    GSE72597 = mk("UC", "Wharton's jelly|Cord Tissue", "Cord-placenta junction", 4, 4,
                  "WJ+CT arm", "CPJ arm",
                  cluster = "G20",
                  note = "same series as GSE76295 (PMID 27107345); merged as 1 independent cluster"),
    GSE76295 = mk("UC", "Wharton's jelly|Cord Tissue", "Cord-placenta junction", 4, 4,
                  "WJ+CT arm (includes Early/Late passage)", "CPJ arm",
                  cluster = "G20",
                  note = "same cluster (see GSE72597); includes Late-passage samples, also listed as stability candidates")
  )
  check(setequal(names(rules), unique(raw$gse)), "rule table and GSM table study sets differ")
  check(length(rules) == 27, "rule table should contain 27 studies")

  # ---- sample-level assignment (keep takes priority; stop on any unmatched sample - fail-safe)---------------------------------
  role_of <- function(title, keep, excl) {
    if (is.na(keep)) return(c("exclude_f4", "F4_exit"))
    if (grepl(keep, title)) return(c("keep", "keep_rule"))
    if (grepl(excl, title)) return(c("exclude", "arm_excl"))
    c("review", "no_rule_match")
  }
  smp <- do.call(rbind, lapply(seq_len(nrow(raw)), function(i) {
    r <- raw[i, ]; ru <- rules[[r$gse]]
    rr <- role_of(r$title, ru$keep, ru$excl)
    data.frame(gse = r$gse, gsm = r$gsm, source_class = ru$cls, platform = r$gpl,
               sample_title = r$title, role = rr[1], rule_id = rr[2],
               stringsAsFactors = FALSE)
  }))
  bad <- smp[smp$role == "review", , drop = FALSE]
  check(nrow(bad) == 0, "samples not matched to any arm rule: ", paste(bad$gsm, collapse = ","))
  check(sum(smp$role == "keep") == 190, "total keep should be 190, actual ", sum(smp$role == "keep"))
  smp <- smp[order(match(smp$gse, names(rules)), smp$gsm), ]

  # ---- study-level summary ---------------------------------------------------------------
  stu <- do.call(rbind, lapply(names(rules), function(g) {
    d <- smp[smp$gse == g, ]; ru <- rules[[g]]
    data.frame(
      gse = g, source_class = ru$cls,
      cluster_group = if (is.null(ru$cluster)) g else ru$cluster,
      platform_gpl = d$platform[1], n_gsm = nrow(d),
      n_keep = sum(d$role == "keep"), n_exclude = sum(d$role == "exclude"),
      n_exclude_f4 = sum(d$role == "exclude_f4"),
      n_donor_final_report = ru$donor_rep, n_donor_recount = ru$donor_recount,
      arm_keep_rule = ru$arm_k, arm_exclude_rule = ru$arm_x,
      evaluable_main = ru$evaluable, v1_check = ru$v1,
      same_cell_exclusion_rule =
        "applied at build time: same-platform samples with Spearman > 0.98 to layer D are treated as the same cells and removed",
      notes = ru$note, stringsAsFactors = FALSE)
  }))
  write_tsv(smp, "curation/c_layer_manifest_v2.tsv")
  write_tsv(stu, "curation/c_layer_studies_v2.tsv")

  # ---- reconciliation output -------------------------------------------------------------------
  log_msg("sample-level ", nrow(smp), " rows: keep ", sum(smp$role == "keep"),
          "; exclude ", sum(smp$role == "exclude"),
          "; exclude_f4 ", sum(smp$role == "exclude_f4"))
  print(as.data.frame(table(smp$source_class[smp$role == "keep"])))
  invisible(list(samples = smp, studies = stu))
})
