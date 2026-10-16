# =============================================================================
# 01b_retriage_candidates.R —— screening of the earlier candidate pool for
# validation-cohort purposes (metadata only; no expression data touched).
# Background: the earlier census adjudicated only the multi-source candidates;
# single-source studies were screened here on the 1,089-candidate pool
# (titles/summaries), producing a shortlist for manual review.
# Rules (plain text, in the spirit of config/archs4_rules.tsv):
#   MSC culture terms x source terms -> single/multi source; treatment designs
#   flagged control_arm_possible; tissue (non-culture)/cell line/iPSC/animal -> out
# Outputs: data/census_retriage_v2.tsv (full classification),
#          docs/retriage_review.md (shortlist for manual review)
# =============================================================================
source("scripts/00_setup.R")

run_step("01b_retriage_candidates", function() {
  v1_dir <- file.path(normalizePath(file.path(getOption("msc.root"), "..", "..")),
                      "\u7814\u7a76\u65b9\u54112")
  src <- file.path(v1_dir, "data", "census_candidates_v1.tsv")
  check(file.exists(src), "v1 candidate pool missing")
  if (!file.exists(proj_path("data", "source_v1", "census_candidates_v1.tsv")))
    file.copy(src, proj_path("data", "source_v1"))
  p <- "data/source_v1/census_candidates_v1.tsv"

  cand <- read_tsv(p)
  geo <- read_tsv("data/census_geo_v2.tsv")
  cand <- cand[!(cand$accession %in% geo$gse) & cand$taxon == "Homo sapiens" &
               grepl("Expression profiling", cand$gdstype), , drop = FALSE]
  txt <- tolower(paste(cand$title, cand$summary, sep = " | "))

  # culture flag: MSC/stromal cell culture (excludes primary tissue, iPSC, cell lines, tumor, animal)
  is_msc <- grepl("mesenchymal (stem|stromal)|\\bmscs?\\b|\\badscs?\\b|\\bascs?\\b|stromal cell", txt)
  not_culture <- grepl("ipsc|induced pluripotent|cell line|\\bhe-la\\b|tumor|cancer|carcinoma|melanoma|adipocyte differentiation of|mouse|mice|rat |porcine|bovine|murine", txt)
  no_msc_word <- !grepl("mesenchymal|stromal|\\bmsc\\b|\\bmscs\\b|\\badsc\\b|\\badscs\\b", txt)

  ad <- grepl("adipose|lipoaspirate|\\badscs?\\b|\\bascs?\\b|subcutaneous fat", txt)
  bm <- grepl("bone marrow|\\bbm-?msc|marrow stromal|iliac", txt)
  uc <- grepl("umbilical|wharton|\\buc-?msc|\\bwj\\b|hucmsc|cord blood", txt)
  ucb_only <- grepl("cord blood|\\bucb\\b", txt) & !uc

  n_src <- (ad & !ucb_only) + (bm & !ucb_only) + (uc & !ucb_only)
  cand$source_hint <- ifelse(n_src >= 2, "multi",
                      ifelse(n_src == 1, ifelse(ucb_only, "UCB", ifelse(ad & !ucb_only, "AD", ifelse(bm & !ucb_only, "BM", "UC"))), "none"))
  cand$triage <- ifelse(!is_msc | no_msc_word | not_culture, "out_not_culture",
                 ifelse(cand$source_hint == "none", "unclear_source",
                 ifelse(cand$source_hint == "UCB", "out_umbilical_cord_blood",
                 ifelse(cand$source_hint == "multi", "candidate_multi", "candidate_single"))))
  # treatment-design flag (control arm possibly recoverable: hypoxia/treatment studies often have normoxia/untreated controls)
  cand$control_arm_possible <- ifelse(cand$triage %in% c("candidate_single", "candidate_multi") &
    grepl("hypox|treat|knock|overexpress|stimul|differentiat|osteogen|adipogen|chondrogen|inflammat|g-csf|ifn|tnf|il-1", txt),
    TRUE, FALSE)

  keep <- cand[cand$triage %in% c("candidate_single", "candidate_multi", "unclear_source"), , drop = FALSE]
  write_tsv(cand[, c("accession", "title", "gdstype", "gpl", "n_samples", "pubmed",
                     "source_hint", "triage", "control_arm_possible")],
            "data/census_retriage_v2.tsv")

  # manual-review shortlist (markdown, grouped by class, titles truncated)
  lines <- c("# C-layer re-screening shortlist (for manual review; metadata only)", "",
             "source: candidates from the 1,089-item pool that were never manually adjudicated; rules in the script header.", "")
  for (cl in c("candidate_single", "candidate_multi", "unclear_source")) {
    sub_k <- keep[keep$triage == cl, , drop = FALSE]
    lines <- c(lines, paste0("## ", cl, " (", nrow(sub_k), ")"), "")
    if (nrow(sub_k) > 0) for (i in seq_len(nrow(sub_k))) {
      lines <- c(lines, sprintf("- **%s** [%s, GPL%s, n=%s, ctrl_arm_possible=%s] %s — PMID:%s",
        sub_k$accession[i], sub_k$source_hint[i], sub_k$gpl[i], sub_k$n_samples[i],
        sub_k$control_arm_possible[i], substr(sub_k$title[i], 1, 110), sub_k$pubmed[i]))
    }
    lines <- c(lines, "")
  }
  writeLines(lines, proj_path("docs", "retriage_review.md"))
  log_msg("re-screening: single-source candidates ", sum(keep$triage == "candidate_single"),
          ", multi-source candidates ", sum(keep$triage == "candidate_multi"),
          ", unclear source ", sum(keep$triage == "unclear_source"),
          "(of which control arm possibly recoverable: ",
          sum(keep$control_arm_possible), ")")
  invisible(keep)
})
