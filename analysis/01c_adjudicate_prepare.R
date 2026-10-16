# =============================================================================
# 01c_adjudicate_prepare.R —— ledger and batch preparation for the item-by-item
# manual adjudication of all 1,089 candidate datasets.
# Discipline: every candidate receives a manual verdict; automated pre-screening
# provides hints only. Verdicts are recorded per batch in
# docs/adjudication_batches/batch_NN_verdicts.tsv (controlled vocabulary below)
# and merged into the ledger by 01d.
# Controlled vocabulary:
#   C_AD/C_BM/C_UC  single-source validation candidates (sample-level check next)
#   MULTI_NEW   new multi-source independent clusters
#   ORTHO_SCRNA single-cell (orthogonal layer)
#   PASSAGE_SERIES  serial-passage series (stability analyses)
#   X_*         exclusion categories (not MSC / not cultured / not human /
#               no expression matrix / cell line or iPSC / treated without
#               control arm / source outside scope / fetal-embryo / no matrix /
#               too few samples / duplicate series / used in the earlier phase)
#   UNCLEAR     insufficient information
# Outputs: data/census_adjudication_v2.tsv (ledger), docs/adjudication_batches/
#          batch_NN.md (reading material, 40 candidates per batch)
# =============================================================================
source("scripts/00_setup.R")

run_step("01c_adjudicate_prepare", function() {
  cand <- read_tsv("data/source_v1/census_candidates_v1.tsv")
  geo <- read_tsv("data/census_geo_v2.tsv")
  check(!anyDuplicated(cand$accession), "duplicate accession in the candidate pool (deduplicate first)")
  cand$taxon <- ifelse(is.na(cand$taxon), "", cand$taxon)
  cand$gdstype <- ifelse(is.na(cand$gdstype), "", cand$gdstype)

  txt <- tolower(paste(cand$title, cand$summary, sep = " | "))
  human <- cand$taxon == "Homo sapiens"
  expr <- grepl("Expression profiling", cand$gdstype)
  is_msc <- grepl("mesenchymal|stromal|\\bmsc\\b|\\bmscs\\b|\\badsc\\b|\\badscs\\b|\\basc\\b|\\bascs\\b", txt)
  not_culture <- grepl("ipsc|induced pluripotent|cell line|\\bhe-la\\b|\\b hela\\b|tumor|cancer|carcinoma|melanoma|glioma|\\bmouse\\b|\\bmice\\b|\\brat\\b|porcine|bovine|murine|canine|equine|ovine|xenopus|zebrafish|drosophila|c\\. elegans", txt)
  ucb <- grepl("cord blood|\\bucb\\b", txt)
  ad <- grepl("adipose|lipoaspirate|\\badscs?\\b|\\bascs?\\b|subcutaneous fat|liposuction", txt) & !ucb
  bm <- grepl("bone marrow|\\bbm-?msc|marrow stromal|iliac", txt) & !ucb
  uc <- grepl("umbilical|wharton|\\buc-?msc|\\bwj\\b|hucmsc", txt) & !ucb
  ooc <- grepl("dental|periodont|apical papilla|tooth|synovi|tendon|meniscus|placenta|amnio|chorion|fetal liver|fetal bone|embry|cord blood|tonsil|trabecular meshwork|gingiv|pulp|hair follicle|dermal papilla|muscle satellite|heart|cardiac|liver|kidney|lung|skin|dermis|dermal fibro|perioste|nasal|umbilical artery", txt) & !ad & !bm & !uc
  treat <- grepl("hypox|treat|knock|overexpress|stimul|differentiat|osteogen|adipogen|chondrogen|inflammat|g-csf|\\bifn|\\btnf|il-1|aging|senesc|passage|irradiat|serum deprivation|3d |scaffold|hydrogel|matrigel|co-?cult|co-culture", txt)

  n_src <- (ad + bm + uc)
  mv <- ifelse(!human, "X_NOT_HUMAN",
        ifelse(!expr, "X_NOT_EXPR",
        ifelse(cand$accession %in% geo$gse, "X_EXPLORATORY_PHASE",
        ifelse(!is_msc, "X_NOT_MSC",
        ifelse(not_culture, "X_CELLLINE_IPSC",
        ifelse(n_src >= 2, "MULTI_NEW?",
        ifelse(n_src == 1, "C_SINGLE?", "X_NOT_MSC")))))))
  mv[mv == "C_SINGLE?" & ucb & !ad & !bm & !uc] <- "X_SOURCE_OOC"
  mv[mv %in% c("MULTI_NEW?", "C_SINGLE?") & ooc] <- "X_SOURCE_OOC"
  mreason <- ifelse(!human, "taxon≠human",
             ifelse(!expr, "not an expression matrix",
             ifelse(!is_msc, "no MSC/stromal term in title/summary",
             ifelse(not_culture, "contains cell line/iPSC/tumor/animal term",
             ifelse(ooc, "source outside the three classes",
             ifelse(ucb & !ad & !bm & !uc, "umbilical cord blood",
             ifelse(treat, "contains treatment/differentiation/co-culture term (control arm pending check)", "")))))))

  led <- data.frame(
    accession = cand$accession, title = cand$title, taxon = cand$taxon,
    gdstype = cand$gdstype, gpl = cand$gpl, n_samples = cand$n_samples,
    pubmed = cand$pubmed, machine_hint = mv, machine_reason = mreason,
    human_verdict = "", human_reason = "", b_idx = NA_integer_,
    stringsAsFactors = FALSE, row.names = NULL)
  # reading material: the 49 earlier-adjudicated items in one batch (quick pass + reconciliation), batches of 40 otherwise
  dir.create(proj_path("docs", "adjudication_batches"), showWarnings = FALSE)
  type_code <- ifelse(grepl("single-cell|single cell|scRNA", cand$gdstype, ignore.case = TRUE), "SCRNA",
               ifelse(grepl("sequencing", cand$gdstype, ignore.case = TRUE), "SEQ",
               ifelse(grepl("array", cand$gdstype, ignore.case = TRUE), "ARR",
               substr(cand$gdstype, 1, 6))))
  names(type_code) <- cand$accession
  sum_short <- substr(cand$summary, 1, 240); names(sum_short) <- cand$accession
  used_earlier <- led$accession %in% geo$gse
  led$b_idx <- ifelse(used_earlier, 0L, NA_integer_)
  rest <- led[!used_earlier, , drop = FALSE]
  n_batch <- ceiling(nrow(rest) / 40)
  led$b_idx[!used_earlier] <- rep(seq_len(n_batch), each = 40, length.out = nrow(rest))
  for (b in sort(unique(led$b_idx))) {
    d <- led[led$b_idx == b, , drop = FALSE]
    if (nrow(d) == 0) next
    fn <- sprintf("docs/adjudication_batches/batch_%02d.md", b)
    lines <- c(sprintf("# adjudication batch %02d (%d items)", b, nrow(d)),
               "columns: index | accession | taxon | type | GPL | n | machine hint (reason) | title | summary (truncated)",
               "verdict vocabulary: C_AD C_BM C_UC MULTI_NEW ORTHO_SCRNA PASSAGE_SERIES X_NOT_MSC X_NOT_CULTURE X_NOT_HUMAN X_NOT_EXPR X_CELLLINE_IPSC X_TREAT_NO_CTRL X_SOURCE_OOC X_FETAL_EMBRYO X_NO_MATRIX X_SMALL_N X_DUP X_EXPLORATORY_PHASE UNCLEAR", "")
    for (i in seq_len(nrow(d))) {
      a <- d$accession[i]
      lines <- c(lines, sprintf("%02d | %s | %s | %s | GPL%s | n=%s | %s(%s) | %s | %s",
        i, a, substr(d$taxon[i], 1, 4), type_code[[a]], d$gpl[i], d$n_samples[i],
        d$machine_hint[i], substr(d$machine_reason[i], 1, 24),
        substr(d$title[i], 1, 105), sum_short[[a]]))
    }
    writeLines(lines, proj_path(fn))
  }
  write_tsv(led, "data/census_adjudication_v2.tsv")   # ledger written last (includes batch index)
  log_msg("ledger ", nrow(led), " rows; batch 00 (earlier adjudications reconciled, ", sum(used_earlier), " items) + 01..", n_batch,
          "(<=40 each, ", nrow(rest), " awaiting manual adjudication)")
  invisible(led)
})
