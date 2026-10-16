# =============================================================================
# 01d_adjudicate_merge.R —— merge the 27 manual-adjudication batches into the
# ledger and write it.
# Inputs: data/census_adjudication_v2.tsv + docs/adjudication_batches/batch_*_verdicts.tsv
# Checks: row counts per batch file; verdicts in the controlled vocabulary; matching accessions
# Outputs: updated ledger (human_verdict/human_reason); docs/adjudication_summary.md
# =============================================================================
source("scripts/00_setup.R")

run_step("01d_adjudicate_merge", function() {
  led <- read_tsv("data/census_adjudication_v2.tsv")
  VALID <- c("C_AD","C_BM","C_UC","MULTI_NEW","ORTHO_SCRNA","PASSAGE_SERIES","X_NOT_MSC",
             "X_NOT_CULTURE","X_NOT_HUMAN","X_NOT_EXPR","X_CELLLINE_IPSC","X_TREAT_NO_CTRL",
             "X_SOURCE_OOC","X_FETAL_EMBRYO","X_NO_MATRIX","X_SMALL_N","X_DUP","X_EXPLORATORY_PHASE","UNCLEAR")
  vs <- list.files(proj_path("docs", "adjudication_batches"), "^batch_[0-9]+_verdicts\\.tsv$")
  vs <- vs[order(as.integer(sub("^batch_([0-9]+)_verdicts\\.tsv$", "\\1", vs)))]
  allv <- do.call(rbind, lapply(vs, function(f) {
    d <- read_tsv(file.path("docs", "adjudication_batches", f))
    d$batch <- sub("^batch_([0-9]+)_verdicts\\.tsv$", "\\1", f)
    d
  }))
  check(!anyDuplicated(allv$accession), "duplicate accession across verdict files")
  check(setequal(allv$accession, led$accession),
        "verdict files and ledger differ: missing ", paste(setdiff(led$accession, allv$accession), collapse = ","),
        "; extra ", paste(setdiff(allv$accession, led$accession), collapse = ","))
  check(all(allv$verdict %in% VALID), "invalid verdict term: ",
        paste(unique(allv$verdict[!allv$verdict %in% VALID]), collapse = ","))

  m <- match(led$accession, allv$accession)
  led$human_verdict <- allv$verdict[m]
  led$human_reason <- allv$reason[m]
  write_tsv(led, "data/census_adjudication_v2.tsv")

  tab <- sort(table(led$human_verdict), decreasing = TRUE)
  lines <- c("# manual adjudication of the 1,089-item pool", "",
             paste0("date: 2026-10-10; method: batchwise reading of titles and summaries; batch files ",
                    length(vs), " files in docs/adjudication_batches/."), "",
             "## Verdict distribution", "", "| Verdict | Count |", "|---|---|",
             vapply(names(tab), function(k) paste0("| ", k, " | ", tab[[k]], " |"), character(1)), "")
  cc <- led[led$human_verdict %in% c("C_AD","C_BM","C_UC"), ]
  lines <- c(lines, "", "## C-layer candidates (Phase 1 sample-level verification required)", "",
             vapply(seq_len(nrow(cc)), function(i)
               paste0("- ", cc$accession[i], " (", cc$human_verdict[i], ") ", cc$human_reason[i]), character(1)), "")
  mu <- led[led$human_verdict == "MULTI_NEW", ]
  if (nrow(mu) > 0) lines <- c(lines, "", "## New multi-source", "",
    vapply(seq_len(nrow(mu)), function(i) paste0("- ", mu$accession[i], ": ", mu$human_reason[i]), character(1)), "")
  un <- led[led$human_verdict == "UNCLEAR", ]
  lines <- c(lines, "", paste0("## UNCLEAR (", nrow(un), " items; control arm / composition pending sample-level check)"), "",
    vapply(seq_len(nrow(un)), function(i) paste0("- ", un$accession[i], ": ", un$human_reason[i]), character(1)), "")
  ps <- led[led$human_verdict == "PASSAGE_SERIES", ]
  lines <- c(lines, "", "## Passage series (Phase 6 candidates)", "",
    vapply(seq_len(nrow(ps)), function(i) paste0("- ", ps$accession[i], ": ", ps$human_reason[i]), character(1)), "")
  writeLines(lines, proj_path("docs", "adjudication_summary.md"))
  log_msg("ledger merged: ", nrow(led), " rows; C-layer candidates ",
          sum(led$human_verdict %in% c("C_AD","C_BM","C_UC")), " (AD ",
          sum(led$human_verdict == "C_AD"), " / BM ", sum(led$human_verdict == "C_BM"),
          " / UC ", sum(led$human_verdict == "C_UC"), "); UNCLEAR ",
          sum(led$human_verdict == "UNCLEAR"))
  invisible(led)
})
