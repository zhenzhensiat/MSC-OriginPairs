# =============================================================================
# 01e_adjudicate_revise.R —— re-check pass over the 1,089 adjudications.
# Method: (1) full-text re-reading of the C/MULTI candidates originally assessed
# from truncated 240-character summaries; (2) re-check of exclusion candidates
# whose text contains control-arm keywords x source terms; (3) superseries /
# subseries consistency checks.
# Correction list (per full-text evidence; certainty labelled):
#   definite: GSE194295 C_BM -> X_NOT_EXPR (MethylCap-seq, not expression)
#             GSE194303 C_BM -> UNCLEAR (superseries without own data)
#             GSE159137 C_AD -> PASSAGE_SERIES (serial passage P0/P1/P2)
#             GSE222263 C_AD -> X_SMALL_N ("2 different cells", 1 donor each)
#   control-arm leads upgraded to UNCLEAR (full texts describe untreated /
#   parental / unstimulated arms): GSE2776, GSE261358, GSE90131, GSE224190,
#   GSE298520, GSE300720
#   superseries reconciliation UNCLEAR -> X_DUP: GSE35959, GSE52845, GSE87798, GSE124862
#   note corrections: GSE118808, GSE100748/49, GSE130393, GSE72597/76295, GSE166327
# Outputs: ledger updated and recorded
# =============================================================================
source("scripts/00_setup.R")

run_step("01e_adjudicate_revise", function() {
  led <- read_tsv("data/census_adjudication_v2.tsv")
  keys <- led$accession
  upd <- function(g, v, r) {
    i <- match(g, keys)
    check(!is.na(i), "not found in the ledger: ", g)
    led$human_verdict[i] <<- v
    led$human_reason[i] <<- paste0(r, "")
  }
  # definite corrections
  upd("GSE194295", "X_NOT_EXPR", "MethylCap-seq methylation subseries, not an expression matrix (expression part: GSE194298)")
  upd("GSE194303", "UNCLEAR", "superseries without own data; expression data belong to subseries GSE194298")
  upd("GSE159137", "PASSAGE_SERIES", "ADMSC serial passage P0/P1/P2 in triplicate (main passage series)")
  upd("GSE222263", "X_SMALL_N", "RNA-seq of 2 different cells: obese and healthy with presumably 1 donor each, below the donor threshold")
  # control-arm leads upgraded to UNCLEAR
  upd("GSE2776",   "UNCLEAR", "chemical neural differentiation vs untreated MSC; untreated control arm pending check, recoverable BM")
  upd("GSE261358", "UNCLEAR", "fibrous dysplasia BMSC; healthy-control BMSC arm pending check")
  upd("GSE90131",  "UNCLEAR", "fucosylation treatment vs untreated BMSC; untreated arm pending check")
  upd("GSE224190", "UNCLEAR", "cytokine pretreatment vs untreated hUC-MSC; untreated arm pending check")
  upd("GSE298520", "UNCLEAR", "engineered promoter vs parental UC-MSC; parental control arm pending check")
  upd("GSE300720", "UNCLEAR", "electrical stimulation vs unstimulated AD-MSC; unstimulated control arm pending check")
  # superseries reconciliation
  upd("GSE35959", "X_DUP", "superseries reconciliation (subseries 35955-58 each adjudicated separately)")
  upd("GSE52845", "X_DUP", "superseries reconciliation (subseries 52841-44 each adjudicated separately)")
  upd("GSE87798", "X_DUP", "superseries reconciliation (subseries 87796/97 each adjudicated separately)")
  upd("GSE124862", "X_DUP", "superseries reconciliation (subseries 124858-60 each adjudicated separately)")
  # note correction (verdicts unchanged)
  note <- function(g, add) {
    i <- match(g, keys)
    led$human_reason[i] <<- paste0(led$human_reason[i], "; ", add)
  }
  note("GSE118808", "all three sources cultured serum-free (medium=xeno_free records)")
  note("GSE100748", "recovered arm = undifferentiated control (undifferentiated controls)")
  note("GSE100749", "recovered arm = undifferentiated control")
  note("GSE130393", "existence of a normoxia control arm pending check (not stated in the abstract)")
  note("GSE72597", "cord regions CPJ/CT/WJ heterogeneous; use requires subregion stratification")
  note("GSE76295", "cord regions CPJ/CT/WJ heterogeneous")
  note("GSE166327", "only 2 donors per class (below the ≥3-donor threshold; cluster registration only)")

  out <- "data/census_adjudication_v2_r2.tsv"   # overwrite guard: revisions stored as r2; v1 kept as history
  write_tsv(led, out)
  cc <- led[led$human_verdict %in% c("C_AD", "C_BM", "C_UC"), ]
  log_msg("re-check corrections complete: C-layer candidates ", nrow(cc), " (AD ", sum(cc$human_verdict == "C_AD"),
          " / BM ", sum(cc$human_verdict == "C_BM"), " / UC ", sum(cc$human_verdict == "C_UC"),
          "); UNCLEAR ", sum(led$human_verdict == "UNCLEAR"),
          "; X_DUP ", sum(led$human_verdict == "X_DUP"))
  invisible(led)
})
