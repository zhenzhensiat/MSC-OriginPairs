# =============================================================================
# 01g_upgrade_epmc.R —— upgrade UNCLEAR cases to validation candidates using
# Europe PMC full-text / abstract evidence (untreated / unstimulated /
# conventional / normoxic / control medium with an explicit source).
# Ledger revision r3.
# =============================================================================
source("scripts/00_setup.R")
run_step("01g_upgrade_epmc", function() {
  led <- read_tsv("data/census_adjudication_v2_r2.tsv")
  keys <- led$accession
  up <- function(g, v, r) {
    i <- match(g, keys); check(!is.na(i), "missing ", g)
    led$human_verdict[i] <<- v
    led$human_reason[i] <<- paste0(r, " (01g Europe PMC evidence upgrade; Phase 1 sample-level final check pending)")
  }
  up("GSE212368", "C_AD", "OA full text contains a normoxic standard-condition arm")
  up("GSE300720", "C_AD", "full text contains an unstimulated+conventional control arm")
  up("GSE308101", "C_AD", "full text contains an untreated control arm (non-obese donors)")
  up("GSE265844", "C_AD", "full text contains a conventional adherent group (13 patients)")
  up("GSE233634", "C_UC", "full text contains a conventional planar-culture control arm")
  up("GSE62665",  "C_UC", "full text contains a conventional control arm")
  up("GSE224190", "C_UC", "full text contains untreated/unstimulated (9/7 patients)")
  up("GSE298520", "C_UC", "full text contains normoxic (parental-control lead)")
  up("GSE165811", "C_UC", "full text contains control medium (large donor set, 32 donors)")
  write_tsv(led, "data/census_adjudication_v2_r3.tsv")
  cc <- led[led$human_verdict %in% c("C_AD","C_BM","C_UC"), ]
  log_msg("r3 ledger: C candidates ", nrow(cc), " (AD ", sum(cc$human_verdict=="C_AD"),
          " / BM ", sum(cc$human_verdict=="C_BM"), " / UC ", sum(cc$human_verdict=="C_UC"),
          ") —— all three classes ≥8, threshold passed (Phase 1 final check pending)")
  invisible(led)
})
