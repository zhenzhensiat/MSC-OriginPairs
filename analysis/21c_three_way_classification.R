# =============================================================================
# 21c_three_way_classification.R —— three-way classification of the published
# markers: (1) replicated in the claimed direction; (2) replicated in the
# opposite direction (two-sided q<0.05 with direction opposite to the claim);
# (3) not replicated; untestable rows reported separately.
# Precondition: direction check of the Lee 2026 claims found the recorded
# direction matches the original text; the marker table stands unchanged.
# This table is a post hoc supplement.
# Inputs (read-only): data/within_study_de_v2.rds, data/pc_precomputed_v2.rds,
#   curation/marker_list.tsv, data/marker_transfer_v2.tsv
# Outputs: data/marker_repro_classification_v2.tsv,
#   data/marker_repro_by_paper_v2.tsv
# =============================================================================
source("scripts/00_setup.R")

run_step("21c_three_way", function() {
  de <- readRDS("data/within_study_de_v2.rds")
  pre <- readRDS("data/pc_precomputed_v2.rds")
  mt <- read_tsv("data/marker_transfer_v2.tsv")          # primary definition, row by row
  claims <- read_tsv("curation/marker_list.tsv")
  mk_t <- expand_claims(claims)

  # two-sided test (direction=0): same data and exclusion logic as the primary definition, only the direction setting changes
  both_t <- mk_t; both_t$direction <- 0L
  both <- family_calls(pc_rows(both_t, de, pre))

  # row-by-row alignment (claim_id x contrast unique; keys verified against the result table)
  key_mt <- paste(mt$claim_id, mt$contrast, sep = "|")
  key_b <- paste(both$claim_id, both$contrast, sep = "|")
  j <- match(key_mt, key_b)
  check(!any(is.na(j)), "two-sided results do not match the primary-table keys")
  check(all(mt$gene == both$gene[j]), "keys match but gene differs")   # both$direction is always 0 (control-arm setting)

  out <- data.frame(
    claim_id = mt$claim_id, gene = mt$gene, contrast = mt$contrast,
    pmid = claims$pmid[match(mt$claim_id, claims$claim_id)],
    direction_claimed = mt$direction,
    K = mt$K, testable = mt$testable,
    pc_q_directed = mt$pc_q, transfer_directed = mt$transfer,
    dir_used_2sided = both$dir_used[j], pc_q_2sided = both$pc_q[j],
    median_logFC = mt$median_logFC,
    stringsAsFactors = FALSE)
  out$class <- ifelse(!out$testable, "not_testable",
    ifelse(out$transfer_directed, "reproduced_claimed",
      ifelse(!is.na(out$pc_q_2sided) & out$pc_q_2sided < CFG$pc_q &
             out$dir_used_2sided == -out$direction_claimed,
             "reproduced_opposite", "not_reproduced")))
  write_tsv(out, "data/marker_repro_classification_v2.tsv")

  # summarize by source paper (pmid; includes the ALL total row)
  agg <- function(d, id) {
    te <- d[d$testable, , drop = FALSE]
    data.frame(pmid = id, n_rows = nrow(d), n_testable = nrow(te),
               reproduced_claimed = sum(te$class == "reproduced_claimed"),
               reproduced_opposite = sum(te$class == "reproduced_opposite"),
               not_reproduced = sum(te$class == "not_reproduced"),
               stringsAsFactors = FALSE)
  }
  byp <- do.call(rbind, c(lapply(split(out, out$pmid), function(d) agg(d, d$pmid[1])),
                          list(agg(out, "ALL"))))
  byp$rate_claimed <- round(byp$reproduced_claimed / byp$n_testable, 4)
  byp$rate_opposite <- round(byp$reproduced_opposite / byp$n_testable, 4)
  write_tsv(byp, "data/marker_repro_by_paper_v2.tsv")

  # ---- number check: 29/99 median logFC opposite to the claim; 52.5-36.4 ~ 16 reversed rows
  te <- out[out$testable, ]
  n_logfc_opp <- sum(sign(te$median_logFC) == -te$direction_claimed)
  log_msg("testable rows ", nrow(te), ": reproduced in the claimed direction ", sum(te$class == "reproduced_claimed"),
          "; reproduced opposite ", sum(te$class == "reproduced_opposite"),
          "; not reproduced ", sum(te$class == "not_reproduced"))
  log_msg("rows with median logFC opposite to the claimed direction: ", n_logfc_opp, " / ", nrow(te))
  a <- byp[byp$pmid == "ALL", ]
  log_msg("per-paper summary written: ", nrow(byp) - 1, " papers + ALL; ALL row = see figures")
  invisible(out)
})
