# =============================================================================
# 21d_audit_extension.R —— extension of the audit to newly curated claims
# (post hoc; the lists were prepared before results were inspected).
# Inputs (read-only): data/within_study_de_v2.rds, data/pc_precomputed_v2.rds,
#   curation/marker_list_ext_v2.tsv (extension claims),
#   data/marker_repro_classification_v2.tsv (for combined statistics)
# Outputs: data/marker_transfer_ext_v2.tsv,
#   data/marker_repro_classification_ext_v2.tsv (same conventions as 21c),
#   data/marker_repro_combined_v2.tsv
# Conventions: reproduced_claimed = directed pc_q < threshold;
#   reproduced_opposite = two-sided pc_q < threshold with direction opposite to
#   the claim; not_reproduced = remaining testable rows; not_testable excluded
#   from proportions. Anti-circularity: discovery GSE excluded per claim.
# =============================================================================
source("scripts/00_setup.R")

run_step("21d_audit_extension", function() {
  de <- readRDS("data/within_study_de_v2.rds")
  pre <- readRDS("data/pc_precomputed_v2.rds")
  claims <- read_tsv("curation/marker_list_ext_v2.tsv")

  mk_t <- expand_claims(claims)
  mt <- family_calls(pc_rows(mk_t, de, pre))
  med_lfc <- stats::aggregate(logFC ~ gene + contrast, data = de, FUN = stats::median)
  names(med_lfc)[3] <- "median_logFC"
  mt <- safe_merge(mt, med_lfc, by = c("gene", "contrast"), type = "m:1")
  write_tsv(mt, "data/marker_transfer_ext_v2.tsv")

  both_t <- mk_t; both_t$direction <- 0L
  both <- family_calls(pc_rows(both_t, de, pre))
  key_mt <- paste(mt$claim_id, mt$contrast, sep = "|")
  key_b <- paste(both$claim_id, both$contrast, sep = "|")
  j <- match(key_mt, key_b)
  check(!any(is.na(j)), "two-sided result keys do not match")
  check(all(mt$gene == both$gene[j]), "keys match but gene differs")

  out <- data.frame(
    claim_id = mt$claim_id, gene = mt$gene, contrast = mt$contrast,
    pmid = claims$pmid[match(mt$claim_id, claims$claim_id)],
    direction_claimed = mt$direction,
    K = mt$K, testable = mt$testable,
    pc_q_directed = mt$pc_q, transfer_directed = mt$transfer,
    dir_used_2sided = both$dir_used[j], pc_q_2sided = both$pc_q[j],
    median_logFC = mt$median_logFC, stringsAsFactors = FALSE)
  out$class <- ifelse(!out$testable, "not_testable",
    ifelse(out$transfer_directed, "reproduced_claimed",
      ifelse(!is.na(out$pc_q_2sided) & out$pc_q_2sided < CFG$pc_q &
             out$dir_used_2sided == -out$direction_claimed,
             "reproduced_opposite", "not_reproduced")))
  write_tsv(out, "data/marker_repro_classification_ext_v2.tsv")

  # per-paper summary (extension part)
  agg <- function(d, id) {
    te <- d[d$testable, , drop = FALSE]
    data.frame(pmid = id, n_rows = nrow(d), n_testable = nrow(te),
               reproduced_claimed = sum(te$class == "reproduced_claimed"),
               reproduced_opposite = sum(te$class == "reproduced_opposite"),
               not_reproduced = sum(te$class == "not_reproduced"),
               stringsAsFactors = FALSE)
  }
  byp <- do.call(rbind, c(lapply(split(out, out$pmid), function(d) agg(d, d$pmid[1])),
                          list(agg(out, "ALL_EXT"))))
  write_tsv(byp, "data/marker_repro_by_paper_ext_v2.tsv")

  # combined statistics (old + new; old table read-only, not recomputed)
  old <- read_tsv("data/marker_repro_classification_v2.tsv")
  old_te <- old[old$class != "not_testable", , drop = FALSE]
  new_te <- out[out$class != "not_testable", , drop = FALSE]
  cnt <- function(d) c(sum(d$class == "reproduced_claimed"),
                       sum(d$class == "reproduced_opposite"),
                       sum(d$class == "not_reproduced"))
  comb <- data.frame(
    scope = c("prior_phase", "extension", "combined"),
    n_testable = c(nrow(old_te), nrow(new_te), nrow(old_te) + nrow(new_te)),
    reproduced_claimed = c(cnt(old_te)[1], cnt(new_te)[1], cnt(old_te)[1] + cnt(new_te)[1]),
    reproduced_opposite = c(cnt(old_te)[2], cnt(new_te)[2], cnt(old_te)[2] + cnt(new_te)[2]),
    not_reproduced = c(cnt(old_te)[3], cnt(new_te)[3], cnt(old_te)[3] + cnt(new_te)[3]),
    stringsAsFactors = FALSE)
  comb$rate_claimed <- round(comb$reproduced_claimed / comb$n_testable, 4)
  comb$rate_opposite <- round(comb$reproduced_opposite / comb$n_testable, 4)
  comb$rate_not <- round(comb$not_reproduced / comb$n_testable, 4)
  write_tsv(comb, "data/marker_repro_combined_v2.tsv")

  te <- out[out$testable, , drop = FALSE]
  log_msg("extension testable ", nrow(te), ": reproduced as claimed ", sum(te$class == "reproduced_claimed"),
          "; reproduced opposite ", sum(te$class == "reproduced_opposite"),
          "; not reproduced ", sum(te$class == "not_reproduced"),
          "(not testable ", sum(out$class == "not_testable"), ")")
  invisible(out)
})
