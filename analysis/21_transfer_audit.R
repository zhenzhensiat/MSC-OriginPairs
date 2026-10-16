# =============================================================================
# 21_transfer_audit.R —— replication audit of published source markers
# (development layer only).
# Inputs: data/studies/*.rds (layer D); curation/marker_list.tsv
#   (claim_id, gene, high_source, low_source, discovery_gse, pmid);
#   curation/hox_list.tsv; data/raw/hgnc_complete_set.txt (gene_group column)
# Outputs: within_study_de_v2.rds (full), within_study_de_v2_q05.tsv (trimmed);
#   hox_transfer_v2.tsv, marker_transfer_v2.tsv, marker_measurability_v2.tsv,
#   null_summary_v2.tsv, family_ranking_v2.tsv (+ raw null-draw RDS)
# =============================================================================
source("scripts/00_setup.R")
need_pkgs(c("limma", "edgeR"))

run_step("21_transfer_audit", function() {

  # ---- 1. Within-study differential expression --------------------------------------------------------
  studies <- load_studies("D")
  de_list <- lapply(studies, function(s) { log_msg("within-study DE: ", s$gse); within_study_de(s) })
  de <- do.call(rbind, Filter(Negate(is.null), de_list))
  check(NROW(de) > 0, "within-study differential expression produced no output")
  saveRDS(de, proj_path("data", "within_study_de_v2.rds"))
  write_tsv(de[de$q_within < 0.05, ], "data/within_study_de_v2_q05.tsv")

  pre <- precompute_pc(de)
  saveRDS(pre, proj_path("data", "pc_precomputed_v2.rds"))
  feat <- pre[, c("gene", "contrast", "K", "expr_pct", "detect")]
  med_lfc <- stats::aggregate(logFC ~ gene + contrast, data = de, FUN = stats::median)
  names(med_lfc)[3] <- "median_logFC"

  # ---- 2. HOX: no fixed direction (two-sided min x 2), within-study BH -------------------------------
  hox <- read_tsv("curation/hox_list.tsv")$gene
  hox_t <- expand.grid(gene = hox, contrast = CFG$contrasts, stringsAsFactors = FALSE)
  hox_t$direction <- 0L; hox_t$exclude <- ""
  hox_res <- family_calls(pc_rows(hox_t, de, pre))
  hox_res <- safe_merge(hox_res, med_lfc, by = c("gene", "contrast"), type = "m:1")
  write_tsv(hox_res, "data/hox_transfer_v2.tsv")

  # ---- 3. published markers: claimed direction (one-sided), excluding the original discovery study (anti-circularity)----------------
  mk_t <- expand_claims(read_tsv("curation/marker_list.tsv"))
  mk_res <- family_calls(pc_rows(mk_t, de, pre))
  mk_res <- safe_merge(mk_res, med_lfc, by = c("gene", "contrast"), type = "m:1")
  write_tsv(mk_res, "data/marker_transfer_v2.tsv")

  # testability reported separately (R12): absent = not detected in any study; low_K = detected but K<2
  meas <- data.frame(
    group = c("published", "HOX"),
    n_rows = c(nrow(mk_res), nrow(hox_res)),
    absent = c(sum(mk_res$K == 0), sum(hox_res$K == 0)),
    low_K = c(sum(mk_res$K > 0 & !mk_res$testable), sum(hox_res$K > 0 & !hox_res$testable)),
    testable = c(sum(mk_res$testable), sum(hox_res$testable)),
    transfer = c(sum(mk_res$transfer), sum(hox_res$transfer)),
    stringsAsFactors = FALSE)
  meas$rate_among_testable <- meas$transfer / meas$testable
  meas$median_abs_logFC_transferred <- c(
    stats::median(abs(mk_res$median_logFC[mk_res$transfer])),
    stats::median(abs(hox_res$median_logFC[hox_res$transfer])))
  write_tsv(meas, "data/marker_measurability_v2.tsv")

  # ---- 4. matched null controls (expression decile x detection quintile x K)-------------------------------
  in_feat <- function(t) t[paste(t$gene, t$contrast) %in% paste(feat$gene, feat$contrast), , drop = FALSE]
  hox_in <- in_feat(hox_t); mk_in <- in_feat(mk_t)
  null_file <- proj_path("data", "null_draws_v2.rds")
  if (!file.exists(null_file)) {   # sample and save to disk first, then compute
    saveRDS(list(hox = sample_matched_null(hox_in, feat, CFG$n_null_sets, derive_seed(21, 0)),
                 published = sample_matched_null(mk_in, feat, CFG$n_null_sets, derive_seed(21, 100000L))),
            null_file)
  }
  draws <- readRDS(null_file)
  obs <- c(hox = transfer_rate(family_calls(pc_rows(hox_in, de, pre))),
           published = transfer_rate(family_calls(pc_rows(mk_in, de, pre))))
  nulls <- list(hox = null_rates(draws$hox, hox_in, de, pre),
                published = null_rates(draws$published, mk_in, de, pre))
  ns <- do.call(rbind, lapply(names(obs), function(g) {
    v <- nulls[[g]]
    data.frame(group = g, observed_rate = obs[[g]], null_mean = mean(v, na.rm = TRUE),
               null_sd = stats::sd(v, na.rm = TRUE), null_q95 = unname(stats::quantile(v, 0.95, na.rm = TRUE)),
               empirical_p = empirical_p(obs[[g]], v),
               max_relax_level = max(draws[[g]]$relax_level), n_sets = length(v), stringsAsFactors = FALSE)
  }))
  write_tsv(ns, "data/null_summary_v2.tsv")

  # ---- 5. gene family ranking (HGNC gene groups; z against each family's own matched null distribution)-----------
  # HGNC site returned HTTP 500 that day: this step deferred (does not affect core audit results)
  hgnc_file <- proj_path("data", "raw", "hgnc_complete_set.txt")
  if (!file.exists(hgnc_file)) {
    log_msg("!! HGNC complete set missing (genenames.org returned 500 that day); family ranking skipped,",
            "rerun 22_family_ranking.R once restored")
    write_tsv(data.frame(message = "HGNC file unavailable; family ranking deferred"),
              "data/family_ranking_v2.tsv")
    return(invisible(list(hox = hox_res, published = mk_res, nulls = ns)))
  }
  fam <- hgnc_families(proj_path("data", "raw", "hgnc_complete_set.txt"))
  fam <- rbind(fam, data.frame(family = "HOX_39", gene = hox, stringsAsFactors = FALSE))
  fam_list <- split(fam$gene, fam$family)
  pre_key <- paste(pre$gene, pre$contrast, sep = "|")
  fr <- do.call(rbind, lapply(seq_along(fam_list), function(i) {
    fn <- names(fam_list)[i]
    tg <- expand.grid(gene = unique(fam_list[[fn]]), contrast = CFG$contrasts, stringsAsFactors = FALSE)
    tg$direction <- 0L; tg$exclude <- ""
    o <- family_calls(pc_either_fast(tg$gene, tg$contrast, pre, pre_key))
    frac <- mean(o$testable)
    base <- data.frame(family = fn, n_genes = length(unique(tg$gene)), frac_testable = frac,
                       rate = NA_real_, null_mean = NA_real_, null_sd = NA_real_, z = NA_real_,
                       empirical_p = NA_real_, included = FALSE, stringsAsFactors = FALSE)
    if (frac < CFG$family_testable_min) return(base)
    tg_in <- in_feat(tg)
    dr <- sample_matched_null(tg_in, feat, CFG$n_family_null_sets, derive_seed(22, i * CFG$n_family_null_sets))
    v <- null_rates(dr, tg_in, de, pre)
    base$rate <- transfer_rate(o)
    base$null_mean <- mean(v, na.rm = TRUE); base$null_sd <- stats::sd(v, na.rm = TRUE)
    base$z <- (base$rate - base$null_mean) / base$null_sd
    base$empirical_p <- empirical_p(base$rate, v)
    base$included <- TRUE
    base
  }))
  fr <- fr[order(-fr$included, -fr$z), ]
  fr$rank_z <- NA_integer_
  fr$rank_z[fr$included] <- seq_len(sum(fr$included))
  write_tsv(fr, "data/family_ranking_v2.tsv")
  log_msg("HOX_39 rank: ", fr$rank_z[fr$family == "HOX_39"], " / ", sum(fr$included))

  # TODO: like-for-like comparison of HOX vs published markers on the same study x contrast sets (R12 sensitivity)
  invisible(list(hox = hox_res, published = mk_res, nulls = ns, families = fr))
})
