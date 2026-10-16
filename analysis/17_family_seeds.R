# =============================================================================
# 17_family_seeds.R —— Monte-Carlo robustness of the gene-family ranking.
# Basis: the family-ranking seeds are kept fixed; robustness is shown by
# repeating the ranking with several seed sets and reporting the rank range.
# Design: same logic as the family-ranking section of 21_transfer_audit.R, with
# the null-seed base swapped to {22,122,222,322,422}; the seed-22 run must match
# the reference family_ranking_v2.tsv bit for bit (regression assertion).
# Inputs (read-only): data/within_study_de_v2.rds, data/pc_precomputed_v2.rds,
#   curation/hox_list.tsv, data/raw/hgnc_complete_set.txt,
#   data/family_ranking_v2.tsv (baseline)
# Outputs: data/family_ranking_seeds_v2.tsv, data/family_seed_summary_v2.tsv 
# =============================================================================
source("scripts/00_setup.R")

run_step("17_family_seeds", function() {
  de <- readRDS(proj_path("data", "within_study_de_v2.rds"))
  pre <- readRDS(proj_path("data", "pc_precomputed_v2.rds"))
  feat <- pre[, c("gene", "contrast", "K", "expr_pct", "detect")]
  hox <- read_tsv("curation/hox_list.tsv")$gene
  fam <- hgnc_families(proj_path("data", "raw", "hgnc_complete_set.txt"))
  fam <- rbind(fam, data.frame(family = "HOX_39", gene = hox, stringsAsFactors = FALSE))
  fam_list <- split(fam$gene, fam$family)
  pre_key <- paste(pre$gene, pre$contrast, sep = "|")
  in_feat <- function(t) t[paste(t$gene, t$contrast) %in% paste(feat$gene, feat$contrast), , drop = FALSE]

  seeds <- c(22L, 122L, 222L, 322L, 422L)
  out <- list()
  for (s in seeds) {
    t0 <- Sys.time()
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
      if (!nrow(tg_in)) return(base)   # guard: skip null draws when the key intersection is empty (not hit on the normal path)
      dr <- sample_matched_null(tg_in, feat, CFG$n_family_null_sets, derive_seed(s, i * CFG$n_family_null_sets))
      v <- null_rates(dr, tg_in, de, pre)
      base$rate <- transfer_rate(o)
      base$null_mean <- mean(v, na.rm = TRUE); base$null_sd <- stats::sd(v, na.rm = TRUE)
      base$z <- (base$rate - base$null_mean) / base$null_sd
      base$empirical_p <- empirical_p(base$rate, v)
      base$included <- TRUE
      base
    }))
    fr$seed <- s
    out[[as.character(s)]] <- fr
    log_msg("seed ", s, " done: ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min; ",
            "included=", sum(fr$included))
  }
  all <- do.call(rbind, out)
  all$rank_z <- NA_integer_
  for (s in seeds) {
    idx <- which(all$seed == s)
    inc <- idx[all$included[idx]]
    ord <- inc[order(-all$z[inc])]
    all$rank_z[ord] <- seq_along(ord)
  }
  write_tsv(all, "data/family_ranking_seeds_v2.tsv")

  # regression assertion: seed=22 matches the reference baseline bit for bit
  base <- read_tsv("data/family_ranking_v2.tsv")
  b22 <- all[all$seed == 22L, ]
  b22 <- b22[match(base$family, b22$family), ]
  same <- isTRUE(all.equal(b22$z, base$z)) && isTRUE(all.equal(b22$rank_z, base$rank_z))
  check(same, "seed=22 did not reproduce the reference baseline (family_ranking_v2.tsv) -- first check that section 5 here matches script 21")

  # summary: HOX_39 rank range across seeds
  hx <- all[all$family == "HOX_39", c("seed", "rate", "null_mean", "null_sd", "z", "empirical_p", "rank_z")]
  n_inc <- vapply(seeds, function(s) sum(all$included[all$seed == s]), 0L)
  hx$n_included <- n_inc
  hx$rank_frac <- hx$rank_z / hx$n_included
  write_tsv(hx, "data/family_seed_summary_v2.tsv")
  log_msg("HOX_39 ranks (5 seeds): ", paste(hx$rank_z, collapse = ","),
          " / ", paste(unique(hx$n_included), collapse = "/"),
          "; z range ", round(min(hx$z), 2), "-", round(max(hx$z), 2))
  invisible(hx)
})
