# =============================================================================
# utils_stats.R —— within-study differential analysis, partial-conjunction tests,
# matched nulls and gene-family ranking (analysis protocol, Phase 2).
#
# Two key design points:
#   1) Units without a claimed direction (HOX, families, their matched nulls):
#      two-sided by "smaller side x 2" (partial conjunction in each direction,
#      the smaller p doubled); published markers with a claimed direction:
#      one-sided in the claimed direction.
#   2) BH correction is applied inside each testing group (observed and each
#      matched-null set follow exactly the same procedure).
# =============================================================================

# ---- direction and contrast mapping -------------------------------------------------------------
canonical_contrast <- function(high, low, contrasts = CFG$contrasts) {
  fwd <- paste(high, low, sep = "-"); rev <- paste(low, high, sep = "-")
  if (fwd %in% contrasts) return(list(contrast = fwd, direction = 1L))
  if (rev %in% contrasts) return(list(contrast = rev, direction = -1L))
  stop("cannot map to a canonical contrast: ", high, " > ", low)
}

# marker table -> one "gene x contrast" test per row
# input columns: claim_id, gene, high_source, low_source (semicolon-separated, only AD/BM/UC allowed), discovery_gse (semicolon-separated, may be empty)
expand_claims <- function(markers) {
  check(all(c("claim_id", "gene", "high_source", "low_source", "discovery_gse") %in% names(markers)),
        "marker table is missing columns; required: claim_id, gene, high_source, low_source, discovery_gse")
  if ("v2_flag" %in% names(markers)) {
    markers <- markers[is.na(markers$v2_flag) | markers$v2_flag != "no_class_projection", , drop = FALSE]
  }
  rows <- lapply(seq_len(nrow(markers)), function(i) {
    m <- markers[i, ]
    highs <- trimws(strsplit(m$high_source, ";", fixed = TRUE)[[1]])
    lows  <- trimws(strsplit(m$low_source,  ";", fixed = TRUE)[[1]])
    check(!anyNA(c(highs, lows)) && length(highs) > 0 && length(lows) > 0 &&
          all(c(highs, lows) %in% CFG$classes) && length(intersect(highs, lows)) == 0,
          "marker ", m$claim_id, " source is not AD/BM/UC or high/low overlap; standardize manually before use")
    grid <- expand.grid(hi = highs, lo = lows, stringsAsFactors = FALSE)
    do.call(rbind, lapply(seq_len(nrow(grid)), function(k) {
      cc <- canonical_contrast(grid$hi[k], grid$lo[k])
      data.frame(claim_id = m$claim_id, gene = m$gene, contrast = cc$contrast,
                 direction = cc$direction,
                 exclude = ifelse(is.na(m$discovery_gse), "", m$discovery_gse),
                 stringsAsFactors = FALSE)
    }))
  })
  do.call(rbind, rows)
}

# ---- within-study differential analysis (limma; voom for counts)-------------------------------------
# study object conventions: see README "data conventions"
within_study_de <- function(s, contrasts = CFG$contrasts, min_per_class = 2L) {
  m <- s$meta[s$meta$arm == "include" & s$meta$source %in% CFG$classes, , drop = FALSE]
  tab <- table(m$source)
  cls <- intersect(CFG$classes, names(tab)[tab >= min_per_class])
  m <- m[m$source %in% cls, , drop = FALSE]
  avail <- contrasts[vapply(strsplit(contrasts, "-", fixed = TRUE), function(p) all(p %in% cls), logical(1))]
  if (!length(avail)) return(NULL)

  grp <- factor(m$source, levels = cls)
  design <- stats::model.matrix(~ 0 + grp)
  colnames(design) <- levels(grp)

  use_counts <- !is.null(s$counts)
  if (use_counts) {
    y <- edgeR::DGEList(s$counts[, m$gsm, drop = FALSE])
    keep <- edgeR::filterByExpr(y, design)
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- edgeR::calcNormFactors(y)
    x <- limma::voom(y, design)
    e <- x$E
  } else {
    x <- s$expr[, m$gsm, drop = FALSE]
    e <- x
  }

  block <- NULL; cor <- NULL
  if (anyDuplicated(m$donor)) {   # multiple samples per donor: use donor as block
    dc <- limma::duplicateCorrelation(x, design, block = m$donor)
    cor <- dc$consensus.correlation; block <- m$donor
  }
  fit <- limma::lmFit(x, design, block = block, correlation = cor)
  cm <- limma::makeContrasts(contrasts = avail, levels = design)
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, cm), trend = (!use_counts && s$scale != "array"))

  det_thr <- if (s$scale == "array") stats::quantile(e, 0.25, na.rm = TRUE) else 0
  detect <- rowMeans(e > det_thr)
  ave <- rowMeans(e)
  expr_pct <- (rank(ave) - 0.5) / length(ave)

  do.call(rbind, lapply(avail, function(cn) {
    tt <- limma::topTable(fit2, coef = cn, number = Inf, sort.by = "none")
    g <- rownames(tt)
    check(identical(g, rownames(fit2$coefficients)), "topTable row order does not match the fitted object")
    data.frame(study = s$gse, contrast = cn, gene = g,
               logFC = tt$logFC, t = tt$t, df = fit2$df.total,
               p = tt$P.Value, q_within = tt$adj.P.Val,
               expr_pct = expr_pct[g], detect = detect[g],
               stringsAsFactors = FALSE, row.names = NULL)
  }))
}

# ---- partial conjunction (Benjamini & Heller 2008; Fisher combines the largest K-r+1 p-values)---------
# must be verified against the original formula; run the H0 calibration tests in tests
pc_fisher <- function(p, r) {
  K <- length(p)
  if (K == 0 || r > K || r < 1) return(NA_real_)
  p <- sort(pmax(p, .Machine$double.xmin))
  tail_p <- p[r:K]
  stats::pchisq(-2 * sum(log(tail_p)), df = 2 * length(tail_p), lower.tail = FALSE)
}

r_majority <- function(K) ceiling(K / 2)

pc_directional <- function(t, df, direction, r_rule = r_majority) {
  pc_fisher(stats::pt(direction * t, df, lower.tail = FALSE), r_rule(length(t)))
}

# precompute: genome-wide partial-conjunction p for each "gene x contrast" in both directions (no excluded studies)
precompute_pc <- function(de, min_K = CFG$min_K) {
  key <- paste(de$gene, de$contrast, sep = "|")
  idx <- split(seq_len(nrow(de)), key)
  res <- vapply(idx, function(r) {
    K <- length(r)
    feat <- c(mean(de$expr_pct[r]), mean(de$detect[r]))   # matching features are computed for all rows
    if (K < min_K) return(c(K, NA_real_, NA_real_, feat))
    c(K, pc_directional(de$t[r], de$df[r], 1L), pc_directional(de$t[r], de$df[r], -1L), feat)
  }, numeric(5))
  kk <- strsplit(names(idx), "|", fixed = TRUE)
  data.frame(gene = vapply(kk, `[`, "", 1), contrast = vapply(kk, `[`, "", 2),
             K = res[1, ], p_plus = res[2, ], p_minus = res[3, ],
             expr_pct = res[4, ], detect = res[5, ],
             stringsAsFactors = FALSE, row.names = NULL)
}

# compute p for a set of test rows (direction: +1/-1 = one-sided claimed direction; 0 = min of both directions x2)
# recompute on the fly when studies are excluded (exclude non-empty); otherwise use the precomputed table
pc_rows <- function(targets, de, pre, min_K = CFG$min_K) {
  check(all(c("gene", "contrast", "direction", "exclude") %in% names(targets)),
        "targets must contain gene, contrast, direction, exclude")
  targets$exclude[is.na(targets$exclude)] <- ""
  pre_key <- paste(pre$gene, pre$contrast, sep = "|")
  de_idx <- if (any(nzchar(targets$exclude))) {
    split(seq_len(nrow(de)), paste(de$gene, de$contrast, sep = "|"))
  } else NULL
  out <- lapply(seq_len(nrow(targets)), function(i) {
    tg <- targets[i, ]
    ex <- tg$exclude
    if (!nzchar(ex)) {
      j <- match(paste(tg$gene, tg$contrast, sep = "|"), pre_key)
      if (is.na(j) || pre$K[j] < min_K) return(c(if (is.na(j)) 0 else pre$K[j], NA, NA))
      pp <- pre$p_plus[j]; pm <- pre$p_minus[j]; K <- pre$K[j]
    } else {
      r <- de_idx[[paste(tg$gene, tg$contrast, sep = "|")]]
      r <- r[!de$study[r] %in% strsplit(ex, ";", fixed = TRUE)[[1]]]
      K <- length(r)
      if (K < min_K) return(c(K, NA, NA))
      pp <- pc_directional(de$t[r], de$df[r], 1L); pm <- pc_directional(de$t[r], de$df[r], -1L)
    }
    if (tg$direction == 0L) {
      d <- if (pp <= pm) 1 else -1
      c(K, min(1, 2 * min(pp, pm)), d)
    } else c(K, if (tg$direction > 0) pp else pm, tg$direction)
  })
  m <- do.call(rbind, out)
  data.frame(targets, K = m[, 1], pc_p = m[, 2], dir_used = m[, 3], stringsAsFactors = FALSE, row.names = NULL)
}

# fast path: no excluded studies and both-directions min x2: vectorized table lookup (for family ranking and HOX null sets)
pc_either_fast <- function(gene, contrast, pre, pre_key = paste(pre$gene, pre$contrast, sep = "|"),
                           min_K = CFG$min_K) {
  j <- match(paste(gene, contrast, sep = "|"), pre_key)
  K <- ifelse(is.na(j), 0, pre$K[j])
  pp <- pre$p_plus[j]; pm <- pre$p_minus[j]
  p <- ifelse(is.na(j) | K < min_K, NA_real_, pmin(1, 2 * pmin(pp, pm)))
  data.frame(gene = gene, contrast = contrast, K = K, pc_p = p,
             dir_used = ifelse(is.na(p), NA_real_, ifelse(pp <= pm, 1, -1)),
             stringsAsFactors = FALSE)
}

# within-group BH and migration calls; migration rate = migrated / testable
family_calls <- function(tbl, q = CFG$pc_q) {
  tbl$testable <- !is.na(tbl$pc_p)
  tbl$pc_q <- NA_real_
  tbl$pc_q[tbl$testable] <- stats::p.adjust(tbl$pc_p[tbl$testable], "BH")
  tbl$transfer <- tbl$testable & !is.na(tbl$pc_q) & tbl$pc_q < q
  tbl
}

transfer_rate <- function(tbl) {
  n <- sum(tbl$testable)
  if (n == 0) return(NA_real_)
  sum(tbl$transfer) / n
}

# ---- matched null (R11)-----------------------------------------------------------
qbin <- function(x, n) {
  br <- unique(stats::quantile(x, probs = seq(0, 1, length.out = n + 1), na.rm = TRUE))
  if (length(br) < 2) return(rep(1L, length(x)))
  as.integer(cut(x, br, include.lowest = TRUE))
}

# feat: gene, contrast, K, expr_pct, detect (from precompute_pc or training-fold features)
add_strata <- function(feat) {
  feat$expr_bin <- stats::ave(feat$expr_pct, feat$contrast, FUN = function(x) qbin(x, CFG$expr_bins))
  feat$det_bin  <- stats::ave(feat$detect,   feat$contrast, FUN = function(x) qbin(x, CFG$detect_bins))
  feat$L0 <- paste(feat$contrast, feat$K, feat$expr_bin, feat$det_bin, sep = "|")
  feat$L1 <- paste(feat$contrast, feat$K, feat$expr_bin, sep = "|")
  feat$L2 <- paste(feat$contrast, feat$expr_bin, sep = "|")
  feat
}

# generate n_sets matched null sets; generate and persist all of them first, then compute
# targets must contain gene, contrast; feat is the genome-wide feature table (including the target genes themselves)
sample_matched_null <- function(targets, feat, n_sets, seed_base) {
  feat <- add_strata(feat)
  key_f <- paste(feat$gene, feat$contrast, sep = "|")
  ti <- match(paste(targets$gene, targets$contrast, sep = "|"), key_f)
  check(!anyNA(ti), "target genes are missing from the feature table: ",
        paste(utils::head(targets$gene[is.na(ti)], 5), collapse = ","))
  tf <- feat[ti, , drop = FALSE]
  cand_feat <- feat[!feat$gene %in% unique(targets$gene), , drop = FALSE]
  pools <- vector("list", nrow(targets)); level <- integer(nrow(targets))
  for (i in seq_len(nrow(targets))) {
    same_c <- cand_feat$contrast == targets$contrast[i]
    for (lv in c("L0", "L1", "L2")) {
      cand <- cand_feat$gene[same_c & cand_feat[[lv]] == tf[[lv]][i]]
      if (length(cand) >= CFG$min_pool) { pools[[i]] <- cand; level[i] <- as.integer(substr(lv, 2, 2)); break }
    }
    if (is.null(pools[[i]])) { pools[[i]] <- cand_feat$gene[same_c]; level[i] <- 3L }
  }
  sets <- lapply(seq_len(n_sets), function(s) {
    set.seed(seed_base + s)
    data.frame(set = s, target_row = seq_len(nrow(targets)),
               target_gene = targets$gene, contrast = targets$contrast,
               gene = vapply(pools, function(p) p[sample.int(length(p), 1L)], character(1)),
               relax_level = level, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, sets)
  if (n_sets >= 2) {
    check(!identical(out$gene[out$set == 1], out$gene[out$set == 2]),
          "the first two null sets are identical: RNG was reset, check the seed logic")
  }
  out
}

# compute the migration rate of each null set with the same pipeline as observed rows (same direction rules, same exclusions, within-group BH)
# align by target_row (the same gene and contrast may have claims in opposite directions; do not merge by gene name)
null_rates <- function(null_long, targets, de, pre) {
  check(all(null_long$target_row %in% seq_len(nrow(targets))), "null set target_row out of range")
  nl <- null_long
  nl$direction <- targets$direction[nl$target_row]
  nl$exclude <- if ("exclude" %in% names(targets)) targets$exclude[nl$target_row] else ""
  nl$exclude[is.na(nl$exclude)] <- ""
  if (all(nl$direction == 0L) && !any(nzchar(nl$exclude))) {   # fast path
    pre_key <- paste(pre$gene, pre$contrast, sep = "|")
    sp <- split(seq_len(nrow(nl)), nl$set)
    return(vapply(sp, function(r) {
      tb <- pc_either_fast(nl$gene[r], nl$contrast[r], pre, pre_key)
      transfer_rate(family_calls(tb))
    }, numeric(1)))
  }
  sp <- split(nl, nl$set)
  vapply(sp, function(d) {
    tb <- pc_rows(d[, c("gene", "contrast", "direction", "exclude")], de, pre)
    transfer_rate(family_calls(tb))
  }, numeric(1))
}

# ---- gene families -------------------------------------------------------------------
# HGNC complete set (hgnc_complete_set.txt); uses the gene_group column ('|'-separated)
hgnc_families <- function(hgnc_file, size = CFG$family_size) {
  h <- utils::read.delim(hgnc_file, quote = "", comment.char = "", stringsAsFactors = FALSE)
  check(all(c("symbol", "status", "gene_group") %in% names(h)), "HGNC file is missing columns; check the downloaded version")
  h <- h[h$status == "Approved" & !is.na(h$gene_group) & nzchar(h$gene_group), c("symbol", "gene_group")]
  # gene_group multi-value field carries double quotes (e.g. "A|B"); strip them before splitting;
  # otherwise family names keep quotes at both ends and one family is split into several "same-name" families (family_ranking parsing bug)
  h$gene_group <- gsub('"', "", h$gene_group, fixed = TRUE)
  lst <- strsplit(h$gene_group, "|", fixed = TRUE)
  fam <- data.frame(family = unlist(lst), gene = rep(h$symbol, lengths(lst)), stringsAsFactors = FALSE)
  check(!any(grepl('"', fam$family, fixed = TRUE)), "family names still contain quotes (quote stripping failed)")
  sz <- table(fam$family)
  fam[fam$family %in% names(sz)[sz >= size[1] & sz <= size[2]], , drop = FALSE]
}
