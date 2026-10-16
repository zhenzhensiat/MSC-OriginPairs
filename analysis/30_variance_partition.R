# =============================================================================
# 30_variance_partition.R —— variance partitioning (exploratory): how much of the
# expression variance is explained by tissue of origin vs study, computed manually
# (partial R^2 on rank-transformed data; the variancePartition package was not
# compatible with the stored objects).
# Output: data/variance_partition_v2.tsv (regenerated).
# =============================================================================
source("scripts/00_setup.R")
need_pkgs("limma")
run_step("30_variance_partition", function() {
  studies <- load_studies("D")
  multi <- Filter(function(s) {
    d <- unique(s$meta$source[s$meta$arm=="include" & s$meta$source %in% CFG$classes])
    length(d) >= 2
  }, studies)
  check(length(multi) >= 10, "multi-source < 10: ", length(multi))
  log_msg("VP: ", length(multi), " multi-source D studies")
  all_genes <- Reduce(intersect, lapply(multi, function(s) rownames(s$expr)))
  check(length(all_genes) >= 2000, "gene overlap: ", length(all_genes))
  parts <- lapply(multi, function(s) {
    m <- s$meta[s$meta$arm=="include" & s$meta$source %in% CFG$classes,,drop=FALSE]
    x <- s$expr[all_genes, m$gsm, drop=FALSE]
    x <- apply(x, 2, function(v) (rank(v, ties.method="average") - 0.5)/length(v))
    meta <- data.frame(gsm=m$gsm, source=m$source, study=s$gse, cluster_id=s$cluster_id,
                       donor=m$donor, passage=m$passage, medium=m$medium,
                       stringsAsFactors=FALSE, row.names=m$gsm)
    list(x=x, meta=meta)
  })
  expr <- do.call(cbind, lapply(parts, `[[`, "x"))
  meta <- do.call(rbind, lapply(parts, `[[`, "meta"))
  check(identical(colnames(expr), meta$gsm), "col mismatch")
  meta$passage[is.na(meta$passage)] <- mean(meta$passage, na.rm=TRUE)
  meta$medium[is.na(meta$medium) | !nzchar(meta$medium)] <- "unknown"
  log_msg("Matrix: ", nrow(expr), " genes x ", ncol(expr), " samples (", length(unique(meta$study)), " studies)")

  # ---- Manual partial R-squared (sequential SS, no lme4 dependency) ----
  # Full model: expr ~ source + study + passage + medium
  # Reduced models: remove each term, compute SS drop / total SS
  design_full <- model.matrix(~ source + study + passage + medium, data=meta)
  # Clean singular columns
  qr_obj <- qr(design_full)
  rank_full <- qr_obj$rank
  if (rank_full < ncol(design_full)) {
    keep_cols <- qr_obj$pivot[seq_len(rank_full)]
    design_full <- design_full[, keep_cols, drop=FALSE]
    log_msg("Design matrix rank deficient: kept ", rank_full, "/", ncol(design_full), " columns")
  }
  # Fit full model with limma (matrix response)
  fit_full <- limma::lmFit(expr, design_full)
  ss_full <- rowSums((expr - fit_full$coefficients %*% t(design_full) )^2, na.rm=TRUE)
  ss_total <- rowSums((expr - rowMeans(expr, na.rm=TRUE))^2, na.rm=TRUE)

  # Partial R-sq for source: remove source from model
  design_no_src <- model.matrix(~ study + passage + medium, data=meta)
  qr_ns <- qr(design_no_src); rk_ns <- qr_ns$rank
  if (rk_ns < ncol(design_no_src)) design_no_src <- design_no_src[, qr_ns$pivot[seq_len(rk_ns)], drop=FALSE]
  fit_ns <- limma::lmFit(expr, design_no_src)
  ss_ns <- rowSums((expr - fit_ns$coefficients %*% t(design_no_src))^2, na.rm=TRUE)

  # Partial R-sq for study: remove study
  design_no_study <- model.matrix(~ source + passage + medium, data=meta)
  qr_nst <- qr(design_no_study); rk_nst <- qr_nst$rank
  if (rk_nst < ncol(design_no_study)) design_no_study <- design_no_study[, qr_nst$pivot[seq_len(rk_nst)], drop=FALSE]
  fit_nst <- limma::lmFit(expr, design_no_study)
  ss_nst <- rowSums((expr - fit_nst$coefficients %*% t(design_no_study))^2, na.rm=TRUE)

  # Assemble
  res <- data.frame(
    gene = rownames(expr),
    R2_source = pmax(0, (ss_ns - ss_full) / pmax(ss_total, 1e-15)),
    R2_study  = pmax(0, (ss_nst - ss_full) / pmax(ss_total, 1e-15)),
    R2_residual = pmax(0, ss_full / pmax(ss_total, 1e-15)),
    stringsAsFactors=FALSE, row.names=NULL)
  # Clip to [0,1] and normalize to sum <= 1
  res$R2_source <- pmin(res$R2_source, 1)
  res$R2_study <- pmin(res$R2_study, 1 - res$R2_source)
  res$R2_residual <- pmin(res$R2_residual, 1 - res$R2_source - res$R2_study)

  hox <- read_tsv("curation/hox_list.tsv")$gene
  mk <- unique(na.omit(read_tsv("curation/marker_list.tsv")$hgnc_symbol))
  res$group <- ifelse(res$gene %in% hox, "HOX", ifelse(res$gene %in% mk, "published", "all"))
  res$analysis <- "exploratory_partial_R2"
  res$collinearity_note <- "single-source studies excluded; sequential SS (source first)"
  write_tsv(res, "data/variance_partition_v2.tsv")
  for (g in c("all","HOX","published")) {
    s <- res[res$group==g,,drop=FALSE]
    if (nrow(s)==0) next
    cat(sprintf("%-10s n=%-5d source=%.4f study=%.4f residual=%.4f\n", g, nrow(s),
      stats::median(s$R2_source, na.rm=TRUE), stats::median(s$R2_study, na.rm=TRUE),
      stats::median(s$R2_residual, na.rm=TRUE)))
  }
  log_msg("Phase 3 done (exploratory, manual partial R-sq)")
  invisible(res)
})
