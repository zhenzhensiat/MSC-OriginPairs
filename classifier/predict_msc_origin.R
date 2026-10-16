# =============================================================================
# predict_msc_origin.R —— MSC-OriginPairs: tissue-of-origin prediction for
# human MSC samples (AD / BM / UC).
#
# Usage:  Rscript predict_msc_origin.R <input.tsv> [model_dir]
#         (model_dir defaults to the folder containing this script)
#
# Input:  a tab-separated file whose first column is named `gene` (HGNC symbols)
#         followed by one column per sample. Values must only preserve the
#         within-sample ordering of expression (higher = higher expression):
#         log-scale values or within-sample ranks. For qPCR data use negative
#         cycle thresholds (−Ct); raw Ct values must not be used directly,
#         because a lower Ct means higher expression and the ordering would be
#         inverted. The bundled example uses the pipeline's log-scale values. Genes that your platform did not measure may be
#         omitted or left empty; samples are scored one at a time, so missing
#         genes in one sample do not affect the others and rules containing a
#         missing gene are skipped automatically.
#
# Output: stdout TSV — sample, call, score_AD, score_BM, score_UC, margin,
#         low_confidence. The margin is the difference between the top two
#         normalised scores; calls with a margin below the low-confidence threshold
#         (model_v2_tau.rds) are flagged low-confidence and should be checked
#         against other information before being relied upon.
#
# Dependencies: R (>= 4.1), multiclassPairs, readr.
# Model: model_v2.rds (rules + weights) and model_v2_tau.rds (threshold) are
#        the final versions, used unchanged for external validation.
# =============================================================================

suppressPackageStartupMessages({
  library(readr)
  library(multiclassPairs)
})

CLASSES <- c("AD", "BM", "UC")

predict_msc_origin <- function(input_path, model_dir = ".") {
  clf <- readRDS(file.path(model_dir, "model_v2.rds"))
  tau <- readRDS(file.path(model_dir, "model_v2_tau.rds"))

  dt <- readr::read_tsv(input_path, show_col_types = FALSE)
  if (ncol(dt) < 2) stop("input needs a gene column plus at least one sample column")
  genes <- as.character(dt[[1]])
  if (anyDuplicated(genes)) stop("duplicated gene symbols in the gene column")
  expr <- as.matrix(dt[, -1, drop = FALSE])
  if (!is.numeric(expr)) stop("sample columns must be numeric")

  out <- vector("list", ncol(expr))
  for (j in seq_len(ncol(expr))) {
    vals <- expr[, j]
    keep <- !is.na(vals)
    sample_id <- colnames(expr)[j]
    if (!any(keep)) {
      message("sample ", sample_id, ": no usable values, skipped")
      out[[j]] <- data.frame(sample = sample_id, call = NA_character_,
                             score_AD = NA_real_, score_BM = NA_real_, score_UC = NA_real_,
                             margin = NA_real_, low_confidence = NA, stringsAsFactors = FALSE)
      next
    }
    sub <- matrix(vals[keep], ncol = 1, dimnames = list(genes[keep], sample_id))
    res <- multiclassPairs::predict_one_vs_rest_TSP(
      classifier = clf, Data = sub, tolerate_missed_genes = TRUE,
      weighted_votes = TRUE, classes = CLASSES, verbose = FALSE)
    if (!all(c("max_score", CLASSES) %in% colnames(res)))
      stop("prediction output columns differ from expectation; check the multiclassPairs version")
    s <- as.numeric(res[1, CLASSES])
    names(s) <- CLASSES
    sn <- s / suppressWarnings(max(sum(s, na.rm = TRUE), .Machine$double.eps))
    top <- sort(sn, decreasing = TRUE)
    margin <- top[1] - top[2]
    out[[j]] <- data.frame(sample = sample_id, call = as.character(res$max_score[1]),
                           score_AD = s["AD"], score_BM = s["BM"], score_UC = s["UC"],
                           margin = margin, low_confidence = margin < tau,
                           stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

# ---- command-line entry point (not run when sourced) ------------------------
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1) stop("usage: Rscript predict_msc_origin.R <input.tsv> [model_dir]")
  model_dir <- if (length(args) >= 2) args[2] else {
    fa <- commandArgs(trailingOnly = FALSE)
    f <- sub("^--file=", "", fa[grep("^--file=", fa)])
    if (length(f)) dirname(normalizePath(f)) else "."
  }
  options(digits = 15)
  out <- predict_msc_origin(args[1], model_dir)
  utils::write.table(out, stdout(), sep = "\t", quote = FALSE, row.names = FALSE)
}
