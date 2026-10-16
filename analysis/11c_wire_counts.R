# =============================================================================
# 11c_wire_counts.R —— wire the counts matrix for GSE283154 (development layer;
# used by voom in Phase 2).
# Mapping: series-matrix sample titles <-> counts column names, keyed on the
# "sample_donor" prefix; only include-arm samples (CT control); the three source
# files (BM/UC/hASC) are merged and aligned by GSM. Supplementary-layer counts
# (GSE291956/GSE275980) are registered as deferred.
# Outputs: updated data/studies/GSE283154.rds; data/counts_wiring_v2.tsv
# =============================================================================
source("scripts/00_setup.R")

run_step("11c_wire_counts", function() {
  g <- "GSE283154"
  s <- readRDS(proj_path("data", "studies", paste0(g, ".rds")))
  check(is.null(s$counts), g, ": counts already wired, not repeating")

  # series matrix: title <-> gsm
  sm_lines <- readLines(gzfile(proj_path("data", "raw", g,
    "GSE283154-GPL24676_series_matrix.txt.gz")))
  get_row <- function(tag) {
    ln <- sm_lines[startsWith(sm_lines, tag)]
    check(length(ln) == 1, g, ": series matrix missing row ", tag)
    trimws(gsub('"', "", strsplit(ln, "\t", fixed = TRUE)[[1]][-1]))
  }
  ttl <- get_row("!Sample_title"); gsm <- get_row("!Sample_geo_accession")
  check(length(ttl) == length(gsm) && !anyDuplicated(gsm), g, ": title/gsm structure abnormal")
  title_key <- vapply(ttl, function(t) paste(strsplit(t, "_", fixed = TRUE)[[1]][1:2],
                                              collapse = "_"), character(1))

  # merge the three counts files; key = column name without arm suffix
  files <- c("GSE283154_raw_counts_BM-MSC_samples.txt.gz",
             "GSE283154_raw_counts_UC-MSC_samples.txt.gz",
             "GSE283154_raw_counts_hASC_samples.txt.gz")
  mats <- lapply(files, function(f) {
    raw <- read.delim(gzfile(proj_path("data", "raw", g, f)), check.names = FALSE,
                      stringsAsFactors = FALSE)
    m <- as.matrix(raw[, -1, drop = FALSE]); storage.mode(m) <- "numeric"
    rownames(m) <- raw[[1]]
    m
  })
  num <- do.call(cbind, mats)
  check(!anyDuplicated(colnames(num)), g, ": duplicate counts column names after merge")
  col_key <- sub("-(CT|Ago)$", "", colnames(num))
  keep <- col_key %in% title_key & grepl("-CT$", colnames(num))   # control arm only
  num <- num[, keep, drop = FALSE]
  colnames(num) <- gsm[match(col_key[keep], title_key)]
  if (anyDuplicated(rownames(num))) num <- rowsum(num, group = rownames(num))

  common <- intersect(colnames(num), s$meta$gsm[s$meta$arm == "include"])
  check(length(common) == sum(s$meta$arm == "include"),
        g, ": include-sample counts misaligned (", length(common), "/",
        sum(s$meta$arm == "include"), ")")
  s$counts <- num[, common, drop = FALSE]
  validate_study(s)
  saveRDS(s, proj_path("data", "studies", paste0(g, ".rds")))

  reg <- data.frame(
    gse = c(g, "GSE291956", "GSE275980"),
    status = c("wired", "deferred_supp", "deferred_supp"),
    n_genes = c(nrow(num), NA, NA),
    n_matched = c(length(common), NA, NA),
    note = c("CT arm, 12 samples, title-prefix mapping", "supp layer: no voom in Phase 2", "same as above"),
    stringsAsFactors = FALSE)
  write_tsv(reg, "data/counts_wiring_v2.tsv")
  log_msg(g, ": counts wired ", nrow(num), " genes x ", length(common), " samples (CT control arm only)")
  invisible(reg)
})
