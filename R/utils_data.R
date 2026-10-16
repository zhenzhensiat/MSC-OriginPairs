# =============================================================================
# utils_data.R —— study objects: reading, validation, training views (conventions in README)
#
# Each study is stored as data/studies/<GSE>.rds, a list with:
#   gse        character, e.g. "GSE199826"
#   cluster_id study-group code (single-source / same-group series share a code)
#   layer      "D" development / "C" validation / "supp" used in an earlier phase / "O" other
#   stability  TRUE/FALSE: used in passage-stability analyses (independent of layer)
#   scale      "logCPM" | "log2TPM" | "array"
#   expr       numeric matrix, rows = HGNC symbols (unique, no NA), cols = GSM; log scale
#   counts     optional raw counts matrix (required for RNA-seq, used by voom)
#   meta       data.frame: gsm, source(AD/BM/UC/other), donor, passage (numeric or NA),
#              medium, arm ("include"/other), sex_check, site
# =============================================================================

validate_study <- function(s) {
  need <- c("gse", "cluster_id", "layer", "scale", "expr", "meta")
  check(all(need %in% names(s)), "study object is missing fields: ", paste(setdiff(need, names(s)), collapse = ","))
  check(s$scale %in% c("logCPM", "log2TPM", "array"), s$gse, ": invalid scale")
  check(s$layer %in% c("D", "C", "supp", "O"), s$gse, ": invalid layer")
  check(is.logical(s$stability) && length(s$stability) == 1, s$gse, ": stability must be TRUE/FALSE")
  x <- s$expr
  check(is.matrix(x) && is.numeric(x), s$gse, ": expr is not a numeric matrix")
  check(!anyNA(rownames(x)) && !anyDuplicated(rownames(x)), s$gse, ": gene names missing or duplicated (merge probes first)")
  check(nrow(x) <= 45000, s$gse, ": gene rows ", nrow(x), " exceed 45000; HGNC mapping (M2) likely not applied")
  check(!anyNA(x), s$gse, ": expr contains NA (imputation not allowed; drop these rows and log them first)")
  check(max(x) < CFG$log_scale_max, s$gse, ": max value ", round(max(x), 1), "; likely not log-transformed")
  m <- s$meta
  check(all(c("gsm", "source", "donor", "passage", "medium", "arm") %in% names(m)), s$gse, ": meta is missing columns")
  check(identical(colnames(x), m$gsm), s$gse, ": expr column names do not match meta$gsm (order must match too)")
  check(!anyDuplicated(m$gsm), s$gse, ": duplicated GSM")
  if (!is.null(s$counts)) check(identical(colnames(s$counts), m$gsm), s$gse, ": counts column names do not match meta")
  invisible(TRUE)
}

# layer selects layers; stability_only = TRUE keeps only stability data (layer ignored)
load_studies <- function(layer = NULL, stability_only = FALSE, dir = proj_path("data", "studies")) {
  files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
  check(length(files) > 0, "no study files under data/studies")
  st <- lapply(files, readRDS)
  if (stability_only) st <- Filter(function(s) isTRUE(s$stability), st)
  if (!is.null(layer)) st <- Filter(function(s) s$layer %in% layer, st)
  check(length(st) > 0, "no studies match the criteria (layer = ", paste(layer, collapse = "/"),
        ", stability_only = ", stability_only, ")")
  for (s in st) validate_study(s)
  st <- lapply(st, function(s) {
    s$meta$study <- s$gse
    s$meta$cluster_id <- s$cluster_id
    s
  })
  names(st) <- vapply(st, `[[`, "", "gse")
  log_msg("loaded ", length(st), " studies (layer ", paste(layer, collapse = "/"), ")")
  st
}

STD_META_COLS <- c("gsm", "study", "cluster_id", "donor", "source", "passage", "medium", "arm")

bind_meta <- function(studies, include_only = TRUE) {
  m <- do.call(rbind, lapply(studies, function(s) s$meta[, STD_META_COLS]))
  if (include_only) m <- m[m$arm == "include" & m$source %in% CFG$classes, , drop = FALSE]
  check(!anyDuplicated(m$gsm), "cross-study duplicated GSM: independence screening not clean")
  rownames(m) <- NULL
  m
}

# training view: exclude held-out clusters; genome-wide matrices only; genes = intersection across training studies (all feature sources share one universe)
# full_genome criterion (2026-10-10): when a study has a full_genome field, that field is authoritative
# (author-prescreened matrices such as GSE38947 can exceed 10000 genes but are not genome-wide); fall back to the gene-count rule when absent
# per-fold universe = training intersection ∩ gene-name intersection of all held-out cluster studies
# (gene names only; expression values/labels are not read, so there is no leakage; matches real deployment, where users know in advance which genes a platform measures)
test_genes <- function(studies, cluster) {
  te <- Filter(function(s) identical(s$cluster_id, cluster), studies)
  check(length(te) >= 1, "held-out cluster ", cluster, " has no studies")
  Reduce(intersect, lapply(te, function(s) rownames(s$expr)))
}

training_view <- function(studies, exclude_cluster, min_genes = CFG$genome_min_genes,
                          restrict_genes = NULL) {
  tr <- Filter(function(s) !identical(s$cluster_id, exclude_cluster) &&
                 (if (is.null(s$full_genome)) nrow(s$expr) >= min_genes else isTRUE(s$full_genome)),
               studies)
  check(length(tr) >= 3, "fewer than 3 training studies")
  genes <- Reduce(intersect, lapply(tr, function(s) rownames(s$expr)))
  if (!is.null(restrict_genes)) genes <- intersect(genes, restrict_genes)
  check(length(genes) >= 2000, "training gene universe has only ", length(genes), " genes")
  parts <- lapply(tr, function(s) {
    m <- s$meta[s$meta$arm == "include" & s$meta$source %in% CFG$classes, STD_META_COLS, drop = FALSE]
    list(x = s$expr[genes, m$gsm, drop = FALSE], m = m)
  })
  list(expr = do.call(cbind, lapply(parts, `[[`, "x")),
       meta = do.call(rbind, lapply(parts, `[[`, "m")),
       studies = tr, genes = genes)
}

test_view <- function(studies, cluster) {
  te <- Filter(function(s) identical(s$cluster_id, cluster), studies)
  check(length(te) >= 1, "held-out cluster ", cluster, " has no studies")
  lapply(te, function(s) {
    m <- s$meta[s$meta$arm == "include" & s$meta$source %in% CFG$classes, STD_META_COLS, drop = FALSE]
    list(expr = s$expr[, m$gsm, drop = FALSE], meta = m)
  })
}

# gene features within the training fold (training data only, leakage-safe), for random anchor matching
gene_features_train <- function(tv) {
  per <- lapply(tv$studies, function(s) {
    x <- s$expr[tv$genes, s$meta$gsm[s$meta$arm == "include" & s$meta$source %in% CFG$classes], drop = FALSE]
    thr <- if (s$scale == "array") stats::quantile(x, 0.25) else 0
    ave <- rowMeans(x)
    cbind(pct = (rank(ave) - 0.5) / length(ave), det = rowMeans(x > thr))
  })
  data.frame(gene = tv$genes, contrast = "all", K = 0L,
             expr_pct = rowMeans(sapply(per, function(p) p[, "pct"])),
             detect = rowMeans(sapply(per, function(p) p[, "det"])),
             stringsAsFactors = FALSE)
}

# ---- series matrix header read-only (for layer C sample-level checks; expression values are never read)----
# local parsing fix: field names split with strsplit, first token, "!" stripped (the original sub capture group returned empty in this environment)
read_series_header <- function(path) {
  con <- gzfile(path, "r"); on.exit(close(con))
  out <- character(0)
  repeat {
    ln <- readLines(con, n = 500, warn = FALSE)
    if (!length(ln)) break
    sa <- grep("^!series_matrix_table_begin", ln)
    if (length(sa)) { out <- c(out, ln[seq_len(sa[1] - 1)]); break }
    out <- c(out, ln)
  }
  smp <- grep("^!Sample_", out, value = TRUE)
  toks <- strsplit(smp, "	", fixed = TRUE)
  fields <- sub("^!", "", vapply(toks, `[`, "", 1))
  vals <- lapply(toks, function(v) gsub('^"|"$', "", v[-1]))
  n <- max(lengths(vals))
  df <- as.data.frame(do.call(cbind, lapply(vals, function(v) c(v, rep(NA, n - length(v))))),
                      stringsAsFactors = FALSE)
  names(df) <- make.unique(fields)
  df
}
