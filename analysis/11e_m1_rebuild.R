# =============================================================================
# 11e_m1_rebuild.R —— probe-level re-annotation rebuild of four studies
# (GSE12274, GSE48022 on GPL570; GSE57151 on GPL16686; GSE19773 on GPL7363).
# Integrity discipline: GPL570 annot.gz checked line by line (>=54,600 data
# rows); platform soft files checked for the end marker and row counts.
# Rules: mapping/deduplication/scale identical to script 11.
# Outputs: data/m1_rebuild_report_v2.tsv (per-gene correlation report);
#   updated expr in the four study objects (expr only); updated coverage audit
#   rows (rewritten). Acceptance: per-gene correlations reported; validate_study
#   passes for all four.
# =============================================================================
source("scripts/00_setup.R")
need_pkgs(c("org.Hs.eg.db", "AnnotationDbi", "hugene20sttranscriptcluster.db", "edgeR"))
library(org.Hs.eg.db)   # select() requires the package object to be attached (same as script 11)

run_step("11e_m1_rebuild", function() {

  # ---- 0. annotation file loading (with integrity checks)------------------------------------------
  ann_dir <- proj_path("data", "raw", "platforms")

  # GPL570.annot.gz: skip ^/!/# metadata; the header row starts with "ID\t"
  load_annot_gz <- function(path, min_rows) {
    con <- gzfile(path, "rt")
    ln <- readLines(con, warn = FALSE)
    close(con)
    ih <- grep("^ID\t", ln)
    check(length(ih) == 1, basename(path), ": unique header row not found")
    dat <- ln[(ih + 1):length(ln)]
    dat <- dat[nzchar(dat)]
    check(length(dat) >= min_rows, basename(path), ": data rows ", length(dat),
          " < expected minimum ", min_rows, " (E7 truncation suspected)")
    df <- read.delim(text = paste(c(ln[ih], dat), collapse = "\n"), header = TRUE,
                     quote = "", check.names = FALSE, stringsAsFactors = FALSE)
    check("Gene symbol" %in% names(df), "missing Gene symbol column")
    m <- setNames(df[["Gene symbol"]], df[["ID"]])
    m <- trimws(sub("\\s*///.*$", "", m))          # multi-value symbols ("A///B", 2,214 rows on this platform): take the first
    m[!nzchar(m) | m == "---"] <- NA_character_
    list(map = m, n_rows = length(dat))
  }

  # GEO CGI platform record (soft.txt): !platform_table_begin ... !platform_table_end
  load_soft_txt <- function(path, sym_col, id_col = "ID", min_rows = 40000) {
    ln <- readLines(path, warn = FALSE)
    ib <- grep("^!platform_table_begin", ln); ie <- grep("^!platform_table_end", ln)
    check(length(ib) == 1 && length(ie) == 1 && ie > ib,
          basename(path), ": platform table markers incomplete (E7 truncation suspected)")
    df <- read.delim(text = paste(ln[ib:ie], collapse = "\n"), header = TRUE,
                     sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE,
                     skip = 1)
    check(nrow(df) >= min_rows, basename(path), ": data rows ", nrow(df), " < ", min_rows)
    check(all(c(id_col, sym_col) %in% names(df)), basename(path), ": missing column ", id_col, "/", sym_col)
    m <- setNames(df[[sym_col]], df[[id_col]])
    m <- trimws(sub("\\s*///.*$", "", m))          # same rule: take the first of multiple values
    m[!nzchar(m) | m == "---"] <- NA_character_
    list(map = m, n_rows = nrow(df))
  }

  maps <- list(
    GPL570  = load_annot_gz(file.path(ann_dir, "GPL570.annot.gz"), 54600),
    GPL7363 = load_soft_txt(file.path(ann_dir, "GPL7363.soft.txt"), sym_col = "Symbol", min_rows = 45000)
  )
  # GPL16686: same source as the v1 route (BioC), key = transcript cluster
  db <- hugene20sttranscriptcluster.db::hugene20sttranscriptcluster.db
  hg <- AnnotationDbi::select(db, keys = AnnotationDbi::keys(db), columns = "SYMBOL")
  hg <- hg[!is.na(hg$SYMBOL), ]
  hg <- hg[!duplicated(hg$PROBEID), ]
  maps$GPL16686 <- list(map = setNames(hg$SYMBOL, hg$PROBEID), n_rows = nrow(hg),
                        sha = paste0("bioc:", as.character(packageVersion("hugene20sttranscriptcluster.db"))))
  for (p in names(maps)) log_msg("annotation source ", p, ": ", maps[[p]]$n_rows, " rows")

  # ---- 1. series matrix loading (assert a unique data table)----------------------------------------
  load_sm <- function(path) {
    ln <- readLines(gzfile(path, "rt"), warn = FALSE)
    ib <- grep("^\"?ID_REF\"?", ln)
    if (length(ib) != 1) ib <- grep("series_matrix_table_begin", ln) + 1
    ie <- grep("series_matrix_table_end", ln)
    check(length(ib) == 1 && length(ie) == 1 && ie > ib, basename(path), ": data table interval abnormal")
    df <- read.delim(text = paste(ln[ib:(ie - 1)], collapse = "\n"), header = TRUE,
                     quote = "\"", check.names = FALSE, stringsAsFactors = FALSE)
    ids <- as.character(df[[1]])   # numeric-only probe IDs (e.g. "16650001" on GPL16686) must not become integers
    m <- as.matrix(df[, -1, drop = FALSE])
    mode(m) <- "numeric"
    na_all <- rowSums(is.na(m)) == ncol(m)         # all-missing probe rows: dropped and recorded (no imputation, E3)
    if (any(na_all)) {
      log_msg(basename(path), ": all-missing probe rows dropped: ", sum(na_all), " rows")
      m <- m[!na_all, , drop = FALSE]; ids <- ids[!na_all]
    }
    check(!anyNA(m), basename(path), ": data contain partial NA/non-numeric values, manual review required")
    list(ids = ids, mat = m)
  }

  collapse_dup <- function(x) {   # same rule as script 11
    if (!anyDuplicated(rownames(x))) return(x)
    mu <- rowMeans(x); iqr <- apply(x, 1, stats::IQR)
    ord <- order(mu, iqr, decreasing = TRUE)
    x[ord[!duplicated(rownames(x)[ord])], , drop = FALSE]
  }

  # ---- 1b. "old matrix" rebuild baseline (from the v1 source, replicating the script-11 pipeline; keeps this script idempotent)--------
  # the first run already replaced expr in the GSE12274/GSE48022 objects with rebuilt values; afterwards reading the objects directly as "old" would
  # contaminate the report into rebuild-vs-rebuild. So the baseline is always rebuilt from the v1 source, and on the studies "not yet replaced"
  # (GSE57151/GSE19773) the baseline is verified to match the original object (fidelity proof).
  geo <- read_tsv("data/census_geo_v2.tsv")
  eg <- AnnotationDbi::select(org.Hs.eg.db, keys = AnnotationDbi::keys(org.Hs.eg.db),
                              columns = c("SYMBOL", "ALIAS"))
  sym_ok <- unique(na.omit(eg$SYMBOL))
  al <- eg[!is.na(eg$ALIAS) & !is.na(eg$SYMBOL), c("ALIAS", "SYMBOL")]
  al_n <- table(al$ALIAS); uniq_alias <- names(al_n[al_n == 1L])
  alias_map <- setNames(al$SYMBOL[match(uniq_alias, al$ALIAS)], uniq_alias)
  map_genes <- function(genes) {
    out <- ifelse(genes %in% sym_ok, genes,
           ifelse(genes %in% names(alias_map), alias_map[genes], NA_character_))
    out[is.na(out)] <- NA_character_
    out
  }
  reconstruct_old <- function(g) {
    x <- readRDS(proj_path("data", "source_v1", "expr", paste0(g, "_symbol_v1.rds")))
    genes <- map_genes(rownames(x$expr))
    e0 <- x$expr[!is.na(genes), , drop = FALSE]
    rownames(e0) <- genes[!is.na(genes)]
    e0 <- collapse_dup(e0)
    is_array <- geo$data_type[geo$gse == g][1] == "array"
    if (max(e0) >= CFG$log_scale_max) {
      is_counts <- all(e0 >= 0) && mean(abs(e0 - round(e0)) < 1e-8) > 0.99
      if (is_counts) e0 <- edgeR::cpm(e0, log = TRUE, prior.count = 0.5)
      else if (is_array) e0 <- log2(e0 + 1)
      else e0 <- log2(e0 + 1)
    }
    e0
  }

  build_new <- function(sm, map) {
    sym <- map[as.character(sm$ids)]
    hit <- mean(!is.na(sym))
    m <- sm$mat[!is.na(sym), , drop = FALSE]
    rownames(m) <- unname(sym[!is.na(sym)])
    g2 <- map_genes(rownames(m))          # HGNC normalisation (same treatment of v1 matrices as script 11)
    m <- m[!is.na(g2), , drop = FALSE]
    rownames(m) <- g2[!is.na(g2)]
    m <- collapse_dup(m)
    if (max(m) >= CFG$log_scale_max) {
      check(all(m >= 0), "rebuild matrix needs a log transform but contains negative values")
      m <- log2(m + 1)
    }
    check(max(m) < CFG$log_scale_max, "rebuild matrix still looks unlogged after transformation")
    check(!anyNA(m), "rebuild matrix contains NA")
    list(expr = m, hit = hit, n_probes_annot = sum(!is.na(sym)))
  }

  spec <- list(
    GSE12274 = list(gpl = "GPL570",  file = "GSE12274/GSE12274_series_matrix.txt.gz"),
    GSE48022 = list(gpl = "GPL570",  file = "GSE48022/GSE48022_series_matrix.txt.gz"),
    GSE57151 = list(gpl = "GPL16686", file = "GSE57151/GSE57151_series_matrix.txt.gz"),
    GSE19773 = list(gpl = "GPL7363", file = "GSE19773/GSE19773_series_matrix.txt.gz")
  )

  # ---- 2. per-study rebuild + per-gene correlation report ------------------------------------------
  rep_rows <- list()
  for (g in names(spec)) {
    old <- readRDS(proj_path("data", "studies", paste0(g, ".rds")))
    old0 <- reconstruct_old(g)
    check(all(old$meta$gsm %in% colnames(old0)), g, ": v1 rebuild baseline is missing samples")
    old0 <- old0[, match(old$meta$gsm, colnames(old0)), drop = FALSE]
    matched <- identical(dim(old$expr), dim(old0)) &&
               isTRUE(all.equal(old$expr, old0, tolerance = 1e-10, check.attributes = FALSE))

    sm <- load_sm(proj_path("data", "raw", spec[[g]]$file))
    check(all(old$meta$gsm %in% colnames(sm$mat)), g, ": matrix is missing study samples: ",
          paste(setdiff(old$meta$gsm, colnames(sm$mat)), collapse = ","))
    sm$mat <- sm$mat[, match(old$meta$gsm, colnames(sm$mat)), drop = FALSE]
    nb <- build_new(sm, maps[[spec[[g]]$gpl]]$map)

    shared <- intersect(rownames(old0), rownames(nb$expr))
    check(length(shared) >= 2000, g, ": shared genes between old and new: ", length(shared))
    rho <- vapply(shared, function(gene)
      suppressWarnings(stats::cor(old0[gene, ], nb$expr[gene, ], method = "pearson")), numeric(1))
    rs <- vapply(shared, function(gene)
      suppressWarnings(stats::cor(old0[gene, ], nb$expr[gene, ], method = "spearman")), numeric(1))

    rep_rows[[g]] <- data.frame(
      gse = g, platform = spec[[g]]$gpl,
      n_probes_total = length(sm$ids), n_probes_annot = nb$n_probes_annot,
      annot_hit = round(nb$hit, 4),
      n_genes_old = nrow(old0), n_genes_new = nrow(nb$expr), n_shared = length(shared),
      n_only_old = nrow(old0) - length(shared), n_only_new = nrow(nb$expr) - length(shared),
      corr_pearson_median = round(stats::median(rho, na.rm = TRUE), 5),
      corr_pearson_q05 = round(unname(stats::quantile(rho, 0.05, na.rm = TRUE)), 5),
      corr_pearson_min = round(min(rho, na.rm = TRUE), 5),
      corr_spearman_median = round(stats::median(rs, na.rm = TRUE), 5),
      n_corr_lt_0.95 = sum(rho < 0.95, na.rm = TRUE),
      old_vs_current_matched = matched,
      samples = ncol(nb$expr), stringsAsFactors = FALSE)

    old$expr <- nb$expr
    validate_study(old)
    saveRDS(old, proj_path("data", "studies", paste0(g, ".rds")))
    log_msg(g, ": rebuilt ", nrow(nb$expr), " genes (v1 baseline ", nrow(old0),
            "; old object = baseline ", matched, "), Pearson median ",
            rep_rows[[g]]$corr_pearson_median, ", annotation rate ", rep_rows[[g]]$annot_hit)
  }
  rep <- do.call(rbind, rep_rows)
  write_tsv(rep, "data/m1_rebuild_report_v2.tsv")

  # ---- 3. update the matching coverage_audit rows (annotation rate/gene count)------------------------------
  aud <- read_tsv("data/coverage_audit_v2.tsv")
  for (g in names(spec)) {
    i <- which(aud$gse == g)
    if (length(i) == 1) {
      aud$n_genes_mapped[i] <- rep$n_genes_new[rep$gse == g]
      aud$annot_rate[i] <- rep$annot_hit[rep$gse == g]
      aud$notes[i] <- paste0(aud$notes[i], "; M1 probe rebuild ", "2026-10-11")
    }
  }
  write_tsv(aud, "data/coverage_audit_v2.tsv")

  # ---- 4. acceptance criteria summary ----------------------------------------------------------
  check(nrow(rep) == 4, "unexpected report row count")
  log_msg("M1 rebuild done: validate_study passed for all 4 studies; report data/m1_rebuild_report_v2.tsv")
  invisible(rep)
})
