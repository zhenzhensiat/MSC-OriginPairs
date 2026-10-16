# =============================================================================
# utils_build.R —— GEO file -> study-object builder (verified code, moved from 46_build_c_layer.R)
# Purpose: scripts 53/54/61/70 reuse the same pipeline logic.
#   read_geo_study(gse, spec, gsm_keep, platform_hint): reads a study from raw GEO
#   downloads ("auto" mode: series matrix -> probe annotation -> cleanup), honouring
#   gsm_keep for column selection; the GPL hint comes from platform_hint, the series
#   matrix file name, or the c_layer_studies table, in that order.
# Dependency: org.Hs.eg.db (symbol mapping) must be available before loading this file.
# =============================================================================

# Dependency: org.Hs.eg.db (symbol mapping) must be available before loading this file
suppressMessages(library(org.Hs.eg.db))

# ---- Shared utilities (migrated from 46_build_c_layer.R, verified)------------------------------
eg <- AnnotationDbi::select(org.Hs.eg.db, keys = AnnotationDbi::keys(org.Hs.eg.db),
                            columns = c("SYMBOL", "ALIAS"))
sym_ok <- unique(na.omit(eg$SYMBOL))
al <- eg[!is.na(eg$ALIAS) & !is.na(eg$SYMBOL), c("ALIAS", "SYMBOL")]
al_n <- table(al$ALIAS)
ua <- names(al_n[al_n == 1L])
alias_map <- setNames(al$SYMBOL[match(ua, al$ALIAS)], ua)
# ENSG (common in RNA-seq expression tables) → unique symbol mapping
en <- AnnotationDbi::select(org.Hs.eg.db, keys = AnnotationDbi::keys(org.Hs.eg.db, "ENSEMBL"),
                            columns = "SYMBOL", keytype = "ENSEMBL")
en <- en[!is.na(en$ENSEMBL) & !is.na(en$SYMBOL), ]
en_n <- table(en$ENSEMBL); en_u <- names(en_n[en_n == 1L])
ens_map <- setNames(en$SYMBOL[match(en_u, en$ENSEMBL)], en_u)
# ENSEMBLTRANS (transcripts, GSE265844 Salmon output) → symbols; strip the version suffix first
et <- AnnotationDbi::select(org.Hs.eg.db, keys = AnnotationDbi::keys(org.Hs.eg.db, "ENSEMBLTRANS"),
                            columns = "SYMBOL", keytype = "ENSEMBLTRANS")
et <- et[!is.na(et$ENSEMBLTRANS) & !is.na(et$SYMBOL), ]
et_n <- table(et$ENSEMBLTRANS); et_u <- names(et_n[et_n == 1L])
etr_map <- setNames(et$SYMBOL[match(et_u, et$ENSEMBLTRANS)], et_u)
map_ensembltrans <- function(x) unname(etr_map[sub("\\.[0-9]+$", "", x)])
# ENTREZID → symbols (the GeneID column of GSE178514 is Entrez)
ez <- AnnotationDbi::select(org.Hs.eg.db, keys = AnnotationDbi::keys(org.Hs.eg.db, "ENTREZID"),
                            columns = "SYMBOL", keytype = "ENTREZID")
ez <- ez[!is.na(ez$ENTREZID) & !is.na(ez$SYMBOL), ]
ez_n <- table(ez$ENTREZID); ez_u <- names(ez_n[ez_n == 1L])
ez_map <- setNames(ez$SYMBOL[match(ez_u, ez$ENTREZID)], ez_u)
# GPL platform annotation: probe ID → symbol (mechanism ported from annot_gz in v1 scripts/02b_bind: NCBI FTP local
# cache + retry backoff). Differences from v1: (1) the completeness gate uses parseability instead of HEAD content-length; (2) no
# annot.gz platforms (7088/14550/16297/17586/28102 etc.) fall back to the SOFT platform table; (3) column names are auto-detected
# (Gene symbol / GeneSymbol / gene_assignment -- the latter takes its symbol segment by " // ").
annot_probe2sym <- function(gpl) {
  gpl <- as.character(gpl)
  check(grepl("^[0-9]+$", gpl), "invalid GPL number: ", gpl)
  pre <- if (as.integer(gpl) >= 1000) paste0(substr(gpl, 1, nchar(gpl) - 3), "nnn") else "nnn"
  base <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL", pre, "/GPL", gpl)
  dir.create(proj_path("data", "gpl_annot"), showWarnings = FALSE)
  curl_bin <- Sys.which("curl")
  old_to <- getOption("timeout"); on.exit(options(timeout = old_to), add = TRUE)
  options(timeout = max(600, old_to))
  # downloads use system curl (NCBI needs --ssl-no-revoke; R's bundled libcurl reports SSL connect error on this site,
  # observed 2026-10-12; --fail: no partial response body on 404); fall back to download.file when curl is unavailable
  dl <- function(url, dest) {
    for (att in 1:4) {
      if (file.exists(dest) && file.size(dest) > 0) return(TRUE)
      if (nzchar(curl_bin)) {
        tryCatch(system2(curl_bin, c("-sS", "--fail", "--ssl-no-revoke", "--max-time", "900", "-o", dest, url),
                         stdout = FALSE, stderr = FALSE), error = function(e) NULL)
      } else {
        tryCatch(utils::download.file(url, dest, quiet = TRUE, mode = "wb"), error = function(e) NULL)
      }
      if (file.exists(dest) && file.size(dest) > 0) return(TRUE)
      Sys.sleep(8 * att)
    }
    FALSE
  }
  # gene_assignment form (Affymetrix HTA): "NM_x // SYMBOL // desc /// …" → SYMBOL
  extract_sym <- function(z) {
    h <- grepl(" // ", z, fixed = TRUE)
    if (any(h)) {
      z[h] <- vapply(strsplit(z[h], " /// ", fixed = TRUE), function(seg) {
        p <- strsplit(seg[1], " // ", fixed = TRUE)[[1]]
        if (length(p) >= 2L) trimws(p[2]) else ""
      }, "")
    }
    trimws(z)
  }
  # (1) primary: annot.gz
  dest_a <- proj_path("data", "gpl_annot", paste0("GPL", gpl, ".annot.gz"))
  if (dl(paste0(base, "/annot/GPL", gpl, ".annot.gz"), dest_a)) {
    p2s <- tryCatch({
      con <- gzfile(dest_a, "rt"); ln <- readLines(con, warn = FALSE); close(con)
      hi <- grep("^ID\t.*Gene symbol", ln)[1]
      if (is.na(hi)) NULL else {
        hd <- strsplit(ln[hi], "\t", fixed = TRUE)[[1]]
        sym_col <- grep("^Gene symbol$", hd)[1]
        if (is.na(sym_col)) NULL else {
          body <- ln[(hi + 1):length(ln)]
          body <- body[grepl("\t", body, fixed = TRUE)]
          m <- do.call(rbind, strsplit(body, "\t", fixed = TRUE))
          keep <- nzchar(m[, sym_col])
          setNames(m[keep, sym_col], m[keep, 1])
        }
      }
    }, error = function(e) NULL)
    if (!is.null(p2s) && length(p2s) > 0) return(p2s)
    unlink(dest_a)   # downloaded but unparseable (bad response body): remove, retry next round
  }
  # (2) fall back to SOFT platform table (!platform_table_begin..end; chunked read to avoid loading large files)
  dest_s <- proj_path("data", "gpl_annot", paste0("GPL", gpl, "_family.soft.gz"))
  if (dl(paste0(base, "/soft/GPL", gpl, "_family.soft.gz"), dest_s)) {
    p2s <- tryCatch({
      con <- gzfile(dest_s, "rt"); on.exit(close(con))
      in_tab <- FALSE; hd <- NULL; sym_col <- NA_integer_
      out_id <- character(0); out_sym <- character(0)
      repeat {
        ln <- readLines(con, n = 200000L, warn = FALSE)
        if (!length(ln)) break
        if (!in_tab) {
          b <- grep("^!platform_table_begin", ln)
          if (!length(b)) next
          in_tab <- TRUE
          ln <- ln[(b[1] + 1):length(ln)]
        }
        done <- FALSE
        e <- grep("^!platform_table_end", ln)
        if (length(e)) { ln <- ln[seq_len(e[1] - 1)]; done <- TRUE }
        if (is.null(hd)) {
          if (!length(ln)) { if (done) break else next }
          hd <- strsplit(ln[1], "\t", fixed = TRUE)[[1]]
          sym_col <- grep("^(gene[ _]?symbol|symbol)$", hd, ignore.case = TRUE)[1]
          if (is.na(sym_col)) sym_col <- grep("^gene[ _]?name$", hd, ignore.case = TRUE)[1]  # Agilent (GPL13607; approved 2026-10-13)
          if (is.na(sym_col)) sym_col <- grep("^gene_assignment$", hd, ignore.case = TRUE)[1]
          if (is.na(sym_col)) sym_col <- grep("^orf$", hd, ignore.case = TRUE)[1]   # GPL7088 (ORF column holds symbols)
          if (is.na(sym_col)) break    # no symbol column in table: give up (p2s=NULL)
          ln <- ln[-1]
        }
        if (length(ln)) {
          sp <- strsplit(ln, "\t", fixed = TRUE)
          okk <- lengths(sp) >= sym_col
          if (any(okk)) {
            out_id <- c(out_id, vapply(sp[okk], `[`, "", 1L))
            out_sym <- c(out_sym, vapply(sp[okk], `[`, "", sym_col))
          }
        }
        if (done) break
      }
      if (!length(out_id)) NULL else {
        sy <- extract_sym(out_sym)
        keep <- nzchar(sy)
        setNames(sy[keep], out_id[keep])
      }
    }, error = function(e) NULL)
    if (!is.null(p2s) && length(p2s) > 0) return(p2s)
    unlink(dest_s)
  }
  NULL
}
map_genes <- function(genes) {
  out <- ifelse(genes %in% sym_ok, genes,
         ifelse(genes %in% names(alias_map), alias_map[genes], NA_character_))
  miss <- is.na(out) & grepl("^ENSG[0-9]+", genes)
  if (any(miss)) out[miss] <- ens_map[genes[miss]]
  out[is.na(out)] <- NA_character_
  out
}
collapse_dup <- function(x) {   # multiple rows per symbol: highest mean wins, ties by IQR (same rule as script 11)
  if (!anyDuplicated(rownames(x))) return(x)
  mu <- rowMeans(x); iqr <- apply(x, 1, stats::IQR)
  ord <- order(mu, iqr, decreasing = TRUE)
  keep <- ord[!duplicated(rownames(x)[ord])]
  x[keep, , drop = FALSE]
}
parse_passage <- function(txt) {
  hit <- grepl("passage[^0-9]*([0-9]+)|\\bP[0-9]{1,2}\\b", txt, ignore.case = TRUE)
  vapply(seq_along(txt), function(i) {
    if (!hit[i]) return(NA_real_)
    up <- toupper(txt[i])
    m <- regmatches(up, regexpr("P[0-9]{1,2}", up))
    if (length(m)) return(as.numeric(sub("P", "", m)))
    as.numeric(sub(".*passage[^0-9]*([0-9]+).*", "\\1", txt[i], ignore.case = TRUE))
  }, numeric(1))
}
parse_medium <- function(txt) {
  t <- tolower(txt)
  ifelse(grepl("platelet lysate|\\bhpl\\b", t), "hPL",
  ifelse(grepl("fetal bovine|\\bfbs\\b", t) & !grepl("platelet", t), "FBS",
  ifelse(grepl("human serum|\\bhs\\b", t) & !grepl("fetal", t), "HS",
  ifelse(grepl("serum.?free|xeno.?free|\\bxf\\b", t), "xeno_free", NA_character_))))
}
parse_sex <- function(txt) {
  t <- tolower(txt)
  f_hit <- grepl("female", t) | grepl("(gender|sex)[ :=]+f($|[^a-z])", t, perl = TRUE)
  m_hit <- (grepl("male", t) & !grepl("female", t)) | grepl("(gender|sex)[ :=]+m($|[^a-z])", t, perl = TRUE)
  ifelse(f_hit & !m_hit, "F", ifelse(m_hit & !f_hit, "M", NA_character_))
}
parse_donor <- function(txt) {
  # layered donor-ID extraction: (1) "donor id: X" (2) donor/patient/subject + identifier (excluding modifiers)
  # (3) known prefixes with numbers (BMSC-3 / HD-4 / L-ASCs-1 …) (4) fallback: first number (TTL-coded studies)
  drop <- c("replicate", "control", "hypoxia", "normal", "high", "low", "glucose",
            "medium", "serum", "free", "passage", "early", "late", "treated", "untreated",
            "primed", "blank", "isolated", "adhered", "donor", "patient", "subject")
  vapply(trimws(txt), function(s) {
    m <- regmatches(s, regexpr("(?i)donor[ _.-]*id[ _.:-]*([A-Za-z0-9_-]+)", s, perl = TRUE))
    if (length(m)) return(sub("(?i)donor[ _.-]*id[ _.:-]*", "", m))
    m <- regmatches(s, regexpr("(?i)(donor|patient|subject|individual)[ _.,:-]+([A-Za-z0-9_-]+)",
                               s, perl = TRUE))
    if (length(m)) {
      v <- sub("(?i)(donor|patient|subject|individual)[ _.,:-]+", "", m)
      if (!tolower(v) %in% drop && !grepl("^[0-9]*$", v)) return(v)
      if (grepl("[A-Za-z]{2,}", v)) return(v)
    }
    m <- regmatches(s, regexpr("(?i)\\b(BMSC|BM|HD|CF|FO|SUB[0-9]+|L-ASCs|O-ASCs|hUC\\.MSC|UC-MSC)[ _.-]*[0-9]+", s, perl = TRUE))
    if (length(m)) return(m)
    m <- regmatches(s, regexpr("[0-9]+", s))
    if (length(m)) return(m)
    NA_character_
  }, character(1))
}
y_genes <- c("DDX3Y", "RPS4Y1", "KDM5D", "EIF1AY", "USP9Y", "ZFY", "UTY")

# ---- File-reading utilities ---------------------------------------------------------------
read_series_table <- function(path) {          # series matrix data table
  ln <- readLines(gzfile(path, "rt"), warn = FALSE)
  b <- grep("^!series_matrix_table_begin", ln); e <- grep("^!series_matrix_table_end", ln)
  check(length(b) == 1 && length(e) == 1 && e > b, basename(path), ": malformed table markers (suspected E7)")
  tb <- ln[(b + 1):(e - 1)]; tb <- tb[nzchar(tb)]
  df <- read.delim(text = paste(tb, collapse = "\n"), header = TRUE, quote = "\"",
                   check.names = FALSE, stringsAsFactors = FALSE)
  id_col <- names(df)[1]
  mat <- as.matrix(df[, -1, drop = FALSE])
  rownames(mat) <- as.character(df[[id_col]])
  mat
}
read_table_auto <- function(path) {   # gz/plain table: auto-locate header + sep
  ln <- readLines(gzfile(path, "rt"), warn = FALSE)
  ln <- ln[nzchar(ln)]
  # header detection: first row satisfying (1) not starting with #/! (2) fields ≥2 and same as the next row (3) almost no numeric fields in this row
  # (<0.1; trait of a real header - the group-title row "\tRaw data (RC)" also matches it but not (4)) (4) the next row contains at least some numeric
  # fields (≥0.1; sparse data rows allowed, e.g. rows 2-3 of GSE130393 with many empty values). GSE200086 two-row header:
  # the row after row 1 is another header (num_frac=0) → keep moving down; row 2 passes.
  sep0 <- function(s) if (grepl("\t", s)) "\t" else if (grepl(",", s)) "," else "\t"
  n_field <- function(s) length(strsplit(s, sep0(s))[[1]])
  num_frac <- function(v) mean(!is.na(suppressWarnings(as.numeric(v))))
  start <- 1L
  for (i in seq_len(max(1L, min(30L, length(ln) - 1L)))) {
    if (grepl("^[#!]", ln[i])) next
    f <- strsplit(ln[i], sep0(ln[i]))[[1]]
    f2 <- strsplit(ln[i + 1], sep0(ln[i]))[[1]]
    if (length(f) >= 2L && length(f2) == length(f) &&
        num_frac(f) < 0.1 && num_frac(f2) >= 0.1) { start <- i; break }
  }
  ln <- ln[start:length(ln)]
  sep <- if (grepl("\t", ln[1])) "\t" else if (grepl(",", ln[1])) "," else ""
  df <- read.delim(text = paste(ln, collapse = "\n"), header = TRUE, quote = "\"",
                   sep = sep, check.names = FALSE, stringsAsFactors = FALSE)
  df
}
read_xlsx_any <- function(path) {              # covers .xlsx / .xls / .xls.gz
  if (grepl("\\.gz$", path)) {
    tmp <- tempfile(fileext = paste0(".", sub("\\.gz$", "", tools::file_ext(sub("\\.gz$", "x", path)))))
    con <- gzfile(path, "rb"); raw <- readBin(con, "raw", n = 5e8); close(con)
    writeBin(raw, tmp); path <- tmp
  }
  as.data.frame(readxl::read_excel(path, col_names = TRUE, .name_repair = "minimal"))
}
find_gene_col <- function(df) {
  # symbol-like columns take precedence over ID-like columns (e.g. tables where GeneName follows GeneId)
  p1 <- grep("^(symbol|gene.?symbol|gene.?name|hugo|hgnc)", names(df), ignore.case = TRUE)
  if (length(p1)) return(p1[1])
  p2 <- grep("^(gene|gene.?id|id|name|ensg|entrez)", names(df), ignore.case = TRUE)
  if (length(p2)) return(p2[1])
  for (j in seq_along(df)) if (is.character(df[[j]]) || is.factor(df[[j]])) return(j)
  NA_integer_
}
to_matrix <- function(df, gene_col) {
  g <- as.character(df[[gene_col]])
  x <- as.matrix(df[, setdiff(seq_along(df), gene_col), drop = FALSE])
  if (!is.numeric(x)) {
    x2 <- apply(x, 2, function(v) suppressWarnings(as.numeric(v)))
    x <- x2
  }
  ok <- !is.na(g) & nzchar(g) & rowSums(is.na(x)) < ncol(x)   # drop all-NA rows
  x <- x[ok, , drop = FALSE]
  rownames(x) <- g[ok]
  x
}
# scale normalisation: counts → cpm(log); linear → log2(x+1); already-logged kept (same rule as script 11)
normalize_scale <- function(x) {
  if (max(x, na.rm = TRUE) >= CFG$log_scale_max) {
    is_counts <- all(x >= 0, na.rm = TRUE) &&
      mean(abs(x - round(x)) < 1e-8, na.rm = TRUE) > 0.99
    if (is_counts) return(list(expr = edgeR::cpm(x, log = TRUE, prior.count = 0.5),
                               scale = "logCPM", counts = x))
    return(list(expr = log2(x + 1), scale = "log2TPM", counts = NULL))
  }
  list(expr = x, scale = "array", counts = NULL)
}
strip_affix <- function(x) {
  x2 <- sub("^(FPKM|TPM|CPM|RPKM|Counts?|Count)[._-]", "", x, ignore.case = TRUE)
  sub("[._-]?(FPKM|TPM|CPM|RPKM|Counts?|Count)$", "", x2, ignore.case = TRUE)
}
match_cols_gsm <- function(cols, hdr) {        # column names → GSM (series-matrix metadata bridge)
  norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  ng <- norm(hdr$Sample_geo_accession); nt <- norm(hdr$Sample_title)
  nd <- if (!is.null(hdr$Sample_description)) norm(hdr$Sample_description) else NULL
  nc <- norm(cols)
  nc2 <- norm(strip_affix(cols))
  out <- rep(NA_character_, length(cols))
  for (i in seq_along(cols)) {
    for (cand in unique(c(nc[i], nc2[i]))) {   # original name → name without quantifier affixes (FPKM.FO1 O → FO1 O)
      if (nchar(cand) < 2) next
      j <- match(cand, ng)
      if (is.na(j)) j <- match(cand, nt)
      if (is.na(j) && !is.null(nd)) j <- match(cand, nd)   # description fallback (GSE165811 X13-X57; approved 2026-10-13, exact layer only)
      if (is.na(j)) {
        h1 <- which(startsWith(nt, cand))       # candidate is a title prefix (unique hit)
        if (length(h1) == 1L) j <- h1
        if (is.na(j)) {                         # title is a prefix of the candidate (take the longest title; e.g. hUC.MSC10_FRAS… vs 10/1)
          h2 <- which(startsWith(cand, nt))
          if (length(h2)) j <- h2[which.max(nchar(nt[h2]))]
        }
      }
      if (is.na(j) && nchar(cand) >= 4L) {      # unique substring (e.g. ctrl1 ∈ …untreatedctrl1)
        h <- which(grepl(cand, nt, fixed = TRUE))
        if (length(h) == 1L) j <- h
      }
      if (!is.na(j)) { out[i] <- hdr$Sample_geo_accession[j]; break }
    }
  }
  out
}


# ---- spec table (migrated from 46; key = GSE, value = read mode)--------------------------------
SPEC <- list(
  GSE144605 = "matrix", GSE90149 = "matrix", GSE111495 = "matrix", GSE14828 = "matrix",
  GSE24871 = "matrix", GSE40033 = "matrix", GSE63099 = "matrix", GSE107214 = "matrix",
  GSE62665 = "matrix", GSE72597 = "matrix", GSE76295 = "matrix",
  GSE130393 = "special_130393",
  GSE217709 = "GSE217709_DKD.vs.controls_RPKM.tsv.gz",
  GSE118808 = "GSE118808_Comparative_analysis_of_MSCs_-Count_Fpkm.txt.gz",
  GSE165811 = "GSE165811_counts.csv.gz",
  GSE169533 = "special_169533",
  GSE224190 = "GSE224190_gene_FPKM.txt.gz",
  GSE269586 = "special_269586",
  GSE222263 = "special_222263",
  GSE308101 = "GSE308101_gene_fpkm_sub.xlsx",
  GSE215041 = "special_215041",
  GSE200086 = "special_tar_gsm",
  GSE265844 = "special_tar_gsm",
  GSE129801 = "special_tar_gsm",
  GSE189524 = "special_189524_txonly",
  GSE276163 = "special_276163",
  GSE178514 = "special_178514",
  GSE123973 = "special_123973",
  GSE269628 = "special_269628",
  GSE272480 = "special_272480",
  GSE310465 = "GSE310465_raw_counts.txt.gz"
)

# ---- read_geo_study -----------------------------------------------------------
# Read one GEO study: raw files -> probe mapping -> cleaning -> columns->GSM -> (optional keep columns) -> scale normalisation
# Returns list(expr, counts, scale, hdr, col_gsm, note). Errors are raised via stop() (with clues).
read_geo_study <- function(gse, spec = NULL, gsm_keep = NULL, platform_hint = NULL) {
  g <- gse
  sp <- if (is.null(spec)) SPEC[[g]] else spec
  if (is.null(sp)) sp <- "auto"
  note <- ""
  hdr <- tryCatch(
    read_series_header(proj_path("data", "raw", g, list.files(proj_path("data", "raw", g),
                                                              pattern = "series_matrix.*gz")[1])),
    error = function(e) NULL)
  raw_expr <- tryCatch({
    if (identical(sp, "auto")) {   # generic mode: read series matrix directly (row-name mapping left to the probe block)
      read_series_table(proj_path("data", "raw", g,
                                  list.files(proj_path("data", "raw", g), pattern = "series_matrix.*gz")[1]))
    } else
if (identical(sp, "matrix")) {
  read_series_table(proj_path("data", "raw", g,
                              list.files(proj_path("data", "raw", g), pattern = "series_matrix.*gz")[1]))
} else if (identical(sp, "special_130393")) {
  # wide table with interleaved columns (SAMPLE k alternating with Detection Pval): take only per-sample signal columns, SAMPLE k ↔ metadata row k
  df <- read_table_auto(proj_path("data", "raw", g, "GSE130393_normalized_data.txt.gz"))
  smp <- grep("^SAMPLE [0-9]+(\\.AVG_Signal)?$", names(df), value = TRUE)
  smp <- smp[order(as.integer(sub("^SAMPLE ([0-9]+).*$", "\\1", smp)))]
  check(length(smp) == nrow(hdr), "GSE130393 SAMPLE k column count does not match metadata row count")
  m <- as.matrix(df[, smp, drop = FALSE])
  rownames(m) <- as.character(df[["SYMBOL"]])
  colnames(m) <- hdr$Sample_geo_accession[as.integer(sub("^SAMPLE ([0-9]+).*$", "\\1", smp))]
  m
} else if (identical(sp, "special_222263")) {
  # Table S1: interleaved FPKM/count columns; take only -FPKM columns; column names O/L-ADSC-k ↔ titles O/L-ASCs-k
  df <- read_xlsx_any(proj_path("data", "raw", g, "GSE222263_processed_data_file.xlsx"))
  fpk <- grep("-FPKM$", names(df), value = TRUE)
  check(length(fpk) >= 3, "GSE222263: no FPKM columns found")
  m <- as.matrix(df[, fpk, drop = FALSE])
  rownames(m) <- as.character(df[["gene_name"]])
  key <- sub("ADSC", "ASCs", sub("-FPKM$", "", fpk))
  norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  j <- match(norm(key), norm(hdr$Sample_title))
  check(!anyNA(j), "GSE222263 column-to-title mapping failed: ", paste(key[is.na(j)], collapse = ","))
  colnames(m) <- hdr$Sample_geo_accession[j]
  m
} else if (identical(sp, "special_169533")) {
  # sample columns "EXP:k_MSC_control_k:FPKM" / "…DEP_k…": control word + number ↔ titles "Control k"/"DEP k"
  df <- read_table_auto(proj_path("data", "raw", g, "GSE169533_expression.genes.xls.gz"))
  smp <- grep(":FPKM$", names(df), value = TRUE)
  key <- tolower(sub("^EXP:[0-9]+_MSC_([A-Za-z]+)_[0-9]+:FPKM$", "\\1", smp))
  num <- sub("^.*_([0-9]+):FPKM$", "\\1", smp)
  key <- ifelse(key == "dep", "treated", key)   # title wording is "…Diesel exhaust particle treated k"
  cand <- paste0(key, num)
  nt <- gsub("[^a-z0-9]", "", tolower(hdr$Sample_title))
  j <- match(cand, nt)
  for (i in which(is.na(j))) {
    h <- which(grepl(cand[i], nt, fixed = TRUE))
    if (length(h) == 1L) j[i] <- h
  }
  check(!anyNA(j), "GSE169533 column mapping failed: ", paste(cand[is.na(j)], collapse = ","))
  m <- as.matrix(df[, smp, drop = FALSE])
  rownames(m) <- as.character(df[["GeneName"]])
  colnames(m) <- hdr$Sample_geo_accession[j]
  m
} else if (identical(sp, "special_269586")) {
  # P0 (master cell bank) samples: columns raw_BM-CEM-GW-<donor> (counts); donor number ↔ title "donor N".
  # column names carry no passage, and P2/P3 donor numbers overlap those of P0 -- only the keep (P0×RNAseq) samples are mapped exactly,
  # assigning keep GSMs directly (P2/P3 columns do not participate in this object).
  keep_g0 <- gsm_keep
  tit0 <- hdr$Sample_title[match(keep_g0, hdr$Sample_geo_accession)]
  check(!anyNA(tit0), "GSE269586 keep GSMs not present in the GPL24676 matrix")
  dn <- sub(".*donor ([0-9]+).*", "\\1", tit0)
  df <- read_table_auto(proj_path("data", "raw", g, "GSE269586_all_counts_merged.txt.gz"))
  j <- match(paste0("raw_BM-CEM-GW-", dn), names(df))
  check(!anyNA(j), "GSE269586 column mapping failed: ", paste(dn[is.na(j)], collapse = ","))
  m <- as.matrix(df[, names(df)[j], drop = FALSE])
  rownames(m) <- as.character(df[["HGNC_symbol"]])
  colnames(m) <- keep_g0
  m
} else if (identical(sp, "special_215041")) {
  # limma output table: row 1 is a group title, row 2 the true header (skip=1). Columns 9-14 = log2 normalized;
  # columns 15-20 = per-sample read counts (integers, spot-checked over 3 rows 2026-10-12); take the counts set → cpm(log), consistent with the layer-D convention.
  # column names repeat from column 9 onward (WJ_1.. appears twice) -- always index by position, never by name.
  df <- as.data.frame(readxl::read_excel(
    proj_path("data", "raw", g, "GSE215041_wj_sg_MSC_RNAseq.xlsx"),
    skip = 1, .name_repair = "minimal"))
  pos <- 15:20
  keep_g0 <- gsm_keep
  tit0 <- hdr$Sample_title[match(keep_g0, hdr$Sample_geo_accession)]
  check(!anyNA(tit0) && all(grepl("wharton", tit0, ignore.case = TRUE)),
        "GSE215041: not all keep samples are wjMSC")
  j <- match(paste0("WJ_", sub(".*-([0-9]+)$", "\\1", tit0)), names(df)[pos])
  check(!anyNA(j), "GSE215041 column mapping failed")
  m <- as.matrix(df[, pos[j], drop = FALSE])
  rownames(m) <- as.character(df[["Gene symbol"]])
  colnames(m) <- keep_g0
  m
} else if (identical(sp, "special_189524_txonly")) {
  stop("GSE189524: gene_name/gene_id are all UCSC transcript IDs (uc*, 182434/182434), no gene-symbol-level annotation -- not built")
} else if (identical(sp, "special_276163")) {
  # supplementary readcount file: gene_name × "Sample k" (column names ↔ !Sample_description exact match, independent of column order)
  check(!is.null(hdr), "GSE276163: series matrix header required for column mapping")
  check("Sample_description" %in% names(hdr), "GSE276163: header lacks Sample_description column")
  df <- read_table_auto(proj_path("data", "raw", g, "GSE276163_ASCs_mRNA_readcount.txt.gz"))
  smp <- grep("^Sample [0-9]+$", names(df), value = TRUE)
  j <- match(smp, hdr$Sample_description)
  check(!anyNA(j) && !anyDuplicated(j), "GSE276163: column mapping failed (", paste(smp[is.na(j)], collapse = ","), ")")
  m <- as.matrix(df[, smp, drop = FALSE])
  rownames(m) <- as.character(df[["gene_name"]])
  colnames(m) <- hdr$Sample_geo_accession[j]
  m
} else if (identical(sp, "special_178514")) {
  # supplementary FPKM file: GeneID(Entrez) × "P<code>-FPKM" (column code ↔ [P<code>] in the TTL;
  # match the number inside brackets by lookahead -- a plain "P[0-9]+" would first hit the passage code of "MSC_P4";
  # [P<code>] of miR_ samples shares numbers with mRNA samples, so it must be excluded before matching, or the first match may land on a miR row)
  check(!is.null(hdr), "GSE178514: series matrix header required for column mapping")
  df <- read_table_auto(proj_path("data", "raw", g, "GSE178514_AllGeneFPKM_mRNA.txt.gz"))
  fc <- grep("-FPKM$", names(df), value = TRUE)
  code <- sub("-FPKM$", "", fc)
  inc <- which(!grepl("^miR_", hdr$Sample_title))
  tc <- regmatches(hdr$Sample_title, regexpr("P[0-9]+(?=\\])", hdr$Sample_title, perl = TRUE))
  j <- match(code, tc[inc])
  check(!anyNA(j), "GSE178514: column mapping failed (", paste(code[is.na(j)], collapse = ","), ")")
  sym <- ez_map[as.character(df[["GeneID"]])]
  ok <- !is.na(sym)
  m <- as.matrix(df[ok, fc, drop = FALSE])
  rownames(m) <- sym[ok]
  colnames(m) <- hdr$Sample_geo_accession[inc[j]]
  m
} else if (identical(sp, "special_123973")) {
  # featureCounts table: Geneid(ENSG) × BAM columns (…US-<number>_…_accepted_hits.bam); a sample split across lanes
  # (multiple columns sharing a US number) → aggregated by summing per US number; US numbers ↔ TTL "US-…: …" (set-equality assertion; approved 2026-10-13)
  check(!is.null(hdr), "GSE123973: series matrix header required for column mapping")
  df <- read_table_auto(proj_path("data", "raw", g, "GSE123973_DPSC_counts_matrix.txt.gz"))
  bam <- grep("accepted_hits\\.bam$", names(df), value = TRUE)
  check(length(bam) >= 3L, "GSE123973: too few BAM columns")
  check(all(vapply(df[, bam, drop = FALSE], is.numeric, logical(1))), "GSE123973: BAM columns are not numeric")
  us <- sub("^.*?(US-[0-9]+).*$", "\\1", bam)
  tu <- sub("^(US-[0-9]+):.*$", "\\1", hdr$Sample_title)
  check(length(unique(us)) == length(unique(tu)) && setequal(unique(us), unique(tu)) && !anyDuplicated(tu),
        "GSE123973: US-number set does not match metadata")
  su <- sort(unique(us))
  vals <- vapply(su, function(u) rowSums(as.matrix(df[, bam[us == u], drop = FALSE])), numeric(nrow(df)))
  m <- vals                          # vapply output is genes(nrow) × samples(length(su))
  rownames(m) <- as.character(df[[1]])
  colnames(m) <- hdr$Sample_geo_accession[match(su, tu)]
  m
} else if (identical(sp, "special_269628")) {
  # headerless first-field file: line 1 = sample column names (I11+NNN, gene column name missing); data rows = ENSG + numeric columns;
  # last three digits of the I number ↔ first TTL segment number (set-equality assertion; approved 2026-10-13)
  check(!is.null(hdr), "GSE269628: series matrix header required for column mapping")
  con <- gzfile(proj_path("data", "raw", g, "GSE269628_raw_counts_OxLDL.txt.gz"), "rt")
  ln <- readLines(con, warn = FALSE); close(con)
  h1 <- gsub('^"|"$', "", strsplit(ln[1], "\t", fixed = TRUE)[[1]])
  tb <- read.delim(text = paste(ln[-1], collapse = "\n"), header = FALSE, quote = "\"",
                   check.names = FALSE, stringsAsFactors = FALSE)
  check(ncol(tb) == length(h1) + 1L, "GSE269628: column count does not match header")
  num <- sub("^I11", "", h1)
  tn <- sub("_.*$", "", hdr$Sample_title)
  check(length(num) == nrow(hdr) && !anyDuplicated(num) && setequal(num, tn),
        "GSE269628: column-number set does not match metadata")
  m <- as.matrix(tb[, -1, drop = FALSE])
  rownames(m) <- as.character(tb[[1]])
  colnames(m) <- hdr$Sample_geo_accession[match(num, tn)]
  m
} else if (identical(sp, "special_272480")) {
  # Partek transposed export: rows = samples (column 1 Sample name + 9 annotation columns), genes from column 11 on;
  # "Sample name" ↔ TTL normalised matching (extra rows discarded); transposed back to genes × samples (approved 2026-10-13)
  check(!is.null(hdr), "GSE272480: series matrix header required for column mapping")
  df <- read_table_auto(proj_path("data", "raw", g, "GSE272480_Partek_bam_files_Normalization_Normalized_counts.txt.gz"))
  gi <- 11L:ncol(df)
  sm <- as.character(df[[1]])
  ns <- gsub("[^a-z0-9]", "", tolower(sm))
  nt <- gsub("[^a-z0-9]", "", tolower(hdr$Sample_title))
  j <- match(ns, nt)
  keep <- !is.na(j)
  check(sum(keep) >= 1L && !anyDuplicated(j[keep]), "GSE272480: sample-name mapping failed")
  mm <- do.call(cbind, lapply(df[keep, gi, drop = FALSE], function(v) suppressWarnings(as.numeric(v))))
  m <- t(mm)
  rownames(m) <- names(df)[gi]
  colnames(m) <- hdr$Sample_geo_accession[j[keep]]
  m
} else if (identical(sp, "special_tar_gsm")) {
  # one file per GSM inside the tar (GSE200086 generic table format / GSE265844 Salmon transcript-level / GSE129801 headerless two-column)
  td <- tempfile(); dir.create(td)
  utils::untar(proj_path("data", "raw", g, paste0(g, "_RAW.tar")), exdir = td)
  fs <- list.files(td, full.names = TRUE)
  mats <- list()
  for (f in fs) {
    gm_f <- sub("_.*$", "", basename(f))
    one <- tryCatch({
      if (identical(g, "GSE265844")) {
        # Salmon quant: Name(ENST) Length EffectiveLength TPM NumReads → sum TPM per gene
        d <- read_table_auto(f)
        v <- as.numeric(d[["TPM"]])
        syn <- map_ensembltrans(as.character(d[["Name"]]))
        okv <- !is.na(syn) & !is.na(v)
        agg <- tapply(v[okv], syn[okv], sum)
        matrix(agg, ncol = 1, dimnames = list(names(agg), gm_f))
      } else if (identical(g, "GSE129801")) {
        # headerless two columns (gene symbol + counts)
        ln2 <- readLines(gzfile(f, "rt"), warn = FALSE)
        ln2 <- ln2[nzchar(ln2)]
        sp2 <- strsplit(ln2, "\t")
        okr <- lengths(sp2) >= 2L
        gs <- vapply(sp2[okr], `[`, "", 1L)
        vs <- suppressWarnings(as.numeric(vapply(sp2[okr], `[`, "", 2L)))
        okv <- !is.na(vs) & nzchar(gs)
        agg <- tapply(vs[okv], gs[okv], sum)
        matrix(agg, ncol = 1, dimnames = list(names(agg), gm_f))
      } else {
        d <- read_table_auto(f)
        to_matrix(d, find_gene_col(d))[, 1, drop = FALSE]
      }
    }, error = function(e) NULL)
    if (is.null(one) || nrow(one) < 100) next
    colnames(one) <- gm_f   # unify column names as GSM (some branches skip it, some set it wrong; always overwrite at the exit)
    mats[[gm_f]] <- one
  }
  check(length(mats) >= 3, "too few valid sample files in tar: ", length(mats))
  common <- Reduce(intersect, lapply(mats, rownames))
  do.call(cbind, lapply(mats, function(mm) mm[common, , drop = FALSE]))
} else if (grepl("\\.tar$", sp)) {
  td <- tempfile(); dir.create(td)
  utils::untar(proj_path("data", "raw", g, sp), exdir = td)
  cand <- list.files(td, recursive = TRUE, full.names = TRUE,
                     pattern = "\\.(txt|csv|tsv|xlsx|xls)(\\.gz)?$")
  if (!length(cand)) stop("no candidate matrix files in tar (raw format only?): ", paste(list.files(td), collapse = ","))
  f <- cand[which.max(file.size(cand))]     # take the largest candidate (the matrix is usually the largest)
  if (grepl("\\.xlsx?$", f)) to_matrix(read_xlsx_any(f), find_gene_col(read_xlsx_any(f)))
  else to_matrix(read_table_auto(f), find_gene_col(read_table_auto(f)))
} else if (grepl("\\.xlsx?(\\.gz)?$", sp)) {
  df <- read_xlsx_any(proj_path("data", "raw", g, sp))
  to_matrix(df, find_gene_col(df))
} else {
  df <- read_table_auto(proj_path("data", "raw", g, sp))
  to_matrix(df, find_gene_col(df))
}
    }, error = function(e) stop("read_fail: ", conditionMessage(e)))
# probe-ID row names (series matrix etc.) → GPL annotation mapping (enabled only when the symbol hit rate is very low; v1 02b mechanism)
if (mean(!is.na(map_genes(rownames(raw_expr)))) < 0.1) {
  smf <- list.files(proj_path("data", "raw", g), pattern = "series_matrix.*gz")
  m0 <- if (length(smf)) regmatches(smf[1], regexpr("GPL[0-9]+", smf[1])) else character(0)
  gpl0 <- if (length(m0)) sub("GPL", "", m0) else if (!is.null(platform_hint)) platform_hint else {
    cs <- tryCatch(read_tsv("curation/c_layer_studies_v2.tsv"), error = function(e) NULL)
    if (!is.null(cs) && any(cs$gse == g)) strsplit(cs$platform_gpl[cs$gse == g], ";")[[1]][1] else NA_character_
  }
  if (is.na(gpl0) && !is.null(hdr) && "Sample_platform_id" %in% names(hdr)) {   # final fallback in the chain (approved 2026-10-13)
    h1 <- unique(sub("^GPL", "", hdr$Sample_platform_id))
    if (length(h1) == 1L && grepl("^[0-9]+$", h1)) gpl0 <- h1
  }
  if (is.na(gpl0)) {
    note <- paste0(note, "no GPL hint, skipping probe annotation; ")
  } else {
    p2s <- tryCatch(annot_probe2sym(gpl0), error = function(e) NULL)
    if (is.null(p2s)) {
      note <- paste0(note, "probe annotation GPL", gpl0, " unavailable; ")
    } else {
      rn1 <- p2s[rownames(raw_expr)]
      hit <- mean(!is.na(rn1) & nzchar(rn1))
      note <- paste0(note, "probe annotation GPL", gpl0, " hit rate ", round(hit, 3), "; ")
      if (hit >= 0.2) {
        okp <- !is.na(rn1) & nzchar(rn1)
        raw_expr <- raw_expr[okp, , drop = FALSE]
        rownames(raw_expr) <- rn1[okp]
      }
    }
  }
}
gm <- map_genes(rownames(raw_expr))
hit_rate <- round(mean(!is.na(gm)), 4)
n_genes_raw <- length(gm)
if (sum(!is.na(gm)) < 100) {   # guard: matrices without usable symbol-level annotation must not reach downstream
  stop("no_symbol_mapping: hit ", sum(!is.na(gm)), " / ", length(gm))
}
e0 <- raw_expr[!is.na(gm), , drop = FALSE]
rownames(e0) <- gm[!is.na(gm)]
# all-NA columns (missing sample data) are removed first (logged) -- otherwise the 'drop NA rows' step below would empty the whole table because a single column is all-empty
col_all_na <- colSums(is.na(e0)) == nrow(e0)
if (any(col_all_na)) {
  note <- paste0(note, "all-NA columns dropped at read-in: ", sum(col_all_na), "; ")
  e0 <- e0[, !col_all_na, drop = FALSE]
}
# NA-containing rows are removed at read-in (logged; validate_study forbids NA and imputation). Must come before collapse:
# the IQR tie-break rule of collapse does not tolerate NA (script 11 feeds it a NA-cleaned v1 matrix; reading raw files directly requires cleaning first)
na_rows <- rowSums(is.na(e0)) > 0
if (any(na_rows)) {
  note <- paste0(note, "NA-containing rows dropped at read-in: ", sum(na_rows), "; ")
  e0 <- e0[!na_rows, , drop = FALSE]
}
if (nrow(e0) < 1000) {   # guard: matrices with <1000 genes after cleaning (panels/sparse data) cannot support genome-scale confirmation
  stop("too_few_genes_after_cleaning: ", nrow(e0), " rows")
}
e0 <- collapse_dup(e0)

col_gsm <- match_cols_gsm(colnames(e0), hdr)
keep_g <- if (is.null(gsm_keep)) col_gsm[!is.na(col_gsm)] else gsm_keep
# one GSM matching several columns (e.g. Count+Fpkm pairs: GSE118808): keep the first column (counts listed first, yielding cpm(log) consistent with layer D)
col_gsm[duplicated(col_gsm)] <- NA_character_
sel <- which(col_gsm %in% keep_g)
if (length(sel) < length(keep_g)) {
  miss <- setdiff(keep_g, col_gsm[sel])
  note <- paste0(note, "missing keep samples: ", length(miss), " (",
                      paste(head(miss, 5), collapse = ","), "); ")
}
if (length(sel) < 3) stop("too_few_keep_cols: ", length(sel))
e1 <- e0[, sel, drop = FALSE]
colnames(e1) <- col_gsm[sel]   # column names → GSM (validate_study requires positional agreement with meta$gsm)
# Drop all-NA / all-zero columns (per sample) and log
col_ok <- colSums(is.na(e1)) == 0
if (!all(col_ok)) {
  note <- paste0(note, "all-NA columns dropped: ", sum(!col_ok), "; ")
  e1 <- e1[, col_ok, drop = FALSE]
}
sc <- normalize_scale(e1)
e2 <- sc$expr
keep_na <- rowSums(is.na(e2)) == 0
if (!all(keep_na)) {
  note <- paste0(note, "NA-containing gene rows dropped: ", sum(!keep_na), "; ")
  e2 <- e2[keep_na, , drop = FALSE]
  if (!is.null(sc$counts)) sc$counts <- sc$counts[keep_na, , drop = FALSE]
}
  list(expr = e2, counts = sc$counts, scale = sc$scale, hdr = hdr,
       col_gsm = col_gsm, note = note, gene_hit_rate = hit_rate,
       n_genes_raw = n_genes_raw, n_genes_mapped = nrow(e0))
}

# ---- build_from_manifest ------------------------------------------------------
# Build a study object from the manifest: read_geo_study (keep = manifest GSMs) → donor prefix → column/row alignment → validate
build_from_manifest <- function(gse, manifest_rows, layer, cluster_id, stability = FALSE,
                                min_genes = 1000L, min_samples = 1L, platform_hint = NULL) {
  check(all(c("gsm", "source", "donor", "arm", "passage", "medium") %in% names(manifest_rows)),
        gse, ": manifest lacks required columns")
  rd <- read_geo_study(gse, gsm_keep = manifest_rows$gsm, platform_hint = platform_hint)
  m <- manifest_rows
  km <- match(colnames(rd$expr), m$gsm)
  check(!anyNA(km), gse, ": expr columns do not match manifest GSMs")
  m <- m[km, , drop = FALSE]
  m$donor <- paste(m$source, m$donor, sep = ":")          # S3 rule: donor prefixed by source
  check(nrow(rd$expr) >= min_genes, gse, ": genes after cleaning ", nrow(rd$expr), " < ", min_genes)
  check(nrow(m) >= min_samples, gse, ": too few samples")
  s <- list(gse = gse, cluster_id = cluster_id, layer = layer, stability = stability,
            full_genome = nrow(rd$expr) >= CFG$genome_min_genes, scale = rd$scale,
            expr = rd$expr, counts = rd$counts, meta = m)
  validate_study(s)
  s
}
