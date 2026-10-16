# =============================================================================
# 11f_build_series.R —— build the four network-series studies into the
# development layer (extended layer).
# Basis: sample selection = standard culture,
# untreated, early passage (<=P5):
#   GSE37065 (GPL6244, BM): early arm 3 samples (P3/4/5; late P10 excluded)
#   GSE9520  (GPL570, BM): all 6 samples (P2; CCM + 17% FBS standard medium)
#   GSE189895 (RNA-seq, BM): primary BMAC arm only, 3 samples (iPSC/MIM arms excluded)
#   GSE194298 (RNA-seq, BM): FBS standard arm, 8 samples (HPL arm excluded; p1, p4/p5 <=P5)
# Rules: probe-annotation mapping + HGNC normalisation + dedup (same as 11/11e);
#   integer counts -> counts slot + cpm(log); non-integer RNA-seq -> log2(x+1)
#   flagged log2TPM ("voom not applicable").
# Outputs: data/studies/{GSE37065,GSE9520,GSE189895,GSE194298}.rds
# =============================================================================
source("scripts/00_setup.R")
need_pkgs(c("org.Hs.eg.db", "AnnotationDbi", "edgeR"))
library(org.Hs.eg.db)

run_step("11f_build_wagner", function() {

  # ---- shared utilities (same rules as 11/11e)---------------------------------------------
  load_annot_gz <- function(path, min_rows) {
    con <- gzfile(path, "rt"); ln <- readLines(con, warn = FALSE); close(con)
    ih <- grep("^ID\t", ln)
    check(length(ih) == 1, basename(path), ": unique header row not found")
    dat <- ln[(ih + 1):length(ln)]; dat <- dat[nzchar(dat)]
    check(length(dat) >= min_rows, basename(path), ": data rows ", length(dat), " < ", min_rows)
    df <- read.delim(text = paste(c(ln[ih], dat), collapse = "\n"), header = TRUE,
                     quote = "", check.names = FALSE, stringsAsFactors = FALSE)
    m <- setNames(df[["Gene symbol"]], df[["ID"]])
    m <- trimws(sub("\\s*///.*$", "", m))
    m[!nzchar(m) | m == "---"] <- NA_character_
    m
  }

  eg <- AnnotationDbi::select(org.Hs.eg.db, keys = AnnotationDbi::keys(org.Hs.eg.db),
                              columns = c("SYMBOL", "ALIAS"))
  sym_ok <- unique(na.omit(eg$SYMBOL))
  al <- eg[!is.na(eg$ALIAS) & !is.na(eg$SYMBOL), c("ALIAS", "SYMBOL")]
  al_n <- table(al$ALIAS); uniq_alias <- names(al_n[al_n == 1L])
  alias_map <- setNames(al$SYMBOL[match(uniq_alias, al$ALIAS)], uniq_alias)
  map_genes <- function(genes) {
    out <- ifelse(genes %in% sym_ok, genes,
           ifelse(genes %in% names(alias_map), alias_map[genes], NA_character_))
    out[is.na(out)] <- NA_character_; out
  }
  collapse_dup <- function(x) {
    if (!anyDuplicated(rownames(x))) return(x)
    mu <- rowMeans(x); iqr <- apply(x, 1, stats::IQR)
    ord <- order(mu, iqr, decreasing = TRUE)
    x[ord[!duplicated(rownames(x)[ord])], , drop = FALSE]
  }
  parse_medium <- function(txt) {
    t <- tolower(txt)
    ifelse(grepl("platelet lysate|\\bhpl\\b", t), "hPL",
    ifelse(grepl("fetal bovine|\\bfbs\\b", t, fixed = FALSE) & !grepl("platelet", t), "FBS",
    ifelse(grepl("human serum|\\bhs\\b", t) & !grepl("fetal", t), "HS",
    ifelse(grepl("serum.?free|xeno.?free", t), "xeno_free", NA_character_))))
  }
  y_genes <- c("DDX3Y", "RPS4Y1", "KDM5D", "EIF1AY", "USP9Y", "ZFY", "UTY")

  field_row <- function(ln, field) {
    i <- grep(paste0("^!", field, "\t"), ln)[1]
    check(!is.na(i), "missing field ", field)
    gsub('^"|"$', '', strsplit(ln[i], "\t", fixed = TRUE)[[1]][-1])
  }
  char_all <- function(ln) {   # merge all characteristics rows by column into text
    idx <- grep("^!Sample_characteristics_ch1\t", ln)
    rows <- lapply(idx, function(i) gsub('^"|"$', '', strsplit(ln[i], "\t", fixed = TRUE)[[1]][-1]))
    do.call(paste, c(rows, list(sep = " | ")))
  }

  load_sm <- function(path) {
    ln <- readLines(gzfile(path, "rt"), warn = FALSE)
    ib <- grep("^\"?ID_REF\"?", ln)
    if (length(ib) != 1) ib <- grep("series_matrix_table_begin", ln) + 1
    ie <- grep("series_matrix_table_end", ln)
    check(length(ib) == 1 && length(ie) == 1 && ie > ib, basename(path), ": data table interval abnormal")
    df <- read.delim(text = paste(ln[ib:(ie - 1)], collapse = "\n"), header = TRUE,
                     quote = "\"", check.names = FALSE, stringsAsFactors = FALSE)
    ids <- as.character(df[[1]])
    m <- as.matrix(df[, -1, drop = FALSE]); mode(m) <- "numeric"
    na_all <- rowSums(is.na(m)) == ncol(m)
    if (any(na_all)) { log_msg(basename(path), ": all-missing probe rows dropped: ", sum(na_all), " rows"); m <- m[!na_all, , drop = FALSE]; ids <- ids[!na_all] }
    check(!anyNA(m), basename(path), ": data contain partial NA, manual review required")
    list(ids = ids, mat = m)
  }

  probe_build <- function(sm, map) {   # array: probe -> symbol -> normalise -> dedup -> scale
    sym <- map[as.character(sm$ids)]
    m <- sm$mat[!is.na(sym), , drop = FALSE]
    rownames(m) <- unname(sym[!is.na(sym)])
    g2 <- map_genes(rownames(m))
    m <- m[!is.na(g2), , drop = FALSE]; rownames(m) <- g2[!is.na(g2)]
    m <- collapse_dup(m)
    if (max(m) >= CFG$log_scale_max) { check(all(m >= 0), "log needed but negative values present"); m <- log2(m + 1) }
    check(max(m) < CFG$log_scale_max && !anyNA(m), "abnormal after scale transformation")
    m
  }
  ens_build <- function(cnt, counts_slot) {   # RNA-seq: ENSG -> symbol -> normalise -> dedup; integers -> cpm(log)
    ens <- rownames(cnt)
    mm <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(ens),
                                columns = "SYMBOL", keytype = "ENSEMBL")
    mm <- mm[!is.na(mm$SYMBOL), ]; mm <- mm[!duplicated(mm$ENSEMBL), ]
    e2s <- setNames(mm$SYMBOL, mm$ENSEMBL)
    sym <- e2s[ens]
    cnt <- cnt[!is.na(sym), , drop = FALSE]
    rownames(cnt) <- unname(sym[!is.na(sym)])
    g2 <- map_genes(rownames(cnt))
    cnt <- cnt[!is.na(g2), , drop = FALSE]; rownames(cnt) <- g2[!is.na(g2)]
    cnt <- collapse_dup(cnt)
    if (counts_slot) {
      check(mean(abs(cnt - round(cnt)) < 1e-8) > 0.99, "integer counts check failed")
      list(expr = edgeR::cpm(cnt, log = TRUE, prior.count = 0.5), counts = cnt, scale = "logCPM")
    } else {
      e <- if (max(cnt) >= CFG$log_scale_max) log2(cnt + 1) else cnt
      list(expr = e, counts = NULL, scale = "log2TPM")
    }
  }
  finish <- function(g, e, counts_m, scale, meta, note) {
    check(identical(colnames(e), meta$gsm), g, ": expr columns do not match meta")
    yg <- rownames(e)[rownames(e) %in% y_genes]
    meta$ygene_mean <- if (length(yg)) colMeans(e[yg, , drop = FALSE]) else NA_real_
    meta$xist <- if ("XIST" %in% rownames(e)) e["XIST", ] else NA_real_
    s <- list(gse = g, cluster_id = "G02", layer = "D", stability = FALSE,
              full_genome = nrow(e) >= CFG$genome_min_genes, scale = scale,
              expr = e, counts = counts_m, meta = meta)
    validate_study(s)
    saveRDS(s, proj_path("data", "studies", paste0(g, ".rds")))
    log_msg(g, " (v2b): ", nrow(e), " genes / ", ncol(e), " samples; scale=", scale,
            "; full_genome=", s$full_genome, "; ", note)
    invisible(s)
  }
  mk_meta <- function(gsm, donor_raw, passage, medium) {
    data.frame(donor_shared = FALSE, gsm = gsm, source = "BM",
               donor = paste("BM", donor_raw, sep = ":"), donor_raw = donor_raw,
               passage = passage, medium = medium, age = NA_real_, sex = NA_character_,
               site = NA_character_, arm = "include", stringsAsFactors = FALSE)
  }

  # ---- 1. GSE37065 (GPL6244; early arm)-----------------------------------------
  ln <- readLines(gzfile(proj_path("data", "raw", "GSE37065", "GSE37065_series_matrix.txt.gz")), warn = FALSE)
  gsm <- field_row(ln, "Sample_geo_accession"); title <- field_row(ln, "Sample_title")
  txt <- char_all(ln)
  keep <- grepl("early", title, ignore.case = TRUE)
  check(sum(keep) == 3, "GSE37065 early arm: unexpected sample count: ", sum(keep))
  passa <- as.numeric(sub(".*passage:\\s*([0-9]+).*", "\\1", txt))
  check(all(passa[keep] <= 5), "GSE37065 early arm passage exceeds 5")
  sm <- load_sm(proj_path("data", "raw", "GSE37065", "GSE37065_series_matrix.txt.gz"))
  check(all(gsm[keep] %in% colnames(sm$mat)), "GSE37065: matrix is missing selected samples")
  sm$mat <- sm$mat[, match(gsm[keep], colnames(sm$mat)), drop = FALSE]
  map6244 <- load_annot_gz(proj_path("data", "raw", "platforms", "GPL6244.annot.gz"), 33000)
  e <- probe_build(sm, map6244)
  donor_raw <- sub("^BM (MSC[0-9]+)_.*$", "\\1", title[keep])
  finish("GSE37065", e, NULL, "array", mk_meta(gsm[keep], donor_raw, passa[keep], parse_medium(txt[keep])),
         "early arm (late P10 excluded)")

  # ---- 2. GSE9520 (GPL570; all samples)---------------------------------------------
  ln <- readLines(gzfile(proj_path("data", "raw", "GSE9520", "GSE9520_series_matrix.txt.gz")), warn = FALSE)
  gsm <- field_row(ln, "Sample_geo_accession"); title <- field_row(ln, "Sample_title")
  txt <- char_all(ln)
  check(length(gsm) == 6 && all(grepl("Passage2", title)), "GSE9520 sample structure abnormal")
  sm <- load_sm(proj_path("data", "raw", "GSE9520", "GSE9520_series_matrix.txt.gz"))
  sm$mat <- sm$mat[, match(gsm, colnames(sm$mat)), drop = FALSE]
  map570 <- load_annot_gz(proj_path("data", "raw", "platforms", "GPL570.annot.gz"), 54600)
  e <- probe_build(sm, map570)
  donor_raw <- sub("^hMSC_?Donor([0-9A-Za-z]+)_Passage.*$", "\\1", title)
  finish("GSE9520", e, NULL, "array", mk_meta(gsm, donor_raw, rep(2, 6), parse_medium(txt)),
         "all 6 samples (P2; Day2/7 same-donor registered)")

  # ---- 3. GSE189895 (RNA-seq; primary BMAC arm)------------------------------------
  ln <- readLines(gzfile(proj_path("data", "raw", "GSE189895", "GSE189895_series_matrix.txt.gz")), n = 300)
  gsm <- field_row(ln, "Sample_geo_accession"); title <- field_row(ln, "Sample_title")
  desc <- field_row(ln, "Sample_description")
  keep <- grepl("^MSPCs BMAC", title)
  check(sum(keep) == 3, "GSE189895 BMAC arm: unexpected sample count: ", sum(keep))
  cnt <- as.matrix(read.csv(gzfile(proj_path("data", "raw", "GSE189895",
                "GSE189895_Normalized_counts_RNAseq_iPSC_diff.csv.gz")),
                check.names = FALSE, row.names = 1))
  cn_num <- sub("^X([0-9]+)$", "\\1", colnames(cnt))
  pick <- match(desc[keep], cn_num)
  check(!anyNA(pick), "GSE189895: sample ID mapping failed")
  csub <- cnt[, pick, drop = FALSE]; colnames(csub) <- gsm[keep]
  nb <- ens_build(csub, counts_slot = TRUE)
  donor_raw <- rep(sub("^MSPCs BMAC ([A-Z]+) .*", "\\1", title[keep])[1], sum(keep))
  finish("GSE189895", nb$expr, nb$counts, nb$scale,
         mk_meta(gsm[keep], donor_raw, rep(NA_real_, 3), rep(NA_character_, 3)),
         "primary BMAC arm only; iPSC/MIM and cord-blood arms excluded (project policy); census notes BMAC A/B but only A measured")

  # ---- 4. GSE194298 (RNA-seq; FBS arm)------------------------------------------
  ln <- readLines(gzfile(proj_path("data", "raw", "GSE194298", "GSE194298_series_matrix.txt.gz")), n = 300)
  gsm <- field_row(ln, "Sample_geo_accession"); title <- field_row(ln, "Sample_title")
  desc <- field_row(ln, "Sample_description")
  all_rows <- grep("^!Sample_characteristics_ch1\t", ln)
  grab_row <- function(prefix) {
    for (i in all_rows) {
      v <- gsub('^"|"$', '', strsplit(ln[i], "\t", fixed = TRUE)[[1]][-1])
      if (any(grepl(paste0("^", prefix, ":"), v))) return(v)
    }
    NULL
  }
  trt_row <- grab_row("treatment"); passa_row <- grab_row("passage"); donor_row <- grab_row("donor")
  check(!is.null(trt_row) && !is.null(passa_row) && !is.null(donor_row), "GSE194298: field row missing")
  keep <- grepl("^treatment:\\s*FBS", trt_row)
  check(sum(keep) == 8, "GSE194298 FBS arm: unexpected sample count: ", sum(keep))
  passa <- as.numeric(sub("^passage:\\s*p([0-9]+).*$", "\\1", passa_row))
  donor_raw <- sub("^donor:\\s*", "", donor_row)
  cnt <- as.matrix(read.csv(gzfile(proj_path("data", "raw", "GSE194298",
                "GSE194298_Norm_counts_RNAseq_hPL_vs_FBS.csv.gz")),
                check.names = FALSE, row.names = 1))
  cn_num <- sub("^X([0-9]+)\\..*$", "\\1", colnames(cnt))
  pick <- match(desc[keep], cn_num)
  check(!anyNA(pick), "GSE194298: sample ID mapping failed")
  csub <- cnt[, pick, drop = FALSE]; colnames(csub) <- gsm[keep]
  nb <- ens_build(csub, counts_slot = FALSE)
  finish("GSE194298", nb$expr, nb$counts, nb$scale,
         mk_meta(gsm[keep], donor_raw[keep], passa[keep], rep("FBS", sum(keep))),
         "FBS standard arm (HPL arm excluded); author-normalised counts (non-integer), voom not applicable")

  log_msg("v2b build complete: validate_study passed for all 4 studies")
})
