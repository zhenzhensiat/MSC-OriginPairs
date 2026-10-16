# =============================================================================
# 11b_build_passage.R —— dedicated build of two supplementary passage series
# (GSE139073 / GSE309666) from raw files.
# GSE139073 (BM, ENSG counts): two GPL subseries; S numbers -> GSM <-> donor
#   <-> passage <-> stimulation; STIM arms excluded (standard untreated arms only).
# GSE309666 (AD, salmon gene counts): counts columns "<donor>_p<N>" matched to
#   donor+passage characteristics.
# Outputs: data/studies/{GSE139073,GSE309666}.rds (with counts slot); appends to
#   coverage/sample-metadata tables.
# =============================================================================
source("scripts/00_setup.R")
need_pkgs(c("org.Hs.eg.db", "AnnotationDbi", "edgeR"))
library(org.Hs.eg.db)

parse_series <- function(path) {   # extract title/gsm/characteristic rows
  ln <- readLines(gzfile(path))
  q <- function(tag) {
    hit <- ln[startsWith(ln, tag)]
    if (length(hit) == 0) return(NULL)
    lapply(hit, function(h) trimws(gsub('"', "", strsplit(h, "\t", fixed = TRUE)[[1]][-1])))
  }
  list(title = q("!Sample_title")[[1]], gsm = q("!Sample_geo_accession")[[1]],
       chars = q("!Sample_characteristics_ch1"))
}
char_lookup <- function(chars, key) {           # "donor: Donor D" -> value by key prefix
  vapply(seq_along(chars[[1]]), function(i) {
    hit <- vapply(chars, function(row) row[i], character(1))
    v <- hit[grepl(paste0("^", key, ":"), hit)]
    if (length(v) == 0) return(NA_character_)
    trimws(sub(paste0("^", key, ":"), "", v[1]))
  }, character(1))
}

run_step("11b_build_passage", function() {
  alloc <- read_tsv("data/layer_allocation_v2.tsv")
  geo <- read_tsv("data/census_geo_v2.tsv")
  dir.create(proj_path("data", "studies"), showWarnings = FALSE)

  # ---- GSE139073: ENSG counts + two GPL series matrices ---------------------------
  g <- "GSE139073"
  ss <- lapply(list.files(proj_path("data", "raw", g), "_series_matrix.*gz$", full.names = TRUE),
               parse_series)
  ttl <- unlist(lapply(ss, `[[`, "title"), use.names = FALSE)
  gsm <- unlist(lapply(ss, `[[`, "gsm"), use.names = FALSE)
  chr <- do.call(c, lapply(ss, function(s) s$chars))
  # characteristic rows aligned by position after per-file concatenation: per-file processing then rbind is safer
  met139 <- do.call(rbind, lapply(ss, function(s) data.frame(
    title = s$title, gsm = s$gsm, donor = char_lookup(s$chars, "donor"),
    passage = char_lookup(s$chars, "passage"), stim = char_lookup(s$chars, "stimulation"),
    stringsAsFactors = FALSE)))
  met139$sid <- sub(".*, *(S[0-9]+).*$", "\\1", met139$title)     # E18: full extraction via sub
  check(all(grepl("^S[0-9]+$", met139$sid)), g, ": S number extraction failed")
  check(!anyDuplicated(met139$sid), g, ": duplicate S numbers")

  raw <- read.delim(gzfile(proj_path("data", "raw", g, "GSE139073_raw_counts.txt.gz")),
                    check.names = FALSE, stringsAsFactors = FALSE)
  num <- as.matrix(raw); storage.mode(num) <- "numeric"
  ens <- sub("\\..*$", "", rownames(num))
  emap <- AnnotationDbi::select(org.Hs.eg.db, keys = ens, columns = "SYMBOL", keytype = "ENSEMBL")
  keep <- !duplicated(emap$ENSEMBL) & !is.na(emap$SYMBOL)
  sym <- ifelse(ens %in% emap$ENSEMBL[keep], emap$SYMBOL[match(ens, emap$ENSEMBL[keep])], NA)
  num <- num[!is.na(sym), , drop = FALSE]; rownames(num) <- sym[!is.na(sym)]
  if (anyDuplicated(rownames(num))) num <- rowsum(num, group = rownames(num))
  colnames(num) <- sub('"', "", colnames(num))
  m139 <- met139[match(colnames(num), met139$sid), , drop = FALSE]
  check(!anyNA(m139$gsm), g, ": counts columns do not align with S numbers")
  keep_arm <- m139$stim %in% c("REST", NA_character_)            # STIM arm excluded
  num_i <- num[, keep_arm, drop = FALSE]; m_i <- m139[keep_arm, ]
  colnames(num_i) <- m_i$gsm
  e <- edgeR::cpm(num_i, log = TRUE, prior.count = 0.5)
  s139 <- list(
    gse = g, cluster_id = alloc$cluster_id[alloc$gse == g], layer = alloc$layer[alloc$gse == g],
    stability = TRUE, full_genome = TRUE, scale = "logCPM", expr = e, counts = num_i,
    meta = data.frame(gsm = m_i$gsm, source = "BM", donor = m_i$donor,
                      passage = as.numeric(sub("P", "", m_i$passage)), medium = NA_character_,
                      age = NA_real_, sex = NA_character_, site = NA_character_,
                      arm = ifelse(is.na(m_i$stim) | m_i$stim == "REST", "include", "stim"),
                      stringsAsFactors = FALSE, row.names = NULL))
  validate_study(s139)
  saveRDS(s139, proj_path("data", "studies", paste0(g, ".rds")))
  log_msg(g, ": ", nrow(e), " genes x ", ncol(e), " samples (REST arm; STIM arm excluded: ",
          sum(!keep_arm), ")")

  # ---- GSE309666: salmon gene counts, column names <donor>_p<N> --------------------------
  g <- "GSE309666"
  ss <- parse_series(proj_path("data", "raw", g, "GSE309666_series_matrix.txt.gz"))
  donor <- char_lookup(ss$chars, "donor"); pas <- char_lookup(ss$chars, "passage")
  tis <- char_lookup(ss$chars, "tissue")
  key <- paste0(sub("Donor ", "", donor, fixed = TRUE), "_", pas)
  raw <- read.delim(gzfile(proj_path("data", "raw", g, "GSE309666_salmon_gene_counts.tsv.gz")),
                    check.names = FALSE, stringsAsFactors = FALSE)
  sym <- raw$gene_name
  num <- as.matrix(raw[, grep("^[A-Z]_p[0-9]+$", names(raw)), drop = FALSE])
  storage.mode(num) <- "numeric"
  sym_ok <- unique(na.omit(AnnotationDbi::select(org.Hs.eg.db, keys = keys(org.Hs.eg.db),
                              columns = "SYMBOL")$SYMBOL))
  ok <- !is.na(sym) & nzchar(sym) & sym %in% sym_ok   # M2: HGNC whitelist filtering
  # gene_id and gene_name are parallel: align num by raw row order (num filtered by column, not by row)
  stopifnot(nrow(num) == length(sym))
  num <- num[ok, , drop = FALSE]; rownames(num) <- sym[ok]
  if (anyDuplicated(rownames(num))) num <- rowsum(num, group = rownames(num))
  m309 <- data.frame(key = colnames(num), gsm = ss$gsm[match(colnames(num), key)],
                     donor = donor[match(colnames(num), key)],
                     passage = as.numeric(sub("p", "", pas[match(colnames(num), key)])),
                     stringsAsFactors = FALSE)
  check(!anyNA(m309$gsm), g, ": counts column names do not align with donor_passage keys")
  colnames(num) <- m309$gsm
  e <- edgeR::cpm(num, log = TRUE, prior.count = 0.5)
  s309 <- list(
    gse = g, cluster_id = alloc$cluster_id[alloc$gse == g], layer = alloc$layer[alloc$gse == g],
    stability = TRUE, full_genome = TRUE, scale = "logCPM", expr = e, counts = num,
    meta = data.frame(gsm = m309$gsm, source = "AD", donor = m309$donor,
                      passage = m309$passage, medium = NA_character_, age = NA_real_,
                      sex = NA_character_, site = NA_character_, arm = "include",
                      stringsAsFactors = FALSE, row.names = NULL))
  validate_study(s309)
  saveRDS(s309, proj_path("data", "studies", paste0(g, ".rds")))
  log_msg(g, ": ", nrow(e), " genes x ", ncol(e), " samples (5 donors p1->p10)")

  # ---- append the two audit tables ------------------------------------------------------------
  for (s in list(s139, s309)) {
    add <- data.frame(gse = s$gse, n_genes_raw = NA_integer_, n_genes_mapped = nrow(s$expr),
                      annot_rate = NA_real_, n_samples = ncol(s$expr),
                      n_include = sum(s$meta$arm == "include"),
                      genome_wide = s$full_genome, scale = s$scale,
                      notes = "11b primary build", stringsAsFactors = FALSE)
    aud <- read_tsv("data/coverage_audit_v2.tsv")
    aud <- rbind(aud, add)
    write_tsv(aud, "data/coverage_audit_v2.tsv")
    mm <- s$meta; mm$study <- s$gse; mm$cluster_id <- s$cluster_id; mm$layer <- s$layer
    sm <- read_tsv("data/sample_metadata_v2.tsv")
    common <- intersect(names(sm), names(mm))
    write_tsv(rbind(sm[, common], mm[, common]), "data/sample_metadata_v2.tsv")
  }
  invisible(TRUE)
})
