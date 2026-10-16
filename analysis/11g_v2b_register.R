# =============================================================================
# 11g_v2b_register.R —— register the four network-series studies into the three
# result tables (completion of the extended layer).
# Updates: data/layer_allocation_v2.tsv (+4 rows, layer D),
#   data/coverage_audit_v2.tsv (+4 rows), data/sample_metadata_v2.tsv (+20 rows)
#   -> all rewritten.
# Annotation-rate definition identical to script 11.
# =============================================================================
source("scripts/00_setup.R")
need_pkgs(c("org.Hs.eg.db", "AnnotationDbi"))
library(org.Hs.eg.db)

run_step("11g_v2b_register", function() {
  new_gse <- c("GSE37065", "GSE9520", "GSE189895", "GSE194298")

  # ---- 1. layer_allocation: in-place C->D layer change (former C candidate rows) + append missing rows ---------------
  alloc <- read_tsv("data/layer_allocation_v2.tsv")
  cov <- which(alloc$gse %in% new_gse)
  old_layers <- paste0(alloc$gse[cov], ":", alloc$layer[cov], collapse = ", ")
  alloc$cluster_id[cov] <- "G02"
  alloc$layer[cov] <- "D"
  miss <- setdiff(new_gse, alloc$gse)
  if (length(miss)) {
    add <- data.frame(cluster_id = "G02", gse = miss, design = "single_source",
                      data_type = ifelse(miss %in% c("GSE189895", "GSE194298"), "rnaseq", "array"),
                      full_genome = TRUE, in_v1_L1 = FALSE, in_v1_L2 = FALSE,
                      stability = FALSE, layer = "D", stringsAsFactors = FALSE)
    alloc <- rbind(alloc, add[, names(alloc)])
  }
  write_tsv(alloc, "data/layer_allocation_v2.tsv")
  log_msg("allocation: in-place layer change for ", length(cov), " rows (was ", old_layers, "); appended ", length(miss), " rows")

  # ---- 2. coverage_audit: append 4 rows (annotation rate recomputed on the spot as in script 11; idempotent)-------------
  aud <- read_tsv("data/coverage_audit_v2.tsv")
  aud <- aud[!(aud$gse %in% new_gse), , drop = FALSE]
  eg <- AnnotationDbi::select(org.Hs.eg.db, keys = AnnotationDbi::keys(org.Hs.eg.db),
                              columns = "SYMBOL")
  sym_ok <- unique(na.omit(eg$SYMBOL))
  probe_rows <- function(sm_path) {
    ln <- readLines(gzfile(sm_path, "rt"), warn = FALSE)
    ib <- grep("^\"?ID_REF\"?", ln)
    if (length(ib) != 1) ib <- grep("series_matrix_table_begin", ln) + 1
    ie <- grep("series_matrix_table_end", ln)
    list(ids = gsub('^"|"$', '', sub("\t.*$", "", ln[(ib + 1):(ie - 1)])))
  }
  hit_rate <- function(g, ids) {
    if (g %in% c("GSE37065", "GSE9520")) {
      plat <- if (g == "GSE37065") "GPL6244.annot.gz" else "GPL570.annot.gz"
      con <- gzfile(proj_path("data", "raw", "platforms", plat), "rt")
      ln <- readLines(con, warn = FALSE); close(con)
      ih <- grep("^ID\t", ln)
      key <- sub("\t.*$", "", ln[(ih + 1):length(ln)])
      sym <- sub("\t.*$", "", sub("^[^\t]*\t[^\t]*\t", "", ln[(ih + 1):length(ln)]))
      sym <- trimws(sub("\\s*///.*$", "", sym))
      m <- setNames(sym, key); m[!nzchar(m) | m == "---"] <- NA
      mean(!is.na(m[ids]))
    } else {
      mm <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(ids),
                                  columns = "SYMBOL", keytype = "ENSEMBL")
      mm <- mm[!is.na(mm$SYMBOL), ]
      mean(unique(ids) %in% mm$ENSEMBL)
    }
  }
  for (g in new_gse) {
    ids <- if (g %in% c("GSE37065", "GSE9520")) {
      probe_rows(proj_path("data", "raw", g, paste0(g, "_series_matrix.txt.gz")))$ids
    } else {
      cfile <- list.files(proj_path("data", "raw", g), "csv.gz$", full.names = TRUE)[1]
      con <- gzfile(cfile, "rt"); hl <- readLines(con, warn = FALSE); close(con)
      sub(',.*$', '', gsub('"', '', hl[-1]))
    }
    s <- readRDS(proj_path("data", "studies", paste0(g, ".rds")))
    hr <- hit_rate(g, ids)
    aud <- rbind(aud, data.frame(
      gse = g, n_genes_raw = length(ids), n_genes_mapped = nrow(s$expr),
      annot_rate = round(hr, 4), n_samples = ncol(s$expr),
      n_include = sum(s$meta$arm == "include"), genome_wide = s$full_genome,
      scale = s$scale, notes = "v2b build", stringsAsFactors = FALSE))
    log_msg(g, ": coverage_audit registered (annotation rate ", round(hr, 4), ")")
  }
  check(sum(aud$gse %in% new_gse) == 4, "unexpected coverage_audit appended row count")
  write_tsv(aud, "data/coverage_audit_v2.tsv")

  # ---- 3. sample_metadata: append sample rows for the 4 studies (projected onto the existing column set; idempotent)--------------
  sm <- read_tsv("data/sample_metadata_v2.tsv")
  sm <- sm[!(sm$study %in% new_gse), , drop = FALSE]
  need_cols <- c("gsm", "source", "donor", "passage", "medium", "age", "sex", "site", "arm")
  add_rows <- do.call(rbind, lapply(new_gse, function(g) {
    s <- readRDS(proj_path("data", "studies", paste0(g, ".rds")))
    mm <- s$meta
    check(all(need_cols %in% names(mm)), g, ": meta is missing columns")
    data.frame(mm[, need_cols], study = g, cluster_id = s$cluster_id, layer = s$layer,
               stringsAsFactors = FALSE)
  }))
  sm2 <- rbind(sm, add_rows)
  write_tsv(sm2, "data/sample_metadata_v2.tsv")
  log_msg("sample_metadata: ", nrow(sm), " -> ", nrow(sm2), " rows")

  # ---- 4. acceptance: re-read the three tables ----------------------------------------------
  a2 <- read_tsv("data/layer_allocation_v2.tsv")
  check(sum(a2$gse %in% new_gse) == 4 && all(a2$layer[a2$gse %in% new_gse] == "D"), "allocation verification failed")
  log_msg("11g done: three tables registered and verified")
})
