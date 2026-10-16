# =============================================================================
# 44a_d4_coverage.R —— platform coverage table for the validation cohort
# (input to the final model).
# Discipline: only platform annotation files and gene names are read; no
# validation-cohort expression file is opened.
# Inputs: validation study platform list; platform annotations under
#   data/raw/platforms/ (annot.gz / soft.txt / Bioc db)
# Output: data/c_platform_coverage_v2.tsv - per-platform gene-list size,
#   coverage of the development training intersection, and the all-platform
#   intersection (candidate "analysis universe").
# Note: sequencing platforms count as genome-wide without restriction;
#   NanoString is a small custom panel (its coverage decides the fixed
#   "<4000 -> drop worst-covered platform" rule).
# =============================================================================
source("scripts/00_setup.R")
need_pkgs(c("org.Hs.eg.db", "AnnotationDbi", "huex10sttranscriptcluster.db",
            "hugene20sttranscriptcluster.db"))
library(org.Hs.eg.db)

run_step("44a_d4_coverage", function() {

  # ---- 0. HGNC normalization (same rules as 11/11e)------------------------------------------
  eg <- AnnotationDbi::select(org.Hs.eg.db, keys = AnnotationDbi::keys(org.Hs.eg.db),
                              columns = c("SYMBOL", "ALIAS"))
  sym_ok <- unique(na.omit(eg$SYMBOL))
  al <- eg[!is.na(eg$ALIAS) & !is.na(eg$SYMBOL), c("ALIAS", "SYMBOL")]
  al_n <- table(al$ALIAS); ua <- names(al_n[al_n == 1L])
  alias_map <- setNames(al$SYMBOL[match(ua, al$ALIAS)], ua)
  norm <- function(x) {
    x <- trimws(x)
    x <- sub("\\s*///.*$", "", x)                 # first value only
    out <- ifelse(x %in% sym_ok, x,
           ifelse(x %in% names(alias_map), alias_map[x], NA_character_))
    unique(na.omit(out))
  }

  # ---- 1. gene lists per C platform -------------------------------------------------------
  load_annot_gz <- function(path) {
    con <- gzfile(path, "rt"); ln <- readLines(con, warn = FALSE); close(con)
    ih <- grep("^ID\t", ln)
    check(length(ih) == 1, basename(path), ": abnormal header row")
    dat <- ln[(ih + 1):length(ln)]; dat <- dat[nzchar(dat)]
    df <- read.delim(text = paste(c(ln[ih], dat), collapse = "\n"), header = TRUE,
                     quote = "", check.names = FALSE, stringsAsFactors = FALSE)
    norm(df[["Gene symbol"]])
  }
  load_soft_col <- function(path, col, min_rows = 100) {
    ln <- readLines(path, warn = FALSE)
    ib <- grep("^!platform_table_begin", ln); ie <- grep("^!platform_table_end", ln)
    check(length(ib) == 1 && length(ie) == 1 && ie > ib, basename(path), ": abnormal table markers (suspected E7)")
    df <- read.delim(text = paste(ln[ib:ie], collapse = "\n"), header = TRUE,
                     sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE,
                     skip = 1)
    check(nrow(df) >= min_rows, basename(path), ": row count ", nrow(df), " too low")
    check(col %in% names(df), basename(path), ": missing column ", col)
    if (col == "gene_assignment") {            # Affy ST format: 'locus // symbol // ... /// next record'
      parts <- strsplit(df[[col]], " /// ", fixed = TRUE)
      syms <- unlist(lapply(parts, function(pp) {
        vapply(strsplit(pp, " // ", fixed = TRUE), function(z) {
          if (length(z) >= 2 && !grepl("^---", z[2])) trimws(z[2]) else NA_character_
        }, character(1))
      }))
      norm(syms)
    } else norm(df[[col]])
  }
  entrez_soft <- function(path, col) {         # Entrez column in the SOFT table -> symbol (GPL9193)
    ln <- readLines(path, warn = FALSE)
    ib <- grep("^!platform_table_begin", ln); ie <- grep("^!platform_table_end", ln)
    check(length(ib) == 1 && length(ie) == 1 && ie > ib, basename(path), ": abnormal table markers")
    df <- read.delim(text = paste(ln[ib:ie], collapse = "\n"), header = TRUE, sep = "\t",
                     quote = "", check.names = FALSE, stringsAsFactors = FALSE, skip = 1)
    check(col %in% names(df), basename(path), ": missing column ", col)
    ek <- AnnotationDbi::select(org.Hs.eg.db, keys = as.character(unique(na.omit(df[[col]]))),
                                columns = "SYMBOL", keytype = "ENTREZID")
    norm(unique(na.omit(ek$SYMBOL)))
  }
  db_genes <- function(db, keytype = "PROBEID") {
    ks <- AnnotationDbi::keys(db)
    m <- AnnotationDbi::select(db, keys = ks, columns = "SYMBOL", keytype = keytype)
    norm(unique(na.omit(m$SYMBOL)))
  }

  plats <- list(
    GPL570  = load_annot_gz(proj_path("data", "raw", "platforms", "GPL570.annot.gz")),
    GPL571  = load_annot_gz(proj_path("data", "raw", "platforms", "GPL571.annot.gz")),
    GPL1708 = load_annot_gz(proj_path("data", "raw", "platforms", "GPL1708.annot.gz")),
    GPL6884 = load_annot_gz(proj_path("data", "raw", "platforms", "GPL6884.annot.gz")),
    GPL10558 = load_annot_gz(proj_path("data", "raw", "platforms", "GPL10558.annot.gz")),
    GPL5188 = db_genes(huex10sttranscriptcluster.db::huex10sttranscriptcluster.db),
    GPL16686 = db_genes(hugene20sttranscriptcluster.db::hugene20sttranscriptcluster.db),
    GPL7088 = load_soft_col(proj_path("data", "raw", "platforms", "GPL7088.soft.txt"), "ORF"),
    GPL9193 = entrez_soft(proj_path("data", "raw", "platforms", "GPL9193.soft.txt"), "Locus Gene ID"),
    GPL17586 = load_soft_col(proj_path("data", "raw", "platforms", "GPL17586.soft.txt"), "gene_assignment"),
    GPL28102 = load_soft_col(proj_path("data", "raw", "platforms", "GPL28102.soft.txt"), "GeneSymbol"),
    GPL14550 = load_soft_col(proj_path("data", "raw", "platforms", "GPL14550.soft.txt"), "GENE_SYMBOL"),
    GPL16297 = load_soft_col(proj_path("data", "raw", "platforms", "GPL16297.soft.txt"), "Gene Symbol"),
    GPL18405 = load_soft_col(proj_path("data", "raw", "platforms", "GPL18405.soft.txt"), "Symbol"),
    GPL32619 = load_soft_col(proj_path("data", "raw", "platforms", "GPL32619.soft.txt"), "ORF", 5)
  )
  seq_platforms <- c("GPL11154", "GPL16791", "GPL18573", "GPL20301", "GPL24014", "GPL24676")

  # ---- 2. layer D training intersection (all genome-wide D studies)-------------------------------------
  d <- load_studies("D")
  fg <- Filter(function(s) isTRUE(s$full_genome), d)
  uni <- Reduce(intersect, lapply(fg, function(s) rownames(s$expr)))
  log_msg("layer D training intersection: ", length(uni), " genes (", length(fg), " genome-wide studies)")

  # ---- 3. coverage table ---------------------------------------------------------------
  cov_rows <- lapply(names(plats), function(gp) {
    ov <- length(intersect(uni, plats[[gp]]))
    data.frame(platform = gp, type = "array", n_genes_platform = length(plats[[gp]]),
               n_overlap = ov, pct_of_D_intersection = round(ov / length(uni), 4),
               stringsAsFactors = FALSE)
  })
  for (gp in seq_platforms) {
    cov_rows[[length(cov_rows) + 1L]] <- data.frame(
      platform = gp, type = "sequencing", n_genes_platform = NA_integer_,
      n_overlap = length(uni), pct_of_D_intersection = 1, stringsAsFactors = FALSE)
  }
  cov <- do.call(rbind, cov_rows)
  cov <- cov[order(cov$pct_of_D_intersection), ]

  # all-C-platform intersection (start from the D intersection and intersect each array platform; sequencing platforms never restricted)
  u <- uni
  for (gp in names(plats)) u <- intersect(u, plats[[gp]])
  attr_note <- paste0("all-C-platform intersection = ", length(u), " genes; ",
                      if (length(u) < 4000) "!! <4000; drop the worst platform per the defined rule" else ">=4000, no exclusion needed")

  # ---- 4. drop-rule simulation (rule draft: drop array platforms in ascending order of coverage of the D intersection,
  #         recompute the joint intersection after each drop; stop once the intersection is >=4000. Sequencing platforms are never dropped.)----------
  remaining <- names(plats)
  dropped <- character(0)
  repeat {
    u2 <- Reduce(intersect, c(list(uni), lapply(remaining, function(gp) plats[[gp]])))
    if (length(u2) >= 4000 || length(remaining) == 0) break
    pct <- cov$pct_of_D_intersection[match(remaining, cov$platform)]
    worst <- remaining[which.min(pct)]
    dropped <- c(dropped, worst)
    remaining <- setdiff(remaining, worst)
  }
  rule_note <- paste0("rule simulation: dropping [", paste(dropped, collapse = ","), "] gives intersection = ",
                      length(u2), " genes (>=4000)")
  cov$note <- ""
  cov <- rbind(cov, data.frame(
    platform = "ALL_C_after_rule", type = "summary", n_genes_platform = NA_integer_,
    n_overlap = length(u2), pct_of_D_intersection = round(length(u2) / length(uni), 4),
    note = rule_note, stringsAsFactors = FALSE))
  write_tsv(cov, "data/c_platform_coverage_v2.tsv")
  # analysis-universe candidate (pending confirmation; script 44 uses the confirmed result)
  write_tsv(data.frame(gene = u2, stringsAsFactors = FALSE),
            "data/universe_candidate_v2.tsv")
  log_msg("coverage table: ", nrow(cov), " rows; full intersection ", length(u), "; ", rule_note)
  invisible(cov)
})
