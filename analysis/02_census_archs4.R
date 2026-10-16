# =============================================================================
# 02_census_archs4.R —— ARCHS4 metadata census (metadata only; expression values
# are never read).
# Inputs: data/raw/archs4/<H5 file> (file name via environment variable ARCHS4_FILE);
#         config/archs4_rules.tsv (reviewed before running);
#         data/census_geo_v2.tsv (census table, needs a gse column)
# Outputs: data/census_archs4_samples_v2.tsv, data/census_archs4_series_v2.tsv,
#          docs/archs4_audit_template.tsv (>=10% manual spot-check template)
# =============================================================================
source("scripts/00_setup.R")
need_pkgs(c("rhdf5"))

run_step("02_census_archs4", function() {
  h5 <- proj_path("data", "raw", "archs4", Sys.getenv("ARCHS4_FILE", "human_gene_v2.latest.h5"))
  check(file.exists(h5), "ARCHS4 file not found: ", h5)

  listing <- rhdf5::h5ls(h5)
  paths <- sub("^//", "/", paste(listing$group, listing$name, sep = "/"))
  pick <- function(cands) {
    hit <- cands[cands %in% paths]
    check(length(hit) >= 1, "field not found in H5: ", paste(cands, collapse = " / "),
          ". Run print(rhdf5::h5ls(h5)) first to check the field names of the current version")
    hit[1]
  }
  rd <- function(p) as.character(rhdf5::h5read(h5, p))
  # VERIFY 2026-10-10: the v11 snapshot uses flat /meta/Sample_* names (only v2-series files use /meta/samples/*),
  # candidates tried in order; once a file is in place, verify it with h5ls before running
  meta <- data.frame(
    gsm             = rd(pick(c("/meta/samples/geo_accession",    "/meta/Sample_geo_accession"))),
    series          = rd(pick(c("/meta/samples/series_id",       "/meta/Sample_series_id"))),
    title           = rd(pick(c("/meta/samples/title",           "/meta/Sample_title"))),
    source_name     = rd(pick(c("/meta/samples/source_name_ch1", "/meta/Sample_source_name_ch1"))),
    characteristics = rd(pick(c("/meta/samples/characteristics_ch1", "/meta/Sample_characteristics_ch1"))),
    stringsAsFactors = FALSE)
  rhdf5::h5closeAll()
  log_msg("ARCHS4 total samples: ", nrow(meta))

  rules <- read_tsv("config/archs4_rules.tsv")
  txt <- tolower(paste(meta$title, meta$source_name, meta$characteristics, sep = " | "))
  rx <- function(type) rules$regex[rules$rule_type == type]

  meta$is_msc <- grepl(rx("msc"), txt, perl = TRUE)
  meta$excluded_arm <- grepl(rx("exclude_arm"), txt, perl = TRUE)

  src <- rules[rules$rule_type == "source", ]
  hits <- vapply(src$regex, function(r) grepl(r, txt, perl = TRUE), logical(nrow(meta)))
  if (!is.matrix(hits)) hits <- matrix(hits, nrow = nrow(meta))
  meta$source_call <- apply(hits, 1, function(h) {
    if (!any(h)) return("unknown")
    pr <- src$priority[h]
    top <- unique(src$label[h][pr == min(pr)])
    if (length(top) == 1) top else "ambiguous"
  })

  msc <- meta[meta$is_msc, , drop = FALSE]
  log_msg("samples called MSC: ", nrow(msc))
  write_tsv(msc, "data/census_archs4_samples_v2.tsv")

  # series aggregation (the series field may contain multiple GSEs, comma-separated)
  sl <- strsplit(msc$series, "[,;[:space:]]+")
  long <- data.frame(gse = unlist(sl), gsm = rep(msc$gsm, lengths(sl)),
                     source = rep(msc$source_call, lengths(sl)),
                     excluded = rep(msc$excluded_arm, lengths(sl)), stringsAsFactors = FALSE)
  long <- long[nzchar(long$gse), ]
  incl <- long[!long$excluded, ]
  by_gse <- split(incl, incl$gse)
  series <- data.frame(
    gse = names(by_gse),
    n_samples = vapply(by_gse, nrow, integer(1)),
    n_AD = vapply(by_gse, function(d) sum(d$source == "AD"), integer(1)),
    n_BM = vapply(by_gse, function(d) sum(d$source == "BM"), integer(1)),
    n_UC = vapply(by_gse, function(d) sum(d$source == "UC"), integer(1)),
    n_other_or_unclear = vapply(by_gse, function(d) sum(!d$source %in% CFG$classes), integer(1)),
    stringsAsFactors = FALSE, row.names = NULL)
  series$n_sources <- (series$n_AD > 0) + (series$n_BM > 0) + (series$n_UC > 0)
  series$design <- ifelse(series$n_sources >= 2, "multi_source", ifelse(series$n_sources == 1, "single_source", "none"))

  geo <- read_tsv("data/census_geo_v2.tsv")
  series$in_geo_census <- series$gse %in% geo$gse
  write_tsv(series, "data/census_archs4_series_v2.tsv")

  # manual spot-check: sample ≥10% per source class (at least 10), cross-checked against the original GEO pages
  set.seed(derive_seed(2, 1))
  aud <- do.call(rbind, lapply(split(msc, msc$source_call), function(d) {
    n <- min(nrow(d), max(10L, ceiling(0.10 * nrow(d))))
    d[sample.int(nrow(d), n), c("gsm", "series", "title", "source_name", "characteristics", "source_call", "excluded_arm")]
  }))
  aud$manual_is_msc <- ""; aud$manual_source <- ""; aud$manual_excluded_arm <- ""; aud$note <- ""
  write_tsv(aud, "docs/archs4_audit_template.tsv")

  log_msg("multi-source series: ", sum(series$design == "multi_source"), "; single-source series: ",
          sum(series$design == "single_source"), " (overlap with the GEO census not removed; independence clustering not performed)")
  invisible(series)
})
