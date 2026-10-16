# =============================================================================
# 47_c_metadata_curation.R - validation-cohort metadata curation (GEO metadata
# only; expression values are not read).
# Purpose: lower the passage/medium missing rate; add age/sex; support strata.
# Inputs: the validation study RDS files (metadata); series-matrix headers
#   (characteristics + growth protocol + title, metadata only); eutils esummary
#   JSON (data/c_layer_meta/*.json).
# Outputs: data/c_metadata_curation_v2.tsv (per-GSM old/new with evidence;
#   read-only); updated study metadata (four columns; originals kept as *_orig;
#   expr/counts untouched — predictions depend on expr only).
# Rules: sample-level key:value parsing; study-level fallback via growth
#   protocol keywords; every change evidenced in the curation table.
# =============================================================================
source("scripts/00_setup.R")

# Helper utilities (must be defined before run_step)
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a
stu_class <- function(g) {
  x <- read_tsv("curation/c_layer_studies_v2.tsv")
  x$source_class[x$gse == g][1]
}

run_step("47_c_metadata_curation", function() {
  # Target studies = the 22 C-layer studies (taken from the cohort build report)
  rep0 <- read_tsv("data/c_layer_build_report_v2.tsv")
  gses <- rep0$gse[rep0$status == "ok"]
  check(length(gses) == 22, "C-layer successful study count is not 22: ", length(gses))

  # Parsing helper: text -> key/value pairs (lowercase keys, stripped spaces)
  parse_kv <- function(txt) {
    out <- list()
    for (seg in strsplit(txt, "[;|]", perl = TRUE)[[1]]) {
      seg <- trimws(seg)
      m <- regmatches(seg, regexpr("^[A-Za-z][A-Za-z0-9 /_.-]{0,40}?\\s*:", seg))
      if (!length(m)) next
      key <- tolower(gsub("[ :]", "", m)); if (nchar(key) > 25) next
      val <- trimws(sub("^[^:]*:", "", seg))
      if (nzchar(val)) out[[key]] <- c(out[[key]], val)
    }
    out
  }
  first_of <- function(kv, keys) {
    for (k in keys) if (!is.null(kv[[k]])) return(kv[[k]][1])
    NA_character_
  }
  # Medium keywords (study-level fallback: full growth_protocol text)
  medium_from_text <- function(txt) {
    t <- tolower(txt)
    if (grepl("platelet lysate|\\bhpl\\b", t)) return("hPL")
    if (grepl("xeno.?free|serum.?free|\\bsf\\b|\\bxf\\b", t)) return("xeno_free")
    if (grepl("fetal (bovine|calf)|\\bfbs\\b|\\bfcs\\b", t)) return("FBS")
    if (grepl("human serum|\\bhs\\b", t)) return("HS")
    NA_character_
  }
  passage_from_text <- function(txt) {
    m <- regmatches(txt, regexpr("(?i)passage[^0-9]{0,8}([0-9]{1,2})", txt, perl = TRUE))
    if (length(m)) return(as.numeric(sub(".*?([0-9]{1,2})$", "\\1", m)))
    m <- regmatches(txt, regexpr("(?i)\\bP[0-9]{1,2}\\b", txt, perl = TRUE))
    if (length(m)) return(as.numeric(sub("(?i)P", "", m)))
    NA_real_
  }
  age_from_text <- function(txt) {
    m <- regmatches(txt, regexpr("(?i)(age|donor.?age)[^0-9]{0,8}([0-9]{1,2})", txt, perl = TRUE))
    if (length(m)) return(as.numeric(sub(".*?([0-9]{1,2})$", "\\1", m)))
    NA_real_
  }
  sex_from_text <- function(txt) {
    t <- tolower(txt)
    f_hit <- grepl("female|\\bwoman\\b|\\bwomen\\b|(gender|sex)[ :=]+f($|[^a-z])", t)
    m_hit <- (grepl("\\bmale\\b|\\bman\\b|\\bmen\\b", t) & !grepl("female", t)) |
             grepl("(gender|sex)[ :=]+m($|[^a-z])", t)
    if (f_hit && !m_hit) "F" else if (m_hit && !f_hit) "M" else NA_character_
  }

  cur <- list(); note_log <- character(0)
  for (g in gses) {
    sfile <- list.files(proj_path("data", "raw", g), pattern = "series_matrix.*gz", full.names = TRUE)[1]
    hdr <- read_series_header(sfile)
    ch_cols <- grep("^Sample_characteristics", names(hdr), value = TRUE)
    gp_cols <- grep("^Sample_growth_protocol", names(hdr), value = TRUE)
    ch_txt <- apply(hdr[, ch_cols, drop = FALSE], 1, function(v) paste(v[!is.na(v)], collapse = " ; "))
    gp_txt <- if (length(gp_cols)) apply(hdr[, gp_cols, drop = FALSE], 1, function(v) paste(v[!is.na(v)], collapse = " ; ")) else rep("", nrow(hdr))
    tt <- if (!is.null(hdr$Sample_title)) hdr$Sample_title else rep("", nrow(hdr))
    # Study-level fallback: merged growth_protocol text across the whole study
    gp_all <- paste(gp_txt, collapse = " ")
    med_study <- medium_from_text(gp_all)
    if (is.na(med_study)) {   # fallback: GSE summary from esummary (a few studies carry medium clues)
      ej <- tryCatch(jsonlite::read_json(proj_path("data", "c_layer_meta", paste0("esummary_", stu_class(g), ".json")), simplifyVector = TRUE),
                     error = function(e) NULL)
      if (!is.null(ej) && !is.null(ej$result[[paste0("200", sub("^GSE", "", g))]])) {
        gs <- ej$result[[paste0("200", sub("^GSE", "", g))]]
        med_study <- medium_from_text(paste(gs$title %||% "", gs$summary %||% "", collapse = " "))
      }
    }
    s <- readRDS(proj_path("data", "studies", paste0(g, ".rds")))
    gm <- hdr$Sample_geo_accession
    for (i in seq_len(nrow(s$meta))) {
      gsm <- s$meta$gsm[i]
      j <- match(gsm, gm)
      if (is.na(j)) { note_log <- c(note_log, paste(g, gsm, "not in series matrix")); next }
      ev <- paste(ch_txt[j], gp_txt[j], tt[j], sep = " ; ")
      kv <- parse_kv(ev)
      new_passage <- suppressWarnings(as.numeric(first_of(kv, c("passage", "passagenumber", "passageno"))))
      if (is.na(new_passage)) new_passage <- passage_from_text(tt[j])
      new_medium <- medium_from_text(first_of(kv, c("medium", "culturemedium", "growthmedium", "serum", "mediamedium")) %||% "")
      if (is.na(new_medium)) new_medium <- medium_from_text(ev)   # no keywords in the key/value pairs (e.g. 'medium: control') -> remap from the full evidence text
      if (is.na(new_medium)) new_medium <- med_study
      new_age <- suppressWarnings(as.numeric(first_of(kv, c("age", "donorage", "donoragenumeric"))))
      if (is.na(new_age)) new_age <- age_from_text(ev)
      new_sex <- first_of(kv, c("sex", "gender"))
      if (!is.na(new_sex)) new_sex <- sex_from_text(new_sex) %||% new_sex
      if (is.na(new_sex)) new_sex <- sex_from_text(ev)
      cur[[length(cur) + 1L]] <- data.frame(
        gse = g, gsm = gsm,
        passage_old = s$meta$passage[i], passage_new = new_passage,
        medium_old = s$meta$medium[i], medium_new = new_medium,
        sex_old = s$meta$sex[i], sex_new = new_sex,
        age_new = new_age,
        evidence = substr(gsub("[\t\n]", " ", ev), 1, 300),
        stringsAsFactors = FALSE)
    }
    log_msg("47 processing ", g, ": ", nrow(s$meta), " samples")
  }
  ct <- do.call(rbind, cur)

  # Update the four meta columns in the rds (originals kept as *_orig; expr/counts untouched)
  for (g in gses) {
    s <- readRDS(proj_path("data", "studies", paste0(g, ".rds")))
    sub <- ct[ct$gse == g, , drop = FALSE]
    k <- match(s$meta$gsm, sub$gsm)
    if (!"passage_orig" %in% names(s$meta)) {   # idempotent: record originals only on the first run (so a re-run does not treat new values as orig)
      s$meta$passage_orig <- s$meta$passage; s$meta$medium_orig <- s$meta$medium; s$meta$sex_orig <- s$meta$sex
    }
    s$meta$passage <- sub$passage_new[k]
    s$meta$medium <- sub$medium_new[k]
    s$meta$sex <- sub$sex_new[k]
    s$meta$age <- sub$age_new[k]
    validate_study(s)
    saveRDS(s, proj_path("data", "studies", paste0(g, ".rds")))
  }

  # Missing-rate before/after comparison + write
  n <- nrow(ct)
  ba <- data.frame(
    field = c("passage", "medium", "sex", "age"),
    na_before = c(sum(is.na(ct$passage_old)), sum(is.na(ct$medium_old)), sum(is.na(ct$sex_old)), NA_integer_),
    na_after = c(sum(is.na(ct$passage_new)), sum(is.na(ct$medium_new)), sum(is.na(ct$sex_new)), sum(is.na(ct$age_new))))
  ct2 <- ct; ct2$post_hoc <- TRUE
  out <- rbind(
    data.frame(gse = "(SUMMARY)", gsm = paste0("n=", n), passage_old = ba$na_before[1], passage_new = ba$na_after[1],
               medium_old = ba$na_before[2], medium_new = ba$na_after[2], sex_old = ba$na_before[3], sex_new = ba$na_after[3],
               age_new = ba$na_after[4], evidence = "NA counts use the before/after definition", stringsAsFactors = FALSE),
    ct[, names(ct)])
  write_tsv(out, "data/c_metadata_curation_v2.tsv")
  log_msg("47 missingness (NA counts): passage ", ba$na_before[1], "->", ba$na_after[1],
          "; medium ", ba$na_before[2], "->", ba$na_after[2],
          "; sex ", ba$na_before[3], "->", ba$na_after[3], "; new NA in age ", ba$na_after[4])
  invisible(ct)
})
