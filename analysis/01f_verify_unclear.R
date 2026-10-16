# =============================================================================
# 01f_verify_unclear.R —— sample-level metadata verification of validation
# candidates and UNCLEAR cases (metadata only; expression matrices untouched).
# Purpose: the eligibility threshold (>=8 studies per source) requires resolving
# UNCLEAR cases - control-arm existence and donor composition.
# Method: NCBI eutils esearch+esummary pull sample titles/characteristics per
# GSE; machine summarisation (arm keywords x donor counts) plus a manual review list.
# Input: data/census_adjudication_v2_r2.tsv (read-only)
# Outputs: data/unclear_verify_eutils.tsv, docs/unclear_verify_review.md
# Rate limiting: >=0.4 s between requests (NCBI limit 3/s); 3 retries on failure
# =============================================================================
source("scripts/00_setup.R")
need_pkgs(c("httr", "jsonlite"))

run_step("01f_verify_unclear", function() {
  led <- read_tsv("data/census_adjudication_v2_r2.tsv")
  tgt <- led[led$human_verdict %in% c("UNCLEAR", "C_AD", "C_BM", "C_UC"), ]
  check(nrow(tgt) > 50, "unexpected target count: ", nrow(tgt))
  log_msg("verification targets: UNCLEAR ", sum(tgt$human_verdict == "UNCLEAR"), " + C candidates ",
          sum(tgt$human_verdict != "UNCLEAR"))

  base <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
  fetch <- function(gse) {
    for (a in 1:3) {
      r <- tryCatch({
        q <- httr::GET(paste0(base, "/esearch.fcgi"), query = list(db = "gds", term = gse, retmode = "json"),
                       timeout(30))
        id <- jsonlite::fromJSON(rawToChar(q$content))$esearchresult$idlist[1]
        if (is.na(id) || !nzchar(id)) stop("no id")
        Sys.sleep(0.4)
        s <- httr::GET(paste0(base, "/esummary.fcgi"), query = list(db = "gds", id = id, retmode = "json"),
                       timeout(30))
        js <- jsonlite::fromJSON(rawToChar(s$content))
        unlist(js$result[[as.character(id)]][c("n_samples", "sampleinfo_title", "sampleinfo_attribute")])
      }, error = function(e) NULL)
      if (!is.null(r)) return(r)
      Sys.sleep(2 * a)
    }
    NULL
  }

  rows <- list()
  for (i in seq_len(nrow(tgt))) {
    g <- tgt$accession[i]
    sm <- fetch(g)
    if (is.null(sm)) {
      rows[[g]] <- data.frame(gse = g, verdict = tgt$human_verdict[i], n_samples_geo = NA,
                              arm_hint = "FETCH_FAIL", donor_hint = "", txt = "")
      next
    }
    ti <- sm["sampleinfo_title"]; at <- sm["sampleinfo_attribute"]
    ti <- if (is.na(ti)) "" else ti
    at <- if (is.na(at)) "" else at
    txt <- tolower(paste(ti, at, sep = " || "))
    arms <- c("untreated", "unstimulated", "normox", "control", "healthy", "non-obese", "normal glucose",
              "1g/l", "fbs", "dmem", "basal", "vehicle", "wildtype", "wt", "parental", "2d", "monolayer",
              "standard", "conventional", "p0", "passage 0", "early passage")
    hit <- arms[ vapply(arms, function(a) grepl(a, txt, fixed = TRUE), logical(1)) ]
    donors <- unique(c(
      regmatches(txt, gregexpr("donor[ _:]*[0-9]{1,2}", txt))[[1]],
      regmatches(txt, gregexpr("\\b(bm|ad|uc|d|hd)[0-9]{2,3}\\b", txt))[[1]]))
    rows[[g]] <- data.frame(
      gse = g, verdict = tgt$human_verdict[i], n_samples_geo = suppressWarnings(as.integer(sm["n_samples"])),
      arm_hint = paste(hit, collapse = ";"), donor_hint = paste(head(donors, 12), collapse = ","),
      txt = substr(txt, 1, 400), stringsAsFactors = FALSE)
    if (i %% 10 == 0) log_msg(i, "/", nrow(tgt), " done")
    Sys.sleep(0.4)
  }
  v <- do.call(rbind, rows)
  write_tsv(v, "data/unclear_verify_eutils.tsv")

  # manual review list: UNCLEAR first (an arm hit upgrades the candidate); donor hints spot-checked for C candidates
  lines <- c("# Sample-level verification hints for UNCLEAR/C candidates (eutils metadata, for manual review)", "")
  for (vd in c("UNCLEAR", "C_AD", "C_BM", "C_UC")) {
    s <- v[v$verdict == vd, ]
    lines <- c(lines, paste0("## ", vd, " (", nrow(s), ")"), "")
    for (i in seq_len(nrow(s))) lines <- c(lines, sprintf(
      "- **%s** n=%s | arm hint [%s] | donor [%s] | %s",
      s$gse[i], s$n_samples_geo[i], substr(s$arm_hint[i], 1, 80),
      substr(s$donor_hint[i], 1, 60), substr(s$txt[i], 1, 220)))
    lines <- c(lines, "")
  }
  writeLines(lines, proj_path("docs", "unclear_verify_review.md"))
  log_msg("verification hints written: ", nrow(v), " rows (FETCH_FAIL ",
          sum(v$arm_hint == "FETCH_FAIL"), ")")
  invisible(v)
})
