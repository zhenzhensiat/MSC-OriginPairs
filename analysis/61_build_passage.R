# =============================================================================
# 61_build_passage.R - build the passage-series studies from eligible
# serial-passage candidates.
# Preconditions: curation/passage_manifest_v2.tsv ready (sample-level: gse,
#   gsm, source, donor, passage, medium, arm); raw matrices downloaded.
# Flow: per GSE call build_from_manifest (utils_build) -> stability=TRUE ->
#   write data/studies/<GSE>.rds. Only arm=="include" samples; passage coerced
#   to numeric (prevents lexicographic early/late reversal); platform hint from
#   the series-matrix header; per-study failures are recorded and skipped (no
#   silent drop; unrecorded failure -> exit).
# Outputs: study objects (stability=TRUE); data/passage_build_report_v2.tsv (written).
# =============================================================================
source("scripts/00_setup.R")
source(proj_path("R", "utils_build.R"), encoding = "UTF-8")
need_pkgs(c("org.Hs.eg.db"))

run_step("61_build_passage", function() {
  man <- read_tsv("curation/passage_manifest_v2.tsv")
  check(all(c("gse", "gsm", "source", "donor", "passage", "medium", "arm") %in% names(man)),
        "manifest is missing required columns")
  gses <- unique(man$gse)
  audit <- list()
  for (g in gses) {
    mr <- man[man$gse == g & man$arm == "include", , drop = FALSE]
    mr$passage <- suppressWarnings(as.numeric(mr$passage))
    check(nrow(mr) >= 3L, g, ": include samples ", nrow(mr), " < 3")
    # Platform hint: Sample_platform_id in the series matrix header (extraction failure degrades to NULL, recorded faithfully by the build layer)
    gpl_hint <- tryCatch({
      fp <- list.files(proj_path("data", "raw", g), pattern = "series_matrix.*gz", full.names = TRUE)[1]
      hd <- read_series_header(fp)
      pu <- unique(sub("^GPL", "", hd$Sample_platform_id))
      if (length(pu) == 1L && grepl("^[0-9]+$", pu)) pu else NULL
    }, error = function(e) NULL)
    info <- data.frame(gse = g, n_samples = nrow(mr), n_donors = length(unique(mr$donor)),
                       n_passage_vals = length(unique(mr$passage[!is.na(mr$passage)])),
                       n_donors_multi_passage = NA_integer_,
                       status = "pending", note = "", stringsAsFactors = FALSE)
    s <- tryCatch(
      build_from_manifest(g, mr, layer = "supp", cluster_id = paste0("passage_", g),
                          stability = TRUE, min_genes = 1000L, min_samples = 3L,
                          platform_hint = gpl_hint),
      error = function(e) conditionMessage(e))
    if (is.character(s)) {
      info$status <- paste0("build_fail: ", substr(s, 1, 200))
      audit[[g]] <- info
      log_msg("61 ", g, ": build failed - ", substr(s, 1, 120))
      next
    }
    # Donor-level passage check (count of donors with >=2 passage points)
    dmulti <- sum(vapply(split(s$meta, s$meta$donor),
                         function(d) length(unique(d$passage[!is.na(d$passage)])) >= 2, logical(1)))
    info$n_donors_multi_passage <- dmulti
    info$status <- if (dmulti >= 3) "ok" else "insufficient_multi_passage_donors"
    info$note <- paste0("scale=", s$scale, "; genes=", nrow(s$expr))
    if (identical(info$status, "ok")) {
      saveRDS(s, proj_path("data", "studies", paste0(g, ".rds")))
      log_msg("61 build ", g, ": ", nrow(s$expr), " genes x ", ncol(s$expr), " samples; >=3 donors with >=2 passage points (", dmulti, ")")
    } else {
      log_msg("61 ", g, ": insufficient donor passage points (", dmulti, " <3), not carried into 60")
    }
    audit[[g]] <- info
  }
  rep_tab <- do.call(rbind, audit)
  write_tsv(rep_tab, "data/passage_build_report_v2.tsv")
  log_msg("61 build audit: ", sum(rep_tab$status == "ok"), " / ", nrow(rep_tab), " passed (stability=TRUE written to disk)")
  invisible(rep_tab)
})
