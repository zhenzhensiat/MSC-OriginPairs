# =============================================================================
# 01i_verify_c.R —— per-study verification of validation candidates
# (eligibility criteria; metadata and publications only, no expression values).
# Europe PMC method sections are checked for: (1) standard culture arm;
# (2) donor numbers; (3) source confirmation.
# =============================================================================
source("scripts/00_setup.R")
need_pkgs(c("httr", "jsonlite"))

run_step("01i_verify_c", function() {
  led <- read_tsv("data/census_adjudication_v2_r3.tsv")
  net <- read_tsv("data/independence_network_v2_r2.tsv")
  geo <- read_tsv("data/census_geo_v2.tsv")
  alloc <- read_tsv("data/layer_allocation_v2.tsv")

  d_cl <- unique(net$cluster_id[net$gse %in% alloc$gse[alloc$layer == "D"]])

  c_gse <- led$accession[led$human_verdict %in% c("C_AD","C_BM","C_UC")]
  c_gse <- setdiff(c_gse, c("GSE194298"))   # this series has been moved out
  geo_c <- setdiff(alloc$gse[alloc$layer == "C"], c("GSE189895","GSE37065","GSE9520"))
  all_c <- sort(unique(c(c_gse, geo_c)))
  cat("C candidate verification: ", length(all_c), " items\n")

  base <- "https://www.ebi.ac.uk/europepmc/webservices/rest"
  results <- list()

  for (g in all_c) {
    vd <- led$human_verdict[led$accession == g]
    if (length(vd) == 0) vd <- "C_GEO"
    idx <- match(g, geo$gse)
    pm <- if (!is.na(idx)) as.character(geo$pmid[idx]) else NA_character_
    pm <- trimws(strsplit(as.character(pm), ";", fixed=TRUE)[[1]][1])
    if (is.na(pm) || pm == "NA") pm <- NA_character_
    # take PMIDs from the candidates table (adjudicated candidates may not be in the geo table)
    if (is.na(pm)) {
      pm <- led$pubmed[led$accession == g]
      if (length(pm) > 0) pm <- trimws(strsplit(as.character(pm), ";", fixed=TRUE)[[1]][1])
      if (is.na(pm) || pm == "NA") pm <- NA_character_
    }
    n_geo <- NA_real_

    # look up the publication in Europe PMC
    arms <- ""; donors_txt <- ""; n_donor <- NA_real_; evidence <- ""
    if (!is.na(pm) && nzchar(pm) && pm != "NA") {
      r <- tryCatch({
        u <- paste0(base, "/search?query=EXT_ID:", pm, "&format=json&resultType=core&pageSize=1")
        js <- jsonlite::fromJSON(readLines(u, warn=FALSE))
        res <- js$resultList$result
        if (is.null(res) || nrow(res)==0) return(NULL)
        abs_txt <- tolower(as.character(res$abstractText[1]))
        oa <- res$isOpenAccess[1] == "Y"
        ft <- ""
        if (oa && !is.na(res$pmcid[1])) {
          ft <- tryCatch(paste(readLines(paste0(base, "/", res$pmcid[1], "/fullTextXML"), warn=FALSE),
                              collapse=" "), error=function(e) "")
        }
        txt <- tolower(paste(abs_txt, ft))
        # (1) standard arm
        arm_kw <- c("healthy (donor|control|subject|volunteer)","untreated","unstimulated","normoxi",
                    "control (medium|condition|group)","conventional (culture|2d|monolayer|plate)",
                    "non-obese","normal (glucose|donor|adult)","standard (culture|condition)",
                    "fbs.*control|control.*fbs","without (treatment|stimul|fbs)","parental","wild.?type")
        hits <- arm_kw[vapply(arm_kw, function(a) grepl(a, txt), logical(1))]
        arms <- paste(hits, collapse=";")
        # (2) donor count (extracted from the Methods section)
        dn <- regmatches(txt, gregexpr(
          "([0-9]+|two|three|four|five|six|seven|eight|nine|ten)[^.]{0,20}(healthy )?(donors?|subjects?|patients?|volunteers?)",
          txt))[[1]]
        donors_txt <- paste(head(unique(dn), 5), collapse="; ")
        if (length(dn) > 0) {
          nums <- c("two"=2,"three"=3,"four"=4,"five"=5,"six"=6,"seven"=7,"eight"=8,"nine"=9,"ten"=10)
          for (d in dn) {
            m <- as.numeric(sub("([0-9]+).*", "\1", d))
            w <- names(nums)[vapply(names(nums), function(n) grepl(n, d), logical(1))]
            v <- if (!is.na(m)) m else if (length(w)>0) nums[w[1]] else NA
            if (!is.na(v) && (is.na(n_donor) || v > n_donor)) n_donor <- v
          }
        }
        evidence <- if (nzchar(ft)) "fulltext" else if (nzchar(abs_txt)) "abstract" else "none"
        Sys.sleep(0.3)
        NULL
      }, error=function(e) NULL)
    }
    # independence
    bad_cl <- if (g %in% net$gse && net$cluster_id[net$gse==g] %in% d_cl) "D_CLUSTER" else ""
    results[[g]] <- data.frame(
      gse = g, verdict = vd, n_geo = n_geo,
      arm_evidence = substr(arms, 1, 120), evidence_source = evidence,
      donor_hint = substr(donors_txt, 1, 80), n_donor_max = n_donor,
      d_cluster = bad_cl, stringsAsFactors = FALSE, row.names = NULL)
    if (length(results) %% 10 == 0) cat(length(results), "/", length(all_c), "\n")
  }

  res <- do.call(rbind, results)
  write_tsv(res, "data/c_verify_detail.tsv")
  # summary
  for (cl in c("C_AD","C_BM","C_UC","C_GEO")) {
    s <- res[res$verdict == cl, , drop=FALSE]
    if (nrow(s) == 0) next
    cat(sprintf("%s: %d candidates | with arm evidence %d | D cluster %d | ≥3 donors %d\n",
      cl, nrow(s), sum(nzchar(s$arm_evidence)), sum(s$d_cluster == "D_CLUSTER"),
      sum(!is.na(s$n_donor_max) & s$n_donor_max >= 3)))
  }
  invisible(res)
})
