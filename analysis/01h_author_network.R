# =============================================================================
# 01h_author_network.R —— author network and independence clusters for the
# development, supplementary and validation-candidate studies.
# Sources: v1 independence_clusters_v1 (contributors column) + census_geo_v2
# (PMIDs) with Europe PMC filling gaps.
# Rule: studies sharing a corresponding/senior author, or >=2 any authors,
# belong to the same group. Network cluster codes are mapped to the published
# neutral numbering at the end of the construction.
# =============================================================================
source("scripts/00_setup.R")
need_pkgs(c("httr", "jsonlite"))

run_step("01h_author_network", function() {
  led <- read_tsv("data/census_adjudication_v2_r3.tsv")
  geo <- read_tsv("data/census_geo_v2.tsv")
  v1_ind <- read_tsv("data/source_v1/independence_clusters_v1.tsv")

  # collect GSEs needing verification: all of layer_allocation + adjudicated C candidates + MULTI_NEW
  alloc <- read_tsv("data/layer_allocation_v2.tsv")
  d_gse <- alloc$gse[alloc$layer == "D"]
  supp_gse <- alloc$gse[alloc$layer == "supp"]
  alloc_c_gse <- alloc$gse[alloc$layer == "C"]       # includes GEO entries assigned to C in the final table (GSE37065 etc.)
  led_c_gse <- led$accession[led$human_verdict %in% c("C_AD","C_BM","C_UC")]
  c_gse <- sort(unique(c(alloc_c_gse, led_c_gse)))    # C layer = union of allocation table and adjudicated candidates
  multi_gse <- led$accession[led$human_verdict == "MULTI_NEW"]
  all_gse <- sort(unique(c(d_gse, supp_gse, c_gse, multi_gse)))
  log_msg("network coverage: D ", length(d_gse), " + supp ", length(supp_gse), " + C ", length(c_gse),
          " + multi ", length(multi_gse), " = ", length(all_gse), " GSEs")

  # author source 1: contributors column of the v1 independence table (split "Last,First" pairs)
  v1_auth <- strsplit(as.character(v1_ind$contributors), "; ")
  names(v1_auth) <- v1_ind$gse

  # author source 2: PMIDs from census_geo_v2 → Europe PMC
  fetch_pmc_authors <- function(pmids) {
    out <- list()
    for (pm in pmids) {
      pm <- trimws(strsplit(as.character(pm), ";", fixed = TRUE)[[1]][1])
      if (is.na(pm) || !nzchar(pm) || pm == "NA") next
      r <- tryCatch({
        u <- paste0("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:", pm,
                    "&format=json&resultType=core&pageSize=1")
        js <- jsonlite::fromJSON(readLines(u, warn = FALSE))
        res <- js$resultList$result
        if (is.null(res) || nrow(res) == 0) return(NULL)
        list(authors = strsplit(as.character(res$authorString[1]), ", ")[[1]],
             last = res$authorString[1])
      }, error = function(e) NULL)
      if (!is.null(r)) out[[pm]] <- r
      Sys.sleep(0.3)
    }
    out
  }
  geo_pmid <- setNames(as.character(geo$pmid), geo$gse)
  pmc_raw <- fetch_pmc_authors(unique(na.omit(geo_pmid[all_gse])))
  pmc_auth <- lapply(as.character(geo_pmid[all_gse]), function(pm) {
    pm <- trimws(strsplit(pm, ";", fixed = TRUE)[[1]][1])
    if (pm %in% names(pmc_raw)) pmc_raw[[pm]]$authors else NULL
  })
  names(pmc_auth) <- all_gse

  # merge the author sets from the two sources
  auth_of <- lapply(all_gse, function(g) {
    a1 <- if (g %in% names(v1_auth)) v1_auth[[g]] else character(0)
    a2 <- if (!is.null(pmc_auth[[g]])) pmc_auth[[g]] else character(0)
    unique(trimws(c(a1, a2)))
  })
  names(auth_of) <- all_gse
  names(auth_of) <- all_gse

  # rule: sharing ≥2 authors of any position → same group (corresponding/last-author roles cannot be told apart a priori from the data structure; a conservative threshold is used)
  shared_count <- function(a1, a2) length(intersect(a1, a2))
  n <- length(all_gse)
  adj <- matrix(FALSE, n, n, dimnames = list(all_gse, all_gse))
  for (i in seq_len(n - 1)) for (j in (i + 1):n) {
    if (j > n) break
    sc <- shared_count(auth_of[[all_gse[i]]], auth_of[[all_gse[j]]])
    if (sc >= 2) adj[i, j] <- adj[j, i] <- TRUE
  }
  # connected components -> groups
  gr <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected")
  comp <- igraph::components(gr)
  cl_id <- paste0("cl_net_", comp$membership)
  names(cl_id) <- all_gse

  # curated group labels (series with no PMID/Europe PMC record)
  g02grp <- c("GSE57151", "GSE12274", "GSE37065", "GSE9520", "GSE189895",
              "GSE194298")  # 194298 = RNA-seq subseries of the same hPL-culture series
  g44grp <- c("GSE291956", "GSE309666")
  for (g in g02grp) if (g %in% all_gse) cl_id[g] <- "cl_net_wagner"
  for (g in g44grp) if (g %in% all_gse) cl_id[g] <- "cl_net_tu_hu"

  out <- data.frame(
    gse = all_gse,
    authors = vapply(auth_of[all_gse], function(a) paste(head(a, 10), collapse = "; "), character(1)),
    n_authors = lengths(auth_of[all_gse]),
    cluster_id = cl_id,
    stringsAsFactors = FALSE, row.names = NULL)
  # linked_gse / shared_authors columns
  out$linked_gse <- vapply(all_gse, function(g) {
    peers <- all_gse[cl_id == cl_id[g] & all_gse != g]
    if (length(peers) == 0) "" else paste(peers, collapse = ";")
  }, character(1))
  write_tsv(out, "data/independence_network_v2_r2.tsv")

  # summary
  tab <- table(cl_id)
  multi_cl <- names(tab)[tab > 1]
  log_msg("author network: ", length(unique(cl_id)), " groups (multi-member groups: ", length(multi_cl), "):",
          paste(multi_cl, tab[multi_cl], sep = "=", collapse = ", "))
  # C candidates in the same group as the D layer
  d_cl <- unique(cl_id[d_gse[d_gse %in% names(cl_id)]])
  c_bad <- c_gse[c_gse %in% names(cl_id) & cl_id[c_gse] %in% d_cl]
  log_msg("!! C candidates in the same group as the D layer (must be moved out of C): ", paste(c_bad, collapse = ","))
  invisible(out)
})
