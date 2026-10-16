# =============================================================================
# 20_marker_list.R —— port and re-verify the marker lists
# (runs before any expression data are inspected).
# (1) marker_list (82 claims): mechanical re-verification of PMID format, current
#     HGNC symbols, high/low in {BM,AD,UC}, contrast direction consistency
# (2) HOX list (39 genes): current HGNC symbols + family membership
# (3) entries failing verification are flagged, never deleted (audit by flag)
# Outputs: curation/marker_list.tsv, curation/hox_list.tsv
# =============================================================================
source("scripts/00_setup.R")
need_pkgs(c("org.Hs.eg.db", "AnnotationDbi"))
library(org.Hs.eg.db)

run_step("20_marker_list", function() {
  v1_dir <- file.path(normalizePath(file.path(getOption("msc.root"), "..", "..")),
                      "\u7814\u7a76\u65b9\u54112")
  src <- file.path(v1_dir, "data", c("marker_list_v1.tsv", "hox_list_v1.tsv"))
  check(all(file.exists(src)), "v1 lists missing: ", paste(src[!file.exists(src)], collapse = ","))
  for (p in src) if (!file.exists(proj_path("data", "source_v1", basename(p))))
    file.copy(p, proj_path("data", "source_v1"))

  ml <- read_tsv("data/source_v1/marker_list_v1.tsv")
  hl <- read_tsv("data/source_v1/hox_list_v1.tsv")
  check(all(c("claim_id", "gene", "high_source", "low_source", "contrast", "pmid") %in% names(ml)),
        "marker list is missing columns")
  check(nrow(ml) == 82, "v1 marker claims should be 82 (81 + final addendum), actual ", nrow(ml))
  check(nrow(hl) == 39, "v1 HOX should be 39, actual ", nrow(hl))

  # ---- HGNC current symbol mapping (same source as script 11)--------------------------------------
  eg <- AnnotationDbi::select(org.Hs.eg.db, keys = keys(org.Hs.eg.db),
                              columns = c("SYMBOL", "ALIAS"))
  sym_ok <- unique(na.omit(eg$SYMBOL))
  al <- eg[!is.na(eg$ALIAS) & !is.na(eg$SYMBOL), c("ALIAS", "SYMBOL")]
  al_n <- table(al$ALIAS); uniq_alias <- names(al_n[al_n == 1L])
  alias_map <- setNames(al$SYMBOL[match(uniq_alias, al$ALIAS)], uniq_alias)

  ml$hgnc_symbol <- ifelse(ml$gene %in% sym_ok, ml$gene,
                    ifelse(ml$gene %in% names(alias_map), alias_map[ml$gene], NA_character_))
  hsym <- hl$symbol                              # v1 HOX table: column 1 is the HGNC ID; symbols are in the symbol column
  hl$hgnc_symbol <- ifelse(hsym %in% sym_ok, hsym,
                    ifelse(hsym %in% names(alias_map), alias_map[hsym], NA_character_))

  # ---- discovery_gse: PMID -> GSE (against the v1 census table; most claims come from training-set papers, which v2 must exclude)----
  geo <- read_tsv("data/census_geo_v2.tsv")
  pmid_to_gse <- function(p) {
    if (is.na(p) || !nzchar(p)) return(NA_character_)
    hits <- geo$gse[vapply(geo$pmid, function(x) !is.na(x) && p %in% strsplit(x, ";", fixed = TRUE)[[1]],
                           logical(1))]
    if (length(hits) == 0) return(NA_character_)
    if (length(hits) > 1) return(paste(hits, collapse = ";"))   # multiple hits: registered for manual adjudication
    hits
  }
  ml$discovery_gse <- vapply(ml$pmid, pmid_to_gse, character(1))

  # ---- source standardisation (manual standardisation required before expand_claims)----
  # rules: (1) strip parenthetical qualifiers (adult / fetal / not-detected style markers in the source list) -> bare source class; the qualifier goes to source_qualifier;
  #       (2) after semicolon split keep only AD/BM/UC; non-class controls (PL/FIB/HSPC etc.) go to the qualifier column;
  #       (3) a side left with no class token after filtering -> no_class_projection (testability only, no pairwise projection)
  # semicolon protection: semicolons inside parentheses become \u0001 so that "adult(BM;AD)" is not torn apart by top-level splitting
  prot_semi <- function(s) {
    s <- gsub("\uff1b", ";", gsub("\uff08", "(", gsub("\uff09", ")", s)))
    ch <- strsplit(s, "", fixed = TRUE)[[1]]; depth <- 0L
    for (i in seq_along(ch)) {
      if (ch[i] == "(") depth <- depth + 1L
      if (ch[i] == ")") depth <- max(depth - 1L, 0L)
      if (ch[i] == ";" && depth > 0L) ch[i] <- "\u0001"
    }
    paste(ch, collapse = "")
  }
  std_side <- function(s) {
    toks <- trimws(strsplit(prot_semi(s), ";", fixed = TRUE)[[1]])
    classes <- c(); quals <- c()
    for (t in toks) {
      pm <- regmatches(t, regexpr("\\([^)]*\\)", t))
      inner <- if (length(pm) == 1 && nzchar(pm))
        trimws(strsplit(gsub("\u0001", ";", substr(pm, 2, nchar(pm) - 1), fixed = TRUE),
                        ";", fixed = TRUE)[[1]]) else character(0)
      inner_cl <- inner[inner %in% CFG$classes]
      bare <- trimws(gsub("\\([^)]*\\)", "", t))
      if (nzchar(bare) && bare %in% CFG$classes) {
        classes <- c(classes, bare)                       # pattern: BM(not detected) -> class BM, qualifier not detected
        quals <- c(quals, inner[!inner %in% CFG$classes])
      } else if (length(inner_cl) > 0) {
        classes <- c(classes, inner_cl)                    # pattern: adult(BM;AD) -> class BM;AD, qualifier adult
        quals <- c(quals, c(bare[nzchar(bare)], inner[!inner %in% CFG$classes]))
      } else {
        quals <- c(quals, t)                               # bare qualifier / non-class control (adult / perinatal / HSPC...)
      }
    }
    list(classes = unique(classes), qual = paste(unique(quals[nzchar(quals)]), collapse = ";"))
  }
  quals <- hi <- lo <- vector("list", nrow(ml))
  for (i in seq_len(nrow(ml))) {
    h <- std_side(ml$high_source[i]); l <- std_side(ml$low_source[i])
    hi[[i]] <- paste(h$classes, collapse = ";"); lo[[i]] <- paste(l$classes, collapse = ";")
    quals[[i]] <- paste(c(h$qual, l$qual), collapse = " | ")
  }
  ml$high_source_std <- unlist(hi); ml$low_source_std <- unlist(lo)
  ml$source_qualifier <- unlist(quals)

  # ---- verification flags (direction semantics = one-vs-rest; pairwise projection is done by expand_claims in script 21)----------
  ml$pmid_ok <- (is.na(ml$pmid) | !nzchar(ml$pmid)) |
                grepl("^[0-9]{6,8}(;[0-9]{6,8})*$", ml$pmid)
  no_proj <- !nzchar(ml$high_source_std) | !nzchar(ml$low_source_std)
  ml$v2_flag <- ifelse(is.na(ml$hgnc_symbol), "symbol_unmapped",
                ifelse(!ml$pmid_ok, "pmid_format",
                ifelse(no_proj, "no_class_projection",
                ifelse(grepl(";", ml$discovery_gse), "discovery_ambiguous", ""))))
  hl$hox_member_ok <- grepl("^HOX[A-D][0-9]{1,2}$", hl$hgnc_symbol)
  hl$v2_flag <- ifelse(is.na(hl$hgnc_symbol), "symbol_unmapped",
                ifelse(!hl$hox_member_ok, "not_hox_pattern", ""))

  out_ml <- ml[, c("claim_id", "gene", "high_source_std", "low_source_std", "discovery_gse", "pmid",
                   "hgnc_symbol", "source_qualifier", "evidence", "flag", "v2_flag")]
  names(out_ml)[match(c("high_source_std", "low_source_std"), names(out_ml))] <- c("high_source", "low_source")
  write_tsv(out_ml, "curation/marker_list.tsv")
  write_tsv(data.frame(gene = hl$hgnc_symbol, hgnc_id = hl$hgnc_id,
                       v2_flag = hl$v2_flag, stringsAsFactors = FALSE),
            "curation/hox_list.tsv")

  log_msg("marker claims: 82 | symbols unmapped ", sum(is.na(ml$hgnc_symbol)),
          " | discovery located ", sum(!is.na(ml$discovery_gse) & !grepl(";", ml$discovery_gse)),
          " (multiple hits pending adjudication ", sum(grepl(";", ml$discovery_gse)), ")",
          " | other flags ", sum(!ml$v2_flag %in% c("", "discovery_ambiguous")))
  log_msg("HOX: 39 | symbols unmapped ", sum(is.na(hl$hgnc_symbol)),
          " | non-HOX naming ", sum(hl$v2_flag == "not_hox_pattern"))
  invisible(ml)
})
