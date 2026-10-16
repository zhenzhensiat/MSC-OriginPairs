# =============================================================================
# 11_build_matrices.R —— port the earlier verified matrices into v2 study
# objects (19 ported studies).
# Route: the earlier symbol-level matrices carry the per-dataset fixes
# (sample mapping, arm filtering); this script ports expr+labels and adds the
# v2 standardisation layer:
#   1) current HGNC symbol mapping (org.Hs.eg.db 3.21.0, versioned offline source)
#   2) metadata standardisation (passage / medium / age / sex / site parsed from
#      the earlier census sample characteristics)
#   3) coverage audit (gene counts; <10,000 flagged as non-genome-wide) and
#      sex checks (XIST + Y genes)
# Inputs (earlier phase, read-only): data/source_v1/expr/<GSE>_symbol_v1.rds,
#   census_samples_v1.tsv, l2_samples_v1.tsv
# Outputs: data/studies/<GSE>.rds (validate_study-compliant),
#   data/coverage_audit_v2.tsv, data/sample_metadata_v2.tsv; the counts slot is
#   wired by 11c (for voom). GSE139073/GSE309666 are built by 11b instead.
# =============================================================================
source("scripts/00_setup.R")
need_pkgs(c("org.Hs.eg.db", "AnnotationDbi", "edgeR"))
library(org.Hs.eg.db)   # select() requires the package object to be attached

run_step("11_build_matrices", function() {
  v1_dir <- file.path(normalizePath(file.path(getOption("msc.root"), "..", "..")),
                      "\u7814\u7a76\u65b9\u54112")

  # ---- 1. v1 input snapshot (matrices and sample tables)---------------------------------
  dir.create(proj_path("data", "source_v1", "expr"), showWarnings = FALSE)
  v1_sample_tab <- file.path(v1_dir, "data", c("census_samples_v1.tsv", "l2_samples_v1.tsv"))
  check(all(file.exists(v1_sample_tab)), "v1 sample tables missing")
  for (p in v1_sample_tab) {
    dst <- file.path("data", "source_v1", basename(p))
    if (!file.exists(proj_path(dst))) file.copy(p, proj_path("data", "source_v1"))
  }
  cs <- read_tsv("data/source_v1/census_samples_v1.tsv")
  l2 <- read_tsv("data/source_v1/l2_samples_v1.tsv")

  alloc <- read_tsv("data/layer_allocation_v2.tsv")
  geo <- read_tsv("data/census_geo_v2.tsv")
  portable <- alloc$gse[alloc$layer %in% c("D", "supp") &
                        file.exists(file.path(v1_dir, "data", "expr",
                                              paste0(alloc$gse, "_symbol_v1.rds")))]
  check(length(portable) == 19, "portable studies should be 19, actual ", length(portable),
        ": ", paste(portable, collapse = ","))
  for (g in portable) {
    dst <- file.path("data", "source_v1", "expr", paste0(g, "_symbol_v1.rds"))
    if (!file.exists(proj_path(dst)))
      file.copy(file.path(v1_dir, "data", "expr", basename(dst)), proj_path(dst))
  }

  # ---- 2. HGNC current symbol mapping (org.Hs.eg.db: SYMBOL direct hit + unique ALIAS mapping)------
  eg <- AnnotationDbi::select(org.Hs.eg.db, keys = keys(org.Hs.eg.db),
                              columns = c("SYMBOL", "ALIAS"))
  sym_ok <- unique(na.omit(eg$SYMBOL))
  al <- eg[!is.na(eg$ALIAS) & !is.na(eg$SYMBOL), c("ALIAS", "SYMBOL")]
  al_n <- table(al$ALIAS)
  uniq_alias <- names(al_n[al_n == 1L])                      # unique-mapping aliases only
  alias_map <- setNames(al$SYMBOL[match(uniq_alias, al$ALIAS)], uniq_alias)
  log_msg("HGNC mapping: current symbols ", length(sym_ok), ", unique aliases ", length(alias_map),
          "(org.Hs.eg.db ", as.character(packageVersion("org.Hs.eg.db")), ")")

  map_genes <- function(genes) {
    out <- ifelse(genes %in% sym_ok, genes,
           ifelse(genes %in% names(alias_map), alias_map[genes], NA_character_))
    out[is.na(out)] <- NA_character_
    out
  }
  collapse_dup <- function(x) {   # duplicate symbols: highest mean kept, ties broken by IQR
    if (!anyDuplicated(rownames(x))) return(x)
    mu <- rowMeans(x); iqr <- apply(x, 1, stats::IQR)
    ord <- order(mu, iqr, decreasing = TRUE)
    keep <- ord[!duplicated(rownames(x)[ord])]
    x[keep, , drop = FALSE]
  }

  # ---- 3. characteristic text parsing (passage / medium / age / sex)---------------------------
  parse_passage <- function(txt) {
    hit <- grepl("passage[^0-9]*([0-9]+)", txt, ignore.case = TRUE)
    vapply(seq_along(txt), function(i) {
      if (!hit[i]) return(NA_real_)
      as.numeric(sub(".*passage[^0-9]*([0-9]+).*", "\\1", txt[i], ignore.case = TRUE))
    }, numeric(1))
  }
  parse_medium <- function(txt) {
    t <- tolower(txt)
    ifelse(grepl("platelet lysate|\\bhpl\\b", t), "hPL",
    ifelse(grepl("fetal bovine|\\bfbs\\b", t) & !grepl("platelet", t), "FBS",
    ifelse(grepl("human serum|\\bhs\\b", t) & !grepl("fetal", t), "HS",
    ifelse(grepl("serum.?free|xeno.?free|\\bxf\\b", t), "xeno_free", NA_character_))))
  }
  # S2 fix (appended): v1 age format has no y suffix ("age 5"), relaxed
  parse_age <- function(txt) {
    vapply(txt, function(s) {
      m <- regmatches(s, regexpr("(?i)(age|donor.?age)[^0-9]{0,8}([0-9]{1,2})", s, perl = TRUE))
      if (length(m) == 0) return(NA_real_)
      as.numeric(sub(".*?([0-9]{1,2})$", "\1", m[1]))
    }, numeric(1))
  }
  # S2 fix (appended): v1 characteristic text often gives sex as a single letter (gender: f / sex: m); abbreviation pattern added
  parse_sex <- function(txt) {
    t <- tolower(txt)
    f_hit <- grepl("female", t) | grepl("(gender|sex)[ :=]+f($|[^a-z])", t, perl = TRUE)
    m_hit <- (grepl("male", t) & !grepl("female", t)) | grepl("(gender|sex)[ :=]+m($|[^a-z])", t, perl = TRUE)
    ifelse(f_hit & !m_hit, "F", ifelse(m_hit & !f_hit, "M", NA_character_))
  }

  y_genes <- c("DDX3Y", "RPS4Y1", "KDM5D", "EIF1AY", "USP9Y", "ZFY", "UTY")

  # ---- 4. per-study porting ------------------------------------------------------------
  dir.create(proj_path("data", "studies"), showWarnings = FALSE)
  audit <- list(); allmeta <- list()
  for (g in portable) {
    x <- readRDS(proj_path("data", "source_v1", "expr", paste0(g, "_symbol_v1.rds")))
    lb <- x$labels
    attr_note <- ""
    check(identical(colnames(x$expr), lb$sample), g, ": v1 expr columns not aligned with labels samples")

    genes <- map_genes(rownames(x$expr))
    hit_rate <- mean(!is.na(genes))
    e0 <- x$expr[!is.na(genes), , drop = FALSE]
    rownames(e0) <- genes[!is.na(genes)]
    e0 <- collapse_dup(e0)

    # scale detection: v1 RDS scales are not uniform (counts / linear array intensity / author normalisation coexist)
    dat_type <- geo$data_type[geo$gse == g]
    is_array <- dat_type == "array"
    counts_m <- NULL
    if (max(e0) >= CFG$log_scale_max) {
      is_counts <- all(e0 >= 0) && mean(abs(e0 - round(e0)) < 1e-8) > 0.99
      if (is_counts) {
        counts_m <- e0
        e <- edgeR::cpm(e0, log = TRUE, prior.count = 0.5)   # analysis protocol Phase 1.3
        scale <- "logCPM"
        attr_note <- paste0(attr_note, ifelse(nzchar(attr_note), "; ", ""), "v1 stored counts; cpm(log) applied here")
      } else if (is_array) {
        e <- log2(e0 + 1)                                     # linear array intensity: log2, cpm forbidden
        scale <- "array"
        attr_note <- paste0(attr_note, ifelse(nzchar(attr_note), "; ", ""), "v1 stored linear intensity; log2(x+1)")
      } else {
        e <- log2(e0 + 1)                                     # non-integer RNA-seq matrix: flagged by the TPM rule
        scale <- "log2TPM"
        attr_note <- paste0(attr_note, ifelse(nzchar(attr_note), "; ", ""),
                           "v1 stored a non-integer matrix (integer rate<0.99); log2(x+1) flagged log2TPM")
      }
    } else {
      e <- e0
      scale <- if (grepl("^scRNA|RNA-seq", dat_type)) "logCPM" else "array"
      if (g == "GSE145728") attr_note <- "author_normCPM (no raw counts, voom not applicable)"
    }
    check(!anyNA(e), g, ": expr contains NA after normalisation")
    check(max(e) < CFG$log_scale_max, g, ": still looks unlogged after transformation, max = ", round(max(e), 1))

    # characteristic text: census_samples joined on gsm with title+characteristics;
    # join-rate assertion >=0.9; the three l2 studies (GSE12274/51869/9451) lack characteristics; use note text and register it
    # S2b: when labels$sample is not a GSM (author-internal IDs), map via series matrix headers (normalisation -> order+suffix check)
    smp_ids <- lb$sample
    if (!all(startsWith(smp_ids, "GSM")) && any(cs$gse == g)) {
      smf <- list.files(proj_path("data", "raw", g), "series_matrix.*gz",
                        recursive = TRUE, full.names = TRUE)
      check(length(smf) >= 1, g, ": sample names are not GSM and no series matrix available for mapping")
      # merge all series matrices of the study (title-union matching across platform subseries)
      tg <- do.call(rbind, lapply(smf, function(f1) {
        hh <- read_series_header(f1)
        data.frame(title = hh[[grep("^Sample_title$", names(hh))[1]]],
                   gsm = hh[[grep("^Sample_geo_accession$", names(hh))[1]]],
                   stringsAsFactors = FALSE)
      }))
      tg <- tg[!duplicated(tg$title) & !is.na(tg$title), , drop = FALSE]
      m <- match(tolower(gsub("[^a-z0-9]", "", smp_ids)),
                 tolower(gsub("[^a-z0-9]", "", tg$title)))
      # third level: prefix-key matching (first two underscore segments, absorbing mid-string title variants like _Lipo/_MA)
      pref <- function(v) tolower(gsub("[^a-z0-9]", "",
                     vapply(strsplit(v, "_", fixed = TRUE), function(z) paste(z[1:2], collapse = "_"), "")))
      if (mean(!is.na(m)) < 0.9) {
        m2 <- match(pref(smp_ids), pref(tg$title))
        m[is.na(m)] <- m2[is.na(m)]
        if (mean(!is.na(m)) >= 0.9) log_msg(g, ": sample names mapped to GSM via prefix keys")
      }
      if (mean(!is.na(m)) < 0.9 && length(smf) == 1) {
        tt <- tg$title
        check(length(tt) == length(smp_ids), g, ": order mapping length mismatch (", length(tt), " vs ", length(smp_ids), ")")
        tail_ok <- mapply(function(a, b) sub(".*_", "", a) == sub(".*_", "", b), smp_ids, tt)
        check(mean(tail_ok) >= 0.9, g, ": order mapping suffix check failed (", round(mean(tail_ok), 2), ")")
        m <- seq_along(smp_ids)
        log_msg(g, ": sample names mapped to GSM via order+suffix check")
      } else check(mean(!is.na(m)) >= 0.9, g, ": title match rate ", round(mean(!is.na(m)), 2), ", mapping failed")
      if (mean(!is.na(match(  # only when title matching is actually used
          tolower(gsub("[^a-z0-9]", "", smp_ids)),
          tolower(gsub("[^a-z0-9]", "", tg$title))))) >= 0.9)
        log_msg(g, ": sample names mapped to GSM via normalised title match")
      smp_ids <- ifelse(is.na(m), smp_ids, tg$gsm[m])
    }
    if (any(cs$gse == g)) {
      tx <- cs$characteristics[match(smp_ids, cs$gsm)]
      ti <- cs$title[match(smp_ids, cs$gsm)]
      check(mean(!is.na(tx)) >= 0.9, g, ": characteristic-text join rate ", round(mean(!is.na(tx)), 2), ", check the join key")
    } else {
      tx <- l2$note[match(lb$sample, l2$sample)]
      ti <- rep(NA_character_, nrow(lb))
      log_msg(g, ": l2 source has no characteristics, using note text")
    }
    txt <- paste(ifelse(is.na(ti), "", ti), ifelse(is.na(tx), "", tx), sep = " | ")

    medium <- parse_medium(txt)
    medium[is.na(medium) & g == "GSE51869"] <- "heterogeneous(FBS/HS/hPL,v1 record)"
    s <- list(
      gse = g,
      cluster_id = alloc$cluster_id[alloc$gse == g],
      layer = alloc$layer[alloc$gse == g],
      stability = isTRUE(geo$stability[geo$gse == g]),
      full_genome = isTRUE(geo$full_genome[geo$gse == g]) && nrow(e) >= CFG$genome_min_genes,
      scale = scale,
      expr = e,
      counts = counts_m,                   # linear counts wired here; the rest filled in by 11c
      meta = data.frame(
        # default donor code = source:index; shared codes kept only where the paper states same-donor simultaneous sampling
        # shared-code whitelist (per v1 census/paper): GSE57151 (same-donor three-tissue series),
        # GSE122778 (census note "14BM+14AD autologous pairs")
        donor_shared = g %in% c("GSE57151", "GSE122778"),
        gsm = lb$sample, source = lb$source,
        donor = if (g %in% c("GSE57151", "GSE122778")) as.character(lb$donor)
                else paste(lb$source, lb$donor, sep = ":"),
        donor_raw = as.character(lb$donor),
        passage = parse_passage(txt), medium = medium,
        age = parse_age(txt), sex = parse_sex(txt),
        site = geo$uc_site[geo$gse == g],
        arm = lb$arm, stringsAsFactors = FALSE))
    # sex check (inferred values are compared in script-12 QC)
    yg <- rownames(e)[rownames(e) %in% y_genes]
    s$meta$ygene_mean <- if (length(yg)) colMeans(e[yg, , drop = FALSE]) else NA_real_
    s$meta$xist <- if ("XIST" %in% rownames(e)) e["XIST", ] else NA_real_
    validate_study(s)
    saveRDS(s, proj_path("data", "studies", paste0(g, ".rds")))

    audit[[length(audit) + 1L]] <- data.frame(
      gse = g, n_genes_raw = nrow(x$expr), n_genes_mapped = nrow(e),
      annot_rate = round(hit_rate, 4), n_samples = ncol(e),
      n_include = sum(s$meta$arm == "include"), genome_wide = s$full_genome,
      scale = scale, notes = attr_note, stringsAsFactors = FALSE)
    mm <- s$meta; mm$study <- g; mm$cluster_id <- s$cluster_id; mm$layer <- s$layer
    allmeta[[length(allmeta) + 1L]] <- mm
    log_msg(g, ": ", nrow(e), " genes / ", ncol(e), " samples (include ",
            sum(s$meta$arm == "include"), "), annotation rate ", round(hit_rate, 3))
  }

  aud <- do.call(rbind, audit)
  low <- aud$annot_rate < CFG$annot_hit_min
  if (any(low)) log_msg("!! annotation hit rate below ", CFG$annot_hit_min, ": ",
                        paste(aud$gse[low], collapse = ","))
  write_tsv(aud, "data/coverage_audit_v2.tsv")
  write_tsv(do.call(rbind, allmeta), "data/sample_metadata_v2.tsv")
  invisible(aud)
})
