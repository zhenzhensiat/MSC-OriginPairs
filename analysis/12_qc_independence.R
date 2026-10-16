# =============================================================================
# 12_qc_independence.R —— sample QC and independence checks.
# (1) sex check: stated vs XIST/Y-inferred (mismatches flagged, never deleted)
# (2) within-study technical replicates: correlation >0.99 with matching
#     metadata -> flagged as merge candidates
# (3) suspected same cells across studies: cross-study correlation >0.98 ->
#     flagged for independence review
# (4) independence review table: v2 clusters + cross-study evidence
# Inputs: data/studies/*.rds
# Outputs: data/sample_qc_v2.tsv, data/independence_clusters_v2.tsv
# Discipline: flag only, never delete.
# =============================================================================
source("scripts/00_setup.R")

run_step("12_qc_independence", function() {
  st <- load_studies()                       # 19 studies (D+supp)
  check(length(st) >= 15, "unexpected study count: ", length(st))

  # ---- 1. sex check --------------------------------------------------------------
  # within-study z scores misclassify single-sex studies; replaced by an absolute Y three-gene threshold,
  # thresholds calibrated per scale (array / logCPM) on declared-sex samples (midpoint of the two-sex Y means)
  infer_sex <- function(ymean, thr) ifelse(is.na(ymean), NA_character_, ifelse(ymean > thr, "M", "F"))
  y3 <- c("RPS4Y1", "DDX3Y", "KDM5D")
  thr_for <- function(st) {
    yv <- vapply(st, function(s) {
      if (all(y3 %in% rownames(s$expr))) mean(s$expr[y3, ]) else NA_real_
    }, numeric(1))
    names(yv) <- vapply(st, `[[`, "", "gse")
    cal <- do.call(rbind, lapply(st, function(s) {
      ok <- !is.na(s$meta$sex) & !is.na(yv[[s$gse]])
      if (!any(ok)) return(NULL)
      data.frame(gse = s$gse, scale = s$scale, sex = s$meta$sex[ok], y = rep(yv[[s$gse]], sum(ok)))
    }))
    if (is.null(cal)) return(list(thr = c(array = 8.5, logCPM = 2), cal = NULL, yv = yv))
    tr <- tapply(cal$y, cal$scale, function(v) v)   # placeholder
    thr <- vapply(split(cal, cal$scale), function(d) {
      m <- tapply(d$y, d$sex, mean)
      if (length(m) == 2) mean(m) else NA_real_
    }, numeric(1))
    thr[is.na(thr)] <- ifelse(names(thr[is.na(thr)]) == "array", 8.5, 2)
    # 2026-10-12 fix (sponsor-approved): values must be fully unnamed -- is.na(thr["array"]) carries the name "array",
    # ifelse inherits that name, so c(array=...) stacks into array.array -> the threshold stays NA -> sex QC fails
    pick_thr <- function(nm, dflt) {
      x <- if (nm %in% names(thr)) unname(thr[[nm]]) else NA_real_
      if (is.na(x)) dflt else x
    }
    list(thr = c(array = pick_thr("array", 8.5), logCPM = pick_thr("logCPM", 2)), cal = cal, yv = yv)
  }
  tf <- thr_for(st)
  log_msg("sex thresholds (Y three-gene means, declared-calibrated): array=", round(tf$thr["array"], 2),
          " logCPM=", round(tf$thr["logCPM"], 2),
          " (calibration samples ", if (is.null(tf$cal)) 0 else nrow(tf$cal), ")")
  qc <- do.call(rbind, lapply(st, function(s) {
    m <- s$meta
    yv <- if (all(y3 %in% rownames(s$expr))) colMeans(s$expr[y3, , drop = FALSE]) else rep(NA_real_, nrow(m))
    inf <- infer_sex(yv, unname(tf$thr[s$scale]))
    data.frame(gsm = m$gsm, study = s$gse, source = m$source, donor = m$donor,
               sex_declared = m$sex, sex_inferred = inf, y_mean = round(yv, 2),
               sex_flag = ifelse(!is.na(m$sex) & !is.na(inf) & m$sex != inf, "mismatch", ""),
               stringsAsFactors = FALSE)
  }))

  # ---- 2. within-study technical replicates (cor>0.99 + same donor and source)-------------------------------
  techrep <- list()
  for (s in st) {
    m <- s$meta[s$meta$arm == "include" & s$meta$source %in% CFG$classes, , drop = FALSE]
    if (nrow(m) < 2) next
    cm <- stats::cor(s$expr[, m$gsm], method = "spearman")
    hi <- which(upper.tri(cm) & cm > 0.99, arr.ind = TRUE)
    if (nrow(hi) == 0) next
    for (k in seq_len(nrow(hi))) {
      i <- hi[k, 1]; j <- hi[k, 2]
      if (m$donor[i] == m$donor[j] && m$source[i] == m$source[j])
        techrep[[length(techrep) + 1L]] <- data.frame(
          study = s$gse, gsm1 = m$gsm[i], gsm2 = m$gsm[j], donor = m$donor[i],
          rho = round(cm[i, j], 5), stringsAsFactors = FALSE)
    }
  }
  tr <- if (length(techrep)) do.call(rbind, techrep) else
    data.frame(study = character(), gsm1 = character(), gsm2 = character(),
               donor = character(), rho = numeric())
  # map merge candidates onto the qc table
  qc$max_within_rho <- NA_real_
  for (i in seq_len(nrow(tr))) {
    qc$max_within_rho[qc$gsm == tr$gsm1[i]] <- pmax(qc$max_within_rho[qc$gsm == tr$gsm1[i]], tr$rho[i], na.rm = TRUE)[1]
    qc$max_within_rho[qc$gsm == tr$gsm2[i]] <- pmax(qc$max_within_rho[qc$gsm == tr$gsm2[i]], tr$rho[i], na.rm = TRUE)[1]
  }
  qc$techrep_flag <- ifelse(!is.na(qc$max_within_rho) & qc$max_within_rho > 0.99, "techrep_candidate", "")

  # ---- 3. suspected same cells across studies (cross-study pairs with cor>0.98)--------------------------------
  common <- Reduce(intersect, lapply(st, function(s) rownames(s$expr)))
  # intersection lower bound aligned to the training_view value 2000 (old Illumina symbol sets shrink the intersection; 2370 measured after the M2 whitelist)
  check(length(common) >= 2000, "cross-study gene intersection only ", length(common))
  # cross-platform Spearman cannot reach 0.98 even for the same cell -- compare within one GPL platform only
  geo0 <- read_tsv("data/census_geo_v2.tsv")
  gpl_of <- setNames(as.character(geo0$gpl), geo0$gse)
  big <- do.call(cbind, lapply(st, function(s) {
    m <- s$meta[s$meta$arm == "include" & s$meta$source %in% CFG$classes, , drop = FALSE]
    s$expr[common, m$gsm, drop = FALSE]
  }))
  lab <- do.call(rbind, lapply(st, function(s) {
    m <- s$meta[s$meta$arm == "include" & s$meta$source %in% CFG$classes, , drop = FALSE]
    data.frame(gsm = m$gsm, study = s$gse, donor = m$donor, source = m$source,
               stringsAsFactors = FALSE, row.names = m$gsm)
  }))
  check(identical(colnames(big), lab$gsm))
  cross <- data.frame(gsm1 = character(), gsm2 = character(), rho = numeric(),
                      study1 = character(), study2 = character(), stringsAsFactors = FALSE)
  for (gp in unique(na.omit(gpl_of[lab$study]))) {
    keep <- lab$gsm[lab$study %in% names(gpl_of)[gpl_of == gp & !is.na(gpl_of)]]
    keep <- intersect(keep, colnames(big))
    st_in <- unique(lab$study[lab$gsm %in% keep])
    if (length(st_in) < 2) next
    cm <- stats::cor(big[, keep, drop = FALSE], method = "spearman")
    hi <- which(upper.tri(cm) & cm > 0.98, arr.ind = TRUE)
    if (nrow(hi) == 0) next
    d <- data.frame(gsm1 = colnames(cm)[hi[, 1]], gsm2 = colnames(cm)[hi[, 2]],
                    rho = round(cm[hi], 5), stringsAsFactors = FALSE)
    d$study1 <- lab$study[match(d$gsm1, lab$gsm)]
    d$study2 <- lab$study[match(d$gsm2, lab$gsm)]
    d <- d[d$study1 != d$study2, , drop = FALSE]
    cross <- rbind(cross, d)
  }
  hi <- which(upper.tri(cm) & cm > 0.98, arr.ind = TRUE)
  qc$cross_study_rho <- NA_real_
  if (nrow(cross) > 0) {
    for (i in seq_len(nrow(cross))) {
      qc$cross_study_rho[qc$gsm == cross$gsm1[i]] <- pmax(qc$cross_study_rho[qc$gsm == cross$gsm1[i]], cross$rho[i], na.rm = TRUE)[1]
      qc$cross_study_rho[qc$gsm == cross$gsm2[i]] <- pmax(qc$cross_study_rho[qc$gsm == cross$gsm2[i]], cross$rho[i], na.rm = TRUE)[1]
    }
  }
  qc$same_cell_flag <- ifelse(!is.na(qc$cross_study_rho) & qc$cross_study_rho > 0.98, "cross_study_suspect", "")
  write_tsv(qc, "data/sample_qc_v2.tsv")

  # ---- 4. independence-cluster review table ----------------------------------------------------------
  geo <- read_tsv("data/census_geo_v2.tsv")
  alloc <- read_tsv("data/layer_allocation_v2.tsv")
  alloc$same_cell_pairs <- vapply(alloc$gse, function(g)
    sum((cross$study1 == g | cross$study2 == g)), integer(1))
  alloc$techrep_pairs <- as.integer(table(tr$study)[alloc$gse])
  alloc$techrep_pairs[is.na(alloc$techrep_pairs)] <- 0L
  colnames(alloc)[which(names(alloc) == "layer")] <- "layer"
  write_tsv(alloc, "data/independence_clusters_v2.tsv")

  log_msg("sex check: declared and inferable ", sum(!is.na(qc$sex_declared) & !is.na(qc$sex_inferred)),
          ", mismatched ", sum(qc$sex_flag == "mismatch"))
  log_msg("technical-replicate candidate pairs: ", nrow(tr), "; same-platform cross-study suspected same-cell pairs (M4 wording): ", nrow(cross))
  if (nrow(cross) > 0) log_msg("studies involved in cross-study pairs: ",
    paste(sort(unique(c(cross$study1, cross$study2))), collapse = ","))
  invisible(qc)
})
