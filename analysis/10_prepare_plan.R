# =============================================================================
# 10_prepare_plan.R —— build the download plan from the earlier plan and
# the v2 layer allocation.
# Rules: only development- and supplementary-layer GSEs are collected
# (chosen = TRUE files, series matrices and selected supplementary files);
# validation-layer data are not downloaded or opened before the analysis plan
# is finalised; earlier files are copied to the v2 layout for 10_download_verify.R.
# Inputs: data/download_plan_v1.tsv (earlier phase); data/layer_allocation_v2.tsv
# Output: data/download_plan_v2.tsv (gse, url, dest_rel) + file copies
# =============================================================================
source("scripts/00_setup.R")

run_step("10_prepare_plan", function() {
  v1_dir <- file.path(normalizePath(file.path(getOption("msc.root"), "..", "..")),
                      "\u7814\u7a76\u65b9\u54112")
  vp_path <- file.path(v1_dir, "data", "download_plan_v1.tsv")
  check(file.exists(vp_path), "v1 download plan not found: ", vp_path)
  vp <- read.delim(vp_path, sep = "\t", quote = "", stringsAsFactors = FALSE,
                   check.names = FALSE, na.strings = "NA")
  check(all(c("gse", "url", "file", "chosen") %in% names(vp)), "v1 plan is missing columns")

  alloc <- read_tsv("data/layer_allocation_v2.tsv")
  want <- alloc$gse[alloc$layer %in% c("D", "supp")]
  check(length(want) == 21, "D+supp should be 21 GSEs, actual ", length(want))

  rows <- vp[vp$gse %in% want & !is.na(vp$chosen) & (vp$chosen == TRUE | vp$chosen == "TRUE"), , drop = FALSE]
  check(nrow(rows) >= length(want), "fewer chosen rows than GSEs; check the chosen column of the v1 plan")
  miss <- setdiff(want, rows$gse)
  check(length(miss) == 0, "no chosen file in the v1 plan for these GSEs: ", paste(miss, collapse = ","))

  plan <- data.frame(gse = rows$gse, url = rows$url,
                     dest_rel = file.path("data", "raw", rows$gse, rows$file),
                     stringsAsFactors = FALSE)
  check(!anyDuplicated(plan$dest_rel), "duplicate destination paths in the download plan")
  # copy v1 files into v2 (no overwrite; missing files are left to script 10 for online download)
  ok <- 0L
  for (i in seq_len(nrow(plan))) {
    src <- file.path(v1_dir, "data_raw", plan$gse[i], basename(plan$dest_rel[i]))
    dst <- proj_path(plan$dest_rel[i])
    if (file.exists(src)) {
      if (!file.exists(dst)) {
        dir.create(dirname(dst), showWarnings = FALSE, recursive = TRUE)
        if (file.copy(src, dst)) ok <- ok + 1L
      } else ok <- ok + 1L
    }
  }
  log_msg("v1 files reusable (pending verification): ", ok, "/", nrow(plan))
  write_tsv(plan, "data/download_plan_v2.tsv")
  invisible(plan)
})
