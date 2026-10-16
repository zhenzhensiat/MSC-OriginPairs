# =============================================================================
# 10_download_verify.R —— integrity-checked downloads.
# Input: data/download_plan_v2.tsv (gse, url, dest_rel; paths relative to root)
# Output: data/raw_manifest.tsv (one row per planned file, failures included)
# Files that already exist are verified first and reused only if they pass,
# otherwise re-downloaded.
# =============================================================================
source("scripts/00_setup.R")
need_pkgs(c("httr", "curl"))

run_step("10_download_verify", function() {
  plan <- read_tsv("data/download_plan_v2.tsv")
  check(all(c("gse", "url", "dest_rel") %in% names(plan)), "download plan is missing columns gse/url/dest_rel")
  check(!anyDuplicated(plan$dest_rel), "duplicate destination paths in the download plan")

  rows <- lapply(seq_len(nrow(plan)), function(i) {
    dest <- proj_path(plan$dest_rel[i])
    if (file.exists(dest)) {
      v <- verify_existing(plan$url[i], dest)
      if (v$status == "ok_reuse") {
        return(data.frame(gse = plan$gse[i], url = plan$url[i], dest = plan$dest_rel[i],
                          expected_bytes = v$expected_bytes, bytes = v$bytes, gzip_ok = v$gzip_ok,
                          tries = 0L, status = "ok_reuse", note = "", stringsAsFactors = FALSE))
      }
      log_msg("existing file failed verification; re-downloading: ", plan$dest_rel[i])
    }
    r <- download_verified(plan$url[i], dest)
    r$dest <- plan$dest_rel[i]
    cbind(gse = plan$gse[i], r[, c("url", "dest", "expected_bytes", "bytes", "gzip_ok", "tries", "status", "note")])
  })
  man <- do.call(rbind, rows)
  expect_rows(man, nrow(plan), "download manifest")
  write_tsv(man, "data/raw_manifest.tsv")
  n_bad <- sum(!man$status %in% c("ok", "ok_reuse"))
  check(n_bad == 0, n_bad, " files failed download or verification; see data/raw_manifest.tsv; fix before proceeding")
  invisible(man)
})
