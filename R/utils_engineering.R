# =============================================================================
# utils_engineering.R —— engineering conventions as code (paths, logging, checks,
# deterministic checks, safe merges). Run tests/run_tests.R before any real-data script.
# =============================================================================

# ---- Paths: locate the project root via a marker file, absolute drive letters forbidden ----------------------------------
find_root <- function(start = getwd(), marker = ".msc_originpairs_root") {
  d <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(d, marker))) return(d)
    up <- dirname(d)
    if (identical(up, d)) stop("project-root marker file not found: ", marker, "; create an empty file at the project root")
    d <- up
  }
}

proj_path <- function(...) {
  root <- getOption("msc.root")
  if (is.null(root)) stop("msc.root is not set; run analysis/00_setup.R first")
  file.path(root, ...)
}

# ---- Logging (timestamps go to logs only, never into result tables)-------------------------------------
log_msg <- function(..., log_file = getOption("msc.log_file")) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
  message(line)
  if (!is.null(log_file)) cat(line, "\n", file = log_file, append = TRUE, sep = "")
  invisible(line)
}

step_log_append <- function(step, status, note = "", file = proj_path("logs", "run_log.tsv")) {
  row <- data.frame(time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), step = step,
                    status = status, note = gsub("[\t\r\n]+", " ", note),
                    stringsAsFactors = FALSE)
  is_new <- !file.exists(file)
  utils::write.table(row, file, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = is_new, append = !is_new)
  invisible(row)
}

# Top-level wrapper: 'finished' is written only on normal completion; the handler just returns a value and the caller records it
run_step <- function(step, fun) {
  step_log_append(step, "started")
  res <- tryCatch(list(ok = TRUE, value = fun()),
                  error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
  if (isTRUE(res$ok)) {
    step_log_append(step, "finished")
    return(invisible(res$value))
  }
  step_log_append(step, "failed", res$msg)
  stop("[", step, "] failed: ", res$msg, call. = FALSE)
}

# ---- Seeds: each random step uses a derived seed-------------------------------------
# steps are spaced by 1e6 to avoid seed-range overlap between different steps; i must be < 1e6; step ≤ 2000
derive_seed <- function(step, i = 0L) {
  check(i >= 0 && i < 1e6, "derive_seed: i out of range")
  as.integer(CFG$master_seed + as.integer(step) * 1000000L + as.integer(i))
}

# ---- Assertions ---------------------------------------------------------------------
check <- function(cond, ...) {
  if (!isTRUE(cond)) stop(paste0(...), call. = FALSE)
  invisible(TRUE)
}

expect_rows <- function(df, n, what = "table") {
  check(NROW(df) == n, what, " rows = ", NROW(df), ", expected = ", n)
}

# Assert key uniqueness and row counts before and after merging
safe_merge <- function(x, y, by, type = c("m:1", "1:1"), all.x = TRUE) {
  type <- match.arg(type)
  check(all(by %in% names(x)) && all(by %in% names(y)), "merge keys not present in both tables: ", paste(by, collapse = ","))
  check(!anyNA(x[by]) && !anyNA(y[by]), "merge keys contain NA")
  kx <- do.call(paste, c(unname(as.list(x[by])), sep = "\r"))
  ky <- do.call(paste, c(unname(as.list(y[by])), sep = "\r"))
  check(!anyDuplicated(ky), "right-table merge keys not unique, e.g.: ", paste(utils::head(unique(ky[duplicated(ky)]), 5), collapse = " / "))
  if (type == "1:1") {
    check(!anyDuplicated(kx), "left-table merge keys not unique")
    miss <- setdiff(kx, ky)
    check(length(miss) == 0, "1:1 merge: ", length(miss), " left-table keys missing from the right table")
  }
  out <- merge(x, y, by = by, all.x = all.x, sort = FALSE)
  if (all.x) check(nrow(out) == nrow(x), "row count changed after merge: ", nrow(x), " -> ", nrow(out))
  out
}

# ---- TSV I/O (uniform format, no timestamps in result tables)---------------------------------------
write_tsv <- function(df, rel_path) {
  path <- proj_path(rel_path)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(df, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  log_msg("wrote ", rel_path, " (", NROW(df), " rows)")
  invisible(path)
}

read_tsv <- function(rel_path) {
  utils::read.delim(proj_path(rel_path), sep = "\t", quote = "", comment.char = "",
                    stringsAsFactors = FALSE, check.names = FALSE, na.strings = "NA")
}

# ---- Download integrity------------------------------
gzip_ok <- function(path) {
  sz <- file.size(path)
  if (is.na(sz) || sz < 20) return(FALSE)
  con <- file(path, "rb"); magic <- readBin(con, "raw", 2L); close(con)
  if (!identical(magic, as.raw(c(0x1f, 0x8b)))) return(FALSE)

  gz <- Sys.which("gzip")
  if (nzchar(gz)) {
    return(identical(as.integer(system2(gz, c("-t", shQuote(path)), stdout = FALSE, stderr = FALSE)), 0L))
  }
  gzip_ok_r(path)
}

# pure-R fallback: stream-decompress and count bytes, compared against the gzip trailer ISIZE (decompressed length mod 2^32)
# Limitation: for multi-member gzip (concatenated gzip streams) only the last member is checked; GEO files are usually single-member
gzip_ok_r <- function(path) {
  sz <- file.size(path)
  if (is.na(sz) || sz < 20) return(FALSE)
  con <- file(path, "rb"); seek(con, where = sz - 4)
  isize <- readBin(con, "integer", n = 1L, size = 4L, endian = "little"); close(con)
  if (isize < 0) isize <- isize + 2^32
  total <- tryCatch({
    zc <- gzcon(file(path, "rb"))
    n <- 0
    repeat {
      b <- readBin(zc, "raw", n = 8e6)
      if (length(b) == 0) break
      n <- n + length(b)
    }
    close(zc)
    n
  }, error = function(e) NA_real_, warning = function(w) NA_real_)
  isTRUE(!is.na(total) && total > 0 && (total %% 2^32) == isize)
}

remote_size <- function(url) {
  tryCatch({
    h <- httr::HEAD(url, httr::timeout(60))
    cl <- httr::headers(h)[["content-length"]]
    if (is.null(cl)) NA_real_ else as.numeric(cl)
  }, error = function(e) NA_real_)
}

# Returns one manifest row; failures also return a row, recorded by the caller
download_verified <- function(url, dest, max_tries = 3L) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  expected <- remote_size(url)
  is_gz <- grepl("\\.gz$", dest, ignore.case = TRUE)
  status <- "failed"; note <- ""; tries_used <- 0L
  for (k in seq_len(max_tries)) {
    tries_used <- k
    tmp <- paste0(dest, ".part")
    if (file.exists(tmp)) file.remove(tmp)
    err <- tryCatch({ curl::curl_download(url, tmp, quiet = TRUE); "" },
                    error = function(e) conditionMessage(e))
    if (nzchar(err)) { note <- err; next }
    size <- file.size(tmp)
    size_ok <- is.na(expected) || isTRUE(size == expected)
    gz_ok <- if (is_gz) gzip_ok(tmp) else NA
    if (size_ok && !isFALSE(gz_ok)) {
      if (file.exists(dest)) file.remove(dest)
      check(file.rename(tmp, dest), "cannot rename ", tmp)
      status <- "ok"
      note <- if (is.na(expected)) "server gave no Content-Length; gzip validation only" else ""
      break
    }
    note <- sprintf("size %s / expected %s; gzip_ok=%s", size, expected, gz_ok)
  }
  data.frame(url = url, dest = dest, expected_bytes = expected,
             bytes = if (file.exists(dest)) file.size(dest) else NA_real_,
             gzip_ok = if (status == "ok" && is_gz) TRUE else NA,
             tries = tries_used, status = status, note = note, stringsAsFactors = FALSE)
}

# Pre-reuse check for files already downloaded by v1 (remote size + gzip integrity)
verify_existing <- function(url, dest) {
  if (!file.exists(dest)) return(data.frame(dest = dest, status = "missing", stringsAsFactors = FALSE))
  expected <- remote_size(url)
  size_ok <- is.na(expected) || isTRUE(file.size(dest) == expected)
  gz_ok <- if (grepl("\\.gz$", dest, ignore.case = TRUE)) gzip_ok(dest) else NA
  data.frame(dest = dest, expected_bytes = expected, bytes = file.size(dest),
             gzip_ok = gz_ok,
             status = if (size_ok && !isFALSE(gz_ok)) "ok_reuse" else "redownload",
             stringsAsFactors = FALSE)
}

# ---- Soft memory check (hard protection via tools/mem_watchdog.ps1)----------------------------
mem_rss_gb <- function() {
  if (!requireNamespace("ps", quietly = TRUE)) return(NA_real_)
  unname(ps::ps_memory_info(ps::ps_handle())[["rss"]]) / 1024^3
}

mem_guard <- function(where = "") {
  g <- mem_rss_gb()
  if (!is.na(g) && g > CFG$mem_limit_gb * 0.9) {
    stop("memory ", round(g, 1), "GB is near the limit (", where, "), stopping voluntarily", call. = FALSE)
  }
  invisible(g)
}
