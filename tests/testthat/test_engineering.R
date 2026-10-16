# Engineering-utility tests (gzip integrity, safe_merge, run_step, derived seeds).

test_that("gzip integrity: intact files pass; truncated and non-gzip files are flagged", {
  set.seed(11)
  tmp <- tempfile(fileext = ".gz")
  con <- gzfile(tmp, "w")
  writeLines(vapply(1:20000, function(i) paste(sample(c("A", "C", "G", "T"), 80, TRUE), collapse = ""), ""), con)
  close(con)
  expect_true(gzip_ok(tmp))
  expect_true(gzip_ok_r(tmp))

  raw <- readBin(tmp, "raw", file.size(tmp))
  bad <- tempfile(fileext = ".gz")
  writeBin(raw[seq_len(floor(length(raw) * 0.6))], bad)
  expect_false(gzip_ok(bad))
  expect_false(gzip_ok_r(bad))

  txt <- tempfile(fileext = ".gz")
  writeLines("not a gzip file at all, just text", txt)
  expect_false(gzip_ok(txt))
})

test_that("safe_merge reports duplicated right keys, missing 1:1 keys and row-count changes", {
  x <- data.frame(id = c("a", "b", "c"), v = 1:3)
  y <- data.frame(id = c("a", "b", "c"), w = 4:6)
  expect_equal(nrow(safe_merge(x, y, "id")), 3)
  y2 <- rbind(y, data.frame(id = "a", w = 9L))
  expect_error(safe_merge(x, y2, "id"), "not unique")
  expect_error(safe_merge(x, y[1:2, ], "id", type = "1:1"), "missing")
  x_na <- data.frame(id = c("a", NA), v = 1:2)
  expect_error(safe_merge(x_na, y, "id"), "NA")
})

test_that("run_step writes failed on error and finished on success", {
  tmp <- tempfile(); dir.create(file.path(tmp, "logs"), recursive = TRUE)
  op <- options(msc.root = tmp, msc.log_file = NULL)
  expect_error(run_step("t_fail", function() stop("boom")), "boom")
  expect_equal(run_step("t_ok", function() 42), 42)
  man <- utils::read.delim(file.path(tmp, "logs", "run_log.tsv"), stringsAsFactors = FALSE)
  expect_equal(man$status, c("started", "failed", "started", "finished"))
  expect_match(man$note[2], "boom")
  options(op)
})

test_that("derived seeds do not overlap", {
  s <- c(derive_seed(40, 0), derive_seed(40, 999999), derive_seed(41, 0), derive_seed(21, 100000))
  expect_equal(length(unique(s)), 4)
  expect_lt(derive_seed(40, 999999), derive_seed(41, 0))
  expect_error(derive_seed(1, 1e6), "out of range")
})

test_that("expect_rows errors on row-count mismatch", {
  expect_error(expect_rows(data.frame(a = 1:3), 4, "test table"), "test table")
  expect_silent(expect_rows(data.frame(a = 1:3), 3))
})
