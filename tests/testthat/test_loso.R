# Leave-one-group-out executor tests (one row per fold; common-successful-fold comparisons).

mk_meta <- function() {
  data.frame(gsm = paste0("s", 1:12), study = rep(c("A", "B", "C"), each = 4),
             cluster_id = rep(c("c1", "c2", "c3"), each = 4), donor = paste0("d", 1:12),
             source = c("AD", "AD", "BM", "BM", "AD", "AD", "BM", "BM", "UC", "UC", "UC", "UC"))
}

test_that("a failed fold still yields one row; predictions come from successful folds only", {
  meta <- mk_meta(); folds <- make_folds(meta)
  fp <- function(cl) {
    if (cl == "c2") stop("simulation failed")
    g <- meta$gsm[meta$cluster_id == cl]
    list(pred = data.frame(gsm = g, pred = meta$source[match(g, meta$gsm)]), n_skipped = 0)
  }
  r <- run_loso(folds, "m1", fp, meta)
  expect_equal(nrow(r$folds), 3)
  expect_equal(r$folds$status, c("ok", "failed", "ok"))
  expect_match(r$folds$note[2], "simulation failed")
  expect_equal(nrow(r$pred), 8)
})

test_that("warning when a multi-class study is predicted as a single class", {
  meta <- mk_meta(); folds <- make_folds(meta)
  fp <- function(cl) { g <- meta$gsm[meta$cluster_id == cl]; list(pred = data.frame(gsm = g, pred = "AD")) }
  r <- run_loso(folds, "m2", fp, meta)
  expect_match(r$folds$note[r$folds$fold == "c1"], "only one class")
  expect_equal(r$folds$note[r$folds$fold == "c3"], "")   # no warning for single-class studies
})

# ---- NA rate and NaN warnings ------------------------

test_that("all-NA predictions and NaN scores trigger warnings", {
  meta <- mk_meta(); folds <- make_folds(meta)
  fp_na <- function(cl) {
    g <- meta$gsm[meta$cluster_id == cl]
    list(pred = data.frame(gsm = g, pred = NA_character_, score_AD = NA_real_,
                           score_BM = NA_real_, score_UC = NA_real_), n_skipped = 0)
  }
  r <- run_loso(folds, "m_na", fp_na, meta)
  expect_match(r$folds$note[r$folds$fold == "c1"], "NA prediction fraction")
  fp_nan <- function(cl) {
    g <- meta$gsm[meta$cluster_id == cl]
    list(pred = data.frame(gsm = g, pred = "AD", score_AD = NaN, score_BM = 0, score_UC = 0),
         n_skipped = 0)
  }
  r2 <- run_loso(folds, "m_nan", fp_nan, meta)
  expect_match(r2$folds$note[r2$folds$fold == "c1"], "NaN")
})

test_that("comparisons use common successful folds only", {
  ft1 <- data.frame(fold = c("c1", "c2", "c3"), status = c("ok", "failed", "ok"))
  ft2 <- data.frame(fold = c("c1", "c2", "c3"), status = c("ok", "ok", "failed"))
  expect_equal(common_ok_folds(list(ft1, ft2)), "c1")
})

test_that("attach_truth joins truth by GSM and errors when no truth is found", {
  meta <- mk_meta()
  pr <- data.frame(gsm = c("s1", "s5"), pred = c("AD", "BM"))
  out <- attach_truth(pr, meta)
  expect_equal(out$true, c("AD", "AD"))
  expect_error(attach_truth(data.frame(gsm = "zzz", pred = "AD"), meta))
})
