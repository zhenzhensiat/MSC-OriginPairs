# Pure model-function tests (no multiclassPairs needed).

test_that("row-wise ANOVA F matches aov", {
  set.seed(51)
  x <- matrix(stats::rnorm(3 * 12), 3)
  g <- rep(c("AD", "BM", "UC"), each = 4)
  f <- row_anova_F(x, g)
  f_ref <- summary(stats::aov(x[1, ] ~ factor(g)))[[1]][["F value"]][1]
  expect_equal(f[1], f_ref, tolerance = 1e-8)
  expect_true(all(is.na(row_anova_F(x, rep("AD", 12)))))
})

mk_tv <- function() {
  set.seed(52)
  genes <- paste0("g", 1:200)
  mk <- function(st) {
    src <- rep(CFG$classes, each = 3)
    x <- matrix(stats::rnorm(200 * 9), 200, dimnames = list(genes, paste0(st, "_", 1:9)))
    x[1:5, src == "AD"] <- x[1:5, src == "AD"] + 5
    x[6:10, src == "BM"] <- x[6:10, src == "BM"] + 5
    x[11:15, src == "UC"] <- x[11:15, src == "UC"] + 5
    list(x = x, m = data.frame(gsm = colnames(x), study = st, cluster_id = paste0("cl_", st), source = src))
  }
  parts <- lapply(c("A", "B", "C"), mk)
  list(expr = do.call(cbind, lapply(parts, `[[`, "x")), meta = do.call(rbind, lapply(parts, `[[`, "m")))
}

test_that("genome anchor: averaging within-study F percentiles recovers the planted genes", {
  tv <- mk_tv()
  top <- genome_anchor_rank(tv, n = 15)
  expect_setequal(top, paste0("g", 1:15))
})

test_that("rank transform and centroid baseline", {
  tv <- mk_tv()
  r <- rank_transform(tv$expr)
  expect_true(all(r > 0 & r < 1))
  expect_equal(order(r[, 1]), order(tv$expr[, 1]))
  pr <- baseline_centroid(tv, tv$expr[, 1:9])
  expect_equal(pr$pred, tv$meta$source[1:9])
})

test_that("mixing: endpoints equal original samples; monotone in fraction", {
  a <- c(5, 8); b <- c(9, 2)
  expect_equal(mix_log2(a, b, 1), a)
  expect_equal(mix_log2(a, b, 0), b)
  m <- vapply(c(0, 0.25, 0.5, 0.75, 1), function(f) mix_log2(a, b, f)[1], numeric(1))
  expect_true(all(diff(m) < 0))   # a[1] < b[1]; the larger f, the closer to a
})

test_that("open set: margins, threshold, rejection", {
  pd <- data.frame(gsm = paste0("s", 1:4), pred = c("AD", "BM", "UC", "AD"),
                   score_AD = c(10, 1, 1, 4), score_BM = c(0, 9, 1, 4), score_UC = c(0, 0, 8, 2),
                   true = c("AD", "BM", "UC", "BM"))
  m <- score_margin(pd)
  expect_equal(unname(m[1]), 1)
  expect_equal(unname(m[4]), 0)
  tau <- choose_tau(pd, coverage = 1)            # smallest margin among correct samples
  ap <- apply_reject(pd, tau)
  expect_true(ap$rejected[4])
  expect_false(any(ap$rejected[1:3]))
  expect_true(is.na(ap$pred_openset[4]))
})
