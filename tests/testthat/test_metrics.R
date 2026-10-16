# Metric tests (balanced accuracy with NA as error, group bootstrap, H2 criterion message).

test_that("balanced accuracy: NA predictions counted as errors", {
  true <- c("AD", "AD", "BM", "BM", "UC", "UC")
  pred <- c("AD", "BM", "BM", "BM", "UC", NA)
  expect_equal(balanced_accuracy(true, pred), mean(c(0.5, 1, 0.5)))
})

test_that("absent classes are not scored as 0 (artefact of fold-level macro-F1)", {
  # studies with only AD/BM classes: the v1 approach scored UC F1 as 0, capping at 0.667
  expect_equal(balanced_accuracy(c("AD", "AD", "BM", "BM"), c("AD", "AD", "BM", "BM")), 1)
  pcm <- per_class_metrics(c("AD", "BM"), c("AD", "BM"))
  expect_true(is.na(pcm$recall[pcm$class == "UC"]))
})

test_that("majority vote: ties recorded as NA", {
  expect_true(is.na(majority_label(c("AD", "BM"))))
  expect_equal(majority_label(c("AD", "AD", "BM")), "AD")
  expect_true(is.na(majority_label(c(NA, NA))))
})

test_that("group bootstrap: stratified by class composition, all three classes in every resample", {
  df <- data.frame(gsm = paste0("s", 1:12), study = rep(paste0("st", 1:6), each = 2),
                   cluster_id = rep(paste0("c", 1:6), each = 2), donor = paste0("d", 1:12),
                   true = rep(c("AD", "BM", "UC"), each = 4), pred = rep(c("AD", "BM", "UC"), each = 4))
  b <- cluster_bootstrap(df, function(d) length(unique(d$true)), R = 200, seed = 1)
  expect_true(all(b$boot == 3))
  expect_true(b$unreliable)                       # only 2 clusters per class < 5
  expect_equal(unname(b$clusters_per_class), c(2L, 2L, 2L))
})

test_that("H2 criterion message", {
  pc <- data.frame(class = CFG$classes, recall = c(0.9, 0.9, 0.9))
  expect_equal(h2_verdict(list(ba = 0.9, per_class = pc), list(lo = 0.8)), "supported")
  expect_equal(h2_verdict(list(ba = 0.9, per_class = pc), list(lo = 0.6)), "not supported")
  pc2 <- pc; pc2$recall[3] <- NA
  expect_equal(h2_verdict(list(ba = 0.9, per_class = pc2), list(lo = 0.8)), "not testable")
  expect_equal(h2_verdict(list(ba = 0.9, per_class = pc), list(lo = 0.8), valid = FALSE), "not testable")
})

test_that("summarize_predictions: complete fields and correct donor-level aggregation", {
  df <- data.frame(gsm = paste0("s", 1:6), study = "A", cluster_id = "c1",
                   donor = c("d1", "d1", "d1", "d2", "d3", "d4"),
                   true = c("AD", "AD", "AD", "BM", "UC", "UC"),
                   pred = c("AD", "AD", "BM", "BM", "UC", "AD"))
  s <- summarize_predictions(df)
  expect_equal(s$n_samples, 6)
  dl <- donor_level(df)
  expect_equal(nrow(dl), 4)
  expect_equal(dl$pred[dl$donor == "d1"], "AD")
  expect_equal(empirical_p(0.9, c(0.1, 0.2, 0.95)), 2 / 4)
})
