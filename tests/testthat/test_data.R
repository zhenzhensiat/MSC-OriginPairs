# Study-object convention tests (structure, dimensions, metadata).

mk_study <- function() {
  x <- matrix(c(5, 6, 7, 8, 3, 4), 3, 2, dimnames = list(c("GAPDH", "HOXA10", "CD200"), c("GSM1", "GSM2")))
  list(gse = "GSE1", cluster_id = "c1", layer = "D", stability = FALSE, scale = "logCPM", expr = x,
       meta = data.frame(gsm = c("GSM1", "GSM2"), source = c("AD", "BM"), donor = c("d1", "d2"),
                         passage = c(3, 3), medium = c("FBS", "FBS"), arm = "include"))
}

test_that("valid study object passes validation", {
  expect_silent(validate_study(mk_study()))
})

test_that("invalid study objects are rejected", {
  s <- mk_study(); s$expr[1, 1] <- NA
  expect_error(validate_study(s), "NA")
  s <- mk_study(); rownames(s$expr)[2] <- "GAPDH"
  expect_error(validate_study(s), "duplicate")
  s <- mk_study(); s$meta$gsm <- rev(s$meta$gsm)
  expect_error(validate_study(s), "column names do not match")
  s <- mk_study(); s$expr <- s$expr * 100
  expect_error(validate_study(s), "not log")
  s <- mk_study(); s$layer <- "S"
  expect_error(validate_study(s), "layer")
  s <- mk_study(); s$stability <- NULL
  expect_error(validate_study(s), "stability")
})

# ---- test-aware universe ----------------------------

mk_big_study <- function(gse, cluster, genes, gsm1, gsm2) {
  x <- matrix(stats::rnorm(length(genes) * 2), length(genes), 2,
              dimnames = list(genes, c(gsm1, gsm2)))
  list(gse = gse, cluster_id = cluster, layer = "D", stability = FALSE, scale = "logCPM",
       full_genome = TRUE, expr = x,
       meta = data.frame(gsm = c(gsm1, gsm2), study = gse, cluster_id = cluster,
                         source = c("AD", "BM"), donor = c("d1", "d2"),
                         passage = c(3, 3), medium = c("FBS", "FBS"), arm = "include"))
}

test_that("test_genes returns the gene intersection of a multi-study cluster", {
  genes <- paste0("G", 1:2500)
  s1 <- mk_big_study("GSE1", "c1", genes, "GSM1", "GSM2")
  s2 <- mk_big_study("GSE2", "c1", genes, "GSM3", "GSM4")
  s2$expr <- s2$expr[-1, , drop = FALSE]         # one gene fewer → intersection 2499
  expect_equal(length(test_genes(list(s1, s2), "c1")), 2499)
  expect_error(test_genes(list(s1, s2), "cX"), "no studies")
})

test_that("training_view restrict_genes takes the intersection and errors when fewer than 2000", {
  genes <- paste0("G", 1:2500)
  st <- lapply(1:4, function(i) mk_big_study(paste0("GSE", i), paste0("c", i), genes,
                                             paste0("GSM", i * 2 - 1), paste0("GSM", i * 2)))
  tv <- training_view(st, "c1", restrict_genes = genes[1:2400])
  expect_equal(length(tv$genes), 2400)
  expect_true(all(tv$genes %in% genes[1:2400]))
  expect_error(training_view(st, "c1", restrict_genes = genes[1:100]), "training gene universe has only 100")
})
