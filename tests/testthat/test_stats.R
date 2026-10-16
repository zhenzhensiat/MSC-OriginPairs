# Statistics tests (partial conjunction calibration, direction mapping, matched nulls).

test_that("partial conjunction: r=1 uniform under H0; r>1 conservative", {
  set.seed(21)
  p1 <- replicate(4000, pc_fisher(stats::runif(6), 1))
  expect_gt(stats::ks.test(p1, "punif")$p.value, 0.001)
  p3 <- replicate(4000, pc_fisher(stats::runif(6), 3))
  expect_lte(mean(p3 < 0.05), 0.06)
})

test_that("partial conjunction: significant only when most studies show an effect", {
  tiny <- 1e-8
  expect_lt(pc_fisher(c(rep(tiny, 4), 0.5, 0.7), 3), 1e-4)       # 4/6 with an effect
  expect_gt(pc_fisher(c(tiny, tiny, 0.3, 0.5, 0.6, 0.8), 3), 0.05) # only 2/6
  expect_true(is.na(pc_fisher(c(0.1, 0.2), 3)))
})

test_that("direction: +1 significant and -1 not when t is positive", {
  expect_lt(pc_directional(rep(6, 4), rep(10, 4), 1L), 1e-4)
  expect_gt(pc_directional(rep(6, 4), rep(10, 4), -1L), 0.9)
})

test_that("contrast mapping and claim expansion", {
  expect_equal(canonical_contrast("BM", "AD"), list(contrast = "BM-AD", direction = 1L))
  expect_equal(canonical_contrast("UC", "BM"), list(contrast = "BM-UC", direction = -1L))
  expect_error(canonical_contrast("UC", "UC"))
  mk <- data.frame(claim_id = "C1", gene = "CD200", high_source = "UC", low_source = "BM;AD",
                   discovery_gse = "GSE108511")
  ex <- expand_claims(mk)
  expect_equal(nrow(ex), 2)
  expect_setequal(ex$contrast, c("BM-UC", "AD-UC"))
  expect_true(all(ex$direction == -1L))
  bad <- mk; bad$low_source <- "BM(not detected)"
  expect_error(expand_claims(bad), "standardize")
})

mk_de <- function() {
  set.seed(31)
  genes <- paste0("g", 1:50)
  do.call(rbind, lapply(paste0("S", 1:4), function(s) {
    data.frame(study = s, contrast = "BM-AD", gene = genes, logFC = stats::rnorm(50),
               t = c(rep(5, 5), stats::rnorm(45)), df = 10, p = NA, q_within = NA,
               expr_pct = stats::runif(50), detect = stats::runif(50))
  }))
}

test_that("precompute, excluded studies and fast path agree", {
  de <- mk_de(); pre <- precompute_pc(de)
  expect_equal(nrow(pre), 50)
  expect_true(all(pre$K == 4))
  tg <- data.frame(gene = c("g1", "g2", "g10"), contrast = "BM-AD",
                   direction = c(1L, 0L, -1L), exclude = c("", "S1", ""))
  r <- pc_rows(tg, de, pre)
  expect_equal(r$K, c(4, 3, 4))
  expect_lt(r$pc_p[1], 1e-3)
  f <- pc_either_fast(c("g1", "g20"), "BM-AD", pre)
  s <- pc_rows(data.frame(gene = c("g1", "g20"), contrast = "BM-AD", direction = 0L, exclude = ""), de, pre)
  expect_equal(f$pc_p, s$pc_p)
  expect_equal(f$K, s$K)
})

test_that("within-group BH and transfer rate", {
  tb <- data.frame(pc_p = c(1e-6, 0.5, NA, 1e-4))
  fc <- family_calls(tb, q = 0.05)
  expect_equal(fc$testable, c(TRUE, TRUE, FALSE, TRUE))
  expect_equal(fc$transfer, c(TRUE, FALSE, FALSE, TRUE))
  expect_equal(transfer_rate(fc), 2 / 3)
})

test_that("matched nulls: never draw target genes, reproducible, distinct sets, K preserved", {
  set.seed(41)
  feat <- data.frame(gene = paste0("g", 1:600), contrast = "BM-AD", K = sample(2:4, 600, TRUE),
                     expr_pct = stats::runif(600), detect = stats::runif(600))
  targets <- data.frame(gene = paste0("g", 1:10), contrast = "BM-AD")
  nl <- sample_matched_null(targets, feat, n_sets = 5, seed_base = 100)
  expect_equal(nrow(nl), 50)
  expect_false(any(nl$gene %in% targets$gene))
  expect_identical(nl, sample_matched_null(targets, feat, n_sets = 5, seed_base = 100))
  expect_false(identical(nl$gene[nl$set == 1], nl$gene[nl$set == 2]))
  k_t <- feat$K[match(nl$target_gene, feat$gene)]
  k_n <- feat$K[match(nl$gene, feat$gene)]
  lv <- nl$relax_level <= 1
  expect_true(all(k_t[lv] == k_n[lv]))
  expect_error(sample_matched_null(data.frame(gene = "nope", contrast = "BM-AD"), feat, 2, 1), "missing from the feature table")
})

test_that("zero-contrast transfer rate: fast path returns one value within 0–1 per set", {
  de <- mk_de(); pre <- precompute_pc(de)
  feat <- pre[, c("gene", "contrast", "K", "expr_pct", "detect")]
  tg <- data.frame(gene = paste0("g", 1:5), contrast = "BM-AD", direction = 0L, exclude = "")
  nl <- sample_matched_null(tg, feat, n_sets = 4, seed_base = 7)
  v <- null_rates(nl, tg, de, pre)
  expect_length(v, 4)
  expect_true(all(v >= 0 & v <= 1))
})

# multi high_source expansion and no_class_projection skipping
test_that("expand_claims splits BM;AD vs UC into two rows with correct direction", {
  mk <- data.frame(claim_id = "T1", gene = "G1", high_source = "BM;AD", low_source = "UC",
                   discovery_gse = NA_character_, v2_flag = "", stringsAsFactors = FALSE)
  out <- expand_claims(mk)
  expect_equal(nrow(out), 2)
  expect_setequal(out$contrast, c("BM-UC", "AD-UC"))
  expect_true(all(out$direction == 1L))
})
test_that("expand_claims skips no_class_projection rows", {
  mk <- data.frame(claim_id = c("T1", "T2"), gene = c("G1", "G2"),
                   high_source = c("BM", "UC"), low_source = c("AD", "BM"),
                   discovery_gse = NA_character_,
                   v2_flag = c("no_class_projection", ""), stringsAsFactors = FALSE)
  out <- expand_claims(mk)
  expect_equal(nrow(out), 1)
  expect_equal(out$claim_id, "T2")
})

# HGNC gene_group multi-value fields are quoted; parsing must leave no residual quotes and merge same-name families
test_that("hgnc_families: gene_group quote stripping and same-name merging without splitting", {
  tmp <- tempfile(fileext = ".tsv")
  writeLines(c(
    "symbol\tstatus\tgene_group",
    paste0("GA", 1:4, "\tApproved\t\"FamA|FamB\""),
    paste0("GB", 1:4, "\tApproved\t\"FamA\""),
    paste0("GC", 1:4, "\tApproved\t\"FamB\""),
    paste0("GD", 1:4, "\tApproved\t\"FamC|FamB\"")
  ), tmp)
  fam <- hgnc_families(tmp, size = c(2L, 100L))
  expect_false(any(grepl('"', fam$family, fixed = TRUE)))
  expect_setequal(unique(fam$family), c("FamA", "FamB", "FamC"))
  expect_setequal(fam$gene[fam$family == "FamA"], c(paste0("GA", 1:4), paste0("GB", 1:4)))
  expect_setequal(fam$gene[fam$family == "FamB"],
                  c(paste0("GA", 1:4), paste0("GC", 1:4), paste0("GD", 1:4)))
  expect_setequal(fam$gene[fam$family == "FamC"], paste0("GD", 1:4))
})
