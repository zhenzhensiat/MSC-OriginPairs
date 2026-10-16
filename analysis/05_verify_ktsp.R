# =============================================================================
# 05_verify_ktsp.R —— end-to-end verification of the multiclassPairs wrapper on
# simulated data (3 classes x 3 studies with planted differences): train (with
# and without anchors) -> predict -> extract rules. Six assertions cover anchor
# membership of every rule, output columns, direction separation, the
# single-rule and zero-rule degenerate cases, and bit-identical agreement of the
# wrapper with a direct package call. str() output saved to logs/05_verify_ktsp_str.txt.
# =============================================================================
source("scripts/00_setup.R")

run_step("05_verify_ktsp", function() {
  need_pkgs(c("multiclassPairs", "switchBox"))
  assert_multiclasspairs_api()
  set.seed(derive_seed(5, 1))
  cls <- CFG$classes
  genes <- paste0("G", sprintf("%04d", 1:3000))
  base <- matrix(rnorm(length(genes) * 90), nrow = length(genes),
                 dimnames = list(genes, paste0("S", 1:90)))   # column names must be supplied first (S1..S90)
  # note: study and source must be interleaved (the platform_wise filter requires every study to contain all three classes)
  meta <- data.frame(gsm = paste0("S", 1:90), source = rep(cls, each = 30),
                     donor = paste0("d", rep(1:9, each = 10)),
                     study = rep(rep(c("ST1", "ST2", "ST3"), each = 10), 3),
                     arm = "include", stringsAsFactors = FALSE, row.names = paste0("S", 1:90))
  # planted marker genes: 15 per class, two-way shift (target class +2.5, other classes -0.5) ——
  # with a one-way +1.5 shift the mathematical upper bound of direction separation for anchor × random pairing is ≈0.5, so assertion 3 has no power (M5 debugging conclusion)
  markers <- lapply(cls, function(cl) sample(genes, 15))
  names(markers) <- cls
  for (cl in cls) {
    base[markers[[cl]], meta$source == cl] <- base[markers[[cl]], meta$source == cl] + 2.5
    base[markers[[cl]], meta$source != cl] <- base[markers[[cl]], meta$source != cl] - 0.5
  }
  tv <- list(expr = base[, meta$gsm, drop = FALSE], meta = meta)

  # ---- without anchors (main model) ----
  clf0 <- ktsp_train(tv, anchors = NULL, seed = derive_seed(5, 2))
  pr0 <- ktsp_predict(clf0, tv$expr)
  check(all(c("max_score", paste0("score_", cls)) %in% names(pr0)), "assertion 2 failed: prediction output columns")
  acc <- mean(pr0$pred == meta$source)
  check(acc > 0.9, "simulation accuracy of the anchor-free model unexpectedly low: ", round(acc, 3))
  rules0 <- ktsp_rules(clf0, tv$expr, meta$source)
  check(all(rules0$separation > 0.5), "assertion 3 failed: rule(s) with direction separation ≤0.5 exist")

  # ---- with anchors ----
  anchors <- unlist(lapply(markers, `[`, 1:5))
  clf1 <- ktsp_train(tv, anchors = anchors, seed = derive_seed(5, 3))
  pr1 <- ktsp_predict(clf1, tv$expr)
  check(all(c("max_score", paste0("score_", cls)) %in% names(pr1)), "assertion 2 failed (anchor version)")
  rules1 <- ktsp_rules(clf1, tv$expr, meta$source)
  in_anchor <- rules1$gene1 %in% anchors | rules1$gene2 %in% anchors
  check(all(in_anchor), "assertion 1 failed: rule(s) with neither gene an anchor (", sum(!in_anchor), " items)")
  check(all(rules1$separation > 0.5), "assertion 3 failed (anchor version)")

  # ---- degenerate-case tolerance (2026-10-11 fixes for assertions 4-6; background: the package drops a dimension and crashes when surviving rules=1, and yields NaN columns when =0) ----
  tsp_of <- function(cl) as.matrix(clf0$classifiers[[cl]]$TSPs)
  n_surv <- function(genes, cl) { t <- tsp_of(cl); sum(t[, 1] %in% genes & t[, 2] %in% genes) }

  # (4) exactly 1 surviving rule: keep only the genes of the first AD rule; BM/UC keep the first two
  keep4 <- unique(c(tsp_of("AD")[1, ], unlist(tsp_of("BM")[1:2, ]), unlist(tsp_of("UC")[1:2, ])))
  n4 <- vapply(cls, function(cl) n_surv(keep4, cl), integer(1)); names(n4) <- cls
  check(n4[["AD"]] == 1 && n4[["BM"]] >= 2 && n4[["UC"]] >= 2,
        "assertion 4 setup failed: surviving rule counts ", paste(n4, collapse = "/"))
  D4 <- tv$expr[keep4, , drop = FALSE]
  pr4 <- ktsp_predict(clf0, D4)
  tspA <- tsp_of("AD"); kA <- tspA[, 1] %in% keep4 & tspA[, 2] %in% keep4
  check(sum(kA) == 1, "assertion 4 setup failed: AD surviving count is not 1")
  check(all(pr4$score_AD == as.numeric(D4[tspA[kA, 1], ] > D4[tspA[kA, 2], ])),
        "assertion 4 failed: single-rule scores do not match the switchBox comparison")
  check(all(!is.na(pr4$pred)), "assertion 4 failed: predictions contain NA")

  # (5) 0 surviving rules: the test matrix contains only the two genes of the first BM rule → AD/UC fully degenerate, BM exactly 1
  # (rule genes are heavily shared across classes —— inherent to one-vs-rest selection —— so the minimal retained set is used instead of gene removal)
  K5 <- unique(c(rules0$gene1[rules0$class == "BM"][1], rules0$gene2[rules0$class == "BM"][1]))
  n5 <- vapply(cls, function(cl) n_surv(K5, cl), integer(1)); names(n5) <- cls
  check(n5[["AD"]] == 0 && n5[["BM"]] == 1, "assertion 5 setup failed: surviving rule counts ", paste(n5, collapse = "/"))
  pr5 <- ktsp_predict(clf0, tv$expr[K5, , drop = FALSE])
  allsc <- as.matrix(pr5[, paste0("score_", cls)])
  check(all(is.na(pr5$score_AD)), "assertion 5 failed: fully degenerate class should be NA")
  check(!any(is.nan(allsc)), "assertion 5 failed: NaN present")
  check(all(!is.na(pr5$pred)) && all(pr5$pred %in% c("BM", "UC")), "assertion 5 failed: predictions should come from the remaining classes")

  # (6) regular case: wrapper output is bit-identical to a direct package call
  pr6 <- ktsp_predict(clf0, tv$expr)
  raw6 <- multiclassPairs::predict_one_vs_rest_TSP(classifier = clf0, Data = tv$expr,
            tolerate_missed_genes = TRUE, weighted_votes = TRUE, classes = cls, verbose = FALSE)
  check(identical(pr6$max_score, as.character(raw6$max_score)), "assertion 6 failed: max_score differs")
  for (cl in cls) check(identical(pr6[[paste0("score_", cl)]], as.numeric(raw6[[cl]])),
                        "assertion 6 failed: ", cl, " score differs from the direct call")

  # ---- str() trace ----
  sink(proj_path("logs", "05_verify_ktsp_str.txt"), split = FALSE)
  cat("=== clf0 (main model) ===\n"); str(clf0, max.level = 2, list.len = 12)
  cat("\n=== clf1 (anchor version) ===\n"); str(clf1, max.level = 2, list.len = 12)
  cat("\n=== pr0 ===\n"); str(pr0)
  cat("\n=== rules1 ===\n"); str(rules1)
  cat("\n=== pr4 (single-rule degenerate) ===\n"); str(pr4)
  cat("\n=== pr5 (fully degenerate class) ===\n"); str(pr5)
  sink()
  log_msg("M5 verification passed: anchor-free accuracy ", round(acc, 3), "; rule counts 0/1 = ",
          nrow(rules0), "/", nrow(rules1), "; rules covered by anchors ", sum(in_anchor), "/", nrow(rules1),
          "; degenerate-case tolerance 4-6 passed (surviving counts ", paste(n4, collapse = "/"), "→1, ",
          paste(n5, collapse = "/"), "→0)")
  invisible(list(clf0 = clf0, clf1 = clf1, rules0 = rules0, rules1 = rules1, pr4 = pr4, pr5 = pr5))
})
