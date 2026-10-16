# =============================================================================
# 94_fig4d_rule_weights.R - extraction of rule weights for the classifier
# figure (read-only).
# Inputs: data/model_v2.rds, data/model_v2_rules.tsv
# Output: data/fig4d_rule_weights_v2.tsv (class, rule_rank, gene1,
#   gene2, weight). Weights are the model-object rule scores, asserted pair by
#   pair against the rule table; no statistics computed.
# =============================================================================
source("scripts/00_setup.R")

run_step("94_fig4d_rule_weights", function() {
  clf <- readRDS(proj_path("data", "model_v2.rds"))
  rules <- read_tsv("data/model_v2_rules.tsv")
  check(nrow(rules) == 30L, "Rule count should be 30, actual ", nrow(rules))

  out <- list()
  for (cl in CFG$classes) {
    rr <- rules[rules$class == cl, , drop = FALSE]
    tsp_full <- as.matrix(clf$classifiers[[cl]]$TSPs)
    score_full <- clf$classifiers[[cl]]$score
    check(length(score_full) == nrow(tsp_full), "score and TSPs row counts differ (class ", cl, ")")
    okrow <- !is.na(tsp_full[, 1]) & !is.na(tsp_full[, 2])   # drop NA placeholder rows (M5)
    g1m <- unname(tsp_full[okrow, 1])
    g2m <- unname(tsp_full[okrow, 2])
    scm <- score_full[okrow]
    check(identical(g1m, rr$gene1) && identical(g2m, rr$gene2),
          "rule table and model object mismatch (class ", cl, ")")
    check(all(is.finite(scm)), "non-finite weight value (class ", cl, ")")
    out[[cl]] <- data.frame(class = cl, rule_rank = seq_along(g1m),
                            gene1 = g1m, gene2 = g2m, weight = scm,
                            stringsAsFactors = FALSE)
  }
  w <- do.call(rbind, out)
  rownames(w) <- NULL
  check(nrow(w) == 30L && !anyNA(w$weight), "weight table structure anomaly")
  write_tsv(w, "data/fig4d_rule_weights_v2.tsv")
  log_msg("94 done: 30 rule weights, range ", sprintf("%.4f", min(w$weight)), "–",
          sprintf("%.4f", max(w$weight)))
  invisible(NULL)
})
