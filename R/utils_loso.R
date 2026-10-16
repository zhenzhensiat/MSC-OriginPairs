# =============================================================================
# utils_loso.R —— generic leave-one-group-out executor
#   * each fold emits exactly one row (success or failure); row count equals fold count
#   * method comparisons use only folds where every method succeeded
#   * automatic warnings: single predicted class in a fold; high rule-skip rate
# =============================================================================

# study object conventions: see README; meta summary must contain at least gsm, study, cluster_id, donor, source
make_folds <- function(meta) {
  cl <- sort(unique(meta$cluster_id))
  stats::setNames(as.list(cl), cl)
}

# fit_predict(test_cluster) must return list(pred = data.frame(gsm, pred, ...), n_skipped = number or NA)
run_loso <- function(folds, method_id, fit_predict, meta) {
  rows <- lapply(names(folds), function(f) {
    n_test <- sum(meta$cluster_id == folds[[f]])
    res <- tryCatch(list(ok = TRUE, value = fit_predict(folds[[f]])),
                    error = function(e) list(ok = FALSE, msg = conditionMessage(e)))
    if (isTRUE(res$ok)) {
      pr <- res$value$pred
      check(all(c("gsm", "pred") %in% names(pr)), "pred returned by fit_predict is missing gsm/pred columns")
      alarm <- character(0)
      na_rate <- mean(is.na(pr$pred))
      if (na_rate > 0) alarm <- c(alarm, sprintf("NA prediction fraction %.3f", na_rate))
      sc <- grep("^score_", names(pr), value = TRUE)
      if (length(sc) && any(vapply(pr[sc], function(v) any(is.nan(v)), logical(1)))) {
        alarm <- c(alarm, "score columns contain NaN")
      }
      np <- pr$pred[!is.na(pr$pred)]
      truth <- meta$source[match(pr$gsm, meta$gsm)]
      if (length(np) >= 4 && length(unique(np)) == 1 && length(unique(truth)) > 1) {
        alarm <- c(alarm, "predictions for a multi-class study contain only one class")
      }
      ns <- res$value$n_skipped
      if (is.null(ns)) ns <- NA
      list(fold = data.frame(method = method_id, fold = f, status = "ok", n_test = n_test,
                             n_pred = sum(!is.na(pr$pred)), n_skipped = ns,
                             note = paste(alarm, collapse = "; "), stringsAsFactors = FALSE),
           pred = data.frame(method = method_id, fold = f, pr, stringsAsFactors = FALSE))
    } else {
      list(fold = data.frame(method = method_id, fold = f, status = "failed", n_test = n_test,
                             n_pred = 0L, n_skipped = NA, note = res$msg, stringsAsFactors = FALSE),
           pred = NULL)
    }
  })
  ft <- do.call(rbind, lapply(rows, `[[`, "fold"))
  expect_rows(ft, length(folds), paste0(method_id, " fold table"))
  if (any(ft$status == "failed")) log_msg(method_id, ": failed folds ", paste(ft$fold[ft$status == "failed"], collapse = ","))
  if (any(nzchar(ft$note) & ft$status == "ok")) log_msg(method_id, " alarms: ", paste(unique(ft$note[nzchar(ft$note)]), collapse = " | "))
  list(folds = ft, pred = do.call(rbind, lapply(rows, `[[`, "pred")))
}

# folds where all methods succeeded
common_ok_folds <- function(fold_tables) {
  ok <- lapply(fold_tables, function(ft) ft$fold[ft$status == "ok"])
  Reduce(intersect, ok)
}

# join predictions with truth, cluster, and donor (1:1 merge on GSM)
attach_truth <- function(pred, meta) {
  m <- meta[, c("gsm", "study", "cluster_id", "donor", "source")]
  names(m)[names(m) == "source"] <- "true"
  out <- safe_merge(pred, m, by = "gsm", type = "m:1")
  check(!anyNA(out$true), "some predicted samples have no ground truth")
  out
}
