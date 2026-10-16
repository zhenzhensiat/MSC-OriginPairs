# =============================================================================
# Run from the repository root:  Rscript tests/run_tests.R
# All tests use simulated data; no real data and no multiclassPairs/limma are
# needed. Everything must pass before any real-data script is run.
# =============================================================================
options(stringsAsFactors = FALSE)
source(file.path("R", "utils_engineering.R"), encoding = "UTF-8")
root <- find_root(getwd())
options(msc.root = root)
source(file.path(root, "config", "config.R"), encoding = "UTF-8")
for (f in c("utils_metrics.R", "utils_stats.R", "utils_loso.R", "utils_data.R", "utils_models.R")) {
  source(file.path(root, "R", f), encoding = "UTF-8")
}
testthat::test_dir(file.path(root, "tests", "testthat"), reporter = "summary", stop_on_failure = TRUE)
