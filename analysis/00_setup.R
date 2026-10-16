# =============================================================================
# 00_setup.R —— first line of every script: source("analysis/00_setup.R")
# Convention: always run from the repository root, e.g.
#   Rscript analysis/21_transfer_audit.R
# The repository root is marked by an empty file .msc_originpairs_root
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1, scipen = 6)

source(file.path("R", "utils_engineering.R"), encoding = "UTF-8")
local({
  root <- find_root(getwd())
  options(msc.root = root)
})
source(proj_path("config", "config.R"), encoding = "UTF-8")

# released data use the neutral group codes (G01-G92); map the config value
CFG$loso_exclude_clusters <- sub("^cl_kuala_lumpur_18$", "G05", CFG$loso_exclude_clusters)
for (f in c("utils_metrics.R", "utils_stats.R", "utils_loso.R", "utils_data.R", "utils_models.R")) {
  source(proj_path("R", f), encoding = "UTF-8")
}

for (d in c("data", "logs", "figures", "curation", "docs")) {
  dir.create(proj_path(d), showWarnings = FALSE, recursive = TRUE)
}

# one log per script; timestamps appear only in log files
local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  nm <- if (length(f)) sub("\\.R$", "", basename(f)) else "interactive"
  options(msc.log_file = proj_path("logs", paste0(nm, "_", format(Sys.Date(), "%Y%m%d"), ".log")))
})

need_pkgs <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  check(length(miss) == 0, "missing R packages: ", paste(miss, collapse = ", "), " (install with renv)")
  for (p in pkgs) log_msg("package ", p, " ", as.character(utils::packageVersion(p)))
  invisible(TRUE)
}

# session info, one file per run
writeLines(utils::capture.output(utils::sessionInfo()),
           proj_path("logs", paste0("sessionInfo_", format(Sys.Date(), "%Y%m%d"), ".txt")))
