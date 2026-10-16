# =============================================================================
# 54_openset_ooc.R - tissues outside the three modelled sources (open set):
# rejection rate and nearest-class assignment (exploratory).
# Questions: (1) are non-BM/AD/UC cells recognisable as outside the three
#   classes? (2) where developmentally related tissues (placenta/amnion ->
#   UC? dental pulp -> ?) are assigned.
# Input: curation/ooc_manifest_v2.tsv (prepared before running): gse, cluster_id,
#   gsm, tissue, donor, passage, medium; tissue vocabulary: placenta / amnion /
#   chorion / dental_pulp / periodontal / synovium / dermal_fibroblast /
#   endometrium / cord_blood / ipsc_derived / other
# Outputs: data/openset_pred_v2.tsv, data/openset_summary_v2.tsv,
#   data/openset_auc_v2.tsv
# =============================================================================
source("scripts/00_setup.R")
source(proj_path("R", "utils_build.R"), encoding = "UTF-8")
need_pkgs(c("multiclassPairs", "switchBox"))

auc_rank <- function(pos, neg) {
  pos <- pos[!is.na(pos)]; neg <- neg[!is.na(neg)]
  if (!length(pos) || !length(neg)) return(NA_real_)
  r <- rank(c(pos, neg))
  (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) / (length(pos) * length(neg))
}

run_step("54_openset_ooc", function() {
  man <- read_tsv("curation/ooc_manifest_v2.tsv")
  man$source <- "other"; man$arm <- "include"
  clf <- readRDS(proj_path("data", "model_v2.rds"))
  tau <- readRDS(proj_path("data", "model_v2_tau.rds"))

  pr <- do.call(rbind, lapply(split(man, man$gse), function(mr) {
    s <- tryCatch(build_from_manifest(mr$gse[1], mr, layer = "O", cluster_id = mr$cluster_id[1]),
                  error = function(e) conditionMessage(e))
    if (is.character(s)) {
      return(data.frame(gsm = mr$gsm, pred = NA_character_, margin = NA_real_, rejected = NA,
                        build_status = s, stringsAsFactors = FALSE))
    }
    p <- apply_reject(ktsp_predict(clf, s$expr), tau)
    data.frame(gsm = p$gsm, pred = p$pred, margin = p$margin, rejected = p$rejected,
               build_status = "ok", stringsAsFactors = FALSE)
  }))
  pr <- safe_merge(pr, man[, c("gsm", "gse", "tissue", "donor")], by = "gsm", type = "1:1")
  write_tsv(pr, "data/openset_pred_v2.tsv")

  sm <- do.call(rbind, lapply(split(pr, pr$tissue), function(d) {
    ok <- d[d$build_status == "ok", , drop = FALSE]
    data.frame(tissue = d$tissue[1], n = nrow(d), n_studies = length(unique(d$gse)),
               rejection_rate = mean(ok$rejected %in% c(TRUE, "TRUE")),
               called_AD = mean(ok$pred == "AD" & !(ok$rejected %in% c(TRUE, "TRUE")), na.rm = TRUE),
               called_BM = mean(ok$pred == "BM" & !(ok$rejected %in% c(TRUE, "TRUE")), na.rm = TRUE),
               called_UC = mean(ok$pred == "UC" & !(ok$rejected %in% c(TRUE, "TRUE")), na.rm = TRUE),
               median_margin = stats::median(ok$margin, na.rm = TRUE),
               post_hoc = TRUE, stringsAsFactors = FALSE)
  }))
  write_tsv(sm, "data/openset_summary_v2.tsv")

  # ability of the margin to distinguish 'in-class (C layer)' from 'out-of-class': AUC (in-class margins should be higher)
  cp <- read_tsv("data/confirm_predictions_v2.tsv")
  auc <- data.frame(comparison = c("C_inclass_vs_all_ooc", paste0("C_inclass_vs_", sort(unique(pr$tissue)))),
                    auc = c(auc_rank(cp$margin, pr$margin),
                            vapply(sort(unique(pr$tissue)), function(t) auc_rank(cp$margin, pr$margin[pr$tissue == t]), numeric(1))),
                    post_hoc = TRUE, stringsAsFactors = FALSE)
  write_tsv(auc, "data/openset_auc_v2.tsv")
  invisible(sm)
})
