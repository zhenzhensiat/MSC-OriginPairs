# =============================================================================
# 11d_rebuild_check.R —— consistency check: raw-file rebuild vs ported matrices.
# Compares probe-level rebuilds against the ported symbol matrices (symbol/alias
# mapping, deduplication, per-gene agreement) before the rebuilt versions replace
# the ported ones.
# =============================================================================
source("scripts/00_setup.R")
need_pkgs(c("org.Hs.eg.db", "AnnotationDbi", "edgeR"))
library(org.Hs.eg.db)
run_step("11d_rebuild_check", function() {
  sym_ok <- unique(na.omit(AnnotationDbi::select(org.Hs.eg.db, keys = keys(org.Hs.eg.db), columns = "SYMBOL")$SYMBOL))
  al <- AnnotationDbi::select(org.Hs.eg.db, keys = keys(org.Hs.eg.db), columns = c("SYMBOL","ALIAS"))
  al <- al[!is.na(al$ALIAS) & !is.na(al$SYMBOL), ]
  an <- table(al$ALIAS); ua <- names(an[an == 1])
  amap <- setNames(al$SYMBOL[match(ua, al$ALIAS)], ua)
  mg <- function(g) ifelse(g %in% sym_ok, g, ifelse(g %in% names(amap), amap[g], NA))
  cdup <- function(x) { if (!anyDuplicated(rownames(x))) return(x)
    o <- order(rowMeans(x), apply(x,1,stats::IQR), decreasing=TRUE)
    x[o[!duplicated(rownames(x)[o])],,drop=FALSE] }
  out <- list()
  g <- "GSE199826"
  cf <- list.files(proj_path("data","raw",g), "raw-counts", full.names=TRUE)[1]
  raw <- read.delim(gzfile(cf), check.names=FALSE, stringsAsFactors=FALSE)
  mat <- as.matrix(raw[,-1,drop=FALSE]); storage.mode(mat) <- "numeric"
  rn <- mg(as.character(raw[[1]]))
  mat <- mat[!is.na(rn),,drop=FALSE]; rownames(mat) <- rn[!is.na(rn)]
  mat <- cdup(mat)
  if (max(mat) >= 30) mat <- edgeR::cpm(mat, log=TRUE, prior.count=0.5)
  old <- readRDS(proj_path("data","studies", paste0(g,".rds")))
  cg <- intersect(rownames(mat), rownames(old$expr))
  k <- min(ncol(mat), ncol(old$expr))
  rhos <- vapply(seq_len(k), function(i) suppressWarnings(stats::cor(mat[cg,i], old$expr[cg,i], method="spearman")), numeric(1))
  out[[g]] <- data.frame(gse=g, type="RNA-seq counts", n_raw=nrow(raw), n_mapped=nrow(mat), n_overlap=length(cg), median_spearman=round(stats::median(rhos),6), note="colnames may differ (author ID vs GSM)", stringsAsFactors=FALSE, row.names=NULL)
  log_msg(g, ": RNA-seq | ", nrow(raw), "->", nrow(mat), " genes | Spearman=", round(stats::median(rhos),6))
  for (g in c("GSE12274","GSE57151")) {
    pl <- ifelse(g=="GSE12274","GPL570","GPL16686")
    out[[g]] <- data.frame(gse=g, type=paste(pl,"array"), n_raw=NA, n_mapped=NA, n_overlap=NA, median_spearman=NA, note=paste0("probe IDs, need ",pl,".annot (pending NCBI)"), stringsAsFactors=FALSE, row.names=NULL)
    log_msg(g, ": ", pl, " | need annot file")
  }
  res <- do.call(rbind, out)
  write_tsv(res, "data/m1_rebuild_check.tsv")
  invisible(res)
})
