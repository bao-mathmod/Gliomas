#!/usr/bin/env bash
set -Euo pipefail
IFS=$'\n\t'

# ======== PATHS (edit if you change layout) ===========
BASE="/mnt/18T/chibao/gliomas/data/upstream/scRNA"
OUT="$BASE/scRNA_clean_2"
MANIFEST="$OUT/scrna_manifest.tsv"    # <-- your provided path
# ======================================================

mkdir -p "$OUT"/{rds,logs,summary,plots,tmp}

log(){ printf "[%s] %s\n" "$(date '+%Y-%m-%d %H:%M:%S%z')" "$*" | tee -a "$OUT/logs/main.log" ; }
err(){ printf "[%s] ERROR: %s\n" "$(date '+%Y-%m-%d %H:%M:%S%z')" "$*" | tee -a "$OUT/logs/error.log" >&2 ; }

log "scRNA ingest+QC started. MANIFEST=$MANIFEST OUT=$OUT"

Rscript --vanilla - "$MANIFEST" "$OUT" <<'RS' || { echo "R failed"; exit 1; }
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(stringr)
  library(dplyr)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(glmGamPoi)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
manifest_path <- args[1]
outdir <- args[2]
log_main  <- file.path(outdir, "logs", "main.log")
log_error <- file.path(outdir, "logs", "error.log")
summary_tsv <- file.path(outdir, "summary", "cohort_summary.tsv")
dir.create(dirname(summary_tsv), recursive=TRUE, showWarnings=FALSE)

log <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"), paste0(..., collapse=" ")),
      file=log_main, append=TRUE)
}
elog <- function(...) {
  cat(sprintf("[%s] ERROR: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"), paste0(..., collapse=" ")),
      file=log_error, append=TRUE)
}

# ---------- Readers ----------
read_10x_h5 <- function(h5_path){
  x <- Read10X_h5(h5_path)
  if (is.list(x)) {
    if (!is.null(x[["Gene Expression"]])) x <- x[["Gene Expression"]] else x <- x[[1]]
  }
  x
}
read_10x_mtx_dir <- function(dirpath){ Read10X(data.dir = dirpath) }

read_10x_mtx_stem <- function(stem){
  pick <- function(sans_ext_gz, sans_ext_plain){
    gz <- paste0(stem, sans_ext_gz); pl <- paste0(stem, sans_ext_plain)
    if (file.exists(gz)) return(gz)
    if (file.exists(pl)) return(pl)
    return(NA_character_)
  }
  bfile <- pick("barcodes.tsv.gz","barcodes.tsv")
  ffile <- { f1 <- pick("features.tsv.gz","features.tsv"); if (!is.na(f1)) f1 else pick("genes.tsv.gz","genes.tsv") }
  mfile <- {
    m1 <- pick("matrix.mtx.gz","matrix.mtx")
    if (!is.na(m1)) m1 else {
      m2 <- { m2a <- paste0(sub("_$", "", stem), "matrix.mtx.gz"); m2b <- paste0(sub("_$", "", stem), "matrix.mtx"); if (file.exists(m2a)) m2a else if (file.exists(m2b)) m2b else NA_character_ }
      if (!is.na(m2)) m2 else { m3 <- paste0(stem, "matrix.gz"); if (file.exists(m3)) m3 else { m4 <- paste0(sub("_$", "", stem), "matrix.gz"); if (file.exists(m4)) m4 else NA_character_ } }
    }
  }
  missing <- setNames(c(is.na(bfile), is.na(ffile), is.na(mfile)), c("barcodes","features/genes","matrix"))
  if (any(missing)) stop(sprintf("Triplet missing for stem=%s :: need %s", stem, paste(names(missing)[missing], collapse=", ")))
  barcodes <- data.table::fread(bfile, header=FALSE, sep="\t", data.table=FALSE)[,1]
  feats_dt <- data.table::fread(ffile, header=FALSE, sep="\t", data.table=FALSE)
  genes <- if (ncol(feats_dt) >= 2) feats_dt[[2]] else feats_dt[[1]]
  M <- Matrix::readMM(mfile); M <- as(M, "dgCMatrix")
  if (nrow(M) != length(genes)) stop("Row mismatch matrix vs features")
  if (ncol(M) != length(barcodes)) stop("Col mismatch matrix vs barcodes")
  rownames(M) <- make.unique(genes)
  colnames(M) <- barcodes
  M
}

read_tsv_wide_sparse <- function(path) {
  dt <- data.table::fread(path, sep="\t", header=TRUE, data.table=FALSE)
  rn <- dt[[1]]; dt[[1]] <- NULL
  unique_rn <- make.unique(as.character(rn))
  if (any(duplicated(rn))) { log("Found and resolved duplicate gene names in", basename(path)) }
  mat <- as.matrix(dt); rownames(mat) <- unique_rn
  Matrix::Matrix(mat, sparse=TRUE)
}

# ---------- scRNA QC helpers ----------
# <-- LATEST MODIFICATION: Robust QC metric calculation
add_qc_metrics <- function(sobj){
  # Calculate mitochondrial percentage safely
  mt_genes <- grep("^MT-", rownames(sobj), value=TRUE)
  if (length(mt_genes) > 0) {
    sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, features=mt_genes)
  } else {
    log("No mitochondrial (MT-) genes found for sample:", sobj$sample_id[1], ". Setting percent.mt to 0.")
    sobj[["percent.mt"]] <- 0
  }
  
  # Calculate ribosomal percentage safely
  rrna_genes <- grep("^RP[SL]", rownames(sobj), value=TRUE)
  if (length(rrna_genes) > 0) {
    sobj[["percent.ribo"]] <- PercentageFeatureSet(sobj, features=rrna_genes)
  } else {
    log("No ribosomal (RP[SL]) genes found for sample:", sobj$sample_id[1], ". Setting percent.ribo to 0.")
    sobj[["percent.ribo"]] <- 0
  }
  sobj
}

flag_ribo_outliers <- function(sobj, k=3){
  rb <- sobj$percent.ribo; rb_med <- stats::median(rb, na.rm=TRUE); rb_mad <- stats::mad(rb, center=rb_med, constant=1, na.rm=TRUE)
  sobj$ribo_high_flag <- rb > (rb_med + k*rb_mad); sobj$ribo_med <- rb_med; sobj$ribo_mad <- rb_mad
  sobj
}

adaptive_cap <- function(v, hard_cap){
  med <- stats::median(v, na.rm=TRUE); madv <- stats::mad(v, center=med, constant=1, na.rm=TRUE)
  upper <- med + 3*madv
  if (!is.finite(upper)) upper <- hard_cap
  pmin(upper, hard_cap)
}

# <-- LATEST MODIFICATION: Add a 'dry run' to prevent filtering all cells
apply_scrna_filters <- function(sobj,
                                min_features=300,
                                min_counts=500,
                                hard_cap_features=8000,
                                hard_cap_pctmt=15){
  cap_feat <- adaptive_cap(sobj$nFeature_RNA, hard_cap_features)
  cap_mt   <- adaptive_cap(sobj$percent.mt,   hard_cap_pctmt)
  
  cells_to_keep <- (
    !is.na(sobj$nFeature_RNA) & sobj$nFeature_RNA >= min_features & sobj$nFeature_RNA <= cap_feat &
    !is.na(sobj$nCount_RNA)   & sobj$nCount_RNA   >= min_counts &
    !is.na(sobj$percent.mt)   & sobj$percent.mt   <= cap_mt
  )
  
  if (sum(cells_to_keep, na.rm = TRUE) == 0) {
    log("WARNING: All cells for sample", sobj$sample_id[1], "would be removed by QC filters. Skipping filtering for this sample.")
    return(sobj)
  }
  subset(sobj, subset = cells_to_keep)
}

save_qc_violin <- function(sobj, outdir, sample_id, stage=c("pre_qc","post_qc")){
  stage <- match.arg(stage); pdir <- file.path(outdir, "plots", sample_id, stage); dir.create(pdir, recursive=TRUE, showWarnings=FALSE)
  feats <- c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo")
  for (feat in feats) if (feat %in% colnames(sobj@meta.data)) {
    p <- VlnPlot(sobj, features=feat, pt.size=0.1) + NoLegend() + ggtitle(paste0(sample_id," • ",feat," (",stage,")"))
    ggplot2::ggsave(file.path(pdir, paste0(feat,".png")), p, width=7, height=5, dpi=150)
  }
}

run_denoise <- function(sobj){
  if (ncol(sobj) < 20) {
    log("Skipping doublet detection for", sobj$sample_id[1], "- not enough cells:", ncol(sobj))
    sobj$doublet_score <- NA; sobj$doublet_call  <- "singlet"
    return(sobj)
  }
  sce <- as.SingleCellExperiment(sobj); dbl <- scDblFinder::scDblFinder(sce, samples = sobj$sample_id)
  sobj$doublet_score <- dbl$scDblFinder.score; sobj$doublet_call  <- dbl$scDblFinder.class
  subset(sobj, subset = !is.na(doublet_call) & doublet_call == "singlet")
}

sct_normalize <- function(sobj){
  suppressWarnings({
    sobj <- SCTransform(sobj, vst.flavor="v2", verbose=FALSE, vars.to.regress = "percent.mt")
    sobj <- RunPCA(sobj, verbose=FALSE); sobj <- RunUMAP(sobj, dims=1:30, verbose=FALSE)
  })
  sobj
}

# ---------- ingest loop ----------
if (!file.exists(summary_tsv)) {
  write.table(data.frame(project_id=character(), sample_id=character(), cells_raw=integer(), cells_postQC=integer(), mean_genes=double(), pct_mt_mean=double(), pct_ribo_mean=double(), ribo_outlier_rate=double(), pct_doublet_removed=double(), thresholds=character(), path_rds=character(), stringsAsFactors=FALSE),
              file=summary_tsv, sep="\t", quote=FALSE, row.names=FALSE)
}

process_one <- function(project, fmt, src, sample){
  if (fmt != "tsv_wide") {
    rds_path_check <- file.path(outdir, "rds", paste0(sample, ".rds"))
    if (file.exists(rds_path_check)) { log(sprintf("Skipping [%s], output already exists: %s", sample, rds_path_check)); return(invisible(NULL)) }
  }
  log(sprintf("Processing [%s] fmt=%s src=%s sample_id=%s", project, fmt, src, sample))
  sobj <- switch(fmt,
    "10x_h5" = { m <- read_10x_h5(src); CreateSeuratObject(m, assay="RNA", project=project) },
    "mtx" = { if (dir.exists(src)) { m <- read_10x_mtx_dir(src) } else { m <- read_10x_mtx_stem(src) }; CreateSeuratObject(m, assay="RNA", project=project) },
    "tsv_wide" = NULL, stop(sprintf("Unknown format: %s", fmt)))
  
  if (!is.null(sobj)) {
    sobj$sample_id <- sample
    sobj <- add_qc_metrics(sobj); raw_n <- ncol(sobj); sobj <- flag_ribo_outliers(sobj, k=3)
    save_qc_violin(sobj, outdir, sobj$sample_id[1], stage="pre_qc")
    cap_feat <- adaptive_cap(sobj$nFeature_RNA, 8000); cap_mt <- adaptive_cap(sobj$percent.mt, 15)
    sobj <- apply_scrna_filters(sobj)
    before_dbl <- ncol(sobj); sobj <- run_denoise(sobj); after_dbl <- ncol(sobj)
    pct_dbl_removed <- ifelse(before_dbl > 0, (before_dbl - after_dbl) / before_dbl * 100, NA_real_)
    save_qc_violin(sobj, outdir, sobj$sample_id[1], stage="post_qc")
    sobj <- sct_normalize(sobj)
    colnames(sobj) <- paste0(colnames(sobj), "-", sobj$sample_id[1])
    rds_path <- file.path(outdir, "rds", paste0(sobj$sample_id[1], ".rds")); saveRDS(sobj, rds_path, compress="xz")
    thr <- sprintf("min_feat=300, min_counts=500, cap_feat<=%d, cap_mt<=%.2f", as.integer(cap_feat[1]), as.numeric(cap_mt[1]))
    df <- data.frame(project_id = project, sample_id  = sobj$sample_id[1], cells_raw  = raw_n, cells_postQC = ncol(sobj), mean_genes = mean(sobj$nFeature_RNA), pct_mt_mean = mean(sobj$percent.mt), pct_ribo_mean = mean(sobj$percent.ribo), ribo_outlier_rate = mean(sobj$ribo_high_flag), pct_doublet_removed = pct_dbl_removed, thresholds = thr, path_rds = rds_path, stringsAsFactors = FALSE)
    write.table(df, file=summary_tsv, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
    log(sprintf("Saved %s (postQC=%d)", rds_path, ncol(sobj)))
    return(invisible(NULL))
  }

  if (fmt == "tsv_wide") {
    sp <- read_tsv_wide_sparse(src); cells <- colnames(sp); gsm_prefix <- stringr::str_extract(cells, "GSM\\d+"); have_gsm <- !is.na(gsm_prefix)
    sample_keys <- if (any(have_gsm)) gsm_prefix else stringr::str_replace(cells, "_.*$", "")
    for (sk in unique(sample_keys)) {
      sid <- if (!is.na(sample) && nzchar(sample)) sample else paste0(project,"__",sk)
      rds_path_check <- file.path(outdir, "rds", paste0(sid, ".rds"))
      if (file.exists(rds_path_check)) { log(sprintf("Skipping [%s], output already exists: %s", sid, rds_path_check)); next }
      log(sprintf("Processing sub-sample [%s] from tsv_wide file", sid))
      idx <- which(sample_keys == sk); sub <- sp[, idx, drop=FALSE]
      sobj <- CreateSeuratObject(sub, assay="RNA", project=project)
      sobj$sample_id <- sid
      sobj <- add_qc_metrics(sobj); raw_n <- ncol(sobj); sobj <- flag_ribo_outliers(sobj, k=3)
      save_qc_violin(sobj, outdir, sid, "pre_qc")
      cap_feat <- adaptive_cap(sobj$nFeature_RNA, 8000); cap_mt <- adaptive_cap(sobj$percent.mt, 15)
      sobj <- apply_scrna_filters(sobj)
      before_dbl <- ncol(sobj); sobj <- run_denoise(sobj); after_dbl <- ncol(sobj)
      pct_dbl_removed <- ifelse(before_dbl > 0, (before_dbl - after_dbl) / before_dbl * 100, NA_real_)
      save_qc_violin(sobj, outdir, sid, "post_qc")
      sobj <- sct_normalize(sobj)
      colnames(sobj) <- paste0(colnames(sobj), "-", sid)
      rds_path <- file.path(outdir, "rds", paste0(sid, ".rds")); saveRDS(sobj, rds_path, compress="xz")
      thr <- sprintf("min_feat=300, min_counts=500, cap_feat<=%d, cap_mt<=%.2f", as.integer(cap_feat[1]), as.numeric(cap_mt[1]))
      df <- data.frame(project_id = project, sample_id  = sid, cells_raw  = raw_n, cells_postQC = ncol(sobj), mean_genes = mean(sobj$nFeature_RNA), pct_mt_mean = mean(sobj$percent.mt), pct_ribo_mean = mean(sobj$percent.ribo), ribo_outlier_rate = mean(sobj$ribo_high_flag), pct_doublet_removed = pct_dbl_removed, thresholds = thr, path_rds = rds_path, stringsAsFactors = FALSE)
      write.table(df, file=summary_tsv, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
      log(sprintf("Saved %s (postQC=%d)", rds_path, ncol(sobj)))
    }
  }
}

# main
man <- data.table::fread(manifest_path, sep="\t", header=TRUE, data.table=FALSE)
for (i in seq_len(nrow(man))) {
  row <- man[i,]; project <- row$project_id; fmt <- row$format; src <- row$path_or_stem; sample  <- row$sample_id
  tryCatch({
    process_one(project, fmt, src, sample)
  }, error=function(e){
    elog(sprintf("Failed: project=%s format=%s src=%s sample=%s :: %s", project, fmt, src, sample, e$message))
  })
}
log(sprintf("All done. Summary at %s", summary_tsv))
RS

log "scRNA ingest+QC finished."