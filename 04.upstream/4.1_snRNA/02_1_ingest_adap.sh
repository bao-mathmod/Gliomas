#!/usr/bin/env bash
set -Euo pipefail
IFS=$'\n\t'

# ======================= CONFIG ==========================
# Adjusted paths for your snRNA data
BASE="/mnt/18T/chibao/gliomas/data/upstream/snRNA/snRNA_adaptive"
OUT="$BASE/snRNA_processed" # New output directory for processed results
MANIFEST="/mnt/18T/chibao/gliomas/data/upstream/snRNA/snRNA_clean_2/snrna_manifest.tsv" # Assuming your manifest is here
# =========================================================

mkdir -p "$OUT"/{rds,logs,summary,plots}

log(){ printf "[%s] %s\n" "$(date '+%Y-%m-%d %H:%M:%S%z')" "$*" | tee -a "$OUT/logs/main.log" ; }
err(){ printf "[%s] ERROR: %s\n" "$(date '+%Y-%m-%d %H:%M:%S%z')" "$*" | tee -a "$OUT/logs/error.log" >&2 ; }

log "snRNA ingest+QC (ADAPTIVE) started. MANIFEST=$MANIFEST OUT=$OUT"

Rscript --vanilla - "$MANIFEST" "$OUT" <<'RS' || { err "Rscript failed"; exit 1; }
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
manifest_path <- args[1]
outdir <- args[2]
log_main  <- file.path(outdir, "logs", "main.log")
log_error <- file.path(outdir, "logs", "error.log")
summary_tsv <- file.path(outdir, "summary", "cohort_summary.tsv")
dir.create(dirname(summary_tsv), recursive=TRUE, showWarnings=FALSE)

# --- Logging helpers ---
log <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"), paste0(..., collapse=" ")),
      file=log_main, append=TRUE)
}
elog <- function(...) {
  cat(sprintf("[%s] ERROR: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"), paste0(..., collapse=" ")),
      file=log_error, append=TRUE)
}

# --- Data Readers (assuming formats from your previous snRNA script) ---
read_10x_h5 <- function(h5_path) Read10X_h5(h5_path)
read_mtx_dir <- function(dirpath) Read10X(data.dir = dirpath)
read_mtx_stem <- function(stem){
    bfile <- list.files(dirname(stem), pattern=paste0(basename(stem),"_barcodes.tsv"), full.names=TRUE)[1]
    ffile <- list.files(dirname(stem), pattern=paste0(basename(stem),"_features.tsv"), full.names=TRUE)[1]
    mfile <- list.files(dirname(stem), pattern=paste0(basename(stem),"_matrix.mtx"), full.names=TRUE)[1]
    if (is.na(bfile)) bfile <- list.files(dirname(stem), pattern=paste0(basename(stem),"_barcodes.tsv.gz"), full.names=TRUE)[1]
    if (is.na(ffile)) ffile <- list.files(dirname(stem), pattern=paste0(basename(stem),"_features.tsv.gz"), full.names=TRUE)[1]
    if (is.na(mfile)) mfile <- list.files(dirname(stem), pattern=paste0(basename(stem),"_matrix.mtx.gz"), full.names=TRUE)[1]
    if (any(is.na(c(bfile, ffile, mfile)))) stop(sprintf("Could not find all MTX files for stem: %s", stem))
    counts <- ReadMtx(mtx = mfile, features = ffile, cells = bfile)
    return(counts)
}
read_tsv_wide_sparse <- function(path) {
  dt <- data.table::fread(path, sep="\t", header=TRUE, data.table=FALSE)
  rn <- dt[[1]]; dt[[1]] <- NULL
  mat <- as.matrix(dt); rownames(mat) <- rn
  Matrix::Matrix(mat, sparse = TRUE)
}

# --- QC and Filtering functions ---
add_qc_metrics <- function(sobj){
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
  return(sobj)
}

# **NEW**: Adaptive filtering function tailored for snRNA-seq data
apply_adaptive_snrna_filters <- function(sobj, mad_k=3) {
  meta <- sobj@meta.data
  
  # 1. Filter by mitochondrial percentage (STRICT FIXED CUTOFF)
  # For nuclei, this is a purity metric. High mito % means cytoplasmic contamination.
  mito_final_thresh <- 1.0
  
  # 2. Filter by feature count (nFeature_RNA) (adaptive upper bound)
  feat_median <- median(meta$nFeature_RNA, na.rm = TRUE)
  feat_mad <- mad(meta$nFeature_RNA, na.rm = TRUE)
  feat_lower_thresh <- 200 # Keep a safe minimum
  feat_upper_thresh <- feat_median + mad_k * feat_mad
  
  # 3. Filter by total counts (nCount_RNA) (adaptive upper bound)
  count_median <- median(meta$nCount_RNA, na.rm = TRUE)
  count_mad <- mad(meta$nCount_RNA, na.rm = TRUE)
  count_upper_thresh <- count_median + mad_k * count_mad
  
  log(sprintf("QC thresholds for %s: nFeature_RNA > %.0f & < %.0f | percent.mt < %.2f | nCount_RNA < %.0f",
      sobj$sample_id[1], feat_lower_thresh, feat_upper_thresh, mito_final_thresh, count_upper_thresh))
  
  subset(sobj, subset = 
    nFeature_RNA > feat_lower_thresh & nFeature_RNA < feat_upper_thresh &
    nCount_RNA < count_upper_thresh &
    percent.mt < mito_final_thresh
  )
}

save_qc_violin <- function(sobj, outdir, sample_id, stage = c("pre_qc","post_qc")){
  pdir <- file.path(outdir, "plots", sample_id)
  dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
  feats <- c("nFeature_RNA","nCount_RNA","percent.mt")
  for (feat in feats) {
    p <- VlnPlot(sobj, features = feat, pt.size = 0) + NoLegend() +
         ggtitle(paste0(sample_id, " • ", feat, " (", stage, ")"))
    ggsave(file.path(pdir, paste0(stage, "_", feat, ".png")), plot=p, width=6, height=5, dpi=100)
  }
}

run_denoise <- function(sobj){
  sce <- as.SingleCellExperiment(sobj)
  dbl_calls <- scDblFinder::scDblFinder(sce, returnType="table")
  sobj$doublet_class  <- dbl_calls[colnames(sobj), "class"]
  subset(sobj, subset = doublet_class == "singlet")
}

sct_normalize <- function(sobj){
  suppressWarnings({
    sobj <- SCTransform(sobj, vst.flavor="v2", verbose=FALSE, vars.to.regress = "percent.mt")
    sobj <- RunPCA(sobj, verbose=FALSE)
    sobj <- RunUMAP(sobj, dims = 1:30, verbose=FALSE)
  })
  return(sobj)
}

# --- Main processing loop ---
manifest <- data.table::fread(manifest_path, sep="\t", header=TRUE, data.table=FALSE)
if (!file.exists(summary_tsv)) {
  write.table(data.frame(project_id=character(), sample_id=character(),
                         cells_raw=integer(), cells_postQC=integer(),
                         median_genes=double(), median_counts=double(), median_pct_mt=double(),
                         path_rds=character(), stringsAsFactors = FALSE),
              file=summary_tsv, sep="\t", quote=FALSE, row.names=FALSE)
}

for (i in seq_len(nrow(manifest))) {
  row <- manifest[i,]
  project <- row$project_id
  fmt     <- row$format
  src     <- row$path_or_stem
  sample  <- row$sample_id
  rds_path <- file.path(outdir, "rds", paste0(sample, ".rds"))
  
  tryCatch({
    log(sprintf("Processing [%s] format=%s sample_id=%s", project, fmt, sample))

    counts_matrix <- switch(fmt,
      "10x_h5"   = read_10x_h5(src),
      "mtx_dir"  = read_mtx_dir(src),
      "mtx_stem" = read_mtx_stem(src),
      "mtx"      = read_mtx_stem(src), # For backward compatibility with old manifest
      "tsv_wide" = read_tsv_wide_sparse(src),
      stop(sprintf("Unknown format: %s", fmt))
    )
    if (is.list(counts_matrix) && !is.data.frame(counts_matrix)) {
        counts_matrix <- counts_matrix[[1]]
    }
    sobj <- CreateSeuratObject(counts_matrix, project = project)
    sobj$sample_id <- sample
    
    sobj <- add_qc_metrics(sobj)
    raw_n <- ncol(sobj)
    save_qc_violin(sobj, outdir, sample, stage = "pre_qc")
    
    sobj <- apply_adaptive_snrna_filters(sobj)
    
    if (ncol(sobj) > 50) {
        sobj <- run_denoise(sobj)
    } else {
        log(sprintf("Skipping doublet removal for %s, too few nuclei post-QC (%d)", sample, ncol(sobj)))
    }
    
    postqc_n <- ncol(sobj)
    if (postqc_n > 50) {
      save_qc_violin(sobj, outdir, sample, stage = "post_qc")
      sobj <- sct_normalize(sobj)
    } else {
      log(sprintf("Skipping normalization for %s, too few nuclei remaining (%d)", sample, postqc_n))
    }
    
    log(sprintf("Saving RDS for %s. Raw nuclei: %d, Post-QC nuclei: %d", sample, raw_n, postqc_n))
    saveRDS(sobj, rds_path, compress="xz")
    
    summary_df <- data.frame(
      project_id = project, sample_id  = sample,
      cells_raw  = raw_n, cells_postQC = postqc_n,
      median_genes = median(sobj$nFeature_RNA, na.rm=TRUE),
      median_counts = median(sobj$nCount_RNA, na.rm=TRUE),
      median_pct_mt = median(sobj$percent.mt, na.rm=TRUE),
      path_rds = rds_path,
      stringsAsFactors = FALSE
    )
    write.table(summary_df, file=summary_tsv, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)

  }, error=function(e){
    msg <- sprintf("Failed on sample %s (src: %s): %s", sample, src, e$message)
    elog(msg)
  })
}

log(sprintf("All processing complete. Summary table at %s", summary_tsv))
RS

log "snRNA ingest+QC (ADAPTIVE) finished."