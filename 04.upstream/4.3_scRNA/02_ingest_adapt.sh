#!/usr/bin/env bash
set -Euo pipefail
IFS=$'\n\t'

# ======================= CONFIG ==========================
BASE="/mnt/18T/chibao/gliomas/data/upstream/scRNA"
OUT="$BASE/scRNA_processed"
MANIFEST="$BASE/scRNA_clean/scrna_manifest.tsv"
# =========================================================

mkdir -p "$OUT"/{rds,logs,summary,plots}

# --- NEW: Define a log for failed samples ---
FAIL_LOG="$OUT/logs/failed_samples.log"
# Clear previous failure logs for a clean run
>"$OUT/logs/error.log"
>"$FAIL_LOG"

log(){ printf "[%s] %s\n" "$(date '+%Y-%m-%d %H:%M:%S%z')" "$*" | tee -a "$OUT/logs/main.log" ; }
err(){ printf "[%s] ERROR: %s\n" "$(date '+%Y-%m-%d %H:%M:%S%z')" "$*" | tee -a "$OUT/logs/error.log" >&2 ; }

log "scRNA ingest+QC (FINAL, RESUMABLE) started. MANIFEST=$MANIFEST OUT=$OUT"

Rscript --vanilla - "$MANIFEST" "$OUT" "$FAIL_LOG" <<'RS' || { err "Rscript failed"; exit 1; }
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
fail_log_path <- args[3] # New argument for the failure log
log_main  <- file.path(outdir, "logs", "main.log")
log_error <- file.path(outdir, "logs", "error.log")
summary_tsv <- file.path(outdir, "summary", "cohort_summary.tsv")
dir.create(dirname(summary_tsv), recursive=TRUE, showWarnings=FALSE)

log <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"), paste0(..., collapse=" ")), file=log_main, append=TRUE)
elog <- function(...) cat(sprintf("[%s] ERROR: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"), paste0(..., collapse=" ")), file=log_error, append=TRUE)
# NEW: Function to log failed sample info for easy retry
log_failure <- function(row) {
    write.table(row, file = fail_log_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
}


# --- Data Readers (No changes here) ---
read_10x_h5 <- function(h5_path) Read10X_h5(h5_path)
read_mtx_dir <- function(dirpath) Read10X(data.dir = dirpath)
read_rds_object <- function(rds_path) {
  obj <- readRDS(rds_path)
  if (is.list(obj) && !is.data.frame(obj)) {
    log("RDS file contains a list; searching for a matrix element.")
    matrix_element <- NULL
    for (i in seq_along(obj)) {
      if (is(obj[[i]], "matrix") || is(obj[[i]], "Matrix")) {
        matrix_element <- obj[[i]]
        log(sprintf("Found matrix in list element %d.", i))
        break
      }
    }
    if (is.null(matrix_element)) stop("Could not find a matrix object within the list in the RDS file.")
    obj <- matrix_element
  }
  if (!is(obj, "dgCMatrix")) obj <- as(obj, "dgCMatrix")
  if (any(is.na(obj@x))) { obj@x[is.na(obj@x)] <- 0 }
  if (any(obj@x %% 1 != 0)) { obj@x <- floor(obj@x) }
  return(CreateSeuratObject(obj))
}
read_mtx_stem <- function(stem){
    find_file <- function(pattern) {
        path <- list.files(dirname(stem), pattern=pattern, full.names=TRUE)
        if (length(path) > 0) return(path[1]) else return(NA_character_)
    }
    bfile <- find_file(paste0(basename(stem),"_barcodes.tsv(.gz)?"))
    ffile <- find_file(paste0(basename(stem),"_features.tsv(.gz)?"))
    mfile <- find_file(paste0(basename(stem),"_matrix.mtx(.gz)?"))
    if (any(is.na(c(bfile, ffile, mfile)))) stop(sprintf("Could not find all MTX files for stem: %s", stem))
    counts <- ReadMtx(mtx = mfile, cells = bfile, features = ffile)
    feature_data <- fread(ffile, header = FALSE)
    gene_names <- if (ncol(feature_data) >= 2) feature_data[[2]] else feature_data[[1]]
    rownames(counts) <- make.unique(gene_names)
    return(counts)
}
read_tsv_wide_sparse <- function(path) {
  dt <- data.table::fread(path, sep="\t", header=TRUE, data.table=FALSE)
  valid_rows <- !is.na(dt[[1]]) & dt[[1]] != ""
  dt <- dt[valid_rows, ]
  rn <- dt[[1]]; dt[[1]] <- NULL
  if (any(duplicated(rn))) {
    log("Found duplicate rownames in TSV file; making them unique.")
    rn <- make.unique(rn)
  }
  mat <- as.matrix(dt)
  rownames(mat) <- rn
  Matrix::Matrix(mat, sparse = TRUE)
}

# --- QC and Filtering functions (No changes here) ---
add_qc_metrics <- function(sobj){
  mt_genes <- grep("^MT-", rownames(sobj), value=TRUE)
  if (length(mt_genes) > 0) {
      sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, features = mt_genes)
  } else {
      log("No mitochondrial genes (starting with 'MT-') found. Setting percent.mt to 0.")
      sobj[["percent.mt"]] <- 0
  }
  return(sobj)
}
apply_adaptive_scrna_filters <- function(sobj, mad_k=3, min_cells_after_filter=20) {
  meta <- sobj@meta.data
  mito_median <- median(meta$percent.mt, na.rm = TRUE)
  mito_mad <- mad(meta$percent.mt, na.rm = TRUE, constant=1.4826)
  mito_upper_thresh <- mito_median + mad_k * mito_mad
  mito_final_thresh <- min(mito_upper_thresh, 20.0)
  feat_median <- median(meta$nFeature_RNA, na.rm = TRUE)
  feat_mad <- mad(meta$nFeature_RNA, na.rm = TRUE, constant=1.4826)
  feat_upper_thresh <- feat_median + mad_k * feat_mad
  count_median <- median(meta$nCount_RNA, na.rm = TRUE)
  count_mad <- mad(meta$nCount_RNA, na.rm = TRUE, constant=1.4826)
  count_upper_thresh <- count_median + mad_k * count_mad
  passing_cells <- which(
    meta$nFeature_RNA > 200 & meta$nFeature_RNA < feat_upper_thresh &
    meta$nCount_RNA < count_upper_thresh &
    meta$percent.mt < mito_final_thresh
  )
  if (length(passing_cells) < min_cells_after_filter) {
    log(sprintf("WARNING: Adaptive filtering would remove all but %d cells. Skipping filtering for sample %s.", 
        length(passing_cells), sobj$sample_id[1]))
    return(sobj)
  }
  log(sprintf("QC thresholds for %s: nFeature_RNA > 200 & < %.0f | percent.mt < %.2f | nCount_RNA < %.0f",
      sobj$sample_id[1], feat_upper_thresh, mito_final_thresh, count_upper_thresh))
  return(subset(sobj, cells = colnames(sobj)[passing_cells]))
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
  return(subset(sobj, subset = doublet_class == "singlet"))
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
  write.table(data.frame(project_id=character(), sample_id=character(), cells_raw=integer(), cells_postQC=integer(),
                         median_genes=double(), median_counts=double(), median_pct_mt=double(),
                         path_rds=character(), stringsAsFactors = FALSE),
              file=summary_tsv, sep="\t", quote=FALSE, row.names=FALSE)
}

for (i in seq_len(nrow(manifest))) {
  row <- manifest[i,]
  project <- row$project_id; fmt <- row$format; src <- row$path_or_stem; sample <- row$sample_id
  rds_path <- file.path(outdir, "rds", paste0(sample, ".rds"))
  
  if (file.exists(rds_path)) {
    log(sprintf("Skipping %s, output RDS already exists.", sample))
    next
  }

  tryCatch({
    log(sprintf("Processing [%s] format=%s sample_id=%s", project, fmt, sample))

    counts_matrix <- switch(fmt,
      "10x_h5"   = read_10x_h5(src), "mtx_dir"  = read_mtx_dir(src),
      "mtx_stem" = read_mtx_stem(src), "rds" = read_rds_object(src),
      "tsv_wide" = read_tsv_wide_sparse(src), stop(sprintf("Unknown format: %s", fmt))
    )

    if (is.list(counts_matrix) && !is.data.frame(counts_matrix) && "Gene Expression" %in% names(counts_matrix)) {
        counts_matrix <- counts_matrix[['Gene Expression']]
    }

    sobj <- if (is(counts_matrix, "Seurat")) counts_matrix else CreateSeuratObject(counts_matrix, project = project)
    sobj$sample_id <- sample
    
    sobj <- add_qc_metrics(sobj)
    raw_n <- ncol(sobj)
    save_qc_violin(sobj, outdir, sample, stage = "pre_qc")
    
    sobj <- apply_adaptive_scrna_filters(sobj)
    
    if (ncol(sobj) > 50) sobj <- run_denoise(sobj)
    
    postqc_n <- ncol(sobj)
    if (postqc_n > 50) {
      save_qc_violin(sobj, outdir, sample, stage = "post_qc")
      sobj <- sct_normalize(sobj)
    }
    
    log(sprintf("Saving RDS for %s. Raw cells: %d, Post-QC cells: %d", sample, raw_n, postqc_n))
    saveRDS(sobj, rds_path, compress="xz")
    
    summary_df <- data.frame(
      project_id=project, sample_id=sample, cells_raw=raw_n, cells_postQC=postqc_n,
      median_genes=median(sobj$nFeature_RNA, na.rm=TRUE), median_counts=median(sobj$nCount_RNA, na.rm=TRUE),
      median_pct_mt=median(sobj$percent.mt, na.rm=TRUE), path_rds=rds_path, stringsAsFactors=FALSE
    )
    write.table(summary_df, file=summary_tsv, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)

  }, error=function(e){
    msg <- sprintf("Failed on sample %s (src: %s): %s", sample, src, e$message)
    elog(msg)
    # **NEW**: Log the manifest row of the failed sample
    log_failure(row)
  })
}

log(sprintf("All processing complete. Summary table at %s", summary_tsv))
RS

# --- NEW: Post-run step to generate a retry manifest ---
RETRY_MANIFEST="$BASE/scRNA_clean/scrna_manifest_retry.tsv"
if [[ -f "$FAIL_LOG" && -s "$FAIL_LOG" ]]; then
    log "Errors were detected. Generating a retry manifest..."
    # Create header for the retry manifest
    head -n 1 "$MANIFEST" > "$RETRY_MANIFEST"
    # Append the failed sample rows
    cat "$FAIL_LOG" >> "$RETRY_MANIFEST"
    log "Retry manifest created at: $RETRY_MANIFEST"
    log "You can now fix the script and re-run it using the retry manifest."
else
    log "Run completed with no errors."
    # Optional: remove the empty retry file if it exists
    rm -f "$RETRY_MANIFEST"
fi

log "scRNA ingest+QC (FINAL, RESUMABLE) finished."