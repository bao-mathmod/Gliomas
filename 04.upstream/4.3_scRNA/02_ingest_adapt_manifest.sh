#!/usr/bin/env bash
set -Euo pipefail
IFS=$'\n\t'

# ======================= CONFIG ==========================
BASE="/mnt/18T/chibao/gliomas/data/upstream/scRNA/set1_1"
OUT="$BASE"
MANIFEST="$OUT/scrna_manifest_1_1.tsv"     # <-- your NEW manifest with sample_uid column
# =========================================================

mkdir -p "$OUT"/{rds,logs,summary,plots}

FAIL_LOG="$OUT/logs/failed_samples.log"
: > "$OUT/logs/error.log"
: > "$FAIL_LOG"

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
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
manifest_path <- args[1]
outdir <- args[2]
fail_log_path <- args[3]

log_main    <- file.path(outdir, "logs", "main.log")
log_error   <- file.path(outdir, "logs", "error.log")
summary_tsv <- file.path(outdir, "summary", "cohort_summary.tsv")
dir.create(dirname(summary_tsv), recursive=TRUE, showWarnings=FALSE)

log  <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"), paste0(..., collapse=" ")), file=log_main,  append=TRUE)
elog <- function(...) cat(sprintf("[%s] ERROR: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"), paste0(..., collapse=" ")), file=log_error, append=TRUE)

# Write raw manifest row so you can regenerate a retry manifest later
log_failure <- function(row) {
  write.table(row, file = fail_log_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
}

# ---------- Readers ----------
read_10x_h5 <- function(h5_path) Read10X_h5(h5_path)
read_mtx_dir <- function(dirpath) Read10X(data.dir = dirpath)
read_rds_object <- function(rds_path) {
  obj <- readRDS(rds_path)
  if (is.list(obj) && !is.data.frame(obj)) {
    log("RDS file contains a list; searching for a matrix element.")
    matrix_element <- NULL
    for (i in seq_along(obj)) {
      if (is(obj[[i]], "matrix") || is(obj[[i]], "Matrix")) {
        matrix_element <- obj[[i]]; log(sprintf("Found matrix in list element %d.", i)); break
      }
    }
    if (is.null(matrix_element)) stop("Could not find a matrix object within the list in the RDS file.")
    obj <- matrix_element
  }
  if (!is(obj, "dgCMatrix")) obj <- as(obj, "dgCMatrix")
  if (any(is.na(obj@x)))      obj@x[is.na(obj@x)] <- 0
  if (any(obj@x %% 1 != 0))   obj@x <- floor(obj@x)
  CreateSeuratObject(obj)
}

read_mtx_stem <- function(stem){
  # helper to find files that may or may not be gzipped
  find_file <- function(pattern) {
    path <- list.files(dirname(stem), pattern = pattern, full.names = TRUE)
    if (length(path) > 0) path[1] else NA_character_
  }

  # accept both prefixed and unprefixed 10x-style names
  bn <- basename(stem)
  bfile <- find_file(paste0("(", bn, "_)?barcodes\\.tsv(\\.gz)?$"))
  ffile <- find_file(paste0("(", bn, "_)?features\\.tsv(\\.gz)?$|(", bn, "_)?genes\\.tsv(\\.gz)?$"))
  mfile <- find_file(paste0("(", bn, "_)?matrix\\.mtx(\\.gz)?$"))
  if (any(is.na(c(bfile, ffile, mfile))))
    stop(sprintf("Could not find all MTX files for stem: %s", stem))

  # detect how many columns the features file has (fast, O(1) I/O)
  feat_dt  <- data.table::fread(ffile, header = FALSE, nrows = 3, data.table = FALSE)
  feat_cols <- ncol(feat_dt)
  feat_col  <- if (feat_cols >= 2) 2 else 1
  message(sprintf("read_mtx_stem(): %s -> features columns=%d, using column %d", ffile, feat_cols, feat_col))

  # read counts; explicitly set feature & cell columns
  counts <- ReadMtx(
    mtx = mfile, cells = bfile, features = ffile,
    feature.column = feat_col, cell.column = 1, unique.features = FALSE
  )

  # re-read only needed feature column to set human-friendly rownames, then de-duplicate
  full_feat <- data.table::fread(ffile, header = FALSE, data.table = FALSE, select = feat_col)
  gene_names <- make.unique(as.character(full_feat[[1]]))

  # clean potential CR line endings (rare but happens on Windows-generated files)
  gene_names <- sub("\r$", "", gene_names)

  rownames(counts) <- gene_names
  counts
}

read_tsv_wide_sparse <- function(path) {
  dt <- data.table::fread(path, sep="\t", header=TRUE, data.table=FALSE)
  valid_rows <- !is.na(dt[[1]]) & dt[[1]] != ""; dt <- dt[valid_rows, ]
  rn <- dt[[1]]; dt[[1]] <- NULL
  if (any(duplicated(rn))) { log("Found duplicate rownames in TSV; making them unique."); rn <- make.unique(rn) }
  mat <- as.matrix(dt); rownames(mat) <- rn
  Matrix::Matrix(mat, sparse = TRUE)
}

# ---------- QC / Filtering ----------
add_qc_metrics <- function(sobj){
  mt_genes <- grep("^MT-", rownames(sobj), value=TRUE)
  if (length(mt_genes) > 0) {
    sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, features = mt_genes)
  } else {
    log("No mitochondrial genes (starting with 'MT-') found. Setting percent.mt to 0.")
    sobj[["percent.mt"]] <- 0
  }
  sobj
}
apply_adaptive_scrna_filters <- function(sobj, mad_k=3, min_cells_after_filter=20) {
  meta <- sobj@meta.data
  mito_median <- median(meta$percent.mt, na.rm = TRUE)
  mito_mad    <- mad(meta$percent.mt, na.rm = TRUE, constant=1.4826)
  mito_upper  <- mito_median + mad_k * mito_mad
  mito_final  <- min(mito_upper, 20.0)

  feat_median <- median(meta$nFeature_RNA, na.rm = TRUE)
  feat_mad    <- mad(meta$nFeature_RNA, na.rm = TRUE, constant=1.4826)
  feat_upper  <- feat_median + mad_k * feat_mad

  count_median <- median(meta$nCount_RNA, na.rm = TRUE)
  count_mad    <- mad(meta$nCount_RNA, na.rm = TRUE, constant=1.4826)
  count_upper  <- count_median + mad_k * count_mad

  passing_cells <- which(
    meta$nFeature_RNA > 200 & meta$nFeature_RNA < feat_upper &
    meta$nCount_RNA  < count_upper &
    meta$percent.mt  < mito_final
  )
  if (length(passing_cells) < min_cells_after_filter) {
    log(sprintf("WARNING: Adaptive filtering would leave %d cells. Skipping filtering for sample %s.",
                length(passing_cells), sobj$sample_uid[1]))
    return(sobj)
  }
  log(sprintf("QC thresholds for %s: nFeature_RNA > 200 & < %.0f | percent.mt < %.2f | nCount_RNA < %.0f",
              sobj$sample_uid[1], feat_upper, mito_final, count_upper))
  subset(sobj, cells = colnames(sobj)[passing_cells])
}
save_qc_violin <- function(sobj, outdir, sample_uid, stage = c("pre_qc","post_qc")){
  pdir <- file.path(outdir, "plots", sample_uid)
  dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
  feats <- c("nFeature_RNA","nCount_RNA","percent.mt")
  for (feat in feats) {
    p <- VlnPlot(sobj, features = feat, pt.size = 0) + NoLegend() +
         ggtitle(paste0(sample_uid, " • ", feat, " (", stage, ")"))
    ggsave(file.path(pdir, paste0(stage, "_", feat, ".png")), plot=p, width=6, height=5, dpi=100)
  }
}
run_denoise <- function(sobj){
  sce <- as.SingleCellExperiment(sobj)
  dbl <- scDblFinder::scDblFinder(sce, returnType="table")
  sobj$doublet_class <- dbl[colnames(sobj), "class"]
  subset(sobj, subset = doublet_class == "singlet")
}
sct_normalize <- function(sobj){
  suppressWarnings({
    sobj <- SCTransform(sobj, vst.flavor="v2", verbose=FALSE, vars.to.regress = "percent.mt")
    sobj <- RunPCA(sobj, verbose=FALSE)
    sobj <- RunUMAP(sobj, dims = 1:30, verbose=FALSE)
  })
  sobj
}

# ---------- Manifest ingest & UID handling ----------
manifest <- data.table::fread(manifest_path, sep="\t", header=TRUE, data.table=FALSE)

required <- c("project_id","format","path_or_stem","sample_id")
missing  <- setdiff(required, names(manifest))
if (length(missing) > 0) stop(sprintf("Manifest missing required columns: %s", paste(missing, collapse=", ")))

sanitize <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
extract_leaf <- function(path) {
  m <- regexpr("SAMN\\d+", path, perl=TRUE)
  if (m[1] != -1) regmatches(path, m)[1] else basename(dirname(path))
}

# Prefer user-provided sample_uid. If absent, build fallback and make.unique.
if (!("sample_uid" %in% names(manifest))) {
  leaf  <- vapply(manifest$path_or_stem, extract_leaf, character(1))
  base  <- paste(manifest$project_id, manifest$sample_id, leaf, sep="__")
  base  <- sanitize(base)
  manifest$sample_uid <- make.unique(base, sep="_v")
} else {
  # If present but some rows empty, fill those only
  need <- which(is.na(manifest$sample_uid) | manifest$sample_uid == "")
  if (length(need) > 0) {
    leaf  <- vapply(manifest$path_or_stem[need], extract_leaf, character(1))
    base  <- paste(manifest$project_id[need], manifest$sample_id[need], leaf, sep="__")
    base  <- sanitize(base)
    manifest$sample_uid[need] <- make.unique(base, sep="_v")
  }
  # Sanitize everything
  manifest$sample_uid <- sanitize(manifest$sample_uid)
  # Enforce global uniqueness just in case
  manifest$sample_uid <- make.unique(manifest$sample_uid, sep="_v")
}

# Ensure per-project RDS folders exist
for (p in unique(manifest$project_id)) {
  dir.create(file.path(outdir, "rds", sanitize(p)), recursive = TRUE, showWarnings = FALSE)
}

# Summary header
if (!file.exists(summary_tsv)) {
  write.table(
    data.frame(project_id=character(), sample_uid=character(), orig_sample_id=character(),
               genome=character(), chemistry=character(),
               cells_raw=integer(), cells_postQC=integer(),
               median_genes=double(), median_counts=double(), median_pct_mt=double(),
               path_rds=character(), stringsAsFactors=FALSE),
    file=summary_tsv, sep="\t", quote=FALSE, row.names=FALSE
  )
}

# ---------- Processing loop ----------
for (i in seq_len(nrow(manifest))) {
  row     <- manifest[i,]
  project <- as.character(row$project_id)
  fmt     <- as.character(row$format)
  src     <- as.character(row$path_or_stem)
  sample  <- as.character(row$sample_id)
  suid    <- as.character(row$sample_uid)
  genome  <- if ("genome"    %in% names(row)) as.character(row$genome)    else NA_character_
  chem    <- if ("chemistry" %in% names(row)) as.character(row$chemistry) else NA_character_

  project_s <- sanitize(project)
  rds_path  <- file.path(outdir, "rds", project_s, paste0(suid, ".rds"))

  if (file.exists(rds_path)) { log(sprintf("Skipping %s (project %s): exists %s", suid, project, rds_path)); next }

  dir.create(dirname(rds_path), recursive = TRUE, showWarnings = FALSE)

  tryCatch({
    log(sprintf("Processing [%s] format=%s sample_uid=%s (orig=%s)", project, fmt, suid, sample))

    counts_matrix <- switch(fmt,
      "10x_h5"   = read_10x_h5(src),
      "mtx_dir"  = read_mtx_dir(src),
      "mtx_stem" = read_mtx_stem(src),
      "rds"      = read_rds_object(src),
      "tsv_wide" = read_tsv_wide_sparse(src),
      stop(sprintf("Unknown format: %s", fmt))
    )

    if (is.list(counts_matrix) && !is.data.frame(counts_matrix) && "Gene Expression" %in% names(counts_matrix)) {
      counts_matrix <- counts_matrix[['Gene Expression']]
    }

    sobj <- if (is(counts_matrix, "Seurat")) counts_matrix else CreateSeuratObject(counts_matrix, project = project_s)
    sobj$project_id     <- project
    sobj$orig_sample_id <- sample
    sobj$sample_uid     <- suid
    if (!is.na(genome)) sobj$genome <- genome
    if (!is.na(chem))   sobj$chemistry <- chem

    sobj <- add_qc_metrics(sobj)
    raw_n <- ncol(sobj)
    save_qc_violin(sobj, outdir, suid, stage = "pre_qc")

    sobj <- apply_adaptive_scrna_filters(sobj)

    if (ncol(sobj) > 50) sobj <- run_denoise(sobj)

    postqc_n <- ncol(sobj)
    if (postqc_n > 50) {
      save_qc_violin(sobj, outdir, suid, stage = "post_qc")
      sobj <- sct_normalize(sobj)
    }

    log(sprintf("Saving RDS for %s. Raw cells: %d, Post-QC cells: %d", suid, raw_n, postqc_n))
    saveRDS(sobj, rds_path, compress="xz")

    summary_df <- data.frame(
      project_id   = project,
      sample_uid   = suid,
      orig_sample_id = sample,
      genome       = genome,
      chemistry    = chem,
      cells_raw    = raw_n,
      cells_postQC = postqc_n,
      median_genes = median(sobj$nFeature_RNA, na.rm=TRUE),
      median_counts= median(sobj$nCount_RNA, na.rm=TRUE),
      median_pct_mt= median(sobj$percent.mt, na.rm=TRUE),
      path_rds     = rds_path,
      stringsAsFactors = FALSE
    )
    write.table(summary_df, file=summary_tsv, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)

  }, error=function(e){
    msg <- sprintf("Failed on sample_uid %s (src: %s): %s", suid, src, e$message)
    elog(msg)
    log_failure(row)
  })
}

log(sprintf("All processing complete. Summary table at %s", summary_tsv))
RS

# --- Post-run: build a retry manifest if any failures occurred ---
RETRY_MANIFEST="$OUT/scrna_manifest_retry.tsv"
if [[ -f "$FAIL_LOG" && -s "$FAIL_LOG" ]]; then
  log "Errors were detected. Generating a retry manifest..."
  head -n 1 "$MANIFEST" > "$RETRY_MANIFEST"
  cat "$FAIL_LOG" >> "$RETRY_MANIFEST"
  log "Retry manifest created at: $RETRY_MANIFEST"
else
  log "Run completed with no errors."
  rm -f "$RETRY_MANIFEST"
fi

log "scRNA ingest+QC (FINAL, RESUMABLE) finished."
