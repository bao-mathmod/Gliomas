#!/usr/bin/env bash
# set -Euo pipefail
# IFS=$'\n\t'

# # ======================= CONFIG ==========================
# # Adjusted paths for your snRNA data
# BASE="/mnt/18T/chibao/gliomas/data/upstream/snRNA/snRNA_adaptive"
# OUT="$BASE/snRNA_processed" # New output directory for processed results
# MANIFEST="/mnt/18T/chibao/gliomas/data/upstream/snRNA/snRNA_clean_2/snrna_manifest.tsv" # Assuming your manifest is here
# # =========================================================

# mkdir -p "$OUT"/{rds,logs,summary,plots}

# log(){ printf "[%s] %s\n" "$(date '+%Y-%m-%d %H:%M:%S%z')" "$*" | tee -a "$OUT/logs/main.log" ; }
# err(){ printf "[%s] ERROR: %s\n" "$(date '+%Y-%m-%d %H:%M:%S%z')" "$*" | tee -a "$OUT/logs/error.log" >&2 ; }

# log "snRNA ingest+QC (ADAPTIVE) started. MANIFEST=$MANIFEST OUT=$OUT"

# Rscript --vanilla - "$MANIFEST" "$OUT" <<'RS' || { err "Rscript failed"; exit 1; }
# #!/usr/bin/env Rscript
# suppressPackageStartupMessages({
#   library(Seurat)
#   library(Matrix)
#   library(data.table)
#   library(dplyr)
#   library(scDblFinder)
#   library(SingleCellExperiment)
#   library(ggplot2)
# })

# args <- commandArgs(trailingOnly = TRUE)
# manifest_path <- args[1]
# outdir <- args[2]
# log_main  <- file.path(outdir, "logs", "main.log")
# log_error <- file.path(outdir, "logs", "error.log")
# summary_tsv <- file.path(outdir, "summary", "cohort_summary.tsv")
# dir.create(dirname(summary_tsv), recursive=TRUE, showWarnings=FALSE)

# # --- Logging helpers ---
# log <- function(...) {
#   cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"), paste0(..., collapse=" ")),
#       file=log_main, append=TRUE)
# }
# elog <- function(...) {
#   cat(sprintf("[%s] ERROR: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"), paste0(..., collapse=" ")),
#       file=log_error, append=TRUE)
# }

# # --- Data Readers (assuming formats from your previous snRNA script) ---
# read_10x_h5 <- function(h5_path) Read10X_h5(h5_path)
# read_mtx_dir <- function(dirpath) Read10X(data.dir = dirpath)
# read_mtx_stem <- function(stem){
#     bfile <- list.files(dirname(stem), pattern=paste0(basename(stem),"_barcodes.tsv"), full.names=TRUE)[1]
#     ffile <- list.files(dirname(stem), pattern=paste0(basename(stem),"_features.tsv"), full.names=TRUE)[1]
#     mfile <- list.files(dirname(stem), pattern=paste0(basename(stem),"_matrix.mtx"), full.names=TRUE)[1]
#     if (is.na(bfile)) bfile <- list.files(dirname(stem), pattern=paste0(basename(stem),"_barcodes.tsv.gz"), full.names=TRUE)[1]
#     if (is.na(ffile)) ffile <- list.files(dirname(stem), pattern=paste0(basename(stem),"_features.tsv.gz"), full.names=TRUE)[1]
#     if (is.na(mfile)) mfile <- list.files(dirname(stem), pattern=paste0(basename(stem),"_matrix.mtx.gz"), full.names=TRUE)[1]
#     if (any(is.na(c(bfile, ffile, mfile)))) stop(sprintf("Could not find all MTX files for stem: %s", stem))
#     counts <- ReadMtx(mtx = mfile, features = ffile, cells = bfile)
#     return(counts)
# }
# read_tsv_wide_sparse <- function(path) {
#   dt <- data.table::fread(path, sep="\t", header=TRUE, data.table=FALSE)
#   rn <- dt[[1]]; dt[[1]] <- NULL
#   mat <- as.matrix(dt); rownames(mat) <- rn
#   Matrix::Matrix(mat, sparse = TRUE)
# }

# # --- QC and Filtering functions ---
# add_qc_metrics <- function(sobj){
#   sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
#   return(sobj)
# }

# # **NEW**: Adaptive filtering function tailored for snRNA-seq data
# apply_adaptive_snrna_filters <- function(sobj, mad_k=3) {
#   meta <- sobj@meta.data
  
#   # 1. Filter by mitochondrial percentage (STRICT FIXED CUTOFF)
#   # For nuclei, this is a purity metric. High mito % means cytoplasmic contamination.
#   mito_final_thresh <- 1.0
  
#   # 2. Filter by feature count (nFeature_RNA) (adaptive upper bound)
#   feat_median <- median(meta$nFeature_RNA, na.rm = TRUE)
#   feat_mad <- mad(meta$nFeature_RNA, na.rm = TRUE)
#   feat_lower_thresh <- 200 # Keep a safe minimum
#   feat_upper_thresh <- feat_median + mad_k * feat_mad
  
#   # 3. Filter by total counts (nCount_RNA) (adaptive upper bound)
#   count_median <- median(meta$nCount_RNA, na.rm = TRUE)
#   count_mad <- mad(meta$nCount_RNA, na.rm = TRUE)
#   count_upper_thresh <- count_median + mad_k * count_mad
  
#   log(sprintf("QC thresholds for %s: nFeature_RNA > %.0f & < %.0f | percent.mt < %.2f | nCount_RNA < %.0f",
#       sobj$sample_id[1], feat_lower_thresh, feat_upper_thresh, mito_final_thresh, count_upper_thresh))
  
#   subset(sobj, subset = 
#     nFeature_RNA > feat_lower_thresh & nFeature_RNA < feat_upper_thresh &
#     nCount_RNA < count_upper_thresh &
#     percent.mt < mito_final_thresh
#   )
# }

# save_qc_violin <- function(sobj, outdir, sample_id, stage = c("pre_qc","post_qc")){
#   pdir <- file.path(outdir, "plots", sample_id)
#   dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
#   feats <- c("nFeature_RNA","nCount_RNA","percent.mt")
#   for (feat in feats) {
#     p <- VlnPlot(sobj, features = feat, pt.size = 0) + NoLegend() +
#          ggtitle(paste0(sample_id, " • ", feat, " (", stage, ")"))
#     ggsave(file.path(pdir, paste0(stage, "_", feat, ".png")), plot=p, width=6, height=5, dpi=100)
#   }
# }

# run_denoise <- function(sobj){
#   sce <- as.SingleCellExperiment(sobj)
#   dbl_calls <- scDblFinder::scDblFinder(sce, returnType="table")
#   sobj$doublet_class  <- dbl_calls[colnames(sobj), "class"]
#   subset(sobj, subset = doublet_class == "singlet")
# }

# sct_normalize <- function(sobj){
#   suppressWarnings({
#     sobj <- SCTransform(sobj, vst.flavor="v2", verbose=FALSE, vars.to.regress = "percent.mt")
#     sobj <- RunPCA(sobj, verbose=FALSE)
#     sobj <- RunUMAP(sobj, dims = 1:30, verbose=FALSE)
#   })
#   return(sobj)
# }

# # --- Main processing loop ---
# manifest <- data.table::fread(manifest_path, sep="\t", header=TRUE, data.table=FALSE)
# if (!file.exists(summary_tsv)) {
#   write.table(data.frame(project_id=character(), sample_id=character(),
#                          cells_raw=integer(), cells_postQC=integer(),
#                          median_genes=double(), median_counts=double(), median_pct_mt=double(),
#                          path_rds=character(), stringsAsFactors = FALSE),
#               file=summary_tsv, sep="\t", quote=FALSE, row.names=FALSE)
# }

# for (i in seq_len(nrow(manifest))) {
#   row <- manifest[i,]
#   project <- row$project_id
#   fmt     <- row$format
#   src     <- row$path_or_stem
#   sample  <- row$sample_id
#   rds_path <- file.path(outdir, "rds", paste0(sample, ".rds"))
  
#   tryCatch({
#     log(sprintf("Processing [%s] format=%s sample_id=%s", project, fmt, sample))

#     counts_matrix <- switch(fmt,
#       "10x_h5"   = read_10x_h5(src),
#       "mtx_dir"  = read_mtx_dir(src),
#       "mtx_stem" = read_mtx_stem(src),
#       "mtx"      = read_mtx_stem(src), # For backward compatibility with old manifest
#       "tsv_wide" = read_tsv_wide_sparse(src),
#       stop(sprintf("Unknown format: %s", fmt))
#     )
#     if (is.list(counts_matrix) && !is.data.frame(counts_matrix)) {
#         counts_matrix <- counts_matrix[[1]]
#     }
#     sobj <- CreateSeuratObject(counts_matrix, project = project)
#     sobj$sample_id <- sample
    
#     sobj <- add_qc_metrics(sobj)
#     raw_n <- ncol(sobj)
#     save_qc_violin(sobj, outdir, sample, stage = "pre_qc")
    
#     sobj <- apply_adaptive_snrna_filters(sobj)
    
#     if (ncol(sobj) > 50) {
#         sobj <- run_denoise(sobj)
#     } else {
#         log(sprintf("Skipping doublet removal for %s, too few nuclei post-QC (%d)", sample, ncol(sobj)))
#     }
    
#     postqc_n <- ncol(sobj)
#     if (postqc_n > 50) {
#       save_qc_violin(sobj, outdir, sample, stage = "post_qc")
#       sobj <- sct_normalize(sobj)
#     } else {
#       log(sprintf("Skipping normalization for %s, too few nuclei remaining (%d)", sample, postqc_n))
#     }
    
#     log(sprintf("Saving RDS for %s. Raw nuclei: %d, Post-QC nuclei: %d", sample, raw_n, postqc_n))
#     saveRDS(sobj, rds_path, compress="xz")
    
#     summary_df <- data.frame(
#       project_id = project, sample_id  = sample,
#       cells_raw  = raw_n, cells_postQC = postqc_n,
#       median_genes = median(sobj$nFeature_RNA, na.rm=TRUE),
#       median_counts = median(sobj$nCount_RNA, na.rm=TRUE),
#       median_pct_mt = median(sobj$percent.mt, na.rm=TRUE),
#       path_rds = rds_path,
#       stringsAsFactors = FALSE
#     )
#     write.table(summary_df, file=summary_tsv, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)

#   }, error=function(e){
#     msg <- sprintf("Failed on sample %s (src: %s): %s", sample, src, e$message)
#     elog(msg)
#   })
# }

# log(sprintf("All processing complete. Summary table at %s", summary_tsv))
# RS

# log "snRNA ingest+QC (ADAPTIVE) finished."

#!/usr/bin/env bash
set -Euo pipefail
IFS=$'\n\t'

# ======================= CONFIG ==========================
# Input manifest (MATCHES A)
MANIFEST="/mnt/18T/chibao/gliomas/data/upstream/snRNA/official/snrna_official_manifest.tsv"

# Output base
BASE="/mnt/18T/chibao/gliomas/data/upstream/snRNA/official"
OUT="$BASE"
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
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
manifest_path <- args[1]
outdir <- args[2]

log_main  <- file.path(outdir, "logs", "main.log")
log_error <- file.path(outdir, "logs", "error.log")
summary_tsv <- file.path(outdir, "summary", "cohort_summary.tsv")
dir.create(dirname(summary_tsv), recursive=TRUE, showWarnings=FALSE)

log <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"), paste0(..., collapse=" ")), file=log_main,  append=TRUE)
elog<- function(...) cat(sprintf("[%s] ERROR: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"), paste0(..., collapse=" ")), file=log_error, append=TRUE)

# ---------- Readers ----------
read_10x_h5   <- function(h5_path) Read10X_h5(h5_path)
read_mtx_dir  <- function(dirpath) Read10X(data.dir = dirpath)

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

  # detect how many columns the features file has (fast)
  feat_dt <- data.table::fread(ffile, header = FALSE, nrows = 3, data.table = FALSE)
  feat_cols <- ncol(feat_dt)
  feat_col  <- if (feat_cols >= 2) 2 else 1

  # read counts; explicitly set feature.column & cell.column to be safe
  counts <- ReadMtx(
    mtx = mfile, cells = bfile, features = ffile,
    feature.column = feat_col, cell.column = 1, unique.features = FALSE
  )

  # enforce unique gene names (and choose human-friendly names)
  full_feat <- data.table::fread(ffile, header = FALSE, data.table = FALSE)
  gene_names <- if (ncol(full_feat) >= 2) full_feat[[2]] else full_feat[[1]]
  rownames(counts) <- make.unique(gene_names)
  counts
}

read_tsv_wide_sparse <- function(path) {
  dt <- data.table::fread(path, sep="\t", header=TRUE, data.table=FALSE)
  valid <- !is.na(dt[[1]]) & nzchar(dt[[1]])
  dt <- dt[valid, ]
  rn <- dt[[1]]; dt[[1]] <- NULL
  if (any(duplicated(rn))) rn <- make.unique(rn)
  mat <- as.matrix(dt); rownames(mat) <- rn
  Matrix::Matrix(mat, sparse = TRUE)
}

# ---------- QC helpers (snRNA friendly) ----------
add_qc_metrics <- function(sobj){
  # 1. Calculate percent.mt
  mt_genes <- grep("^MT-", rownames(sobj), value=TRUE)
  if (length(mt_genes) > 0) {
    sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, features = mt_genes)
  } else {
    # snRNA often lacks MT- genes; avoid NA/NaN flows
    sobj[["percent.mt"]] <- 0
    log(sprintf("No MT- genes found for %s; setting percent.mt=0.", sobj$sample_uid[1]))
  }
  
  # 2. Calculate percent.ribo (MOVED INSIDE)
  # Use ignore.case=TRUE for robustness (e.g., Rps, rpl)
  ribo_genes <- grep("^RP[SL][0-9]", rownames(sobj), value = TRUE, ignore.case = TRUE)
  if (length(ribo_genes) > 0) {
    sobj[["percent.ribo"]] <- PercentageFeatureSet(sobj, features = ribo_genes)
  } else {
    sobj[["percent.ribo"]] <- 0
    log(sprintf("No Ribo genes (RP[SL]) found for %s; setting percent.ribo=0.", sobj$sample_uid[1]))
  }
  
  # 3. Return the object
  sobj
}

# snRNA rule-of-thumb: strict mito; adaptive high-end gates for features/counts
apply_adaptive_snrna_filters <- function(sobj, mad_k=3, min_cells_after_filter=20) {
  md <- sobj@meta.data
  mito_median <- median(md$percent.mt, na.rm=TRUE)
  # keep strict cap at 1%; if MT all zeros, this remains 0–1%
  mito_final  <- 1.0

  feat_median <- median(md$nFeature_RNA, na.rm=TRUE)
  feat_mad    <- mad(md$nFeature_RNA, na.rm=TRUE, constant=1.4826)
  feat_lower  <- 200
  feat_upper  <- feat_median + mad_k * feat_mad

  cnt_median  <- median(md$nCount_RNA, na.rm=TRUE)
  cnt_mad     <- mad(md$nCount_RNA, na.rm=TRUE, constant=1.4826)
  cnt_upper   <- cnt_median + mad_k * cnt_mad

  passing <- which(
    md$nFeature_RNA > feat_lower & md$nFeature_RNA < feat_upper &
    md$nCount_RNA   < cnt_upper  &
    md$percent.mt   < mito_final
  )

  if (length(passing) < min_cells_after_filter) {
    log(sprintf("WARNING: Filtering would leave %d nuclei in %s; skipping filtering.",
                length(passing), sobj$sample_uid[1]))
    return(sobj)
  }

  log(sprintf("QC %s: nFeature_RNA > %d & < %.0f | nCount_RNA < %.0f | percent.mt < %.2f (median mt=%.2f)",
              sobj$sample_uid[1], feat_lower, feat_upper, cnt_upper, mito_final, mito_median))
  subset(sobj, cells = colnames(sobj)[passing])
}

save_qc_violin <- function(sobj, outdir, suid, stage=c("pre_qc","post_qc")){
  pdir <- file.path(outdir, "plots", suid)
  dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
  for (feat in c("nFeature_RNA","nCount_RNA","percent.mt", "percent.ribo")) {
    p <- VlnPlot(sobj, features=feat, pt.size=0) + NoLegend() + ggtitle(paste0(suid," • ",feat," (",stage,")"))
    ggsave(file.path(pdir, paste0(stage,"_",feat,".png")), plot=p, width=6, height=5, dpi=100)
  }
}

run_denoise <- function(sobj){
  sce <- as.SingleCellExperiment(sobj)
  tbl <- scDblFinder::scDblFinder(sce, returnType="table")
  sobj$doublet_class <- tbl[colnames(sobj), "class"]
  subset(sobj, subset = doublet_class == "singlet")
}

sct_normalize <- function(sobj){
  suppressWarnings({
    sobj <- SCTransform(sobj, vst.flavor="v2", verbose=FALSE, vars.to.regress=c("percent.mt", "percent.ribo"))
    sobj <- RunPCA(sobj, verbose=FALSE)
    sobj <- RunUMAP(sobj, dims=1:30, verbose=FALSE)
  })
  sobj
}

# ---------- Manifest & ID handling ----------
manifest <- data.table::fread(manifest_path, sep="\t", header=TRUE, data.table=FALSE)
req <- c("project_id","format","path_or_stem","sample_id")
miss <- setdiff(req, names(manifest))
if (length(miss) > 0) stop(sprintf("Manifest missing required columns: %s", paste(miss, collapse=", ")))

sanitize <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
extract_leaf <- function(path) {
  m <- regexpr("SAMN\\d+", path, perl=TRUE)
  if (m[1] != -1) regmatches(path, m)[1] else basename(dirname(path))
}

# Prefer provided sample_uid; otherwise build & make.unique
if (!("sample_uid" %in% names(manifest))) {
  leaf <- vapply(manifest$path_or_stem, extract_leaf, character(1))
  built <- sanitize(paste(manifest$project_id, manifest$sample_id, leaf, sep="__"))
  manifest$sample_uid <- make.unique(built, sep="_v")
} else {
  # Fill empties, sanitize, ensure uniqueness
  need <- which(is.na(manifest$sample_uid) | manifest$sample_uid == "")
  if (length(need) > 0) {
    leaf <- vapply(manifest$path_or_stem[need], extract_leaf, character(1))
    built <- sanitize(paste(manifest$project_id[need], manifest$sample_id[need], leaf, sep="__"))
    manifest$sample_uid[need] <- make.unique(built, sep="_v")
  }
  manifest$sample_uid <- make.unique(sanitize(manifest$sample_uid), sep="_v")
}

# Ensure per-project RDS folders
for (p in unique(manifest$project_id)) {
  dir.create(file.path(outdir, "rds", sanitize(p)), recursive=TRUE, showWarnings=FALSE)
}

# Write summary header if missing
if (!file.exists(summary_tsv)) {
  write.table(
    data.frame(project_id=character(), sample_uid=character(), orig_sample_id=character(),
               genome=character(), chemistry=character(),
               cells_raw=integer(), cells_postQC=integer(),
               median_genes=double(), median_counts=double(), median_pct_mt=double(),
               median_pct_ribo=double(), # <-- ADD THIS COLUMN
               path_rds=character(), stringsAsFactors=FALSE),
    file=summary_tsv, sep="\t", quote=FALSE, row.names=FALSE
  )
}

# ---------- Main loop ----------
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i,]
  project <- as.character(row$project_id)
  fmt     <- as.character(row$format)
  src     <- as.character(row$path_or_stem)
  sample  <- as.character(row$sample_id)
  suid    <- as.character(row$sample_uid %||% sample)
  genome  <- if ("genome"    %in% names(row)) as.character(row$genome)    else NA_character_
  chem    <- if ("chemistry" %in% names(row)) as.character(row$chemistry) else NA_character_

  project_s <- sanitize(project)
  rds_path  <- file.path(outdir, "rds", project_s, paste0(suid, ".rds"))

  if (file.exists(rds_path)) { log(sprintf("Skip %s (exists: %s)", suid, rds_path)); next }

  dir.create(dirname(rds_path), recursive=TRUE, showWarnings=FALSE)

  tryCatch({
    log(sprintf("Processing [%s] fmt=%s uid=%s (orig=%s)", project, fmt, suid, sample))

    counts <- switch(fmt,
      "10x_h5"   = read_10x_h5(src),
      "mtx_dir"  = read_mtx_dir(src),
      "mtx_stem" = read_mtx_stem(src),
      "tsv_wide" = read_tsv_wide_sparse(src),
      # backward compatibility if an old manifest said "mtx"
      "mtx"      = read_mtx_stem(src),
      stop(sprintf("Unknown format: %s", fmt))
    )
    if (is.list(counts) && !is.data.frame(counts)) {
      # 10x Read10X can return a list; prefer Gene Expression if present else first
      counts <- if ("Gene Expression" %in% names(counts)) counts[["Gene Expression"]] else counts[[1]]
    }

    sobj <- CreateSeuratObject(counts, project=project_s)
    sobj$project_id     <- project
    sobj$orig_sample_id <- sample
    sobj$sample_uid     <- suid
    if (!is.na(genome)) sobj$genome <- genome
    if (!is.na(chem))   sobj$chemistry <- chem

    sobj <- add_qc_metrics(sobj)
    raw_n <- ncol(sobj)
    save_qc_violin(sobj, outdir, suid, stage="pre_qc")

    sobj <- apply_adaptive_snrna_filters(sobj)

    if (ncol(sobj) > 50) {
      sobj <- run_denoise(sobj)
    } else {
      log(sprintf("Skip doublet removal: too few nuclei in %s (%d)", suid, ncol(sobj)))
    }

    postqc_n <- ncol(sobj)
    if (postqc_n > 50) {
      save_qc_violin(sobj, outdir, suid, stage="post_qc")
      sobj <- sct_normalize(sobj)
    } else {
      log(sprintf("Skip normalization: too few nuclei in %s (%d)", suid, postqc_n))
    }

    log(sprintf("Saving %s | raw=%d postQC=%d", suid, raw_n, postqc_n))
    saveRDS(sobj, rds_path, compress="xz")

    summary_df <- data.frame(
      project_id   = project,
      sample_uid   = suid,
      orig_sample_id = sample,
      genome       = genome,
      chem         = chem,
      cells_raw    = raw_n,
      cells_postQC = postqc_n,
      median_genes = median(sobj$nFeature_RNA, na.rm=TRUE),
      median_counts= median(sobj$nCount_RNA, na.rm=TRUE),
      median_pct_mt= median(sobj$percent.mt, na.rm=TRUE),
      median_pct_ribo= median(sobj$percent.ribo, na.rm=TRUE), # <-- ADD THIS VALUE
      path_rds     = rds_path,
      stringsAsFactors=FALSE
    )
    write.table(summary_df, file=summary_tsv, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)

  }, error=function(e){
    elog(sprintf("Failed uid %s (src=%s): %s", suid, src, e$message))
  })
}

log(sprintf("All processing complete. Summary at %s", summary_tsv))
RS

log "snRNA ingest+QC (ADAPTIVE) finished."

