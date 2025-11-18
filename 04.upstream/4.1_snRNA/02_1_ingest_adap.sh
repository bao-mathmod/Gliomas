# #!/usr/bin/env bash
# # set -Euo pipefail
# # IFS=$'\n\t'

# # # ======================= CONFIG ==========================
# # # Adjusted paths for your snRNA data
# # BASE="/mnt/18T/chibao/gliomas/data/upstream/snRNA/snRNA_adaptive"
# # OUT="$BASE/snRNA_processed" # New output directory for processed results
# # MANIFEST="/mnt/18T/chibao/gliomas/data/upstream/snRNA/snRNA_clean_2/snrna_manifest.tsv" # Assuming your manifest is here
# # # =========================================================

# # mkdir -p "$OUT"/{rds,logs,summary,plots}

# # log(){ printf "[%s] %s\n" "$(date '+%Y-%m-%d %H:%M:%S%z')" "$*" | tee -a "$OUT/logs/main.log" ; }
# # err(){ printf "[%s] ERROR: %s\n" "$(date '+%Y-%m-%d %H:%M:%S%z')" "$*" | tee -a "$OUT/logs/error.log" >&2 ; }

# # log "snRNA ingest+QC (ADAPTIVE) started. MANIFEST=$MANIFEST OUT=$OUT"

# # Rscript --vanilla - "$MANIFEST" "$OUT" <<'RS' || { err "Rscript failed"; exit 1; }
# # #!/usr/bin/env Rscript
# # suppressPackageStartupMessages({
# #   library(Seurat)
# #   library(Matrix)
# #   library(data.table)
# #   library(dplyr)
# #   library(scDblFinder)
# #   library(SingleCellExperiment)
# #   library(ggplot2)
# # })

# # args <- commandArgs(trailingOnly = TRUE)
# # manifest_path <- args[1]
# # outdir <- args[2]
# # log_main  <- file.path(outdir, "logs", "main.log")
# # log_error <- file.path(outdir, "logs", "error.log")
# # summary_tsv <- file.path(outdir, "summary", "cohort_summary.tsv")
# # dir.create(dirname(summary_tsv), recursive=TRUE, showWarnings=FALSE)

# # # --- Logging helpers ---
# # log <- function(...) {
# #   cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"), paste0(..., collapse=" ")),
# #       file=log_main, append=TRUE)
# # }
# # elog <- function(...) {
# #   cat(sprintf("[%s] ERROR: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"), paste0(..., collapse=" ")),
# #       file=log_error, append=TRUE)
# # }

# # # --- Data Readers (assuming formats from your previous snRNA script) ---
# # read_10x_h5 <- function(h5_path) Read10X_h5(h5_path)
# # read_mtx_dir <- function(dirpath) Read10X(data.dir = dirpath)
# # read_mtx_stem <- function(stem){
# #     bfile <- list.files(dirname(stem), pattern=paste0(basename(stem),"_barcodes.tsv"), full.names=TRUE)[1]
# #     ffile <- list.files(dirname(stem), pattern=paste0(basename(stem),"_features.tsv"), full.names=TRUE)[1]
# #     mfile <- list.files(dirname(stem), pattern=paste0(basename(stem),"_matrix.mtx"), full.names=TRUE)[1]
# #     if (is.na(bfile)) bfile <- list.files(dirname(stem), pattern=paste0(basename(stem),"_barcodes.tsv.gz"), full.names=TRUE)[1]
# #     if (is.na(ffile)) ffile <- list.files(dirname(stem), pattern=paste0(basename(stem),"_features.tsv.gz"), full.names=TRUE)[1]
# #     if (is.na(mfile)) mfile <- list.files(dirname(stem), pattern=paste0(basename(stem),"_matrix.mtx.gz"), full.names=TRUE)[1]
# #     if (any(is.na(c(bfile, ffile, mfile)))) stop(sprintf("Could not find all MTX files for stem: %s", stem))
# #     counts <- ReadMtx(mtx = mfile, features = ffile, cells = bfile)
# #     return(counts)
# # }
# # read_tsv_wide_sparse <- function(path) {
# #   dt <- data.table::fread(path, sep="\t", header=TRUE, data.table=FALSE)
# #   rn <- dt[[1]]; dt[[1]] <- NULL
# #   mat <- as.matrix(dt); rownames(mat) <- rn
# #   Matrix::Matrix(mat, sparse = TRUE)
# # }

# # # --- QC and Filtering functions ---
# # add_qc_metrics <- function(sobj){
# #   sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
# #   return(sobj)
# # }

# # # **NEW**: Adaptive filtering function tailored for snRNA-seq data
# # apply_adaptive_snrna_filters <- function(sobj, mad_k=3) {
# #   meta <- sobj@meta.data
  
# #   # 1. Filter by mitochondrial percentage (STRICT FIXED CUTOFF)
# #   # For nuclei, this is a purity metric. High mito % means cytoplasmic contamination.
# #   mito_final_thresh <- 1.0
  
# #   # 2. Filter by feature count (nFeature_RNA) (adaptive upper bound)
# #   feat_median <- median(meta$nFeature_RNA, na.rm = TRUE)
# #   feat_mad <- mad(meta$nFeature_RNA, na.rm = TRUE)
# #   feat_lower_thresh <- 200 # Keep a safe minimum
# #   feat_upper_thresh <- feat_median + mad_k * feat_mad
  
# #   # 3. Filter by total counts (nCount_RNA) (adaptive upper bound)
# #   count_median <- median(meta$nCount_RNA, na.rm = TRUE)
# #   count_mad <- mad(meta$nCount_RNA, na.rm = TRUE)
# #   count_upper_thresh <- count_median + mad_k * count_mad
  
# #   log(sprintf("QC thresholds for %s: nFeature_RNA > %.0f & < %.0f | percent.mt < %.2f | nCount_RNA < %.0f",
# #       sobj$sample_id[1], feat_lower_thresh, feat_upper_thresh, mito_final_thresh, count_upper_thresh))
  
# #   subset(sobj, subset = 
# #     nFeature_RNA > feat_lower_thresh & nFeature_RNA < feat_upper_thresh &
# #     nCount_RNA < count_upper_thresh &
# #     percent.mt < mito_final_thresh
# #   )
# # }

# # save_qc_violin <- function(sobj, outdir, sample_id, stage = c("pre_qc","post_qc")){
# #   pdir <- file.path(outdir, "plots", sample_id)
# #   dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
# #   feats <- c("nFeature_RNA","nCount_RNA","percent.mt")
# #   for (feat in feats) {
# #     p <- VlnPlot(sobj, features = feat, pt.size = 0) + NoLegend() +
# #          ggtitle(paste0(sample_id, " • ", feat, " (", stage, ")"))
# #     ggsave(file.path(pdir, paste0(stage, "_", feat, ".png")), plot=p, width=6, height=5, dpi=100)
# #   }
# # }

# # run_denoise <- function(sobj){
# #   sce <- as.SingleCellExperiment(sobj)
# #   dbl_calls <- scDblFinder::scDblFinder(sce, returnType="table")
# #   sobj$doublet_class  <- dbl_calls[colnames(sobj), "class"]
# #   subset(sobj, subset = doublet_class == "singlet")
# # }

# # sct_normalize <- function(sobj){
# #   suppressWarnings({
# #     sobj <- SCTransform(sobj, vst.flavor="v2", verbose=FALSE, vars.to.regress = "percent.mt")
# #     sobj <- RunPCA(sobj, verbose=FALSE)
# #     sobj <- RunUMAP(sobj, dims = 1:30, verbose=FALSE)
# #   })
# #   return(sobj)
# # }

# # # --- Main processing loop ---
# # manifest <- data.table::fread(manifest_path, sep="\t", header=TRUE, data.table=FALSE)
# # if (!file.exists(summary_tsv)) {
# #   write.table(data.frame(project_id=character(), sample_id=character(),
# #                          cells_raw=integer(), cells_postQC=integer(),
# #                          median_genes=double(), median_counts=double(), median_pct_mt=double(),
# #                          path_rds=character(), stringsAsFactors = FALSE),
# #               file=summary_tsv, sep="\t", quote=FALSE, row.names=FALSE)
# # }

# # for (i in seq_len(nrow(manifest))) {
# #   row <- manifest[i,]
# #   project <- row$project_id
# #   fmt     <- row$format
# #   src     <- row$path_or_stem
# #   sample  <- row$sample_id
# #   rds_path <- file.path(outdir, "rds", paste0(sample, ".rds"))
  
# #   tryCatch({
# #     log(sprintf("Processing [%s] format=%s sample_id=%s", project, fmt, sample))

# #     counts_matrix <- switch(fmt,
# #       "10x_h5"   = read_10x_h5(src),
# #       "mtx_dir"  = read_mtx_dir(src),
# #       "mtx_stem" = read_mtx_stem(src),
# #       "mtx"      = read_mtx_stem(src), # For backward compatibility with old manifest
# #       "tsv_wide" = read_tsv_wide_sparse(src),
# #       stop(sprintf("Unknown format: %s", fmt))
# #     )
# #     if (is.list(counts_matrix) && !is.data.frame(counts_matrix)) {
# #         counts_matrix <- counts_matrix[[1]]
# #     }
# #     sobj <- CreateSeuratObject(counts_matrix, project = project)
# #     sobj$sample_id <- sample
    
# #     sobj <- add_qc_metrics(sobj)
# #     raw_n <- ncol(sobj)
# #     save_qc_violin(sobj, outdir, sample, stage = "pre_qc")
    
# #     sobj <- apply_adaptive_snrna_filters(sobj)
    
# #     if (ncol(sobj) > 50) {
# #         sobj <- run_denoise(sobj)
# #     } else {
# #         log(sprintf("Skipping doublet removal for %s, too few nuclei post-QC (%d)", sample, ncol(sobj)))
# #     }
    
# #     postqc_n <- ncol(sobj)
# #     if (postqc_n > 50) {
# #       save_qc_violin(sobj, outdir, sample, stage = "post_qc")
# #       sobj <- sct_normalize(sobj)
# #     } else {
# #       log(sprintf("Skipping normalization for %s, too few nuclei remaining (%d)", sample, postqc_n))
# #     }
    
# #     log(sprintf("Saving RDS for %s. Raw nuclei: %d, Post-QC nuclei: %d", sample, raw_n, postqc_n))
# #     saveRDS(sobj, rds_path, compress="xz")
    
# #     summary_df <- data.frame(
# #       project_id = project, sample_id  = sample,
# #       cells_raw  = raw_n, cells_postQC = postqc_n,
# #       median_genes = median(sobj$nFeature_RNA, na.rm=TRUE),
# #       median_counts = median(sobj$nCount_RNA, na.rm=TRUE),
# #       median_pct_mt = median(sobj$percent.mt, na.rm=TRUE),
# #       path_rds = rds_path,
# #       stringsAsFactors = FALSE
# #     )
# #     write.table(summary_df, file=summary_tsv, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)

# #   }, error=function(e){
# #     msg <- sprintf("Failed on sample %s (src: %s): %s", sample, src, e$message)
# #     elog(msg)
# #   })
# # }

# # log(sprintf("All processing complete. Summary table at %s", summary_tsv))
# # RS

# # log "snRNA ingest+QC (ADAPTIVE) finished."

# #!/usr/bin/env bash
# set -Euo pipefail
# IFS=$'\n\t'

# # ======================= CONFIG ==========================
# # Input manifest (MATCHES A)
# MANIFEST="/mnt/18T/chibao/gliomas/data/upstream/snRNA/official/snrna_official_manifest.tsv"

# # Output base
# BASE="/mnt/18T/chibao/gliomas/data/upstream/snRNA/official"
# OUT="$BASE"
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
#   library(stringr)
# })

# args <- commandArgs(trailingOnly = TRUE)
# manifest_path <- args[1]
# outdir <- args[2]

# log_main  <- file.path(outdir, "logs", "main.log")
# log_error <- file.path(outdir, "logs", "error.log")
# summary_tsv <- file.path(outdir, "summary", "cohort_summary.tsv")
# dir.create(dirname(summary_tsv), recursive=TRUE, showWarnings=FALSE)

# log <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"), paste0(..., collapse=" ")), file=log_main,  append=TRUE)
# elog<- function(...) cat(sprintf("[%s] ERROR: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"), paste0(..., collapse=" ")), file=log_error, append=TRUE)

# # ---------- Readers ----------
# read_10x_h5   <- function(h5_path) Read10X_h5(h5_path)
# read_mtx_dir  <- function(dirpath) Read10X(data.dir = dirpath)

# read_mtx_stem <- function(stem){
#   # helper to find files that may or may not be gzipped
#   find_file <- function(pattern) {
#     path <- list.files(dirname(stem), pattern = pattern, full.names = TRUE)
#     if (length(path) > 0) path[1] else NA_character_
#   }
#   # accept both prefixed and unprefixed 10x-style names
#   bn <- basename(stem)
#   bfile <- find_file(paste0("(", bn, "_)?barcodes\\.tsv(\\.gz)?$"))
#   ffile <- find_file(paste0("(", bn, "_)?features\\.tsv(\\.gz)?$|(", bn, "_)?genes\\.tsv(\\.gz)?$"))
#   mfile <- find_file(paste0("(", bn, "_)?matrix\\.mtx(\\.gz)?$"))
#   if (any(is.na(c(bfile, ffile, mfile))))
#     stop(sprintf("Could not find all MTX files for stem: %s", stem))

#   # detect how many columns the features file has (fast)
#   feat_dt <- data.table::fread(ffile, header = FALSE, nrows = 3, data.table = FALSE)
#   feat_cols <- ncol(feat_dt)
#   feat_col  <- if (feat_cols >= 2) 2 else 1

#   # read counts; explicitly set feature.column & cell.column to be safe
#   counts <- ReadMtx(
#     mtx = mfile, cells = bfile, features = ffile,
#     feature.column = feat_col, cell.column = 1, unique.features = FALSE
#   )

#   # enforce unique gene names (and choose human-friendly names)
#   full_feat <- data.table::fread(ffile, header = FALSE, data.table = FALSE)
#   gene_names <- if (ncol(full_feat) >= 2) full_feat[[2]] else full_feat[[1]]
#   rownames(counts) <- make.unique(gene_names)
#   counts
# }

# read_tsv_wide_sparse <- function(path) {
#   dt <- data.table::fread(path, sep="\t", header=TRUE, data.table=FALSE)
#   valid <- !is.na(dt[[1]]) & nzchar(dt[[1]])
#   dt <- dt[valid, ]
#   rn <- dt[[1]]; dt[[1]] <- NULL
#   if (any(duplicated(rn))) rn <- make.unique(rn)
#   mat <- as.matrix(dt); rownames(mat) <- rn
#   Matrix::Matrix(mat, sparse = TRUE)
# }

# # ---------- QC helpers (snRNA friendly) ----------
# add_qc_metrics <- function(sobj){
#   # 1. Calculate percent.mt
#   mt_genes <- grep("^MT-", rownames(sobj), value=TRUE)
#   if (length(mt_genes) > 0) {
#     sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, features = mt_genes)
#   } else {
#     # snRNA often lacks MT- genes; avoid NA/NaN flows
#     sobj[["percent.mt"]] <- 0
#     log(sprintf("No MT- genes found for %s; setting percent.mt=0.", sobj$sample_uid[1]))
#   }
  
#   # 2. Calculate percent.ribo (MOVED INSIDE)
#   # Use ignore.case=TRUE for robustness (e.g., Rps, rpl)
#   ribo_genes <- grep("^RP[SL][0-9]", rownames(sobj), value = TRUE, ignore.case = TRUE)
#   if (length(ribo_genes) > 0) {
#     sobj[["percent.ribo"]] <- PercentageFeatureSet(sobj, features = ribo_genes)
#   } else {
#     sobj[["percent.ribo"]] <- 0
#     log(sprintf("No Ribo genes (RP[SL]) found for %s; setting percent.ribo=0.", sobj$sample_uid[1]))
#   }
  
#   # 3. Return the object
#   sobj
# }

# # snRNA rule-of-thumb: strict mito; adaptive high-end gates for features/counts
# apply_adaptive_snrna_filters <- function(sobj, mad_k=3, min_cells_after_filter=20) {
#   md <- sobj@meta.data
#   mito_median <- median(md$percent.mt, na.rm=TRUE)
#   # keep strict cap at 1%; if MT all zeros, this remains 0–1%
#   mito_final  <- 1.0

#   feat_median <- median(md$nFeature_RNA, na.rm=TRUE)
#   feat_mad    <- mad(md$nFeature_RNA, na.rm=TRUE, constant=1.4826)
#   feat_lower  <- 200
#   feat_upper  <- feat_median + mad_k * feat_mad

#   cnt_median  <- median(md$nCount_RNA, na.rm=TRUE)
#   cnt_mad     <- mad(md$nCount_RNA, na.rm=TRUE, constant=1.4826)
#   cnt_upper   <- cnt_median + mad_k * cnt_mad

#   passing <- which(
#     md$nFeature_RNA > feat_lower & md$nFeature_RNA < feat_upper &
#     md$nCount_RNA   < cnt_upper  &
#     md$percent.mt   < mito_final
#   )

#   if (length(passing) < min_cells_after_filter) {
#     log(sprintf("WARNING: Filtering would leave %d nuclei in %s; skipping filtering.",
#                 length(passing), sobj$sample_uid[1]))
#     return(sobj)
#   }

#   log(sprintf("QC %s: nFeature_RNA > %d & < %.0f | nCount_RNA < %.0f | percent.mt < %.2f (median mt=%.2f)",
#               sobj$sample_uid[1], feat_lower, feat_upper, cnt_upper, mito_final, mito_median))
#   subset(sobj, cells = colnames(sobj)[passing])
# }

# save_qc_violin <- function(sobj, outdir, suid, stage=c("pre_qc","post_qc")){
#   pdir <- file.path(outdir, "plots", suid)
#   dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
#   for (feat in c("nFeature_RNA","nCount_RNA","percent.mt", "percent.ribo")) {
#     p <- VlnPlot(sobj, features=feat, pt.size=0) + NoLegend() + ggtitle(paste0(suid," • ",feat," (",stage,")"))
#     ggsave(file.path(pdir, paste0(stage,"_",feat,".png")), plot=p, width=6, height=5, dpi=100)
#   }
# }

# run_denoise <- function(sobj){
#   sce <- as.SingleCellExperiment(sobj)
#   tbl <- scDblFinder::scDblFinder(sce, returnType="table")
#   sobj$doublet_class <- tbl[colnames(sobj), "class"]
#   subset(sobj, subset = doublet_class == "singlet")
# }

# sct_normalize <- function(sobj){
#   suppressWarnings({
#     sobj <- SCTransform(sobj, vst.flavor="v2", verbose=FALSE, vars.to.regress=c("percent.mt", "percent.ribo"))
#     sobj <- RunPCA(sobj, verbose=FALSE)
#     sobj <- RunUMAP(sobj, dims=1:30, verbose=FALSE)
#   })
#   sobj
# }

# # ---------- Manifest & ID handling ----------
# manifest <- data.table::fread(manifest_path, sep="\t", header=TRUE, data.table=FALSE)
# req <- c("project_id","format","path_or_stem","sample_id")
# miss <- setdiff(req, names(manifest))
# if (length(miss) > 0) stop(sprintf("Manifest missing required columns: %s", paste(miss, collapse=", ")))

# sanitize <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)
# extract_leaf <- function(path) {
#   m <- regexpr("SAMN\\d+", path, perl=TRUE)
#   if (m[1] != -1) regmatches(path, m)[1] else basename(dirname(path))
# }

# # Prefer provided sample_uid; otherwise build & make.unique
# if (!("sample_uid" %in% names(manifest))) {
#   leaf <- vapply(manifest$path_or_stem, extract_leaf, character(1))
#   built <- sanitize(paste(manifest$project_id, manifest$sample_id, leaf, sep="__"))
#   manifest$sample_uid <- make.unique(built, sep="_v")
# } else {
#   # Fill empties, sanitize, ensure uniqueness
#   need <- which(is.na(manifest$sample_uid) | manifest$sample_uid == "")
#   if (length(need) > 0) {
#     leaf <- vapply(manifest$path_or_stem[need], extract_leaf, character(1))
#     built <- sanitize(paste(manifest$project_id[need], manifest$sample_id[need], leaf, sep="__"))
#     manifest$sample_uid[need] <- make.unique(built, sep="_v")
#   }
#   manifest$sample_uid <- make.unique(sanitize(manifest$sample_uid), sep="_v")
# }

# # Ensure per-project RDS folders
# for (p in unique(manifest$project_id)) {
#   dir.create(file.path(outdir, "rds", sanitize(p)), recursive=TRUE, showWarnings=FALSE)
# }

# # Write summary header if missing
# if (!file.exists(summary_tsv)) {
#   write.table(
#     data.frame(project_id=character(), sample_uid=character(), orig_sample_id=character(),
#                genome=character(), chemistry=character(),
#                cells_raw=integer(), cells_postQC=integer(),
#                median_genes=double(), median_counts=double(), median_pct_mt=double(),
#                median_pct_ribo=double(), # <-- ADD THIS COLUMN
#                path_rds=character(), stringsAsFactors=FALSE),
#     file=summary_tsv, sep="\t", quote=FALSE, row.names=FALSE
#   )
# }

# # ---------- Main loop ----------
# for (i in seq_len(nrow(manifest))) {
#   row <- manifest[i,]
#   project <- as.character(row$project_id)
#   fmt     <- as.character(row$format)
#   src     <- as.character(row$path_or_stem)
#   sample  <- as.character(row$sample_id)
#   suid    <- as.character(row$sample_uid %||% sample)
#   genome  <- if ("genome"    %in% names(row)) as.character(row$genome)    else NA_character_
#   chem    <- if ("chemistry" %in% names(row)) as.character(row$chemistry) else NA_character_

#   project_s <- sanitize(project)
#   rds_path  <- file.path(outdir, "rds", project_s, paste0(suid, ".rds"))

#   if (file.exists(rds_path)) { log(sprintf("Skip %s (exists: %s)", suid, rds_path)); next }

#   dir.create(dirname(rds_path), recursive=TRUE, showWarnings=FALSE)

#   tryCatch({
#     log(sprintf("Processing [%s] fmt=%s uid=%s (orig=%s)", project, fmt, suid, sample))

#     counts <- switch(fmt,
#       "10x_h5"   = read_10x_h5(src),
#       "mtx_dir"  = read_mtx_dir(src),
#       "mtx_stem" = read_mtx_stem(src),
#       "tsv_wide" = read_tsv_wide_sparse(src),
#       # backward compatibility if an old manifest said "mtx"
#       "mtx"      = read_mtx_stem(src),
#       stop(sprintf("Unknown format: %s", fmt))
#     )
#     if (is.list(counts) && !is.data.frame(counts)) {
#       # 10x Read10X can return a list; prefer Gene Expression if present else first
#       counts <- if ("Gene Expression" %in% names(counts)) counts[["Gene Expression"]] else counts[[1]]
#     }

#     sobj <- CreateSeuratObject(counts, project=project_s)
#     sobj$project_id     <- project
#     sobj$orig_sample_id <- sample
#     sobj$sample_uid     <- suid
#     if (!is.na(genome)) sobj$genome <- genome
#     if (!is.na(chem))   sobj$chemistry <- chem

#     sobj <- add_qc_metrics(sobj)
#     raw_n <- ncol(sobj)
#     save_qc_violin(sobj, outdir, suid, stage="pre_qc")

#     sobj <- apply_adaptive_snrna_filters(sobj)

#     if (ncol(sobj) > 50) {
#       sobj <- run_denoise(sobj)
#     } else {
#       log(sprintf("Skip doublet removal: too few nuclei in %s (%d)", suid, ncol(sobj)))
#     }

#     postqc_n <- ncol(sobj)
#     if (postqc_n > 50) {
#       save_qc_violin(sobj, outdir, suid, stage="post_qc")
#       sobj <- sct_normalize(sobj)
#     } else {
#       log(sprintf("Skip normalization: too few nuclei in %s (%d)", suid, postqc_n))
#     }

#     log(sprintf("Saving %s | raw=%d postQC=%d", suid, raw_n, postqc_n))
#     saveRDS(sobj, rds_path, compress="xz")

#     summary_df <- data.frame(
#       project_id   = project,
#       sample_uid   = suid,
#       orig_sample_id = sample,
#       genome       = genome,
#       chem         = chem,
#       cells_raw    = raw_n,
#       cells_postQC = postqc_n,
#       median_genes = median(sobj$nFeature_RNA, na.rm=TRUE),
#       median_counts= median(sobj$nCount_RNA, na.rm=TRUE),
#       median_pct_mt= median(sobj$percent.mt, na.rm=TRUE),
#       median_pct_ribo= median(sobj$percent.ribo, na.rm=TRUE), # <-- ADD THIS VALUE
#       path_rds     = rds_path,
#       stringsAsFactors=FALSE
#     )
#     write.table(summary_df, file=summary_tsv, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)

#   }, error=function(e){
#     elog(sprintf("Failed uid %s (src=%s): %s", suid, src, e$message))
#   })
# }

# log(sprintf("All processing complete. Summary at %s", summary_tsv))
# RS

# log "snRNA ingest+QC (ADAPTIVE) finished."

#!/usr/bin/env bash
set -Euo pipefail
IFS=$'\n\t'

# ======================= CONFIG ==========================
BASE="/mnt/18T/chibao/gliomas/data/upstream/snRNA/official_dry_run"
OUT="$BASE"
MANIFEST="$OUT/snrna_official_manifest.tsv"   # manifest for snRNA only
# =========================================================

mkdir -p "$OUT"/{rds,logs,summary,plots}

FAIL_LOG="$OUT/logs/failed_samples.log"
: > "$OUT/logs/error.log"
: > "$FAIL_LOG"

log(){ printf "[%s] %s\n" "$(date '+%Y-%m-%d %H:%M:%S%z')" "$*" | tee -a "$OUT/logs/main.log" ; }
err(){ printf "[%s] ERROR: %s\n" "$(date '+%Y-%m-%d %H:%M:%S%z')" "$*" | tee -a "$OUT/logs/error.log" >&2 ; }

if [[ ! -f "$MANIFEST" ]]; then
  err "Manifest not found: $MANIFEST"
  exit 1
fi

if ! command -v Rscript >/dev/null 2>&1; then
  err "Rscript not found in PATH"
  exit 1
fi

log "snRNA ingest+QC (DRY-RUN SUBSET) started. MANIFEST=$MANIFEST OUT=$OUT"

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
summary_tsv <- file.path(outdir, "summary", "snrna_cohort_summary.tsv")
dir.create(dirname(summary_tsv), recursive=TRUE, showWarnings=FALSE)

log  <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"), paste0(..., collapse=" ")), file=log_main,  append=TRUE)
elog <- function(...) cat(sprintf("[%s] ERROR: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S%z"), paste0(..., collapse=" ")), file=log_error, append=TRUE)

# Write raw manifest row so you can regenerate a retry manifest later
log_failure <- function(row) {
  write.table(row, file = fail_log_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
}

# ---------- Readers ----------
read_10x_h5 <- function(h5_path) {
  if (!file.exists(h5_path)) {
    stop(sprintf("10x_h5 file not found: %s", h5_path))
  }
  Read10X_h5(h5_path)
}

read_mtx_dir <- function(dirpath) {
  if (!dir.exists(dirpath)) {
    stop(sprintf("mtx_dir path does not exist: %s", dirpath))
  }
  Read10X(data.dir = dirpath)
}

read_rds_object <- function(rds_path) {
  if (!file.exists(rds_path)) stop(sprintf("RDS file not found: %s", rds_path))
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
  if (!dir.exists(dirname(stem))) {
    stop(sprintf("Directory for stem does not exist: %s", dirname(stem)))
  }

  # helper to find files that may or may not be gzipped
  find_file <- function(pattern) {
    path <- list.files(dirname(stem), pattern = pattern, full.names = TRUE)
    if (length(path) > 0) path[1] else NA_character_
  }

  bn <- basename(stem)
  bfile <- find_file(paste0("(", bn, "_)?barcodes\\.tsv(\\.gz)?$"))
  ffile <- find_file(paste0("(", bn, "_)?features\\.tsv(\\.gz)?$|(", bn, "_)?genes\\.tsv(\\.gz)?$"))
  mfile <- find_file(paste0("(", bn, "_)?matrix\\.mtx(\\.gz)?$"))

  if (any(is.na(c(bfile, ffile, mfile))))
    stop(sprintf("Could not find all MTX files for stem: %s (bfile=%s, ffile=%s, mfile=%s)", stem, bfile, ffile, mfile))

  feat_dt  <- data.table::fread(ffile, header = FALSE, nrows = 3, data.table = FALSE)
  feat_cols <- ncol(feat_dt)
  feat_col  <- if (feat_cols >= 2) 2 else 1
  message(sprintf("read_mtx_stem(): %s -> features columns=%d, using column %d", ffile, feat_cols, feat_col))

  counts <- ReadMtx(
    mtx = mfile, cells = bfile, features = ffile,
    feature.column = feat_col, cell.column = 1, unique.features = FALSE
  )

  full_feat <- data.table::fread(ffile, header = FALSE, data.table = FALSE, select = feat_col)
  gene_names <- make.unique(as.character(full_feat[[1]]))
  gene_names <- sub("\r$", "", gene_names)

  rownames(counts) <- gene_names
  counts
}

read_tsv_wide_sparse <- function(path) {
  if (!file.exists(path)) stop(sprintf("TSV file not found: %s", path))
  dt <- data.table::fread(path, sep="\t", header=TRUE, data.table=FALSE)
  if (ncol(dt) < 2) stop(sprintf("TSV file %s has <2 columns; expecting gene + at least 1 cell column.", path))
  valid_rows <- !is.na(dt[[1]]) & dt[[1]] != ""; dt <- dt[valid_rows, ]
  rn <- dt[[1]]; dt[[1]] <- NULL
  if (any(duplicated(rn))) { log("Found duplicate rownames in TSV; making them unique."); rn <- make.unique(rn) }
  mat <- as.matrix(dt); rownames(mat) <- rn
  Matrix::Matrix(mat, sparse = TRUE)
}

# ---------- QC / Filtering for snRNA ----------
add_qc_metrics_snrna <- function(sobj){
  # Mitochondrial percentage (mostly near 0% in nuclei, for diagnostics only)
  mt_genes <- grep("^MT-", rownames(sobj), value=TRUE)
  if (length(mt_genes) > 0) {
    sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, features = mt_genes)
  } else {
    log("No mitochondrial genes (starting with 'MT-') found. Setting percent.mt to 0.")
    sobj[["percent.mt"]] <- 0
  }
  sobj
}

apply_adaptive_snrna_filters <- function(sobj, mad_k=3, min_cells_after_filter=20) {
  meta <- sobj@meta.data

  # Basic sanity checks
  if (!all(c("nFeature_RNA","nCount_RNA") %in% colnames(meta))) {
    stop("Metadata missing nFeature_RNA or nCount_RNA after CreateSeuratObject.")
  }

  feat_median <- median(meta$nFeature_RNA, na.rm = TRUE)
  feat_mad    <- mad(meta$nFeature_RNA, na.rm = TRUE, constant=1.4826)
  feat_upper  <- feat_median + mad_k * feat_mad

  count_median <- median(meta$nCount_RNA, na.rm = TRUE)
  count_mad    <- mad(meta$nCount_RNA, na.rm = TRUE, constant=1.4826)
  count_upper  <- count_median + mad_k * count_mad

  # For snRNA: do NOT filter on percent.mt (usually close to 0)
  passing_cells <- which(
    meta$nFeature_RNA > 200 &
    meta$nFeature_RNA < feat_upper &
    meta$nCount_RNA  < count_upper
  )

  if (length(passing_cells) < min_cells_after_filter) {
    log(sprintf("WARNING (snRNA): Adaptive filtering would leave %d cells. Skipping filtering for sample %s.",
                length(passing_cells), sobj$sample_uid[1]))
    return(sobj)
  }

  log(sprintf("snRNA QC thresholds for %s: nFeature_RNA > 200 & < %.0f | nCount_RNA < %.0f (no percent.mt filter applied)",
              sobj$sample_uid[1], feat_upper, count_upper))

  subset(sobj, cells = colnames(sobj)[passing_cells])
}

save_qc_violin <- function(sobj, outdir, sample_uid, stage = c("pre_qc","post_qc")){
  pdir <- file.path(outdir, "plots", sample_uid)
  dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
  feats <- c("nFeature_RNA","nCount_RNA","percent.mt")
  for (feat in feats) {
    if (!feat %in% colnames(sobj@meta.data)) next
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

sct_normalize_snrna <- function(sobj){
  suppressWarnings({
    # For snRNA: often no need to regress percent.mt; keep it simple & stable
    sobj <- SCTransform(sobj, vst.flavor="v2", verbose=FALSE, vars.to.regress = NULL)
    sobj <- RunPCA(sobj, npcs = 40, verbose=FALSE)
    sobj <- RunUMAP(sobj, dims = 1:40, verbose=FALSE)
  })
  sobj
}

# ---------- Manifest ingest & UID handling ----------
manifest <- data.table::fread(manifest_path, sep="\t", header=TRUE, data.table=FALSE)

required <- c("project_id","format","path_or_stem","sample_id")
missing  <- setdiff(required, names(manifest))
if (length(missing) > 0) stop(sprintf("Manifest missing required columns: %s", paste(missing, collapse=", ")))

# For this dry-run, restrict to the two snRNA projects you named
allowed_projects <- c("PRJNA1081384","PRJNA1155637")
manifest <- manifest[manifest$project_id %in% allowed_projects, , drop=FALSE]

if (nrow(manifest) == 0) {
  stop("After filtering for PRJNA1081384 and PRJNA1155637, manifest has 0 rows.")
}

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
  need <- which(is.na(manifest$sample_uid) | manifest$sample_uid == "")
  if (length(need) > 0) {
    leaf  <- vapply(manifest$path_or_stem[need], extract_leaf, character(1))
    base  <- paste(manifest$project_id[need], manifest$sample_id[need], leaf, sep="__")
    base  <- sanitize(base)
    manifest$sample_uid[need] <- make.unique(base, sep="_v")
  }
  manifest$sample_uid <- sanitize(manifest$sample_uid)
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

total_samples <- nrow(manifest)
log(sprintf("Total snRNA samples to process (after filter): %d", total_samples))

# ---------- Processing loop ----------
for (i in seq_len(total_samples)) {
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

  progress_str <- sprintf("Sample %d/%d (%.1f%%)", i, total_samples, 100 * i / total_samples)

  if (file.exists(rds_path)) {
    log(sprintf("%s: Skipping %s (project %s): RDS exists: %s", progress_str, suid, project, rds_path))
    next
  }

  dir.create(dirname(rds_path), recursive = TRUE, showWarnings = FALSE)

  tryCatch({
    total_steps <- 6
    step <- 1

    log(sprintf("%s: [Step %d/%d] Starting processing for format=%s sample_uid=%s (orig=%s)",
                progress_str, step, total_steps, fmt, suid, sample))

    # ---- Step 1: read counts ----
    log(sprintf("%s: [Step %d/%d] Reading counts from: %s",
                progress_str, step, total_steps, src))

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

    if (!inherits(counts_matrix, "dgCMatrix")) {
      counts_matrix <- as(counts_matrix, "dgCMatrix")
    }

    if (nrow(counts_matrix) == 0 || ncol(counts_matrix) == 0) {
      stop(sprintf("Counts matrix for %s is empty: %d genes x %d cells", suid, nrow(counts_matrix), ncol(counts_matrix)))
    }

    step <- step + 1
    log(sprintf("%s: [Step %d/%d] Creating Seurat object and adding snRNA QC metrics",
                progress_str, step, total_steps))

    sobj <- CreateSeuratObject(counts_matrix, project = project_s)
    sobj$project_id     <- project
    sobj$orig_sample_id <- sample
    sobj$sample_uid     <- suid
    if (!is.na(genome)) sobj$genome <- genome
    if (!is.na(chem))   sobj$chemistry <- chem

    sobj <- add_qc_metrics_snrna(sobj)
    raw_n <- ncol(sobj)
    save_qc_violin(sobj, outdir, suid, stage = "pre_qc")

    step <- step + 1
    log(sprintf("%s: [Step %d/%d] Applying snRNA adaptive QC filters",
                progress_str, step, total_steps))

    sobj <- apply_adaptive_snrna_filters(sobj)

    step <- step + 1
    if (ncol(sobj) > 50) {
      log(sprintf("%s: [Step %d/%d] Running doublet detection (scDblFinder) for snRNA",
                  progress_str, step, total_steps))
      sobj <- run_denoise(sobj)
    } else {
      log(sprintf("%s: [Step %d/%d] Skipping doublet detection (only %d cells after QC)",
                  progress_str, step, total_steps, ncol(sobj)))
    }

    postqc_n <- ncol(sobj)

    step <- step + 1
    if (postqc_n > 50) {
      save_qc_violin(sobj, outdir, suid, stage = "post_qc")
      log(sprintf("%s: [Step %d/%d] Running SCTransform + PCA(40) + UMAP(40 dims) for snRNA",
                  progress_str, step, total_steps))
      sobj <- sct_normalize_snrna(sobj)
    } else {
      log(sprintf("%s: [Step %d/%d] Skipping SCTransform/UMAP (only %d cells post-QC)",
                  progress_str, step, total_steps, postqc_n))
    }

    step <- step + 1
    log(sprintf("%s: [Step %d/%d] Saving RDS and updating summary",
                progress_str, step, total_steps))

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

    log(sprintf("%s: Completed sample %s. Raw cells: %d, Post-QC cells: %d",
                progress_str, suid, raw_n, postqc_n))

  }, error=function(e){
    msg <- sprintf("Failed on sample_uid %s (src: %s): %s", suid, src, e$message)
    elog(msg)
    log_failure(row)
  })
}

log(sprintf("All snRNA processing complete. Summary table at %s", summary_tsv))
RS

# --- Post-run: build a retry manifest if any failures occurred ---
RETRY_MANIFEST="$OUT/snrna_manifest_retry.tsv"
if [[ -f "$FAIL_LOG" && -s "$FAIL_LOG" ]]; then
  log "Errors were detected. Generating a retry manifest..."
  head -n 1 "$MANIFEST" > "$RETRY_MANIFEST"
  cat "$FAIL_LOG" >> "$RETRY_MANIFEST"
  log "Retry manifest created at: $RETRY_MANIFEST"
else
  log "Run completed with no errors."
  rm -f "$RETRY_MANIFEST"
fi

log "snRNA ingest+QC (DRY-RUN SUBSET) finished."
