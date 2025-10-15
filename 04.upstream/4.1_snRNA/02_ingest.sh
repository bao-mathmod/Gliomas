#!/usr/bin/env bash
set -Euo pipefail
IFS=$'\n\t'

BASE="/mnt/18T/chibao/gliomas/data/upstream/snRNA"
OUT="/mnt/18T/chibao/gliomas/data/upstream/snRNA/set2"
MANIFEST="/mnt/18T/chibao/gliomas/data/upstream/snRNA/set2/snrna_manifest_2.tsv"

mkdir -p "$OUT"/{rds,logs,summary,tmp}

log(){ printf "[%s] %s\n" "$(date '+%Y-%m-%d %H:%M:%S%z')" "$*" | tee -a "$OUT/logs/main.log" ; }
err(){ printf "[%s] ERROR: %s\n" "$(date '+%Y-%m-%d %H:%M:%S%z')" "$*" | tee -a "$OUT/logs/error.log" >&2 ; }

log "snRNA ingest+QC started. MANIFEST=$MANIFEST OUT=$OUT"

# NOTE: arguments ($MANIFEST, $OUT) MUST appear *before* the here-doc marker.
Rscript --vanilla - "$MANIFEST" "$OUT" <<'RS' || { err "Rscript failed"; exit 1; }
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(stringr)
  library(dplyr)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(glmGamPoi)
  library(ggplot2)  # for plots
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

safe_fread <- function(p) data.table::fread(p, sep="\t", header=TRUE, data.table=FALSE)

# --- readers ---
read_10x_h5 <- function(h5_path){
  Read10X_h5(h5_path)
}
read_10x_mtx_dir <- function(dirpath){
  Read10X(data.dir = dirpath)
}
read_10x_mtx_stem <- function(stem){
  # Resolve actual paths for the three triplet files (gz or plain)
  pick <- function(sans_ext_gz, sans_ext_plain){
    gz <- paste0(stem, sans_ext_gz)
    pl <- paste0(stem, sans_ext_plain)
    if (file.exists(gz)) return(gz)
    if (file.exists(pl)) return(pl)
    return(NA_character_)
  }

  bfile <- pick("barcodes.tsv.gz", "barcodes.tsv")
  ffile <- pick("features.tsv.gz", "features.tsv")
  mfile <- pick("matrix.mtx.gz",  "matrix.mtx")

  missing <- setNames(c(is.na(bfile), is.na(ffile), is.na(mfile)),
                      c("barcodes","features","matrix"))
  if (any(missing)) {
    stop(sprintf(
      "Triplet missing for stem=%s :: need %s",
      stem,
      paste(names(missing)[missing], collapse=", ")
    ))
  }

  # Read with data.table/readMM
  # barcodes
  barcodes <- data.table::fread(bfile, header=FALSE, sep="\t", data.table=FALSE)[,1]
  # features (accept 2 or 3 column 10x ‘features.tsv’)
  feats_dt <- data.table::fread(ffile, header=FALSE, sep="\t", data.table=FALSE)
  if (ncol(feats_dt) >= 2) {
    genes <- feats_dt[[2]]
  } else {
    genes <- feats_dt[[1]]
  }
  # matrix
  M <- Matrix::readMM(mfile)
  M <- as(M, "dgCMatrix")  # ensure sparse and silence dgT->dgC warnings

  # Sanity checks
  if (nrow(M) != length(genes)) stop("Row mismatch matrix vs features")
  if (ncol(M) != length(barcodes)) stop("Col mismatch matrix vs barcodes")

  rownames(M) <- make.unique(genes)
  colnames(M) <- barcodes
  return(M)  # Return a matrix compatible with CreateSeuratObject()
}

# --- reader for wide TSV matrices (single all-samples table) ---
read_tsv_wide_sparse <- function(path) {
  dt <- data.table::fread(path, sep="\t", header=TRUE, data.table=FALSE)
  rn <- dt[[1]]; dt[[1]] <- NULL
  mat <- as.matrix(dt)
  rownames(mat) <- rn
  Matrix::Matrix(mat, sparse = TRUE)
}

# --- nuclei QC helpers ---
add_qc_metrics <- function(sobj){
  mt_genes <- grep("^MT-", rownames(sobj), value=TRUE, ignore.case=FALSE)
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, features = mt_genes)

  # ribosomal kept as DIAGNOSTIC only (no hard filtering)
  rrna_genes <- grep("^RP[SL]", rownames(sobj), value=TRUE)
  sobj[["percent.ribo"]] <- PercentageFeatureSet(sobj, features = rrna_genes)
  sobj
}

# Flag (do not filter) ribo outliers: stored in metadata for diagnostics
flag_ribo_outliers <- function(sobj, k = 3){
  rb <- sobj$percent.ribo
  rb_med <- stats::median(rb, na.rm=TRUE)
  rb_mad <- stats::mad(rb, center = rb_med, constant = 1, na.rm=TRUE)
  sobj$ribo_high_flag <- rb > (rb_med + k * rb_mad)
  sobj$ribo_med <- rb_med
  sobj$ribo_mad <- rb_mad
  sobj
}

# CONSTANT thresholds (your original), with mito tightened to 5% and no ribo in filter
apply_snrna_filters <- function(sobj,
                                min_features=400, max_features=7000,
                                min_counts=800, max_mt=5){
  # Add !is.na() checks to all numeric filters for safety
  subset(sobj, subset = !is.na(nFeature_RNA) & nFeature_RNA >= min_features &
                   !is.na(nFeature_RNA) & nFeature_RNA <= max_features &
                   !is.na(nCount_RNA)   & nCount_RNA   >= min_counts &
                   !is.na(percent.mt)   & percent.mt   <= max_mt )
}

# Save violin plots for QC metrics at a given stage ("pre_qc" / "post_qc")
save_qc_violin <- function(sobj, outdir, sample_id, stage = c("pre_qc","post_qc")){
  stage <- match.arg(stage)
  pdir <- file.path(outdir, "plots", sample_id, stage)
  dir.create(pdir, recursive = TRUE, showWarnings = FALSE)

  feats <- c("nFeature_RNA","nCount_RNA","percent.mt","percent.ribo")
  meta <- sobj@meta.data

  for (feat in feats) {
    if (!feat %in% colnames(meta)) next
    p <- VlnPlot(sobj, features = feat, pt.size = 0.1) +
         NoLegend() +
         ggtitle(paste0(sample_id, " • ", feat, " (", stage, ")"))
    ggplot2::ggsave(filename = file.path(pdir, paste0(feat, ".png")),
                    plot = p, width = 7, height = 5, dpi = 150)
  }

  if ("ribo_high_flag" %in% colnames(meta)) {
    df <- data.frame(flag = meta$ribo_high_flag)
    p2 <- ggplot(df, aes(x = flag)) + geom_bar() +
          ggtitle(paste0(sample_id, " • ribo_high_flag (", stage, ")")) +
          xlab("ribo_high_flag") + ylab("count")
    ggplot2::ggsave(filename = file.path(pdir, "ribo_high_flag.png"),
                    plot = p2, width = 6, height = 4, dpi = 150)
  }
}

run_denoise <- function(sobj){
  sce <- as.SingleCellExperiment(sobj)
  dbl <- scDblFinder::scDblFinder(sce, samples = sobj$sample_id)
  sobj$doublet_score <- dbl$scDblFinder.score
  sobj$doublet_call  <- dbl$scDblFinder.class
  # Add !is.na() check before comparing the string
  subset(sobj, subset = !is.na(doublet_call) & doublet_call == "singlet")
}
sct_normalize <- function(sobj){
  suppressWarnings({
    sobj <- SCTransform(sobj, vst.flavor="v2", verbose=FALSE, vars.to.regress = "percent.ribo")
    sobj <- RunPCA(sobj, verbose=FALSE)
    sobj <- RunUMAP(sobj, dims = 1:30, verbose=FALSE)
  })
  sobj
}

manifest <- data.table::fread(manifest_path, sep="\t", header=TRUE, data.table=FALSE)
if (!file.exists(summary_tsv)) {
  write.table(data.frame(project_id=character(), sample_id=character(),
                         cells_raw=integer(), cells_postQC=integer(),
                         mean_genes=double(), pct_mt_mean=double(), pct_ribo_mean=double(),
                         ribo_outlier_rate=double(),
                         path_rds=character(), stringsAsFactors = FALSE),
              file=summary_tsv, sep="\t", quote=FALSE, row.names=FALSE)
}

for (i in seq_len(nrow(manifest))) {
  row <- manifest[i,]
  project <- row$project_id
  fmt     <- row$format
  src     <- row$path_or_stem
  sample  <- row$sample_id
  tryCatch({
    log(sprintf("Processing [%s] format=%s src=%s sample_id=%s", project, fmt, src, sample))

    if (fmt == "10x_h5") {
      m <- read_10x_h5(src)
      if (is.list(m)) {
        m <- if (!is.null(m[["Gene Expression"]])) m[["Gene Expression"]] else m[[1]]
      }
      sobj <- CreateSeuratObject(m, assay="RNA", project=project)
      sobj$sample_id <- sample

    } else if (fmt == "mtx_dir") {
      m <- read_10x_mtx_dir(src)
      sobj <- CreateSeuratObject(m, assay="RNA", project=project)
      sobj$sample_id <- sample

    } else if (fmt == "mtx") {
      m <- read_10x_mtx_stem(src)
      sobj <- CreateSeuratObject(m, assay="RNA", project=project)
      sobj$sample_id <- sample

    } else if (fmt == "tsv_wide") {
      sp <- read_tsv_wide_sparse(src)
      all_cells <- colnames(sp)
      gsm_prefix <- stringr::str_extract(all_cells, "GSM\\d+")
      have_gsm <- !is.na(gsm_prefix)
      sample_keys <- if (any(have_gsm)) gsm_prefix else stringr::str_replace(all_cells, "_.*$", "")
      for (sk in unique(sample_keys)) {
        cell_idx <- which(sample_keys == sk)
        sub <- sp[, cell_idx, drop=FALSE]
        sid <- if (!is.na(sample) && nzchar(sample)) sample else paste0(project, "__", sk)
        sobj <- CreateSeuratObject(sub, assay="RNA", project=project)
        sobj$sample_id <- sid

        # --- PRE-QC ---
        sobj <- add_qc_metrics(sobj); raw_n <- ncol(sobj)
        sobj <- flag_ribo_outliers(sobj, k = 3)
        save_qc_violin(sobj, outdir, sobj$sample_id[1], stage = "pre_qc")

        # --- Hard filter (constants), no ribo in filter, tighter mito ---
        sobj <- apply_snrna_filters(sobj)

        # --- Doublet removal ---
        sobj <- run_denoise(sobj)

        # --- POST-QC ---
        save_qc_violin(sobj, outdir, sobj$sample_id[1], stage = "post_qc")

        # --- Normalize + reduce ---
        sobj <- sct_normalize(sobj)
        colnames(sobj) <- paste0(colnames(sobj), "-", sobj$sample_id[1])

        rds_path <- file.path(outdir, "rds", paste0(sobj$sample_id[1], ".rds"))
        saveRDS(sobj, rds_path, compress="xz")

        df <- data.frame(
          project_id = project,
          sample_id  = sobj$sample_id[1],
          cells_raw  = raw_n,
          cells_postQC = ncol(sobj),
          mean_genes = mean(sobj$nFeature_RNA),
          pct_mt_mean = mean(sobj$percent.mt),
          pct_ribo_mean = mean(sobj$percent.ribo),
          ribo_outlier_rate = mean(sobj$ribo_high_flag),
          path_rds = rds_path,
          stringsAsFactors = FALSE
        )
        write.table(df, file=summary_tsv, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
        log(sprintf("Saved %s (postQC=%d)", rds_path, ncol(sobj)))
      }
      next

    } else {
      stop(sprintf("Unknown format: %s", fmt))
    }

    # --- PRE-QC ---
    sobj <- add_qc_metrics(sobj); raw_n <- ncol(sobj)
    sobj <- flag_ribo_outliers(sobj, k = 3)
    save_qc_violin(sobj, outdir, sobj$sample_id[1], stage = "pre_qc")

    # --- Hard filter (constants), no ribo in filter, tighter mito ---
    sobj <- apply_snrna_filters(sobj)

    # --- Doublet removal ---
    sobj <- run_denoise(sobj)

    # --- POST-QC ---
    save_qc_violin(sobj, outdir, sobj$sample_id[1], stage = "post_qc")

    # --- Normalize + reduce ---
    sobj <- sct_normalize(sobj)
    colnames(sobj) <- paste0(colnames(sobj), "-", sobj$sample_id[1])

    rds_path <- file.path(outdir, "rds", paste0(sobj$sample_id[1], ".rds"))
    saveRDS(sobj, rds_path, compress="xz")

    df <- data.frame(
      project_id = project,
      sample_id  = sobj$sample_id[1],
      cells_raw  = raw_n,
      cells_postQC = ncol(sobj),
      mean_genes = mean(sobj$nFeature_RNA),
      pct_mt_mean = mean(sobj$percent.mt),
      pct_ribo_mean = mean(sobj$percent.ribo),
      ribo_outlier_rate = mean(sobj$ribo_high_flag),
      path_rds = rds_path,
      stringsAsFactors = FALSE
    )
    write.table(df, file=summary_tsv, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
    log(sprintf("Saved %s (postQC=%d)", rds_path, ncol(sobj)))

  }, error=function(e){
    msg <- sprintf("Failed: project=%s format=%s src=%s sample=%s :: %s", project, fmt, src, sample, e$message)
    elog(msg)
  })
}

log(sprintf("All done. Summary at %s", summary_tsv))
RS

log "snRNA ingest+QC finished."
