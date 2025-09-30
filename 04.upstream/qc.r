#!/usr/bin/env Rscript
# Step 2: Per-cell QC & ambient RNA correction (SoupX), with per-sample thresholds
# Outputs:
#   - 02_qc_filtered.rds
#   - qc_summary.tsv
#   - QC_plots.pdf

suppressPackageStartupMessages({
  library(Seurat)         # v5 recommended
  library(Matrix)
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(cli)
})

# ------------------------------ CONFIG ---------------------------------------
RAW_RDS           <- "/mnt/12T/chibao/data/stuff_data/special/seurat_objects/01_raw_seurat.rds"
SAMPLESHEET       <- "/mnt/12T/chibao/data/stuff_data/special/sample_metadata.tsv"
OUTDIR            <- "/mnt/12T/chibao/data/stuff_data/special/seurat_objects"
PLOT_PDF          <- file.path(OUTDIR, "QC_plots.pdf")
OUT_RDS           <- file.path(OUTDIR, "02_qc_filtered.rds")
QC_SUMMARY_TSV    <- file.path(OUTDIR, "qc_summary.tsv")

# If TRUE, attempt SoupX ambient correction per Run_ID from Cell Ranger outs
RUN_SOUPX         <- TRUE

# Species heuristics for MT/ribo/Hb gene patterns (human default)
# You can flip to mouse patterns automatically if Reference_Genome suggests mm10.
HUMAN_PATTERNS <- list(
  mt   = "^MT-",
  ribo = "^RPL|^RPS",
  hbb  = "^HB[AB][A-Z]*"   # HBA*, HBB*
)
MOUSE_PATTERNS <- list(
  mt   = "^mt-",
  ribo = "^Rpl|^Rps",
  hbb  = "^Hb[ab][a-z]*"
)

# Thresholding strategy parameters (per-sample)
# - minGenes: infer via knee + protect by quantile floor
# - maxMT: set data-driven by upper whisker (Q3 + 1.5*IQR) capped to 30%; never below 5%
# - maxCounts: cap extreme library sizes at Q3 + 4*IQR (rarely binding)
MIN_GENES_QUANTILE_FLOOR <- 0.05
MAX_MT_CAP               <- 30   # %
MAX_MT_FLOOR             <- 5    # %
MAX_COUNT_IQR_MULT       <- 4

# Optional: Doublet detection (disabled by default; uncomment later if needed)
RUN_DOUBLETS <- FALSE
# -----------------------------------------------------------------------------

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

message(">> Loading inputs")
so_all <- readRDS(RAW_RDS)
meta   <- read.delim(SAMPLESHEET, check.names = FALSE, sep = "\t")

# Basic checks
stopifnot(all(c("Sample_ID","Run_ID","Reference_Genome","CellRanger_Output_Dir") %in% colnames(so_all@meta.data)))
stopifnot(all(c("Run_ID","CellRanger_Output_Dir","Reference_Genome") %in% colnames(meta)))

# ------------------------- helpers -------------------------------------------
species_patterns_for <- function(ref) {
  ref_l <- tolower(ref %||% "")
  # crude heuristic
  if (grepl("mm10|mouse|murine|mm39", ref_l)) return(MOUSE_PATTERNS)
  return(HUMAN_PATTERNS)
}

# Safe coalesce
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Compute knee (barcode-rank) gene cutoff:
# Sort cells by nFeature_RNA desc; find elbow via a simple second-derivative heuristic.
# Returns suggested minGenes threshold.
suggest_min_genes_knee <- function(nfeat_vec) {
  v <- sort(as.numeric(nfeat_vec), decreasing = TRUE)
  if (length(v) < 100) return(round(quantile(v, 0.05)))   # fallback for tiny samples
  x <- seq_along(v)
  # smooth with rolling median to reduce noise
  k <- pmin(51, floor(length(v)/20) * 2 + 1) # odd window ~5%
  v_s <- stats::runmed(v, k = k)
  # approximate elbow as index of max curvature (discrete second derivative)
  d1 <- c(0, diff(v_s))
  d2 <- c(0, diff(d1))
  elbow_idx <- which.min(d2) # most negative curvature
  cutoff <- max(10, round(v_s[elbow_idx] * 0.5))  # conservative (half of elbow height)
  return(cutoff)
}

# Boxplot-upper whisker rule (Q3 + 1.5*IQR)
upper_whisker <- function(x) {
  q <- quantile(x, probs = c(.25,.75), na.rm = TRUE)
  iqr <- q[2] - q[1]
  as.numeric(q[2] + 1.5*iqr)
}

# ---------------------- SoupX ambient correction -----------------------------
# We correct per Run_ID (most stable) using each run's Cell Ranger outs.
# Requires SoupX; install if needed: remotes::install_github("constantAmateur/SoupX")
# Helper: read a 10x matrix from a directory with matrix.mtx(.gz)/features.tsv(.gz)/barcodes.tsv(.gz)
# Returns a dgCMatrix. If multiple feature types exist, prefer "Gene Expression".
.read10x_mat <- function(dir) {
  m <- Seurat::Read10X(dir)
  if (is.list(m)) {
    if ("Gene Expression" %in% names(m)) {
      m <- m[["Gene Expression"]]
    } else {
      m <- m[[1]]
    }
  }
  m
}

# Preferred: use SoupX::load10X on the *outs* directory (it auto-detects raw/filtered)
# Fallback: build a SoupChannel from matrices read via Seurat::Read10X
run_soupx_on_run <- function(run_dir, fallback_contamination = 0.05) {
  if (!dir.exists(run_dir)) stop("Cell Ranger outs dir not found: ", run_dir)

  raw_dir  <- file.path(run_dir, "raw_feature_bc_matrix")
  filt_dir <- file.path(run_dir, "filtered_feature_bc_matrix")
  if (!dir.exists(filt_dir)) stop("filtered_feature_bc_matrix missing under: ", run_dir)

  # Helper: read dgCMatrix from a 10x directory
  .read10x_mat <- function(dir) {
    m <- Seurat::Read10X(dir)
    if (is.list(m)) {
      if ("Gene Expression" %in% names(m)) m <- m[["Gene Expression"]] else m <- m[[1]]
    }
    m
  }

  # 1) Try the robust "one-liner" path
  sc <- tryCatch({
    SoupX::load10X(run_dir, keepDroplets = TRUE)  # uses raw+filtered if available
  }, error = function(e) {
    message("   [SoupX] load10X() failed: ", conditionMessage(e))
    NULL
  })

  # 2) Manual construction if load10X failed (or raw missing)
  if (is.null(sc)) {
    if (dir.exists(raw_dir)) {
      tod <- .read10x_mat(raw_dir)   # raw/droplets
      toc <- .read10x_mat(filt_dir)  # filtered/cells
    } else {
      message("   [SoupX] raw_feature_bc_matrix missing; using filtered for both tod/toc (less ideal).")
      tod <- .read10x_mat(filt_dir); toc <- tod
    }
    sc <- SoupX::SoupChannel(tod = tod, toc = toc)
  }

  # 3) Try autoEstCont with a ladder of increasingly permissive params
  tries <- list(
    list(tfidfMin = 1.0, soupQuantile = 0.90),
    list(tfidfMin = 0.5, soupQuantile = 0.80),
    list(tfidfMin = 0.25, soupQuantile = 0.70),
    list(tfidfMin = 0.10, soupQuantile = 0.50)
  )

  est_ok <- FALSE
  last_err <- NULL
  for (p in tries) {
    message(sprintf(
      "   [SoupX] autoEstCont() try: tfidfMin=%.2f, soupQuantile=%.2f",
      p$tfidfMin, p$soupQuantile
    ))
    sc_try <- tryCatch({
      SoupX::autoEstCont(sc, doPlot = FALSE,
                         tfidfMin = p$tfidfMin,
                         soupQuantile = p$soupQuantile)
    }, error = function(e) e)

    if (inherits(sc_try, "error")) {
      last_err <- sc_try
      message("   [SoupX] autoEstCont() failed: ", conditionMessage(sc_try))
    } else {
      sc <- sc_try
      est_ok <- TRUE
      break
    }
  }

    # 4) Adjust counts (either using estimated rho, or a conservative fixed rho)
  adj <- NULL
  if (est_ok) {
    message("   [SoupX] autoEstCont() succeeded, adjusting counts with estimated contamination.")
    adj <- SoupX::adjustCounts(sc, roundToInt = TRUE)
  } else {
    # NEW: explicitly set a fallback rho before adjustCounts()
    message("   [SoupX][WARN] autoEstCont() failed for this run even after relax. ",
            "Falling back to fixed contamination = ", fallback_contamination, 
            ". Setting rho via setContaminationFraction() and proceeding.")
    sc  <- SoupX::setContaminationFraction(sc, contFrac = fallback_contamination, forceAccept = TRUE)
    adj <- SoupX::adjustCounts(sc, roundToInt = TRUE)

  }
}


# Apply SoupX per Run_ID and rebuild a corrected Seurat for each Run_ID
soupx_correct_per_run <- function(so_subset) {
  run_id <- unique(so_subset$Run_ID)
  stopifnot(length(run_id) == 1)
  cr_dir <- unique(so_subset$CellRanger_Output_Dir)
  stopifnot(length(cr_dir) == 1)

  message("   [SoupX] Run_ID=", run_id, " | dir=", cr_dir)
  adj <- run_soupx_on_run(cr_dir)

  # Subset adjusted matrix to the barcodes present in this object.
  # Our cell barcodes were renamed in Step 1 to "<Sample_ID>__<Run_ID>__<barcode>".
  # We need to map back: take raw barcode tail after the last "__".
  old_bc  <- colnames(so_subset)
  raw_bc  <- sub(".*__", "", old_bc)  # original 10x barcode
  keep_bc <- intersect(colnames(adj), raw_bc)
  if (length(keep_bc) == 0) {
    stop("No overlapping barcodes between SoupX matrix and Seurat colnames for run: ", run_id)
  }
  adj_sub <- adj[, keep_bc, drop = FALSE]

  # Put back the prefixed cell names order
  adj_sub <- adj_sub[, raw_bc, drop = FALSE]  # align columns to raw_bc order (NA will drop)
  colnames(adj_sub) <- old_bc                 # restore prefixed names

  # Recreate Seurat object with corrected counts; keep metadata
  so_corr <- CreateSeuratObject(counts = adj_sub,
                                project = unique(so_subset$Project_ID),
                                min.cells = 0, min.features = 0)
  so_corr@meta.data <- so_subset@meta.data[match(colnames(so_corr), rownames(so_subset@meta.data)), , drop = FALSE]
  so_corr$CellID <- colnames(so_corr)
  return(so_corr)
}

if (RUN_SOUPX) {
  message(">> Running SoupX ambient RNA correction per Run_ID ...")
  runs <- split(colnames(so_all), so_all$Run_ID)
  so_runs <- lapply(names(runs), function(rid) {
    so_r <- subset(so_all, cells = runs[[rid]])
    soupx_correct_per_run(so_r)
  })
  # Merge corrected runs back together
  so_all <- Reduce(function(a,b) merge(a, b), so_runs)
  DefaultAssay(so_all) <- "RNA"
}

# -------------------------- QC metrics ---------------------------------------
message(">> Computing QC metrics")

# choose gene patterns by sample (if mixed species, we compute per-cell using its Reference_Genome)
meta_ref <- so_all$Reference_Genome
patterns_per_cell <- vapply(meta_ref, function(ref) {
  ref_l <- tolower(ref %||% "")
  if (grepl("mm10|mouse|murine|mm39", ref_l)) "mouse" else "human"
}, character(1))

gene_names <- rownames(so_all)
is_human <- patterns_per_cell == "human"

# Precompute feature sets (vectorized once for speed)
mt_human   <- grepl(HUMAN_PATTERNS$mt, gene_names, perl = TRUE)
ribo_human <- grepl(HUMAN_PATTERNS$ribo, gene_names, perl = TRUE)
hbb_human  <- grepl(HUMAN_PATTERNS$hbb, gene_names, perl = TRUE)

mt_mouse   <- grepl(MOUSE_PATTERNS$mt, gene_names, perl = TRUE)
ribo_mouse <- grepl(MOUSE_PATTERNS$ribo, gene_names, perl = TRUE)
hbb_mouse  <- grepl(MOUSE_PATTERNS$hbb, gene_names, perl = TRUE)

# Compute percentages per cell by switching index sets per genome
# calc_pct <- function(ind) Matrix::colSums(GetAssayData(so_all, slot = "counts")[ind, , drop = FALSE]) /
#                             Matrix::colSums(GetAssayData(so_all, slot = "counts")) * 100

calc_pct <- function(ind) Matrix::colSums(so_all[["RNA"]]$counts[ind, , drop = FALSE]) /
                            Matrix::colSums(so_all[["RNA"]]$counts) * 100

pct_mt   <- numeric(ncol(so_all))
pct_ribo <- numeric(ncol(so_all))
pct_hbb  <- numeric(ncol(so_all))

pct_mt[is_human]    <- calc_pct(mt_human)[is_human]
pct_ribo[is_human]  <- calc_pct(ribo_human)[is_human]
pct_hbb[is_human]   <- calc_pct(hbb_human)[is_human]
pct_mt[!is_human]   <- calc_pct(mt_mouse)[!is_human]
pct_ribo[!is_human] <- calc_pct(ribo_mouse)[!is_human]
pct_hbb[!is_human]  <- calc_pct(hbb_mouse)[!is_human]

so_all$percent.mt    <- pct_mt
so_all$percent.ribo  <- pct_ribo
so_all$percent.hbb   <- pct_hbb

# Library size & features
# so_all$nCount_RNA  <- Matrix::colSums(GetAssayData(so_all, slot = "counts"))
# so_all$nFeature_RNA<- Matrix::colSums(GetAssayData(so_all, slot = "counts") > 0)
# Library size & features
so_all$nCount_RNA   <- Matrix::colSums(so_all[["RNA"]]$counts)
so_all$nFeature_RNA <- Matrix::colSums(so_all[["RNA"]]$counts > 0)

# -------------------- Per-sample thresholds (data-driven) ---------------------
library(pkgconfig)
message(">> Deriving per-sample QC thresholds")

samples <- unique(so_all$Sample_ID)
qc_rows <- list()

# one big PDF for plots
pdf(PLOT_PDF, width = 8.5, height = 6)

for (sid in samples) {
  so_s <- subset(so_all, subset = Sample_ID == sid)

  # Violin plots per Run_ID, per-sample summary
  p1 <- VlnPlot(so_s, features = c("nCount_RNA","nFeature_RNA","percent.mt","percent.ribo","percent.hbb"),
                group.by = "Run_ID", pt.size = 0, ncol = 3) + ggtitle(paste0("QC violins: ", sid))
  print(p1)

  # Barcode rank "knee" curve (per Run_ID shown together)
  df_rank <- so_s@meta.data %>%
    dplyr::select(nFeature_RNA, Run_ID) %>%
    mutate(nFeature_RNA = as.numeric(nFeature_RNA)) %>%
    group_by(Run_ID) %>%
    arrange(desc(nFeature_RNA), .by_group = TRUE) %>%
    mutate(rank = row_number()) %>%
    ungroup()
  p2 <- ggplot(df_rank, aes(x = rank, y = nFeature_RNA, color = Run_ID)) +
    geom_line() +
    scale_y_log10() + scale_x_log10() +
    labs(title = paste0("Barcode-rank (knee) per Run_ID: ", sid),
         x = "Rank (log10)", y = "Detected genes per cell (log10)") +
    theme_bw()
  print(p2)

  # Suggest thresholds
  # 1) minGenes from knee + quantile floor
  knee_cut <- suggest_min_genes_knee(so_s$nFeature_RNA)
  q_floor  <- as.numeric(quantile(so_s$nFeature_RNA, MIN_GENES_QUANTILE_FLOOR, na.rm = TRUE))
  minGenes <- max(100, min(knee_cut, max(q_floor, 100)))  # never below 100; conservative

  # 2) maxMT from upper whisker, floored at 5% and capped at 30%
  uw_mt    <- upper_whisker(so_s$percent.mt)
  maxMT    <- min(MAX_MT_CAP, max(MAX_MT_FLOOR, round(uw_mt, 1)))

  # 3) maxCounts to remove extreme libraries
  uw_cnt   <- upper_whisker(so_s$nCount_RNA)
  iqr_cnt  <- IQR(so_s$nCount_RNA, na.rm = TRUE)
  q3_cnt   <- quantile(so_s$nCount_RNA, .75, na.rm = TRUE)
  maxCounts<- as.numeric(q3_cnt + MAX_COUNT_IQR_MULT * iqr_cnt)

  # Show thresholds on plots
  p3 <- ggplot(so_s@meta.data, aes(x = nFeature_RNA)) +
    geom_histogram(bins = 60) + geom_vline(xintercept = minGenes, color = "red") +
    ggtitle(paste0(sid, " | nFeature_RNA (red=minGenes=", minGenes, ")"))
  print(p3)

  p4 <- ggplot(so_s@meta.data, aes(x = percent.mt)) +
    geom_histogram(bins = 60) + geom_vline(xintercept = maxMT, color = "red") +
    ggtitle(paste0(sid, " | percent.mt (red=maxMT=", maxMT, "%)"))
  print(p4)

  p5 <- ggplot(so_s@meta.data, aes(x = nCount_RNA)) +
    geom_histogram(bins = 60) + geom_vline(xintercept = maxCounts, color = "red") +
    ggtitle(paste0(sid, " | nCount_RNA (red=maxCounts≈", round(maxCounts), ")"))
  print(p5)

  # Record thresholds
  qc_rows[[sid]] <- data.frame(
    Sample_ID = sid,
    threshold_min_genes = minGenes,
    threshold_max_percent_mt = maxMT,
    threshold_max_counts = round(maxCounts),
    stringsAsFactors = FALSE
  )
}
dev.off()

qc_tbl <- bind_rows(qc_rows)
# -------------------- Apply per-sample QC thresholds --------------------------
message(">> Applying per-sample QC filters")

keep_cells <- logical(ncol(so_all))
names(keep_cells) <- colnames(so_all)

for (i in seq_len(nrow(qc_tbl))) {
  sid <- qc_tbl$Sample_ID[i]
  minG <- qc_tbl$threshold_min_genes[i]
  maxMT<- qc_tbl$threshold_max_percent_mt[i]
  maxC <- qc_tbl$threshold_max_counts[i]

  idx <- which(so_all$Sample_ID == sid)
  pass <- (so_all$nFeature_RNA[idx] >= minG) &
          (so_all$percent.mt[idx] <= maxMT) &
          (so_all$nCount_RNA[idx] <= maxC)
  keep_cells[idx] <- pass
}

retained <- names(keep_cells)[keep_cells]
so_filt  <- subset(so_all, cells = retained)

# Summaries per Sample_ID
summ <- so_filt@meta.data %>%
  group_by(Sample_ID) %>%
  summarise(retained_cells = n(),
            median_nFeature = median(nFeature_RNA),
            median_nCount   = median(nCount_RNA),
            median_mt       = median(percent.mt)) %>%
  right_join(qc_tbl, by = "Sample_ID") %>%
  arrange(Sample_ID)

write_tsv(summ, QC_SUMMARY_TSV)

message(">> Saving filtered object: ", OUT_RDS)
saveRDS(so_filt, OUT_RDS)

message(">> Done. Retained cells: ", ncol(so_filt), " / ", ncol(so_all))
