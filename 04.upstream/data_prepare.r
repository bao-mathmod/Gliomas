#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(stringr)
})

# ------------------------- CONFIG ---------------------------------
SAMPLESHEET <- "/mnt/12T/chibao/data/stuff_data/special/sample_metadata.tsv"   # from your bash builder
OUTDIR       <- "/mnt/12T/chibao/data/stuff_data/special/seurat_objects"
SAVE_BY_SAMPLE <- TRUE     # if TRUE: also save one RDS per Sample_ID (merged across runs)
CELLRANGER_VERSION_EXPECTED <- "9.0.1"   # optional sanity check; set "" to skip
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
# ------------------------------------------------------------------

message(">> Loading samplesheet: ", SAMPLESHEET)
meta <- read.delim(SAMPLESHEET, check.names = FALSE, sep = "\t")

required_cols <- c(
  "Project_ID","Sample_ID","Experiment_ID","Run_ID",
  "CellRanger_Output_Dir","CellRanger_Version","Reference_Genome"
)
missing_req <- setdiff(required_cols, colnames(meta))
if (length(missing_req) > 0) {
  stop("Missing required column(s) in samplesheet: ", paste(missing_req, collapse = ", "))
}

# Optional: warn if CR version differs from what you expect
if (nzchar(CELLRANGER_VERSION_EXPECTED)) {
  bad_ver <- meta %>% filter(CellRanger_Version != CELLRANGER_VERSION_EXPECTED)
  if (nrow(bad_ver) > 0) {
    warning("Some rows have CellRanger_Version != ", CELLRANGER_VERSION_EXPECTED,
            "\nRows (first 5):\n",
            paste(utils::capture.output(print(head(bad_ver[,c('Run_ID','CellRanger_Version')]))), collapse="\n"))
  }
}

# Utility: make a safe prefix for cell barcodes
safe_id <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)   # keep simple chars
  x
}

# Read a single Cell Ranger "outs/filtered_feature_bc_matrix" and return a Seurat object
read_one_run <- function(row) {
  cr_dir <- row[["CellRanger_Output_Dir"]]
  if (!dir.exists(cr_dir)) {
    stop("Directory not found: ", cr_dir)
  }

  # Cell Ranger layout:
  # - filtered_feature_bc_matrix/ (matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz)
  # - metrics_summary.csv (optional sanity check)
  matrix_dir <- file.path(cr_dir, "filtered_feature_bc_matrix")
  if (!dir.exists(matrix_dir)) {
    stop("filtered_feature_bc_matrix not found under: ", cr_dir)
  }

  # Read counts
  counts <- Read10X(matrix_dir)
  # Create a Seurat object
  so <- CreateSeuratObject(counts = counts,
                           project = as.character(row[["Project_ID"]]),
                           min.cells = 0, min.features = 0)

  # Make unique barcodes: <Sample_ID>__<Run_ID>__<barcode>
  pref <- paste0(safe_id(row[["Sample_ID"]]), "__", safe_id(row[["Run_ID"]]), "__")
  colnames(so) <- paste0(pref, colnames(so))
  so$CellID <- colnames(so)  # add explicit cell id column

  # Attach *all* samplesheet columns to each cell
  for (cn in colnames(meta)) {
    so[[cn]] <- as.character(row[[cn]])
  }

  # Optional: pull a few QC numbers from metrics_summary.csv (if present)
  mfile <- file.path(cr_dir, "metrics_summary.csv")
  if (file.exists(mfile)) {
    # a very tolerant reader (doesn't assume order)
    ms <- tryCatch({
      suppressWarnings(readr::read_csv(mfile, show_col_types = FALSE, n_max = 2))
    }, error = function(e) NULL)

    if (!is.null(ms) && nrow(ms) >= 1) {
      # standardize names
      nm <- tolower(gsub("[^a-z0-9]+","_", names(ms)))
      names(ms) <- nm

      # a few useful fields if present
      pick <- c("estimated_number_of_cells","mean_reads_per_cell","median_genes_per_cell",
                "number_of_reads","valid_barcodes","sequencing_saturation",
                "fraction_reads_in_cells")
      for (p in pick) {
        if (p %in% names(ms)) {
          val <- ms[[p]][1]
          # strip % and commas
          val <- gsub(",", "", as.character(val))
          val <- gsub("%$", "", val)
          so[[paste0("ms_", p)]] <- val
        }
      }
    }
  }

  so
}

# ---- Read all runs ------------------------------------------------
message(">> Reading all runs listed in samplesheet...")
objs <- vector("list", nrow(meta))
for (i in seq_len(nrow(meta))) {
  row <- meta[i, , drop = FALSE]
  msg <- paste0("[", i, "/", nrow(meta), "] ", row$Run_ID, " (", row$CellRanger_Output_Dir, ")")
  message("   - ", msg)
  objs[[i]] <- read_one_run(row)
}

# ---- Merge runs belonging to the same biological Sample_ID --------
# Rationale: a patient/sample may have been sequenced in multiple SRR runs.
# We keep Run_ID in metadata, but at the cell level we often want a single object per Sample_ID.
message(">> Merging runs per Sample_ID ...")
by_sample <- split(objs, sapply(objs, function(x) x$Sample_ID[1]))

merge_list <- lapply(names(by_sample), function(sid) {
  lst <- by_sample[[sid]]
  if (length(lst) == 1) {
    return(lst[[1]])
  } else {
    # Reduce merge: keep all genes (union). Seurat handles sparse matrices efficiently.
    m <- Reduce(function(a, b) merge(a, y = b, project = a@project.name),
                lst)
    return(m)
  }
})
names(merge_list) <- names(by_sample)

# ---- Save each Sample_ID object (optional) ------------------------
if (SAVE_BY_SAMPLE) {
  message(">> Saving per-sample Seurat objects ...")
  for (sid in names(merge_list)) {
    fn <- file.path(OUTDIR, paste0("sample_", safe_id(sid), "_raw.rds"))
    saveRDS(merge_list[[sid]], fn)
  }
}

# ---- Merge all samples into one project object --------------------
message(">> Merging all samples into one Seurat object ...")
all_obj <- Reduce(function(a, b) merge(a, y = b, project = a@project.name),
                  merge_list)

# basic sanity columns
stopifnot("CellID" %in% colnames(all_obj@meta.data))
stopifnot("Sample_ID" %in% colnames(all_obj@meta.data))
stopifnot("Run_ID" %in% colnames(all_obj@meta.data))

# ---- Write outputs -----------------------------------------------
out_rds   <- file.path(OUTDIR, "01_raw_seurat.rds")
out_meta  <- file.path(OUTDIR, "01_cell_metadata.tsv")
out_cells <- file.path(OUTDIR, "01_cell_ids.txt")

message(">> Saving combined object: ", out_rds)
saveRDS(all_obj, out_rds)

message(">> Exporting per-cell metadata: ", out_meta)
mtd <- all_obj@meta.data %>% tibble::rownames_to_column("CellID_raw_index")
write_tsv(mtd, out_meta)

message(">> Exporting CellIDs only: ", out_cells)
writeLines(colnames(all_obj), con = out_cells)

message(">> Done. Cells: ", ncol(all_obj), " | Genes: ", nrow(all_obj))
