suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(copykat)      # after installation, see below
  library(ggplot2)
})

# Config
# Path where CopyKAT writes its outputs
copykat_outdir <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/copykat/output"

dir.create(copykat_outdir, recursive = TRUE, showWarnings = FALSE)

# Helper: check if CopyKAT has finished for a sample
is_copykat_done <- function(sample_id, out_dir = copykat_outdir) {
  # Use the prediction file as the "flag" that the run completed
  pred_file <- file.path(
    out_dir,
    paste0(sample_id, "_copykat_prediction.txt")
  )
  file.exists(pred_file)
}


# Load object
adult_obj <- readRDS(
  "/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/harmony_cleaned_annotated_v2.rds"
)

adult_obj

# Check assays
Assays(adult_obj)
DefaultAssay(adult_obj)

adult_obj[["RNA"]]
adult_obj[["SCT"]]

# Use raw counts from RNA assay
DefaultAssay(adult_obj) <- "RNA"
rna_counts <- LayerData(adult_obj, assay = "RNA", layer = "counts")
dim(rna_counts)
head(rownames(rna_counts))
head(colnames(rna_counts))

# Run per sample
table(adult_obj$sample_uid)  # or other column you used

# Normal Reference
table(adult_obj$general_cell_type)
normal_types <- c("TILs", "Myeloid")

norm_cells_global <- colnames(adult_obj)[adult_obj$general_cell_type %in% normal_types]
length(norm_cells_global)

# Run CopyKat
run_copykat_for_sample <- function(
  obj,
  sample_id,
  sample_col = "sample_uid",
  assay = "RNA",
  genome = "hg20",             # hg20 in CopyKAT = human hg19 style annotation
  id.type = "S",               # "S" = SYMBOL; "E" = ENSEMBL
  norm_cells_global = NULL,    # vector of colnames(obj) to treat as normal
  n.cores = 8,
  out_dir = copykat_outdir
) {
  message("Running CopyKAT for sample: ", sample_id)

  # Subset to that sample
  obj_s <- subset(obj, subset = (!!as.name(sample_col)) == sample_id)
  message("Cells in sample ", sample_id, ": ", ncol(obj_s))

  # Extract raw counts
  counts_s <- GetAssayData(obj_s, assay = assay, slot = "counts")
  message("Counts matrix dims: ", paste(dim(counts_s), collapse = " x "))

  # Optionally select normal cells for baseline
  if (!is.null(norm_cells_global)) {
    norm_cells_s <- intersect(colnames(counts_s), norm_cells_global)
  } else {
    norm_cells_s <- character(0)  # no known normals
  }

  # Convert to matrix (CopyKAT expects a standard matrix)
  rawmat <- as.matrix(counts_s)

  # Ensure output dir exists and temporarily setwd there
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(out_dir)

  ck_res <- copykat(
    rawmat = rawmat,
    id.type = id.type,
    ngene.chr = 5,
    win.size = 25,
    KS.cut = 0.1,
    sam.name = sample_id,
    distance = "euclidean",
    norm.cell.names = norm_cells_s,
    output.seg = "FALSE",
    plot.genes = TRUE,
    genome = genome,
    n.cores = n.cores
  )

  # Optionally also save the R object
  saveRDS(
    ck_res,
    file = file.path(out_dir, paste0(sample_id, "_copykat_clustering_results.rds"))
  )

  return(ck_res)
}


# Run CopyKat on all samples and collect results
sample_ids <- sort(unique(adult_obj$sample_uid))

# Filter to only unfinished samples
to_run <- sample_ids[!vapply(sample_ids, is_copykat_done, logical(1))]

message("Total samples: ", length(sample_ids))
message("Already completed: ", length(sample_ids) - length(to_run))
message("To (re)run: ", length(to_run))

copykat_results <- list()
library(future)
plan("multisession", workers = 10)

for (sid in to_run) {
  message("==== Processing sample: ", sid, " ====")
  copykat_results[[sid]] <- run_copykat_for_sample(
    obj = adult_obj,
    sample_id = sid,
    sample_col = "sample_uid",
    genome = "hg20",
    id.type = "S",
    norm_cells_global = norm_cells_global,
    n.cores = 8,
    out_dir = copykat_outdir
  )
}

# Map CopyKat results back to Seurat object
pred_list <- lapply(names(copykat_results), function(sid) {
  pred <- copykat_results[[sid]]$prediction
  pred$sample_uid <- sid
  pred
})

pred_df <- do.call(rbind, pred_list)
head(pred_df)
table(pred_df$copykat.pred)

# Align by cell names and add to metadata
# Make sure rownames = cell names
rownames(pred_df) <- pred_df$cell.names

# Create a vector aligned to adult_obj cells (NA if not in a prediction set)
copykat_pred_vec <- pred_df[colnames(adult_obj), "copykat.pred"]

# Add to metadata
adult_obj$copykat.pred <- copykat_pred_vec

# Optionally: a simpler binary status "tumor" vs "normal"
adult_obj$copykat.status <- dplyr::case_when(
  adult_obj$copykat.pred == "aneuploid" ~ "tumor_like",
  adult_obj$copykat.pred == "diploid"  ~ "normal_like",
  TRUE                                  ~ "undefined"
)
table(adult_obj$copykat.status, useNA = "ifany")

# Visualize 
DefaultAssay(adult_obj) <- "SCT"

DimPlot(
  adult_obj,
  reduction = "umap.harmony",
  group.by = "copykat.status",
  label = TRUE
) + ggtitle("CopyKAT CNV status on Harmony UMAP")
