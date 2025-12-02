# 1. install JAGS (if not already, via conda outside R):
# conda install -c conda-forge jags

# 2. Load packages
# suppressPackageStartupMessages({
#   library(Seurat)
#   library(infercnv)
#   library(dplyr)
#   library(ggplot2)
#   library(future)
# })

# options(future.globals.maxSize = 100 * 1024^3)  # 100 GB
# # 3. Load object
# adult_obj <- readRDS(
#   "/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/harmony_cleaned_annotated_v2.rds"
# )

# adult_obj
# Assays(adult_obj)
# DefaultAssay(adult_obj)

# DefaultAssay(adult_obj) <- "RNA"
# rna_counts <- LayerData(adult_obj, assay = "RNA", layer = "counts")
# dim(rna_counts)

# # 4. Decide reference and per-sample running
# # 4.1 Sample IDs
# table(adult_obj$sample_uid)
# sample_ids <- sort(unique(adult_obj$sample_uid))
# sample_ids

# # 4.2 Normal reference cell types
# table(adult_obj$general_cell_type)

# normal_types <- c("Myeloid", "TILs")

# # All reference cells across the whole object:
# reference_cells_global <- colnames(adult_obj)[adult_obj$general_cell_type %in% normal_types]
# length(reference_cells_global)

# # 5. Prepare input
# # 5.1 Create output
# infercnv_root <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/infercnv"
# dir.create(infercnv_root, showWarnings = FALSE, recursive = TRUE)
# # 5.2 Gene annotation file
# suppressPackageStartupMessages(library(data.table))

# gtf_file <- "/mnt/18T/chibao/gliomas/book/ref/gencode.v38.annotation.gtf"
# gene_order_file <- file.path(infercnv_root, "gene_order_hg38_from_gencode.txt")

# library(data.table)

# # Read GTF without comment lines
# gtf <- fread(
#   cmd = paste("grep -v '^#' ", gtf_file),
#   header = FALSE,
#   sep = "\t"
# )

# # Name the columns we care about
# colnames(gtf)[c(1,3,4,5,9)] <- c("chr", "feature", "start", "end", "attr")

# # Keep only gene rows
# gene_gtf <- gtf[feature == "gene"]

# # Extract gene_name from attributes
# gene_gtf[, gene := sub('.*gene_name "([^"]+)".*', '\\1', attr)]

# # Filter to genes present in your Seurat object
# genes_in_exp <- rownames(adult_obj)
# gene_gtf <- gene_gtf[gene %in% genes_in_exp]

# # Deduplicate: one row per gene
# # Strategy: take the FIRST chromosome, and min/max of coordinates
# gene_gtf_unique <- gene_gtf[, .(
#   chr   = chr[1],
#   start = min(start),
#   end   = max(end)
# ), by = gene]

# # Remove unwanted chromosomes
# gene_gtf_unique <- gene_gtf_unique[!chr %in% c("chrM", "chrX", "chrY")]

# # (Optional but recommended) Keep only canonical autosomes
# # gene_gtf_unique <- gene_gtf_unique[grepl("^chr[0-9]+$", chr)]

# # Order by chromosome & start
# gene_gtf_unique <- gene_gtf_unique[order(chr, start)]

# # Final sanity check: no duplicated gene names
# stopifnot(!anyDuplicated(gene_gtf_unique$gene))

# # Write gene order file: gene\tchr\tstart\tend
# fwrite(
#   gene_gtf_unique,
#   file = gene_order_file,
#   sep = "\t",
#   col.names = FALSE,
#   quote = FALSE
# )

# # Double-check after writing
# go2 <- read.table(gene_order_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
# colnames(go2) <- c("gene", "chr", "start", "end")
# sum(duplicated(go2$gene))  # should be 0


# # 6. Prepare infercnv input matrix
# prepare_infercnv_input_for_sample <- function(
#   obj,
#   sample_id,
#   sample_col = "sample_uid",
#   group_col = "general_cell_type",
#   out_root = infercnv_root,
#   min_cells_per_group = 5 # Increased safety margin
# ) {
#   message("Preparing inferCNV input for sample: ", sample_id)

#   # 1. Subset to the specific sample
#   obj_s <- subset(obj, subset = (!!as.name(sample_col)) == sample_id)
  
#   # [SEURAT V5] Ensure layers are joined
#   obj_s <- JoinLayers(obj_s)

#   # 2. IDENTIFY AND REMOVE TINY GROUPS
#   # Check group sizes
#   group_counts <- table(obj_s[[group_col]])
#   message("  Original Group Sizes:")
#   print(group_counts)
  
#   # Find groups to keep (must have >= min_cells_per_group)
#   groups_to_keep <- names(group_counts)[group_counts >= min_cells_per_group]
#   groups_to_drop <- names(group_counts)[group_counts < min_cells_per_group]
  
#   if (length(groups_to_drop) > 0) {
#     message("  !!! DROPPING groups with < ", min_cells_per_group, " cells: ", 
#             paste(groups_to_drop, collapse = ", "))
    
#     # Subset the object to keep only valid groups
#     # We use the cell names that belong to the valid groups
#     cells_to_keep <- colnames(obj_s)[obj_s[[group_col]][,1] %in% groups_to_keep]
#     obj_s <- subset(obj_s, cells = cells_to_keep)
#   }

#   message("  Final Cell Count for Analysis: ", ncol(obj_s))
  
#   if (ncol(obj_s) < 10) {
#     message("  --> Too few cells remaining. Skipping.")
#     return(NULL)
#   }

#   # 3. Extract Counts & SANITIZE NAMES
#   counts_s <- GetAssayData(obj_s, assay = "RNA", layer = "counts")
  
#   # Replace hyphens with dots to prevent mismatch errors
#   clean_names <- gsub("-", ".", colnames(counts_s))
#   colnames(counts_s) <- clean_names
  
#   # 4. Create Output Directory
#   out_dir <- file.path(out_root, sample_id)
#   dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

#   # 5. Write Counts Matrix
#   counts_file <- file.path(out_dir, paste0(sample_id, ".counts.matrix.txt"))
#   write.table(
#     as.matrix(counts_s),
#     file = counts_file,
#     sep = "\t",
#     quote = FALSE,
#     col.names = NA
#   )

#   # 6. Write Annotations
#   annots <- data.frame(
#     cell_id = clean_names,
#     group   = obj_s[[group_col]][,1], # Ensure we get the vector
#     stringsAsFactors = FALSE
#   )
  
#   annots_file <- file.path(out_dir, paste0(sample_id, ".annotations.txt"))
#   write.table(
#     annots,
#     file = annots_file,
#     sep = "\t",
#     quote = FALSE,
#     row.names = FALSE,
#     col.names = FALSE
#   )

#   return(list(counts_file = counts_file, annots_file = annots_file, out_dir = out_dir))
# }

# run_infercnv_for_sample <- function(
#   sample_id,
#   counts_file,
#   annots_file,
#   gene_order_file,
#   normal_types,
#   out_dir,
#   min_ref_cells = 10,
#   min_genes = 50,
#   safe_min_group_size = 5,
#   num_threads = 60  # <--- CHANGED TO 10 FOR STABILITY
# ) {
#   message("Processing sample: ", sample_id)

#   # 1. Read annotations
#   annots <- read.table(annots_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
#   colnames(annots) <- c("cell_id", "group")

#   # 2. Check counts for each potential reference group
#   ref_counts <- table(annots$group[annots$group %in% normal_types])
  
#   message("  Reference group sizes in annotation:")
#   print(ref_counts)

#   # 3. Filter small groups
#   valid_ref_groups <- names(ref_counts)[ref_counts >= safe_min_group_size]
#   total_valid_ref <- sum(ref_counts[valid_ref_groups])
  
#   message("  Valid reference groups (>= ", safe_min_group_size, " cells): ", 
#           paste(valid_ref_groups, collapse = ", "))

#   if (length(valid_ref_groups) == 0) {
#     message("  --> No valid reference groups large enough. Skipping sample.")
#     return(NULL)
#   }

#   if (total_valid_ref < min_ref_cells) {
#     message("  --> Total reference cells (", total_valid_ref, 
#             ") < min_ref_cells (", min_ref_cells, "). Skipping sample.")
#     return(NULL)
#   }

#   # 4. Create Object
#   infercnv_obj <- CreateInfercnvObject(
#     raw_counts_matrix = counts_file,
#     annotations_file  = annots_file,
#     delim             = "\t",
#     gene_order_file   = gene_order_file,
#     ref_group_names   = valid_ref_groups,
#     chr_exclude       = c("chrX", "chrY", "chrM")
#   )

#   expr <- infercnv_obj@expr.data
#   if (is.null(dim(expr)) || any(dim(expr) < 2)) {
#     message("  --> Matrix collapsed to vector or empty. Skipping.")
#     return(NULL)
#   }

#   message("Running inferCNV for sample: ", sample_id)
#   message("Using num_threads: ", num_threads)

#   # 5. Run inferCNV
#   # We use the passed 'num_threads' here.
#   # If this still fails, you can set HMM=FALSE, but try threads=1 first.
#   infercnv_obj <- infercnv::run(
#     infercnv_obj,
#     cutoff            = 0.1,
#     out_dir           = file.path(out_dir, "infercnv_out"),
#     cluster_by_groups = TRUE,
#     denoise           = TRUE,
#     HMM               = TRUE,
#     num_threads       = num_threads 
#   )

#   return(infercnv_obj)
# }

# # Loop over all samples
# infercnv_results <- list()

# # Keep your existing run_infercnv_for_sample (the one with safe_min_group_size)
# # And run the loop:

# for (sid in sample_ids) {
#   message("------------------------------------------------------")
#   message("Processing sample: ", sid)

#   out_dir_sid  <- file.path(infercnv_root, sid)
#   done_file    <- file.path(out_dir_sid, "infercnv_out", "infercnv.observation_groupings.txt")

#   if (file.exists(done_file)) {
#     message("  --> inferCNV output already exists. Skipping.")
#     next
#   }

#   # 1. Prepare Files (Removes 1-cell groups automatically)
#   prep <- prepare_infercnv_input_for_sample(
#     obj = adult_obj,
#     sample_id = sid,
#     sample_col = "sample_uid",
#     group_col = "general_cell_type",
#     out_root = infercnv_root,
#     min_cells_per_group = 5  # Strict filter
#   )

#   # If prep returns NULL (too few cells or all groups dropped), skip
#   if (is.null(prep)) {
#     message("  --> Skipping execution for ", sid)
#     next
#   }

#   # 2. Run inferCNV
#   infercnv_results[[sid]] <- run_infercnv_for_sample(
#     sample_id       = sid,
#     counts_file     = prep$counts_file,
#     annots_file     = prep$annots_file,
#     gene_order_file = gene_order_file,
#     normal_types    = normal_types,
#     out_dir         = prep$out_dir,
#     num_threads     = 60,  # Use 10 threads for stability
#     safe_min_group_size = 5 # Matches the prep filter
#   )
# }

# # Map back to Seurat object 
# # For one 
# adult_obj <- infercnv::add_to_seurat(
#   seurat_obj            = adult_obj,
#   assay_name            = "RNA",
#   infercnv_output_path  = "/mnt/18T/chibao/gliomas/analysis/infercnv_adult/P1/infercnv_out",
#   top_n                 = 10,
#   column_prefix         = "infercnv_P1_"
# )

# # For all
# for (sid in sample_ids) {
#   out_path_sid <- file.path(infercnv_root, sid, "infercnv_out")

#   adult_obj <- infercnv::add_to_seurat(
#     seurat_obj           = adult_obj,
#     assay_name           = "RNA",
#     infercnv_output_path = out_path_sid,
#     top_n                = 10,
#     column_prefix        = paste0("infercnv_", sid, "_")
#   )
# }

# head(colnames(adult_obj@meta.data))

# # Summary file
# # Identify all inferCNV columns
# inf_cols <- grep("^infercnv_", colnames(adult_obj@meta.data), value = TRUE)

# # Simple CNV burden: number of CNV segments per cell
# adult_obj$infercnv_cnv_count <- rowSums(adult_obj@meta.data[, inf_cols, drop = FALSE] != 0)

# # Or a weighted sum (if the columns store magnitude)
# adult_obj$infercnv_cnv_sum <- rowSums(abs(adult_obj@meta.data[, inf_cols, drop = FALSE]))

# summary(adult_obj$infercnv_cnv_count)
# summary(adult_obj$infercnv_cnv_sum)

# # Visualize
# DefaultAssay(adult_obj) <- "SCT"

# p_inf <- DimPlot(
#   adult_obj,
#   reduction = "umap.harmony",
#   group.by = "infercnv_status",
#   label = TRUE
# ) + ggtitle("inferCNV CNV burden on Harmony UMAP")

# p_inf

###############################################################
library(infercnv)
library(Seurat)
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(future)

options(future.globals.maxSize = 460 * 1024^2) # 80 GB RAM limit for parallel processes

# Load the seurat object
os_ec_subset <- readRDS("/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/harmony_cleaned_annotated_v2.rds")

# Set the default assay to "RNA" for raw counts
DefaultAssay(os_ec_subset) <- "RNA"
os_ec_subset <- JoinLayers(os_ec_subset, assay = "RNA")
# Step 1: Prepare Your Input Data

# Assuming you already have your Seurat object 'os_subset'
# Extract gene expression data and metadata
# Extract raw counts from the RNA assay
raw_counts <- LayerData(os_ec_subset, assay = "RNA", layer = "counts")
#
cell_annotations <- os_ec_subset@meta.data$general_cell_type

# Create annotation file (cells as rows, one column for cell type)
cell_anno_df <- data.frame(
  cell_type = cell_annotations,
  row.names = colnames(raw_counts)
)

# Define reference cells (normal cells)
reference_cells <- rownames(cell_anno_df)[cell_anno_df$cell_type %in% c("Myeloid", 'TILs')]

# > reference_cells |> length()
# [1] 4637


# Step 2: Run inferCNV Analysis

# Define constants for inferCNV
# wget https://data.broadinstitute.org/Trinity/CTAT/cnv/hg38_gencode_v27.txt -O /mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/infercnv_2/gene_ordering_hg38_genecode_v27.txt

infercnv_root <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/infercnv"
gene_order_file <- file.path(infercnv_root, "gene_order_hg38_from_gencode.txt")
# gene_order_file <- '/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/infercnv_2/gene_ordering_hg38_genecode_v27.txt'
base_out_dir <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/infercnv_2/"

# Create inferCNV object
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = raw_counts,
  annotations_file = cell_anno_df,
  delim = "\t",
  gene_order_file = gene_order_file,  # Gene positions file
  ref_group_names = c('Myeloid', 'TILs')  # Reference cell type
)

options(scipen = 100)
# Run inferCNV
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = base_out_dir,
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = TRUE,
  num_threads = 10,
  
  # === KEY SETTINGS ===
  #tumor_subcluster_partition_method = "random_trees", # Avoids the SNN graph error
  resume_mode = TRUE,    # Picks up from Step 17 automatically
  
  # === DISABLE PLOTTING (To fix memory crash) ===
  no_plot = TRUE,        # Disables final heatmap
  plot_steps = FALSE,    # Disables intermediate heatmaps
  no_prelim_plot = TRUE, # Disables preliminary heatmaps
  write_expr_matrix = TRUE # Ensures you get the text data since we aren't plotting
)

# Option A
library(infercnv)

# 1. Re-create the initial object (takes 2 seconds)
# You need your raw_counts and annotations from your setup script
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = raw_counts,
  annotations_file = cell_anno_df,
  delim = "\t",
  gene_order_file = gene_order_file,
  ref_group_names = c('Myeloid', 'TILs') 
)

# 2. Run the command with Plotting DISABLED
# infercnv will scan 'base_out_dir', see Step 17 exists, load it, 
# and automatically resume at Step 18 without plotting.
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = base_out_dir,
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = TRUE,
  num_threads = 10,
  
  # You can keep 'random_trees' here, but since Step 17 (Leiden) is already 
  # on disk, it will likely continue using the Leiden clusters it already found.
  tumor_subcluster_partition_method = "random_trees", 
  
  # === CRITICAL: STOP THE CRASH ===
  no_plot = TRUE,       # Do not try to draw heatmaps
  plot_steps = FALSE    # Do not draw intermediate steps
)

# Option B
library(infercnv)
options(scipen = 100) # Prevents scientific notation errors in filenames/indices

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = 0.1,
  out_dir = base_out_dir,
  cluster_by_groups = TRUE,
  denoise = TRUE,
  HMM = TRUE,
  num_threads = 50,
  
  # === KEY SETTINGS ===
  #tumor_subcluster_partition_method = "random_trees", # Avoids the SNN graph error
  resume_mode = TRUE,    # Picks up from Step 17 automatically
  
  # === DISABLE PLOTTING (To fix memory crash) ===
  no_plot = TRUE,        # Disables final heatmap
  plot_steps = FALSE,    # Disables intermediate heatmaps
  no_prelim_plot = TRUE, # Disables preliminary heatmaps
  write_expr_matrix = TRUE # Ensures you get the text data since we aren't plotting
)

# Define constants for inferCNV
gene_order_file <- file.path(infercnv_root, "gene_order_hg38_from_gencode.txt")
base_out_dir <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/infercnv_2/"

# Run 2
# Alternative with different HMM settings
infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff=0.1,
    out_dir= base_out_dir,
    cluster_by_groups=TRUE,
    denoise=TRUE,
    HMM=TRUE,
    HMM_type="i3",  # Try a simpler HMM model (i3 instead of i6)
    num_threads=48
)

saveRDS(infercnv_obj, file = "/mnt/d12b/SC_BONE_RAW/data/out_infercnv_os_ec_hmm_i3/infercnv_obj/infercnv_obj.rds")


####### Official ########
# ==============================================================================
# PART 1: SETUP & LIBRARIES
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(infercnv)
  library(dplyr)
  library(ggplot2)
  library(future)
  library(data.table)
})

# Set massive memory limit for Seurat/Future
options(future.globals.maxSize = 100 * 1024^3)

# Define Paths
seurat_path   <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/harmony_cleaned_annotated_v2.rds"
infercnv_root <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/infercnv_3"
gtf_file      <- "/mnt/18T/chibao/gliomas/book/ref/gencode.v38.annotation.gtf"
final_out_path <- file.path(infercnv_root, "adult_obj_with_cnv_final.rds")

# Create output directory
dir.create(infercnv_root, showWarnings = FALSE, recursive = TRUE)

# Load Object
message("Loading Seurat Object...")
adult_obj <- readRDS(seurat_path)
DefaultAssay(adult_obj) <- "RNA"

# Remove 
# Remove specified cell type annotations from the Seurat object
remove_types <- c("Melanoma_Metas", "Lung_Metas")
message("Removing annotations: ", paste(remove_types, collapse = ", "))

message("Counts before removal:")
print(table(adult_obj$general_cell_type))

adult_obj <- subset(adult_obj, subset = !(general_cell_type %in% remove_types))

# Drop unused factor levels in metadata
# obj@meta.data$general_cell_type <- droplevels(obj@meta.data$general_cell_type)

# Ensure Idents reflect the current metadata column (optional but helpful for downstream code)
if ("general_cell_type" %in% colnames(adult_obj@meta.data)) {
    Idents(adult_obj) <- "general_cell_type"
}

message("Removal complete. Counts after removal:")
print(table(adult_obj$general_cell_type))

# Check cell types
adult_obj$general_cell_type |> unique()

# Save RDS
saveRDS(adult_obj, '/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/harmony_cleaned_annotated_v3.rds')

# Define Parameters
normal_types <- c("Myeloid", "TILs", 'B_cell', 'Stromal/Endothelial')
sample_ids   <- sort(unique(adult_obj$sample_uid))

# ==============================================================================
# PART 2: GENE ORDER FILE GENERATION
# ==============================================================================
message("Generating Gene Order File...")
gene_order_file <- file.path(infercnv_root, "gene_order_hg38_from_gencode.txt")

# Read GTF (Filtering comments)
gtf <- fread(cmd = paste("grep -v '^#' ", gtf_file), header = FALSE, sep = "\t")
colnames(gtf)[c(1,3,4,5,9)] <- c("chr", "feature", "start", "end", "attr")

# Filter for genes present in Seurat object
gene_gtf <- gtf[feature == "gene"]
gene_gtf[, gene := sub('.*gene_name "([^"]+)".*', '\\1', attr)]
gene_gtf <- gene_gtf[gene %in% rownames(adult_obj)]

# Deduplicate (One row per gene)
gene_gtf_unique <- gene_gtf[, .(chr = chr[1], start = min(start), end = max(end)), by = gene]
gene_gtf_unique <- gene_gtf_unique[!chr %in% c("chrM", "chrX", "chrY")] # Exclude sex/mito
gene_gtf_unique <- gene_gtf_unique[order(chr, start)]

# Write File
fwrite(gene_gtf_unique, file = gene_order_file, sep = "\t", col.names = FALSE, quote = FALSE)
message("Gene order file created.")

# ==============================================================================
# PART 3: INFERCNV EXECUTION LOOPS
# ==============================================================================

# --- Helper Function 1: Prepare Input ---
prepare_infercnv_input_for_sample <- function(obj, sample_id, group_col = "general_cell_type", out_root, min_cells = 5) {
  obj_s <- subset(obj, subset = sample_uid == sample_id)
  obj_s <- JoinLayers(obj_s)
  
  # Remove tiny groups
  group_counts <- table(obj_s[[group_col]])
  keep_groups <- names(group_counts)[group_counts >= min_cells]
  
  if (length(keep_groups) == 0) return(NULL)
  
  # Filter object
  cells_keep <- colnames(obj_s)[obj_s[[group_col]][,1] %in% keep_groups]
  obj_s <- subset(obj_s, cells = cells_keep)
  
  if (ncol(obj_s) < 10) return(NULL)
  
  # Sanitize Names (Hyphen to Dot for File System Safety)
  counts_s <- GetAssayData(obj_s, assay = "RNA", layer = "counts")
  clean_names <- gsub("-", ".", colnames(counts_s))
  colnames(counts_s) <- clean_names
  
  # Write files
  out_dir <- file.path(out_root, sample_id)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  counts_file <- file.path(out_dir, paste0(sample_id, ".counts.matrix.txt"))
  write.table(as.matrix(counts_s), file = counts_file, sep = "\t", quote = FALSE, col.names = NA)
  
  annots <- data.frame(cell_id = clean_names, group = obj_s[[group_col]][,1], stringsAsFactors = FALSE)
  annots_file <- file.path(out_dir, paste0(sample_id, ".annotations.txt"))
  write.table(annots, file = annots_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  return(list(counts_file = counts_file, annots_file = annots_file, out_dir = out_dir))
}

# --- Helper Function 2: Run InferCNV ---
run_infercnv_for_sample <- function(sid, prep, gene_file, refs, min_ref = 10, threads = 40) {
  annots <- read.table(prep$annots_file, header = FALSE, sep = "\t")
  colnames(annots) <- c("cell", "group")
  
  # Check Reference Cells
  ref_counts <- table(annots$group[annots$group %in% refs])
  valid_refs <- names(ref_counts)[ref_counts >= 5]
  
  if (sum(ref_counts[valid_refs]) < min_ref) {
    message("  -> Not enough reference cells. Skipping.")
    return(NULL)
  }
  
  # Create & Run
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = prep$counts_file,
    annotations_file  = prep$annots_file,
    delim = "\t",
    gene_order_file   = gene_file,
    ref_group_names   = valid_refs,
    chr_exclude       = c("chrX", "chrY", "chrM")
  )
  
  infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff = 0.1,
    out_dir = file.path(prep$out_dir, "infercnv_out"),
    cluster_by_groups = TRUE,
    denoise = TRUE,
    HMM = TRUE,
    num_threads = threads
  )
  return(infercnv_obj)
}

# --- Execution Loop ---
message("Starting InferCNV Loop...")
for (sid in sample_ids) {
  message("Processing: ", sid)
  
  # Check if done
  if (file.exists(file.path(infercnv_root, sid, "infercnv_out", "infercnv.observation_groupings.txt"))) {
    message("  -> Output exists. Skipping.")
    next
  }
  
  # Run Prep
  prep <- prepare_infercnv_input_for_sample(adult_obj, sid, out_root = infercnv_root)
  if (is.null(prep)) next
  
  # Run InferCNV
  try({
    run_infercnv_for_sample(sid, prep, gene_order_file, normal_types)
  })
}

# ==============================================================================
# PART 4: HYBRID RENAMING (FIXING HYPHEN/DOT MISMATCH)
# ==============================================================================
# This ensures Seurat cell names match the files generated by InferCNV
message("------------------------------------------------")
message("PHASE 2: Harmonizing Cell Names (Hybrid Renaming)")

current_cells <- colnames(adult_obj)
new_names_vec <- setNames(current_cells, current_cells) # Default: No change

count_renamed <- 0

for (sid in sample_ids) {
  # Check the actual output file
  obs_file <- file.path(infercnv_root, sid, "infercnv_out", "infercnv.observation_groupings.txt")
  if (!file.exists(obs_file)) next
  
  # Read header to see if InferCNV used dots
  header <- read.table(obs_file, header = TRUE, nrows = 5, stringsAsFactors = FALSE)
  sample_cell <- if(!is.null(rownames(header))) rownames(header)[1] else header[1,1]
  
  # If file uses dots, we must rename Seurat cells to dots
  if (grepl("\\.", sample_cell) && !grepl("-", sample_cell)) {
    cells_in_sample <- colnames(adult_obj)[adult_obj$sample_uid == sid]
    if (length(cells_in_sample) > 0) {
      new_names_vec[cells_in_sample] <- gsub("-", ".", cells_in_sample)
      count_renamed <- count_renamed + 1
    }
  }
}

message("Renaming samples to match InferCNV output format: ", count_renamed)
adult_obj <- RenameCells(adult_obj, new.names = new_names_vec)

# ==============================================================================
# PART 5: ROBUST METADATA EXTRACTION (PARSING HMM)
# ==============================================================================
message("------------------------------------------------")
message("PHASE 3: Extracting CNV Metrics from HMM")

# Define Robust Extraction Function (Handles Cluster Name Prefixes)
extract_cnv_metrics_v2 <- function(sample_id, root_dir) {
  out_dir <- file.path(root_dir, sample_id, "infercnv_out")
  obs_file <- file.path(out_dir, "infercnv.observation_groupings.txt")
  if (!file.exists(obs_file)) return(NULL)
  
  # 1. Read Cell-to-Cluster Map
  groups_df <- read.table(obs_file, header = TRUE, stringsAsFactors = FALSE)
  cnv_data <- data.frame(
    cell_id = rownames(groups_df), 
    infercnv_cluster = groups_df[,1], 
    stringsAsFactors = FALSE
  )
  
  # 2. Read HMM Predictions
  hmm_files <- list.files(out_dir, pattern = "pred_cnv_regions.dat", full.names = TRUE)
  if (length(hmm_files) > 0) {
    hmm_raw <- tryCatch(read.table(hmm_files[1], header = TRUE, stringsAsFactors = FALSE), error = function(e) NULL)
    
    if (!is.null(hmm_raw)) {
      # Clean Cluster Names (Remove "Annotation." prefix if present)
      target_clusters <- unique(cnv_data$infercnv_cluster)
      hmm_names <- unique(hmm_raw$cell_group_name)
      
      if (length(intersect(hmm_names, target_clusters)) == 0) {
        hmm_raw$clean_name <- sub("^[^.]+\\.", "", hmm_raw$cell_group_name)
      } else {
        hmm_raw$clean_name <- hmm_raw$cell_group_name
      }
      
      # Aggregate Metrics
      metrics <- hmm_raw %>%
        mutate(len_mb = (end - start) / 1e6) %>%
        group_by(clean_name) %>%
        summarise(
          total_cnv_mb  = sum(len_mb[state != 3]),
          total_loss_mb = sum(len_mb[state < 3]),
          total_gain_mb = sum(len_mb[state > 3]),
          has_cnv       = if_else(total_cnv_mb > 0, "Yes", "No")
        )
      
      cnv_data <- left_join(cnv_data, metrics, by = c("infercnv_cluster" = "clean_name"))
    }
  }
  
  # Fill NAs
  if ("total_cnv_mb" %in% colnames(cnv_data)) {
    cnv_data[is.na(cnv_data)] <- 0
    cnv_data$has_cnv[cnv_data$has_cnv == 0] <- "No"
  }
  return(cnv_data)
}

# Run Extraction Loop
results_list <- list()
for (sid in sample_ids) {
  res <- extract_cnv_metrics_v2(sid, infercnv_root)
  if (!is.null(res)) results_list[[sid]] <- res
}

# Merge Results
full_cnv_meta <- dplyr::bind_rows(results_list)
rownames(full_cnv_meta) <- full_cnv_meta$cell_id
full_cnv_meta$cell_id <- NULL

# ==============================================================================
# PART 6: MERGE & SAVE
# ==============================================================================
message("------------------------------------------------")
message("PHASE 4: Merging and Saving")

# Add Metadata (Safe Check)
common_cells <- intersect(rownames(full_cnv_meta), colnames(adult_obj))
message("Cells with CNV data found in Seurat object: ", length(common_cells))

if (length(common_cells) > 0) {
  adult_obj <- AddMetaData(adult_obj, metadata = full_cnv_meta)
  
  # Validation Print
  print(head(adult_obj@meta.data[, c("infercnv_cluster", "total_cnv_mb", "has_cnv")]))
  
  # Save Final Object
  message("Saving final object to: ", final_out_path)
  saveRDS(adult_obj, file = final_out_path)
  
  if (file.exists(final_out_path)) {
    message("SUCCESS. Pipeline complete.")
  } else {
    message("ERROR: File save failed.")
  }
} else {
  message("CRITICAL ERROR: No name overlap found. Metadata not added.")
}

DimPlot(
  adult_obj,
  group.by = "has_cnv",
  reduction = "umap.harmony",
  cols = c("No" = "lightgrey", "Yes" = "red")
) + ggtitle("Cells with Detected CNVs")