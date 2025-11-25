# 1. install JAGS (if not already, via conda outside R):
# conda install -c conda-forge jags

# 2. Load packages
suppressPackageStartupMessages({
  library(Seurat)
  library(infercnv)
  library(dplyr)
  library(ggplot2)
  library(future)
})

options(future.globals.maxSize = 100 * 1024^3)  # 100 GB
# 3. Load object
adult_obj <- readRDS(
  "/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/harmony_cleaned_annotated_v2.rds"
)

adult_obj
Assays(adult_obj)
DefaultAssay(adult_obj)

DefaultAssay(adult_obj) <- "RNA"
rna_counts <- LayerData(adult_obj, assay = "RNA", layer = "counts")
dim(rna_counts)

# 4. Decide reference and per-sample running
# 4.1 Sample IDs
table(adult_obj$sample_uid)
sample_ids <- sort(unique(adult_obj$sample_uid))
sample_ids

# 4.2 Normal reference cell types
table(adult_obj$general_cell_type)

normal_types <- c("Myeloid", "TILs")

# All reference cells across the whole object:
reference_cells_global <- colnames(adult_obj)[adult_obj$general_cell_type %in% normal_types]
length(reference_cells_global)

# 5. Prepare input
# 5.1 Create output
infercnv_root <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/infercnv"
dir.create(infercnv_root, showWarnings = FALSE, recursive = TRUE)
# 5.2 Gene annotation file
suppressPackageStartupMessages(library(data.table))

gtf_file <- "/mnt/18T/chibao/gliomas/book/ref/gencode.v38.annotation.gtf"
gene_order_file <- file.path(infercnv_root, "gene_order_hg38_from_gencode.txt")

library(data.table)

# Read GTF without comment lines
gtf <- fread(
  cmd = paste("grep -v '^#' ", gtf_file),
  header = FALSE,
  sep = "\t"
)

# Name the columns we care about
colnames(gtf)[c(1,3,4,5,9)] <- c("chr", "feature", "start", "end", "attr")

# Keep only gene rows
gene_gtf <- gtf[feature == "gene"]

# Extract gene_name from attributes
gene_gtf[, gene := sub('.*gene_name "([^"]+)".*', '\\1', attr)]

# Filter to genes present in your Seurat object
genes_in_exp <- rownames(adult_obj)
gene_gtf <- gene_gtf[gene %in% genes_in_exp]

# Deduplicate: one row per gene
# Strategy: take the FIRST chromosome, and min/max of coordinates
gene_gtf_unique <- gene_gtf[, .(
  chr   = chr[1],
  start = min(start),
  end   = max(end)
), by = gene]

# Remove unwanted chromosomes
gene_gtf_unique <- gene_gtf_unique[!chr %in% c("chrM", "chrX", "chrY")]

# (Optional but recommended) Keep only canonical autosomes
# gene_gtf_unique <- gene_gtf_unique[grepl("^chr[0-9]+$", chr)]

# Order by chromosome & start
gene_gtf_unique <- gene_gtf_unique[order(chr, start)]

# Final sanity check: no duplicated gene names
stopifnot(!anyDuplicated(gene_gtf_unique$gene))

# Write gene order file: gene\tchr\tstart\tend
fwrite(
  gene_gtf_unique,
  file = gene_order_file,
  sep = "\t",
  col.names = FALSE,
  quote = FALSE
)

# Double-check after writing
go2 <- read.table(gene_order_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(go2) <- c("gene", "chr", "start", "end")
sum(duplicated(go2$gene))  # should be 0


# 6. Prepare infercnv input matrix
prepare_infercnv_input_for_sample <- function(
  obj,
  sample_id,
  sample_col = "sample_uid",
  group_col = "general_cell_type",
  out_root = infercnv_root,
  min_cells_per_group = 5 # Increased safety margin
) {
  message("Preparing inferCNV input for sample: ", sample_id)

  # 1. Subset to the specific sample
  obj_s <- subset(obj, subset = (!!as.name(sample_col)) == sample_id)
  
  # [SEURAT V5] Ensure layers are joined
  obj_s <- JoinLayers(obj_s)

  # 2. IDENTIFY AND REMOVE TINY GROUPS
  # Check group sizes
  group_counts <- table(obj_s[[group_col]])
  message("  Original Group Sizes:")
  print(group_counts)
  
  # Find groups to keep (must have >= min_cells_per_group)
  groups_to_keep <- names(group_counts)[group_counts >= min_cells_per_group]
  groups_to_drop <- names(group_counts)[group_counts < min_cells_per_group]
  
  if (length(groups_to_drop) > 0) {
    message("  !!! DROPPING groups with < ", min_cells_per_group, " cells: ", 
            paste(groups_to_drop, collapse = ", "))
    
    # Subset the object to keep only valid groups
    # We use the cell names that belong to the valid groups
    cells_to_keep <- colnames(obj_s)[obj_s[[group_col]][,1] %in% groups_to_keep]
    obj_s <- subset(obj_s, cells = cells_to_keep)
  }

  message("  Final Cell Count for Analysis: ", ncol(obj_s))
  
  if (ncol(obj_s) < 10) {
    message("  --> Too few cells remaining. Skipping.")
    return(NULL)
  }

  # 3. Extract Counts & SANITIZE NAMES
  counts_s <- GetAssayData(obj_s, assay = "RNA", layer = "counts")
  
  # Replace hyphens with dots to prevent mismatch errors
  clean_names <- gsub("-", ".", colnames(counts_s))
  colnames(counts_s) <- clean_names
  
  # 4. Create Output Directory
  out_dir <- file.path(out_root, sample_id)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # 5. Write Counts Matrix
  counts_file <- file.path(out_dir, paste0(sample_id, ".counts.matrix.txt"))
  write.table(
    as.matrix(counts_s),
    file = counts_file,
    sep = "\t",
    quote = FALSE,
    col.names = NA
  )

  # 6. Write Annotations
  annots <- data.frame(
    cell_id = clean_names,
    group   = obj_s[[group_col]][,1], # Ensure we get the vector
    stringsAsFactors = FALSE
  )
  
  annots_file <- file.path(out_dir, paste0(sample_id, ".annotations.txt"))
  write.table(
    annots,
    file = annots_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  return(list(counts_file = counts_file, annots_file = annots_file, out_dir = out_dir))
}

run_infercnv_for_sample <- function(
  sample_id,
  counts_file,
  annots_file,
  gene_order_file,
  normal_types,
  out_dir,
  min_ref_cells = 10,
  min_genes = 50,
  safe_min_group_size = 5,
  num_threads = 10  # <--- CHANGED TO 10 FOR STABILITY
) {
  message("Processing sample: ", sample_id)

  # 1. Read annotations
  annots <- read.table(annots_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(annots) <- c("cell_id", "group")

  # 2. Check counts for each potential reference group
  ref_counts <- table(annots$group[annots$group %in% normal_types])
  
  message("  Reference group sizes in annotation:")
  print(ref_counts)

  # 3. Filter small groups
  valid_ref_groups <- names(ref_counts)[ref_counts >= safe_min_group_size]
  total_valid_ref <- sum(ref_counts[valid_ref_groups])
  
  message("  Valid reference groups (>= ", safe_min_group_size, " cells): ", 
          paste(valid_ref_groups, collapse = ", "))

  if (length(valid_ref_groups) == 0) {
    message("  --> No valid reference groups large enough. Skipping sample.")
    return(NULL)
  }

  if (total_valid_ref < min_ref_cells) {
    message("  --> Total reference cells (", total_valid_ref, 
            ") < min_ref_cells (", min_ref_cells, "). Skipping sample.")
    return(NULL)
  }

  # 4. Create Object
  infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = counts_file,
    annotations_file  = annots_file,
    delim             = "\t",
    gene_order_file   = gene_order_file,
    ref_group_names   = valid_ref_groups,
    chr_exclude       = c("chrX", "chrY", "chrM")
  )

  expr <- infercnv_obj@expr.data
  if (is.null(dim(expr)) || any(dim(expr) < 2)) {
    message("  --> Matrix collapsed to vector or empty. Skipping.")
    return(NULL)
  }

  message("Running inferCNV for sample: ", sample_id)
  message("Using num_threads: ", num_threads)

  # 5. Run inferCNV
  # We use the passed 'num_threads' here.
  # If this still fails, you can set HMM=FALSE, but try threads=1 first.
  infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff            = 0.1,
    out_dir           = file.path(out_dir, "infercnv_out"),
    cluster_by_groups = TRUE,
    denoise           = TRUE,
    HMM               = TRUE,
    num_threads       = num_threads 
  )

  return(infercnv_obj)
}

# Loop over all samples
infercnv_results <- list()

# Keep your existing run_infercnv_for_sample (the one with safe_min_group_size)
# And run the loop:

for (sid in sample_ids) {
  message("------------------------------------------------------")
  message("Processing sample: ", sid)

  out_dir_sid  <- file.path(infercnv_root, sid)
  done_file    <- file.path(out_dir_sid, "infercnv_out", "infercnv.observation_groupings.txt")

  if (file.exists(done_file)) {
    message("  --> inferCNV output already exists. Skipping.")
    next
  }

  # 1. Prepare Files (Removes 1-cell groups automatically)
  prep <- prepare_infercnv_input_for_sample(
    obj = adult_obj,
    sample_id = sid,
    sample_col = "sample_uid",
    group_col = "general_cell_type",
    out_root = infercnv_root,
    min_cells_per_group = 5  # Strict filter
  )

  # If prep returns NULL (too few cells or all groups dropped), skip
  if (is.null(prep)) {
    message("  --> Skipping execution for ", sid)
    next
  }

  # 2. Run inferCNV
  infercnv_results[[sid]] <- run_infercnv_for_sample(
    sample_id       = sid,
    counts_file     = prep$counts_file,
    annots_file     = prep$annots_file,
    gene_order_file = gene_order_file,
    normal_types    = normal_types,
    out_dir         = prep$out_dir,
    num_threads     = 10,  # Use 10 threads for stability
    safe_min_group_size = 5 # Matches the prep filter
  )
}



# Map back to Seurat object 
# For one 
adult_obj <- infercnv::add_to_seurat(
  seurat_obj            = adult_obj,
  assay_name            = "RNA",
  infercnv_output_path  = "/mnt/18T/chibao/gliomas/analysis/infercnv_adult/P1/infercnv_out",
  top_n                 = 10,
  column_prefix         = "infercnv_P1_"
)

# For all
for (sid in sample_ids) {
  out_path_sid <- file.path(infercnv_root, sid, "infercnv_out")

  adult_obj <- infercnv::add_to_seurat(
    seurat_obj           = adult_obj,
    assay_name           = "RNA",
    infercnv_output_path = out_path_sid,
    top_n                = 10,
    column_prefix        = paste0("infercnv_", sid, "_")
  )
}

head(colnames(adult_obj@meta.data))

# Summary file
# Identify all inferCNV columns
inf_cols <- grep("^infercnv_", colnames(adult_obj@meta.data), value = TRUE)

# Simple CNV burden: number of CNV segments per cell
adult_obj$infercnv_cnv_count <- rowSums(adult_obj@meta.data[, inf_cols, drop = FALSE] != 0)

# Or a weighted sum (if the columns store magnitude)
adult_obj$infercnv_cnv_sum <- rowSums(abs(adult_obj@meta.data[, inf_cols, drop = FALSE]))

summary(adult_obj$infercnv_cnv_count)
summary(adult_obj$infercnv_cnv_sum)

# Visualize
DefaultAssay(adult_obj) <- "SCT"

p_inf <- DimPlot(
  adult_obj,
  reduction = "umap.harmony",
  group.by = "infercnv_status",
  label = TRUE
) + ggtitle("inferCNV CNV burden on Harmony UMAP")

p_inf
