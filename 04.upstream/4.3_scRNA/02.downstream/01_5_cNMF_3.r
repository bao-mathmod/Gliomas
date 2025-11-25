# ==============================================================================
# Step 3: Import cNMF Annotations into Seurat
# ==============================================================================

library(Seurat)
library(ggplot2)
library(dplyr)

# 1. Load Seurat Object (if not already loaded)
obj <- readRDS("/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/harmony_cleaned_annotated_v2.rds")

# Remove specified cell type annotations from the Seurat object
remove_types <- c("Melanoma_Metas", "Lung_Metas")
message("Removing annotations: ", paste(remove_types, collapse = ", "))

message("Counts before removal:")
print(table(obj$general_cell_type))

obj <- subset(obj, subset = !(general_cell_type %in% remove_types))

# Drop unused factor levels in metadata
# obj@meta.data$general_cell_type <- droplevels(obj@meta.data$general_cell_type)

# Ensure Idents reflect the current metadata column (optional but helpful for downstream code)
if ("general_cell_type" %in% colnames(obj@meta.data)) {
    Idents(obj) <- "general_cell_type"
}

message("Removal complete. Counts after removal:")
print(table(obj$general_cell_type))
# ------------------------------------------------------------------------------
# 2. Config Paths
cnmf_dir <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/cNMF"
run_name <- "Glioma_Adult_Myeloid"
k_val <- 18
plot_dir <- file.path(cnmf_dir, "plot.png")

# ------------------------------------------------------------------------------
# 3. Load the Usage Matrix (CRITICAL: Selecting the 0.02 file)
# ------------------------------------------------------------------------------
# Note the ".dt_0_02" in the filename matching your chosen threshold
usage_file <- file.path(cnmf_dir, run_name, 
                        paste0(run_name, ".usages.k_", k_val, ".dt_0_02.consensus.txt"))

if (!file.exists(usage_file)) {
  stop("File not found! Check if the density threshold in the filename matches your run.")
}

# Read the usage matrix
usages <- read.table(usage_file, header = TRUE, row.names = 1, sep = "\t")

# 4. Normalize Usages
# cNMF usages are raw values. We normalize them so each cell sums to 1 (100%).
usages_norm <- as.data.frame(t(apply(usages, 1, function(x) x / sum(x))))

# Rename columns to be descriptive (e.g., cNMF_1, cNMF_2...)
colnames(usages_norm) <- paste0("cNMF_", colnames(usages_norm))

# 5. Add Metadata to Seurat
# This matches cell barcodes and adds the 18 new columns to your object
obj <- AddMetaData(obj, metadata = usages_norm)

# ------------------------------------------------------------------------------
# 6. Quick Validation Plots
# ------------------------------------------------------------------------------

# A. Visualize the programs on your UMAP
# This will show you which programs correspond to T-cells, Tumor, etc.
# We plot the first 4 as an example.
p <- FeaturePlot(obj, features = colnames(usages_norm)[1:4], ncol = 2, min.cutoff = "q10", max.cutoff = "q90")
ggsave(filename = plot_dir, plot = p)
# B. Identify Gene Markers for the Programs
# Load the spectra file to see which genes define "cNMF_1", "cNMF_2", etc.
spectra_file <- file.path(cnmf_dir, run_name, 
                          paste0(run_name, ".gene_spectra_score.k_", k_val, ".dt_0_02.txt"))
spectra <- read.table(spectra_file, header = TRUE, row.names = 1, sep = "\t")
spectra <- as.matrix(spectra)
# Function to print top 10 genes for a specific program (e.g., Program 1)
print_program_genes <- function(prog_num) {
  top_genes <- names(sort(spectra[prog_num,], decreasing = TRUE))[1:20]
  cat(paste0("Program ", prog_num, ": ", paste(top_genes, collapse=", "), "\n"))
}

# Print markers for all 18 programs
for (i in 1:18) {
  print_program_genes(i)
}

# Save the annotated object
saveRDS(obj, "/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/harmony_annotated_cNMF.rds")

#######
# 1. Define your cNMF program columns (assuming you added them to metadata)
cnmf_features <- paste0("cNMF_X", 1:18)

# 2. DotPlot Visualization
# This shows the average usage of each program in each cell type
p1 <- DotPlot(obj, features = cnmf_features, group.by = "general_cell_type") + 
    RotatedAxis() +
    scale_color_gradientn(colours = c("white", "blue", "red")) +
    ggtitle("Mapping cNMF Programs to Cell Types") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10)) +
    theme(plot.margin = unit(c(-1, -2.5, -0.5, -1.5), "cm"))

ggsave(filename = plot_dir, plot = p1)