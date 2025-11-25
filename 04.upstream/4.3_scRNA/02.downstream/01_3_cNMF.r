# Load libraries
library(Seurat)
library(Matrix)
library(data.table)

# 1. Load your Final Object
# (Update path if needed)
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

# 2. Handle Seurat V5 Layers (CRITICAL STEP)
# cNMF needs a single matrix, not split layers. We join them back together.
# We focus on the "RNA" assay (Raw counts), not the "SCT" assay.
DefaultAssay(obj) <- "RNA"
obj
obj <- JoinLayers(obj)
obj

# 3. Feature Selection (Replicating Script 05 from Paper)
# The paper identifies variable genes on raw data to define biological signals.
# We select the top 3,000 genes (as per the paper's parameters for discovery).
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
variable_genes <- VariableFeatures(obj)

# 4. Define Output Directory
# Create a dedicated folder for cNMF to keep things clean
cnmf_dir <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/cNMF"
#if (!dir.exists(cnmf_dir)) { dir.create(cnmf_dir, recursive = TRUE) }

# 5. Export Data for cNMF
# We export the RAW counts (sparse matrix) and the list of variable genes.

# A. Export Variable Genes List
write.table(variable_genes, file = file.path(cnmf_dir, "variable_genes.txt"), 
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# B. Export the Count Matrix (.mtx format)
# This is compatible with the Python cNMF tool and memory efficient.
counts_matrix <- LayerData(obj, assay = "RNA", layer = "counts")
writeMM(counts_matrix, file = file.path(cnmf_dir, "matrix.mtx"))

# C. Export Barcodes and Genes (Required for .mtx inputs)
write.table(colnames(counts_matrix), file = file.path(cnmf_dir, "barcodes.tsv"),
            quote = FALSE, row.names = FALSE, col.names = FALSE)

genes <- rownames(LayerData(obj, assay = "RNA", layer = "counts"))

# Create a 2-column Data Frame (GeneID <tab> GeneSymbol)
# We just duplicate the names for both columns
genes_df <- data.frame(GeneID = genes, GeneSymbol = genes)

# Overwrite the genes.tsv file
# Critical: Use sep = "\t" to make it a tab-separated file
write.table(genes_df, file = file.path(cnmf_dir, "genes.tsv"),
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Check:
cnmf_dir <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/cNMF"

# Look at file headers
system(paste("head -n 5", file.path(cnmf_dir, "matrix.mtx")))
system(paste("head -n 5", file.path(cnmf_dir, "genes.tsv")))
system(paste("head -n 5", file.path(cnmf_dir, "barcodes.tsv")))


message("Data export complete. Proceed to cNMF execution.")