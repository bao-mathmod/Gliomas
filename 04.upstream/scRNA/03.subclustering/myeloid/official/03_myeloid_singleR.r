suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(ggplot2)
  library(pheatmap)
  library(BiocParallel)
})

# ==============================================================================
# USER SETTINGS
# ==============================================================================
sample_id <- "Glioma_Myeloid"
base_dir <- "/mnt/18T/chibao/gliomas/data_official/02_myeloid/new/rpca/obj"
out <- '/mnt/18T/chibao/gliomas/data_official/stuff'
# Input: The integrated myeloid Seurat object
in_path <- file.path(base_dir, "01_myeloid_integrated_rpca_obj.rds")

# Output: Where to save the annotated results
out_dir <- file.path(out, "02_Annotation_outputs", sample_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Processing Power
num_cores <- 8
bpparam <- MulticoreParam(workers = num_cores)

message("\n=== Starting Annotation for: ", sample_id, " ===")

# ==============================================================================
# STEP 1: Load Seurat Data & Prepare Matrix
# ==============================================================================
message("Loading Seurat Object...")
seu <- readRDS(in_path)
seu

# Extract the SCT data layer (or RNA data layer if SCT isn't present).
if ("SCT" %in% names(seu@assays)) {
  mat <- GetAssayData(seu, assay = "SCT", layer = "data")
} else {
  mat <- GetAssayData(seu, assay = "RNA", layer = "data")
}

message("Matrix extracted. Cells: ", ncol(mat), " | Genes: ", nrow(mat))

# ==============================================================================
# STEP 2: Load Reference Datasets (celldex)
# ==============================================================================
message("Downloading/Loading celldex references...")
# 1. Blueprint Encode (Excellent for deep myeloid lineage profiles)
ref.blueprint <- celldex::BlueprintEncodeData()

# 2. Human Primary Cell Atlas (Crucial because it contains Microglia and tissue Macrophages)
ref.hpca <- celldex::HumanPrimaryCellAtlasData()

# ==============================================================================
# STEP 3: Run SingleR (Multi-Reference Mode)
# ==============================================================================
message("Running SingleR Annotation (This may take a few minutes)...")

# Use "label.fine" for detailed sub-typing (e.g., separating myeloid subsets)
pred <- SingleR(
  test = mat,
  ref = list(Blueprint = ref.blueprint, HPCA = ref.hpca),
  labels = list(ref.blueprint$label.fine, ref.hpca$label.fine),
  BPPARAM = bpparam
)

message("SingleR Annotation Complete!")

# ==============================================================================
# STEP 4: Diagnostics & Quality Control
# ==============================================================================
message("Generating Diagnostic Plots...")

# 1. Score Heatmap
png(file.path(out_dir, paste0(sample_id, "_SingleR_Score_Heatmap.png")), width = 1200, height = 800, res = 150)
plotScoreHeatmap(pred, show.labels = TRUE, max.labels = 40)
dev.off()

# 2. Delta Distribution
png(file.path(out_dir, paste0(sample_id, "_SingleR_Delta_Scores.png")), width = 800, height = 600, res = 150)
plotDeltaDistribution(pred)
dev.off()

# ==============================================================================
# STEP 5: Integrate Back into Seurat
# ==============================================================================
message("Adding labels to Seurat metadata...")

# Add the final predicted fine labels
fine_labels <- pred$pruned.labels
fine_labels[is.na(fine_labels)] <- "Unknown"
seu$SingleR_Fine <- fine_labels

# Add the main (broad) labels for high-level overviews
pred_main <- SingleR(test = mat, ref = ref.hpca, labels = ref.hpca$label.main, BPPARAM = bpparam)
main_labels <- pred_main$pruned.labels
main_labels[is.na(main_labels)] <- "Unknown"
seu$SingleR_Main <- main_labels

# ==============================================================================
# STEP 6: Domain-Specific Glioma Profiling (Microglia vs. TAMs)
# ==============================================================================
message("Scoring cells for Microglia and blood-derived TAM signatures...")

# Classic Microglia resting/homeostatic markers
microglia_genes <- list(c("CX3CR1", "P2RY12", "TMEM119", "HEXB", "TREM2"))

# Classic Blood-derived Macrophage/TAM markers in glioma
tam_genes <- list(c("CD163", "CD68", "MRC1", "ITGAM", "TGFB1"))

# Score Microglia
seu <- AddModuleScore(object = seu, features = microglia_genes, name = "Microglia_Score")
colnames(seu@meta.data)[colnames(seu@meta.data) == "Microglia_Score1"] <- "Microglia_Score"

# Score TAMs
seu <- AddModuleScore(object = seu, features = tam_genes, name = "TAM_Score")
colnames(seu@meta.data)[colnames(seu@meta.data) == "TAM_Score1"] <- "TAM_Score"

# ==============================================================================
# STEP 7: Visualization
# ==============================================================================
message("Generating UMAPs...")

# Note: Since this is an integrated object, ensure your reduction name is correct 
# (e.g., it might be "umap.harmony" or "integrated_umap" depending on your pipeline).
umap_reduction <- "umap.rpca" 

# 1. UMAP by Broad Cell Type
p1 <- DimPlot(seu, reduction = 'umap.rpca', group.by = "SingleR_Main", label = TRUE, repel = TRUE) + 
  ggtitle(paste0(sample_id, " Broad Cell Types")) + NoLegend()
ggsave(file.path(out_dir, paste0(sample_id, "_UMAP_SingleR_Main.png")), p1, width = 7, height = 6)

# 2. UMAP by Fine Cell Type
p2 <- DimPlot(seu, reduction = 'umap.rpca', group.by = "SingleR_Fine", label = TRUE, repel = TRUE, label.size = 3) + 
  ggtitle(paste0(sample_id, " Detailed Myeloid Subtypes")) + NoLegend()
ggsave(file.path(out_dir, paste0(sample_id, "_UMAP_SingleR_Fine.png")), p2, width = 8, height = 7)

# 3. FeaturePlots of Microglia vs TAM Scores
p3 <- FeaturePlot(seu, features = c("Microglia_Score", "TAM_Score"), reduction = 'umap.harmony', ncol = 2) & 
  scale_color_viridis_c(option = "magma")
ggsave(file.path(out_dir, paste0(sample_id, "_UMAP_Microglia_vs_TAM.png")), p3, width = 12, height = 5)

# ==============================================================================
# STEP 8: Save Outputs
# ==============================================================================
message("Saving final annotated object...")

# Save RDS
saveRDS(seu, file.path(out_dir, paste0(sample_id, "_Annotated_seurat.rds")))

# Save CSV Summary of Cell Counts
cell_counts <- as.data.frame(table(seu$SingleR_Fine))
colnames(cell_counts) <- c("Cell_Type", "Count")
cell_counts <- cell_counts[order(-cell_counts$Count), ]
write.csv(cell_counts, file.path(out_dir, paste0(sample_id, "_CellType_Counts.csv")), row.names = FALSE)

message("✅ DONE! Outputs saved to: ", out_dir)


###############################################
suppressPackageStartupMessages({
  library(Seurat)
  library(Azimuth)
  library(SeuratData)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

# ==============================================================================
# SETTINGS
# ==============================================================================
sample_id <- "Glioma_Myeloid_Azimuth"
base_dir <- "/mnt/18T/chibao/gliomas/data_official/02_myeloid"

# Input: The integrated myeloid Seurat object
in_path <- file.path(base_dir, "01_myeloid_integrated_obj.rds")

# Output: Directory for Azimuth results
out_dir <- file.path(base_dir, "03_Annotation_outputs_Azimuth", sample_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("\n=== Starting Azimuth Annotation for: ", sample_id, " ===")

# ==============================================================================
# STEP 1: Load Data & Format for Azimuth
# ==============================================================================
message("Loading Seurat Object...")
seu <- readRDS(in_path)
seu

# Azimuth explicitly requires the 'RNA' assay (raw counts) to run its mapping.
DefaultAssay(seu) <- "RNA"
seu
backup <- seu

message("Data loaded. Cells: ", ncol(seu))

# ==============================================================================
# STEP 2: Run Azimuth (Using Human Motor Cortex Reference)
# ==============================================================================
message("Running Azimuth Brain Reference Mapping...")

# Using the human motor cortex reference to capture Microglia and CNS Macrophages
seu <- RunAzimuth(
  query = seu, 
  reference = "pbmcref",
  k.weight = 20,
  verbose = TRUE
)

message("Azimuth mapping complete!")

# ==============================================================================
# STEP 3: Quality Control & Diagnostics
# ==============================================================================
message("Generating QC Plots...")

# Note: The brain reference uses different metadata columns than the PBMC reference.
# It uses 'predicted.class', 'predicted.subclass', and 'predicted.cluster'.

p_score <- VlnPlot(seu, features = "mapping.score", group.by = "predicted.subclass", pt.size = 0.1) +
  ggtitle("Azimuth Mapping Confidence by Cell Subclass") +
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "red") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_dir, paste0(sample_id, "_Azimuth_Mapping_Scores.png")), p_score, width = 10, height = 6)

# ==============================================================================
# STEP 4: Visualization
# ==============================================================================
message("Generating UMAPs...")

# 1. Broad Class (e.g., Non-Neuronal vs Neuronal)
p_c1 <- DimPlot(seu, reduction = "ref.umap", group.by = "predicted.class", label = TRUE, label.size = 4, repel = TRUE) +
  ggtitle(paste0(sample_id, " Broad Class")) + NoLegend()
ggsave(file.path(out_dir, paste0(sample_id, "_UMAP_Azimuth_Class.png")), p_c1, width = 7, height = 6)

# 2. Detailed Subclass (e.g., Microglia, Perivascular Macrophages)
p_c2 <- DimPlot(seu, reduction = "ref.umap", group.by = "predicted.subclass", label = TRUE, label.size = 3, repel = TRUE) +
  ggtitle(paste0(sample_id, " Subclass (Microglia/Macrophages)")) + NoLegend()
ggsave(file.path(out_dir, paste0(sample_id, "_UMAP_Azimuth_Subclass.png")), p_c2, width = 9, height = 7)

# ==============================================================================
# STEP 5: Save Outputs
# ==============================================================================
message("Saving annotated object and summaries...")

# Save RDS
saveRDS(seu, file.path(out_dir, paste0(sample_id, "_Azimuth_annotated_seurat.rds")))

# Create a master CSV containing Barcode predictions and mapping scores
anno_df <- data.frame(
  barcode = colnames(seu),
  Class = seu$predicted.class,
  Subclass = seu$predicted.subclass,
  Cluster = seu$predicted.cluster,
  Mapping_Score = seu$mapping.score
)
write.csv(anno_df, file.path(out_dir, paste0(sample_id, "_Azimuth_Barcodes.csv")), row.names = FALSE)

# Create a summary count of Subclass cell types
cell_counts <- as.data.frame(table(seu$predicted.subclass))
colnames(cell_counts) <- c("Cell_Subclass", "Count")
cell_counts <- cell_counts %>% arrange(desc(Count))
write.csv(cell_counts, file.path(out_dir, paste0(sample_id, "_Azimuth_CellCounts_Subclass.csv")), row.names = FALSE)

message("✅ DONE! Azimuth outputs saved to: ", out_dir)