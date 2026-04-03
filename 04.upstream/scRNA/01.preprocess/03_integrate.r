#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(future)
  library(presto)
  library(harmony) # Recommended alternative method
  library(glmGamPoi)
  library(SeuratData)
  library(SeuratWrappers)
  library(Azimuth)
})

# == Configuration ==========================================================
plan("multicore", workers = 1)
options(future.globals.maxSize = 500 * 1024^3) # 450 GB
options(Seurat.object.assay.version = "v5")
set.seed(1234) # for reproducibility

# == Paths ==================================================================
in_dir  <- "/mnt/18T/chibao/gliomas/data_official/00_raw_data_adult_GBM/2_QC_output/post_QC_newparam/official/cohort_official_GBM/rds"
out_dir <- "/mnt/18T/chibao/gliomas/data_official/00_raw_data_adult_GBM/03_integrated/official/harmony/obj"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# == 1. Load, Merge, and Pre-process Data ====================================
message("Step 1: Loading individual Seurat objects...")
rds_files <- list.files(in_dir, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
objs_list <- lapply(rds_files, readRDS)

message("Merging objects into a single Seurat object with layers...")
# This creates a single object where each original object is a layer in the 'RNA' assay
# We use the 'sample_uid' from your metadata to name the layers
layer_names <- sapply(objs_list, function(x) unique(x$sample_uid))
merged_obj <- merge(x = objs_list[[1]], y = objs_list[2:length(objs_list)], add.cell.ids = layer_names)
merged_obj
# Join 
DefaultAssay(merged_obj) <- "RNA"
merged_obj
merged_obj <- JoinLayers(merged_obj)
merged_obj

# Split the RNA assay by the 'orig.ident' or 'sample_uid' which now corresponds to your layers
# merged_obj[["RNA"]] <- split(merged_obj[["RNA"]], f = merged_obj$sample_uid)
merged_obj[["RNA"]] <- split(merged_obj[["RNA"]], f = merged_obj$orig.ident)
merged_obj
# Clean up memory
rm(objs_list)
gc()

DefaultAssay(merged_obj) <- "RNA"
message("Performing Cell Cycle Scoring...")
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
merged_obj <- NormalizeData(merged_obj, verbose = FALSE)
merged_obj <- CellCycleScoring(merged_obj, s.features = s.genes, g2m.features = g2m.genes)

# == 2. SCTransform Normalization ============================================
message("Step 2: Running SCTransform on each layer...")
# SCTransform is run on the merged object, and it will automatically process each layer independently.
# We regress out cell cycle scores, a critical step for tumor data.
merged_obj <- SCTransform(merged_obj,
                          method = "glmGamPoi",
                          vars.to.regress = c("S.Score", "G2M.Score", "percent.mt"),
                          verbose = FALSE)

message("Step 2.1: Save file for Backup...")
# Save intermediate object for backup
saveRDS(merged_obj, file.path(out_dir, "merge_cohort_new_backup_SCT.rds"))

# For safety, reload the backup
# merged_obj <- readRDS('/mnt/18T/chibao/gliomas/data_official/00_raw_data_adult/03_integrated/merge_backup_SCT.rds')
merged_obj
backup <- merged_obj
# Run PCA on the SCT assay. Seurat will automatically perform this for each layer.
message("Running PCA and UMAP on each layer...")
merged_obj <- RunPCA(merged_obj, assay = "SCT", npcs = 100, verbose = FALSE)
elbow_plot <- ElbowPlot(merged_obj, ndims = 100, reduction = "pca")
ggsave(file.path(out_dir, "elbow_plot.png"), plot = elbow_plot, width = 6, height = 4)

# For UMAP
merged_obj <- RunUMAP(merged_obj, dims = 1:40, verbose = FALSE)

harmony_obj <- merged_obj

rpca_obj <- merged_obj


message("Step 2.2: Integrating layers using Harmony...")
harmony_obj <- IntegrateLayers(
  object = harmony_obj,
  method = HarmonyIntegration,
  normalization.method = "SCT",
  # reference = reference_datasets,
  orig.reduction = "pca",
  #k.weight = 50,
  new.reduction = "harmony", # Name of the new integrated reduction
  dims = 1:40, # Using more PCs can be beneficial for complex datasets
  verbose = TRUE
)
# == 4. Downstream Analysis and Clustering ===================================
message("Step 4: UMAP, Neighbors, and Clustering on integrated data...")
# Harmony
# harmony_obj <- RunPCA(harmony_obj, assay = "SCT", npcs = 100, verbose = FALSE)
harmony_obj <- RunUMAP(harmony_obj, reduction = "harmony", dims = 1:40, reduction.name = "umap.harmony")
harmony_obj <- FindNeighbors(harmony_obj, reduction = "harmony", dims = 1:40)
harmony_obj <- FindClusters(
  harmony_obj,
  #graph.name = "SCT_snn",
  resolution = c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.4, 0.5, 0.6,  0.7, 0.8, 1.0, 1.2),
  # algorithm = 1,     # Louvain (stable), switch to 2 (SLM) if desired
  verbose = FALSE
)
# Viz for Hamrony 
p1 <- DimPlot(harmony_obj, reduction = 'umap.harmony', group.by = 'SCT_snn_res.0.04', label = TRUE)
p2 <- DimPlot(harmony_obj, reduction = 'umap.harmony', group.by = 'SCT_snn_res.0.06', label = TRUE)
p3 <- DimPlot(harmony_obj, reduction = 'umap.harmony', group.by = 'SCT_snn_res.0.08', label = TRUE)

p <- p1 + p2 + p3


ggsave(file.path(out_dir, 'res.png'), plot = p, width = 16, height = 8)
ggsave(file.path(out_dir, 'dot.png'), plot = p, width = 12, height = 10)

# Prep markers
harmony_obj <- PrepSCTFindMarkers(harmony_obj, assay = "SCT", verbose = TRUE)
Idents(harmony_obj) <- 'SCT_snn_res.0.08'
harmony_markers_0.08 <- FindAllMarkers(harmony_obj,
                              assay = "SCT",
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25,
                              test.use = "wilcox")
harmony_markers_0.08 <- harmony_markers_0.08 |> arrange(desc(avg_log2FC))
write.csv(harmony_markers_0.08, file = '/mnt/18T/chibao/gliomas/data_official/00_raw_data_adult_GBM/03_integrated/official/harmony/markers/harmony_markers_GBM_0.08.csv')
p <- DimPlot(harmony_obj, reduction = 'umap.harmony', group.by = 'general_cell_type', label = TRUE)
ggsave(file.path(out_dir, 'annote.png'), plot = p)

saveRDS(harmony_obj, file.path(out_dir, "harmony_integrated_cohort_GBM_orig_ident.rds"))




# == 4. Downstream Analysis and Clustering ===================================
# RPCA
message("Step 2.2: Integrating layers using RPCA...")
rpca_obj <- IntegrateLayers(
  object = rpca_obj, method = RPCAIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca", new.reduction = "rpca", k.weight = 40,
  verbose = FALSE
)
# rpca_obj <- RunPCA(rpca_obj, assay = "SCT", npcs = 100, verbose = FALSE)
rpca_obj <- RunUMAP(rpca_obj, reduction = "rpca", dims = 1:40, reduction.name = "umap.rpca")
rpca_obj <- FindNeighbors(rpca_obj, reduction = "rpca", dims = 1:40)
rpca_obj <- FindClusters(
  rpca_obj,
  #graph.name = "SCT_snn",
  resolution = c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.4, 0.5, 0.6,  0.7, 0.8, 1.0, 1.2),
  # algorithm = 1,     # Louvain (stable), switch to 2 (SLM) if desired
  verbose = FALSE
)

# saveRDS(harmony_obj, file.path(out_dir, "harmony_integrated_orig.ident.rds"))
saveRDS(rpca_obj, "/mnt/18T/chibao/gliomas/data_official/00_raw_data_adult/03_integrated/02_rpca/orig.ident/rpca_integrated_orig.ident.rds")
