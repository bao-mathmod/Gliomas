library(Seurat)
library(SeuratObject)
library(tidyverse)
library(ggplot2)
library(future)
library(glmGamPoi)

# Ensure object convert to V5 format
options(Seurat.object.assay.version = 'v5')

# Load the integrated cleaned object
obj <- readRDS('/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/obj_full/harmony_cleaned_annotated_v3.rds')
obj
obj$general_cell_type |> table()

# Check the clean
Idents(obj) |> table()

# Subset tils cells
Idents(obj) <- obj$general_cell_type
tils <- subset(obj, idents = c('TILs'))
tils
tils@meta.data$general_cell_type |> unique()

# Check NA values in tils
sum(is.na(tils))

# Change to RNA assay
DefaultAssay(tils) <- "RNA"
tils
# Join Layers
tils <- JoinLayers(tils)
tils

# Split the RNA assay by the orig.ident 
tils[["RNA"]] <- split(tils[["RNA"]], f = tils$orig.ident)
tils

# Perform Cell Cycle Scoring
DefaultAssay(tils) <- "RNA"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
tils <- NormalizeData(tils, verbose = FALSE)
tils
tils <- CellCycleScoring(tils, s.features = s.genes, g2m.features = g2m.genes)

backup_obj <- tils
# Remove existence of SCT assays and any reductions related
if ("SCT" %in% Assays(tils)) tils[["SCT"]] <- NULL
for (r in c("pca","harmony","umap","umap.harmony")) {
  if (r %in% names(tils@reductions)) tils@reductions[[r]] <- NULL
}

tils

# Re-create SCT assay for tils
plan("multicore", workers = 1)
options(future.globals.maxSize = 350 * 1024^3) # 450 GB
options(Seurat.object.assay.version = "v5")
set.seed(1234) # for reproducibility

tils <- SCTransform(tils,
                       method = 'glmGamPoi',
                       vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
                       verbose = TRUE)
tils
backup_SCT <- tils
# Save file for backup
saveRDS(backup_SCT,'/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/subclusters/tils/tils_SCT.rds')

# Run the Harmony Integration
tils <- RunPCA(tils, assay = "SCT", verbose = FALSE)
tils <- IntegrateLayers(
  object = tils,
  method = HarmonyIntegration,
  normalization.method = "SCT",
  # reference = reference_datasets,
  orig.reduction = "pca",
  #k.weight = 50,
  new.reduction = "harmony", # Name of the new integrated reduction
  dims = 1:30, # Using more PCs can be beneficial for complex datasets
  verbose = TRUE
)
tils <- RunUMAP(tils, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
tils <- FindNeighbors(tils, reduction = "harmony", dims = 1:30)
tils <- FindClusters(
  tils,
  #graph.name = "SCT_snn",
  resolution = c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.4, 0.5, 0.6,  0.7, 0.8, 1.0, 1.2),
  algorithm = 1,     # Louvain (stable), switch to 2 (SLM) if desired
  verbose = FALSE
)

# Save the final object
saveRDS(tils, '/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/subclusters/tils/tils_final_integrated.rds')

# Visualize
plot_dir <- '/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/subclusters/tils/tils_harmony.png'
p <- DimPlot(tils, reduction = 'umap.harmony', group.by = 'SCT_snn_res.0.05', label = TRUE)
ggsave(filename = plot_dir, plot = p)



