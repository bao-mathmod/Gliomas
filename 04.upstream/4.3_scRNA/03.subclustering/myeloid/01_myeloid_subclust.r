library(Seurat)
library(SeuratObject)
library(tidyverse)
library(ggplot2)
library(future)
library(glmGamPoi)

# Ensure object convert to V5 format
options(Seurat.object.assay.version = 'v5')

# Load the integrated cleaned object
obj <- readRDS('/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/harmony_cleaned_annotated.rds')
obj

# Check the clean
Idents(obj) |> table()

# Subset myeloid cells
Idents(obj) <- obj$general_cell_type
myeloid <- subset(obj, idents = c('Myeloid'))
myeloid
myeloid@meta.data$general_cell_type |> unique()

# Check NA values in myeloid
sum(is.na(myeloid))

# Change to RNA assay
DefaultAssay(myeloid) <- "RNA"
myeloid

# Join Layers
myeloid <- JoinLayers(myeloid)
myeloid

# Split the RNA assay by the orig.ident 
myeloid[["RNA"]] <- split(myeloid[["RNA"]], f = myeloid$orig.ident)
myeloid

# Perform Cell Cycle Scoring
DefaultAssay(myeloid) <- "RNA"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
myeloid <- NormalizeData(myeloid, verbose = FALSE)
myeloid
myeloid <- CellCycleScoring(myeloid, s.features = s.genes, g2m.features = g2m.genes)

backup_obj <- myeloid
# Remove existence of SCT assays and any reductions related
if ("SCT" %in% Assays(myeloid)) myeloid[["SCT"]] <- NULL
for (r in c("pca","harmony","umap","umap.harmony")) {
  if (r %in% names(myeloid@reductions)) myeloid@reductions[[r]] <- NULL
}

myeloid

# Re-create SCT assay for myeloid
plan("multicore", workers = 20)
options(future.globals.maxSize = 400 * 1024^3) # 450 GB
options(Seurat.object.assay.version = "v5")
set.seed(1234) # for reproducibility

myeloid <- SCTransform(myeloid,
                       method = 'glmGamPoi',
                       vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
                       verbose = TRUE)
myeloid
# Save file for backup
saveRDS(myeloid,'/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/subclusters/myeloid/myeloid_SCT.rds')

# Run the Harmony Integration
myeloid <- RunPCA(myeloid, assay = "SCT", verbose = FALSE)
myeloid <- IntegrateLayers(
  object = myeloid,
  method = HarmonyIntegration,
  normalization.method = "SCT",
  # reference = reference_datasets,
  orig.reduction = "pca",
  #k.weight = 50,
  new.reduction = "harmony", # Name of the new integrated reduction
  dims = 1:30, # Using more PCs can be beneficial for complex datasets
  verbose = TRUE
)
myeloid <- RunUMAP(myeloid, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
myeloid <- FindNeighbors(myeloid, reduction = "harmony", dims = 1:30)
myeloid <- FindClusters(
  myeloid,
  #graph.name = "SCT_snn",
  resolution = c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.4, 0.5, 0.6,  0.7, 0.8, 1.0, 1.2),
  algorithm = 1,     # Louvain (stable), switch to 2 (SLM) if desired
  verbose = FALSE
)

# Save the final object
saveRDS(myeloid, '/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/subclusters/myeloid/myeloid_final_integrated.rds')

# Visualize
plot_dir <- '/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/subclusters/myeloid/myeloid_harmony.png'
p <- DimPlot(myeloid, reduction = 'umap.harmony', group.by = 'SCT_snn_res.0.05', label = TRUE)
ggsave(filename = plot_dir, plot = p)

