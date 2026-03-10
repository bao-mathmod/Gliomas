library(Seurat)
library(SeuratObject)
library(tidyverse)
library(ggplot2)
library(future)
library(glmGamPoi)

# Ensure object convert to V5 format
options(Seurat.object.assay.version = 'v5')

# Load the integrated cleaned object
# obj <- readRDS('/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/harmony_cleaned_annotated.rds')
obj <- readRDS('/mnt/18T/chibao/gliomas/data_official/01_integrated_obj/new/obj/02_rpca_integrated_orig.ident_annotated.rds')
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
for (r in c("pca","harmony","umap","umap.rpca")) {
  if (r %in% names(myeloid@reductions)) myeloid@reductions[[r]] <- NULL
}

myeloid

# Re-create SCT assay for myeloid
plan("multicore", workers = 10)
options(future.globals.maxSize = 400 * 1024^3) # 450 GB
options(Seurat.object.assay.version = "v5")
set.seed(1234) # for reproducibility

myeloid <- SCTransform(myeloid,
                       method = 'glmGamPoi',
                       vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"),
                       verbose = TRUE)
myeloid
backup_SCT <- myeloid


# Save file for backup
saveRDS(backup_SCT,'/mnt/18T/chibao/gliomas/data_official/02_myeloid/new/01_myeloid_SCT.rds')

out_dir <- '/mnt/18T/chibao/gliomas/data_official/02_myeloid/new'

# Check the npcs suit for downstream step 
myeloid <- RunPCA(myeloid, assay = "SCT", npcs = 50, verbose = FALSE)
elbow_plot <- ElbowPlot(myeloid, ndims = 50, reduction = "pca")
ggsave(file.path(out_dir, "elbow_plot.png"), plot = elbow_plot, width = 6, height = 4)

# For UMAP
myeloid <- RunUMAP(myeloid, dims = 1:20, verbose = FALSE)


# Set separate objects for harmony and rpca 
harmony_obj <- myeloid
rpca_obj <- myeloid

# Run the Harmony Integration
harmony_obj <- IntegrateLayers(
  object = harmony_obj,
  method = HarmonyIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca",
  new.reduction = "harmony", # Name of the new integrated reduction
  dims = 1:20, # Check with elbow output
  verbose = TRUE)

harmony_obj <- RunUMAP(harmony_obj, reduction = "harmony", dims = 1:20, reduction.name = "umap.harmony")
harmony_obj <- FindNeighbors(harmony_obj, reduction = "harmony", dims = 1:20)
harmony_obj <- FindClusters(
  harmony_obj,
  #graph.name = "SCT_snn",
  resolution = c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.4, 0.5, 0.6,  0.7, 0.8, 1.0, 1.2),
  # algorithm = 1,     # Louvain (stable), switch to 2 (SLM) if desired
  verbose = FALSE
)

# Save the final object
# saveRDS(myeloid, '/mnt/18T/chibao/gliomas/data_official/02_myeloid/01_myeloid_integrated_obj.rds')
saveRDS(harmony_obj, '/mnt/18T/chibao/gliomas/data_official/02_myeloid/new/harmony/01_myeloid_integrated_harmony_obj.rds')

# Visualize
# plot_dir <- '/mnt/18T/chibao/gliomas/data_official/02_myeloid/new/harmony'
# p <- DimPlot(harmony_obj, reduction = 'umap.harmony', group.by = 'SCT_snn_res.0.05', label = TRUE)
# ggsave(file.path(plot_dir, "umap_harmony.png"), plot = p)

############### Test with RPCA Integration ###############
# Run the RPCA Integration
rpca_obj <- IntegrateLayers(
  object = rpca_obj,
  method = RPCAIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca",
  new.reduction = "rpca", # Name of the new integrated reduction
  dims = 1:20,
  verbose = FALSE
)

rpca_obj <- RunUMAP(rpca_obj, reduction = "rpca", dims = 1:20, reduction.name = "umap.rpca")
rpca_obj <- FindNeighbors(rpca_obj, reduction = "rpca", dims = 1:20)
rpca_obj <- FindClusters(
  rpca_obj,
  #graph.name = "SCT_snn",
  resolution = c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.4, 0.5, 0.6,  0.7, 0.8, 1.0, 1.2),
  # algorithm = 1,     # Louvain (stable), switch to 2 (SLM) if desired
  verbose = FALSE
)

# Save the final object
saveRDS(rpca_obj, '/mnt/18T/chibao/gliomas/data_official/02_myeloid/new/rpca/01_myeloid_integrated_rpca_obj.rds')


