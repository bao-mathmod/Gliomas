library(Seurat)
library(SeuratObject)
library(tidyverse)
library(ggplot2)
library(future)
library(glmGamPoi)
library(presto)
library(viridis)
# library(ggpubr)
library(scales)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(purrr)
library(tidyr)

# == Configuration ==========================================================
plan("multicore", workers = 5)
options(future.globals.maxSize = 450 * 1024^3) # 450 GB
options(Seurat.object.assay.version = "v5")
set.seed(1234) # for reproducibility

# Load the integrated cleaned object
# obj <- readRDS('/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/harmony_cleaned_annotated.rds')
obj <- readRDS('/mnt/18T/chibao/gliomas/data_official/00_raw_data_adult_GBM/03_integrated/official/harmony/obj/harmony_integrated_cohort_GBM_orig_ident.rds')
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
for (r in c("pca","harmony","umap","umap.rpca", 'umap.harmony')) {
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
saveRDS(backup_SCT,'/mnt/18T/chibao/gliomas/data_official/02_myeloid/new/official/object/01_myeloid_cohort_GBM_SCT.rds')

myeloid <- readRDS('/mnt/18T/chibao/gliomas/data_official/02_myeloid/new/official/object/01_myeloid_cohort_GBM_SCT.rds')
out_dir <- '/mnt/18T/chibao/gliomas/data_official/02_myeloid/new/official'

# Check the npcs suit for downstream step 
myeloid <- RunPCA(myeloid, assay = "SCT", npcs = 50, verbose = FALSE)
elbow_plot <- ElbowPlot(myeloid, ndims = 50, reduction = "pca")
ggsave(file.path(out_dir, "elbow_plot.png"), plot = elbow_plot, width = 6, height = 4)

# For UMAP
myeloid <- RunUMAP(myeloid, dims = 1:30, verbose = FALSE)

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
  dims = 1:30, # Check with elbow output
  verbose = TRUE)

harmony_obj <- RunUMAP(harmony_obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
harmony_obj <- FindNeighbors(harmony_obj, reduction = "harmony", dims = 1:30)
harmony_obj <- FindClusters(
  harmony_obj,
  #graph.name = "SCT_snn",
  resolution = c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.4, 0.5, 0.6,  0.7, 0.8, 1.0, 1.2),
  # algorithm = 1,     # Louvain (stable), switch to 2 (SLM) if desired
  verbose = FALSE
)

p1 <- DimPlot(harmony_obj, reduction = 'umap.harmony', group.by = 'SCT_snn_res.0.08', label = TRUE)
p2 <- DimPlot(harmony_obj, reduction = 'umap.harmony', group.by = 'SCT_snn_res.0.1', label = TRUE)
p3 <- DimPlot(harmony_obj, reduction = 'umap.harmony', group.by = 'SCT_snn_res.0.2', label = TRUE)
p4 <- DimPlot(harmony_obj, reduction = 'umap.harmony', group.by = 'SCT_snn_res.0.4', label = TRUE)

p <- p1 + p2 + p3 + p4


ggsave(file.path(fig_dir, 'res.png'), plot = p, width = 14, height = 10)


# Save the final object
# saveRDS(myeloid, '/mnt/18T/chibao/gliomas/data_official/02_myeloid/01_myeloid_integrated_obj.rds')
# saveRDS(harmony_obj, '/mnt/18T/chibao/gliomas/data_official/02_myeloid/new/harmony/obj/01_myeloid_integrated_cohort_new_harmony_obj.rds')

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
  dims = 1:35,
  verbose = FALSE
)

rpca_obj <- RunUMAP(rpca_obj, reduction = "rpca", dims = 1:35, reduction.name = "umap.rpca")
rpca_obj <- FindNeighbors(rpca_obj, reduction = "rpca", dims = 1:35)
rpca_obj <- FindClusters(
  rpca_obj,
  #graph.name = "SCT_snn",
  resolution = c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.4, 0.5, 0.6,  0.7, 0.8, 1.0, 1.2),
  # algorithm = 1,     # Louvain (stable), switch to 2 (SLM) if desired
  verbose = FALSE
)

# Save the final object
# saveRDS(rpca_obj, '/mnt/18T/chibao/gliomas/data_official/02_myeloid/new/rpca/01_myeloid_integrated_rpca_obj.rds')

####################### GO PATHWAY ANALYSIS #######################
# Add markers 
cluster0_genes <- harmony_markers_0.2 |> filter(cluster == 0 & p_val_adj < 0.05 & avg_log2FC > 1) |> pull(gene)
cluster1_genes <- harmony_markers_0.2 |> filter(cluster == 1 & p_val_adj < 0.05 & avg_log2FC > 1) |> pull(gene)
cluster2_genes <- harmony_markers_0.2 |> filter(cluster == 2 & p_val_adj < 0.05 & avg_log2FC > 1) |> pull(gene)
cluster3_genes <- harmony_markers_0.2 |> filter(cluster == 3 & p_val_adj < 0.05 & avg_log2FC > 1) |> pull(gene)
cluster4_genes <- harmony_markers_0.2 |> filter(cluster == 4 & p_val_adj < 0.05 & avg_log2FC > 1) |> pull(gene)
cluster5_genes <- harmony_markers_0.2 |> filter(cluster == 5 & p_val_adj < 0.05 & avg_log2FC > 1) |> pull(gene)
cluster6_genes <- harmony_markers_0.2 |> filter(cluster == 6 & p_val_adj < 0.05 & avg_log2FC > 1) |> pull(gene)
cluster7_genes <- harmony_markers_0.2 |> filter(cluster == 7 & p_val_adj < 0.05 & avg_log2FC > 1) |> pull(gene)
cluster8_genes <- harmony_markers_0.2 |> filter(cluster == 8 & p_val_adj < 0.05 & avg_log2FC > 1) |> pull(gene)
cluster9_genes <- harmony_markers_0.2 |> filter(cluster == 9 & p_val_adj < 0.05 & avg_log2FC > 1) |> pull(gene)
cluster10_genes <- harmony_markers_0.2 |> filter(cluster == 10 & p_val_adj < 0.05 & avg_log2FC > 1) |> pull(gene)
cluster11_genes <- harmony_markers_0.2 |> filter(cluster == 11 & p_val_adj < 0.05 & avg_log2FC > 1) |> pull(gene)
cluster12_genes <- harmony_markers_0.2 |> filter(cluster == 12 & p_val_adj < 0.05 & avg_log2FC > 1) |> pull(gene)
cluster13_genes <- harmony_markers_0.2 |> filter(cluster == 13 & p_val_adj < 0.05 & avg_log2FC > 1) |> pull(gene)
cluster14_genes <- harmony_markers_0.2 |> filter(cluster == 14 & p_val_adj < 0.05 & avg_log2FC > 1) |> pull(gene)
cluster15_genes <- harmony_markers_0.2 |> filter(cluster == 15 & p_val_adj < 0.05 & avg_log2FC > 1) |> pull(gene)
cluster16_genes <- harmony_markers_0.2 |> filter(cluster == 16 & p_val_adj < 0.05 & avg_log2FC > 1) |> pull(gene)
cluster17_genes <- harmony_markers_0.2 |> filter(cluster == 17 & p_val_adj < 0.05 & avg_log2FC > 1) |> pull(gene)
cluster18_genes <- harmony_markers_0.2 |> filter(cluster == 18 & p_val_adj < 0.05 & avg_log2FC > 1) |> pull(gene)
cluster19_genes <- harmony_markers_0.2 |> filter(cluster == 19 & p_val_adj < 0.05 & avg_log2FC > 1) |> pull(gene)

# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(
    cluster8_genes,
    fromType = "SYMBOL", 
    toType = "ENTREZID", 
    OrgDb = org.Hs.eg.db)

# GO enrichment analysis
go_result <- enrichGO(gene = entrez_ids$ENTREZID,
                                         OrgDb = org.Hs.eg.db,
                                         ont = "BP",  # Biological Process
                                         pAdjustMethod = "BH",
                                         pvalueCutoff = 0.05,
                                         qvalueCutoff = 0.05
                                         )

# Visualize GO results
p <- dotplot(go_result, showCategory = 10) + 
    ggtitle("GO Enrichment Analysis for Cluster 8") +
    theme(plot.title = element_text(size = 16, face = "bold"),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          axis.title = element_text(size = 14)) +
    scale_size(range = c(3, 8))

# Save
fig_dir <- "/mnt/18T/chibao/gliomas/data_official/02_myeloid/new/official/figures" 
ggsave(file.path(fig_dir, "go_enrichment_cluster8.png"), plot = p, width = 8, height = 8)

#################### Marker Visualization ##################
# Feature Plot
p <- FeaturePlot(harmony_obj, features = c('CXCL3', 'CXCL2', 'FN1', 'VCAN'), reduction = 'umap.harmony', ncol = 2, cols = c('lightgrey', 'red'))
ggsave(file.path(fig_dir, "feature_plot_marker_8.png"), plot = p, width = 14, height = 8)

# Dot Plot
# top_markers <- c('CD74', 'LYZ', 'VCAN', "JAML", 'ITGAM', 'CD68', 'CD163', 'CD14', 'FCGRB3A', 'FCGR3B', 'TREM2', 'TMEM119', "P2RY12")
# top_markers <- c('SPP1', 'NUPR1', 'MIF', 'LGALS1', 'CSTB', 'MT1G', 'MT1E', 'MT1X', 'MT1F')
top_markers <- c('CXCL3', 'CXCL2', 'FN1', 'VCAN')
# Flatten the markers for DotPlot
markers_genes <- unlist(top_markers)
markers_genes <- make.unique(markers_genes)

# Set factor levels in the order of your original markers vector
markers_genes <- factor(markers_genes, levels = markers_genes)

# Create a mapping from gene -> "Cluster0:Gene" label
label_map <- setNames(
  paste0("ClusterX:", markers_genes),
  markers_genes
)

Idents(harmony_obj) <- 'SCT_snn_res.0.2' # Ensure the correct cluster identities are set
p <- DotPlot(harmony_obj, features = markers_genes, cluster.idents = TRUE) +
  scale_color_gradient(low = "lightgrey", high = "red") +
  theme_minimal() +
  xlab("Marker Genes") +
  ylab("Cluster ID") +
  ggtitle("Dot Plot of All Marker Genes Across Clusters") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    plot.title  = element_text(hjust = 0.5)
  ) +
  # Override labels so "CD3D" is displayed as "T_cells:CD3D", etc.
  scale_x_discrete(labels = label_map)

ggsave(file.path(fig_dir, "dot_plot_8.png"), plot = p, width = 8, height = 8)


################### Dot Plot Enhance ##################
markers <- c(
    # Cluster 5
    #'SPP1', 'NUPR1', 'MIF', 'LGALS1', 'CSTB',
    # Cluster 17
    #'MT1G', 'MT1E', 'MT1X', 'MT1F'

    # X
    'HMOX1', 'GSTP1', 'UQCRB'
)

p <- DotPlot(
  object = harmony_obj,
  features = markers,
  group.by = "SCT_snn_res.0.2"
) +
  scale_color_gradientn(
    colors = c("#2166AC", "white", "#B2182B"), 
    limits = c(-2.5, 2.5)
  ) +
  scale_size(range = c(3, 10)) +
  RotatedAxis() +
  labs(
    title = "Expression of Markers Across Cell Types",
    x = "Marker",
    y = "Cell Types",
    color = "Scaled Expression",
    size = "Percent Expressing"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 24, hjust = 0.5),
    axis.title.x = element_text(face = "bold", size = 20),
    axis.title.y = element_text(face = "bold", size = 20),
    axis.text.x = element_text(color = "black", size = 16, angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black", size = 16),
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 10),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1)
)

ggsave(file.path(fig_dir, "dot_plot_3.png"), plot = p, width = 8, height = 8)





# # Remove cluster contaminated 
# # clusters to remove
# rem <- c("9", "14", "16", "18", "19")

# # keep everything else
# keep <- setdiff(unique(Idents(harmony_obj)), rem)

# harmony_obj_cleaned <- subset(harmony_obj, idents = keep)

# # verify
# table(Idents(harmony_obj_cleaned))

# backup_harmony <- harmony_obj
# harmony_obj <- harmony_obj_cleaned