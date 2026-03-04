#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(stringr); library(future)
})
plan("multicore", workers = 8)
options(future.globals.maxSize = 60120 * 1024^2) # Sets limit to 80 GB
library(glmGamPoi)
# Ensure the object convert to V5
options(Seurat.object.assay.version = "v5") 

in_dir  <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/set1/rds"
out_dir <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/set1/integrated"; dir.create(out_dir, FALSE, TRUE)

objs <- lapply(list.files(in_dir, pattern="\\.rds$", full.names=TRUE, recursive = TRUE), readRDS)

# Proactively rename cells to ensure global uniqueness before integration
message("Renaming cells to be globally unique using 'sample_uid'...")
for (i in seq_along(objs)) {
  # It's good practice to ensure there's only one unique ID per object
  sample_identifier <- unique(objs[[i]]$sample_uid)
  if (length(sample_identifier) != 1) {
    stop("Object at index ", i, " has zero or more than one sample_uid.")
  }
  # Rename the cells by adding the sample_uid as a prefix
  objs[[i]] <- RenameCells(objs[[i]], add.cell.id = sample_identifier)
}

# Validate 
all_cell_names <- unlist(lapply(objs, colnames))

# 2. Check if there are any duplicated names
if (any(duplicated(all_cell_names))) {
  message("❌ VALIDATION FAILED: Duplicate cell names still exist after renaming.")
} else {
  message("✅ VALIDATION SUCCESS: All cell names are now globally unique. You are ready for the next step!")
}

# Updated to Seurat v5 integration workflow (RPCA + SCT)
objs <- lapply(objs, function(x) {
  # Update
  x <- UpdateSeuratObject(x) 
  
  # Check if RNA assay is updated to Assay5
   if (class(x@assays$RNA) != "Assay5") {
     message(paste0("Warning: ", unique(as.character(x$sample_id))[1]))
   }

  return(x)
})

# 1) DO NOT hard-intersect the full gene sets. Let Seurat pick shared features.
#    (SelectIntegrationFeatures ensures chosen features exist across all objects.)
#    If you truly want to be conservative, keep your intersect step,
#    but it's not required for SCT integration.

# 2) Ensure SCT exists (SCT v2 + glmGamPoi recommended for speed)
objs <- lapply(objs, function(x) {
  if (!"SCT" %in% names(x@assays)) {
    x <- SCTransform(x, vst.flavor="v2", method="glmGamPoi", verbose=FALSE)
  }
  x
})

# Rerun SCT to ensure consistency
## 0) Fresh SCT on all objects — consistent, post-HGNC
objs <- lapply(objs, function(x) {
  DefaultAssay(x) <- "RNA"
  x <- SCTransform(x, vst.flavor = "v2", method = "glmGamPoi",
                   return.only.var.genes = FALSE, verbose = FALSE)
  DefaultAssay(x) <- "SCT"
  x
})

# 1) Pick integration features from SCT and drop MT/ribo
features <- SelectIntegrationFeatures(object.list=objs, nfeatures=4000)
features <- features[!grepl("^MT-", features) & !grepl("^RP[SL]", features)]
features <- head(features, 3000)

objs <- lapply(objs, function(x) { VariableFeatures(x) <- features; x })

## (safety) If any stale scale.data exists with wrong rows, wipe it so Prep recomputes
objs <- lapply(objs, function(x) {
  sd <- LayerData(x, assay = "SCT", layer = "scale.data")
  if (!is.null(sd) && !identical(rownames(sd), features))
    LayerData(x, assay = "SCT", layer = "scale.data") <- NULL
  x
})

## 2) (optional but helpful) warn on very small samples
sizes <- vapply(objs, ncol, integer(1))
if (any(sizes < 200)) {
  message("Tiny objects (<200 cells) at indices: ",
          paste(which(sizes < 200), collapse = ", "),
          " — consider excluding from anchors & mapping later.")
}

# 3. Exclude tiny objects from anchor building
tiny_idx <- c(5) # change as needed
objs_main <- objs[-tiny_idx]
objs_tiny <- objs[tiny_idx]

# RPCA integration with SCT normalization
# 1) Feature selection + Prep (SCT path)
# 2. Integrate only the main cohort
features <- SelectIntegrationFeatures(object.list = objs, nfeatures = 3000)
objs <- PrepSCTIntegration(object.list = objs, anchor.features = features)

# In case: Have tiny objects only
features <- SelectIntegrationFeatures(object.list = objs_main, nfeatures = 3000)
objs_main <- PrepSCTIntegration(object.list = objs_main, anchor.features = features)


# 2) IMPORTANT in RPCA: RunPCA on EACH object using the same features (v5 vignette)
objs <- lapply(objs, function(x) RunPCA(x, features = features))

# In case: Have tiny objects only
objs_main <- lapply(objs_main, function(x) RunPCA(x, features = features))

# 3) Anchors (RPCA + SCT). Consider k.anchor tweak & reference-based integration for big cohorts
anchors <- FindIntegrationAnchors(object.list = objs_main, # Change to objs_main if tiny objects excluded
                                  normalization.method = "SCT",
                                  reduction = "rpca",
                                  anchor.features = features,
                                  dims = 1:30)

integrated <- IntegrateData(anchorset = anchors,
                            normalization.method = "SCT",
                            dims = 1:30)

DefaultAssay(integrated) <- "integrated"

# 4) Downstream on integrated assay
integrated <- RunPCA(integrated, npcs = 40, verbose = FALSE)
integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:30, return.model = TRUE, verbose = FALSE)

# Save for backup
saveRDS(integrated, file.path(out_dir, "scrna_integrated.rds"), compress = "xz")
# saveRDS(objs_tiny, file.path(out_dir, "scrna_tiny_objects.rds"), compress = "xz")
## 4.1) Map each tiny object separately
## Pre-req (reference):
## integrated already has PCA; if you want UMAP projection, run with return.model=TRUE earlier:
## integrated <- RunUMAP(integrated, reduction="pca", dims=1:30, return.model=TRUE)

map_one <- function(q, ref, dims = 1:30, min_feats = 500) {
  stopifnot("SCT" %in% names(q@assays), "SCT" %in% names(ref@assays))
  DefaultAssay(q)   <- "SCT"
  # keep PCA/UMAP tied to integrated assay in ref:
  DefaultAssay(ref) <- "integrated"

  # Build a shared feature set from SCT assays
  ref_sct_rows <- rownames(LayerData(ref, assay = "SCT", layer = "counts"))
  qry_sct_rows <- rownames(LayerData(q,   assay = "SCT", layer = "counts"))
  feats <- intersect(ref_sct_rows, qry_sct_rows)

  # (snRNA tip) drop MT/ribo to avoid noise
  feats <- feats[!grepl("^MT-", feats) & !grepl("^RP[SL]", feats)]

  # Fallback: if too small, intersect with ref SCT variable features
  if (length(feats) < min_feats) {
    ref_var <- VariableFeatures(ref[["SCT"]])
    feats <- intersect(feats, ref_var)
  }

  if (length(feats) < min_feats) {
    stop(sprintf("Shared SCT features too few for transfer (n=%d). Consider merging this tiny query or mapping after merging).", length(feats)))
  }

  tr.anchors <- FindTransferAnchors(
    reference             = ref,
    query                 = q,
    normalization.method  = "SCT",
    reference.assay       = "SCT",
    query.assay           = "SCT",
    reference.reduction   = "pca",   # uses PCA on ref (integrated assay)
    dims                  = dims,
    features              = feats
  )

  q_mapped <- MapQuery(
    anchorset            = tr.anchors,
    reference            = ref,
    query                = q,
    reference.reduction  = "pca",
    reduction.model      = "umap"    # works if ref UMAP was computed with return.model=TRUE
    # , refdata = list(celltype = ref$celltype) # add if you have labels to transfer
  )
  q_mapped
}

## Run for the tiny objects (list):
mapped_list <- lapply(objs_tiny, map_one, ref = integrated, dims = 1:30)

# 4.2) Merge mapped tiny objects into integrated
integrated_plus <- merge(integrated, y = mapped_list)
integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:30, verbose = FALSE)
integrated <- FindClusters(integrated, resolution = c(0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4), verbose = FALSE)

# integrated_plus <- FindNeighbors(integrated_plus, reduction = "pca", dims = 1:30, verbose = FALSE)
# integrated_plus <- FindClusters(integrated_plus, resolution = c(0.2, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4), verbose = FALSE)

saveRDS(integrated, file.path(out_dir, "scrna_integrated_enrich.rds"), compress = "xz")

# Minimal audit export (as you prefer)
md <- integrated_plus@meta.data
emb <- Embeddings(integrated_plus, "umap")
md$umap_1 <- emb[,1]; md$umap_2 <- emb[,2]
write.table(md, file.path(out_dir, "snrna_integrated_plus_metadata.tsv"),
            sep="\t", quote=FALSE, row.names=TRUE)


integrated <- PrepSCTFindMarkers(integrated, assay = "SCT")
# 1. ❗ IMPORTANT: Before running FindAllMarkers, ensure the SCT assay is properly prepared.
# 2. ✅ Now, run FindAllMarkers again on the SCT assay.
# The function will now use the properly prepared data.
all_markers <- FindAllMarkers(
  object = integrated, 
  assay = "SCT", 
  only.pos = TRUE, 
  min.pct = 0.25, 
  logfc.threshold = 0.25,
  test.use = "wilcox" 
)

#################################################################################
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
plan("multicore", workers = 20)
options(future.globals.maxSize = 400 * 1024^3) # 450 GB
options(Seurat.object.assay.version = "v5")
set.seed(1234) # for reproducibility

# == Paths ==================================================================
in_dir  <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/rds_adult"
out_dir <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult"
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

# Split the RNA assay by the 'orig.ident' which now corresponds to your layers
# merged_obj[["RNA"]] <- split(merged_obj[["RNA"]], f = merged_obj$orig.ident)
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
                          vars.to.regress = c("S.Score", "G2M.Score"),
                          verbose = FALSE)

message("Step 2.1: Save file for Backup...")
# Save intermediate object for backup
saveRDS(merged_obj, file.path(out_dir, "merge_backup_SCT_2samp.rds"))

# For safety, reload the backup
merged_obj <- readRDS('/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/merge_backup_SCT.rds')

# merged_obj <- readRDS('/mnt/18T/chibao/gliomas/data/upstream/sc_sn/official/integrated_v5_optimized/merge_backup_SCT.rds')
merged_obj
# Run PCA on the SCT assay. Seurat will automatically perform this for each layer.
message("Running PCA and UMAP on each layer...")
merged_obj <- RunPCA(merged_obj, assay = "SCT", verbose = FALSE)
merged_obj <- RunUMAP(merged_obj, dims = 1:30, verbose = FALSE)

harmony_obj <- merged_obj
scvi_obj <- merged_obj
message("Step 2.2: Integrating layers using Harmony...")
# backup <- readRDS('/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/merge_backup.rds')
# backup <- RunPCA(backup, assay = "SCT", verbose = FALSE)
# backup <- RunUMAP(backup, dims = 1:30, verbose = FALSE)
harmony_obj <- IntegrateLayers(
  object = harmony_obj,
  method = HarmonyIntegration,
  normalization.method = "SCT",
  # reference = reference_datasets,
  orig.reduction = "pca",
  #k.weight = 50,
  new.reduction = "harmony", # Name of the new integrated reduction
  dims = 1:30, # Using more PCs can be beneficial for complex datasets
  verbose = TRUE
)

# == 3. Integration using Seurat v5 `IntegrateLayers` ========================
message("Step 3: Integrating layers using reference-based RPCA...")
# RPCA Integration
merged_obj <- IntegrateLayers(
  object = merged_obj,
  method = RPCAIntegration,
  normalization.method = "SCT",
  # reference = reference_datasets,
  orig.reduction = "pca",
  k.weight = 50,
  new.reduction = "integrated.rpca", # Name of the new integrated reduction
  dims = 1:30, # Using more PCs can be beneficial for complex datasets
  verbose = TRUE
)

# Deep integration using scVI
Sys.setenv(
  RETICULATE_PYTHON = "/opt/miniforge/envs/scvi/bin/python",
  # Hide all GPUs from CUDA & Lightning (must be set before Python initializes!)
  CUDA_VISIBLE_DEVICES = "",
  PL_ACCELERATOR = "cpu",
  PL_DEVICES = "1"
)
library(reticulate)
py_config()
reticulate::py_run_string("
import torch
print('CUDA available:', torch.cuda.is_available())
")
# Expect: CUDA available: False

DefaultAssay(scvi_obj) <- "RNA"
features_to_integrate <- VariableFeatures(scvi_obj, assay = "SCT")
scvi_obj <- IntegrateLayers(
  object = scvi_obj,
  method = scVIIntegration,
  new.reduction = "integrated.scvi",
  features = features_to_integrate,  # <-- This is the main fix
  conda_env = "/opt/miniforge/envs/scvi",
  verbose = TRUE
)

# Harmony 

# Save
saveRDS(scvi_obj, file.path('/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized', "sc_integrated_scvi.rds"))


# == 4. Downstream Analysis and Clustering ===================================
message("Step 4: UMAP, Neighbors, and Clustering on integrated data...")
# Note: All downstream steps now use the new 'integrated.rpca' reduction
# Harmony
harmony_obj <- RunPCA(harmony_obj, assay = "SCT", npcs = 100, verbose = FALSE)
harmony_obj <- RunUMAP(harmony_obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
harmony_obj <- FindNeighbors(harmony_obj, reduction = "harmony", dims = 1:30)
harmony_obj <- FindClusters(
  harmony_obj,
  #graph.name = "SCT_snn",
  resolution = c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.4, 0.5, 0.6,  0.7, 0.8, 1.0, 1.2),
  algorithm = 1,     # Louvain (stable), switch to 2 (SLM) if desired
  verbose = FALSE
)


# Batchelor
library(batchelor)
set.seed(1101)
# Convert to SCE
DefaultAssay(merged_obj) <- "SCT"
# as.SingleCellExperiment will put RNA@counts into `counts` (if present)
# and the current assay's data (SCT) into `logcounts`.
sce <- as.SingleCellExperiment(merged_obj, assay = "SCT")
assayNames(sce)      # should include "counts" (raw) and "logcounts" (SCT)
table(sce$orig.ident)
# Split by batch
sce_norm <- batchelor::multiBatchNorm(sce, batch = sce$orig.ident)
sce_norm <- correctExperiments(sce_norm, batch=sce_norm$orig.ident, PARAM=FastMnnParam())
seu <- sce_norm
library(scater)
sce_norm <- runUMAP(sce_norm, dimred="corrected")
library(bluster)
colLabels(sce_norm) <- clusterRows(reducedDim(sce_norm, "corrected"), NNGraphParam())

# 1. Specify the filename and open the PDF device
pdf(file = "/mnt/18T/chibao/gliomas/data/upstream/snRNA/official/integrated_v5_optimized/adult/my_combined_plot.pdf", width = 10, height = 5) 

# 2. Re-run your plotting command, but print the result directly
p <- gridExtra::grid.arrange(
    plotUMAP(sce_norm, colour_by="label", text_by="label"),
    plotUMAP(sce_norm, colour_by="batch"),
    ncol=2
)

# Explicitly print the gtable object 'p' to the open device
grid::grid.draw(p)

# 3. Close the device to finalize the file
dev.off()

ggsave(filename = plot_dir, plot = p1)
# RPCA
merged_obj <- RunUMAP(merged_obj, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
# merged_obj <- RunTSNE(merged_obj, 
#                       reduction = "integrated.rpca", 
#                       dims = 1:30, 
#                       reduction.name = "tsne.rpca",
#                       verbose = FALSE)
merged_obj <- FindNeighbors(merged_obj, reduction = "integrated.rpca", dims = 1:30)
merged_obj <- FindClusters(
  merged_obj,
  #graph.name = "SCT_snn",
  resolution = c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.4, 0.5, 0.6,  0.7, 0.8, 1.0, 1.2),
  algorithm = 1,     # Louvain (stable), switch to 2 (SLM) if desired
  verbose = FALSE
)

# scVI
# scvi_obj <- RunPCA(scvi_obj, assay = "SCT", npcs = 100, verbose = FALSE)
# scvi_obj <- RunUMAP(scvi_obj, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.scvi")

scvi_obj <- readRDS('/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/sc_integrated_scvi.rds')
scvi_obj
DefaultAssay(scvi_obj) <- "SCT"
scvi_obj <- FindNeighbors(
  scvi_obj,
  reduction   = "integrated.scvi",
  dims        = 1:30,
  k.param     = 100,
  nn.method   = "rann",
  graph.name  = c("scvi_nn", "scvi_snn")   # <-- first = kNN, second = SNN
)


scvi_obj <- FindClusters(
  scvi_obj,
  graph.name = "scvi_snn",
  resolution = c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.4, 0.5, 0.6,  0.7, 0.8, 1.0, 1.2),
  algorithm = 1,     # Louvain (stable), switch to 2 (SLM) if desired
  verbose = FALSE
)
scvi_obj <- RunUMAP(
  scvi_obj,
  reduction      = "integrated.scvi",
  dims           = 1:30,
  nn.name        = "scvi_nn",
  reduction.name = "umap.scvi",
  verbose        = FALSE
)


scvi_obj <- FindNeighbors(
  scvi_obj,
  reduction = "integrated.scvi",
  dims      = 1:30,
  k.param   = 100,
  nn.method = "annoy",         # optional: for very large objects; needs RcppAnnoy installed
  annoy.metric = "cosine"     # optional: if you prefer cosine in neighbor graph
  #graph.name = sprintf("SNN_scvi_k%d_d%d", best_k, latent_k)  # keeps graphs identifiable
)

scvi_obj <- FindClusters(
  scvi_obj,
  #graph.name = "SCT_snn",
  resolution = c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.2, 0.4, 0.5, 0.6,  0.7, 0.8, 1.0, 1.2),
  algorithm = 1,     # Louvain (stable), switch to 2 (SLM) if desired
  verbose = FALSE
)

# 1) Pull the scVI latent (30D unless you re-integrated with ndims)
latent_name <- "integrated.scvi"
X <- Embeddings(scvi_obj, latent_name)
stopifnot(is.matrix(X) && is.finite(X[1,1]))

# 2) Set your tuned params (you can tweak these)
umap_metric     <- "cosine"
umap_neighbors  <- 200
umap_min_dist   <- 0.50
umap_epochs     <- 500
annoy_trees     <- 75
annoy_search_k  <- 6000
set.seed(1234)

# 3) Run uwot::umap directly with PCA init (bypasses spectral init)
#    NOTE: n_threads controls both Annoy search and SGD threads here.
U <- uwot::umap(
  X,
  n_neighbors   = umap_neighbors,
  metric        = umap_metric,
  min_dist      = umap_min_dist,
  n_epochs      = umap_epochs,
  init          = "pca",        # <-- crucial to avoid spectral
  n_trees       = annoy_trees,  # Annoy index
  search_k      = annoy_search_k,
  n_threads     = parallel::detectCores(),  # or set an integer
  verbose       = TRUE
)

# 4) Attach as a Seurat reduction
umap_name <- sprintf("umap.scvi.nd%d.k%d.md%02d.%s",
                     ncol(X), umap_neighbors, as.integer(umap_min_dist*100), umap_metric)

scvi_obj[[umap_name]] <- CreateDimReducObject(
  embeddings = U,
  key       = "UMAP_",
  assay     = DefaultAssay(scvi_obj)
)

# 5) Quick checks
cat("UMAP dims:", ncol(Embeddings(scvi_obj, umap_name)), " (expect 2)\n")
cat("UMAP rows:", nrow(Embeddings(scvi_obj, umap_name)), " (expect ~", nrow(X), ")\n")

# == 5. VALIDATION AND VISUALIZATION ========================================
message("Step 5: Validating and visualizing integration results...")

# 5.1. Visual Inspection for Batch Correction
# A good integration shows a "salt-and-pepper" mix of datasets in each cluster.
# Bad integration shows clusters dominated by one dataset.
p1 <- DimPlot(merged_obj, reduction = "umap.rpca", group.by = "orig.ident", shuffle = TRUE, label.size = 2) +
      ggtitle("UMAP Colored by Dataset (Batch)")
p2 <- DimPlot(merged_obj, reduction = "umap.rpca", group.by = "rpca_clusters_res.0.8", label = TRUE, repel = TRUE) +
      ggtitle("UMAP Colored by Cluster (Resolution 0.8)")

# Save the validation plots
ggsave(file.path(out_dir, "umap_validation_plots.png"), plot = p1 + p2, width = 16, height = 7, dpi = 300)

# 5.2. Quantitative Assessment of Cluster Stability (Conceptual)
# At this stage, you should choose the best clustering resolution.
# A "good" resolution produces stable clusters that are not just artifacts of parameter choice.
# The 'scclusteval' package is excellent for this. The workflow involves:
#   1. Subsampling your data (e.g., 80% of cells) repeatedly.
#   2. Re-running the clustering on each subsample.
#   3. Calculating the Jaccard similarity index between original and re-clustered cells.
#   4. Choosing a resolution that maximizes the number of stable clusters (Jaccard index > 0.75 is a good sign).
# This is a deeper analysis, and you can apply it using the 'scclusteval' GitHub tutorials.
message("VALIDATION NOTE: Review the UMAP plots. For quantitative analysis, consider using 'scclusteval' to find the most stable clustering resolution.")


# == 6. Marker Gene Identification ===========================================
message("Step 6: Finding marker genes for a chosen resolution...")
# After choosing a stable resolution (e.g., 0.8), set it as the main identity.
Idents(merged_obj) <- "SCT_snn_res.0.04"

# Prepare the SCT assay for DE analysis
# This is required when using multiple SCT models
merged_obj <- PrepSCTFindMarkers(merged_obj, assay = "SCT", verbose = TRUE)

all_markers <- FindAllMarkers(
  object = harmony_obj,
  assay = "SCT",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

cluster0 <- all_markers %>% filter(cluster == 0) %>% arrange(desc(avg_log2FC))
cluster1 <- all_markers %>% filter(cluster == 1) %>% arrange(desc(avg_log2FC))
cluster2 <- all_markers %>% filter(cluster == 2) %>% arrange(desc(avg_log2FC))
cluster3 <- all_markers %>% filter(cluster == 3) %>% arrange(desc(avg_log2FC))
cluster4 <- all_markers %>% filter(cluster == 4) %>% arrange(desc(avg_log2FC))
cluster5 <- all_markers %>% filter(cluster == 5) %>% arrange(desc(avg_log2FC))
cluster6 <- all_markers %>% filter(cluster == 6) %>% arrange(desc(avg_log2FC))
cluster7 <- all_markers %>% filter(cluster == 7) %>% arrange(desc(avg_log2FC))
cluster8 <- all_markers %>% filter(cluster == 8) %>% arrange(desc(avg_log2FC))
cluster9 <- all_markers %>% filter(cluster == 9) %>% arrange(desc(avg_log2FC))
cluster10 <- all_markers %>% filter(cluster == 10) %>% arrange(desc(avg_log2FC))
cluster11 <- all_markers %>% filter(cluster == 11) %>% arrange(desc(avg_log2FC))
cluster12 <- all_markers %>% filter(cluster == 12) %>% arrange(desc(avg_log2FC))
cluster13 <- all_markers %>% filter(cluster == 13) %>% arrange(desc(avg_log2FC))
cluster14 <- all_markers %>% filter(cluster == 14) %>% arrange(desc(avg_log2FC))
cluster15 <- all_markers %>% filter(cluster == 15) %>% arrange(desc(avg_log2FC))
cluster16 <- all_markers %>% filter(cluster == 16) %>% arrange(desc(avg_log2FC))
cluster17 <- all_markers %>% filter(cluster == 17) %>% arrange(desc(avg_log2FC))
cluster18 <- all_markers %>% filter(cluster == 18) %>% arrange(desc(avg_log2FC))
cluster19 <- all_markers %>% filter(cluster == 19) %>% arrange(desc(avg_log2FC))
cluster20 <- all_markers %>% filter(cluster == 20) %>% arrange(desc(avg_log2FC)) 
cluster21 <- all_markers %>% filter(cluster == 21) %>% arrange(desc(avg_log2FC))
cluster22 <- all_markers %>% filter(cluster == 22) %>% arrange(desc(avg_log2FC))
cluster23 <- all_markers %>% filter(cluster == 23) %>% arrange(desc(avg_log2FC))
cluster24 <- all_markers %>% filter(cluster == 24) %>% arrange(desc(avg_log2FC))
cluster25 <- all_markers %>% filter(cluster == 25) %>% arrange(desc(avg_log2FC))
cluster26 <- all_markers %>% filter(cluster == 26) %>% arrange(desc(avg_log2FC))
cluster27 <- all_markers %>% filter(cluster == 27) %>% arrange(desc(avg_log2FC))
cluster28 <- all_markers %>% filter(cluster == 28) %>% arrange(desc(avg_log2FC))

cluster0 |> head(10)
cluster1 |> head(10)
cluster2 |> head(10)
cluster3 |> head(10)
cluster4 |> head(10)
cluster5 |> head(10)
cluster6 |> head(10)
cluster7 |> head(10)
cluster8 |> head(10)
cluster9 |> head(10)
cluster10 |> head(10)
cluster11 |> head(10)
cluster12 |> head(10)
cluster13 |> head(10)
cluster14 |> head(10)
cluster15 |> head(10)
cluster16 |> head(10)
cluster17 |> head(10)
cluster18 |> head(10)
cluster19 |> head(10)
cluster20 |> head(10)
cluster21 |> head(10)
cluster22 |> head(10)
cluster23 |> head(10)
cluster24 |> head(10)
cluster25 |> head(10)
cluster26 |> head(10)
cluster27 |> head(10)
cluster28 |> head(10)

# Save results
saveRDS(merged_obj, file.path(out_dir, "scrna_integrated_harmony_final.rds"))
saveRDS(harmony_obj, file.path(out_dir, "scrna_integrated_cca.rds"))

write.csv(all_markers, file.path(out_dir, "all_cluster_markers_res0.8.csv"))

message("✅ Workflow Complete!")

###############################################################################################
#!/usr/bin/env Rscript

# == 0. SETUP ===============================================================
# This script integrates multiple single-cell datasets using the Seurat v5
# reference-based RPCA workflow. It is optimized for large, diverse projects
# and includes steps for robust QC, validation, and handling of small samples.

# -- Load Libraries ---------------------------------------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(future)
  library(glmGamPoi)
})

# -- Configure Environment --------------------------------------------------
plan("multicore", workers = 8) # Set to the number of cores you want to use
options(future.globals.maxSize = 80 * 1024^2) # 80 GB RAM limit for parallel processes
options(Seurat.object.assay.version = "v5")
set.seed(1234) # for reproducibility

# -- Define Paths -----------------------------------------------------------
in_dir  <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/set1_manifest/rds"
out_dir <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/set1_manifest/integrated_v5_optmi"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# == 1. LOAD & PREPARE DATA =================================================
# Load all individual Seurat objects and separate them into 'main' and 'tiny'
# based on cell count. This prevents small, potentially low-quality samples
# from negatively influencing the main integration.

message("Step 1: Loading and separating objects by size...")
rds_files <- list.files(in_dir, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
all_objs_list <- lapply(rds_files, readRDS)

cell_counts <- sapply(all_objs_list, ncol)
tiny_indices <- which(cell_counts < 200)
main_indices <- which(cell_counts >= 200)

main_objs_list <- all_objs_list[main_indices]
tiny_objs_list <- all_objs_list[tiny_indices]

message(paste("Found", length(main_objs_list), "main objects for integration."))
message(paste("Found", length(tiny_objs_list), "tiny objects (<200 cells) for mapping."))

# Clean up memory
rm(all_objs_list)
gc()


# == 2. MERGE & PRE-PROCESS MAIN OBJECTS ====================================
# The modern Seurat v5 workflow uses a single object with 'layers' for integration,
# which is more efficient than a list of objects.

message("Step 2: Merging main objects and running pre-processing...")
layer_names <- sapply(main_objs_list, function(x) unique(x$sample_uid))
merged_main_obj <- merge(x = main_objs_list[[1]], y = main_objs_list[2:length(main_objs_list)], add.cell.ids = layer_names)

# Split the RNA assay by sample; this creates the 'layers' for integration.
# merged_main_obj[["RNA"]] <- split(merged_main_obj[["RNA"]], f = merged_main_obj$orig.ident)

rm(main_objs_list)
gc()

# Score cell cycle phases. This is a critical step to mitigate a major source
# of technical variation, especially in cancer datasets.
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
merged_main_obj <- CellCycleScoring(merged_main_obj, s.features = s.genes, g2m.features = g2m.genes)


# == 3. SCT NORMALIZE MAIN OBJECTS ==========================================
message("Step 3: Running SCTransform on each layer of the main object...")
# SCTransform is run on the merged object and automatically processes each layer.
# We regress out cell cycle scores to improve downstream biological clustering.
merged_main_obj <- SCTransform(merged_main_obj,
                               method = "glmGamPoi",
                               vars.to.regress = c("S.Score", "G2M.Score"),
                               verbose = FALSE)


# == 4. INTEGRATE MAIN OBJECTS ==============================================
message("Step 4: Integrating layers using reference-based RPCA...")
# We use the streamlined `IntegrateLayers` function.
# Reference-based integration is faster and more scalable for large projects
#.

# --- IMPORTANT: CHOOSE YOUR REFERENCE DATASET(S) ---
# Select one or more high-quality datasets (e.g., high cell/feature counts)
# to serve as the integration reference. Here, we assume the first two are best.
# reference_datasets <- c(1, 2)
# if (length(layer_names) < 2) {
#     reference_datasets <- 1
# }

integrated_main <- IntegrateLayers(
  object = merged_main_obj,
  method = RPCAIntegration, # RPCA is fast and robust for similar datasets.
  normalization.method = "SCT",
  # reference = reference_datasets,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca", # Name for the new integrated reduction
  dims = 1:50, # Using more PCs can help capture complex biology.
  verbose = TRUE
)

# Clean up the pre-integration object to save memory
rm(merged_main_obj)
gc()


# == 5. CLUSTERING & VISUALIZATION OF MAIN OBJECT ===========================
message("Step 5: Clustering and UMAP visualization...")
# All downstream analysis on the integrated data uses the 'integrated.rpca' reduction.
integrated_main <- RunUMAP(integrated_main, reduction = "integrated.rpca", dims = 1:50, reduction.name = "umap.rpca", return.model = TRUE)
integrated_main <- FindNeighbors(integrated_main, reduction = "integrated.rpca", dims = 1:50)
integrated_main <- FindClusters(integrated_main, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2), cluster.name = "rpca_clusters")

# -- Visual Validation of Batch Correction --
# A good integration shows a "salt-and-pepper" mix of datasets in each cluster.
p1 <- DimPlot(integrated_main, reduction = "umap.rpca", group.by = "orig.ident", shuffle = TRUE, label.size = 2) +
      ggtitle("UMAP Colored by Dataset (Batch)")
p2 <- DimPlot(integrated_main, reduction = "umap.rpca", group.by = "rpca_clusters_res.0.8", label = TRUE, repel = TRUE) +
      ggtitle("UMAP Colored by Cluster (Resolution 0.8)")

ggsave(file.path(out_dir, "umap_main_integration_validation.png"), plot = p1 | p2, width = 16, height = 7, dpi = 300)

message("VALIDATION NOTE: Review UMAP plots. Consider using 'scclusteval' to find the most stable clustering resolution.")
# For a deeper analysis, choose a resolution that maximizes the number of stable clusters.


# == 6. MAP & MERGE TINY SAMPLES ============================================
message("Step 6: Mapping tiny objects onto the integrated reference...")
final_integrated_obj <- integrated_main

if (length(tiny_objs_list) > 0) {
  # Normalize the tiny objects using the same SCT method
  tiny_objs_list <- lapply(tiny_objs_list, FUN = function(x) {
    SCTransform(x, method = "glmGamPoi", vars.to.regress = c("S.Score", "G2M.Score"), verbose = FALSE)
  })

  # Find transfer anchors and map each tiny object to the main integrated reference
  mapped_tiny_list <- lapply(tiny_objs_list, FUN = function(query_obj) {
    anchors <- FindTransferAnchors(
      reference = integrated_main,
      query = query_obj,
      normalization.method = "SCT",
      reference.reduction = "integrated.rpca",
      dims = 1:50
    )
    MapQuery(
      anchorset = anchors,
      reference = integrated_main,
      query = query_obj,
      refdata = list(predicted.cluster = "rpca_clusters_res.0.8"), # Transfer cluster labels
      reference.reduction = "integrated.rpca",
      reduction.model = "umap.rpca" # Project onto the reference UMAP
    )
  })

  message("Merging mapped tiny objects into the final object...")
  final_integrated_obj <- merge(integrated_main, y = mapped_tiny_list)
} else {
  message("No tiny objects to map.")
}


# == 7. FIND MARKERS & SAVE RESULTS =========================================
message("Step 7: Finding marker genes...")
# After choosing a stable resolution (e.g., 0.8), set it as the main identity.
Idents(final_integrated_obj) <- "rpca_clusters_res.0.8"

# Prepare the SCT assay for DE analysis. This is required when using multiple
# SCT models from different batches/layers.
final_integrated_obj <- PrepSCTFindMarkers(final_integrated_obj, assay = "SCT", verbose = TRUE)

all_markers <- FindAllMarkers(
  object = final_integrated_obj,
  assay = "SCT",
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

# -- Save Final Outputs --
message("Saving final integrated object and marker list...")
saveRDS(final_integrated_obj, file.path(out_dir, "glioma_integrated_final_object.rds"))
write.csv(all_markers, file.path(out_dir, "all_cluster_markers_res0.8.csv"), row.names = FALSE)

message("✅ Workflow Complete!")