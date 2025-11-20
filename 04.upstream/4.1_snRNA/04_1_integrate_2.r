# #!/usr/bin/env Rscript
# suppressPackageStartupMessages({
#   library(Seurat); library(dplyr); library(stringr); library(future)
#   library(glmGamPoi)
# })
# plan("multicore", workers = 8)
# options(future.globals.maxSize = 133120 * 1024^2) # ~130 GB
# options(Seurat.object.assay.version = "v5")
# library(digest)

# in_dir  <- "/mnt/18T/chibao/gliomas/data/upstream/snRNA/snRNA_clean_2/rds_hgnc_test"
# out_dir <- "/mnt/18T/chibao/gliomas/data/upstream/snRNA/snRNA_clean_2/integrated_test"; dir.create(out_dir, FALSE, TRUE)

# objs <- lapply(list.files(in_dir, pattern="\\.rds$", full.names=TRUE), readRDS)

# # 0) Update & sanity for RNA assay version
# objs <- lapply(objs, UpdateSeuratObject)

# # 1) Create a stable per-object ID to prefix cell names
# make_id <- function(o, fallback_prefix = "sample") {
#   # priority: sample_id field -> project.name -> basename of file not available here
#   sid <- NA
#   if ("sample_id" %in% colnames(o@meta.data)) sid <- unique(as.character(o$sample_id))
#   sid <- sid[!is.na(sid) & sid != ""]
#   if (length(sid) == 1) return(sid)
#   # else fallback to project name + a random short token
#   paste0(fallback_prefix, "_", gsub("[^A-Za-z0-9]+","", o@project.name), "_", substr(digest::digest(colnames(o)[1]), 1, 6))
# }

# ids <- vapply(objs, make_id, character(1))

# # 2) Enforce unique cell names per object and across all objects
# for (i in seq_along(objs)) {
#   new_names <- paste0(ids[i], "_", colnames(objs[[i]]))
#   # If already prefixed, avoid double prefixing
#   if (!all(startsWith(colnames(objs[[i]]), paste0(ids[i], "_")))) {
#     objs[[i]] <- RenameCells(objs[[i]], new.names = new_names)
#   }
#   # Align meta.data rownames to cell names
#   md <- objs[[i]]@meta.data
#   if (!identical(rownames(md), colnames(objs[[i]]))) {
#     md <- md[match(colnames(objs[[i]]), rownames(md)), , drop=FALSE]
#     rownames(md) <- colnames(objs[[i]])
#     objs[[i]]@meta.data <- md
#   }
# }

# stopifnot(!any(duplicated(unlist(lapply(objs, colnames)))))  # hard fail if still duplicates

# # 3) Fresh SCT on all objects (post-HGNC), no var-gene-only restriction
# objs <- lapply(objs, function(x) {
#   DefaultAssay(x) <- "RNA"
#   x <- SCTransform(x, vst.flavor = "v2", method = "glmGamPoi",
#                    return.only.var.genes = FALSE, verbose = FALSE)
#   DefaultAssay(x) <- "SCT"
#   x
# })

# # 4) Optionally drop tiny objects from anchor-building (auto-detect)
# sizes <- vapply(objs, ncol, integer(1))
# tiny_thresh <- 200L
# tiny_idx <- which(sizes < tiny_thresh)
# objs_main <- if (length(tiny_idx)) objs[-tiny_idx] else objs
# objs_tiny <- if (length(tiny_idx)) objs[tiny_idx] else list()

# if (length(objs_main) < 2) stop("Not enough main objects for integration after excluding tiny sets.")

# # 5) Select integration features (SCT) and remove MT/Ribo
# features <- SelectIntegrationFeatures(object.list = objs_main, nfeatures = 4000)
# features <- features[!grepl("^MT-", features) & !grepl("^RP[SL]", features)]
# features <- head(unique(features), 3000L)

# # 6) Ensure VariableFeatures match
# objs_main <- lapply(objs_main, function(x) { VariableFeatures(x) <- features; x })

# # 7) Wipe stale scale.data layers if features don’t match
# objs_main <- lapply(objs_main, function(x) {
#   sd <- tryCatch(LayerData(x, assay = "SCT", layer = "scale.data"), error = function(e) NULL)
#   if (!is.null(sd)) {
#     if (!identical(rownames(sd), features)) {
#       LayerData(x, assay = "SCT", layer = "scale.data") <- NULL
#     }
#   }
#   x
# })

# # 8) Prep + RPCA per object
# objs_main <- PrepSCTIntegration(object.list = objs_main, anchor.features = features)
# objs_main <- lapply(objs_main, function(x) RunPCA(x, features = features, verbose = FALSE))

# # 9) Anchors & Integrate
# anchors <- FindIntegrationAnchors(object.list = objs_main,
#                                   normalization.method = "SCT",
#                                   reduction = "rpca",
#                                   anchor.features = features,
#                                   dims = 1:30)

# integrated <- IntegrateData(anchorset = anchors,
#                             normalization.method = "SCT",
#                             dims = 1:30)

# # 10) (Optional) Map tiny objects onto the integrated reference
# if (length(objs_tiny)) {
#   # MapQuery in v5 world uses reference reductions
#   for (i in seq_along(objs_tiny)) {
#     objs_tiny[[i]] <- SCTransform(objs_tiny[[i]], vst.flavor="v2", method="glmGamPoi",
#                                   return.only.var.genes = FALSE, verbose = FALSE)
#   }
#   ref <- integrated
#   ref <- RunPCA(ref, npcs = 50, verbose = FALSE)
#   ref <- RunUMAP(ref, dims = 1:30, verbose = FALSE)
#   # Map each tiny object (simplified; adapt if you need transfer of labels)
#   mapped <- lapply(objs_tiny, function(q) {
#     q <- PrepSCTIntegration(object.list = list(q), anchor.features = features)[[1]]
#     q  <- RunPCA(q, features = features, verbose = FALSE)
#     anchors_q <- FindTransferAnchors(reference = ref, query = q,
#                                      normalization.method = "SCT",
#                                      reduction = "rpca", dims = 1:30,
#                                      reference.reduction = "pca")
#     MapQuery(anchorset = anchors_q, reference = ref, query = q,
#              refdata = NULL, # add label transfer here if you have
#              reference.reduction = "pca", reduction.model = "umap")
#   })
#   # Merge mapped with integrated if desired (optional)
# }

# saveRDS(integrated, file.path(out_dir, "integrated_sct_rpca.rds"), compress = "xz")
# cat("Integrated object saved to:", file.path(out_dir, "integrated_sct_rpca.rds"), "\n")

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(future)
  library(presto)
  library(harmony)
  library(glmGamPoi)
  library(SeuratData)
  library(SeuratWrappers)
  library(Azimuth)
})

# ========================= CONFIG ==========================================
plan("multicore", workers = 20)
options(future.globals.maxSize = 400 * 1024^3)  # 400 GB
options(Seurat.object.assay.version = "v5")
set.seed(1234)

# ---- Paths: ADAPT TO YOUR snRNA LAYOUT ----
in_dir  <- "/mnt/18T/chibao/gliomas/data/upstream/snRNA/dry_run/rds"  # <-- snRNA RDS from ingest/QC
out_dir <- "/mnt/18T/chibao/gliomas/data/upstream/snRNA/dry_run/integrated_v5_optimized"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

backup_rds <- file.path(out_dir, "merge_backup_SCT_snrna.rds")
final_rds  <- file.path(out_dir, "snrna_harmony_integrated_dry_run.rds")

message("=== snRNA Harmony Integration (dry-run for PRJNA1081384 + PRJNA1155637) ===")
message("Input dir : ", in_dir)
message("Output dir: ", out_dir)

# ========================= 1. LOAD RDS OBJECTS ==============================
message("Step 1: Loading individual snRNA Seurat objects...")

rds_files <- list.files(in_dir, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
if (length(rds_files) == 0) {
  stop("No .rds files found in in_dir: ", in_dir)
}
message("Found ", length(rds_files), " RDS files.")

objs_list <- lapply(rds_files, readRDS)

# Basic sanity: check for sample_uid metadata
if (!all(vapply(objs_list, function(x) "sample_uid" %in% colnames(x@meta.data), logical(1)))) {
  stop("Not all objects have a 'sample_uid' column in meta.data. Please ensure ingest/QC step added it.")
}

# Extract sample_uids per object (should be unique per RDS)
sample_uids <- vapply(objs_list, function(x) {
  su <- unique(x$sample_uid)
  if (length(su) != 1) {
    stop("Object has multiple sample_uid values. Expected exactly one per RDS. sample_uid: ", paste(su, collapse = ","))
  }
  su
}, character(1))

message("Sample UIDs detected: ", paste(sample_uids, collapse = ", "))

# ========================= 2. MERGE + LAYERS ================================
message("Step 2: Merging objects into a single Seurat object with layers (per sample_uid)...")

# Use sample_uid as add.cell.ids for clarity
merged_obj <- merge(
  x  = objs_list[[1]],
  y  = objs_list[2:length(objs_list)],
  add.cell.ids = sample_uids
)

merged_obj

rm(objs_list); gc()

DefaultAssay(merged_obj) <- "RNA"

# Join layers so the RNA assay is layered (Seurat v5 style)
merged_obj <- JoinLayers(merged_obj)
merged_obj
# JoinLayers expects each layer to correspond to a grouping of cells. Here we want per-sample layers.
backup_obj <- merged_obj
# We use 'sample_uid' instead of orig.ident, because orig.ident is project-level.
message("Step 2.1: Splitting RNA assay into layers by sample_uid/orig.ident...")
merged_obj[["RNA"]] <- split(merged_obj[["RNA"]], f = merged_obj$orig.ident) # use orig.ident or sample_uid here
merged_obj
DefaultAssay(merged_obj) <- "RNA"

message("Merged object summary:")
print(merged_obj)

# Optional sanity: number of cells per sample_uid
cell_counts <- table(merged_obj$sample_uid)
message("Cells per sample_uid:")
print(cell_counts)

# ========================= 3. CELL CYCLE SCORING (snRNA) ====================
message("Step 3: Cell cycle scoring on snRNA (RNA assay)...")
# You can disable this entire block if you decide not to regress cell cycle for nuclei.

s.genes  <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Normalize RNA for CC scoring only; SCT will use raw counts underneath
merged_obj <- NormalizeData(merged_obj, verbose = FALSE)
merged_obj <- CellCycleScoring(merged_obj, s.features = s.genes, g2m.features = g2m.genes)

# You now have S.Score and G2M.Score in meta.data

# ========================= 4. SCT NORMALIZATION (snRNA) =====================
message("Step 4: Running SCTransform (glmGamPoi) on snRNA per layer...")
message("  • No percent.mt regression (snRNA has very low mt%)")
message("  • Regressing cell cycle scores (S.Score, G2M.Score) – optional, but useful for tumor data")

# If you want to *not* regress cell cycle, set vars.to.regress = NULL
merged_obj <- SCTransform(
  merged_obj,
  method         = "glmGamPoi",
  vars.to.regress = c("S.Score", "G2M.Score"),
  verbose        = FALSE
)

message("Step 4.1: Saving backup post-SCT...")
saveRDS(merged_obj, backup_rds)
message("Backup saved to: ", backup_rds)

# If wanted, you can re-load (e.g. to confirm reproducibility or after crash)
merged_obj <- readRDS(backup_rds)

# ========================= 5. PCA ON SCT (snRNA) ============================
message("Step 5: Running PCA on SCT (snRNA, 40 PCs)...")
DefaultAssay(merged_obj) <- "SCT"

merged_obj <- RunPCA(merged_obj, npcs = 40, verbose = FALSE)

message("PCA complete. Variance explained by first few PCs:")
print(merged_obj[["pca"]]@stdev[1:10]^2 / sum(merged_obj[["pca"]]@stdev^2))

# Initial UMAP purely for QC (pre-integration)
message("Step 5.1: Running pre-integration UMAP (dims 1:40) for QC...")
merged_obj <- RunUMAP(merged_obj, reduction = "pca", dims = 1:40, reduction.name = "umap", verbose = FALSE)

# ========================= 6. HARMONY INTEGRATION ===========================
message("Step 6: Harmony integration across snRNA samples...")

harmony_obj <- merged_obj

# Determine which grouping variables to use in Harmony
# Minimum: sample_uid; optionally also project_id and processing (cellbender vs cellranger)
group_vars <- c("sample_uid")
if ("orig.ident" %in% colnames(harmony_obj@meta.data)) {
  group_vars <- c(group_vars, "orig.ident")
}
if ("processing" %in% colnames(harmony_obj@meta.data)) {
  group_vars <- c(group_vars, "processing")
}

message("Harmony will use group.by.vars = ", paste(group_vars, collapse = ", "))

harmony_obj <- IntegrateLayers(
  object               = harmony_obj,
  method               = HarmonyIntegration,
  normalization.method = "SCT",
  orig.reduction       = "pca",
  new.reduction        = "harmony",
  dims                 = 1:40,
  verbose              = TRUE
  #group.by.vars        = group_vars # Test result (v1: with this, v2: without)
)

# After Harmony, recompute PCA/UMAP/Neighbors/Clusters based on integrated embedding
message("Step 6.1: Running PCA (100 PCs) on SCT and UMAP on Harmony dims 1:40...")
harmony_obj <- RunPCA(harmony_obj, assay = "SCT", npcs = 100, verbose = FALSE)
harmony_obj <- RunUMAP(
  harmony_obj,
  reduction      = "harmony",
  dims           = 1:40,
  reduction.name = "umap.harmony",
  verbose        = FALSE
)

message("Step 6.2: Finding neighbors and clusters on Harmony embedding...")
harmony_obj <- FindNeighbors(harmony_obj, reduction = "harmony", dims = 1:40)

harmony_obj <- FindClusters(
  harmony_obj,
  resolution = c(0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09,
                 0.1, 0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.2),
  algorithm  = 1,
  verbose    = FALSE
)

plot_dir <- file.path('/mnt/18T/chibao/gliomas/data/upstream/snRNA/dry_run/integrated_v5_optimized', "harmony_2.png")
p <- DimPlot(harmony_obj, reduction = 'umap.harmony', group.by = 'SCT_snn_res.0.02', label = TRUE)
# library(ggplot2)
ggsave(filename = plot_dir, plot = p)

message("Clustering complete. Example: table at resolution 0.2 (SCT_snn_res.0.2):")
if ("SCT_snn_res.0.2" %in% colnames(harmony_obj@meta.data)) {
  print(table(harmony_obj$SCT_snn_res.0.2))
} else {
  message("Resolution 0.2 not found in meta.data (naming may differ by Seurat version).")
}

# ========================= 7. SAVE FINAL OBJECT ============================
message("Step 7: Saving final snRNA Harmony-integrated object...")
saveRDS(harmony_obj, final_rds)
message("Final integrated object saved to: ", final_rds)

message("=== snRNA Harmony Integration complete ===")
