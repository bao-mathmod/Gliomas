#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(stringr); library(future)
})
plan("multicore", workers = 8)
options(future.globals.maxSize = 133120 * 1024^2) # Sets limit to 130 GB
library(glmGamPoi)
# Ensure the object convert to V5
options(Seurat.object.assay.version = "v5") 

# in_dir  <- "/mnt/18T/chibao/gliomas/data/upstream/snRNA/snRNA_clean_2/rds"
# out_dir <- "/mnt/18T/chibao/gliomas/data/upstream/snRNA/snRNA_clean_2/integrated"; dir.create(out_dir, FALSE, TRUE)

in_dir  <- "/mnt/18T/chibao/gliomas/data/upstream/snRNA/snRNA_clean_2/rds_hgnc_test"
out_dir <- "/mnt/18T/chibao/gliomas/data/upstream/snRNA/snRNA_clean_2/integrated_test"; dir.create(out_dir, FALSE, TRUE)

objs <- lapply(list.files(in_dir, pattern="\\.rds$", full.names=TRUE), readRDS)

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

# Rerun SCT after HGNC harmonization to ensure consistency
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
tiny_idx <- c(5, 34, 35) # change as needed
objs_main <- objs[-tiny_idx]
objs_tiny <- objs[tiny_idx]

# RPCA integration with SCT normalization
# 1) Feature selection + Prep (SCT path)
# 2. Integrate only the main cohort
features <- SelectIntegrationFeatures(object.list = objs_main, nfeatures = 3000)
objs_main <- PrepSCTIntegration(object.list = objs_main, anchor.features = features)


# 2) IMPORTANT in RPCA: RunPCA on EACH object using the same features (v5 vignette)
objs_main <- lapply(objs_main, function(x) RunPCA(x, features = features))

# 3) Anchors (RPCA + SCT). Consider k.anchor tweak & reference-based integration for big cohorts
anchors <- FindIntegrationAnchors(object.list = objs_main,
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
saveRDS(integrated, file.path(out_dir, "snrna_integrated.rds"), compress = "xz")

## 4.1) Map each tiny object separately
## Pre-req (reference):
## integrated already has PCA; if you want UMAP projection, run with return.model=TRUE earlier:
## integrated <- RunUMAP(integrated, reduction="pca", dims=1:30, return.model=TRUE)

map_one <- function(q, ref, dims = 1:20, min_feats = 500) {
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
mapped_list <- lapply(objs_tiny, map_one, ref = integrated, dims = 1:20)

# 4.2) Merge mapped tiny objects into integrated
integrated_plus <- merge(integrated, y = mapped_list)

integrated_plus <- FindNeighbors(integrated_plus, reduction = "pca", dims = 1:20, verbose = FALSE)
integrated_plus <- FindClusters(integrated_plus, resolution = c(0.2, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4), verbose = FALSE)

saveRDS(integrated_plus, file.path(out_dir, "snrna_integrated_plus.rds"), compress = "xz")

# Minimal audit export (as you prefer)
md <- integrated_plus@meta.data
emb <- Embeddings(integrated_plus, "umap")
md$umap_1 <- emb[,1]; md$umap_2 <- emb[,2]
write.table(md, file.path(out_dir, "snrna_integrated_plus_metadata.tsv"),
            sep="\t", quote=FALSE, row.names=TRUE)
