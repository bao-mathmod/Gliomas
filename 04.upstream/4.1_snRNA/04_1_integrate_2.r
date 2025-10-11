#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(stringr); library(future)
  library(glmGamPoi)
})
plan("multicore", workers = 8)
options(future.globals.maxSize = 133120 * 1024^2) # ~130 GB
options(Seurat.object.assay.version = "v5")
library(digest)

in_dir  <- "/mnt/18T/chibao/gliomas/data/upstream/snRNA/snRNA_clean_2/rds_hgnc_test"
out_dir <- "/mnt/18T/chibao/gliomas/data/upstream/snRNA/snRNA_clean_2/integrated_test"; dir.create(out_dir, FALSE, TRUE)

objs <- lapply(list.files(in_dir, pattern="\\.rds$", full.names=TRUE), readRDS)

# 0) Update & sanity for RNA assay version
objs <- lapply(objs, UpdateSeuratObject)

# 1) Create a stable per-object ID to prefix cell names
make_id <- function(o, fallback_prefix = "sample") {
  # priority: sample_id field -> project.name -> basename of file not available here
  sid <- NA
  if ("sample_id" %in% colnames(o@meta.data)) sid <- unique(as.character(o$sample_id))
  sid <- sid[!is.na(sid) & sid != ""]
  if (length(sid) == 1) return(sid)
  # else fallback to project name + a random short token
  paste0(fallback_prefix, "_", gsub("[^A-Za-z0-9]+","", o@project.name), "_", substr(digest::digest(colnames(o)[1]), 1, 6))
}

ids <- vapply(objs, make_id, character(1))

# 2) Enforce unique cell names per object and across all objects
for (i in seq_along(objs)) {
  new_names <- paste0(ids[i], "_", colnames(objs[[i]]))
  # If already prefixed, avoid double prefixing
  if (!all(startsWith(colnames(objs[[i]]), paste0(ids[i], "_")))) {
    objs[[i]] <- RenameCells(objs[[i]], new.names = new_names)
  }
  # Align meta.data rownames to cell names
  md <- objs[[i]]@meta.data
  if (!identical(rownames(md), colnames(objs[[i]]))) {
    md <- md[match(colnames(objs[[i]]), rownames(md)), , drop=FALSE]
    rownames(md) <- colnames(objs[[i]])
    objs[[i]]@meta.data <- md
  }
}

stopifnot(!any(duplicated(unlist(lapply(objs, colnames)))))  # hard fail if still duplicates

# 3) Fresh SCT on all objects (post-HGNC), no var-gene-only restriction
objs <- lapply(objs, function(x) {
  DefaultAssay(x) <- "RNA"
  x <- SCTransform(x, vst.flavor = "v2", method = "glmGamPoi",
                   return.only.var.genes = FALSE, verbose = FALSE)
  DefaultAssay(x) <- "SCT"
  x
})

# 4) Optionally drop tiny objects from anchor-building (auto-detect)
sizes <- vapply(objs, ncol, integer(1))
tiny_thresh <- 200L
tiny_idx <- which(sizes < tiny_thresh)
objs_main <- if (length(tiny_idx)) objs[-tiny_idx] else objs
objs_tiny <- if (length(tiny_idx)) objs[tiny_idx] else list()

if (length(objs_main) < 2) stop("Not enough main objects for integration after excluding tiny sets.")

# 5) Select integration features (SCT) and remove MT/Ribo
features <- SelectIntegrationFeatures(object.list = objs_main, nfeatures = 4000)
features <- features[!grepl("^MT-", features) & !grepl("^RP[SL]", features)]
features <- head(unique(features), 3000L)

# 6) Ensure VariableFeatures match
objs_main <- lapply(objs_main, function(x) { VariableFeatures(x) <- features; x })

# 7) Wipe stale scale.data layers if features don’t match
objs_main <- lapply(objs_main, function(x) {
  sd <- tryCatch(LayerData(x, assay = "SCT", layer = "scale.data"), error = function(e) NULL)
  if (!is.null(sd)) {
    if (!identical(rownames(sd), features)) {
      LayerData(x, assay = "SCT", layer = "scale.data") <- NULL
    }
  }
  x
})

# 8) Prep + RPCA per object
objs_main <- PrepSCTIntegration(object.list = objs_main, anchor.features = features)
objs_main <- lapply(objs_main, function(x) RunPCA(x, features = features, verbose = FALSE))

# 9) Anchors & Integrate
anchors <- FindIntegrationAnchors(object.list = objs_main,
                                  normalization.method = "SCT",
                                  reduction = "rpca",
                                  anchor.features = features,
                                  dims = 1:30)

integrated <- IntegrateData(anchorset = anchors,
                            normalization.method = "SCT",
                            dims = 1:30)

# 10) (Optional) Map tiny objects onto the integrated reference
if (length(objs_tiny)) {
  # MapQuery in v5 world uses reference reductions
  for (i in seq_along(objs_tiny)) {
    objs_tiny[[i]] <- SCTransform(objs_tiny[[i]], vst.flavor="v2", method="glmGamPoi",
                                  return.only.var.genes = FALSE, verbose = FALSE)
  }
  ref <- integrated
  ref <- RunPCA(ref, npcs = 50, verbose = FALSE)
  ref <- RunUMAP(ref, dims = 1:30, verbose = FALSE)
  # Map each tiny object (simplified; adapt if you need transfer of labels)
  mapped <- lapply(objs_tiny, function(q) {
    q <- PrepSCTIntegration(object.list = list(q), anchor.features = features)[[1]]
    q  <- RunPCA(q, features = features, verbose = FALSE)
    anchors_q <- FindTransferAnchors(reference = ref, query = q,
                                     normalization.method = "SCT",
                                     reduction = "rpca", dims = 1:30,
                                     reference.reduction = "pca")
    MapQuery(anchorset = anchors_q, reference = ref, query = q,
             refdata = NULL, # add label transfer here if you have
             reference.reduction = "pca", reduction.model = "umap")
  })
  # Merge mapped with integrated if desired (optional)
}

saveRDS(integrated, file.path(out_dir, "integrated_sct_rpca.rds"), compress = "xz")
cat("Integrated object saved to:", file.path(out_dir, "integrated_sct_rpca.rds"), "\n")
