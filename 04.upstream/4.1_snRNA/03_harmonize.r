# Different authors mix Ensembl IDs (±version) and gene symbols. Standardize row names to HGNC symbols (drop Ensembl versions), then deduplicate.

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(stringr)
  library(AnnotationDbi); library(org.Hs.eg.db)
})

in_dir  <- "/mnt/18T/chibao/gliomas/data/upstream/snRNA/snRNA_clean_2/rds"
out_dir <- "/mnt/18T/chibao/gliomas/data/upstream/snRNA/snRNA_clean_2/rds_hgnc_2"; dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

map_to_symbol <- function(ids){
  # If already looks like symbols (A-Z, -, digits), leave as-is
  looks_ens <- grepl("^ENSG\\d+", ids)
  clean_ids <- sub("\\.\\d+$", "", ids)        # drop .version
  symbols <- ifelse(looks_ens,
                    mapIds(org.Hs.eg.db, keys=clean_ids, keytype="ENSEMBL", column="SYMBOL"),
                    ids)
  symbols <- ifelse(is.na(symbols) | symbols=="", ids, symbols)
  make.unique(symbols)
}

rds <- list.files(in_dir, pattern="\\.rds$", full.names=TRUE)
for (f in rds) {
  obj <- readRDS(f)
  sym <- map_to_symbol(rownames(obj))
  # collapse duplicates by symbol (sum)
  if (any(duplicated(sym))) {
    m <- LayerData(obj, assay="SCT", slot="counts")
    # If SCT not present yet (depending on your earlier script), fall back to RNA
    if (nrow(m)==0) m <- LayerData(obj, assay="RNA", slot="counts")
    m <- as.matrix(m)
    rownames(m) <- sym
    m <- rowsum(m, group=rownames(m))
    obj <- CreateSeuratObject(m, project=obj@project.name, meta.data=obj@meta.data)
  } else {
    rownames(obj) <- sym
  }
  out <- file.path(out_dir, basename(f))
  saveRDS(obj, out, compress="xz")
  cat("Saved:", out, "\n")
}


#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(stringr)
  library(AnnotationDbi); library(org.Hs.eg.db)
})

in_dir  <- "/mnt/18T/chibao/gliomas/data/upstream/snRNA/snRNA_clean_2/rds"
out_dir <- "/mnt/18T/chibao/gliomas/data/upstream/snRNA/snRNA_clean_2/rds_hgnc"; dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# Map Ensembl(±version) -> HGNC symbol
map_to_symbol <- function(ids){
  looks_ens <- grepl("^ENSG\\d+", ids)
  clean_ids <- sub("\\.\\d+$", "", ids)  # drop .version
  syms <- ifelse(looks_ens,
                 mapIds(org.Hs.eg.db, keys=clean_ids, keytype="ENSEMBL", column="SYMBOL", multiVals="first"),
                 ids)
  # If mapping fails, fall back to original (post-version-stripping for Ensembl)
  fallback <- ifelse(looks_ens, clean_ids, ids)
  syms <- ifelse(is.na(syms) | syms == "", fallback, syms)
  syms
}

# Collapse duplicate symbols by summing counts (sparse-aware)
rowsum_sparse <- function(mat, groups){
  # Matrix::rowsum is for dense; use tapply-like grouping for sparse
  grp_levels <- unique(groups)
  idx_list <- split(seq_along(groups), groups)
  res_list <- lapply(idx_list, function(ix){
    if (length(ix) == 1) mat[ix,,drop=FALSE] else Matrix::colSums(mat[ix,,drop=FALSE])
  })
  # Bind efficiently
  out <- do.call(rbind, res_list)
  rownames(out) <- names(idx_list)
  if (!inherits(out, "dgCMatrix")) out <- as(out, "dgCMatrix")
  out
}

rds_files <- list.files(in_dir, pattern="\\.rds$", full.names=TRUE)

for (f in rds_files) {
  obj <- readRDS(f)

  # 1) Always pull the RNA counts (most universal & safe)
  if (!"RNA" %in% names(obj@assays)) {
    stop("Object lacks an RNA assay: ", f)
  }
  counts <- GetAssayData(obj, assay="RNA", slot="counts")
  if (is.null(counts) || nrow(counts) == 0) {
    stop("Empty RNA counts in: ", f)
  }

  # 2) Map IDs to HGNC
  old_ids <- rownames(counts)
  new_syms <- map_to_symbol(old_ids)

  # 3) Collapse duplicates by sum
  if (anyDuplicated(new_syms)) {
    counts <- rowsum_sparse(counts, new_syms)
  } else {
    rownames(counts) <- new_syms
  }

  # 4) Rebuild Seurat object (preserve metadata)
  meta <- obj@meta.data
  # ensure consistent cell order
  meta <- meta[colnames(counts), , drop = FALSE]
  new_obj <- CreateSeuratObject(counts = counts,
                                project = obj@project.name %||% tools::file_path_sans_ext(basename(f)),
                                meta.data = meta,
                                min.cells = 0, min.features = 0)

  # # Optional: re-compute mito/ribo % now that symbols are HGNC
  # # (adjust patterns for your data)
  # mito_genes <- grep("^MT-", rownames(new_obj), value=TRUE)
  # if (length(mito_genes) > 0) {
  #   new_obj[["percent.mt"]] <- PercentageFeatureSet(new_obj, features = mito_genes)
  # }

  # # 5) Re-normalize (choose one):
  # # (A) LogNormalize pipeline
  # new_obj <- NormalizeData(new_obj)
  # new_obj <- FindVariableFeatures(new_obj, nfeatures = 3000)
  # new_obj <- ScaleData(new_obj, features = rownames(new_obj), verbose=FALSE)
  # If you prefer SCT, comment the three lines above and do:
  # new_obj <- SCTransform(new_obj, verbose=FALSE)

  # 6) Save
  out <- file.path(out_dir, basename(f))
  saveRDS(new_obj, out, compress = "xz")
  cat("Saved:", out, "| genes:", nrow(new_obj), "| cells:", ncol(new_obj), "\n")
}


########################## Approach 2 ###########################
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(stringr)
  library(AnnotationDbi); library(org.Hs.eg.db)
})

in_dir  <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/set1/rds"
out_dir <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/set1/rds_hgnc"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Map Ensembl(±version) -> HGNC symbol, with robust fallback
map_to_symbol <- function(ids){
  looks_ens  <- grepl("^ENSG\\d+", ids)
  clean_ids  <- sub("\\.\\d+$", "", ids)  # drop .version
  syms <- ifelse(
    looks_ens,
    mapIds(org.Hs.eg.db,
           keys     = clean_ids,
           keytype  = "ENSEMBL",
           column   = "SYMBOL",
           multiVals= "first"),
    ids
  )
  fallback <- ifelse(looks_ens, clean_ids, ids)
  syms <- ifelse(is.na(syms) | syms == "", fallback, syms)
  syms
}

# Fast sparse row-collapsing by groups (sum)
rowsum_sparse <- function(mat, groups){
  f <- factor(groups)                                 # length = nrow(mat)
  G <- Matrix::sparse.model.matrix(~ 0 + f)           # nrow(mat) x n_groups
  out <- Matrix::t(G) %*% mat                         # n_groups x ncol(mat)
  rownames(out) <- levels(f)
  out
}

rds_files <- list.files(in_dir, pattern = "\\.rds$", full.names = TRUE)
for (f in rds_files) {
  obj <- readRDS(f)

  if (!"RNA" %in% names(obj@assays)) stop("Object lacks an RNA assay: ", f)

  counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
  if (is.null(counts) || nrow(counts) == 0) stop("Empty RNA counts in: ", f)

  old_ids  <- rownames(counts)
  new_syms <- map_to_symbol(old_ids)

  # Collapse duplicate symbols by sum (sparse)
  if (anyDuplicated(new_syms)) {
    counts <- rowsum_sparse(counts, new_syms)
  } else {
    rownames(counts) <- new_syms
  }

  # Rebuild Seurat object and preserve metadata (align by cells)
  meta <- obj@meta.data
  meta <- meta[colnames(counts), , drop = FALSE]

  new_obj <- CreateSeuratObject(
    counts  = counts,
    project = obj@project.name %||% tools::file_path_sans_ext(basename(f)),
    meta.data = meta,
    min.cells = 0, min.features = 0
  )

  # Optional (recommended): recompute percent.mt now that symbols are standardized
  mito <- grep("(?i)^MT[-.]", rownames(new_obj), value = TRUE, perl = TRUE)
  if (length(mito) > 0) new_obj[["percent.mt"]] <- PercentageFeatureSet(new_obj, features = mito)

  # Save
  out <- file.path(out_dir, basename(f))
  saveRDS(new_obj, out, compress = "xz")
  cat(sprintf("Saved: %s | genes: %d | cells: %d\n", out, nrow(new_obj), ncol(new_obj)))
}

