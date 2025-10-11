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


###################
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(stringr)
  library(AnnotationDbi); library(org.Hs.eg.db); library(Matrix)
})

in_dir  <- "/mnt/18T/chibao/gliomas/data/upstream/snRNA/snRNA_clean_2/rds"
out_dir <- "/mnt/18T/chibao/gliomas/data/upstream/snRNA/snRNA_clean_2/rds_hgnc_test"; dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# Map Ensembl (±version) -> HGNC symbol, keep original if no map, ensure uniqueness
map_to_symbol <- function(ids){
  ids        <- as.character(ids)
  looks_ens  <- grepl("^ENSG\\d+", ids)
  clean_ids  <- sub("\\.\\d+$", "", ids)  # drop .version
  symbols    <- ifelse(looks_ens,
                       mapIds(org.Hs.eg.db, keys=clean_ids, keytype="ENSEMBL", column="SYMBOL"),
                       ids)
  # fallback to original if NA/empty
  symbols[is.na(symbols) | symbols==""] <- ids[is.na(symbols) | symbols==""]
  # trim whitespace and normalize weird unicode
  symbols <- str_trim(symbols)
  # Make syntactically valid but preserve case (avoid changing biology)
  # Ensure uniqueness deterministically
  symbols <- make.unique(symbols, sep = "_")
  return(symbols)
}

safe_counts <- function(obj) {
  # Prefer RNA counts explicitly; if missing, try SCT counts; else error
  m <- tryCatch({
    LayerData(obj, assay = "RNA", layer = "counts")
  }, error = function(e) NULL)
  if (is.null(m) || nrow(m) == 0) {
    m <- tryCatch(LayerData(obj, assay = "SCT", layer = "counts"), error = function(e) NULL)
  }
  if (is.null(m) || nrow(m) == 0) stop("No counts layer found in RNA or SCT for object: ", obj@project.name)
  as(m, "dgCMatrix")
}

rds <- list.files(in_dir, pattern="\\.rds$", full.names=TRUE)
for (f in rds) {
  obj <- readRDS(f)

  m_counts <- safe_counts(obj)
  new_genes <- map_to_symbol(rownames(m_counts))

  # Collapse duplicates by **symbol** (sum) using rowsum on dense columns in chunks to save RAM
  rownames(m_counts) <- new_genes
  if (any(duplicated(new_genes))) {
    # rowsum requires numeric matrix; do it sparsely:
    # strategy: split by gene, sum sparse rows
    idx <- split(seq_len(nrow(m_counts)), rownames(m_counts))
    uniq_genes <- names(idx)
    summed <- lapply(idx, function(ix) {
      if (length(ix) == 1) m_counts[ix, , drop=FALSE]
      else Matrix::colSums(m_counts[ix, , drop=FALSE,])
    })
    # Bind back (ensure sparse)
    # Convert vectors to 1-row sparse matrices consistently
    as_dgc <- function(x) {
      if (is.matrix(x)) return(as(x, "dgCMatrix"))
      i <- which(x != 0)
      if (length(i) == 0) return(Matrix(0, nrow=1, ncol=length(x), sparse=TRUE))
      Matrix::sparseMatrix(i = rep(1, length(i)), j = i, x = x[i],
                           dims = c(1, length(x)))
    }
    summed_list <- lapply(summed, as_dgc)
    m_counts <- do.call(rbind, summed_list)
    rownames(m_counts) <- uniq_genes
    # Ensure dimnames preserved
    colnames(m_counts) <- colnames(obj)
  }

  # Rebuild Seurat with aligned meta.data (rownames must equal colnames(counts))
  md <- obj@meta.data
  # Force meta rownames == cell names exactly
  if (!identical(rownames(md), colnames(m_counts))) {
    md <- md[match(colnames(m_counts), rownames(md)), , drop=FALSE]
    rownames(md) <- colnames(m_counts)
  }

  obj2 <- CreateSeuratObject(counts = m_counts,
                             project = obj@project.name,
                             meta.data = md)
  # carry over sample_id if present
  if ("sample_id" %in% colnames(obj@meta.data)) {
    obj2$sample_id <- obj@meta.data[colnames(obj2), "sample_id", drop=TRUE]
  }

  out <- file.path(out_dir, basename(f))
  saveRDS(obj2, out, compress="xz")
  cat("Saved:", out, " | genes:", nrow(obj2), " cells:", ncol(obj2), "\n")
}
