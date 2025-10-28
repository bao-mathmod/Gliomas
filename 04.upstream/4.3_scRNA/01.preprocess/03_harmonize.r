# Different authors mix Ensembl IDs (±version) and gene symbols. Standardize row names to HGNC symbols (drop Ensembl versions), then deduplicate.
# Approach 1
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(stringr)
  library(AnnotationDbi); library(org.Hs.eg.db)
})

in_dir  <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/set1/rds"
out_dir <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/set1/rds_hgnc"; dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

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

########################## Approach 2 ###########################
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(stringr)
  library(AnnotationDbi); library(org.Hs.eg.db)
})

in_dir  <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/set1_manifest/rds"
out_dir <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/set1_manifest/rds_hgnc"
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

rds_files <- list.files(in_dir, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
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

########################## Biomart ############################
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(stringr); library(data.table)
  library(biomaRt)
  suppressWarnings({suppressMessages(library(AnnotationDbi)); suppressMessages(library(org.Hs.eg.db))})
})

# ===================== CONFIG =====================
in_dir   <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/set1/rds"        # input RDS (unmapped / mixed)
out_dir  <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/set1/rds_hgnc_biomart"   # output mapped RDS
cache_dir <- file.path(out_dir, "_biomart_cache")                         # ENSG->HGNC cache
chunk_size <- 5000     # biomaRt query chunk
retries    <- 3        # biomaRt retries per chunk
sleep_sec  <- 3        # wait between retries
recalc_percent_mt <- TRUE
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
# ==================================================

`%||%` <- function(a,b) if (is.null(a) || length(a)==0) b else a
strip_ver <- function(x) sub("\\.\\d+$","", x)

rowsum_sparse <- function(mat, groups){
  f <- factor(groups)                                  # length = nrow(mat)
  G <- Matrix::sparse.model.matrix(~ 0 + f)            # nrow x n_groups
  out <- Matrix::t(G) %*% mat                           # n_groups x ncol
  rownames(out) <- levels(f)
  out
}

bm_map_cached <- function(ens_ids, mart, chunk=5000, retries=3, sleep=3, cache_dir){
  ens_ids <- unique(strip_ver(ens_ids))
  if (!length(ens_ids)) return(setNames(character(0), character(0)))
  key <- file.path(cache_dir, "cache.tsv.gz")
  cache <- if (file.exists(key)) fread(key, sep="\t", header=TRUE) else
            data.frame(ensembl_gene_id=character(), hgnc_symbol=character())
  out <- setNames(rep(NA_character_, length(ens_ids)), ens_ids)

  # fill from cache
  m <- match(ens_ids, cache$ensembl_gene_id)
  hit <- which(!is.na(m))
  if (length(hit)) out[hit] <- cache$hgnc_symbol[m[hit]]

  need <- names(out)[is.na(out)]
  if (length(need)) {
    pieces <- split(need, ceiling(seq_along(need)/chunk))
    fetched <- lapply(pieces, function(vec){
      attempt <- 1
      repeat {
        res <- try({
          getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
                filters    = "ensembl_gene_id",
                values     = vec,
                mart       = mart)
        }, silent = TRUE)
        if (!inherits(res, "try-error")) return(res)
        if (attempt >= retries) {
          warning(sprintf("biomaRt failed after %d retries; returning NAs for %d IDs.", retries, length(vec)))
          return(data.frame(ensembl_gene_id = vec, hgnc_symbol = NA_character_))
        }
        Sys.sleep(sleep); attempt <- attempt + 1
      }
    })
    fetched <- rbindlist(fetched, use.names=TRUE, fill=TRUE)
    # update outs
    m2 <- match(names(out), fetched$ensembl_gene_id)
    upd <- which(!is.na(m2))
    out[upd] <- fetched$hgnc_symbol[m2[upd]]
    # update cache (dedup by last)
    cache <- rbind(cache, fetched)
    cache <- cache[!duplicated(cache$ensembl_gene_id, fromLast = TRUE), ]
    fwrite(cache, key, sep="\t")
  }
  out
}

# connect once
message("Connecting to Ensembl via biomaRt ...")
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "useast")

files <- list.files(in_dir, pattern="\\.rds$", full.names=TRUE, recursive = TRUE)
for (i in seq_along(files)) {
  f <- files[i]
  cat(sprintf("\n[%d/%d] %s\n", i, length(files), basename(f)))

  obj <- readRDS(f)
  stopifnot("RNA" %in% names(obj@assays))

  counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
  stopifnot(inherits(counts, "dgCMatrix"))

  # --- before-stats (for logging)
  pre_genes <- nrow(counts)
  pre_cells <- ncol(counts)
  pre_ensg  <- sum(grepl("^ENSG\\d+(\\.\\d+)?$", rownames(counts)))

  # strip Ensembl version
  base_ids <- strip_ver(rownames(counts))  # length == nrow(counts)

  # 1) biomaRt mapping (returns unique ENSG -> SYMBOL)
  bm_map <- bm_map_cached(base_ids, mart,
                          chunk = chunk_size, retries = retries,
                          sleep = sleep_sec, cache_dir = cache_dir)
  # 2) RE-INDEX to per-row order (critical!)
  new_syms <- bm_map[match(base_ids, names(bm_map))]

  # 3) fallback to org.Hs.eg.db for remaining NA/blank
  need_idx <- which(is.na(new_syms) | new_syms == "")
  if (length(need_idx) && requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    need_keys <- unique(base_ids[need_idx])
    db_map <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                    keys      = need_keys,
                                    keytype   = "ENSEMBL",
                                    column    = "SYMBOL",
                                    multiVals = "first")
    m <- match(base_ids[need_idx], names(db_map))
    hit <- which(!is.na(m))
    if (length(hit)) new_syms[need_idx[hit]] <- db_map[m[hit]]
  }

  # 4) final fallback: keep base ENSG for any still NA/blank
  still_na <- which(is.na(new_syms) | new_syms == "")
  if (length(still_na)) new_syms[still_na] <- base_ids[still_na]

  # --- sanity checks before collapsing
  if (length(new_syms) != nrow(counts)) {
    stop(sprintf("Grouping labels length (%d) != nrow(counts) (%d)",
                 length(new_syms), nrow(counts)))
  }
  if (anyNA(new_syms)) {
    bad <- sum(is.na(new_syms))
    stop(sprintf("Found %d NA group labels after mapping; check fallback logic.", bad))
  }

  # 5) apply names & collapse duplicates sparsely
  if (anyDuplicated(new_syms)) {
    counts <- rowsum_sparse(counts, new_syms)  # requires rowsum_sparse defined earlier
  } else {
    rownames(counts) <- new_syms
  }

  # 6) rebuild object, preserve metadata (align by cells)
  meta <- obj@meta.data
  meta <- meta[colnames(counts), , drop = FALSE]

  project_name <- if (!is.null(obj@project.name) && length(obj@project.name) > 0) {
    obj@project.name
  } else {
    tools::file_path_sans_ext(basename(f))
  }

  new_obj <- CreateSeuratObject(
    counts    = counts,
    project   = project_name,
    meta.data = meta,
    min.cells = 0,
    min.features = 0
  )

  if (recalc_percent_mt) {
    mito <- grep("(?i)^MT[-.]", rownames(new_obj), value = TRUE, perl = TRUE)
    if (length(mito)) new_obj[["percent.mt"]] <- PercentageFeatureSet(new_obj, features = mito)
  }

  new_obj@misc$mapping <- list(
    method = "biomaRt-first; fallback org.Hs.eg.db; fallback ENSG",
    cache  = cache_dir,
    date   = as.character(Sys.time())
  )

  out <- file.path(out_dir, basename(f))
  saveRDS(new_obj, out, compress = "xz")

  # --- after-stats (log)
  post_counts <- GetAssayData(new_obj, assay = "RNA", slot = "counts")
  post_genes  <- nrow(post_counts)
  post_ensg   <- sum(grepl("^ENSG\\d+(\\.\\d+)?$", rownames(post_counts)))
  cat(sprintf("Saved: %s | genes: %d -> %d | cells: %d | remaining_ENSG: %d -> %d\n",
              out, pre_genes, post_genes, pre_cells, pre_ensg, post_ensg))
}
