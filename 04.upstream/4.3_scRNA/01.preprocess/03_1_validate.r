library(Seurat); library(Matrix); library(data.table)

in_dir  <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/set1_1/rds_hgnc"
files   <- list.files(in_dir, pattern="\\.rds$", full.names=TRUE)

audit_one <- function(path){
  obj <- readRDS(path)
  counts <- GetAssayData(obj, assay="RNA", slot="counts")

  n_genes <- nrow(counts); n_cells <- ncol(counts)
  n_ensg  <- sum(grepl("^ENSG\\d+(\\.\\d+)?$", rownames(counts)))
  pct_ensg <- if (n_genes>0) 100*n_ensg/n_genes else NA_real_

  dup_n   <- sum(duplicated(rownames(counts)))
  n_mt    <- sum(grepl("(?i)^MT[-.]", rownames(counts), perl=TRUE))
  dens    <- round(100 * length(counts@x) / (n_genes * n_cells), 4)

  # capture a few example remaining ENSG if any
  ex_ensg <- paste(head(rownames(counts)[grepl("^ENSG", rownames(counts))], 5), collapse=";")

  data.frame(
    file = basename(path),
    genes = n_genes, cells = n_cells,
    remaining_ENSG = n_ensg, pct_remaining_ENSG = round(pct_ensg,3),
    duplicates = dup_n, mito_genes = n_mt, sparsity_pct = 100 - dens,
    examples_remaining_ENSG = ex_ensg,
    stringsAsFactors = FALSE
  )
}

aud <- rbindlist(lapply(files, audit_one), fill=TRUE)
out <- file.path(in_dir, "mapping_audit.tsv")
fwrite(aud, out, sep="\t")
print(aud[order(-pct_remaining_ENSG), ][1:min(10, nrow(aud)), ])
cat("Wrote: ", out, "\n")


pre  <- readRDS("/mnt/18T/chibao/gliomas/data/upstream/scRNA/set1/rds/PRJNA683876/PRJNA683876__PRJNA683876__filtered_feature_bc_matrix__SAMN17039397.rds")
post <- readRDS("/mnt/18T/chibao/gliomas/data/upstream/scRNA/set1/rds_hgnc/PRJNA683876__PRJNA683876__filtered_feature_bc_matrix__SAMN17039397.rds")
pre_counts  <- GetAssayData(pre,  assay="RNA", slot="counts")
post_counts <- GetAssayData(post, assay="RNA", slot="counts")

# number of ENSG rows removed via merging
pre_ensg  <- sum(grepl("^ENSG", rownames(pre_counts)))
post_ensg <- sum(grepl("^ENSG", rownames(post_counts)))

cat(sprintf("ENSG rows: pre=%d → post=%d\n", pre_ensg, post_ensg))
cat(sprintf("Genes: pre=%d → post=%d (Δ=%d)\n",
            nrow(pre_counts), nrow(post_counts), nrow(post_counts)-nrow(pre_counts)))

# Optional: prove counts are conserved after merge (within a small tolerance)
total_pre  <- sum(pre_counts@x)
total_post <- sum(post_counts@x)
cat(sprintf("Total UMI counts: pre=%.0f vs post=%.0f (should be equal)\n", total_pre, total_post))


##################################3
#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

# ===================== CONFIG =====================
in_dir  <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/set1_manifest/rds"
out <- "/mnt/18T/chibao/gliomas/data/stuff"
out_tsv <- file.path(out, "raw_ENSG_audit_set1_mani.tsv")
# ==================================================

audit_one <- function(path){
  obj <- readRDS(path)
  stopifnot("RNA" %in% names(obj@assays))
  counts <- GetAssayData(obj, assay="RNA", slot="counts")

  n_genes <- nrow(counts)
  n_cells <- ncol(counts)
  n_ensg  <- sum(grepl("^ENSG\\d+(\\.\\d+)?$", rownames(counts)))
  pct_ensg <- round(100 * n_ensg / n_genes, 3)
  ex_ensg <- paste(head(rownames(counts)[grepl("^ENSG", rownames(counts))], 5), collapse=";")

  data.frame(
    file = basename(path),
    genes = n_genes,
    cells = n_cells,
    remaining_ENSG = n_ensg,
    pct_remaining_ENSG = pct_ensg,
    examples_ENSG = ex_ensg,
    stringsAsFactors = FALSE
  )
}

files <- list.files(in_dir, pattern="\\.rds$", full.names=TRUE, recursive=TRUE)
res <- rbindlist(lapply(files, audit_one), use.names=TRUE, fill=TRUE)
fwrite(res, out_tsv, sep="\t")

cat("Wrote:", out_tsv, "\n\n")
print(res[order(-pct_remaining_ENSG), ])
