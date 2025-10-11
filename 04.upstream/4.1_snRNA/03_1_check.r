#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(dplyr); library(readr); library(stringr)
})

in_dir  <- "/mnt/18T/chibao/gliomas/data/upstream/snRNA/snRNA_clean_2/rds_hgnc_test"
out_dir <- file.path(in_dir, "_validation"); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

rds_files <- list.files(in_dir, pattern="\\.rds$", full.names=TRUE)
objs <- lapply(rds_files, readRDS)
names(objs) <- basename(rds_files)

get_counts <- function(o) {
  m <- tryCatch(LayerData(o, assay="RNA", layer="counts"), error=function(e) NULL)
  if (is.null(m) || nrow(m)==0) {
    m <- tryCatch(LayerData(o, assay="SCT", layer="counts"), error=function(e) NULL)
  }
  if (is.null(m) || nrow(m)==0) stop("No counts layer found for: ", o@project.name)
  as(m, "dgCMatrix")
}

summaries <- lapply(names(objs), function(nm){
  o <- objs[[nm]]
  m <- get_counts(o)

  # Basic
  genes  <- nrow(m); cells <- ncol(m)
  nnz    <- length(m@x)
  spars  <- 1 - (nnz / (genes * cells))

  # Gene name checks
  rn <- rownames(m)
  has_na_rn   <- any(is.na(rn))
  blank_rn    <- any(rn == "")
  dups_genes  <- any(duplicated(rn))

  # Metadata alignment
  md <- o@meta.data
  md_ok <- identical(rownames(md), colnames(o))

  # QC-ish summaries
  libsize_med <- median(Matrix::colSums(m))
  nfeat_med   <- median(Matrix::colSums(m > 0))

  # MT presence (optional)
  mt_genes <- sum(grepl("^MT-", rn))
  rpl_genes <- sum(grepl("^RPL", rn))
  rps_genes <- sum(grepl("^RPS", rn))

  tibble(
    file = nm,
    genes = genes, cells = cells,
    sparsity = round(spars, 4),
    duplicated_genes = dups_genes,
    na_gene_names = has_na_rn,
    blank_gene_names = blank_rn,
    meta_aligned = md_ok,
    libsize_median = libsize_med,
    nfeatures_median = nfeat_med,
    n_MT_genes = mt_genes,
    n_RPL_genes = rpl_genes,
    n_RPS_genes = rps_genes,
    tiny_flag = cells < 200
  )
})

summary_df <- bind_rows(summaries)
write_tsv(summary_df, file.path(out_dir, "per_object_summary.tsv"))
cat("Wrote:", file.path(out_dir, "per_object_summary.tsv"), "\n")

# Pairwise gene-set Jaccard (feature overlap)
gene_sets <- lapply(objs, function(o) rownames(get_counts(o)))
obj_names <- names(objs)
J <- matrix(NA_real_, nrow=length(objs), ncol=length(objs),
            dimnames=list(obj_names, obj_names))
for (i in seq_along(objs)) for (j in seq_along(objs)) {
  a <- gene_sets[[i]]; b <- gene_sets[[j]]
  J[i,j] <- length(intersect(a,b)) / length(union(a,b))
}
J_df <- as.data.frame(J) %>% tibble::rownames_to_column("object")
write_tsv(J_df, file.path(out_dir, "pairwise_jaccard_features.tsv"))
cat("Wrote:", file.path(out_dir, "pairwise_jaccard_features.tsv"), "\n")

# Hard stop hints
if (any(!summary_df$meta_aligned)) {
  bad <- summary_df$file[!summary_df$meta_aligned]
  stop("meta.data not aligned to cell names for: ", paste(bad, collapse=", "))
}
if (any(summary_df$duplicated_genes)) {
  bad <- summary_df$file[summary_df$duplicated_genes]
  stop("Duplicated gene names remain after Script 1 in: ", paste(bad, collapse=", "),
       "\nInvestigate symbol mapping/collapse for these files.")
}
