base_dir <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/subclusters/myeloid/annotated"
in_rds   <- file.path(base_dir, "myeloid_clean_annotated.rds")

out_dir  <- file.path('/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/subclusters/myeloid/GSVA')
# dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE)
dir.create(file.path(out_dir, "rds"), showWarnings = FALSE)

set.seed(1)
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
})

obj <- readRDS(in_rds)

# Ensure you only keep annotated cells
obj <- subset(obj, subset = !is.na(general_cell_type))

# Set identity to your myeloid state label
Idents(obj) <- "general_cell_type"

# Check what you have available
DefaultAssay(obj)
Assays(obj)
head(obj@meta.data[, c("general_cell_type")], 3)

# 2. Build Pseudobulk
suppressPackageStartupMessages({
  library(Seurat)
})

# Pick your sample column (adjust if you use sample_uid/patient_id instead)
sample_col <- "sample_uid"

stopifnot(sample_col %in% colnames(obj@meta.data))

# Pseudobulk counts for each sample × myeloid state
pb_counts <- AggregateExpression(
  obj,
  assays = "RNA",
  slot = "counts",
  group.by = c(sample_col, "general_cell_type"),
  return.seurat = FALSE
)$RNA

dim(pb_counts)  # genes x (sample_state)
pb_counts[1:3, 1:3]

# 3. Convert Pseudobulk
suppressPackageStartupMessages({
  library(edgeR)
})

y <- DGEList(counts = pb_counts)
y <- calcNormFactors(y)
pb_logcpm <- cpm(y, log = TRUE, prior.count = 1)

dim(pb_logcpm)

# 4. GSVA loading
suppressPackageStartupMessages({
  library(msigdbr)
  library(tibble)
})

# Hallmark (50 sets)
h_df <- msigdbr(species = "Homo sapiens", category = "H")
hallmark <- split(h_df$gene_symbol, h_df$gs_name)

length(hallmark)  # should be ~50

c2_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")

# Keyword filter focused on myeloid biology
keep_kw <- c(
  "INTERFERON", "ANTIGEN", "MHC", "TOLL", "NF_KB", "CYTOKINE", "CHEMOKINE",
  "COMPLEMENT", "PHAGOCYT", "NEUTROPHIL", "MONOCYTE", "MACROPHAGE",
  "HYPOX", "GLYCOLYSIS", "OXIDATIVE", "RESPIRATORY",
  "APOPTOS", "ANGIOGEN", "VEGF"
)

c2_df2 <- c2_df %>%
  filter(str_detect(gs_name, str_c(keep_kw, collapse="|")))

reactome_myeloid <- split(c2_df2$gene_symbol, c2_df2$gs_name)
length(reactome_myeloid)

# Merge gene sets
all_sets <- c(hallmark, reactome_myeloid)

# Clean gene symbols (strip weird whitespace)
all_sets <- lapply(all_sets, \(g) unique(str_trim(g)))

# Keep only genes present in pseudobulk matrix
genes_present <- rownames(pb_logcpm)
all_sets <- lapply(all_sets, \(g) intersect(g, genes_present))

# Size filter (good defaults for scRNA pseudobulk)
minSize <- 10
maxSize <- 300
all_sets <- all_sets[lengths(all_sets) >= minSize & lengths(all_sets) <= maxSize]

length(all_sets)
summary(lengths(all_sets))

# Run 
suppressPackageStartupMessages({
  library(GSVA)
})

# GSVA method (pathway-centric variation across samples)
gp <- gsvaParam(
  exprData  = pb_logcpm,
  geneSets  = all_sets,
  kcdf      = "Gaussian",
  minSize   = minSize,
  maxSize   = maxSize,
  maxDiff   = TRUE
)

gsva_mat <- gsva(gp)  # pathways x pseudo-samples
dim(gsva_mat)

# Save
saveRDS(gsva_mat, file.path(out_dir, "rds/gsva_pseudobulk_matrix.rds"))
write.csv(gsva_mat, file.path(out_dir, "tables/gsva_pseudobulk_matrix.csv"))

# Build clean sample annotation
pb_cols <- colnames(gsva_mat)

annot <- tibble(pb_id = pb_cols) %>%
  tidyr::separate(pb_id, into = c("sample", "general_cell_type"), sep = "_", extra="merge", remove = FALSE)

# If your sample names themselves contain underscores, you may need a safer delimiter;
# another approach is to create a combined key yourself before AggregateExpression.

annot$general_cell_type <- factor(annot$general_cell_type)
head(annot)
saveRDS(annot, file.path(out_dir, "rds/pseudobulk_annotation.rds"))

# Viz
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
})

# Z-score per pathway for heatmap readability
gsva_z <- t(scale(t(gsva_mat)))

# Show only Hallmark as an initial clean figure
hallmark_names <- names(hallmark)
keep_rows <- intersect(rownames(gsva_z), hallmark_names)

mat_h <- gsva_z[keep_rows, , drop=FALSE]

ha <- HeatmapAnnotation(
  cell_state = annot$general_cell_type,
  show_annotation_name = TRUE
)

pdf(file.path(out_dir, "plots/heatmap_hallmark_gsva.pdf"), width = 10, height = 8)
Heatmap(
  mat_h,
  top_annotation = ha,
  name = "GSVA (z)",
  show_column_names = FALSE,
  show_row_names = TRUE,
  cluster_columns = TRUE,
  cluster_rows = TRUE
)
dev.off()

# Option 2: 
library(Seurat)
library(edgeR)
library(GSVA)
library(msigdbr)
library(dplyr)
library(stringr)
library(ComplexHeatmap)

obj <- readRDS(in_rds)
obj <- subset(obj, subset = !is.na(general_cell_type))
Idents(obj) <- "general_cell_type"

# 1) Pseudobulk counts per cell type (columns = cell types only)
pb_counts_ct <- AggregateExpression(
  obj,
  assays = "RNA",
  slot = "counts",
  group.by = "general_cell_type",
  return.seurat = FALSE
)$RNA

# 2) logCPM
y <- edgeR::DGEList(counts = pb_counts_ct)
y <- edgeR::calcNormFactors(y)
pb_logcpm_ct <- edgeR::cpm(y, log = TRUE, prior.count = 1)

# 3) Gene sets: C2 / CP:REACTOME, keyword-filtered for myeloid biology
c2_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")

keep_kw <- c(
  "INTERFERON", "ANTIGEN", "MHC", "TOLL",
  "NF_KB", "NFKB", "CYTOKINE", "CHEMOKINE",
  "COMPLEMENT", "PHAGOCYT", "NEUTROPHIL", "MONOCYTE", "MACROPHAGE",
  "HYPOX", "GLYCOLYSIS", "OXIDATIVE", "RESPIRATORY",
  "APOPTOS", "ANGIOGEN", "VEGF"
)

c2_df2 <- c2_df %>%
  filter(str_detect(gs_name, regex(paste(keep_kw, collapse = "|"), ignore_case = TRUE)))

reactome_sets <- split(c2_df2$gene_symbol, c2_df2$gs_name)

# Filter sets to genes present + size
genes_present <- rownames(pb_logcpm_ct)
reactome_sets <- lapply(reactome_sets, function(g) intersect(unique(g), genes_present))
reactome_sets <- reactome_sets[lengths(reactome_sets) >= 10 & lengths(reactome_sets) <= 300]

if (length(reactome_sets) == 0) {
  stop("No Reactome gene sets left after keyword + size filtering. Check subcategory='CP:REACTOME' and your keep_kw.")
}
message("Reactome sets kept: ", length(reactome_sets))

# 4) GSVA
gp <- gsvaParam(exprData = pb_logcpm_ct, geneSets = reactome_sets, kcdf = "Gaussian")
gsva_ct <- gsva(gp)  # pathways x cell_types

# Optional: keep only top variable pathways for a readable heatmap
top_k <- 80  # set NULL to disable
if (!is.null(top_k) && nrow(gsva_ct) > top_k) {
  sds <- apply(gsva_ct, 1, sd)
  keep <- names(sort(sds, decreasing = TRUE))[1:top_k]
  gsva_ct <- gsva_ct[keep, , drop = FALSE]
  message("Heatmap pathways shown (top_k): ", nrow(gsva_ct))
}

# 5) Heatmap (z-score per pathway to highlight "high in this cell type")
gsva_ct_z <- t(scale(t(gsva_ct)))

pdf(file.path(out_dir, "plots/heatmap_REACTOME_gsva_by_celltype.pdf"), width = 7, height = 9)
Heatmap(gsva_ct_z, name = "GSVA (z)", cluster_rows = TRUE, cluster_columns = TRUE)
dev.off()

# Option 3
library(stringr)
library(tidyr)
library(dplyr)
library(edgeR)
library(GSVA)
library(msigdbr)
library(ComplexHeatmap)

obj <- readRDS(in_rds)
obj <- subset(obj, subset = !is.na(general_cell_type))
stopifnot("sample_uid" %in% colnames(obj@meta.data))

# Make a safe combined group id
obj$pb_group <- paste(obj$sample_uid, obj$general_cell_type, sep = "||")

# 1) Pseudobulk per (sample_uid × cell_type)
pb_counts <- AggregateExpression(
  obj,
  assays = "RNA",
  slot = "counts",
  group.by = "pb_group",
  return.seurat = FALSE
)$RNA

# 2) logCPM
y <- edgeR::DGEList(counts = pb_counts)
y <- edgeR::calcNormFactors(y)
pb_logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 1)

# 3) Gene sets: C2 / CP:REACTOME, keyword-filtered for myeloid biology
c2_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")

keep_kw <- c(
  "INTERFERON", "ANTIGEN", "MHC", "TOLL",
  "NF_KB", "NFKB", "CYTOKINE", "CHEMOKINE",
  "COMPLEMENT", "PHAGOCYT", "NEUTROPHIL", "MONOCYTE", "MACROPHAGE",
  "HYPOX", "GLYCOLYSIS", "OXIDATIVE", "RESPIRATORY",
  "APOPTOS", "ANGIOGEN", "VEGF"
)

c2_df2 <- c2_df %>%
  filter(str_detect(gs_name, regex(paste(keep_kw, collapse = "|"), ignore_case = TRUE)))

reactome_sets <- split(c2_df2$gene_symbol, c2_df2$gs_name)

# Filter sets to genes present + size
genes_present <- rownames(pb_logcpm)
reactome_sets <- lapply(reactome_sets, function(g) intersect(unique(g), genes_present))
reactome_sets <- reactome_sets[lengths(reactome_sets) >= 10 & lengths(reactome_sets) <= 300]

if (length(reactome_sets) == 0) {
  stop("No Reactome gene sets left after keyword + size filtering. Check subcategory='CP:REACTOME' and your keep_kw.")
}
message("Reactome sets kept: ", length(reactome_sets))

# 4) GSVA
gp <- gsvaParam(exprData = pb_logcpm, geneSets = reactome_sets, kcdf = "Gaussian")
gsva_mat <- gsva(gp)   # pathways x (sample_uid||cell_type)

# 5) Annotation: split pb_group back into sample + cell_type
annot <- tibble(pb_group = colnames(gsva_mat)) %>%
  tidyr::separate(pb_group, into = c("sample_uid", "general_cell_type"), sep = "\\|\\|", remove = FALSE)

annot$general_cell_type <- factor(annot$general_cell_type)

# 6) Collapse to cell-type mean GSVA (columns = cell types only)
cell_types <- levels(annot$general_cell_type)

gsva_mean_by_ct <- sapply(cell_types, function(ct) {
  cols <- annot$pb_group[annot$general_cell_type == ct]
  rowMeans(gsva_mat[, cols, drop = FALSE])
})

# (optional) "% positive across samples" per pathway per cell type
gsva_posrate_by_ct <- sapply(cell_types, function(ct) {
  cols <- annot$pb_group[annot$general_cell_type == ct]
  rowMeans(gsva_mat[, cols, drop = FALSE] > 0)
})

saveRDS(gsva_mean_by_ct, file.path(out_dir, "rds/gsva_REACTOME_mean_by_celltype.rds"))
write.csv(gsva_mean_by_ct, file.path(out_dir, "tables/gsva_REACTOME_mean_by_celltype.csv"))

# Optional: keep only top variable pathways for a readable heatmap
top_k <- 80  # set NULL to disable
mat_show <- gsva_mean_by_ct
if (!is.null(top_k) && nrow(mat_show) > top_k) {
  sds <- apply(mat_show, 1, sd)
  keep <- names(sort(sds, decreasing = TRUE))[1:top_k]
  mat_show <- mat_show[keep, , drop = FALSE]
  message("Heatmap pathways shown (top_k): ", nrow(mat_show))
}

# 7) Heatmap for figure (z-score per pathway across cell types)
gsva_mean_z <- t(scale(t(mat_show)))

pdf(file.path(out_dir, "plots/heatmap_REACTOME_gsva_celltype_MEAN.pdf"), width = 7, height = 9)
Heatmap(gsva_mean_z, name = "GSVA mean (z)", cluster_rows = TRUE, cluster_columns = TRUE)
dev.off()
