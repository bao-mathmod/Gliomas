library(Seurat)
library(SeuratDisk)

# 1) Load your Seurat object
seu_path <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/harmony_cleaned_annotated_v3.rds"
seu <- readRDS(seu_path)
backup <- seu
DefaultAssay(seu)
DefaultAssay(seu) <- 'SCT'

# 2) Decide where to write output
out_dir <- "/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult"
setwd(out_dir)

# 3) Save as .h5Seurat
h5s_file <- "harmony_cleaned_annotated_v3.h5Seurat"
SaveH5Seurat(
  seu,
  filename = h5s_file,
  overwrite = TRUE
)

# 4) Convert .h5Seurat → .h5ad
Convert(
  h5s_file,
  dest = "h5ad",
  assay = DefaultAssay(seu),  # often "RNA", "SCT", or "integrated"
  overwrite = TRUE
)
