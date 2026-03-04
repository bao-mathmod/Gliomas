# Load libraries
library(Seurat)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(BiocParallel)

# 1. Load your Seurat v5 object
myeloid_obj <- readRDS('/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/subclusters/myeloid/myeloid_clean.rds')
myeloid_obj
# 2. Verify the Default Assay
# SingleR should be run on the actual expression data (RNA), NOT the integrated data.
# Integrated data is for clustering/visualization; raw/normalized data is for annotation.
DefaultAssay(myeloid_obj) <- "RNA"
myeloid_obj
myeloid_obj <- JoinLayers(myeloid_obj)
myeloid_obj
# 3. Handling Seurat v5 Compatibility
# We convert the Seurat object to a SingleCellExperiment (SCE) object.
# This bridges the gap between Seurat v5 layers and SingleR's expectations.
sce <- as.SingleCellExperiment(myeloid_obj, assay = "RNA")

# Alternative manual extraction if conversion fails:
# counts_matrix <- LayerData(myeloid_obj, assay = "RNA", layer = "data")

# FOR HUMAN DATA (Likely your case for Gliomas):
# Blueprint/ENCODE is cited in the paper for annotating human data.
ref <- celldex::BlueprintEncodeData()

# OR Human Primary Cell Atlas (covers a broad range of cell types)
# ref <- celldex::HumanPrimaryCellAtlasData()

# FOR MOUSE DATA (If applicable):
# The paper used ImmGen to annotate mouse lung macrophages[cite: 193, 289].
# ref <- celldex::ImmGenData()

# Run SingleR
# test = your single cell data (sce object)
# ref = the reference dataset loaded above
# labels = the specific annotation column in the reference (usually label.main or label.fine)
bpp <- MulticoreParam(workers = 8)  # use SerialParam() if needed

pred.myeloid <- SingleR(test = sce, 
                        ref = ref, 
                        labels = ref$label.fine,
                        BPPARAM = bpp)

# View the structure of the results
head(pred.myeloid)

# Check the quality of predictions (scores)
# The paper mentions SingleR uses the 80th percentile of correlation values to assign scores [cite: 956]
plotScoreHeatmap(pred.myeloid)

