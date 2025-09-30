# 06_atac_config.sh
# Central configuration for ATAC preprocessing and Cell Ranger ATAC runs

# Root produced by step 05 (symlinked fastqs ready for Cell Ranger ATAC)
LINKED_ATAC_ROOT="/mnt/12T/chibao/data/official_data/atac_prep/fastq_cellranger_atac"

# Path where the manifest will be written (overwritten each run of 06a)
MANIFEST="/mnt/12T/chibao/data/cellranger_data/fastq_cellranger/_cr_atac_samples.tsv"

# Path to Cell Ranger ATAC reference (already downloaded)
REF_PATH="/mnt/12T/chibao/env_tool/cellranger_atac/refgenome_atac/refdata-cellranger-arc-GRCh38-2024-A"

# Where to write Cell Ranger ATAC output (per-project subfolders created)
OUT_ROOT="/mnt/12T/chibao/data/cellranger_data/cell_ranger_atac_out"

# Parallelism and resource controls
MAX_JOBS=2        # number of samples to process at once
LOCAL_CORES="20"  # cores per job
LOCAL_MEM="48"    # GB memory per job
