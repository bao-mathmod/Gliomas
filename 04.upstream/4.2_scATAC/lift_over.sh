#!/bin/bash
#
# SCRIPT: run_liftover_PRJNA578617.sh
# PURPOSE: Convert peak coordinates for all samples in PRJNA578617 from hg19 to hg38.
# VERSION: Conda-compatible
#

# --- (USER) PLEASE MODIFY THESE PATHS ---
# -----------------------------------------------------------------------------------
INPUT_BASE_DIR="/mnt/18T/chibao/gliomas/data/output_cell/snATAC"
OUTPUT_BASE_DIR="/mnt/18T/chibao/gliomas/data/upstream/scATAC"
# -----------------------------------------------------------------------------------

# --- Automatic Path Configuration ---
set -e
PROJECT_DIR="${INPUT_BASE_DIR}/PRJNA578617"
CONVERTED_DATA_DIR="${OUTPUT_BASE_DIR}/converted_data/PRJNA578617_hg38"
LOG_DIR="${OUTPUT_BASE_DIR}/logs"
TOOLS_DIR="${OUTPUT_BASE_DIR}/tools"

mkdir -p "${CONVERTED_DATA_DIR}"
mkdir -p "${LOG_DIR}"
mkdir -p "${TOOLS_DIR}"

LOG_FILE="${LOG_DIR}/liftover_summary_$(date +%F).log"
CHAIN_FILE="${TOOLS_DIR}/hg19ToHg38.over.chain.gz"

# --- 1. Download Chain File (liftOver executable is now managed by Conda) ---
echo "--- Step 1: Checking for and downloading chain file ---"
if [ ! -f "$CHAIN_FILE" ]; then
    echo "Downloading hg19 to hg38 chain file to ${CHAIN_FILE}..."
    wget -qO "$CHAIN_FILE" http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
else
    echo "Chain file already exists."
fi
echo "--- Tools are ready. ---"
echo ""

# --- 2. Initialize Log File ---
echo "LiftOver Summary for Project ${PROJECT_DIR}" > ${LOG_FILE}
# ... (rest of the script is identical) ...
echo "Run started on: $(date)" >> ${LOG_FILE}
echo "Output will be saved in: ${CONVERTED_DATA_DIR}" >> ${LOG_FILE}
echo "--------------------------------------------------" >> ${LOG_FILE}

# --- 3. Process Each Sample ---
echo "--- Step 2: Processing samples in ${PROJECT_DIR} ---"
SAMPLES=$(find "${PROJECT_DIR}" -name "*_peaks.bed.gz" | sed -E "s|(.*)_peaks.bed.gz|\1|")

for sample_prefix in ${SAMPLES}; do
    sample_base=$(basename "${sample_prefix}")
    echo "Processing sample: ${sample_base}"

    original_peaks_gz="${sample_prefix}_peaks.bed.gz"
    original_matrix_gz="${sample_prefix}_matrix.mtx.gz"
    original_barcodes_gz="${sample_prefix}_barcodes.tsv.gz"
    
    output_dir_sample="${CONVERTED_DATA_DIR}/${sample_base}"
    mkdir -p "${output_dir_sample}"

    temp_peaks_hg19="${output_dir_sample}/${sample_base}_hg19_peaks.bed"
    temp_peaks_hg38="${output_dir_sample}/${sample_base}_hg38_peaks.bed"
    unmapped_log="${output_dir_sample}/${sample_base}_unmapped.log"
    final_peaks_gz="${output_dir_sample}/${sample_base}_peaks.bed.gz"

    gunzip -c "${original_peaks_gz}" > "${temp_peaks_hg19}"

    # MODIFICATION #2: Call liftOver directly from the Conda environment
    liftOver "${temp_peaks_hg19}" "${CHAIN_FILE}" "${temp_peaks_hg38}" "${unmapped_log}"

    gzip -c "${temp_peaks_hg38}" > "${final_peaks_gz}"
    cp "${original_matrix_gz}" "${output_dir_sample}/"
    cp "${original_barcodes_gz}" "${output_dir_sample}/"

    # --- 4. Logging and Cleanup ---
    original_count=$(wc -l < "${temp_peaks_hg19}")
    lifted_count=$(wc -l < "${temp_peaks_hg38}")
    unmapped_count=$(wc -l < "${unmapped_log}")

    echo "Sample: ${sample_base}" >> ${LOG_FILE}
    echo "  - Original hg19 peaks: ${original_count}" >> ${LOG_FILE}
    echo "  - Successfully lifted to hg38: ${lifted_count}" >> ${LOG_FILE}
    echo "  - Unmapped regions: ${unmapped_count}" >> ${LOG_FILE}
    echo "  - Output directory: ${output_dir_sample}" >> ${LOG_FILE}
    echo "--------------------------------------------------" >> ${LOG_FILE}
    
    rm "${temp_peaks_hg19}" "${temp_peaks_hg38}"
    
    echo "Finished processing ${sample_base}."
    echo ""
done

echo "--- LiftOver process complete. See ${LOG_FILE} for summary. ---"