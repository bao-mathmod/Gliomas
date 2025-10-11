#!/bin/bash
#
# SCRIPT: create_metadata.sh
# PURPOSE: To automatically scan data directories and generate a metadata.csv file
#          for importing various scATAC/snATAC samples into R/Seurat.
#

# --- (USER) PLEASE MODIFY THESE PATHS ---
# -----------------------------------------------------------------------------------
# Base directory containing the 'ATAC_multiome' and 'snATAC' folders
DATA_BASE_DIR="/mnt/18T/chibao/gliomas/data/output_cell"

# Base directory where the lifted hg19->hg38 data was saved
LIFTED_DATA_DIR="/mnt/18T/chibao/gliomas/data/upstream/scATAC/converted_data/PRJNA578617_hg38"

# Desired location for the final metadata.csv file
OUTPUT_FILE="/mnt/18T/chibao/gliomas/data/upstream/scATAC/metadata.csv"
# -----------------------------------------------------------------------------------

# --- Script Initialization ---
set -e
echo "--- Starting metadata generation ---"
# Create the file and write the header
echo "sample_id,project_id,data_type,fragments_path,peaks_path,barcodes_path,matrix_path" > "${OUTPUT_FILE}"

# --- 1. Process Multiome Fragment Files ---
echo "[1/3] Finding multiome fragment samples..."
MULTIOME_DIR="${DATA_BASE_DIR}/ATAC_multiome"
for frag_file in $(find "${MULTIOME_DIR}" -name "*-atac_fragments.tsv.gz" -o -name "*_atac_fragments.tsv.gz"); do
    
    # Extract Project ID (e.g., PRJNA1155637)
    project_id=$(basename $(dirname "${frag_file}"))
    
    # Extract Sample ID (e.g., GSM8492625)
    sample_id=$(basename "${frag_file}" | sed -E 's/(_|\-).+//')

    # Get the full, absolute path to the file
    full_path=$(realpath "${frag_file}")

    # Write to CSV
    echo "${sample_id},${project_id},multiome_fragments,${full_path},NA,NA,NA" >> "${OUTPUT_FILE}"
done

# --- 2. Process Native hg38 snATAC Matrix Files ---
echo "[2/3] Finding native hg38 snATAC matrix samples..."
SNATAC_DIR="${DATA_BASE_DIR}/snATAC"
# Find all peak files, but exclude the original hg19 project and the special RData project
for peak_file in $(find "${SNATAC_DIR}" -name "*_peaks.bed.gz" | grep -v "PRJNA578617/" | grep -v "PRJNA941288/"); do

    # Get the common file prefix (e.g. /path/to/GSM4119513_SF11964_snATAC)
    prefix=$(echo "${peak_file}" | sed 's/_peaks.bed.gz$//')
    
    project_id=$(basename $(dirname "${peak_file}"))
    sample_id=$(basename "${peak_file}" | cut -d'_' -f1)

    # Get absolute paths
    peaks_path=$(realpath "${peak_file}")
    barcodes_path=$(realpath "${prefix}_barcodes.tsv.gz")
    matrix_path=$(realpath "${prefix}_matrix.mtx.gz")

    echo "${sample_id},${project_id},snatac_matrix,NA,${peaks_path},${barcodes_path},${matrix_path}" >> "${OUTPUT_FILE}"
done

# --- 3. Process Lifted hg38 snATAC Matrix Files ---
echo "[3/3] Finding lifted hg38 snATAC matrix samples..."
for peak_file in $(find "${LIFTED_DATA_DIR}" -name "*_peaks.bed.gz"); do

    prefix=$(echo "${peak_file}" | sed 's/_peaks.bed.gz$//')

    # Project ID is known for this batch
    project_id="PRJNA578617" 
    sample_id=$(basename "${peak_file}" | cut -d'_' -f1)

    peaks_path=$(realpath "${peak_file}")
    barcodes_path=$(realpath "${prefix}_barcodes.tsv.gz")
    matrix_path=$(realpath "${prefix}_matrix.mtx.gz")

    echo "${sample_id},${project_id},snatac_matrix,NA,${peaks_path},${barcodes_path},${matrix_path}" >> "${OUTPUT_FILE}"
done

echo "--- Metadata generation complete. File saved to: ${OUTPUT_FILE} ---"
echo "Note: Special case project PRJNA941288 was skipped and must be added manually."