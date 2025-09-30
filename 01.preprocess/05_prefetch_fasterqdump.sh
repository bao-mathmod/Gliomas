# prefetch SRR21832248
# fasterq-dump --split-3 --threads 8 SRR21832248

# fasterq-dump --split-files --include-technical --threads 8 SRR21832249

#!/bin/bash

# # Define an array of SRA accession numbers
# SRA_ACCESSIONS=("SRR10315837" "SRR10315835" "SRR10315838")

# # Define log and summary file names
# LOG_FILE="fastq_dump_progress.log"
# SUMMARY_FILE="fastq_dump_summary.txt"

# # Clear previous logs and summary files
# > "$LOG_FILE"
# > "$SUMMARY_FILE"

# # Log script start time
# echo "Script started at $(date)" | tee -a "$LOG_FILE"

# # Loop through each accession number
# for accession in "${SRA_ACCESSIONS[@]}"; do
#     start_time=$(date +%s)
#     echo "Processing $accession..." | tee -a "$LOG_FILE"
    
#     # Run prefetch and log its output
#     echo "Running prefetch for $accession..." | tee -a "$LOG_FILE"
#     prefetch "$accession" >> "$LOG_FILE" 2>&1
    
#     # Run fasterq-dump and log its output
#     echo "Running fasterq-dump for $accession..." | tee -a "$LOG_FILE"
#     fasterq-dump --split-files --include-technical --threads 8 "$accession" >> "$LOG_FILE" 2>&1
    
#     end_time=$(date +%s)
#     duration=$((end_time - start_time))
    
#     # Create a summary entry
#     echo "----------------------------------------" | tee -a "$LOG_FILE"
#     echo "Summary for $accession:" | tee -a "$SUMMARY_FILE"
#     echo "Status: COMPLETED" | tee -a "$SUMMARY_FILE"
#     echo "Duration: $duration seconds" | tee -a "$SUMMARY_FILE"
#     echo "Log file: $LOG_FILE" | tee -a "$SUMMARY_FILE"
#     echo "----------------------------------------" | tee -a "$LOG_FILE"
# done

# # Log script completion time
# echo "Script finished at $(date)" | tee -a "$LOG_FILE"
# echo "All processing complete. Check $SUMMARY_FILE for a summary and $LOG_FILE for detailed output."


# proj=PRJNA578617
# run=SRR10315838
# meta="/mnt/12T/chibao/data/official_data/metadata_enriched/$proj/${proj}.extra.tsv"

# awk -F'\t' -v r="$run" '
#   NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
#   $h["run_accession"]==r {
#     printf "run=%s\tlibrary_strategy=%s\tlibrary_layout=%s\tinstrument_platform=%s\n",
#            r, $h["library_strategy"], $h["library_layout"], $h["instrument_platform"]
#   }
# ' "$meta"

# fdir="/mnt/12T/chibao/data/official_data/fastq/PRJNA578617/SAMN13066825/SRX7026892/SRR10315838"

# # Length of the first sequence of each file
# for f in "$fdir"/SRR10315838*.fastq.gz; do
#   echo ">>> $(basename "$f")"
#   zcat "$f" | awk 'NR%4==2{print "read_length=" length; exit}'
# done

# # Look at the first read header of the unlabeled file
# zcat "$fdir/SRR10315838.fastq.gz" | head -n 4


# BASE_DIR_SCAN="/mnt/12T/chibao/data/official_data/fastq"
# run_dir="$BASE_DIR_SCAN/PRJNA683876/SAMN17039397/SRX9660939/SRR13228579"

# rel="${run_dir#$BASE_DIR_SCAN/}"
# IFS='/' read -r study sample exp run <<< "$rel"
# printf "study=%s\nsample=%s\nexp=%s\nrun=%s\n" "$study" "$sample" "$exp" "$run"
# Expect:
# study=PRJNA683876
# sample=SAMN17039397
# exp=SRX9660939
# run=SRR13228579

#!/bin/bash

# Define the base directory to search for SRA accessions
BASE_PATH="/mnt/12T/chibao/data/official_data/fastq_by_tech/scATAC/PRJNA961045"

# --- Define the output directory for FASTQ and log files ---
OUTPUT_PATH="/mnt/12T/chibao/data/official_data/fastq_by_tech/bulk_or_plate_RNA/official_result/PRJNA961045"
mkdir -p "$OUTPUT_PATH"

# --- Control Dry-Run Mode ---
# Set to "true" for a simulation run (no files downloaded)
# Set to "false" to perform the actual download and processing
DRY_RUN="false"

# --- Create the output directory if it doesn't exist ---
# This ensures log files can be saved there from the start.
mkdir -p "$OUTPUT_PATH"

# Define log and summary file names within the output directory
LOG_FILE="$OUTPUT_PATH/fastq_dump_progress.log"
SUMMARY_FILE="$OUTPUT_PATH/fastq_dump_summary.txt"

# Log script start time
echo "Script started at $(date)" | tee -a "$LOG_FILE"

# Find all directories that start with "SRR" within the base path.
echo "Searching for SRA accessions in $BASE_PATH..." | tee -a "$LOG_FILE"
find "$BASE_PATH" -type d -name "SRR*" | while read -r accession_path; do
    # Extract the SRA accession number (e.g., "SRR13228579") from the path.
    accession=$(basename "$accession_path")
    
    start_time=$(date +%s)
    echo "Processing $accession..." | tee -a "$LOG_FILE"

    if [[ "$DRY_RUN" == "true" ]]; then
        # Dry-run mode: Only print the commands without executing them
        echo "[DRY-RUN] Simulating: prefetch $accession" | tee -a "$LOG_FILE"
        echo "[DRY-RUN] Simulating: fasterq-dump --split-files --include-technical --threads 10 --outdir $OUTPUT_PATH $accession" | tee -a "$LOG_FILE"
    else
        # Full run mode: Execute the actual commands
        echo "Running prefetch for $accession..." | tee -a "$LOG_FILE"
        prefetch "$accession" >> "$LOG_FILE" 2>&1
        
        echo "Running fasterq-dump for $accession to $OUTPUT_PATH..." | tee -a "$LOG_FILE"
        fasterq-dump --threads 18 --split-files --include-technical --outdir "$OUTPUT_PATH" "$accession" >> "$LOG_FILE" 2>&1
    fi
    
    end_time=$(date +%s)
    duration=$((end_time - start_time))
    
    # Create a summary entry
    echo "----------------------------------------" | tee -a "$LOG_FILE"
    echo "Summary for $accession:" | tee -a "$SUMMARY_FILE"
    echo "Status: COMPLETED" | tee -a "$SUMMARY_FILE"
    echo "Duration: $duration seconds" | tee -a "$SUMMARY_FILE"
    echo "Log file: $LOG_FILE" | tee -a "$SUMMARY_FILE"
    echo "Output files saved to: $OUTPUT_PATH" | tee -a "$SUMMARY_FILE"
    echo "----------------------------------------" | tee -a "$LOG_FILE"
done

# Log script completion time
echo "Script finished at $(date)" | tee -a "$LOG_FILE"
echo "All processing complete. Check $SUMMARY_FILE for a summary and $LOG_FILE for detailed output."