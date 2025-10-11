# #!/bin/bash

# # Define an array of directories to process
# DIRECTORIES=(
#     "/mnt/12T/chibao/data/official_data/fastq_by_tech/scATAC/PRJNA578617/SAMN13066825/SRX7026892/SRR10315838"
#     "/mnt/12T/chibao/data/official_data/fastq_by_tech/scATAC/PRJNA578617/SAMN13066826/SRX7026891/SRR10315837"
#     "/mnt/12T/chibao/data/official_data/fastq_by_tech/scATAC/PRJNA578617/SAMN13066827/SRX7026890/SRR10315836"
#     "/mnt/12T/chibao/data/official_data/fastq_by_tech/scATAC/PRJNA578617/SAMN13066828/SRX7026889/SRR10315835"
# )

# # Define the output file name
# OUTPUT_FILE="fastq_headers.txt"

# # Clear the output file from previous runs
# > "$OUTPUT_FILE"

# # Loop through each directory
# for dir in "${DIRECTORIES[@]}"; do
#     echo "========================================" >> "$OUTPUT_FILE"
#     echo "Processing directory: $dir" >> "$OUTPUT_FILE"
#     echo "========================================" >> "$OUTPUT_FILE"
    
#     # Loop through common FASTQ suffixes
#     for suffix in "_1.fastq.gz" "_2.fastq.gz" "_3.fastq.gz" ".fastq.gz"; do
#         file="$dir/$(basename "$dir")$suffix"

#         # Check if the file exists before trying to process it
#         if [[ -f "$file" ]]; then
#             echo "--- Displaying head of: $(basename "$file") ---" >> "$OUTPUT_FILE"
#             zcat "$file" | sed -n '1,8p' >> "$OUTPUT_FILE"
#             echo "" >> "$OUTPUT_FILE"
#         fi
#     done
# done

# echo "Script complete. All output is in $OUTPUT_FILE."

# This script will gunzip and rename files after spliting by fasterq-dump for Cell Ranger processing 
process_one() {
  local RUN="$1"
  local SRC_DIR="/mnt/18T/chibao/gliomas/data/fastq/sra/PRJNA1213849"
  local DEST_BASE="/mnt/18T/chibao/gliomas/data/fastq/sra/PRJNA1213849"
  local DEST_DIR="${DEST_BASE}/${RUN}"
  mkdir -p "${DEST_DIR}"

  # R1
  pigz -p 14 "${SRC_DIR}/${RUN}_3.fastq"
  mv "${SRC_DIR}/${RUN}_3.fastq.gz" "${DEST_DIR}/${RUN}_S1_L001_R1_001.fastq.gz"

  # R2
  pigz -p 14 "${SRC_DIR}/${RUN}_4.fastq"
  mv "${SRC_DIR}/${RUN}_4.fastq.gz" "${DEST_DIR}/${RUN}_S1_L001_R2_001.fastq.gz"
}

export -f process_one

# Kick off two runs in parallel:
bash -lc 'process_one SRR32178321' &
bash -lc 'process_one SRR32178322' &
bash -lc 'process_one SRR32178323' &
bash -lc 'process_one SRR32178324' &
bash -lc 'process_one SRR32178325' &
bash -lc 'process_one SRR32178326' &
bash -lc 'process_one SRR32178327' &
bash -lc 'process_one SRR32178328' &
bash -lc 'process_one SRR32178329' &
bash -lc 'process_one SRR32178330' &
bash -lc 'process_one SRR32178331' &
bash -lc 'process_one SRR32178332' &
bash -lc 'process_one SRR32178333' &
bash -lc 'process_one SRR32178334' &


wait
echo "All done."
