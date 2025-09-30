# This script will gunzip and rename files after spliting by fasterq-dump for Cell Ranger processing 
process_one() {
  local RUN="$1"
  local SRC_DIR="/mnt/12T/chibao/data/official_data/fastq_by_tech/bulk_or_plate_RNA/official_result/PRJNA961045"
  local DEST_BASE="/mnt/12T/chibao/data/official_data/fastq_by_tech/bulk_or_plate_RNA/official_cell/PRJNA961045"
  local DEST_DIR="${DEST_BASE}/${RUN}"
  mkdir -p "${DEST_DIR}"

  # R1
  pigz -p 9 "${SRC_DIR}/${RUN}_2.fastq"
  mv "${SRC_DIR}/${RUN}_2.fastq.gz" "${DEST_DIR}/${RUN}_S1_L001_R1_001.fastq.gz"

  # R3 (I1)
  pigz -p 9 "${SRC_DIR}/${RUN}_3.fastq"
  mv "${SRC_DIR}/${RUN}_3.fastq.gz" "${DEST_DIR}/${RUN}_S1_L001_R3_001.fastq.gz"

  # R2
  pigz -p 9 "${SRC_DIR}/${RUN}_4.fastq"
  mv "${SRC_DIR}/${RUN}_4.fastq.gz" "${DEST_DIR}/${RUN}_S1_L001_R2_001.fastq.gz"

  # I2
  pigz -p 9 "${SRC_DIR}/${RUN}_1.fastq"
  mv "${SRC_DIR}/${RUN}_1.fastq.gz" "${DEST_DIR}/${RUN}_S1_L001_I2_001.fastq.gz"
}

export -f process_one

# Kick off two runs in parallel:
bash -lc 'process_one SRR24283909' &
bash -lc 'process_one SRR24283910' &
bash -lc 'process_one SRR24283911' &
bash -lc 'process_one SRR24283912' &
bash -lc 'process_one SRR24283917' &
bash -lc 'process_one SRR24283920' &
bash -lc 'process_one SRR24283921' &
bash -lc 'process_one SRR24283923' &


wait
echo "All done."
