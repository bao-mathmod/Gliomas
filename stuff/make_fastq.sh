#!/usr/bin/env bash
set -euo pipefail

# === INPUT: path to the big project folder (containing SRR* subfolders) ===
# Default to your PRJNA869964; you can pass another path as $1 if needed.
PROJECT_DIR="${1:-/mnt/12T/chibao/data/official_data/fastq_by_tech/bulk_or_plate_RNA/official_cell/PRJNA869964}"

# === OUTPUT: TSV file path ===
OUT_TSV="${2:-${PROJECT_DIR}/project_fastq.tsv}"

# Basic checks
if [[ ! -d "$PROJECT_DIR" ]]; then
  echo "ERROR: Project directory not found: $PROJECT_DIR" >&2
  exit 1
fi

PROJECT_NAME="$(basename "$PROJECT_DIR")"

# Write header
echo -e "project\tsample\tfastq_dir" > "$OUT_TSV"

# Iterate only immediate subfolders (each should be an SRR run)
# If your runs aren’t strictly SRR*, remove the -name filter.
while IFS= read -r -d '' RUN_DIR; do
  SAMPLE="$(basename "$RUN_DIR")"

  # Confirm there are FASTQs in this run folder
  # (adjust pattern if your names differ)
  shopt -s nullglob
  FASTQS=("$RUN_DIR"/*.fastq "$RUN_DIR"/*.fastq.gz)
  shopt -u nullglob

  if (( ${#FASTQS[@]} > 0 )); then
    # Record the row: project, sample, fastq_dir
    echo -e "${PROJECT_NAME}\t${SAMPLE}\t${RUN_DIR}" >> "$OUT_TSV"
  fi
done < <(find "$PROJECT_DIR" -mindepth 1 -maxdepth 1 -type d -name 'SRR*' -print0)

echo "Wrote $(($(wc -l < "$OUT_TSV")-1)) rows to: $OUT_TSV"