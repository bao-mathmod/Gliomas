#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

DIR="/mnt/12T/chibao/data/official_data/fastq_by_tech/bulk_or_plate_RNA/official_result/PRJNA683876"
OUT="read_length_report.txt"

# Create/overwrite the output file with a header
echo -e "filename\tfirst_read_length" > "$OUT"

# Loop through all .fastq or .fastq.gz files in that directory
for f in "$DIR"/*.fastq "$DIR"/*.fastq.gz; do
  # Skip if glob didn't match
  [[ -e "$f" ]] || continue

  # Get first read length
  if [[ "$f" == *.gz ]]; then
    len=$(zcat -- "$f" | sed -n '2p' | awk '{print length($0)}')
  else
    len=$(sed -n '2p' "$f" | awk '{print length($0)}')
  fi

  echo -e "$(basename "$f")\t$len" >> "$OUT"
done

echo "Done. Results saved to $OUT"
