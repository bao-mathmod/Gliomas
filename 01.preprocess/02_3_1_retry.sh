#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

# Adjust only if your layout differs
SPLIT_ID="split4.1"
FASTQ_ROOT="/mnt/18T/chibao/gliomas/data/fastq/ena/${SPLIT_ID}"
QUAR="$FASTQ_ROOT/_quarantine"

OFFICIAL_META_ROOT="/mnt/18T/chibao/gliomas/data/metadata/official"
EXPECTED="$OFFICIAL_META_ROOT/${SPLIT_ID}/_expected_files.tsv"

# Your downloader script (with process_one_file)
SCRIPT="/mnt/18T/chibao/gliomas/code/01.preprocess/02_fastq_split4.sh"
REQUEUE="$FASTQ_ROOT/requeue_quarantine.tsv" # Define REQUEUE here

# Check if the existing manifest is valid before proceeding
if [[ ! -s "$REQUEUE" || $(wc -l < "$REQUEUE") -le 1 ]]; then
  echo "ERROR: Requeue manifest not found or is empty: $REQUEUE" >&2
  echo "Please run the full script first, or ensure the file exists and has data rows." >&2
  exit 1
fi

echo "Using existing requeue manifest: $REQUEUE"

# Load only functions
export ENA_LIB_ONLY=1
source "$SCRIPT"

# Safer transfer profile
export CONNS=1
export SPLITS=1
export SPLIT_MIN="16M"
export MAX_RETRIES=6
export BACKOFF_INITIAL=5
export MISMATCH_ACTION="quarantine"  # Consider "delete" if you want a stricter cleanup

# Sequential retries
tail -n +2 "$REQUEUE" | while IFS=$'\t' read -r project study sample exp run file url md5exp target; do
  echo "Retrying: $file  ($study/$sample/$exp/$run)"
  # The original script does not use $target, but process_one_file will determine it
  process_one_file "$SPLIT_ID" "$project" "$study" "$sample" "$exp" "$run" "$url" "$md5exp"
done

SUMMARY="$FASTQ_ROOT/_download_summary.tsv"

# Success count (project..file)
# ... (rest of the summary logic from the original script) ...
OK_NOW=$(awk -F'\t' 'NR>1{k=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6; want[k]=1} END{for(k in want) print k}' "$REQUEUE" \
  | grep -F -f - "$SUMMARY" \
  | awk -F'\t' '$9=="success" || $9=="skipped_existing"{c++} END{print c+0}')

TOTAL=$(( $(wc -l < "$REQUEUE") - 1 ))
echo "Requeue results: OK=${OK_NOW} / ${TOTAL}"

# Remaining mismatches (project..file)
grep -F -f <(awk -F'\t' 'NR>1{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' "$REQUEUE") \
  "$SUMMARY" | awk -F'\t' '$9=="md5_mismatch"'

# Optional cleanup (This part is less useful now since the files were removed earlier)
# find "$QUAR" -type f -name "*.fastq.gz" -delete