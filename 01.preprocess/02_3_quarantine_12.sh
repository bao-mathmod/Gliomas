#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

# Adjust only if your layout differs
SPLIT_ID="split12"
FASTQ_ROOT="/mnt/18T/chibao/gliomas/data/fastq/ena/${SPLIT_ID}"
QUAR="$FASTQ_ROOT/_quarantine"

OFFICIAL_META_ROOT="/mnt/18T/chibao/gliomas/data/metadata/official"
# Per-split enhanced index (appended by your enhancer)
EXPECTED="$OFFICIAL_META_ROOT/${SPLIT_ID}/_expected_files.tsv"

# Your downloader script (with process_one_file)
SCRIPT="/mnt/18T/chibao/gliomas/code/01.preprocess/02_fastq_split9.sh"

# Sanity: EXPECTED must exist
if [[ ! -s "$EXPECTED" ]]; then
  echo "ERROR: Expected index not found or empty: $EXPECTED" >&2
  exit 1
fi

# If no quarantine files, exit early
if ! find "$QUAR" -type f -name "*.fastq.gz" -print -quit | grep -q .; then
  echo "No quarantined FASTQ files found under $QUAR — nothing to do."
  exit 0
fi

# Count files and list a few
find "$QUAR" -type f -name "*.fastq.gz" | tee "$FASTQ_ROOT/_quarantine.list" | wc -l
head -n 10 "$FASTQ_ROOT/_quarantine.list" || true

REQUEUE="$FASTQ_ROOT/requeue_quarantine.tsv"
echo -e "project\tstudy\tsample\texperiment\trun\tfile\turl\tmd5_expected\ttarget_path" > "$REQUEUE"

# Build key map from EXPECTED
# EXPECTED: split_id(1) project_id(2) study(3) sample(4) experiment(5) run(6) file(7) url(8) md5(9) target(10)
awk -F'\t' 'NR>1{print $3"|"$4"|"$5"|"$6"|"$7"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' "$EXPECTED" \
  > "$FASTQ_ROOT/_expected_index.map"

# Build requeue TSV
while read -r f; do
  rel="${f#$QUAR/}"  # study/sample/exp/run/file
  study="${rel%%/*}"; rest="${rel#*/}"
  sample="${rest%%/*}"; rest="${rest#*/}"
  exp="${rest%%/*}";  rest="${rest#*/}"
  run="${rest%%/*}";  file="${rest#*/}"
  key="${study}|${sample}|${exp}|${run}|${file}"

  match=$(awk -F'\t' -v k="$key" '$1==k{print $0; found=1} END{if(!found) exit 1}' "$FASTQ_ROOT/_expected_index.map" || true)
  if [[ -n "$match" ]]; then
    echo "$match" | awk -F'\t' 'BEGIN{OFS="\t"}{print $2,$3,$4,$5,$6,$7,$8,$9,$10}' >> "$REQUEUE"
  else
    echo "WARN: No EXPECTED row for quarantined file: $f" >&2
  fi
done < "$FASTQ_ROOT/_quarantine.list"

echo "Built requeue manifest: $REQUEUE"
if [[ $(wc -l < "$REQUEUE") -le 1 ]]; then
  echo "No rows to retry; exiting."
  exit 0
fi

# Remove quarantined files
cp "$FASTQ_ROOT/_quarantine.list" "$FASTQ_ROOT/_quarantine.list.bak"
if [[ -s "$FASTQ_ROOT/_quarantine.list" ]]; then
  xargs -a "$FASTQ_ROOT/_quarantine.list" -I{} rm -f -- "{}"
fi

# Load only functions
export ENA_LIB_ONLY=1
source "$SCRIPT"

# Safer transfer profile
export CONNS=1
export SPLITS=1
export SPLIT_MIN="16M"
export MAX_RETRIES=6
export BACKOFF_INITIAL=5
export MISMATCH_ACTION="quarantine"  # consider "delete" if you want a stricter cleanup

# Sequential retries
tail -n +2 "$REQUEUE" | while IFS=$'\t' read -r project study sample exp run file url md5exp target; do
  echo "Retrying: $file  ($study/$sample/$exp/$run)"
  process_one_file "$SPLIT_ID" "$project" "$study" "$sample" "$exp" "$run" "$url" "$md5exp"
done

SUMMARY="$FASTQ_ROOT/_download_summary.tsv"

# Success count (project..file)
OK_NOW=$(awk -F'\t' 'NR>1{k=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6; want[k]=1} END{for(k in want) print k}' "$REQUEUE" \
  | grep -F -f - "$SUMMARY" \
  | awk -F'\t' '$9=="success" || $9=="skipped_existing"{c++} END{print c+0}')

TOTAL=$(( $(wc -l < "$REQUEUE") - 1 ))
echo "Requeue results: OK=${OK_NOW} / ${TOTAL}"

# Remaining mismatches (project..file)
grep -F -f <(awk -F'\t' 'NR>1{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' "$REQUEUE") \
  "$SUMMARY" | awk -F'\t' '$9=="md5_mismatch"'

# Optional cleanup
find "$QUAR" -type f -name "*.fastq.gz" -delete
