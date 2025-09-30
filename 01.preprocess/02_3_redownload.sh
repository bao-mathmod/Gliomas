#!/bin/bash
set -Eeuo pipefail
IFS=$'\n\t'

SPLIT_ID="split4"
FASTQ_ROOT="/mnt/12T/chibao/data/official_data/fastq/${SPLIT_ID}"
EXPECTED="/mnt/12T/chibao/data/official_data/metadata/official/${SPLIT_ID}/_expected_files.tsv"
REQUEUE="$FASTQ_ROOT/requeue_specific.tsv"

echo -e "project\tstudy\tsample\texperiment\trun\tfile\turl\tmd5_expected\ttarget_path" > "$REQUEUE"

# helper to append a row from EXPECTED by (study,sample,exp,run,file)
append_from_expected () {
  local study="$1" sample="$2" exp="$3" run="$4" file="$5"
  awk -F'\t' -v s="$study" -v sa="$sample" -v e="$exp" -v r="$run" -v f="$file" \
    'NR>1 && $3==s && $4==sa && $5==e && $6==r && $7==f {print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' \
    "$EXPECTED" >> "$REQUEUE"
}

append_from_expected "PRJNA1081384" "SAMN40176213" "SRX23769605" "SRR28124868" "SRR28124868_1.fastq.gz"
# append_from_expected "PRJNA1081384" "SAMN40176224" "SRX23769594" "SRR28124881" "SRR28124881_2.fastq.gz"
# append_from_expected "PRJNA1081384" "SAMN40176227" "SRX23769591" "SRR28124882" "SRR28124882_2.fastq.gz"
# append_from_expected "PRJNA1081384" "SAMN40176228" "SRX23769590" "SRR28124883" "SRR28124883_1.fastq.gz"

echo "Wrote $(($(wc -l < "$REQUEUE")-1)) rows to $REQUEUE"

# remove partials and old target files for these rows
tail -n +2 "$REQUEUE" | while IFS=$'\t' read -r project study sample exp run file url md5exp target; do
  rm -f -- "$target" "$target.aria2" 2>/dev/null || true
done

SCRIPT="/mnt/12T/chibao/code/01.preprocess/02.fastq_split4.sh"
export ENA_LIB_ONLY=1
source "$SCRIPT"

# safer network profile
export CONNS=1
export SPLITS=1
export SPLIT_MIN="16M"
export MAX_RETRIES=6
export BACKOFF_INITIAL=5
export MISMATCH_ACTION="quarantine"   # or "delete" if you prefer not to keep bad copies

# (optional) if URLs are ftp://, convert to https:// just in REQUEUE
sed -i -E 's#\tftp://#\thttps://#' "$REQUEUE"

# sequential retry
tail -n +2 "$REQUEUE" | while IFS=$'\t' read -r project study sample exp run file url md5exp target; do
  echo "Retrying: $file  ($study/$sample/$exp/$run)"
  process_one_file "$SPLIT_ID" "$project" "$study" "$sample" "$exp" "$run" "$url" "$md5exp"
done

SUMMARY="$FASTQ_ROOT/_download_summary.tsv"

# keys for these requeued rows (project..file)
OK_NOW=$(awk -F'\t' 'NR>1{k=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6; want[k]=1} END{for(k in want) print k}' "$REQUEUE" \
  | grep -F -f - "$SUMMARY" \
  | awk -F'\t' '$9=="success" || $9=="skipped_existing"{c++} END{print c+0}')

TOTAL=$(( $(wc -l < "$REQUEUE") - 1 ))
echo "Requeue results: OK=${OK_NOW} / ${TOTAL}"

echo "Still mismatched:"
grep -F -f <(awk -F'\t' 'NR>1{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' "$REQUEUE") \
  "$SUMMARY" | awk -F'\t' '$9=="md5_mismatch"'
