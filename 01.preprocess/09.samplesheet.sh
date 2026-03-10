#!/usr/bin/env bash

# bash /mnt/18T/chibao/gliomas/code/01.preprocess/09.samplesheet.sh \
#   --root /mnt/18T/chibao/gliomas/data_official/00_raw_data_adult/0_fastq/1_fastq_fetch \
#   --projects "PRJNA647809" \
#   --sample-from samn \
#   --output "/mnt/18T/chibao/gliomas/data_official/00_raw_data_adult/1_raw_cell_ranger/new_cohort/new_cohort_samplesheet.csv"

set -Eeuo pipefail

# === CONFIG (defaults) ===
ROOT="/mnt/18T/chibao/gliomas/data/fastq/official"
OUTPUT=""                 # if empty & --split-by-project=false -> prints to stdout
SAMPLE_FROM="samn"        # samn | srr | prefix | dir
SPLIT_BY_PROJECT="false"  # true -> write one CSV per PRJNA under ROOT
PROJECTS=()               # optional filter: space-separated PRJNA IDs

usage() {
  cat <<EOF
Usage:
  $(basename "$0") [--root DIR] [--output FILE] [--sample-from samn|srr|prefix|dir]
                   [--split-by-project true|false]
                   [--projects "PRJNA869964 PRJNA683876"]

Notes:
  * Emits CSV with header: sample,fastq_1,fastq_2
  * Handles multiple fastq naming schemes (_1.fastq.gz, _R1_001.fastq.gz, _R1.fastq.gz)
EOF
}

# --- parse args ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    --root) ROOT="$2"; shift 2 ;;
    --output) OUTPUT="$2"; shift 2 ;;
    --sample-from) SAMPLE_FROM="$2"; shift 2 ;;
    --split-by-project) SPLIT_BY_PROJECT="$2"; shift 2 ;;
    --projects) IFS=' ' read -r -a PROJECTS <<< "$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1"; usage; exit 1 ;;
  esac
done

[[ -d "$ROOT" ]] || { echo "[ERR] ROOT not found: $ROOT" >&2; exit 1; }

emit_header() { echo "sample,fastq_1,fastq_2"; }

# derive sample name from a full path to R1
derive_sample() {
  local r1="$1"
  case "$SAMPLE_FROM" in
    samn)
      # Regex extract SAMN followed by numbers, regardless of folder depth
      if echo "$r1" | grep -qEo 'SAMN[0-9]+'; then
        echo "$r1" | grep -Eo 'SAMN[0-9]+' | head -n 1
      else
        # Fallback to 3 directories up if regex fails
        basename "$(dirname "$(dirname "$(dirname "$r1")")")"
      fi
      ;;
    srr)
      # Extract SRR followed by numbers
      if echo "$r1" | grep -qEo 'SRR[0-9]+'; then
        echo "$r1" | grep -Eo 'SRR[0-9]+' | head -n 1
      else
        basename "$(dirname "$r1")"
      fi
      ;;
    prefix)
      local base; base="$(basename "$r1")"
      printf "%s" "${base%%_*}"
      ;;
    dir)
      # immediate parent directory
      printf "%s" "$(basename "$(dirname "$r1")")"
      ;;
    *)
      local base; base="$(basename "$r1")"
      printf "%s" "${base%%_*}"
      ;;
  esac
}

# build rows for one PRJNA; prints to stdout
build_rows_for_project() {
  local prj="$1"
  local prjdir="$ROOT/$prj"
  [[ -d "$prjdir" ]] || return 0

  # Find multiple variations of R1 files
  find "$prjdir" -type f \( -name "*_R1_001.fastq.gz" -o -name "*_1.fastq.gz" -o -name "*_R1.fastq.gz" \) -print0 \
  | while IFS= read -r -d '' r1; do
      
      # Dynamically determine the R2 file based on the R1 suffix found
      local r2=""
      if [[ "$r1" == *_R1_001.fastq.gz ]]; then
          r2="${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
      elif [[ "$r1" == *_1.fastq.gz ]]; then
          r2="${r1/_1.fastq.gz/_2.fastq.gz}"
      elif [[ "$r1" == *_R1.fastq.gz ]]; then
          r2="${r1/_R1.fastq.gz/_R2.fastq.gz}"
      fi

      if [[ -f "$r2" ]]; then
        sample="$(derive_sample "$r1")"
        printf "%s,%s,%s\n" "$sample" "$r1" "$r2"
      else
        echo "[WARN] Missing R2 for: $r1" >&2
      fi
    done
}

# --- discover projects if not provided ---
if [[ ${#PROJECTS[@]} -eq 0 ]]; then
  mapfile -t PROJECTS < <(find "$ROOT" -mindepth 1 -maxdepth 1 -type d -name "PRJNA*" -printf "%f\n" | sort)
fi
[[ ${#PROJECTS[@]} -gt 0 ]] || { echo "[ERR] No PRJNA* projects under $ROOT" >&2; exit 1; }

if [[ "$SPLIT_BY_PROJECT" == "true" ]]; then
  outdir="${OUTPUT:-.}"
  # if output ends in .csv, take the directory name instead
  [[ "$outdir" == *.csv ]] && outdir=$(dirname "$outdir")
  mkdir -p "$outdir"
  
  for prj in "${PROJECTS[@]}"; do
    outfile="$outdir/${prj}.csv"
    echo "[INFO] Writing $outfile"
    {
        emit_header
        build_rows_for_project "$prj"
    } | awk '!seen[$0]++' > "$outfile"
  done
else
  # single big CSV
  {
      emit_header
      for prj in "${PROJECTS[@]}"; do
        build_rows_for_project "$prj"
      done
  } | awk '!seen[$0]++' > "${OUTPUT:-/dev/stdout}"

  if [[ -n "$OUTPUT" ]]; then
      echo "[INFO] Wrote $OUTPUT"
  fi
fi