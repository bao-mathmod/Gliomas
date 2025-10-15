#!/usr/bin/env bash

# bash /mnt/18T/chibao/gliomas/code/01.preprocess/09.samplesheet.sh \
#   --root /mnt/18T/chibao/gliomas/data/fastq/official \
#   --output /mnt/18T/chibao/gliomas/data/fastq/official/all_projects.csv \
#   --sample-from samn \
#    --projects "PRJNA1125010 PRJNA1134206 PRJNA1141154 PRJNA797449 PRJNA968165" \
#   --split-by-project false

# PRJNA1241125 PRJNA1213849 PRJNA1098411 

set -Eeuo pipefail

# === CONFIG (defaults) ===
ROOT="/mnt/18T/chibao/gliomas/data/fastq/official"
OUTPUT=""                 # if empty & --split-by-project=false -> prints to stdout
SAMPLE_FROM="samn"        # samn | prefix | dir
SPLIT_BY_PROJECT="false"  # true -> write one CSV per PRJNA under ROOT
PROJECTS=()               # optional filter: space-separated PRJNA IDs

usage() {
  cat <<EOF
Usage:
  $(basename "$0") [--root DIR] [--output FILE] [--sample-from samn|prefix|dir]
                   [--split-by-project true|false]
                   [--projects "PRJNA869964 PRJNA683876"]

Notes:
  * Assumes flat structure: ROOT/PRJNA*/SAMN*/SRX*/SRR*/*.fastq.gz
  * Emits CSV with header: sample,fastq_1,fastq_2
  * Groups are handled by nf-core/scrnaseq automatically via 'sample' column
    (so multiple rows per sample = multiple runs/lanes merged by cellranger)
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
      # .../PRJNAxxxxxx/SAMNxxxxxx/SRX.../SRR.../file
      # go up: file -> SRR -> SRX -> SAMN
      local samn
      samn="$(basename "$(dirname "$(dirname "$(dirname "$r1")")")")"
      printf "%s" "$samn"
      ;;
    prefix)
      local base; base="$(basename "$r1")"
      printf "%s" "${base%%_*}"
      ;;
    dir)
      # immediate parent directory (SRR)
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

  # Find every *_R1_001.fastq.gz and pair to its corresponding R2 file
  # Use -print0 for safety; then unique & sorted
  find "$prjdir" -type f -name "*_R1_001.fastq.gz" -print0 \
  | while IFS= read -r -d '' r1; do
      # More specific replacement to guarantee finding the R2 file
      r2="${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
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
  # one CSV per project, named <PRJNA>.csv next to the script (or set --output as a directory)
  outdir="${OUTPUT:-.}"
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
  # single big CSV (stdout or --output FILE)
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