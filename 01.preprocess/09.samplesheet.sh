#!/usr/bin/env bash

# bash /mnt/12T/chibao/code/01.preprocess/09.samplesheet.sh \
#   --root /mnt/12T/chibao/data/official_data/fastq_by_tech/bulk_or_plate_RNA/official_cell \
#   --output /mnt/10T2/chibao/gliomas/data_cell/all_projects.csv \
#   --sample-from samn \
#    --projects "PRJNA869964 PRJNA887804 PRJNA887805 PRJNA995768 PRJNA1081384 PRJNA1126248 PRJNA1131103 PRJNA1212512 PRJNA683876" \
#   --split-by-project false



set -Eeuo pipefail

# === CONFIG (defaults) ===
ROOT="/mnt/12T/chibao/data/official_data/fastq_by_tech/bulk_or_plate_RNA/official_cell"
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

  # Find every *_R1_*.fastq.gz and pair to R2
  # Use -print0 for safety; then unique & sorted
  emit_header
  find "$prjdir" -type f -name "*_R1_*.fastq.gz" -print0 \
  | while IFS= read -r -d '' r1; do
      r2="${r1/_R1_/_R2_}"
      if [[ -f "$r2" ]]; then
        sample="$(derive_sample "$r1")"
        printf "%s,%s,%s\n" "$sample" "$r1" "$r2"
      else
        echo "[WARN] Missing R2 for: $r1" >&2
      fi
    done \
  | awk 'NR==1 || !seen[$0]++'  # de-dup rows, keep header once
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
    build_rows_for_project "$prj" > "$outfile"
  done
else
  # single big CSV (stdout or --output FILE)
  if [[ -n "$OUTPUT" ]]; then
    tmp="$(mktemp)"
    # Write header once
    emit_header > "$tmp"
    for prj in "${PROJECTS[@]}"; do
      # append rows skipping header
      build_rows_for_project "$prj" | awk 'NR>1'
    done >> "$tmp"
    # De-dup again, but keep header
    {
      head -n1 "$tmp"
      tail -n +2 "$tmp" | awk '!seen[$0]++'
    } > "$OUTPUT"
    rm -f "$tmp"
    echo "[INFO] Wrote $OUTPUT"
  else
    # print to stdout
    emit_header
    for prj in "${PROJECTS[@]}"; do
      build_rows_for_project "$prj" | awk 'NR>1'
    done | awk 'NR==1 || !seen[$0]++'
  fi
fi
