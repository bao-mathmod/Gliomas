#!/usr/bin/env bash
set -euo pipefail

# Usage: Use for file downloading manual not from SRA tools
#   ./rename_to_cellranger.sh /mnt/10T/chibao/gliomas/fastq/split4
#   bash /mnt/12T/chibao/code/01.preprocess/06_4_naming.sh --dry-run /mnt/12T/chibao/data/official_data/fastq/split4/PRJNA1081384
#
# Notes:
# - Renames files to: <SAMN_ID>_S1_L001_R1_001.fastq.gz etc.
# - If no SAMN ID found in path, falls back to SRR_ID from filename.
# - Lane is fixed to L001 and chunk is fixed to 001 (typical for public SRA runs).
# - Writes/append to: rename_cellranger_summary.tsv at the base directory.

DRY_RUN=0
if [[ "${1:-}" == "--dry-run" ]]; then
  DRY_RUN=1
  shift
fi

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 [--dry-run] <BASE_DIR>" >&2
  exit 1
fi

BASE_DIR="$1"
if [[ ! -d "$BASE_DIR" ]]; then
  echo "Error: BASE_DIR not found: $BASE_DIR" >&2
  exit 1
fi

LOG_FILE="$BASE_DIR/rename_cellranger_summary.tsv"
# Create log header once if file not exists
if [[ ! -f "$LOG_FILE" ]]; then
  echo -e "timestamp\tstatus\told_path\tnew_path\tmessage" > "$LOG_FILE"
fi

timestamp() { date +"%Y-%m-%d %H:%M:%S"; }

# Helper: write a log line
log() {
  local status="$1"; shift
  local old="$1"; shift
  local new="$1"; shift
  local msg="$1"; shift || true
  echo -e "$(timestamp)\t${status}\t${old}\t${new}\t${msg}" >> "$LOG_FILE"
}

# A regex that looks like a Cell Ranger name already (to avoid double renaming):
#   <SAMPLE>_S<NUM>_L<3 digits>_(R1|R2|I1)_<3 digits>.fastq.gz
is_cellranger_like() {
  local name="$1"
  [[ "$name" =~ _S[0-9]+_L[0-9]{3}_(R1|R2|I1)_[0-9]{3}\.fastq\.gz$ ]]
}

# Map raw read suffix to Cell Ranger read tag
read_tag() {
  case "$1" in
    1) echo "R1" ;;
    2) echo "R2" ;;
    3) echo "I1" ;; # optional index read if provided as _3
    *) echo "" ;;
  esac
}

# NEW: Helper function to extract SAMN ID from a path
get_samn_id() {
  local path="$1"
  # Search for a directory component that matches SAMN<digits>
  if [[ "$path" =~ /(SAMN[0-9]+)/ ]]; then
    echo "${BASH_REMATCH[1]}"
  else
    echo ""
  fi
}

# Walk and process only *.fastq.gz under BASE_DIR
# We specifically match files ending with _1.fastq.gz, _2.fastq.gz, or _3.fastq.gz
# and skip anything already in CR naming.
while IFS= read -r -d '' f; do
  bn="$(basename -- "$f")"
  dir="$(dirname -- "$f")"

  # Skip if already cellranger-like
  if is_cellranger_like "$bn"; then
    log "SKIP-ALREADY-FORMATTED" "$f" "" "Already Cell Ranger style"
    continue
  fi

  # Expect form: PREFIX_1.fastq.gz (e.g., SRR12345678_1.fastq.gz)
  if [[ "$bn" =~ ^(.+)_([123])\.fastq\.gz$ ]]; then
    file_prefix="${BASH_REMATCH[1]}"
    read_num="${BASH_REMATCH[2]}"
  else
    # Not matching expected pattern; skip
    log "SKIP-NONMATCH" "$f" "" "Name does not match '<prefix>_[1-3].fastq.gz'"
    continue
  fi

  # Determine the sample prefix to use for the new filename
  # NEW: Prioritize SAMN ID from path
  sample_prefix=$(get_samn_id "$f")
  if [[ -z "$sample_prefix" ]]; then
    # Fallback: use the original filename prefix (e.g., SRR ID)
    sample_prefix="$file_prefix"
  fi

  tag="$(read_tag "$read_num")"
  if [[ -z "$tag" ]]; then
    log "SKIP-READTAG" "$f" "" "Unknown read index: $read_num"
    continue
  fi

  # Build new name: <PREFIX>_S1_L001_<TAG>_001.fastq.gz
  new_bn="${sample_prefix}_S1_L001_${tag}_001.fastq.gz"
  new_path="${dir}/${new_bn}"

  # If target exists, do not overwrite
  if [[ -e "$new_path" ]]; then
    log "SKIP-EXISTS" "$f" "$new_path" "Target already exists"
    continue
  fi

  if [[ $DRY_RUN -eq 1 ]]; then
    echo "[DRY-RUN] mv -- \"$f\" \"$new_path\""
    log "DRY-RUN" "$f" "$new_path" "Would rename"
  else
    mv -- "$f" "$new_path"
    log "RENAMED" "$f" "$new_path" "OK"
  fi

done < <(find "$BASE_DIR" -type f -name "*.fastq.gz" -print0)