#!/usr/bin/env bash
set -Eeuo pipefail

# This script will rename in 2 scenario 
# + From the output of pigz (step 06 - gunzip) into the logical cell ranger structure 
# + From the output of directly fastq download

# Default config (can be overridden by CLI flags)
ROOT="/mnt/12T/chibao/data/official_data/fastq_by_tech/bulk_or_plate_RNA/official_cell/PRJNA1212512"
DRY_RUN=false                 # set false to actually mv
RENAME_SINGLE_RUN=true      # <-- NEW: rename even if a sample has only 1 SRR?
GLOB_EXT="{fastq,fq}.gz"     # support *.fastq.gz and *.fq.gz

usage() {
  cat <<EOF
Usage:
  $(basename "$0") [--root DIR] [--dry-run true|false] [--rename-single-run true|false]

What it does:
- Renames FASTQs under ROOT/PRJNA.../SAMN.../SRX.../SRR.../
  to Cell Ranger style: <SAMN>_S1_L00X_RY_001.fastq.gz
- X increments per SRR (run) per SAMN.
- Y is R1/R2/R3 (and optionally I1/I2 if present).

Notes:
- If --rename-single-run=false, samples with only 1 SRR are left unchanged.
- If --rename-single-run=true, single-run samples are renamed with lane L001.
EOF
}

# -------- arg parsing --------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --root) ROOT="$2"; shift 2 ;;
    --dry-run) DRY_RUN="$2"; shift 2 ;;
    --rename-single-run) RENAME_SINGLE_RUN="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1"; usage; exit 1 ;;
  esac
done

[[ -d "$ROOT" ]] || { echo "[ERR] ROOT not found: $ROOT" >&2; exit 1; }
shopt -s nullglob

# -------- Helpers --------
read_tag() {
  # Infer read tag: R1/R2/R3 or I1/I2 (index). Fallback R1.
  # Args: basename only
  local b="$1"
  shopt -s nocasematch
  if   [[ "$b" =~ (^|[_-])I1([^0-9]|_|$) ]] || [[ "$b" =~ _I1_ ]] || [[ "$b" =~ _I1\. ]]; then
    echo "I1"
  elif [[ "$b" =~ (^|[_-])I2([^0-9]|_|$) ]] || [[ "$b" =~ _I2_ ]] || [[ "$b" =~ _I2\. ]]; then
    echo "I2"
  elif [[ "$b" =~ (^|[_-])R1([^0-9]|_|$) ]] || [[ "$b" =~ _R1_ ]] || [[ "$b" =~ _R1\. ]] || [[ "$b" =~ _1\.f(ast)?q(\.gz)?$ ]]; then
    echo "R1"
  elif [[ "$b" =~ (^|[_-])R2([^0-9]|_|$) ]] || [[ "$b" =~ _R2_ ]] || [[ "$b" =~ _R2\. ]] || [[ "$b" =~ _2\.f(ast)?q(\.gz)?$ ]]; then
    echo "R2"
  elif [[ "$b" =~ (^|[_-])R3([^0-9]|_|$) ]] || [[ "$b" =~ _R3_ ]] || [[ "$b" =~ _R3\. ]] || [[ "$b" =~ _3\.f(ast)?q(\.gz)?$ ]]; then
    echo "R3"
  else
    echo "R1"
  fi
}

safe_mv() {
  local src="$1" dst="$2"
  if [[ "$DRY_RUN" == true ]]; then
    echo "    [DRY] mv \"$src\" \"$dst\""
  else
    if [[ -e "$dst" ]]; then
      if cmp -s "$src" "$dst"; then
        echo "    [SKIP] exists (same content): $dst"
      else
        echo "    [WARN] destination exists, skipping: $dst" >&2
      fi
    else
      mv -- "$src" "$dst"
      echo "    [OK] $src -> $(basename "$dst")"
    fi
  fi
}

rename_one_srr_dir() {
  local samn="$1" srr_dir="$2" lane="$3"
  printf "  SRR=%s → lane L%03d\n" "$(basename "$srr_dir")" "$lane"

  # enable nullglob so empty globs don't yield the literal pattern
  shopt -s nullglob
  local found_any=false

  # iterate both *.fastq.gz and *.fq.gz
  for fq in "$srr_dir"/*.fastq.gz "$srr_dir"/*.fq.gz; do
    [[ -e "$fq" ]] || continue
    found_any=true

    local base tag ext stem ext2 newname newpath
    base="$(basename "$fq")"
    tag="$(read_tag "$base")"   # R1/R2/R3/I1/I2

    # preserve extension (.fastq.gz or .fq.gz)
    ext="${fq##*.}"             # gz
    stem="${fq%.*}"             # drop .gz
    ext2="${stem##*.}"          # fastq or fq

    newname="${samn}_S1_L$(printf '%03d' "$lane")_${tag}_001.${ext2}.${ext}"
    newpath="${srr_dir}/${newname}"

    # avoid no-op (already named)
    if [[ "$fq" == "$newpath" ]]; then
      echo "    [SKIP] already named: $base"
      continue
    fi

    safe_mv "$fq" "$newpath"
  done

  if [[ "$found_any" == false ]]; then
    echo "    [WARN] no *.fastq.gz or *.fq.gz files found in $srr_dir"
  fi

  # optional: restore default nullglob behavior
  shopt -u nullglob || true
}


# -------- Main --------
# Iterate all SAMN under ROOT (supports many PRJNA trees if ROOT is set above them)
while IFS= read -r -d '' samn_dir; do
  [[ -d "$samn_dir" ]] || continue
  samn="$(basename "$samn_dir")"
  echo ">>> Processing $samn"

  # discover SRR level: ROOT/.../SAMN/*/SRR*
  mapfile -t srr_dirs < <(find "$samn_dir" -mindepth 2 -maxdepth 2 -type d -name 'SRR*' | sort)
  num_srr=${#srr_dirs[@]}

  if (( num_srr == 0 )); then
    echo "  (no SRR folders found)"
    continue
  fi

  if (( num_srr == 1 )); then
    if [[ "$RENAME_SINGLE_RUN" == true ]]; then
      echo "  Single SRR (renaming enabled) — assigning L001"
      rename_one_srr_dir "$samn" "${srr_dirs[0]}" 1
    else
      echo "  Only one SRR run for $samn — skipping rename (use --rename-single-run true to force)"
    fi
    continue
  fi

  # Multiple SRRs: assign lanes L001, L002, ...
  lane=1
  for srr_dir in "${srr_dirs[@]}"; do
    rename_one_srr_dir "$samn" "$srr_dir" "$lane"
    lane=$((lane+1))
  done
done < <(find "$ROOT" -type d -name 'SAMN*' -print0)
