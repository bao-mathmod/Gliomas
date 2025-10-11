# This script will create a samplesheet for processing snRNA data in the tsv format
# This Bash script detects all three formats:
      # 10x HDF5 (*filtered_feature_bc_matrix.h5) → format=10x_h5
      # Matrix Market triplets (prefix + _barcodes.tsv[.gz], _features.tsv[.gz], _matrix.mtx[.gz]) → format=mtx with a stem that our R loader expects
      # Wide TSV (the big all-samples counts table) → format=tsv_wide

#!/usr/bin/env bash
set -Euo pipefail
IFS=$'\n\t'

# ============ CONFIG (edit as needed) ==========================
# Root folder that contains PRJNA* subfolders you showed:
BASE_DIR="/mnt/18T/chibao/gliomas/data/output_cell/snRNA"

# Where to write the manifest:
OUT_DIR="$BASE_DIR/../snRNA_clean"
OUT_FILE="$OUT_DIR/snrna_manifest.tsv"

# Default annotation (you can revise later)
DEF_GENOME="GRCh38"
DEF_CHEM="10x"
DEF_NOTE="nuclei"
# ===============================================================

mkdir -p "$OUT_DIR"

# Header
printf "project_id\tformat\tpath_or_stem\tsample_id\tgenome\tchemistry\tnote\n" > "$OUT_FILE"

# Track stems already written to avoid duplicates if both .gz and plain exist
declare -A SEEN_STEM

# --- helper: print one TSV row ---
emit_row() {
  local project="$1" fmt="$2" pathstem="$3" sample="$4"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$project" "$fmt" "$pathstem" "$sample" "$DEF_GENOME" "$DEF_CHEM" "$DEF_NOTE" >> "$OUT_FILE"
}

# --- helper: derive project id from absolute path under BASE_DIR ---
project_from_path() {
  local path="$1"
  local rel="${path#$BASE_DIR/}"
  echo "${rel%%/*}"  # first path component
}

# --- 1) 10x HDF5 (.h5 filtered matrices) ---
# Matches both "-filtered_feature_bc_matrix.h5" and "_filtered_feature_bc_matrix.h5"
while IFS= read -r h5; do
  [[ -z "$h5" ]] && continue
  project="$(project_from_path "$h5")"
  base="$(basename "$h5")"
  # sample id: strip the common suffix; also remove trailing ".h5"
  sample="${base%filtered_feature_bc_matrix.h5}"
  # remove any trailing "-" or "_" left behind
  sample="${sample%-}"
  sample="${sample%_}"
  # prepend project
  sample_id="${project}__${sample}"
  emit_row "$project" "10x_h5" "$h5" "$sample_id"
done < <(find "$BASE_DIR" -type f -name "*filtered_feature_bc_matrix.h5" | sort)

# --- 2) Wide TSV counts (one big matrix) ---
while IFS= read -r tsv; do
  [[ -z "$tsv" ]] && continue
  project="$(project_from_path "$tsv")"
  # sample_id left blank; we will split by column prefixes later in R
  emit_row "$project" "tsv_wide" "$tsv" ""
done < <(find "$BASE_DIR" -type f -regex ".*_snRNAseq_counts_mtx_allsamples\.tsv(\.gz)?$" | sort)

# --- 3) Matrix Market triplets (prefix_ + barcodes/features/matrix) ---
# Strategy:
#   - Find every *barcodes.tsv[.gz]
#   - Build its "stem" by stripping the trailing 'barcodes.tsv' or 'barcodes.tsv.gz'
#   - Check that matching features + matrix exist (plain OR gz)
#   - Emit one row per stem (de-duplicated)
while IFS= read -r bfile; do
  [[ -z "$bfile" ]] && continue
  dir="$(dirname "$bfile")"
  base="$(basename "$bfile")"

  # Accept both suffixes
  if [[ "$base" == *_barcodes.tsv.gz ]]; then
    stem="${bfile%barcodes.tsv.gz}"
  elif [[ "$base" == *_barcodes.tsv ]]; then
    stem="${bfile%barcodes.tsv}"
  else
    continue
  fi

  # Check corresponding features + matrix exist in either gz or plain
  feats_gz="${stem}features.tsv.gz"
  feats_pl="${stem}features.tsv"
  mtx_gz="${stem}matrix.mtx.gz"
  mtx_pl="${stem}matrix.mtx"

  if { [[ -f "$feats_gz" ]] || [[ -f "$feats_pl" ]]; } && { [[ -f "$mtx_gz" ]] || [[ -f "$mtx_pl" ]]; }; then
    # Avoid duplicates if we see both gz and plain
    key="$stem"
    if [[ -n "${SEEN_STEM[$key]:-}" ]]; then
      continue
    fi
    SEEN_STEM[$key]=1

    project="$(project_from_path "$bfile")"

    # sample_id is project + basename(stem) with trailing '_' trimmed
    stem_base="$(basename "$stem")"
    sample_core="${stem_base%_}"        # drop trailing underscore if present
    sample_id="${project}__${sample_core}"

    # path_or_stem should be EXACTLY the stem we just computed
    emit_row "$project" "mtx" "$stem" "$sample_id"
  fi
done < <(find "$BASE_DIR" -type f -name "*_barcodes.tsv" -o -name "*_barcodes.tsv.gz" | sort)

echo "[OK] Wrote manifest: $OUT_FILE"