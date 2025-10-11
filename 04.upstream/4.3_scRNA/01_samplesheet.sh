#!/usr/bin/env bash
set -Euo pipefail
IFS=$'\n\t'

# ============ CONFIG (edit as needed) ==========================
BASE_DIR="/mnt/18T/chibao/gliomas/data/output_cell/scRNA"
OUT_DIR="/mnt/18T/chibao/gliomas/data/upstream/scRNA/scRNA_clean"
OUT_FILE="$OUT_DIR/scrna_manifest.tsv"

DEF_GENOME="GRCh38"
DEF_CHEM="10x"
DEF_NOTE="cells"
# ===============================================================

mkdir -p "$OUT_DIR"
printf "project_id\tformat\tpath_or_stem\tsample_id\tgenome\tchemistry\tnote\n" > "$OUT_FILE"

declare -A SEEN_STEM
declare -A SEEN_DIR

emit_row() {
  local project="$1" fmt="$2" pathstem="$3" sample="$4"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$project" "$fmt" "$pathstem" "$sample" "$DEF_GENOME" "$DEF_CHEM" "$DEF_NOTE" >> "$OUT_FILE"
}

project_from_path() {
  local path="$1"
  local rel="${path#$BASE_DIR/}"
  echo "${rel%%/*}"
}

echo "🔎 Searching for scRNA-seq data in '$BASE_DIR'..."

# --- 1) 10x HDF5 (.h5) - CORRECTED LOGIC ---
while IFS= read -r h5; do
  [[ -z "$h5" ]] && continue
  project="$(project_from_path "$h5")"
  base_fname="$(basename "$h5")"

  # Check if the filename is generic or has a prefix
  if [[ "$base_fname" == "filtered_feature_bc_matrix.h5" ]]; then
    # Case 1: Generic filename, sample ID is the parent directory name
    sample=$(basename "$(dirname "$h5")")
  else
    # Case 2: Sample ID is a prefix in the filename
    sample="${base_fname%_filtered_feature_bc_matrix.h5}"
    sample="${sample%filtered_feature_bc_matrix.h5}"
  fi
  
  sample_id="${project}__${sample}"
  emit_row "$project" "10x_h5" "$h5" "$sample_id"
done < <(find "$BASE_DIR" -type f -name "*filtered_feature_bc_matrix.h5" | sort)

# --- 2) R Data Structure (.rds count files) ---
while IFS= read -r rds; do
  [[ -z "$rds" ]] && continue
  project="$(project_from_path "$rds")"
  sample="$(basename "$rds" .counts.rds)"
  sample_id="${project}__${sample}"
  emit_row "$project" "rds" "$rds" "$sample_id"
done < <(find "$BASE_DIR" -type f -name "*.counts.rds" | sort)

# --- 3) Wide TSV counts (priority cases like PRJNA955813) ---
while IFS= read -r tsv; do
  [[ -z "$tsv" ]] && continue
  project="$(project_from_path "$tsv")"
  sample="$(basename "$tsv" .tsv.gz)"
  sample_id="${project}__${sample}"
  emit_row "$project" "tsv_wide" "$tsv" "$sample_id"
done < <(find "$BASE_DIR" -type f \( -name "*countsMatrix.tsv.gz" -o -name "*_raw_counts_*_filtered_cells.tsv.gz" \) | sort)

# --- 4) Matrix Market (MTX) triplets - handles two types ---
while IFS= read -r bfile; do
  [[ -z "$bfile" ]] && continue
  dir="$(dirname "$bfile")"
  base="$(basename "$bfile")"
  project="$(project_from_path "$bfile")"

  if [[ "$base" == "barcodes.tsv.gz" ]]; then
    if [[ -f "$dir/features.tsv.gz" && -f "$dir/matrix.mtx.gz" ]]; then
      [[ -n "${SEEN_DIR[$dir]:-}" ]] && continue
      SEEN_DIR[$dir]=1
      sample="$(basename "$dir")"
      sample_id="${project}__${sample}"
      emit_row "$project" "mtx_dir" "$dir" "$sample_id"
    fi
  else
    if [[ "$base" == *_barcodes.tsv.gz ]]; then
      stem="${bfile%_barcodes.tsv.gz}"
    elif [[ "$base" == *_barcodes.tsv ]]; then
      stem="${bfile%_barcodes.tsv}"
    else
      continue
    fi
    if { [[ -f "${stem}_features.tsv.gz" ]] || [[ -f "${stem}_features.tsv" ]]; } && \
       { [[ -f "${stem}_matrix.mtx.gz" ]] || [[ -f "${stem}_matrix.mtx" ]]; }; then
      [[ -n "${SEEN_STEM[$stem]:-}" ]] && continue
      SEEN_STEM[$stem]=1
      sample="$(basename "$stem")"
      sample_id="${project}__${sample}"
      emit_row "$project" "mtx_stem" "$stem" "$sample_id"
    fi
  fi
done < <(find "$BASE_DIR" -type f \( -name "barcodes.tsv.gz" -o -name "*_barcodes.tsv.gz" -o -name "*_barcodes.tsv" \) | sort)

echo "✅ Done! Wrote manifest to: $OUT_FILE"