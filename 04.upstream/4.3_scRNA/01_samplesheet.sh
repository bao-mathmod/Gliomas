#!/usr/bin/env bash
set -Euo pipefail
IFS=$'\n\t'

# ============ CONFIG (edit as needed) ==========================
BASE_DIR="/mnt/18T/chibao/gliomas/data/output_cell/scRNA_2/set1"
OUT_DIR="/mnt/18T/chibao/gliomas/data/upstream/scRNA/set1"
OUT_FILE="$OUT_DIR/scrna_manifest_1_enrich_2.tsv"

DEF_GENOME="GRCh38"
DEF_CHEM="10x"
DEF_NOTE="cells"
# ===============================================================

mkdir -p "$OUT_DIR"

# Temporary file before we enforce global uniqueness
TMP_FILE="$(mktemp)"
trap 'rm -f "$TMP_FILE"' EXIT

printf "project_id\tformat\tpath_or_stem\tsample_id\tsample_uid\tgenome\tchemistry\tnote\n" > "$TMP_FILE"

declare -A SEEN_STEM
declare -A SEEN_DIR
declare -A UID_COUNT  # to make.unique globally

sanitize() { sed -E 's/[^A-Za-z0-9._-]+/_/g' <<<"$1"; }

project_from_path() {
  local path="$1"
  local rel="${path#$BASE_DIR/}"
  echo "${rel%%/*}"
}

leaf_token() {
  local p="$1"
  if grep -Eo 'SAMN[0-9]+' <<<"$p" >/dev/null; then
    grep -Eo 'SAMN[0-9]+' <<<"$p" | head -n1
  else
    basename "$(dirname "$p")"
  fi
}

make_unique_uid() {
  local uid="$1"
  if [[ -z "${UID_COUNT[$uid]:-}" ]]; then
    UID_COUNT[$uid]=1
    echo "$uid"
  else
    local n="${UID_COUNT[$uid]}"
    (( n++ ))
    UID_COUNT[$uid]=$n
    echo "${uid}_v${n}"
  fi
}

emit_row() {
  local project="$1" fmt="$2" pathstem="$3" sample="$4"
  local leaf base_uid uniq_uid
  leaf="$(leaf_token "$pathstem")"
  base_uid="$(sanitize "${project}__${sample}__${leaf}")"
  uniq_uid="$(make_unique_uid "$base_uid")"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$project" "$fmt" "$pathstem" "$sample" "$uniq_uid" "$DEF_GENOME" "$DEF_CHEM" "$DEF_NOTE" >> "$TMP_FILE"
}

echo "🔎 Searching for scRNA-seq data in '$BASE_DIR'..."

# --- 1) 10x HDF5 (.h5) ---
while IFS= read -r h5; do
  [[ -z "$h5" ]] && continue
  project="$(project_from_path "$h5")"
  base_fname="$(basename "$h5")"
  if [[ "$base_fname" == "filtered_feature_bc_matrix.h5" ]]; then
    sample=$(basename "$(dirname "$h5")")
  else
    sample="${base_fname%_filtered_feature_bc_matrix.h5}"
    sample="${sample%filtered_feature_bc_matrix.h5}"
  fi
  sample_id="${project}__${sample}"
  emit_row "$project" "10x_h5" "$h5" "$sample_id"
done < <(find "$BASE_DIR" -type f -name "*filtered_feature_bc_matrix.h5" | sort)

# --- 2) .counts.rds ---
while IFS= read -r rds; do
  [[ -z "$rds" ]] && continue
  project="$(project_from_path "$rds")"
  sample="$(basename "$rds" .counts.rds)"
  sample_id="${project}__${sample}"
  emit_row "$project" "rds" "$rds" "$sample_id"
done < <(find "$BASE_DIR" -type f -name "*.counts.rds" | sort)

# --- 3) Wide TSV counts ---
while IFS= read -r tsv; do
  [[ -z "$tsv" ]] && continue
  project="$(project_from_path "$tsv")"
  sample="$(basename "$tsv" .tsv.gz)"
  sample_id="${project}__${sample}"
  emit_row "$project" "tsv_wide" "$tsv" "$sample_id"
done < <(find "$BASE_DIR" -type f \( -name "*countsMatrix.tsv.gz" -o -name "*_raw_counts_*_filtered_cells.tsv.gz" \) | sort)

# --- 4) MTX triplets ---
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

# Finalize
mv "$TMP_FILE" "$OUT_FILE"
echo "✅ Done! Wrote manifest to: $OUT_FILE"

# Quick sanity stats (example numbers printed to terminal)
rows=$(($(wc -l < "$OUT_FILE")-1))
dups_by_sampleid=$(cut -f4 "$OUT_FILE" | tail -n +2 | sort | uniq -d | wc -l)
echo "📊 Rows: $rows | Duplicate sample_id groups: $dups_by_sampleid | All sample_uid are unique."
