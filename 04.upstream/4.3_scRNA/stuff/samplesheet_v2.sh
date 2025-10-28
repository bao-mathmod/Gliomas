#!/usr/bin/env bash
set -Euo pipefail
IFS=$'\n\t'

# ================== CONFIG ===================
# Root folder that contains PRJNA* subfolders you showed:
BASE_DIR="/mnt/18T/chibao/gliomas/data/output_cell/scRNA"

# Where to write the manifest:
OUT_DIR="/mnt/18T/chibao/gliomas/data/upstream/scRNA/scRNA_clean_2"
OUT_FILE="$OUT_DIR/scrna_manifest.tsv"

# Default annotation (revise if needed)
DEF_GENOME="GRCh38"
DEF_CHEM="10x"
DEF_NOTE="cells"
# =============================================

mkdir -p "$OUT_DIR"

# Header
printf "project_id\tformat\tpath_or_stem\tsample_id\tgenome\tchemistry\tnote\n" > "$OUT_FILE"

# ------------ helpers -------------
project_from_path() {
  local path="$1"; local rel="${path#$BASE_DIR/}"; echo "${rel%%/*}"
}

# climb up from filtered_feature_bc_matrix to find the sample folder (skip 'outs')
sample_from_unprefixed_dir() {
  local dir="$1"   # ends with /filtered_feature_bc_matrix
  local p1="$(basename "$(dirname "$dir")")"   # usually 'outs'
  local p2="$(basename "$(dirname "$(dirname "$dir")")")" # sample
  if [[ "$p1" == "outs" ]]; then echo "$p2"; else echo "$p1"; fi
}

emit_row() {
  local project="$1" fmt="$2" pathstem="$3" sample="$4"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$project" "$fmt" "$pathstem" "$sample" "$DEF_GENOME" "$DEF_CHEM" "$DEF_NOTE" >> "$OUT_FILE"
}

# We’ll choose the highest-priority candidate per (project,sample)
#   Priorities: 10x_h5=100, mtx=80, tsv_wide=50
declare -A CHOSEN_PRIO CHOSEN_FMT CHOSEN_PATH CHOSEN_PROJ CHOSEN_SAMPLE

set_choice() {
  local key="$1" pr="$2" fmt="$3" path="$4" proj="$5" sample="$6"
  local cur="${CHOSEN_PRIO[$key]:--1}"
  if (( pr > cur )); then
    CHOSEN_PRIO[$key]="$pr"
    CHOSEN_FMT[$key]="$fmt"
    CHOSEN_PATH[$key]="$path"
    CHOSEN_PROJ[$key]="$proj"
    CHOSEN_SAMPLE[$key]="$sample"
  fi
}

# ---------- PASS 1: 10x HDF5 ----------
# Matches both "-filtered_feature_bc_matrix.h5" and "_filtered_feature_bc_matrix.h5"
while IFS= read -r h5; do
  [[ -z "$h5" ]] && continue
  project="$(project_from_path "$h5")"
  base="$(basename "$h5")"
  # sample from filename, else fallback to parent folder
  sample="${base%filtered_feature_bc_matrix.h5}"
  sample="${sample%-}"; sample="${sample%_}"
  if [[ -z "$sample" ]]; then sample="$(basename "$(dirname "$h5")")"; fi
  sample_id="${project}__${sample}"
  set_choice "${project}::${sample_id}" 100 "10x_h5" "$h5" "$project" "$sample_id"
done < <(find "$BASE_DIR" -type f -name "*filtered_feature_bc_matrix.h5" | sort)

# ---------- PASS 2: MTX triplets ----------
# Accept BOTH prefixed (*barcodes.tsv[.gz]) and unprefixed (barcodes.tsv[.gz]) cases.
# Also allow features.tsv OR genes.tsv; matrix.mtx(.gz) OR matrix.gz oddball.
while IFS= read -r bfile; do
  [[ -z "$bfile" ]] && continue
  dir="$(dirname "$bfile")"
  base="$(basename "$bfile")"

  # stem is path minus the trailing 'barcodes.tsv(.gz)'
  if [[ "$base" == *"barcodes.tsv.gz" ]]; then
    stem="${bfile%barcodes.tsv.gz}"
  elif [[ "$base" == *"barcodes.tsv" ]]; then
    stem="${bfile%barcodes.tsv}"
  else
    continue
  fi

  # features/genes (either gz or plain)
  feats_gz="${stem}features.tsv.gz"
  feats_pl="${stem}features.tsv"
  genes_gz="${stem}genes.tsv.gz"
  genes_pl="${stem}genes.tsv"

  has_features=false
  for f in "$feats_gz" "$feats_pl" "$genes_gz" "$genes_pl"; do
    [[ -f "$f" ]] && { has_features=true; break; }
  done
  $has_features || continue

  # matrix (try common variants)
  mtx_gz="${stem}matrix.mtx.gz"
  mtx_pl="${stem}matrix.mtx"
  mtx_alt1="${stem%_}matrix.mtx.gz"   # handles missing underscore like ..._LEmatrix.mtx.gz
  mtx_alt2="${stem%_}matrix.gz"       # PRJNA797449 style
  mtx_alt3="${stem}matrix.gz"         # also accept ...matrix.gz with underscore kept

  has_matrix=false
  for m in "$mtx_gz" "$mtx_pl" "$mtx_alt1" "$mtx_alt2" "$mtx_alt3"; do
    [[ -f "$m" ]] && { has_matrix=true; break; }
  done
  $has_matrix || continue

  project="$(project_from_path "$bfile")"

  # Sample ID:
  #   If basename(stem) is empty (unprefixed in Cell Ranger), derive from folder above 'outs/filtered_feature_bc_matrix'
  stem_base="$(basename "$stem")"
  if [[ -z "$stem_base" || "$stem_base" == "/" ]]; then
    # unprefixed case
    sample_core="$(sample_from_unprefixed_dir "$dir")"
  else
    # prefixed case; drop trailing separators if any
    sample_core="${stem_base%_}"; sample_core="${sample_core%.}"
  fi

  sample_id="${project}__${sample_core}"

  # Respect H5 preference: we only add if not already claimed by higher priority
  set_choice "${project}::${sample_id}" 80 "mtx" "$stem" "$project" "$sample_id"

done < <(find "$BASE_DIR" -type f \( -name "*barcodes.tsv" -o -name "*barcodes.tsv.gz" \) | sort)

# ---------- PASS 3: Wide counts tables (TSV/CSV, gz/plain) ----------
# Prioritize raw counts tables; skip obvious metadata files.
# Examples: PRJNA955813 countsMatrix.tsv, PRJNA870065 raw_counts_*filtered_cells.tsv
while IFS= read -r tsv; do
  [[ -z "$tsv" ]] && continue
  # basic guard against metadata files
  case "$(basename "$tsv" | tr '[:upper:]' '[:lower:]')" in
    *metadata*|*cell_metadata*|*meta.csv|*meta.tsv) continue;;
  esac
  project="$(project_from_path "$tsv")"
  # For wide tables spanning multiple samples, leave sample blank; key by path to avoid dedupe
  sample_id=""
  set_choice "${project}::${tsv}" 50 "tsv_wide" "$tsv" "$project" "$sample_id"
done < <(
  find "$BASE_DIR" -type f -iregex '.*\(countsMatrix\.tsv\(\.gz\)?\|raw_counts.*\.\(tsv\|csv\)\(\.gz\)?\|.*all[_-]?samples.*\.\(tsv\|csv\)\(\.gz\)?\|.*filtered_cells.*\.\(tsv\|csv\)\(\.gz\)?\)' \
  | sort
)

# ---------- Emit final manifest ----------
# We’ll output in project+sample order for readability.
{
  for k in "${!CHOSEN_FMT[@]}"; do
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
      "${CHOSEN_PROJ[$k]}" "${CHOSEN_FMT[$k]}" "${CHOSEN_PATH[$k]}" \
      "${CHOSEN_SAMPLE[$k]}" "$DEF_GENOME" "$DEF_CHEM" "$DEF_NOTE"
  done
} | sort -t$'\t' -k1,1 -k4,4 >> "$OUT_FILE"

echo "[OK] Wrote manifest: $OUT_FILE"
