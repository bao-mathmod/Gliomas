#!/usr/bin/env bash
# ============================================================
# 06a_atac_build_manifest.sh
# Scan linked scATAC fastqs and build a manifest (project, sample, fastq_dir).
# Run this AFTER:
#   bash 05_atac_prep_detect_and_link.sh --do-link
# ============================================================

set -Euo pipefail
IFS=$'\n\t'

# --- Optional shared config (will be skipped if missing) ---
CONFIG="/mnt/12T/chibao/code/03.cellranger_ATAC/config.sh"
[[ -s "$CONFIG" ]] && source "$CONFIG" || true

# -------------------- Defaults (edit if no config file) --------------------
: "${LINKED_ATAC_ROOT:=/mnt/12T/chibao/data/official_data/atac_prep/fastq_cellranger_atac}"
: "${MANIFEST:=/mnt/12T/chibao/data/cellranger_data/fastq_cellranger/_cr_atac_samples.tsv}"

ts(){ date '+%Y-%m-%d %H:%M:%S%z'; }
log(){ printf '[%s] %s\n' "$(ts)" "$*" >&2; }

[[ -d "$LINKED_ATAC_ROOT" ]] || { log "ERROR: linked root not found: $LINKED_ATAC_ROOT"; exit 1; }
mkdir -p "$(dirname "$MANIFEST")"

log "=== BUILD MANIFEST ==="
log "LINKED_ATAC_ROOT=$LINKED_ATAC_ROOT"
log "MANIFEST=$MANIFEST"

# Header (overwrite each time)
echo -e "project\tsample\tfastq_dir" > "$MANIFEST"

# We expect leaves like: <project>/<sample>/<run>/ with Cell Ranger–style names
while IFS= read -r -d '' leaf; do
  rel="${leaf#${LINKED_ATAC_ROOT}/}"
  IFS='/' read -r project sample run <<< "$rel"

  # Require R1/R2/R3; I1 optional
  if ! compgen -G "$leaf/${sample}_S1_L001_R1_001.fastq.gz" >/dev/null || \
     ! compgen -G "$leaf/${sample}_S1_L001_R2_001.fastq.gz" >/dev/null || \
     ! compgen -G "$leaf/${sample}_S1_L001_R3_001.fastq.gz" >/dev/null ; then
    log "SKIP (missing R1/R2/R3): $leaf"
    continue
  fi

  echo -e "${project}\t${sample}\t${leaf}" >> "$MANIFEST"
done < <(find -L "$LINKED_ATAC_ROOT" -name '*_R1_001.fastq.gz' -printf '%h\0' | sort -zu | uniq -z)


rows=$(( $(wc -l <"$MANIFEST") - 1 ))
if (( rows == 0 )); then
  log "ERROR: Manifest empty (no complete R1/R2/R3 sets)."
  exit 2
fi

log "DONE. Manifest rows: $rows → $MANIFEST"
