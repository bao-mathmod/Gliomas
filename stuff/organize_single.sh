#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# Config / Arguments
# Usage:
#   ./stage_fastq_for_cellranger.sh \
#      --config pairs.tsv \
#      --dest /mnt/12T/chibao/data/official_data/fastq_by_tech/bulk_or_plate_RNA/official_cell \
#      --threads 8 \
#      --summary /mnt/12T/chibao/data/official_data/fastq_by_tech/bulk_or_plate_RNA/run_stage_summary.tsv \
#      --log /mnt/12T/chibao/data/official_data/fastq_by_tech/bulk_or_plate_RNA/run_stage.log
#
# Required:
#   --config   TSV with columns: run, r1_path, r2_path
#   --dest     Base destination directory where <RUN>/ will be created
#
# Optional (recommended):
#   --threads  Threads for pigz (default: 8)
#   --summary  Path to a TSV summary file (appended; created if not exists)
#   --log      Path to a detailed log file (appended; created if not exists)
###############################################################################

CONFIG=""
DEST_BASE=""
THREADS="8"
SUMMARY=""
LOGFILE=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)  CONFIG="$2"; shift 2;;
    --dest)    DEST_BASE="$2"; shift 2;;
    --threads) THREADS="${2}"; shift 2;;
    --summary) SUMMARY="$2"; shift 2;;
    --log)     LOGFILE="$2"; shift 2;;
    *) echo "Unknown argument: $1"; exit 1;;
  esac
done

if [[ -z "${CONFIG}" || -z "${DEST_BASE}" ]]; then
  echo "ERROR: --config and --dest are required."
  exit 1
fi

if [[ ! -f "${CONFIG}" ]]; then
  echo "ERROR: Config TSV not found: ${CONFIG}"
  exit 1
fi

mkdir -p "${DEST_BASE}"

# Pick compressor
if command -v pigz >/dev/null 2>&1; then
  COMPRESS_CMD=(pigz -p "${THREADS}")
  COMPRESS_NAME="pigz"
else
  COMPRESS_CMD=(gzip)
  COMPRESS_NAME="gzip"
fi

timestamp() { date +"%Y-%m-%d %H:%M:%S%z"; }

log() {
  local msg="[$(timestamp)] $*"
  echo "${msg}"
  if [[ -n "${LOGFILE}" ]]; then
    mkdir -p "$(dirname "${LOGFILE}")"
    echo "${msg}" >> "${LOGFILE}"
  fi
}

# Prepare summary header
if [[ -n "${SUMMARY}" && ! -f "${SUMMARY}" ]]; then
  mkdir -p "$(dirname "${SUMMARY}")"
  echo -e "ts\trun\tr1_src\tr2_src\tdest_dir\tr1_out\tr2_out\tstatus\treason" > "${SUMMARY}"
fi

# Read CONFIG (skip header if present)
# Expected columns: run \t r1_path \t r2_path
tail -n +1 "${CONFIG}" | while IFS=$'\t' read -r RUN R1_SRC R2_SRC; do
  # Skip header line if it contains "run"
  if [[ "${RUN}" == "run" ]]; then
    continue
  fi

  STATUS="OK"
  REASON=""
  DEST_DIR="${DEST_BASE}/${RUN}"
  mkdir -p "${DEST_DIR}"

  # Target output names (Cell Ranger convention)
  R1_OUT="${DEST_DIR}/${RUN}_S1_L001_R1_001.fastq.gz"
  R2_OUT="${DEST_DIR}/${RUN}_S1_L001_R2_001.fastq.gz"

  log "Processing ${RUN} using ${COMPRESS_NAME} -> ${DEST_DIR}"

  # Validate sources
  if [[ ! -f "${R1_SRC}" ]]; then
    STATUS="ERROR"; REASON="R1 source not found"
  fi
  if [[ ! -f "${R2_SRC}" ]]; then
    STATUS="ERROR"; REASON="${REASON:+${REASON}; }R2 source not found"
  fi

  if [[ "${STATUS}" == "ERROR" ]]; then
    log "ERROR: ${RUN} ${REASON}"
    if [[ -n "${SUMMARY}" ]]; then
      echo -e "$(timestamp)\t${RUN}\t${R1_SRC}\t${R2_SRC}\t${DEST_DIR}\t${R1_OUT}\t${R2_OUT}\t${STATUS}\t${REASON}" >> "${SUMMARY}"
    fi
    continue
  fi

  # Compress R1 if needed
  if [[ -f "${R1_OUT}" ]]; then
    log "Skip R1: already exists -> ${R1_OUT}"
  else
    if [[ "${R1_SRC}" == *.gz ]]; then
      # Already compressed -> just copy/rename into place
      log "R1 is already gzipped; copying to ${R1_OUT}"
      cp -n "${R1_SRC}" "${R1_OUT}" || { STATUS="ERROR"; REASON="Failed to copy gz R1"; }
    else
      log "Compressing R1: ${R1_SRC}"
      "${COMPRESS_CMD[@]}" -c "${R1_SRC}" > "${R1_OUT}" || { STATUS="ERROR"; REASON="Compression failed R1"; }
    fi
  fi

  # Compress R2 if needed
  if [[ -f "${R2_OUT}" ]]; then
    log "Skip R2: already exists -> ${R2_OUT}"
  else
    if [[ "${R2_SRC}" == *.gz ]]; then
      log "R2 is already gzipped; copying to ${R2_OUT}"
      cp -n "${R2_SRC}" "${R2_OUT}" || { STATUS="ERROR"; REASON="${REASON:+${REASON}; }Failed to copy gz R2"; }
    else
      log "Compressing R2: ${R2_SRC}"
      "${COMPRESS_CMD[@]}" -c "${R2_SRC}" > "${R2_OUT}" || { STATUS="ERROR"; REASON="${REASON:+${REASON}; }Compression failed R2"; }
    fi
  fi

  if [[ "${STATUS}" == "OK" ]]; then
    log "DONE: ${RUN}"
  else
    log "FAILED: ${RUN} :: ${REASON}"
  fi

  # Append to summary
  if [[ -n "${SUMMARY}" ]]; then
    echo -e "$(timestamp)\t${RUN}\t${R1_SRC}\t${R2_SRC}\t${DEST_DIR}\t${R1_OUT}\t${R2_OUT}\t${STATUS}\t${REASON:-NA}" >> "${SUMMARY}"
  fi

done
