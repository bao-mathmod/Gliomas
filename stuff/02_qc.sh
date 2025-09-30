#!/usr/bin/env bash
set -Euo pipefail
IFS=$'\n\t'

CR_FASTQ_ROOT="/mnt/12T/chibao/cell_ranger/fastq_cellranger"
MANIFEST="$CR_FASTQ_ROOT/_cr_samples.tsv"
QC_ROOT="$CR_FASTQ_ROOT/_qc"
LOG="$QC_ROOT/_qc.log"

mkdir -p "$QC_ROOT"
: > "$LOG"
ts(){ date '+%Y-%m-%d %H:%M:%S%z'; }
log(){ printf '[%s] %s\n' "$(ts)" "$*" | tee -a "$LOG" >&2; }

command -v fastqc >/dev/null || { log "ERROR: fastqc not found"; exit 1; }
command -v multiqc >/dev/null || { log "ERROR: multiqc not found"; exit 1; }

[[ -s "$MANIFEST" ]] || { log "ERROR: manifest not found: $MANIFEST"; exit 1; }

# Run FastQC for each FASTQ in the symlink workspace
# Find FASTQs (follow symlinks) and run FastQC
find -L "$CR_FASTQ_ROOT" -type f -name "*.fastq.gz" -not -path "$QC_ROOT/*" -print0 \
| xargs -0 -r -n 8 -P "$(nproc)" fastqc -o "$QC_ROOT" --noextract 2>&1 | tee -a "$LOG"

# Aggregate with MultiQC
multiqc "$QC_ROOT" -o "$QC_ROOT" 2>&1 | tee -a "$LOG"

log "FastQC HTML + MultiQC summary in $QC_ROOT"