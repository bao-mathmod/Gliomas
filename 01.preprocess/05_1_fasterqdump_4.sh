#!/usr/bin/env bash
set -euo pipefail

# --- Output directory for FASTQ files ---
OUTPUT_PATH="/mnt/18T/chibao/gliomas/data/fastq/sra/PRJNA1213849"

# --- Log directory ---
LOG_DIR="/mnt/18T/chibao/gliomas/data/fastq/sra/logs/PRJNA1213849"
mkdir -p "$OUTPUT_PATH" "$LOG_DIR"

LOG_FILE="$LOG_DIR/fastq_dump_progress.log"
SUMMARY_FILE="$LOG_DIR/fastq_dump_summary.txt"

# --- Control Dry-Run Mode ---
DRY_RUN="false"
THREADS=18

echo "Script started at $(date)" | tee -a "$LOG_FILE"

# --- Only this specific accession ---
ACCESSIONS=(
SRR32178326
SRR32178327
SRR32178328
SRR32178323
SRR32178334
SRR32178330
SRR32178331
SRR32178322
SRR32178321
SRR32178324
SRR32178325
SRR32178329
SRR32178332
SRR32178333
)

TOTAL=${#ACCESSIONS[@]}
SUCCESS=0
FAIL=0
INDEX=0

for accession in "${ACCESSIONS[@]}"; do
  INDEX=$((INDEX+1))
  start_time=$(date +%s)
  echo "[${INDEX}/${TOTAL}] Processing $accession..." | tee -a "$LOG_FILE"

  if [[ "$DRY_RUN" == "true" ]]; then
    echo "[DRY-RUN] prefetch $accession" | tee -a "$LOG_FILE"
    echo "[DRY-RUN] fasterq-dump --threads $THREADS --split-files --include-technical --outdir \"$OUTPUT_PATH\" $accession" | tee -a "$LOG_FILE"
    run_status="COMPLETED (DRY-RUN)"
    exit_code=0
  else
    set +e
    prefetch "$accession" >> "$LOG_FILE" 2>&1
    rc_prefetch=$?

    if [[ $rc_prefetch -eq 0 ]]; then
      fasterq-dump --threads "$THREADS" --split-files --include-technical \
                   --outdir "$OUTPUT_PATH" "$accession" >> "$LOG_FILE" 2>&1
      rc_dump=$?
    else
      rc_dump=99
    fi
    set -e

    if [[ $rc_prefetch -eq 0 && $rc_dump -eq 0 ]]; then
      run_status="COMPLETED"
      exit_code=0
    else
      run_status="FAILED (prefetch=$rc_prefetch, fasterq-dump=$rc_dump)"
      exit_code=1
    fi
  fi

  end_time=$(date +%s)
  duration=$((end_time - start_time))

  {
    echo "----------------------------------------"
    echo "Summary for $accession:"
    echo "Status: $run_status"
    echo "Duration: $duration seconds"
    echo "Log file: $LOG_FILE"
    echo "Output files saved to: $OUTPUT_PATH"
  } >> "$SUMMARY_FILE"

  if [[ $exit_code -eq 0 ]]; then
    SUCCESS=$((SUCCESS+1))
  else
    FAIL=$((FAIL+1))
  fi

  echo "----------------------------------------" | tee -a "$LOG_FILE"
done

echo "Script finished at $(date)" | tee -a "$LOG_FILE"
echo "Totals: SUCCESS=${SUCCESS}, FAIL=${FAIL}, TOTAL=${TOTAL}" | tee -a "$LOG_FILE"
