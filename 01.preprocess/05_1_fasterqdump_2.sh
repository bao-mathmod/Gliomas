#!/usr/bin/env bash
set -euo pipefail

# --- Output directory for FASTQ files ---
OUTPUT_PATH="/mnt/18T/chibao/gliomas/data/fastq/sra/PRJNA1011262"

# --- Log directory ---
LOG_DIR="/mnt/18T/chibao/gliomas/data/fastq/sra/logs/PRJNA1011262"
mkdir -p "$OUTPUT_PATH" "$LOG_DIR"

LOG_FILE="$LOG_DIR/fastq_dump_progress.log"
SUMMARY_FILE="$LOG_DIR/fastq_dump_summary.txt"

# --- Control Dry-Run Mode ---
DRY_RUN="false"
THREADS=18

echo "Script started at $(date)" | tee -a "$LOG_FILE"

# --- Only this specific accession ---
ACCESSIONS=(
SRR25835550
SRR25835551
SRR25835552
SRR25835553
SRR25835554
SRR25835555
SRR25835556
SRR25835557
SRR25835558
SRR25835559
SRR25835560
SRR25835561
SRR25835562
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
