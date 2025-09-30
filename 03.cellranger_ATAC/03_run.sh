#!/usr/bin/env bash
# ============================================================
# 06b_atac_run_from_manifest.sh
# Run cellranger-atac count for each row in the manifest produced by 06a.
# ============================================================

set -Euo pipefail
IFS=$'\n\t'

# --- Optional shared config (will be skipped if missing) ---
CONFIG="/mnt/12T/chibao/code/03.cellranger_ATAC/config.sh"
[[ -s "$CONFIG" ]] && source "$CONFIG" || true

# -------------------- Defaults (edit if no config file) --------------------
: "${MANIFEST:=/mnt/12T/chibao/data/cellranger_data/fastq_cellranger/_cr_atac_samples.tsv}"
: "${REF_PATH:=/mnt/12T/chibao/env_tool/cellranger_atac/refgenome_atac/refdata-cellranger-arc-GRCh38-2024-A}"
: "${OUT_ROOT:=/mnt/12T/chibao/data/cellranger_data/cell_ranger_atac_out}"
: "${MAX_JOBS:=2}"
: "${LOCAL_CORES:=20}"
: "${LOCAL_MEM:=48}"   # GB

ts(){ date '+%Y-%m-%d %H:%M:%S%z'; }
log(){ printf '[%s] %s\n' "$(ts)" "$*" >&2; }

command -v cellranger-atac >/dev/null || { log "ERROR: cellranger-atac not in PATH"; exit 1; }
command -v flock           >/dev/null || { log "ERROR: flock missing (install util-linux)"; exit 1; }
[[ -s "$MANIFEST" ]] || { log "ERROR: manifest not found or empty: $MANIFEST"; exit 1; }
[[ -d "$REF_PATH"  ]] || { log "ERROR: reference not found: $REF_PATH"; exit 1; }
mkdir -p "$OUT_ROOT"

MASTER_LOG="$OUT_ROOT/cellranger_atac_master.log"
exec > >(tee -a "$MASTER_LOG") 2>&1
log "=== START ATAC RUN ==="
log "MANIFEST=$MANIFEST"
log "REF_PATH=$REF_PATH"
log "OUT_ROOT=$OUT_ROOT"
log "MAX_JOBS=$MAX_JOBS localcores=$LOCAL_CORES localmem=$LOCAL_MEM"

RUNS_SUMMARY="$OUT_ROOT/runs_summary.tsv"
FAIL_SUMMARY="$OUT_ROOT/failures.tsv"
_init(){ [[ -f "$1" ]] || printf '%s\n' "$2" > "$1"; }
_init "$RUNS_SUMMARY" "ts\tproject\tsample\tfastq_dir\tstatus\texit_code\tout_dir\tweb_summary\tlog_file"
_init "$FAIL_SUMMARY" "ts\tproject\tsample\tfastq_dir\texit_code\tout_dir\tlog_file\treason"

safe_append(){
  local file="$1"; shift
  local line="$*"
  (
    flock -w 30 9 || { echo "WARN: lock timeout for $file" >&2; exit 0; }
    printf '%s\n' "$line" >> "$file"
  ) 9>>"${file}.lock"
}

extract_reason(){
  local logfile="$1"
  local reason=""
  if [[ -s "$logfile" ]]; then
    reason=$(awk '/Pipestance failed|Log message:/,0 {print}' "$logfile" | head -n 8 | tr '\t' ' ' | paste -sd ' ' -)
    [[ -z "$reason" ]] && reason=$(grep -E "Pipestance failed|ERROR|CRITICAL|No such file|reference|FASTQ|barcode|peaks" "$logfile" | tail -n 1 || true)
  fi
  reason="$(echo "${reason:-unknown}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' | cut -c1-400)"
  printf '%s' "$reason"
}

run_one(){
  local project="$1" sample="$2" fastq_dir="$3"
  local id="$sample"               # change to "${sample}_$(basename "$fastq_dir")" if you want per-run IDs
  local outdir="$OUT_ROOT/$project"
  mkdir -p "$outdir"

  local sample_log="$outdir/${id}.log"
  printf '==== %s START sample=%s ====\n' "$(ts)" "$id" >> "$sample_log"

  # Defensive checks
  if [[ ! -d "$fastq_dir" ]]; then
    log "SKIP (fastq_dir missing): $sample"
    safe_append "$RUNS_SUMMARY" "$(ts)"$'\t'"$project"$'\t'"$sample"$'\t'"$fastq_dir"$'\t'"SKIP_NO_FASTQDIR"$'\t'0$'\t'"$outdir"$'\t'""$'\t'"$sample_log"
    return 0
  fi
  if ! compgen -G "$fastq_dir/${sample}_S1_L001_R1_001.fastq.gz" >/dev/null || \
     ! compgen -G "$fastq_dir/${sample}_S1_L001_R2_001.fastq.gz" >/dev/null || \
     ! compgen -G "$fastq_dir/${sample}_S1_L001_R3_001.fastq.gz" >/dev/null ; then
    log "SKIP (missing R1/R2/R3): $sample -> $fastq_dir"
    safe_append "$RUNS_SUMMARY" "$(ts)"$'\t'"$project"$'\t'"$sample"$'\t'"$fastq_dir"$'\t'"SKIP_MISSING_R1R2R3"$'\t'0$'\t'"$outdir"$'\t'""$'\t'"$sample_log"
    return 0
  fi

  # Already done?
  local ws="$outdir/$id/outs/web_summary.html"
  if [[ -s "$ws" ]]; then
    log "OK (exists): $id (skip)"
    safe_append "$RUNS_SUMMARY" "$(ts)"$'\t'"$project"$'\t'"$sample"$'\t'"$fastq_dir"$'\t'"OK_EXISTS"$'\t'0$'\t'"$outdir/$id/outs"$'\t'"$ws"$'\t'"$sample_log"
    return 0
  fi

  log "RUN: $id  fastqs=$fastq_dir  out=$outdir"

  cmd=( cellranger-atac count
        --id="$id"
        --reference="$REF_PATH"
        --fastqs="$fastq_dir"
        --sample="$sample"
      )
  [[ -n "$LOCAL_CORES" ]] && cmd+=( --localcores="$LOCAL_CORES" )
  [[ -n "$LOCAL_MEM"   ]] && cmd+=( --localmem="$LOCAL_MEM" )

  ( cd "$outdir" && "${cmd[@]}" ) >>"$sample_log" 2>&1
  local rc=$?

  if [[ $rc -eq 0 && -s "$ws" ]]; then
    log "DONE: $id"
    safe_append "$RUNS_SUMMARY" "$(ts)"$'\t'"$project"$'\t'"$sample"$'\t'"$fastq_dir"$'\t'"DONE"$'\t'"$rc"$'\t'"$outdir/$id/outs"$'\t'"$ws"$'\t'"$sample_log"
    printf '==== %s DONE sample=%s (exit=%s) ====\n' "$(ts)" "$id" "$rc" >> "$sample_log"
  else
    log "FAIL ($rc): $id"
    local reason; reason="$(extract_reason "$sample_log")"
    printf '%s\n' "$reason" > "$outdir/$id/FAIL_REASON.txt" || true
    safe_append "$RUNS_SUMMARY" "$(ts)"$'\t'"$project"$'\t'"$sample"$'\t'"$fastq_dir"$'\t'"FAIL"$'\t'"$rc"$'\t'"$outdir/$id/outs"$'\t'"$ws"$'\t'"$sample_log"
    safe_append "$FAIL_SUMMARY" "$(ts)"$'\t'"$project"$'\t'"$sample"$'\t'"$fastq_dir"$'\t'"$rc"$'\t'"$outdir/$id/outs"$'\t'"$sample_log"$'\t'"$reason"
    printf '==== %s FAIL sample=%s (exit=%s) ====\n' "$(ts)" "$id" "$rc" >> "$sample_log"
  fi
}

export -f run_one log ts safe_append extract_reason
export REF_PATH OUT_ROOT RUNS_SUMMARY FAIL_SUMMARY LOCAL_CORES LOCAL_MEM

log "Driving runs from manifest… (MAX_JOBS=$MAX_JOBS)"
if (( MAX_JOBS > 1 )); then
  active=0
  awk -F'\t' 'NR>1{print $1"\t"$2"\t"$3}' "$MANIFEST" \
  | while IFS=$'\t' read -r project sample fastq_dir; do
      [[ -z "${project:-}" || -z "${sample:-}" || -z "${fastq_dir:-}" ]] && continue
      bash -lc 'run_one "$@"' _ "$project" "$sample" "$fastq_dir" &
      (( active+=1 ))
      if (( active >= MAX_JOBS )); then
        wait -n; (( active-=1 ))
      fi
    done
  wait
else
  awk -F'\t' 'NR>1{print $1"\t"$2"\t"$3}' "$MANIFEST" \
  | while IFS=$'\t' read -r project sample fastq_dir; do
      [[ -z "${project:-}" || -z "${sample:-}" || -z "${fastq_dir:-}" ]] && continue
      run_one "$project" "$sample" "$fastq_dir" || true
    done
fi

log "=== ALL DONE ==="
log "Master log:   $MASTER_LOG"
log "Summary TSV:  $RUNS_SUMMARY"
log "Failures TSV: $FAIL_SUMMARY"
