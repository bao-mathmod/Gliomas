# #!/usr/bin/env bash
# set -euo pipefail
# IFS=$'\n\t'

# # ===== Required paths (edit to your real paths) =====
# MANIFEST="/mnt/rdisk/gliomas/data/official_data/fastq_cellranger/_cr_samples.tsv"  # project_id  sample  fastq_dir
# TRANSCRIPTOME="/mnt/rdisk/gliomas/data/cellranger/ref_genome/refdata-gex-GRCh38-2020-A"
# OUT_ROOT="/mnt/rdisk/gliomas/data/cellranger_out"   # where per-sample ID folders will be created

# # ===== Options (simple) =====
# CREATE_BAM="false"          # set to "true" if you also want BAMs (v8+); CRAM is default
# INCLUDE_INTRONS=""          # set to "true" for snRNA/pre-mRNA, else leave empty
# CHEMISTRY=""                # e.g. "SC3Pv3", "auto"; leave empty to let cellranger decide
# EXPECT_CELLS=""             # e.g. "6000"; leave empty to let cellranger estimate
# MAX_JOBS=1                  # 1 = sequential; increase (e.g., 2, 3) for light parallel

# # ===== Tiny helpers =====
# ts(){ date '+%Y-%m-%d %H:%M:%S%z'; }
# log(){ printf '[%s] %s\n' "$(ts)" "$*" >&2; }

# # ===== Sanity checks =====
# command -v cellranger >/dev/null || { log "ERROR: cellranger not found in PATH"; exit 1; }
# [[ -d "$TRANSCRIPTOME" ]] || { log "ERROR: transcriptome not found: $TRANSCRIPTOME"; exit 1; }
# [[ -s "$MANIFEST" ]] || { log "ERROR: manifest not found or empty: $MANIFEST"; exit 1; }
# mkdir -p "$OUT_ROOT"

# # ===== Runner for one line from manifest =====
# run_one() {
#   local project="$1" sample="$2" fastq_dir="$3"
#   local id="$sample"                           # keep identical to your working single-sample command
#   local outdir="$OUT_ROOT/$project"
#   mkdir -p "$outdir"

#   # quick check the fastqs exist
#   if [[ ! -d "$fastq_dir" ]]; then
#     log "SKIP (fastq_dir missing): $sample -> $fastq_dir"
#     return 0
#   fi
#   if ! compgen -G "$fastq_dir/*_R1_001.fastq.gz" >/dev/null || \
#      ! compgen -G "$fastq_dir/*_R2_001.fastq.gz" >/dev/null; then
#     log "SKIP (no R1/R2 pattern in dir): $sample -> $fastq_dir"
#     return 0
#   fi

#   # skip if already done (web_summary exists)
#   if [[ -s "$outdir/$id/outs/web_summary.html" ]]; then
#     log "OK (exists, skipping): $id"
#     return 0
#   fi

#   log "START: $id  fastqs=$fastq_dir  out=$outdir"

#   # build the exact command you used, plus optional knobs if set
#   cmd=( cellranger count
#         --id="$id"
#         --transcriptome="$TRANSCRIPTOME"
#         --fastqs="$fastq_dir"
#         --sample="$sample"
#         --create-bam="$CREATE_BAM"
#         --localcores=45
#         --localmem=440
#       )
#   [[ -n "$INCLUDE_INTRONS" ]] && cmd+=( --include-introns="$INCLUDE_INTRONS" )
#   [[ -n "$CHEMISTRY"      ]] && cmd+=( --chemistry="$CHEMISTRY" )
#   [[ -n "$EXPECT_CELLS"   ]] && cmd+=( --expect-cells="$EXPECT_CELLS" )

#   ( cd "$outdir" && "${cmd[@]}" )
#   local rc=$?

#   if [[ $rc -eq 0 ]]; then
#     log "DONE: $id → $outdir/$id/outs"
#   else
#     log "FAIL ($rc): $id"
#   fi
#   return $rc
# }

# export -f run_one log ts
# export TRANSCRIPTOME OUT_ROOT CREATE_BAM INCLUDE_INTRONS CHEMISTRY EXPECT_CELLS

# # ===== Read manifest and run =====
# if (( MAX_JOBS > 1 )); then
#   # Light parallel using bash job control
#   active=0
#   while IFS=$'\t' read -r project sample fastq_dir; do
#     [[ -z "${project:-}" || -z "${sample:-}" || -z "${fastq_dir:-}" ]] && continue
#     bash -lc 'run_one "$@"' _ "$project" "$sample" "$fastq_dir" &
#     (( active+=1 ))
#     if (( active >= MAX_JOBS )); then
#       wait -n
#       (( active-=1 ))
#     fi
#   done < <(awk -F'\t' 'NR>1{print $1"\t"$2"\t"$3}' "$MANIFEST")
#   wait
# else
#   # Sequential (safest)
#   awk -F'\t' 'NR>1{print $1"\t"$2"\t"$3}' "$MANIFEST" \
#   | while IFS=$'\t' read -r project sample fastq_dir; do
#       [[ -z "${project:-}" || -z "${sample:-}" || -z "${fastq_dir:-}" ]] && continue
#       run_one "$project" "$sample" "$fastq_dir"
#     done
# fi

# log "All done."

# #!/usr/bin/env bash
# set -euo pipefail
# IFS=$'\n\t'

# # ===== Required paths (edit to your real paths) =====
# MANIFEST="/mnt/rdisk/gliomas/data/official_data/fastq_cellranger/_cr_samples.tsv"  # project_id  sample  fastq_dir
# TRANSCRIPTOME="/mnt/rdisk/gliomas/data/cellranger/ref_genome/refdata-gex-GRCh38-2020-A"
# OUT_ROOT="/mnt/rdisk/gliomas/data/cellranger_out"   # where per-sample ID folders will be created

# # ===== Options (simple) =====
# CREATE_BAM="false"          # set to "true" if you also want BAMs (v8+); CRAM is default
# INCLUDE_INTRONS=""          # set to "true" for snRNA/pre-mRNA, else leave empty
# CHEMISTRY=""                # e.g. "SC3Pv3", "auto"; leave empty to let cellranger decide
# EXPECT_CELLS=""             # e.g. "6000"; leave empty to let cellranger estimate
# MAX_JOBS=1                  # 1 = sequential; increase (e.g., 2, 3) for light parallel

# # (optional) local cores/mem; set to "" to omit
# LOCAL_CORES="45"
# LOCAL_MEM="450"   # GB

# # ===== Tiny helpers =====
# ts(){ date '+%Y-%m-%d %H:%M:%S%z'; }
# log(){ printf '[%s] %s\n' "$(ts)" "$*" >&2; }

# # ===== Sanity checks =====
# command -v cellranger >/dev/null || { log "ERROR: cellranger not found in PATH"; exit 1; }
# command -v flock       >/dev/null || { log "ERROR: flock not found (install util-linux)"; exit 1; }
# [[ -d "$TRANSCRIPTOME" ]] || { log "ERROR: transcriptome not found: $TRANSCRIPTOME"; exit 1; }
# [[ -s "$MANIFEST" ]] || { log "ERROR: manifest not found or empty: $MANIFEST"; exit 1; }
# mkdir -p "$OUT_ROOT"

# # ===== Global master log (high-level) =====
# LOGFILE="$OUT_ROOT/cellranger_master.log"
# # tee the whole script's stdout+stderr to master log + console
# exec > >(tee -a "$LOGFILE") 2>&1
# log "=== START run === MANIFEST=$MANIFEST  OUT_ROOT=$OUT_ROOT  TRANSCRIPTOME=$TRANSCRIPTOME"

# # ===== Summary tables (parallel-safe with flock) =====
# RUNS_SUMMARY="$OUT_ROOT/runs_summary.tsv"
# FAIL_SUMMARY="$OUT_ROOT/failures.tsv"

# # init headers once
# _init_table(){
#   local path="$1" header="$2"
#   if [[ ! -f "$path" ]]; then
#     printf '%s\n' "$header" > "$path"
#   fi
# }
# _init_table "$RUNS_SUMMARY" "ts\tproject\tsample\tfastq_dir\tstatus\texit_code\tout_dir\tweb_summary\tlog_file"
# _init_table "$FAIL_SUMMARY" "ts\tproject\tsample\tfastq_dir\texit_code\tout_dir\tlog_file"

# # atomic append
# safe_append(){
#   local file="$1"; shift
#   local tmp line
#   line="$*"
#   (
#     flock -w 30 9 || { echo "WARN: lock timeout for $file" >&2; exit 0; }
#     printf '%s\n' "$line" >> "$file"
#   ) 9>>"${file}.lock"
# }

# # ===== Runner for one line from manifest =====
# run_one() {
#   local project="$1" sample="$2" fastq_dir="$3"
#   local id="$sample"
#   local outdir="$OUT_ROOT/$project"
#   mkdir -p "$outdir"

#   # per-sample log collects full cellranger output
#   local sample_log="$outdir/${id}.log"
#   printf '==== %s START sample=%s ====\n' "$(ts)" "$id" >> "$sample_log"

#   # quick checks
#   if [[ ! -d "$fastq_dir" ]]; then
#     log "SKIP (fastq_dir missing): $sample -> $fastq_dir"
#     safe_append "$RUNS_SUMMARY" "$(ts)"$'\t'"$project"$'\t'"$sample"$'\t'"$fastq_dir"$'\t'"SKIP_NO_FASTQDIR"$'\t'0$'\t'"$outdir"$'\t'""$'\t'"$sample_log"
#     return 0
#   fi
#   if ! compgen -G "$fastq_dir/*_R1_001.fastq.gz" >/dev/null || \
#      ! compgen -G "$fastq_dir/*_R2_001.fastq.gz" >/dev/null; then
#     log "SKIP (no R1/R2 pattern in dir): $sample -> $fastq_dir"
#     safe_append "$RUNS_SUMMARY" "$(ts)"$'\t'"$project"$'\t'"$sample"$'\t'"$fastq_dir"$'\t'"SKIP_NO_R1R2"$'\t'0$'\t'"$outdir"$'\t'""$'\t'"$sample_log"
#     return 0
#   fi

#   # skip if already done
#   if [[ -s "$outdir/$id/outs/web_summary.html" ]]; then
#     log "OK (exists, skipping): $id"
#     safe_append "$RUNS_SUMMARY" "$(ts)"$'\t'"$project"$'\t'"$sample"$'\t'"$fastq_dir"$'\t'"OK_EXISTS"$'\t'0$'\t'"$outdir/$id/outs"$'\t'"$outdir/$id/outs/web_summary.html"$'\t'"$sample_log"
#     return 0
#   fi

#   log "START: $id  fastqs=$fastq_dir  out=$outdir"

#   # build command
#   cmd=( cellranger count
#         --id="$id"
#         --transcriptome="$TRANSCRIPTOME"
#         --fastqs="$fastq_dir"
#         --sample="$sample"
#         --create-bam="$CREATE_BAM"
#       )
#   [[ -n "$INCLUDE_INTRONS" ]] && cmd+=( --include-introns="$INCLUDE_INTRONS" )
#   [[ -n "$CHEMISTRY"      ]] && cmd+=( --chemistry="$CHEMISTRY" )
#   [[ -n "$EXPECT_CELLS"   ]] && cmd+=( --expect-cells="$EXPECT_CELLS" )
#   [[ -n "$LOCAL_CORES"    ]] && cmd+=( --localcores="$LOCAL_CORES" )
#   [[ -n "$LOCAL_MEM"      ]] && cmd+=( --localmem="$LOCAL_MEM" )

#   # run: send full output to per-sample log (master still gets high-level via log())
#   ( cd "$outdir" && "${cmd[@]}" ) >>"$sample_log" 2>&1
#   local rc=$?

#   local ws="$outdir/$id/outs/web_summary.html"
#   if [[ $rc -eq 0 ]]; then
#     log "DONE: $id → $outdir/$id/outs"
#     safe_append "$RUNS_SUMMARY" "$(ts)"$'\t'"$project"$'\t'"$sample"$'\t'"$fastq_dir"$'\t'"DONE"$'\t'"$rc"$'\t'"$outdir/$id/outs"$'\t'"$ws"$'\t'"$sample_log"
#     printf '==== %s DONE sample=%s (exit=%s) ====\n' "$(ts)" "$id" "$rc" >> "$sample_log"
#   else
#     log "FAIL ($rc): $id"
#     safe_append "$RUNS_SUMMARY" "$(ts)"$'\t'"$project"$'\t'"$sample"$'\t'"$fastq_dir"$'\t'"FAIL"$'\t'"$rc"$'\t'"$outdir/$id/outs"$'\t'"$ws"$'\t'"$sample_log"
#     safe_append "$FAIL_SUMMARY" "$(ts)"$'\t'"$project"$'\t'"$sample"$'\t'"$fastq_dir"$'\t'"$rc"$'\t'"$outdir/$id/outs"$'\t'"$sample_log"
#     printf '==== %s FAIL sample=%s (exit=%s) ====\n' "$(ts)" "$id" "$rc" >> "$sample_log"
#   fi
#   return $rc
# }

# export -f run_one log ts safe_append
# export TRANSCRIPTOME OUT_ROOT CREATE_BAM INCLUDE_INTRONS CHEMISTRY EXPECT_CELLS LOCAL_CORES LOCAL_MEM RUNS_SUMMARY FAIL_SUMMARY

# # ===== Read manifest and run =====
# if (( MAX_JOBS > 1 )); then
#   active=0
#   while IFS=$'\t' read -r project sample fastq_dir; do
#     [[ -z "${project:-}" || -z "${sample:-}" || -z "${fastq_dir:-}" ]] && continue
#     bash -lc 'run_one "$@"' _ "$project" "$sample" "$fastq_dir" &
#     (( active+=1 ))
#     if (( active >= MAX_JOBS )); then
#       wait -n
#       (( active-=1 ))
#     fi
#   done < <(awk -F'\t' 'NR>1{print $1"\t"$2"\t"$3}' "$MANIFEST")
#   wait
# else
#   awk -F'\t' 'NR>1{print $1"\t"$2"\t"$3}' "$MANIFEST" \
#   | while IFS=$'\t' read -r project sample fastq_dir; do
#       [[ -z "${project:-}" || -z "${sample:-}" || -z "${fastq_dir:-}" ]] && continue
#       run_one "$project" "$sample" "$fastq_dir"
#     done
# fi

# log "=== ALL DONE ==="
# log "Master log:   $LOGFILE"
# log "Summary TSV:  $RUNS_SUMMARY"
# log "Failures TSV: $FAIL_SUMMARY"

#!/usr/bin/env bash
set -Euo pipefail
IFS=$'\n\t'

# ===== Required paths (edit to your real paths) =====
MANIFEST="/mnt/12T/chibao/data/official_data/fastq_by_tech/bulk_or_plate_RNA/official_cell/_cr_samples.tsv"  # project_id  sample  fastq_dir
TRANSCRIPTOME="/mnt/12T/chibao/env_tool/cellranger/ref_genome/refdata-gex-GRCh38-2020-A"
OUT_ROOT="/mnt/12T/chibao/data/cellranger_data/cell_ranger_out"   # where per-sample ID folders will be created

# ===== Options (simple) =====
CREATE_BAM="true"          # set to "true" if you also want BAMs (v8+); CRAM is default
INCLUDE_INTRONS=""          # set to "true" for snRNA/pre-mRNA, else leave empty
CHEMISTRY=""                # e.g. "SC3Pv3", "SC3Pv2", "SC5P-PE"; leave empty to auto-detect
EXPECT_CELLS=""             # e.g. "6000"; leave empty to let cellranger estimate
MAX_JOBS=3                 # 1 = sequential; increase (e.g., 2, 3) for light parallel

# (optional) local cores/mem; set to "" to omit
LOCAL_CORES="16"
LOCAL_MEM="26"   # GB

# ===== Tiny helpers =====
ts(){ date '+%Y-%m-%d %H:%M:%S%z'; }
log(){ printf '[%s] %s\n' "$(ts)" "$*" >&2; }

# ===== Sanity checks =====
command -v cellranger >/dev/null || { log "ERROR: cellranger not found in PATH"; exit 1; }
command -v flock       >/dev/null || { log "ERROR: flock not found (install util-linux)"; exit 1; }
[[ -d "$TRANSCRIPTOME" ]] || { log "ERROR: transcriptome not found: $TRANSCRIPTOME"; exit 1; }
[[ -s "$MANIFEST" ]] || { log "ERROR: manifest not found or empty: $MANIFEST"; exit 1; }
mkdir -p "$OUT_ROOT"

# ===== Global master log (high-level) =====
LOGFILE="$OUT_ROOT/cellranger_master.log"
exec > >(tee -a "$LOGFILE") 2>&1
log "=== START run === MANIFEST=$MANIFEST  OUT_ROOT=$OUT_ROOT  TRANSCRIPTOME=$TRANSCRIPTOME"

# ===== Summary tables (parallel-safe with flock) =====
RUNS_SUMMARY="$OUT_ROOT/runs_summary.tsv"
FAIL_SUMMARY="$OUT_ROOT/failures.tsv"

_init_table(){
  local path="$1" header="$2"
  if [[ ! -f "$path" ]]; then
    printf '%s\n' "$header" > "$path"
  fi
}
_init_table "$RUNS_SUMMARY" "ts\tproject\tsample\tfastq_dir\tstatus\texit_code\tout_dir\tweb_summary\tlog_file"
# NOTE: add a reason column
_init_table "$FAIL_SUMMARY" "ts\tproject\tsample\tfastq_dir\texit_code\tout_dir\tlog_file\treason"

# atomic append
safe_append(){
  local file="$1"; shift
  local line="$*"
  (
    flock -w 30 9 || { echo "WARN: lock timeout for $file" >&2; exit 0; }
    printf '%s\n' "$line" >> "$file"
  ) 9>>"${file}.lock"
}

# pull a short failure reason from the sample log
extract_reason(){
  local logfile="$1"
  local reason=""
  if [[ -s "$logfile" ]]; then
    # Prefer Cell Ranger "Log message:" block if present
    reason=$(awk '/Log message:/,0 {print}' "$logfile" | head -n 5 | tr '\t' ' ' | paste -sd ' ' -)
    # Fallback to last error-like line
    if [[ -z "$reason" ]]; then
      reason=$(grep -E "Pipestance failed|ERROR|CRITICAL|chemistr|reference|fastq|No such file" "$logfile" | tail -n 1 || true)
    fi
  fi
  # Keep it single-line and not too long
  reason="$(echo "${reason:-unknown}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' | cut -c1-400)"
  printf '%s' "$reason"
}

run_one() {
  local project="$1" sample="$2" fastq_dir="$3"
  local id="$sample"
  local outdir="$OUT_ROOT/$project"
  mkdir -p "$outdir"

  local sample_log="$outdir/${id}.log"
  printf '==== %s START sample=%s ====\n' "$(ts)" "$id" >> "$sample_log"

  # quick checks
  if [[ ! -d "$fastq_dir" ]]; then
    log "SKIP (fastq_dir missing): $sample -> $fastq_dir"
    safe_append "$RUNS_SUMMARY" "$(ts)"$'\t'"$project"$'\t'"$sample"$'\t'"$fastq_dir"$'\t'"SKIP_NO_FASTQDIR"$'\t'0$'\t'"$outdir"$'\t'""$'\t'"$sample_log"
    return 0
  fi
  if ! compgen -G "$fastq_dir/*_R1_001.fastq.gz" >/dev/null || \
     ! compgen -G "$fastq_dir/*_R2_001.fastq.gz" >/dev/null; then
    log "SKIP (no R1/R2 pattern in dir): $sample -> $fastq_dir"
    safe_append "$RUNS_SUMMARY" "$(ts)"$'\t'"$project"$'\t'"$sample"$'\t'"$fastq_dir"$'\t'"SKIP_NO_R1R2"$'\t'0$'\t'"$outdir"$'\t'""$'\t'"$sample_log"
    return 0
  fi

  # skip if already done
  local ws="$outdir/$id/outs/web_summary.html"
  if [[ -s "$ws" ]]; then
    log "OK (exists, skipping): $id"
    safe_append "$RUNS_SUMMARY" "$(ts)"$'\t'"$project"$'\t'"$sample"$'\t'"$fastq_dir"$'\t'"OK_EXISTS"$'\t'0$'\t'"$outdir/$id/outs"$'\t'"$ws"$'\t'"$sample_log"
    return 0
  fi

  log "START: $id  fastqs=$fastq_dir  out=$outdir"

  # build command
  cmd=( cellranger count
        --id="$id"
        --transcriptome="$TRANSCRIPTOME"
        --fastqs="$fastq_dir"
        --sample="$sample"
        --create-bam="$CREATE_BAM"
      )
  [[ -n "$INCLUDE_INTRONS" ]] && cmd+=( --include-introns="$INCLUDE_INTRONS" )
  [[ -n "$CHEMISTRY"      ]] && cmd+=( --chemistry="$CHEMISTRY" )
  [[ -n "$EXPECT_CELLS"   ]] && cmd+=( --expect-cells="$EXPECT_CELLS" )
  [[ -n "$LOCAL_CORES"    ]] && cmd+=( --localcores="$LOCAL_CORES" )
  [[ -n "$LOCAL_MEM"      ]] && cmd+=( --localmem="$LOCAL_MEM" )

  # run Cell Ranger; capture exit code
  ( cd "$outdir" && "${cmd[@]}" ) >>"$sample_log" 2>&1
  local rc=$?

  if [[ $rc -eq 0 && -s "$ws" ]]; then
    log "DONE: $id → $outdir/$id/outs"
    safe_append "$RUNS_SUMMARY" "$(ts)"$'\t'"$project"$'\t'"$sample"$'\t'"$fastq_dir"$'\t'"DONE"$'\t'"$rc"$'\t'"$outdir/$id/outs"$'\t'"$ws"$'\t'"$sample_log"
    printf '==== %s DONE sample=%s (exit=%s) ====\n' "$(ts)" "$id" "$rc" >> "$sample_log"
  else
    log "FAIL ($rc): $id"
    local reason; reason="$(extract_reason "$sample_log")"
    # write per-sample reason file
    printf '%s\n' "$reason" > "$outdir/$id/FAIL_REASON.txt" || true
    # summaries
    safe_append "$RUNS_SUMMARY" "$(ts)"$'\t'"$project"$'\t'"$sample"$'\t'"$fastq_dir"$'\t'"FAIL"$'\t'"$rc"$'\t'"$outdir/$id/outs"$'\t'"$ws"$'\t'"$sample_log"
    safe_append "$FAIL_SUMMARY" "$(ts)"$'\t'"$project"$'\t'"$sample"$'\t'"$fastq_dir"$'\t'"$rc"$'\t'"$outdir/$id/outs"$'\t'"$sample_log"$'\t'"$reason"
    printf '==== %s FAIL sample=%s (exit=%s) ====\n' "$(ts)" "$id" "$rc" >> "$sample_log"
  fi

  return 0   # IMPORTANT: never propagate failure; keep the batch going
}

export -f run_one log ts safe_append extract_reason
export TRANSCRIPTOME OUT_ROOT CREATE_BAM INCLUDE_INTRONS CHEMISTRY EXPECT_CELLS LOCAL_CORES LOCAL_MEM RUNS_SUMMARY FAIL_SUMMARY

# ===== Read manifest and run =====
if (( MAX_JOBS > 1 )); then
  active=0
  while IFS=$'\t' read -r project sample fastq_dir; do
    [[ -z "${project:-}" || -z "${sample:-}" || -z "${fastq_dir:-}" ]] && continue
    # background; each run captures its own rc/summary
    bash -lc 'run_one "$@"' _ "$project" "$sample" "$fastq_dir" &
    (( active+=1 ))
    if (( active >= MAX_JOBS )); then
      wait -n
      (( active-=1 ))
    fi
  done < <(awk -F'\t' 'NR>1{print $1"\t"$2"\t"$3}' "$MANIFEST")
  wait
else
  # Ensure a failure doesn't abort the while-loop (|| true)
  awk -F'\t' 'NR>1{print $1"\t"$2"\t"$3}' "$MANIFEST" \
  | while IFS=$'\t' read -r project sample fastq_dir; do
      [[ -z "${project:-}" || -z "${sample:-}" || -z "${fastq_dir:-}" ]] && continue
      run_one "$project" "$sample" "$fastq_dir" || true
    done
fi

log "=== ALL DONE ==="
log "Master log:   $LOGFILE"
log "Summary TSV:  $RUNS_SUMMARY"
log "Failures TSV: $FAIL_SUMMARY"

