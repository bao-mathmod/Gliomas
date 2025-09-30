#!/usr/bin/env bash
set -Euo pipefail
IFS=$'\n\t'

# =========================
# CONFIG (edit these)
# =========================
SHEET="/mnt/10T/chibao/gliomas/data/PRJNA961045/PRJNA961045_v2.csv"   # CSV: [sample?],fastq1,fastq2,index1,index2
REFERENCE="/mnt/12T/chibao/env_tool/cellranger_atac/refgenome_atac/refdata-cellranger-arc-GRCh38-2024-A"
OUT_ROOT="/mnt/12T/chibao/data/cellranger_data"

# Optional runtime tuning
LOCAL_CORES="16"         # "" to omit
LOCAL_MEM="32"           # GB, "" to omit
# Parallelism
MAX_JOBS=2    # 1 = sequential; increase to run samples in parallel

# Logging model (same shape as your GEX script)
LOGFILE="$OUT_ROOT/cellranger_atac_master.log"
RUNS_SUMMARY="$OUT_ROOT/atac_runs_summary.tsv"
FAIL_SUMMARY="$OUT_ROOT/atac_failures.tsv"

# =========================
# HELPERS
# =========================
ts(){ date '+%Y-%m-%d %H:%M:%S%z'; }
log(){ printf '[%s] %s\n' "$(ts)" "$*" >&2; }

_init_table(){
  local path="$1" header="$2"
  [[ -f "$path" ]] || printf '%s\n' "$header" > "$path"
}

safe_append(){
  local file="$1"; shift
  local line="$*"
  (
    flock -w 30 9 || { echo "WARN: lock timeout for $file" >&2; exit 0; }
    printf '%s\n' "$line" >> "$file"
  ) 9>>"${file}.lock"
}

extract_reason(){
  local logfile="$1" reason=""
  if [[ -s "$logfile" ]]; then
    reason=$(awk '/Log message:/,0 {print}' "$logfile" | head -n 5 | tr '\t' ' ' | paste -sd ' ' -)
    [[ -z "$reason" ]] && reason=$(grep -E "Pipestance failed|ERROR|CRITICAL|reference|fastq|No such file|index|barcode|I2" "$logfile" | tail -n 1 || true)
  fi
  reason="$(echo "${reason:-unknown}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//' | cut -c1-400)"
  printf '%s' "$reason"
}

# Return dirname of a file path
dir_of(){ python3 - "$1" <<'PY'
import os,sys; print(os.path.dirname(sys.argv[1]) or ".")
PY
}

# Extract the FASTQ "sample prefix" (text before "_S" in Illumina-style file name).
# Falls back to basename without first underscore chunk if pattern missing.
prefix_of_fastq(){
  local base ; base="$(basename "$1")"
  if [[ "$base" =~ ^([^_]+)_S[0-9]+_L[0-9]{3}_[IR][12]_001\.fastq(\.gz)?$ ]]; then
    printf '%s' "${BASH_REMATCH[1]}"
  else
    # Try common SRR-style: SRRxxxxx_1.fastq.gz -> SRRxxxxx
    printf '%s' "${base%%_*}"
  fi
}

# =========================
# SANITY CHECKS
# =========================
command -v cellranger-atac >/dev/null || { echo "ERROR: cellranger-atac not found in PATH"; exit 1; }
command -v flock >/dev/null || { echo "ERROR: flock not found (install util-linux)"; exit 1; }
mkdir -p "$OUT_ROOT"

# Reference guard: warn if looks like ARC
if [[ "$REFERENCE" =~ cellranger-arc ]]; then
  log "WARNING: Reference path looks like ARC (multiome). cellranger-atac needs an ATAC bundle (e.g., refdata-cellranger-atac-...)."
fi
[[ -d "$REFERENCE" ]] || { echo "ERROR: reference not found: $REFERENCE"; exit 1; }
[[ -s "$SHEET"    ]] || { echo "ERROR: sheet not found or empty: $SHEET"; exit 1; }

# Init logs
exec > >(tee -a "$LOGFILE") 2>&1
log "=== START ATAC === SHEET=$SHEET  OUT_ROOT=$OUT_ROOT  REFERENCE=$REFERENCE"
_init_table "$RUNS_SUMMARY" "ts\tsample_logic\tprefixes\tfastq_dirs\tstatus\texit_code\touts_dir\tweb_summary\tlog_file"
_init_table "$FAIL_SUMMARY" "ts\tsample_logic\tprefixes\tfastq_dirs\texit_code\touts_dir\tlog_file\treason"

# =========================
# PARSE CSV & GROUP
# =========================
# We accept either:
#   sample,fastq1,fastq2,index1,index2
# or (if sample omitted):
#   fastq1,fastq2,index1,index2      (we will derive 'sample_logic' from directory ancestry)
#
# We will:
#  - Verify that index2 exists (I5, required). index1 optional.
#  - Group rows by 'sample_logic' (prefer CSV 'sample'; else derive from the directory 3-level basename).
#  - For each group, accumulate:
#       - FASTQ directories (unique) for --fastqs
#       - Sample prefixes from R1 basenames (unique) for --sample
#

python3 - "$SHEET" <<'PY' >"$OUT_ROOT/.atac_plan.tsv"
import csv, os, sys, re
sheet = sys.argv[1]
with open(sheet, newline='') as f:
    reader = csv.DictReader(f)
    headers = [h.strip().lower() for h in reader.fieldnames]
    need = {'fastq1','fastq2','index2'}
    if not need.issubset(set(headers)):
        sys.exit("CSV must have at least: fastq1, fastq2, index2 (I7 index1 optional, sample optional)")
    has_sample = 'sample' in headers

    # very lightweight heuristic to derive a logical sample if not provided
    def derive_logic(p):
        # SAMN if present in path; else the SRX/parent; else the immediate parent dir
        parts = p.split('/')
        for key in ('SAMN','SAM','biosample','SAMPLE'):
            for s in parts:
                if key in s.upper():
                    return s
        # fallback: take 3-level parent to avoid SRR folder
        return '/'.join(parts[-4:-1]) if len(parts) >= 4 else os.path.basename(os.path.dirname(p)) or "SAMPLE_UNKNOWN"

    plan = {}  # logic_sample -> {dirs:set(), prefixes:set()}
    def prefix_of_fastq(p):
        b = os.path.basename(p)
        m = re.match(r'^([^_]+)_S\d+_L\d{3}_[IR]\d_001\.fastq(?:\.gz)?$', b)
        if m: return m.group(1)
        return b.split('_')[0]

    for row in reader:
        fastq1 = row.get('fastq1','').strip()
        fastq2 = row.get('fastq2','').strip()
        i1     = row.get('index1','').strip()
        i2     = row.get('index2','').strip()
        sample = row.get('sample','').strip() if has_sample else ''
        # sanity
        if not (fastq1 and fastq2 and i2):
            print(f"# SKIP row (missing fastq1/fastq2/index2): {row}", file=sys.stderr)
            continue
        # derive logical sample
        logic = sample if sample else derive_logic(fastq1)
        d1 = os.path.dirname(fastq1); d2 = os.path.dirname(fastq2)
        di2 = os.path.dirname(i2);   di1 = os.path.dirname(i1) if i1 else ""
        # prefixes from fastq1 filename
        pref = prefix_of_fastq(fastq1)
        # collect
        entry = plan.setdefault(logic, {"dirs": set(), "prefixes": set(), "missing": []})
        for d in (d1,d2,di2,di1):
            if d: entry["dirs"].add(d)
        entry["prefixes"].add(pref)
        # record missing optional I1
        if not i1:
            entry["missing"].append("no_I1")

    # emit plan table: logic_sample \t fastqs_csv \t prefixes_csv \t flags
    print("sample_logic\tfastq_dirs\tprefixes\tflags")
    for logic, v in plan.items():
        flags = ",".join(sorted(set(v["missing"]))) if v["missing"] else ""
        print(f"{logic}\t{','.join(sorted(v['dirs']))}\t{','.join(sorted(v['prefixes']))}\t{flags}")
PY

if [[ ! -s "$OUT_ROOT/.atac_plan.tsv" ]]; then
  log "ERROR: no runnable rows parsed from CSV. Check headers/paths."
  exit 1
fi

# =========================
# RUN EACH LOGICAL SAMPLE
# =========================

run_one_sample() {
  local logic_sample="$1" fastq_dirs="$2" prefixes="$3" flags="$4"

  # Per-sample out + log
  local outdir="$OUT_ROOT/$logic_sample"
  mkdir -p "$outdir"
  local sample_log="$outdir/${logic_sample}.log"
  printf '==== %s START sample=%s ====\n' "$(ts)" "$logic_sample" >> "$sample_log"

  # Soft warnings
  [[ "$flags" == *"no_I1"* ]] && log "WARN: $logic_sample has no I1 (index1). This is OK (I1 optional), but ensure your run truly lacks i7."

  # Build command
  local cmd=( cellranger-atac count
              --id="$logic_sample"
              --reference="$REFERENCE"
              --fastqs="$fastq_dirs"
              --sample="$prefixes" )
  [[ -n "$LOCAL_CORES" ]] && cmd+=( --localcores="$LOCAL_CORES" )
  [[ -n "$LOCAL_MEM"   ]] && cmd+=( --localmem="$LOCAL_MEM" )

  log "START: $logic_sample  fastqs=[$fastq_dirs]  prefixes=[$prefixes]  out=$outdir"
  ( cd "$outdir" && "${cmd[@]}" ) >>"$sample_log" 2>&1
  local rc=$?

  local ws="$outdir/$logic_sample/outs/web_summary.html"
  if [[ $rc -eq 0 && -s "$ws" ]]; then
    log "DONE: $logic_sample → $outdir/$logic_sample/outs"
    safe_append "$RUNS_SUMMARY" \
      "$(ts)"$'\t'"$logic_sample"$'\t'"$prefixes"$'\t'"$fastq_dirs"$'\t'"DONE"$'\t'"$rc"$'\t'"$outdir/$logic_sample/outs"$'\t'"$ws"$'\t'"$sample_log"
    printf '==== %s DONE sample=%s (exit=%s) ====\n' "$(ts)" "$logic_sample" "$rc" >> "$sample_log"
  else
    log "FAIL ($rc): $logic_sample"
    local reason; reason="$(extract_reason "$sample_log")"
    printf '%s\n' "$reason" > "$outdir/$logic_sample/FAIL_REASON.txt" || true
    safe_append "$RUNS_SUMMARY" \
      "$(ts)"$'\t'"$logic_sample"$'\t'"$prefixes"$'\t'"$fastq_dirs"$'\t'"FAIL"$'\t'"$rc"$'\t'"$outdir/$logic_sample/outs"$'\t'"$ws"$'\t'"$sample_log"
    safe_append "$FAIL_SUMMARY" \
      "$(ts)"$'\t'"$logic_sample"$'\t'"$prefixes"$'\t'"$fastq_dirs"$'\t'"$rc"$'\t'"$outdir/$logic_sample/outs"$'\t'"$sample_log"$'\t'"$reason"
    printf '==== %s FAIL sample=%s (exit=%s) ====\n' "$(ts)" "$logic_sample" "$rc" >> "$sample_log"
  fi
}

# Wait-for-slot helper (works with/without 'wait -n')
wait_for_slot() {
  local max="$1"
  # If wait -n exists, use it; else count jobs and wait for any to finish
  if help wait 2>/dev/null | grep -q -- 'wait [-fn]'; then
    while (( $(jobs -rp | wc -l) >= max )); do
      wait -n || wait || true
    done
  else
    while (( $(jobs -rp | wc -l) >= max )); do
      local pid; pid=$(jobs -rp | head -n1)
      [[ -n "$pid" ]] && wait "$pid" || sleep 1
    done
  fi
}

export -f run_one_sample ts log safe_append extract_reason
export OUT_ROOT REFERENCE RUNS_SUMMARY FAIL_SUMMARY LOCAL_CORES LOCAL_MEM

# Header skip + parallel dispatch
if (( MAX_JOBS > 1 )); then
  while IFS=$'\t' read -r logic_sample fastq_dirs prefixes flags; do
    [[ "$logic_sample" == "sample_logic" ]] && continue
    wait_for_slot "$MAX_JOBS"
    # Launch in background
    bash -lc 'run_one_sample "$@"' _ "$logic_sample" "$fastq_dirs" "$prefixes" "$flags" &
  done < "$OUT_ROOT/.atac_plan.tsv"
  # Wait for remaining jobs
  wait
else
  # Sequential
  while IFS=$'\t' read -r logic_sample fastq_dirs prefixes flags; do
    [[ "$logic_sample" == "sample_logic" ]] && continue
    run_one_sample "$logic_sample" "$fastq_dirs" "$prefixes" "$flags"
  done < "$OUT_ROOT/.atac_plan.tsv"
fi

log "=== ALL DONE (ATAC) ==="
log "Master log:   $LOGFILE"
log "Summary TSV:  $RUNS_SUMMARY"
log "Failures TSV: $FAIL_SUMMARY"
