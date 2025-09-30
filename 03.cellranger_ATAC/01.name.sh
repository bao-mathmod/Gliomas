#!/usr/bin/env bash
# ============================================================
# 05_atac_prep_detect_and_link.sh
#
# Detect per-run read roles (R1/R2/R3[/I1]) for (putative) scATAC datasets,
# summarize, and (optionally) create Cell Ranger ATAC–style symlinks.
#
# WHY: Public datasets mix naming styles. For ATAC you *must* map the
# short barcode read to R2 and long genomic reads to R1/R3. This script
# detects by measuring read lengths and applies robust rules.
#
# USAGE:
#   bash 05_atac_prep_detect_and_link.sh \
#     [--fastq-root /mnt/12T/chibao/data/official_data/fastq] \
#     [--enriched-root /mnt/12T/chibao/data/official_data/metadata_enriched] \
#     [--expected-index /mnt/12T/chibao/data/official_data/fastq/_expected_files.tsv] \
#     [--out-root /mnt/12T/chibao/data/official_data/atac_prep] \
#     [--min-short 10] [--max-short 32] \
#     [--do-link] [--dry-run]
#
# Inputs:
#   - FASTQs in a PRJ/Study/Sample/Experiment/Run tree (your downloader layout)
#   - (Optional but recommended) enriched metadata to check library_strategy
#   - The expected index (to discover run directories reliably)
#
# Outputs (under --out-root):
#   _prep.log                         # human log (append-safe)
#   runs_detect.tsv                   # per file: measured length, inferred role
#   runs_plan.tsv                     # per run: chosen mapping, decision, reason
#   runs_actions.tsv                  # symlink actions (if --do-link)
#   problems.tsv                      # incomplete/ambiguous runs with guidance
#
# Symlink layout (if --do-link):
#   <out-root>/fastq_cellranger_atac/<project>/<sample>/<run>/
#       <sample>_S1_L001_R1_001.fastq.gz
#       <sample>_S1_L001_R2_001.fastq.gz  (short barcode)
#       <sample>_S1_L001_R3_001.fastq.gz
#       [<sample>_S1_L001_I1_001.fastq.gz] (if present)
#
# Re-runnable: appends TSVs, updates links idempotently, never edits source.

#In 10x scATAC, the expected reads are:

# R1 = long genomic (≈50 bp)

# R2 = short barcode/index read (often ~16–28 cycles, platform-dependent)

# R3 = long genomic (≈50 bp)

# I1 = sample index (present only if multiplexed)
# ============================================================

#!/usr/bin/env bash
set -Euo pipefail
IFS=$'\n\t'

# ---------- Defaults ----------
# NOTE: FASTQ_ROOT now points to the *scATAC view* produced by 04_organize_fastq_by_tech.sh
FASTQ_ROOT="/mnt/12T/chibao/data/official_data/fastq_by_tech/scATAC"
ENRICHED_ROOT="/mnt/12T/chibao/data/official_data/metadata_enriched"
# We consume the organizer’s symlink manifest as the “expected index”
EXPECTED_INDEX="/mnt/12T/chibao/data/official_data/fastq_by_tech/_symlink_manifest.tsv"
OUT_ROOT="/mnt/12T/chibao/data/official_data/atac_prep"

# length heuristics (inclusive)
MIN_SHORT=10
MAX_SHORT=32

DO_LINK="false"
DRY_RUN="false"

# ---------- CLI ----------
while (( $# )); do
  case "$1" in
    --fastq-root)     FASTQ_ROOT="${2:-}"; shift 2;;
    --enriched-root)  ENRICHED_ROOT="${2:-}"; shift 2;;
    --expected-index) EXPECTED_INDEX="${2:-}"; shift 2;;
    --out-root)       OUT_ROOT="${2:-}"; shift 2;;
    --min-short)      MIN_SHORT="${2:-}"; shift 2;;
    --max-short)      MAX_SHORT="${2:-}"; shift 2;;
    --do-link)        DO_LINK="true"; shift 1;;
    --dry-run)        DRY_RUN="true"; shift 1;;
    *) echo "WARN: unknown arg $1" >&2; shift 1;;
  esac
done

# ---------- Paths & logs ----------
mkdir -p "$OUT_ROOT"
LOG="$OUT_ROOT/_prep.log"
DETECT="$OUT_ROOT/runs_detect.tsv"
PLAN="$OUT_ROOT/runs_plan.tsv"
ACTIONS="$OUT_ROOT/runs_actions.tsv"
PROBLEMS="$OUT_ROOT/problems.tsv"

ts(){ date '+%Y-%m-%d %H:%M:%S%z'; }
log(){ printf '[%s] %s\n' "$(ts)" "$*" | tee -a "$LOG" >&2; }

[[ -f "$DETECT"  ]] || echo -e "ts\tproject\tsample\texperiment\trun\tfile\tpath\tfirst_seq_len\trole_hint" > "$DETECT"
[[ -f "$PLAN"    ]] || echo -e "ts\tproject\tsample\trun\tmapping\tdecision\treason" > "$PLAN"
[[ -f "$ACTIONS" ]] || echo -e "ts\tproject\tsample\trun\tsrc\tdst\taction\tnote" > "$ACTIONS"
[[ -f "$PROBLEMS" ]] || echo -e "ts\tproject\tsample\trun\tproblem\tdetail" > "$PROBLEMS"

# ---------- Sanity ----------
[[ -s "$EXPECTED_INDEX" ]] || { log "ERROR: missing expected index: $EXPECTED_INDEX"; exit 1; }
command -v zcat >/dev/null || { log "ERROR: zcat not found"; exit 1; }

log "=== START ATAC PREP ==="
log "FASTQ_ROOT=$FASTQ_ROOT"
log "EXPECTED_INDEX=$EXPECTED_INDEX"
log "ENRICHED_ROOT=$ENRICHED_ROOT"
log "OUT_ROOT=$OUT_ROOT"
log "MIN_SHORT=$MIN_SHORT  MAX_SHORT=$MAX_SHORT  DO_LINK=$DO_LINK  DRY_RUN=$DRY_RUN"

# ---------- Helpers ----------
trim(){ awk '{$1=$1;print}' <<<"$*"; }
project_meta_extra(){
  local proj="$1"
  local t="$ENRICHED_ROOT/$proj/${proj}.extra.tsv"
  [[ -s "$t" && $(wc -l <"$t") -gt 1 ]] && echo "$t" || echo ""
}
lib_strategy_of(){
  local proj="$1" run="$2" meta; meta="$(project_meta_extra "$proj")"
  [[ -z "$meta" ]] && { echo "-"; return 0; }
  awk -F'\t' -v r="$run" '
    NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
    $h["run_accession"]==r {print ($h["library_strategy"]==""?"-":$h["library_strategy"]); exit}
    END{if(NR==1) print "-"}
  ' "$meta"
}
first_len(){
  # Echo exactly ONE integer; 0 on any error/corruption/odd file.
  local f="$1" n
  # Don’t use “|| echo 0” here; capture whatever prints, then sanitize.
  n=$( zcat -- "$f" 2>/dev/null | awk 'NR==2{print length; exit}' 2>/dev/null || true )
  # Keep only a pure integer; anything else becomes 0
  if [[ "$n" =~ ^[0-9]+$ ]]; then
    echo "$n"
  else
    echo 0
  fi
}

role_hint_by_name(){
  local bn="$1"
  if   [[ "$bn" =~ (_I1|\.I1) ]]; then echo "I1"
  elif [[ "$bn" =~ (_R1|\.R1|_1\.f) ]]; then echo "R1?"
  elif [[ "$bn" =~ (_R2|\.R2|_2\.f) ]]; then echo "R2?/R3?"
  else echo "-"
  fi
}
dst_dir_for(){ echo "$OUT_ROOT/fastq_cellranger_atac/$1/$2/$3"; }
mk_link(){
  local proj="$1" sample="$2" run="$3" src="$4" dst="$5"
  local action note="-"
  if [[ ! -e "$src" ]]; then
    action="skip_missing_src"
  else
    mkdir -p "$(dirname "$dst")"
    if [[ -L "$dst" ]]; then
      local cur; cur="$(readlink "$dst" || true)"
      if [[ "$cur" == "$src" ]]; then action="exists_ok"
      else
        if [[ "$DRY_RUN" == "true" ]]; then action="would_update_link"; note="from:$cur"
        else rm -f "$dst"; ln -s "$src" "$dst"; action="updated_link"; note="from:$cur"
        fi
      fi
    elif [[ -e "$dst" ]]; then
      action="conflict_exists"; note="not_symlink"
    else
      if [[ "$DRY_RUN" == "true" ]]; then action="would_create"
      else ln -s "$src" "$dst"; action="created"
      fi
    fi
  fi
  echo -e "$(ts)\t$proj\t$sample\t$run\t$src\t$dst\t$action\t$note" >> "$ACTIONS"
}

# ---------- Build run list from _symlink_manifest.tsv ----------
# Manifest columns: ts(1) project(2) tech_bucket(3) src_path(4) dst_path(5) action(6) note(7)
mapfile -t RUN_KEYS < <(
  awk -F'\t' -v root="$FASTQ_ROOT" '
    NR==1 { next }
    $3!="scATAC" { next }
    $6 !~ /(created|exists_ok|updated_link|would_create|would_update_link)/ { next }

    {
      dst=$5
      # must be inside FASTQ_ROOT
      if (index(dst, root"/")!=1) next

      # .../scATAC/<PRJNA>/<SAMN>/<SRX>/<SRR>/<file>
      n=split(dst, a, "/")
      if (n < 7) next  # safety

      file=a[n]
      run=a[n-1]
      srx=a[n-2]        # <- renamed from "exp" to avoid awk exp()
      samp=a[n-3]
      proj=a[n-4]

      if (proj=="" || samp=="" || run=="") next

      print proj "\t" samp "\t" srx "\t" run
    }
  ' "$EXPECTED_INDEX" | sort -u
)

# ---------- Per run detection ----------
for key in "${RUN_KEYS[@]}"; do
IFS=$'\t' read -r proj sample exp run <<< "$key"   # "exp" here is just a Bash var, fine


  # New layout from organizer:
  # <FASTQ_ROOT>/<project>/<sample>/<run>
rundir="$FASTQ_ROOT/$proj/$sample/$exp/$run"

  [[ -d "$rundir" ]] || {
    echo -e "$(ts)\t$proj\t$sample\t$run\trun_dir_not_on_disk\t$rundir" >> "$PROBLEMS"
    continue
  }

  shopt -s nullglob
  files=( "$rundir"/*.fastq.gz "$rundir"/*.fq.gz )
  shopt -u nullglob
  if (( ${#files[@]} == 0 )); then
    echo -e "$(ts)\t$proj\t$sample\t$run\tno_fastqs\t$rundir" >> "$PROBLEMS"
    continue
  fi

  declare -A LEN=() HINT=()
  for f in "${files[@]}"; do
    bn="$(basename "$f")"
    l="$(first_len "$f")"; l="${l:-0}"
    h="$(role_hint_by_name "$bn")"
    LEN["$bn"]="$l"; HINT["$bn"]="$h"
    echo -e "$(ts)\t$proj\t$sample\t$exp\t$run\t$bn\t$f\t$l\t$h" >> "$DETECT"
  done

short=() long=()
for bn in "${!LEN[@]}"; do
  l="${LEN[$bn]}"
  if [[ "$l" =~ ^[0-9]+$ ]]; then
    if (( l >= MIN_SHORT && l <= MAX_SHORT )); then
      short+=( "$bn" )
    elif (( l > 0 )); then
      long+=( "$bn" )
    fi
  else
    # Optional: log a hint once per file
    log "WARN: non-numeric length for $bn (got: '$l'); treating as 0"
  fi
done

  libstr="$(lib_strategy_of "$proj" "$run")"
  lowlib="$(echo "$libstr" | tr '[:upper:]' '[:lower:]')"

  mapping="" reason="" decision=""

  if [[ "$lowlib" == *atac* ]]; then
    if (( ${#short[@]} == 1 && ${#long[@]} >= 2 )); then
      IFS=$'\n' read -r -d '' -a sorted_long < <(printf '%s\n' "${long[@]}" | sort && printf '\0')
      r2="${short[0]}"
      r1=""
      r3=""
      for bn in "${sorted_long[@]}"; do
        if [[ "$bn" =~ _1\.f(ast)?q\.gz$ || "$bn" =~ _R1 ]]; then r1="$bn"; break; fi
      done
      for bn in "${sorted_long[@]}"; do
        if [[ "$bn" =~ _2\.f(ast)?q\.gz$ || "$bn" =~ _R3 ]]; then r3="$bn"; fi
      done
      if [[ -z "$r1" || -z "$r3" ]]; then
        r1="${sorted_long[0]}"
        r3="${sorted_long[1]:-${sorted_long[0]}}"
      fi

      i1=""
      for bn in "${!LEN[@]}"; do
        if [[ "$bn" =~ (_I1|\.I1) ]]; then i1="$bn"; break; fi
      done

      mapping="R1=$r1;R2=$r2;R3=$r3;I1=${i1:-none}"
      decision="scATAC"
      reason="1_short_${#short[@]}_[${MIN_SHORT}-${MAX_SHORT}], ${#long[@]}_long; lib=${libstr}"

    elif (( ${#short[@]} == 0 && ${#long[@]} == 2 )); then
      mapping="long1=${long[0]};long2=${long[1]}"
      decision="bulk_like_ATAC"
      reason="no_short_barcode_read_found; lib=${libstr}"
    else
      decision="ambiguous_ATAC"
      reason="short=${#short[@]} long=${#long[@]}; lib=${libstr}"
    fi
  else
    decision="not_ATAC_in_metadata"
    reason="library_strategy=${libstr}"
  fi

  echo -e "$(ts)\t$proj\t$sample\t$run\t$mapping\t$decision\t$reason" >> "$PLAN"

  if [[ "$DO_LINK" == "true" && "$decision" == "scATAC" ]]; then
    dstdir="$(dst_dir_for "$proj" "$sample" "$run")"
    # parse mapping → filenames
    src_r1="$rundir/${mapping#*R1=}"; src_r1="${src_r1%%;*}"
    src_r2="$rundir/${mapping#*R2=}"; src_r2="${src_r2%%;*}"
    src_r3="$rundir/${mapping#*R3=}"; src_r3="${src_r3%%;*}"

    dst_r1="$dstdir/${sample}_S1_L001_R1_001.fastq.gz"
    dst_r2="$dstdir/${sample}_S1_L001_R2_001.fastq.gz"
    dst_r3="$dstdir/${sample}_S1_L001_R3_001.fastq.gz"

    mk_link "$proj" "$sample" "$run" "$src_r1" "$dst_r1"
    mk_link "$proj" "$sample" "$run" "$src_r2" "$dst_r2"
    mk_link "$proj" "$sample" "$run" "$src_r3" "$dst_r3"

    if [[ "$mapping" == *"I1="* && "$mapping" != *"I1=none"* ]]; then
      i1bn="${mapping##*I1=}"; i1bn="${i1bn%%;*}"
      src_i1="$rundir/$i1bn"
      dst_i1="$dstdir/${sample}_S1_L001_I1_001.fastq.gz"
      mk_link "$proj" "$sample" "$run" "$src_i1" "$dst_i1"
    fi
  fi

  if [[ "$decision" != "scATAC" ]]; then
    echo -e "$(ts)\t$proj\t$sample\t$run\t$decision\t$reason" >> "$PROBLEMS"
  fi
done

log "=== DONE ATAC PREP ==="
log "Detect TSV:   $DETECT"
log "Plan TSV:     $PLAN"
log "Actions TSV:  $ACTIONS"
log "Problems TSV: $PROBLEMS"

