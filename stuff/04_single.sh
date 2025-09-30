#!/usr/bin/env bash
set -Euo pipefail
IFS=$'\n\t'

# ================================================================
# Post-download audit & (optional) recovery for ENA FASTQs
# Compatible with your downloader’s OUTPUT_DIR layout and globals
# ================================================================

# -------- Config (edit as needed) --------------------------------
# OUTPUT_DIR="/mnt/12T/chibao/data/official_data/fastq"     # same as in downloader
# BASE_DIR_SCAN="$OUTPUT_DIR"                               # scan inside OUTPUT_DIR
# GLOBAL_SUMMARY="$OUTPUT_DIR/_download_summary.tsv"
# EXPECTED_INDEX="$OUTPUT_DIR/_expected_files.tsv"

# Where the FASTQs live (input)
BASE_DIR_SCAN="/mnt/12T/chibao/data/official_data/fastq_by_tech/bulk_or_plate_RNA"

# Where to write logs/summaries (output)
OUTPUT_DIR="/mnt/12T/chibao/data/stuff_data/single"

# These indexes still point to FASTQ tree unless you also moved them
GLOBAL_SUMMARY="$BASE_DIR_SCAN/_download_summary.tsv"
EXPECTED_INDEX="$BASE_DIR_SCAN/_expected_files.tsv"

# Recovery behavior
DO_RECOVER="true"           # if true, attempt SRA prefetch + fasterq-dump when needed
INCLUDE_TECH="true"         # include technical reads (10x / barcodes)
THREADS=8                   # fasterq-dump threads
PREFETCH_ROOT="${OUTPUT_DIR}/_sra_cache"   # where .sra files go
RECOVERY_TMP="${OUTPUT_DIR}/_recovery_tmp" # staging for split files

# Heuristic: treat single-file Illumina "PAIRED" with /3 or no /[12] as 10x-style
ASSUME_10X_IF_SINGLE_AND_PAIRED="true"

# Tools (fail-soft check)
NEED_TOOLS=(curl column awk grep sed wc prefetch fasterq-dump pigz zcat)
# prefetch/fasterq-dump from SRA Toolkit; pigz optional (gzip otherwise)

# -------- Logging -------------------------------------------------
ts(){ date '+%Y-%m-%d %H:%M:%S%z'; }
log(){ printf '[%s] %s\n' "$(ts)" "$*" >&2; }
warn(){ log "WARN: $*"; }
err(){ log "ERROR: $*"; }

# -------- Init outputs (GLOBAL, append-safe) ---------------------
mkdir -p "$RECOVERY_TMP" "$PREFETCH_ROOT"

GLOBAL_AUDIT="$OUTPUT_DIR/_postcheck_audit.tsv"           # one line per run observed on disk
GLOBAL_ENA_CHECK="$OUTPUT_DIR/_postcheck_ena.tsv"         # ENA API echo per checked run
GLOBAL_INTERLEAVE="$OUTPUT_DIR/_postcheck_interleave.tsv" # interleaving probe results
GLOBAL_SRA_LAYOUT="$OUTPUT_DIR/_postcheck_sra_layout.tsv" # reads/spot + technical counts
GLOBAL_ACTIONS="$OUTPUT_DIR/_postcheck_actions.tsv"       # attempted recoveries, moves, renames
GLOBAL_ERRORS="$OUTPUT_DIR/_postcheck_errors.tsv"         # all errors (non-fatal)
GLOBAL_SUMMARY2="$OUTPUT_DIR/_postcheck_summary_by_project.tsv" # project roll-up

init_file(){
  local f="$1" header="$2"
  [[ -f "$f" ]] || { echo -e "$header" > "$f"; }
}
init_file "$GLOBAL_AUDIT"      "project_id\tstudy\tsample\texperiment\trun\tfastq_count\tfiles_list"
init_file "$GLOBAL_ENA_CHECK"  "run\tlibrary_layout\tlibrary_strategy\tinstrument_platform\tinstrument_model\tread_count\tbase_count\tfastq_ftp\tsubmitted_ftp\tsubmitted_format"
init_file "$GLOBAL_INTERLEAVE" "study\tsample\texperiment\trun\tfile\tprobe\tresult\tmsg"
init_file "$GLOBAL_SRA_LAYOUT" "run\tspots\treads_read\treads_written\ttechnical_reads"
init_file "$GLOBAL_ACTIONS"    "ts\tproject_id\tstudy\tsample\texperiment\trun\taction\tstatus\tmessage"
init_file "$GLOBAL_ERRORS"     "ts\tproject_id\tstudy\tsample\texperiment\trun\twhere\tmessage"
init_file "$GLOBAL_SUMMARY2"   "project_id\truns\truns_ok2\truns_only1\truns_zero\trecovered_runs"

# -------- Helpers -------------------------------------------------
safe_join_semicolons(){ local IFS=';'; echo "$*"; }

have(){ command -v "$1" >/dev/null 2>&1; }

need_tools(){
  local missing=0
  for t in "${NEED_TOOLS[@]}"; do
    if ! have "$t"; then
      warn "Missing tool: $t"
      missing=1
    fi
  done
  if (( missing )); then
    warn "Some tools missing. Script will skip related steps and just log."
  fi
}
need_tools

# ENA API for a RUN (SRR/ERR/DRR)
ena_query(){
  local run="$1"
  curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=run_accession,library_layout,library_strategy,instrument_platform,instrument_model,read_count,base_count,fastq_ftp,submitted_ftp,submitted_format&download=false"
}

# Basic interleave probe: look for /1 or /2 alternation or newer Illumina tags
probe_interleave(){
  local fq="$1"
  # take first 400 records (~1600 lines) to keep it fast
  local headz; headz="$(zcat "$fq" 2>/dev/null | sed -n '1,1600p' || true)"
  local c1; c1="$(grep -cE '/1(\s|$)' <<<"$headz" || true)"
  local c2; c2="$(grep -cE '/2(\s|$)' <<<"$headz" || true)"
  if (( c1>0 && c2>0 )); then
    echo "yes"
    return 0
  fi
  # try Illumina CASAVA 1.8+ space-delimited read number (e.g., " 1:N:" / " 2:N:")
  local c1b; c1b="$(grep -cE '^[^[:space:]]+ [12]:[YN]:' <<<"$headz" || true)"
  if (( c1b>0 )); then
    echo "maybe"  # cannot confirm alternation, but headers carry read numbers
    return 0
  fi
  echo "no"
}

record_error(){
  local project="$1" study="$2" sample="$3" exp="$4" run="$5" where="$6" msg="$7"
  echo -e "$(ts)\t${project}\t${study}\t${sample}\t${exp}\t${run}\t${where}\t${msg}" >> "$GLOBAL_ERRORS"
}

record_action(){
  local project="$1" study="$2" sample="$3" exp="$4" run="$5" action="$6" status="$7" msg="$8"
  echo -e "$(ts)\t${project}\t${study}\t${sample}\t${exp}\t${run}\t${action}\t${status}\t${msg}" >> "$GLOBAL_ACTIONS"
}

summarize_by_project(){
  # Build per-project roll-up using GLOBAL_AUDIT and GLOBAL_ACTIONS
  # runs: total runs seen; ok2: fastq_count>=2; only1: ==1; zero:==0; recovered: count of actions=RECOVER_OK
  awk -F'\t' 'NR>1{key=$1; runs[key]++; 
                   if($6+0>=2) ok2[key]++;
                   else if($6+0==1) one[key]++;
                   else zero[key]++}
               END{for(k in runs){
                 print k, (runs[k]+0), (ok2[k]+0), (one[k]+0), (zero[k]+0)
               }}' OFS='\t' "$GLOBAL_AUDIT" | sort > "${OUTPUT_DIR}/._tmp_proj_counts.tsv" || true

  awk -F'\t' 'NR>1 && $7=="recover" && $8=="OK"{rec[$2]++} END{for(k in rec) print k,rec[k]}' OFS='\t' "$GLOBAL_ACTIONS" | sort > "${OUTPUT_DIR}/._tmp_proj_rec.tsv" || true

  join -t $'\t' -a1 -e0 -o '1.1,1.2,1.3,1.4,1.5,2.2' "${OUTPUT_DIR}/._tmp_proj_counts.tsv" "${OUTPUT_DIR}/._tmp_proj_rec.tsv" \
    | awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6}' \
    > "$GLOBAL_SUMMARY2"

  rm -f "${OUTPUT_DIR}/._tmp_proj_counts.tsv" "${OUTPUT_DIR}/._tmp_proj_rec.tsv" 2>/dev/null || true
}

# -------- Scan runs present on disk --------------------------------
# Layout: OUTPUT_DIR/study/sample/experiment/run/*.fastq.gz
runs_found=0

while IFS= read -r -d '' run_dir; do
  runs_found=$((runs_found+1))
  #rel="${run_dir#$OUTPUT_DIR/}"
  rel="${run_dir#$BASE_DIR_SCAN/}"
  # rel = study/sample/exp/run
  IFS='/' read -r study sample exp run <<< "$rel"
  if [[ ! "$run" =~ ^[SED]RR[0-9]+$ ]]; then
  record_error "UNKNOWN" "$study" "$sample" "$exp" "$run" "parse" "Run not parsed as accession from path: $rel"
  continue
  fi
  project_id="UNKNOWN"
  # try to find a matching project_id from EXPECTED_INDEX
  if [[ -f "$EXPECTED_INDEX" ]]; then
    project_id="$(awk -F'\t' -v s="$study" -v sa="$sample" -v e="$exp" -v r="$run" \
      'NR>1 && $3==s && $4==sa && $5==e && $6==r {print $2; exit}' "$EXPECTED_INDEX" 2>/dev/null || echo "UNKNOWN")"
  fi

  # count fastqs
  #mapfile -t fq_files < <(find "$run_dir" -maxdepth 1 -type f -name '*.fastq.gz' -printf '%f\n' | LC_ALL=C sort) uncomment when do with not symlink
  mapfile -t fq_files < <(find -L "$run_dir" -maxdepth 1 -type f -name '*.fastq.gz' -printf '%f\n' | LC_ALL=C sort)
  fq_count="${#fq_files[@]}"
  files_joined="$(safe_join_semicolons "${fq_files[@]:-}")"
  echo -e "${project_id}\t${study}\t${sample}\t${exp}\t${run}\t${fq_count}\t${files_joined}" >> "$GLOBAL_AUDIT"

  # If 2+ FASTQs, this run is “OK” for most pipelines; continue
  if (( fq_count >= 2 )); then
    continue
  fi

  # For 0 or 1 files: inspect more deeply
  # 1) ENA metadata
  if have curl; then
    ena_line="$(ena_query "$run" | tail -n +2 || true)"
    # If nothing comes back, still proceed with recovery heuristics
    if [[ -n "$ena_line" ]]; then
      IFS=$'\t' read -r r_acc layout strategy plat model rcount bcount fastq_ftp submitted_ftp submitted_fmt <<< "$ena_line"
      echo -e "${r_acc}\t${layout}\t${strategy}\t${plat}\t${model}\t${rcount}\t${bcount}\t${fastq_ftp}\t${submitted_ftp}\t${submitted_fmt}" >> "$GLOBAL_ENA_CHECK"
    else
      layout="NA"; strategy="NA"; plat="NA"; model="NA"; fastq_ftp=""
      record_error "$project_id" "$study" "$sample" "$exp" "$run" "ENA" "No ENA record returned"
    fi
  else
    layout="NA"; strategy="NA"; plat="NA"; model="NA"; fastq_ftp=""
    record_error "$project_id" "$study" "$sample" "$exp" "$run" "ENA" "curl missing; skipped ENA check"
  fi

  # 2) If there is exactly 1 FASTQ, probe interleaving cheaply
  if (( fq_count == 1 )); then
    fq="${run_dir}/${fq_files[0]}"
    if [[ -s "$fq" ]] && have zcat; then
      inter="$(probe_interleave "$fq")"
      case "$inter" in
        yes)  echo -e "${study}\t${sample}\t${exp}\t${run}\t${fq_files[0]}\t'header'\tyes\t'headers contain /1 and /2'" >> "$GLOBAL_INTERLEAVE" ;;
        maybe)echo -e "${study}\t${sample}\t${exp}\t${run}\t${fq_files[0]}\t'header'\tmaybe\t'CASAVA read numbers present'" >> "$GLOBAL_INTERLEAVE" ;;
        *)    echo -e "${study}\t${sample}\t${exp}\t${run}\t${fq_files[0]}\t'header'\tno\t'no /1 or /2 detected in headers'" >> "$GLOBAL_INTERLEAVE" ;;
      esac
    else
      record_error "$project_id" "$study" "$sample" "$exp" "$run" "probe" "Cannot probe interleaving (zcat missing or file empty)"
    fi
  fi

  # 3) Decide whether to attempt recovery via SRA
  attempt_recover="false"
  reason=""

  if [[ "${DO_RECOVER}" == "true" ]]; then
    # If Illumina & layout says PAIRED but fq_count < 2, very likely ENA exposed only read #3
    if (( fq_count < 2 )) && [[ "${layout:-}" == "PAIRED" ]] && [[ "${plat:-}" == ILLUMINA || "${model:-}" == *Illumina* ]]; then
      if [[ "${ASSUME_10X_IF_SINGLE_AND_PAIRED}" == "true" ]]; then
        attempt_recover="true"; reason="Illumina PAIRED but only ${fq_count} FASTQ on disk"
      fi
    fi
    # If zero files at all → try to recover
    if (( fq_count == 0 )); then
      attempt_recover="true"; reason="No FASTQ present"
    fi
  fi

# Recovery Block
if [[ "$attempt_recover" == "true" ]] && have prefetch && have fasterq-dump; then
  record_action "$project_id" "$study" "$sample" "$exp" "$run" "recover" "START" "$reason"

  mkdir -p "$PREFETCH_ROOT" "$RECOVERY_TMP"

  if prefetch --output-directory "$PREFETCH_ROOT" "$run" >/dev/null 2>&1; then
    # Build fasterq-dump args
    dump_args=( "--split-files" "--threads" "$THREADS" "-O" "$RECOVERY_TMP" )
    [[ "$INCLUDE_TECH" == "true" ]] && dump_args=( "--include-technical" "${dump_args[@]}" )

    # Optional: point fasterq temp to RECOVERY_TMP to avoid /tmp pressure
    # dump_args+=( "--temp" "$RECOVERY_TMP" )

    if fasterq-dump "${dump_args[@]}" "$run" >/dev/null 2>&1; then
      # Compress outputs (pigz if available). Be careful: don't quote the glob.
      shopt -s nullglob
      to_zip=( ${RECOVERY_TMP}/${run}_*.fastq )
      if (( ${#to_zip[@]} )); then
        if have pigz; then pigz -p "$THREADS" "${to_zip[@]}"; else gzip "${to_zip[@]}"; fi
      fi
      shopt -u nullglob

      # Move into run_dir without clobbering existing files/symlinks
      moved=0
      shopt -s nullglob
      for f in ${RECOVERY_TMP}/${run}_*.fastq.gz; do
        [[ -f "$f" ]] || continue
        base="$(basename "$f")"            # SRRxxxxx_1.fastq.gz, _2, _3 …
        new="$base"                         # keep SRR-style names (pipeline-friendly)
        dest="$run_dir/$new"
        if [[ -L "$dest" || -e "$dest" ]]; then
          record_action "$project_id" "$study" "$sample" "$exp" "$run" "recover" "SKIP" "dest exists ($dest); not overwriting"
        else
          if mv -n "$f" "$dest"; then
            moved=$((moved+1))
          fi
        fi
      done
      shopt -u nullglob

      if (( moved > 0 )); then
        record_action "$project_id" "$study" "$sample" "$exp" "$run" "recover" "OK" "moved ${moved} file(s) into run_dir"
        # Optional: capture stats if sra-stat is present
        if have sra-stat; then
          layout_line="$(sra-stat -s "$run" 2>/dev/null | tr -s ' ' | sed -n '1,5p' | paste -sd ';' -)"
          # You could append layout_line to GLOBAL_SRA_LAYOUT if you want.
        fi
      else
        record_action "$project_id" "$study" "$sample" "$exp" "$run" "recover" "PARTIAL" "dumped but could not move/locate .fastq.gz"
      fi
    else
      record_action "$project_id" "$study" "$sample" "$exp" "$run" "recover" "FAIL" "fasterq-dump failed"
      record_error  "$project_id" "$study" "$sample" "$exp" "$run" "fasterq-dump" "failed"
    fi
  else
    record_action "$project_id" "$study" "$sample" "$exp" "$run" "recover" "FAIL" "prefetch failed"
    record_error  "$project_id" "$study" "$sample" "$exp" "$run" "prefetch" "failed"
  fi

  # Clean temp of this run (ignore if none)
  rm -f ${RECOVERY_TMP}/${run}_*.fastq.gz 2>/dev/null || true

  # Re-count in run_dir (symlink-aware) and probe again if changed
  mapfile -t fq_files2 < <(find -L "$run_dir" -maxdepth 1 -type f -name '*.fastq.gz' -printf '%f\n' | LC_ALL=C sort)
  fq_count2="${#fq_files2[@]}"
  if (( fq_count2 != fq_count )); then
    echo -e "${project_id}\t${study}\t${sample}\t${exp}\t${run}\t${fq_count2}\t$(safe_join_semicolons "${fq_files2[@]}")" >> "$GLOBAL_AUDIT"
    for f in "${fq_files2[@]}"; do
      if [[ -s "$run_dir/$f" ]]; then
        inter2="$(probe_interleave "$run_dir/$f")"
        echo -e "${study}\t${sample}\t${exp}\t${run}\t${f}\t'header'\t${inter2}\t'post-recovery probe'" >> "$GLOBAL_INTERLEAVE"
      fi
    done
  fi
elif [[ "$attempt_recover" == "true" ]]; then
  record_action "$project_id" "$study" "$sample" "$exp" "$run" "recover" "SKIP" "prefetch/fasterq-dump not available"
  record_error  "$project_id" "$study" "$sample" "$exp" "$run" "recover" "tools missing"
fi

done < <(find "$BASE_DIR_SCAN" -mindepth 4 -maxdepth 4 -type d -print0 | LC_ALL=C sort -z)

# -------- Also catch runs that are "expected" but missing their directory ----
if [[ -f "$EXPECTED_INDEX" ]]; then
  # Build a unique list of expected (study/sample/exp/run)
  awk -F'\t' 'NR>1{print $3"/"$4"/"$5"/"$6}' "$EXPECTED_INDEX" | sort -u \
  | while read -r rel; do
      [[ -z "$rel" ]] && continue
      #run_dir="$OUTPUT_DIR/$rel" uncomment when do with not symlink
      run_dir="$BASE_DIR_SCAN/$rel"
      if [[ ! -d "$run_dir" ]]; then
        IFS='/' read -r study sample exp run <<< "$rel"
        project_id="$(awk -F'\t' -v s="$study" -v sa="$sample" -v e="$exp" -v r="$run" \
          'NR>1 && $3==s && $4==sa && $5==e && $6==r {print $2; exit}' "$EXPECTED_INDEX" 2>/dev/null || echo "UNKNOWN")"
        echo -e "${project_id}\t${study}\t${sample}\t${exp}\t${run}\t0\t" >> "$GLOBAL_AUDIT"
      fi
    done
fi

# -------- Build per-project roll-up --------------------------------
summarize_by_project

log "Postcheck complete."
log "Wrote:"
log "  GLOBAL_AUDIT:      $GLOBAL_AUDIT"
log "  GLOBAL_ENA_CHECK:  $GLOBAL_ENA_CHECK"
log "  GLOBAL_INTERLEAVE: $GLOBAL_INTERLEAVE"
log "  GLOBAL_SRA_LAYOUT: $GLOBAL_SRA_LAYOUT (best-effort)"
log "  GLOBAL_ACTIONS:    $GLOBAL_ACTIONS"
log "  GLOBAL_ERRORS:     $GLOBAL_ERRORS"
log "  SUMMARY (project): $GLOBAL_SUMMARY2"
