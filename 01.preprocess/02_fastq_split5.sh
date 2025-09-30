#!/usr/bin/env bash
set -Euo pipefail
IFS=$'\n\t'

# ========================
# Config (edit as needed)
# ========================
BASE_DIR="/mnt/18T/chibao/gliomas/data/metadata/official/split5.1"   # Step-1 output root
OUTPUT_DIR="/mnt/18T/chibao/gliomas/data/fastq/ena/split5.1"         # where FASTQs are saved

# Logs & reports (GLOBAL files accumulate across all splits)
GLOBAL_SUMMARY="$OUTPUT_DIR/_download_summary.tsv"          # per-file event log (append-only)
PROJECT_AUDIT="$OUTPUT_DIR/_project_audit.tsv"              # READY/TSV usage per project (append-only)
METADATA_INV="$OUTPUT_DIR/_metadata_inventory.tsv"          # all TSVs seen & whether used (append-only)
EXPECTED_INDEX="$OUTPUT_DIR/_expected_files.tsv"            # flattened expected files (append-only)
COVERAGE_REPORT="$OUTPUT_DIR/_coverage_report.tsv"          # per-project coverage stats (regenerated each run)
MISSING_REPORT="$OUTPUT_DIR/_missing_files.tsv"             # expected - covered (regenerated each run)
MISMATCH_REPORT="$OUTPUT_DIR/_md5_mismatches.tsv"           # md5 mismatches (regenerated each run)
ORPHANS_REPORT="$OUTPUT_DIR/_orphans.tsv"                   # files in OUTPUT with no metadata (regenerated each run)

USER_AGENT="ena-fastq-downloader/2.2"

# aria2c tuning (used if available)
CONNS=4               # was 6; fewer per-server conns is gentler
SPLITS=4              # was 6; keep <= CONNS
SPLIT_MIN="8M"        # was 1M; larger pieces = fewer tiny segments
ARIA_MAX_TRIES=8      # aria2's own tries (on top of our outer loop)
ARIA_RETRY_WAIT=5     # seconds between aria2 retries
ARIA_TIMEOUT=60       # socket timeout
ARIA_CONN_TIMEOUT=15  # connect timeout
ARIA_FILE_ALLOC="prealloc"  # better on HDDs; use "falloc" on XFS/ext4+SSD


MAX_RETRIES=6
BACKOFF_INITIAL=3

# === MD5 & mismatch handling ===
MISMATCH_ACTION="quarantine"   # one of: keep | delete | quarantine
QUARANTINE_ROOT="$OUTPUT_DIR/_quarantine"

# Resume mode: when true, skip work already marked done and avoid duplicate appends
RESUME=false

# ========================
# Split ID handling
# ========================
SPLIT_ID="${SPLIT_ID:-}"  # allow env var
parse_args() {
  while (( "$#" )); do
    case "$1" in
      --split-id)
        SPLIT_ID="${2:-}"; shift 2 ;;
      --base-dir)
        BASE_DIR="${2:-}"; shift 2 ;;
      --output-dir)
        OUTPUT_DIR="${2:-}"; shift 2 ;;
      --resume)
        RESUME=true; shift ;;
      *)
        echo "WARN: Unknown arg: $1" >&2; shift ;;
    esac
  done
}

deduce_split_id() {
  if [[ -n "${SPLIT_ID:-}" ]]; then
    printf '%s' "$SPLIT_ID"; return
  fi
  local b; b="$(basename "$BASE_DIR" || true)"
  if [[ "$b" =~ ^split[0-9A-Za-z._-]*$ ]]; then
    printf '%s' "$b"
  else
    printf 'unknown'
  fi
}

# ========================
# Logging
# ========================
timestamp(){ date '+%Y-%m-%d %H:%M:%S%z'; }
log(){ printf '[%s] %s\n' "$(timestamp)" "$*" >&2; }

# Non-fatal trap: log unhandled errors but don’t exit whole job
trap 'log "ERROR: Unhandled error occurred."' ERR

# ========================
# Init (append-safe)
# ========================
init_files(){
  mkdir -p "$OUTPUT_DIR"

  # Add headers ONLY if files do not exist yet (so multiple splits append safely)
  if [[ ! -f "$GLOBAL_SUMMARY" ]]; then
    echo -e "split_id\tts\tproject_id\tstudy\tsample\texperiment\trun\tfile\tstatus\tmd5_expected\tmd5_actual\tbytes\tmessage" > "$GLOBAL_SUMMARY"
  fi
  if [[ ! -f "$PROJECT_AUDIT" ]]; then
    echo -e "split_id\tts\tproject_id\tready\tused_tsv\tstatus\tmessage" > "$PROJECT_AUDIT"
  fi
  if [[ ! -f "$METADATA_INV" ]]; then
    echo -e "split_id\tts\tproject_id\ttsv_path\tusage\tmessage" > "$METADATA_INV"
  fi
  if [[ ! -f "$EXPECTED_INDEX" ]]; then
    echo -e "split_id\tproject_id\tstudy\tsample\texperiment\trun\tfile\turl\tmd5_expected\ttarget_path" > "$EXPECTED_INDEX"
  fi
}

# Append helpers (GLOBAL)
append_summary(){
  local split_id="$1" ts="$2" project="$3" study="$4" sample="$5" exp="$6" run="$7" file="$8" status="$9" md5_exp="${10}" md5_act="${11}" bytes="${12}" msg="${13}"
  echo -e "${split_id}\t${ts}\t${project}\t${study}\t${sample}\t${exp}\t${run}\t${file}\t${status}\t${md5_exp}\t${md5_act}\t${bytes}\t${msg}" >> "$GLOBAL_SUMMARY"
}

append_project_audit(){
  local split_id="$1" ts="$2" project_id="$3" ready="$4" used_tsv="$5" status="$6" msg="$7"
  echo -e "${split_id}\t${ts}\t${project_id}\t${ready}\t${used_tsv}\t${status}\t${msg}" >> "$PROJECT_AUDIT"
}

append_metadata_inv(){
  local split_id="$1" ts="$2" project_id="$3" tsv_path="$4" usage="$5" message="$6"
  echo -e "${split_id}\t${ts}\t${project_id}\t${tsv_path}\t${usage}\t${message}" >> "$METADATA_INV"
}

append_expected_index(){
  local split_id="$1" project_id="$2" study="$3" sample="$4" exp="$5" run="$6" file="$7" url="$8" md5_exp="$9" target="${10}"
  echo -e "${split_id}\t${project_id}\t${study}\t${sample}\t${exp}\t${run}\t${file}\t${url}\t${md5_exp}\t${target}" >> "$EXPECTED_INDEX"
}

# Per-run log (in each run folder)
init_run_log(){
  local rundir="$1"
  local runlog="$rundir/_run_log.tsv"
  if [[ ! -f "$runlog" ]]; then
    echo -e "split_id\tts\tstatus\tfile\tmd5_expected\tmd5_actual\tbytes\tmessage" > "$runlog"
  fi
}
append_run_log(){
  local rundir="$1" split_id="$2" ts="$3" status="$4" file="$5" md5_exp="$6" md5_act="$7" bytes="$8" msg="$9"
  echo -e "${split_id}\t${ts}\t${status}\t${file}\t${md5_exp}\t${md5_act}\t${bytes}\t${msg}" >> "$rundir/_run_log.tsv"
}

# ========================
# Helpers
# ========================
normalize_url(){
  local raw="$1"
  if [[ "$raw" =~ ^https?:// || "$raw" =~ ^ftp:// ]]; then
    printf '%s' "$raw"
  else
    printf 'https://%s' "$raw"
  fi
}

md5_of(){ [[ -f "$1" ]] && md5sum "$1" | awk '{print $1}' || echo ""; }
bytes_of(){ [[ -f "$1" ]] && { stat -c '%s' "$1" 2>/dev/null || wc -c <"$1"; } || echo 0; }

sleep_backoff(){
  local attempt="$1"
  local delay=$(( BACKOFF_INITIAL * (attempt - 1) ))
  sleep "$delay"
}

handle_mismatch_file(){
  local path="$1"
  case "$MISMATCH_ACTION" in
    keep) : ;;
    delete)
      rm -f -- "$path" "$path.aria2" 2>/dev/null || true
      ;;
    quarantine)
      mkdir -p "$QUARANTINE_ROOT"
      local rel="${path#$OUTPUT_DIR/}"
      local qpath="$QUARANTINE_ROOT/$rel"
      mkdir -p "$(dirname "$qpath")"
      mv -f -- "$path" "$qpath" 2>/dev/null || true
      rm -f -- "$path.aria2" 2>/dev/null || true
      echo "$qpath"
      ;;
  esac
}

# ========================
# Resume helpers
# ========================
# Build a composite key for a file row
key_of(){ # project study sample exp run file
  printf '%s|%s|%s|%s|%s|%s' "$1" "$2" "$3" "$4" "$5" "$6"
}

# DONE set from GLOBAL_SUMMARY for statuses that mean the file is completed
declare -A DONE_KEYS=()
build_done_set(){
  if [[ -f "$GLOBAL_SUMMARY" ]]; then
    while IFS= read -r k; do
      DONE_KEYS["$k"]=1
    done < <(
      awk -F'\t' 'NR>1 && ($9=="success" || $9=="skipped_existing"){print $3"|"$4"|"$5"|"$6"|"$7"|"$8}' "$GLOBAL_SUMMARY" \
      | sort -u
    )
  fi
}

# Unique append helpers (avoid duplicates when RESUME=true)
append_unique(){
  local f="$1" L="$2"
  if [[ -f "$f" ]] && grep -Fqx -- "$L" "$f"; then
    return 0
  fi
  echo "$L" >> "$f"
}

append_expected_index_unique(){
  local split_id="$1" project_id="$2" study="$3" sample="$4" exp="$5" run="$6" file="$7" url="$8" md5_exp="$9" target="${10}"
  local line="${split_id}\t${project_id}\t${study}\t${sample}\t${exp}\t${run}\t${file}\t${url}\t${md5_exp}\t${target}"
  append_unique "$EXPECTED_INDEX" "$line"
}

append_project_audit_unique(){
  local split_id="$1" ts="$2" project_id="$3" ready="$4" used_tsv="$5" status="$6" msg="$7"
  local line="${split_id}\t${ts}\t${project_id}\t${ready}\t${used_tsv}\t${status}\t${msg}"
  $RESUME && append_unique "$PROJECT_AUDIT" "$line" || echo -e "$line" >> "$PROJECT_AUDIT"
}

append_metadata_inv_unique(){
  local split_id="$1" ts="$2" project_id="$3" tsv_path="$4" usage="$5" message="$6"
  local line="${split_id}\t${ts}\t${project_id}\t${tsv_path}\t${usage}\t${message}"
  $RESUME && append_unique "$METADATA_INV" "$line" || echo -e "$line" >> "$METADATA_INV"
}

# ========================
# Downloaders
# ========================
download_file(){
  # url, outdir, outfile, md5_expected (optional)
  local url="$1" outdir="$2" outfile="$3" md5_exp="${4:-}"

  mkdir -p "$outdir"
  local dest="$outdir/$outfile"
  local tmp="$dest.part"

  if command -v aria2c >/dev/null 2>&1; then
    # Build common options
    local aopts=(
      --user-agent="$USER_AGENT"
      -x "$CONNS" -s "$SPLITS" -k "$SPLIT_MIN"
      -j 1 -c
      --max-tries="$ARIA_MAX_TRIES"
      --retry-wait="$ARIA_RETRY_WAIT"
      --timeout="$ARIA_TIMEOUT"
      --connect-timeout="$ARIA_CONN_TIMEOUT"
      --file-allocation="$ARIA_FILE_ALLOC"
      --continue=true
      --check-integrity=true
      --auto-file-renaming=false
      --http-no-cache=true
    )
    # If we know MD5, let aria2 verify it (and retry if it fails)
    if [[ -n "$md5_exp" ]]; then
      aopts+=( "--checksum=md5=${md5_exp}" )
    fi

    aria2c "${aopts[@]}" -d "$outdir" -o "$outfile" "$url"
    return $?
  fi

  if command -v curl >/dev/null 2>&1; then
    # Curl fallback: resume + retries. md5 is checked after download in our script.
    curl -A "$USER_AGENT" -L --fail --retry 3 --retry-delay 3 -C - -o "$tmp" "$url" && mv -f "$tmp" "$dest"
    return $?
  fi

  if command -v wget >/dev/null 2>&1; then
    # Wget fallback: resume. md5 checked later by script.
    wget --user-agent="$USER_AGENT" -c -O "$dest" "$url"
    return $?
  fi

  log "ERROR: No downloader found (aria2c|curl|wget)."
  return 127
}


check_md5(){
  local file="$1" expected="$2"
  [[ -z "$expected" ]] && return 1
  local actual; actual="$(md5_of "$file")"
  [[ -n "$actual" && "$actual" == "$expected" ]]
}

# ========================
# Core: process one FASTQ
# ========================
process_one_file(){
  local split_id="$1" project="$2" study="$3" sample="$4" exp="$5" run="$6" url_raw="$7" md5_exp="$8"

  local url; url="$(normalize_url "$url_raw")"
  local fname; fname="$(basename "$url")"
  local rundir="$OUTPUT_DIR/$study/$sample/$exp/$run"
  local dest="$rundir/$fname"

  mkdir -p "$rundir"
  init_run_log "$rundir"

  # Skip if already good on disk
  if [[ -f "$dest" ]] && check_md5 "$dest" "$md5_exp"; then
    local bytes; bytes="$(bytes_of "$dest")"
    append_run_log "$rundir" "$split_id" "$(timestamp)" "skipped_existing" "$fname" "$md5_exp" "$md5_exp" "$bytes" "Already present and checksum OK"
    append_summary   "$split_id" "$(timestamp)" "$project" "$study" "$sample" "$exp" "$run" "$fname" "skipped_existing" "$md5_exp" "$md5_exp" "$bytes" "Already present"
    return 0
  fi

  # Be quieter in resume mode for "started/downloading" spam
  if ! $RESUME; then
    append_run_log "$rundir" "$split_id" "$(timestamp)" "started" "$fname" "$md5_exp" "" "0" "Begin download"
    append_summary   "$split_id" "$(timestamp)" "$project" "$study" "$sample" "$exp" "$run" "$fname" "started" "$md5_exp" "" "0" "Begin download"
  fi

  local attempt=1
  while (( attempt <= MAX_RETRIES )); do
    if ! $RESUME; then
      append_run_log "$rundir" "$split_id" "$(timestamp)" "downloading" "$fname" "$md5_exp" "" "0" "Attempt ${attempt}/${MAX_RETRIES}: $url"
    fi
    if download_file "$url" "$rundir" "$fname" "$md5_exp"; then
      local actual; actual="$(md5_of "$dest")"
      local bytes;  bytes="$(bytes_of "$dest")"
      if [[ -n "$md5_exp" && -n "$actual" && "$actual" == "$md5_exp" ]]; then
        append_run_log "$rundir" "$split_id" "$(timestamp)" "checksum_ok" "$fname" "$md5_exp" "$actual" "$bytes" "Download and MD5 verified"
        append_summary   "$split_id" "$(timestamp)" "$project" "$study" "$sample" "$exp" "$run" "$fname" "success" "$md5_exp" "$actual" "$bytes" "OK"
        return 0
      else
        append_run_log "$rundir" "$split_id" "$(timestamp)" "checksum_fail" "$fname" "$md5_exp" "$actual" "$bytes" "MD5 mismatch on attempt ${attempt}/${MAX_RETRIES}"
        append_summary   "$split_id" "$(timestamp)" "$project" "$study" "$sample" "$exp" "$run" "$fname" "md5_mismatch" "$md5_exp" "$actual" "$bytes" "Attempt ${attempt}/${MAX_RETRIES}: mismatch"
      fi
    else
      append_run_log "$rundir" "$split_id" "$(timestamp)" "error" "$fname" "$md5_exp" "" "0" "Downloader failed; will retry if attempts remain"
    fi
    (( attempt++ ))
    if (( attempt <= MAX_RETRIES )); then sleep_backoff $((attempt-1)); fi
  done

  # Final state after retries
  local final_md5; final_md5="$(md5_of "$dest")"
  local bytes; bytes="$(bytes_of "$dest")"

  local qmsg=""
  if [[ -f "$dest" ]]; then
    local qpath; qpath="$(handle_mismatch_file "$dest")" || true
    if [[ -n "$qpath" ]]; then
      qmsg=" | quarantined: $qpath"
    elif [[ "$MISMATCH_ACTION" == "delete" ]]; then
      qmsg=" | deleted"
    fi
  fi

  append_run_log "$rundir" "$split_id" "$(timestamp)" "failed" "$fname" "$md5_exp" "$final_md5" "$bytes" "Exhausted retries; MD5 mismatch or download error${qmsg}"
  append_summary   "$split_id" "$(timestamp)" "$project" "$study" "$sample" "$exp" "$run" "$fname" "md5_mismatch" "$md5_exp" "$final_md5" "$bytes" "Exhausted retries; MD5 mismatch${qmsg}"
  return 1
}

# ========================
# Parse one project TSV and download all rows
# (also builds EXPECTED_INDEX and metadata audits)
# ========================
process_project(){
  local split_id="$1" project_dir="$2"
  local project_id; project_id="$(basename "$project_dir")"
  local FAIL_COUNT=0

  local ready="no"
  [[ -f "$project_dir/READY" ]] && ready="yes"

  # prefer latest.tsv; fallback to ${project_id}.tsv
  local used_tsv="" status="" msg=""
  if [[ -f "$project_dir/latest.tsv" ]]; then
    used_tsv="$project_dir/latest.tsv"; status="used"; msg="OK"
  elif [[ -f "$project_dir/${project_id}.tsv" ]]; then
    used_tsv="$project_dir/${project_id}.tsv"; status="used"; msg="OK (no symlink)"
  else
    used_tsv=""; status="skipped"; msg="No metadata TSV"
  fi

  append_project_audit_unique "$split_id" "$(timestamp)" "$project_id" "$ready" "$used_tsv" "$status" "$msg"

  # inventory all *.tsv in the folder and mark usage
  while IFS= read -r -d '' t; do
    local usage="ignored" reason="not selected"
    [[ "$t" == "$used_tsv" ]] && { usage="used"; reason="primary metadata"; }
    [[ "$(basename "$t")" == "_run_log.tsv" ]] && { usage="ignored"; reason="run log"; }
    append_metadata_inv_unique "$split_id" "$(timestamp)" "$project_id" "$t" "$usage" "$reason"
  done < <(find "$project_dir" -maxdepth 1 -type f -name "*.tsv" -print0 | sort -z)

  # skip if no READY or no TSV
  if [[ "$ready" != "yes" || -z "$used_tsv" ]]; then
    log "Project $project_id: ready=$ready used_tsv='$used_tsv' — skipping downloads."
    return 0
  fi

  log "Project: $project_id | TSV: $used_tsv"

  # Header -> index map
  local header; header="$(head -n1 "$used_tsv")"
  IFS=$'\t' read -r -a cols <<< "$header"
  declare -A H=()
  for i in "${!cols[@]}"; do H["${cols[$i]}"]="$i"; done

  for must in study_accession sample_accession experiment_accession run_accession fastq_ftp fastq_md5; do
    if [[ -z "${H[$must]+x}" ]]; then
      log "WARNING: Missing required column '$must' in $used_tsv. Skipping this TSV."
      return 0
    fi
  done

  # Iterate rows, populate EXPECTED_INDEX, then download
  while IFS=$'\t' read -r -a fields; do
    [[ "${#fields[@]}" -eq 0 ]] && continue

    local study="${fields[${H[study_accession]}]:-NA}"
    local sample="${fields[${H[sample_accession]}]:-NA}"
    local exp="${fields[${H[experiment_accession]}]:-NA}"
    local run="${fields[${H[run_accession]}]:-NA}"
    local url_list="${fields[${H[fastq_ftp]}]:-}"
    local md5_list="${fields[${H[fastq_md5]}]:-}"
    [[ -z "$url_list" ]] && continue

    IFS=';' read -r -a urls <<< "$url_list"
    IFS=';' read -r -a md5s <<< "$md5_list"
    if (( ${#md5s[@]} < ${#urls[@]} )); then
      for ((k=${#md5s[@]}; k<${#urls[@]}; k++)); do md5s+=(""); done
      log "NOTE: $project_id ${study}/${sample}/${exp}/${run} has ${#urls[@]} URLs but fewer MD5s."
    fi

    local rundir="$OUTPUT_DIR/$study/$sample/$exp/$run"
    mkdir -p "$rundir"
    init_run_log "$rundir"
    append_run_log "$rundir" "$split_id" "$(timestamp)" "queued" "-" "-" "-" "0" "Row queued for processing"

    for idx in "${!urls[@]}"; do
      local url_raw="${urls[$idx]}"
      local md5_exp="${md5s[$idx]}"
      local url="$(normalize_url "$url_raw")"
      local fname="$(basename "$url")"
      local target="$rundir/$fname"

      append_expected_index_unique "$split_id" "$project_id" "$study" "$sample" "$exp" "$run" "$fname" "$url" "$md5_exp" "$target"

      # If resume is ON and this key is already DONE, skip entirely
      if $RESUME; then
        local comp_key
        comp_key="$(key_of "$project_id" "$study" "$sample" "$exp" "$run" "$fname")"
        if [[ -n "${DONE_KEYS[$comp_key]:-}" ]]; then
          append_run_log "$rundir" "$split_id" "$(timestamp)" "resume_skip_done" "$fname" "$md5_exp" "$md5_exp" "$(bytes_of "$target")" "Already marked done in summary"
          continue
        fi
      fi

      if process_one_file "$split_id" "$project_id" "$study" "$sample" "$exp" "$run" "$url_raw" "$md5_exp"; then
        :
      else
        FAIL_COUNT=$((FAIL_COUNT+1))
      fi
    done
  done < <(tail -n +2 "$used_tsv")

  log "Project $project_id finished with ${FAIL_COUNT} file failure(s) recorded for split_id=${split_id}."
}

# ========================
# Coverage calculators (recomputed fresh each run)
# ========================
build_mismatch_report(){
  if grep -q $'\tmd5_mismatch\t' "$GLOBAL_SUMMARY" 2>/dev/null; then
    awk -F'\t' 'BEGIN{OFS="\t"}
      NR==1{print "split_id","ts","project_id","study","sample","experiment","run","file","md5_expected","md5_actual","bytes","message"; next}
      $9=="md5_mismatch"{print $1,$2,$3,$4,$5,$6,$7,$8,$10,$11,$12,$13}' "$GLOBAL_SUMMARY" > "$MISMATCH_REPORT"
    log "MD5 mismatch report written: $MISMATCH_REPORT"
  fi
}

build_coverage_reports(){
  # Per-project coverage stats
  echo -e "project_id\texpected\tdownloaded_ok\tskipped_existing\tmismatched\tmissing" > "$COVERAGE_REPORT"

  # keys: project|study|sample|exp|run|file  (EXPECTED_INDEX now has split_id in col1)
  awk -F'\t' 'NR>1{print $2"|"$3"|"$4"|"$5"|"$6"|"$7}' "$EXPECTED_INDEX" | sort -u > "$OUTPUT_DIR/.expected.keys"
  awk -F'\t' 'NR>1 && ($9=="success" || $9=="skipped_existing"){print $3"|"$4"|"$5"|"$6"|"$7"|"$8}' "$GLOBAL_SUMMARY" | sort -u > "$OUTPUT_DIR/.covered.keys"
  awk -F'\t' 'NR>1 && $9=="md5_mismatch"{print $3"|"$4"|"$5"|"$6"|"$7"|"$8}' "$GLOBAL_SUMMARY" | sort -u > "$OUTPUT_DIR/.mismatch.keys" || true

  # Missing = expected - covered
  comm -23 "$OUTPUT_DIR/.expected.keys" "$OUTPUT_DIR/.covered.keys" > "$OUTPUT_DIR/.missing.keys"

  # Build _missing_files.tsv  (EXPECTED_INDEX columns shifted by split_id)
  echo -e "project_id\tstudy\tsample\texperiment\trun\tfile\turl\tmd5_expected\ttarget_path" > "$MISSING_REPORT"
  while IFS= read -r k; do
    proj="${k%%|*}"; rest="${k#*|}"
    IFS='|' read -r study sample exp run file <<< "$rest"
    awk -F'\t' -v p="$proj" -v s="$study" -v sa="$sample" -v e="$exp" -v r="$run" -v f="$file" \
      'NR>1 && $2==p && $3==s && $4==sa && $5==e && $6==r && $7==f {print $2,$3,$4,$5,$6,$7,$8,$9,$10}' "$EXPECTED_INDEX"
  done < "$OUTPUT_DIR/.missing.keys" >> "$MISSING_REPORT"

  # Summaries
  awk -F'\t' 'NR>1{exp[$2]++} END{for(p in exp) print p,exp[p]}' OFS='\t' "$EXPECTED_INDEX" | sort > "$OUTPUT_DIR/.exp.counts"
  awk -F'\t' 'NR>1 && ($9=="success" || $9=="skipped_existing"){ok[$3]++} END{for(p in ok) print p,ok[p]}' OFS='\t' "$GLOBAL_SUMMARY" | sort > "$OUTPUT_DIR/.ok.counts"
  awk -F'\t' 'NR>1 && $9=="skipped_existing"{sk[$3]++} END{for(p in sk) print p,sk[p]}' OFS='\t' "$GLOBAL_SUMMARY" | sort > "$OUTPUT_DIR/.sk.counts"
  awk -F'\t' 'NR>1{mm[$1]++} END{for(p in mm) print p,mm[p]}' OFS='\t' "$MISSING_REPORT" | sort > "$OUTPUT_DIR/.miss.counts"
  awk -F'|' '{mm[$1]++} END{for(p in mm) print p,mm[p]}' "$OUTPUT_DIR/.mismatch.keys" 2>/dev/null | sort > "$OUTPUT_DIR/.md5m.counts" || true

  join -a1 -e0 -o '1.1,1.2,2.2' -t $'\t' "$OUTPUT_DIR/.exp.counts" "$OUTPUT_DIR/.ok.counts" \
  | join -a1 -e0 -o '1.1,1.2,1.3,2.2' -t $'\t' - "$OUTPUT_DIR/.sk.counts" \
  | join -a1 -e0 -o '1.1,1.2,1.3,1.4,2.2' -t $'\t' - "$OUTPUT_DIR/.md5m.counts" \
  | join -a1 -e0 -o '1.1,1.2,1.3,1.4,1.5,2.2' -t $'\t' - "$OUTPUT_DIR/.miss.counts" \
  | awk -F'\t' 'BEGIN{OFS="\t"; print "project_id","expected","downloaded_ok","skipped_existing","mismatched","missing"} {print}' \
  > "$COVERAGE_REPORT"

  # ORPHANS (on disk without metadata)
  echo -e "study\tsample\texperiment\trun\tfile\tfull_path" > "$ORPHANS_REPORT"
  awk -F'\t' 'NR>1{print $3"|"$4"|"$5"|"$6"|"$7}' "$EXPECTED_INDEX" | sort -u > "$OUTPUT_DIR/.expected.ssersu.keys"
  find "$OUTPUT_DIR" -type f -name "*.fastq.gz" \
    | sed -E "s|^$OUTPUT_DIR/||" \
    | awk -F'/' 'NF>=5{print $1"|"$2"|"$3"|"$4"|"$NF"\t"$0}' \
    | sort -u \
    | while IFS=$'\t' read -r key rel; do
        if ! grep -qx "$key" "$OUTPUT_DIR/.expected.ssersu.keys"; then
          IFS='|' read -r s sa e r f <<< "$key"
          echo -e "$s\t$sa\t$e\t$r\t$f\t$OUTPUT_DIR/$rel" >> "$ORPHANS_REPORT"
        fi
      done
}

# ========================
# Main
# ========================
main() {
  parse_args "$@"
  SPLIT_ID="$(deduce_split_id)"
  init_files

  $RESUME && build_done_set
  $RESUME && log "Resume mode: found ${#DONE_KEYS[@]} completed file keys; will skip them."

  log "Starting ENA download job. split_id=$SPLIT_ID  BASE_DIR=$BASE_DIR  OUTPUT_DIR=$OUTPUT_DIR  resume=$RESUME"

  # Process each project directory under BASE_DIR
  shopt -s nullglob
  project_dirs=( "$BASE_DIR"/*/ )
  shopt -u nullglob

  if (( ${#project_dirs[@]} == 0 )); then
    log "No project directories found in $BASE_DIR"
    # still rebuild reports from whatever exists so far
    build_mismatch_report
    build_coverage_reports
    return 0
  fi

  for project_dir in "${project_dirs[@]}"; do
    process_project "$SPLIT_ID" "$project_dir"
  done

  build_mismatch_report
  build_coverage_reports

  TOTAL_MISMATCH=$(awk -F'\t' 'NR>1 && $9=="md5_mismatch"{c++} END{print c+0}' "$GLOBAL_SUMMARY")
  if [[ -s "$MISSING_REPORT" ]]; then
    TOTAL_MISSING=$(awk -F'\t' 'NR>1{c++} END{print c+0}' "$MISSING_REPORT")
  else
    TOTAL_MISSING=0
  fi
  log "Final: split_id=$SPLIT_ID md5_mismatch=$TOTAL_MISMATCH, missing=$TOTAL_MISSING"

  log "Done."
  log "Global summary:   $GLOBAL_SUMMARY"
  log "Project audit:    $PROJECT_AUDIT"
  log "Metadata index:   $METADATA_INV"
  log "Expected files:   $EXPECTED_INDEX"
  log "Coverage report:  $COVERAGE_REPORT"
  log "Missing files:    $MISSING_REPORT"
  log "Orphans report:   $ORPHANS_REPORT"
}

# Only run main if not in library mode
if [[ -z "${ENA_LIB_ONLY:-}" ]]; then
  main "$@"
fi
