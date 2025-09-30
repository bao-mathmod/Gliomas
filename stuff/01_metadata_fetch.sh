#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

# =========================
# Config (edit as needed)
# =========================
CSV_FILE="/mnt/12T/chibao/data/official_data/glioma_accession.csv"
META_ROOT="/mnt/12T/chibao/data/official_data/metadata"
GLOBAL_MANIFEST="$META_ROOT/_manifest_projects.tsv"
USER_AGENT="ena-meta-fetch/1.0"
MAX_RETRIES=5
RETRY_DELAY=2   # seconds

# ENA fields expected by your Step 2 downloader:
REQUIRED_FIELDS=(
  study_accession
  sample_accession
  experiment_accession
  run_accession
  scientific_name
  fastq_md5
  fastq_ftp
  sra_ftp
)

# =========================
# Logging helpers
# =========================
ts(){ date '+%Y-%m-%d %H:%M:%S%z'; }
log(){ printf '[%s] %s\n' "$(ts)" "$*" >&2; }

init_global_manifest(){
  mkdir -p "$META_ROOT"
  if [[ ! -f "$GLOBAL_MANIFEST" ]]; then
    echo -e "ts\tproject_id\tstatus\tbytes\trows\tmd5\ttsv_path\tmessage" > "$GLOBAL_MANIFEST"
  fi
}

append_global_manifest(){
  local ts_="$1" proj="$2" status="$3" bytes="$4" rows="$5" md5="$6" path="$7" msg="$8"
  echo -e "${ts_}\t${proj}\t${status}\t${bytes}\t${rows}\t${md5}\t${path}\t${msg}" >> "$GLOBAL_MANIFEST"
}

init_run_log(){
  local proj_dir="$1"
  local logf="$proj_dir/_run_log.tsv"
  if [[ ! -f "$logf" ]]; then
    echo -e "ts\tstage\tstatus\tfield\tvalue\tmessage" > "$logf"
  fi
}

append_run_log(){
  local proj_dir="$1" stage="$2" status="$3" field="$4" value="$5" msg="$6"
  echo -e "$(ts)\t${stage}\t${status}\t${field}\t${value}\t${msg}" >> "$proj_dir/_run_log.tsv"
}

# =========================
# Utilities
# =========================
trim(){ awk '{$1=$1;print}' <<<"$*"; }

csv_find_col_idx(){
  # Find the 1-based column index in a CSV header for any of these names
  local header="$1"; shift
  local candidates=("$@")
  local IFS=,
  read -r -a cols <<<"$header"
  for name in "${candidates[@]}"; do
    for i in "${!cols[@]}"; do
      if [[ "$(trim "${cols[$i]}")" == "$name" ]]; then
        echo $((i+1))
        return 0
      fi
    done
  done
  echo 0
}

bytes_of(){ [[ -f "$1" ]] && { stat -c '%s' "$1" 2>/dev/null || wc -c <"$1"; } || echo 0; }
md5_of(){ [[ -f "$1" ]] && md5sum "$1" | awk '{print $1}' || echo ""; }

tsv_has_required_fields(){
  local header="$1"; shift
  local missing=()
  for f in "$@"; do
    if ! awk -F'\t' -v want="$f" 'NR==1{for(i=1;i<=NF;i++) if($i==want) found=1} END{exit(found?0:1)}' <<<"$header"; then
      missing+=("$f")
    fi
  done
  if (( ${#missing[@]} )); then
    printf '%s' "${missing[*]}"
    return 1
  fi
  return 0
}

# =========================
# ENA downloader (with retry)
# =========================
fetch_project_tsv(){
  local project_id="$1"
  local proj_dir="$META_ROOT/$project_id"
  local final_tsv="$proj_dir/${project_id}.tsv"
  local tmp_tsv="$proj_dir/${project_id}.tsv.part"
  local url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${project_id}&result=read_run&fields=$(IFS=,; echo "${REQUIRED_FIELDS[*]}")&format=tsv&download=true&limit=0"

  mkdir -p "$proj_dir"
  init_run_log "$proj_dir"

  append_run_log "$proj_dir" "prepare" "info" "url" "$url" "Built ENA request URL"

  local attempt=1
  while (( attempt <= MAX_RETRIES )); do
    append_run_log "$proj_dir" "download" "start" "attempt" "$attempt" "Fetching TSV"
    if curl -fL -A "$USER_AGENT" --retry "$MAX_RETRIES" --retry-delay "$RETRY_DELAY" \
           --retry-connrefused -s "$url" -o "$tmp_tsv"; then
      append_run_log "$proj_dir" "download" "ok" "attempt" "$attempt" "Downloaded"
      break
    else
      append_run_log "$proj_dir" "download" "error" "attempt" "$attempt" "curl failed; will retry"
      sleep "$RETRY_DELAY"
    fi
    ((attempt++))
  done

  if [[ ! -s "$tmp_tsv" ]]; then
    append_run_log "$proj_dir" "download" "failed" "file" "$tmp_tsv" "Empty or missing file after retries"
    append_global_manifest "$(ts)" "$project_id" "failed" "0" "0" "" "" "Download failed"
    return 1
  fi

  # Validate header + rows
  local header; header="$(head -n1 "$tmp_tsv")"
  if [[ -z "$header" ]]; then
    append_run_log "$proj_dir" "validate" "error" "header" "-" "Missing header"
    append_global_manifest "$(ts)" "$project_id" "failed" "$(bytes_of "$tmp_tsv")" "0" "" "$tmp_tsv" "Missing header"
    return 1
  fi

  if ! tsv_has_required_fields "$header" "${REQUIRED_FIELDS[@]}"; then
    local missing; missing="$(tsv_has_required_fields "$header" "${REQUIRED_FIELDS[@]}"; true || true)"
    append_run_log "$proj_dir" "validate" "error" "missing_fields" "$missing" "Header missing required fields"
    append_global_manifest "$(ts)" "$project_id" "failed" "$(bytes_of "$tmp_tsv")" "0" "" "$tmp_tsv" "Missing fields: $missing"
    return 1
  fi

  # Count data rows (excluding header)
  local rows; rows="$(($(wc -l <"$tmp_tsv") - 1))"
  if (( rows < 0 )); then rows=0; fi
  append_run_log "$proj_dir" "validate" "ok" "rows" "$rows" "TSV looks valid"

  # Atomically move into place and write markers
  mv -f "$tmp_tsv" "$final_tsv"
  ln -sfn "$final_tsv" "$proj_dir/latest.tsv"  # handy pointer for Step 2
  printf '%s\n' "READY $(ts)" > "$proj_dir/READY"

  local md5; md5="$(md5_of "$final_tsv")"
  local bytes; bytes="$(bytes_of "$final_tsv")"
  append_run_log "$proj_dir" "finalize" "ok" "md5" "$md5" "Saved TSV"
  append_global_manifest "$(ts)" "$project_id" "success" "$bytes" "$rows" "$md5" "$final_tsv" "OK"

  log "Project $project_id: ${rows} rows, ${bytes} bytes, md5=${md5}"
}

# =========================
# Read CSV → unique project IDs
# =========================
read_projects_from_csv(){
  local file="$1"
  if [[ ! -s "$file" ]]; then
    log "ERROR: CSV not found or empty: $file"
    exit 1
  fi

  # Read header and find the column index for project/study accession
  local header; header="$(head -n1 "$file")"
  local idx
  idx=$(csv_find_col_idx "$header" "study_accession" "project_accession" "project_id" "study_id")
  if (( idx == 0 )); then
    log "WARNING: Could not detect project column in header. Falling back to column 2."
    idx=2
  fi

  # Extract column -> dedupe -> strip quotes/spaces
  awk -F',' -v i="$idx" 'NR>1 {gsub(/^[ \t"]+|[ \t"]+$/, "", $i); if($i!="") print $i}' "$file" \
    | sort -u
}

# =========================
# Main
# =========================
main(){
  init_global_manifest
  log "Reading projects from CSV: $CSV_FILE"
  mapfile -t projects < <(read_projects_from_csv "$CSV_FILE")

  if (( ${#projects[@]} == 0 )); then
    log "No projects found in CSV. Nothing to do."
    exit 0
  fi

  log "Found ${#projects[@]} unique project IDs."
  for pid in "${projects[@]}"; do
    log "==> Processing $pid"
    fetch_project_tsv "$pid" || log "Project $pid failed."
  done

  log "Done. Global manifest: $GLOBAL_MANIFEST"
}

main "$@"
