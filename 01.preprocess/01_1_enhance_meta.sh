#######################################
# Section 1: Fetch ENA metadata TSVs
#######################################

# #!/usr/bin/env bash
# set -Eeuo pipefail
# IFS=$'\n\t'

# # =========================
# # Config (edit as needed)
# # =========================
# CSV_FILE="/mnt/12T/chibao/data/official_data/glioma_accession.csv"
# META_ROOT="/mnt/12T/chibao/data/stuff_data/special/official"
# GLOBAL_MANIFEST="$META_ROOT/_manifest_projects.tsv"
# USER_AGENT="ena-meta-fetch/1.0"
# MAX_RETRIES=5
# RETRY_DELAY=2   # seconds

# # ENA fields (expanded)
# REQUIRED_FIELDS=(
#   study_accession
#   sample_accession
#   experiment_accession
#   run_accession
#   scientific_name
#   fastq_md5
#   fastq_ftp
#   sra_ftp
#   library_name
#   library_strategy
#   library_layout
#   instrument_platform
#   instrument_model
#   run_alias
#   read_count
# )

# # =========================
# # Logging helpers
# # =========================
# ts(){ date '+%Y-%m-%d %H:%M:%S%z'; }
# log(){ printf '[%s] %s\n' "$(ts)" "$*" >&2; }

# init_global_manifest(){
#   mkdir -p "$META_ROOT"
#   if [[ ! -f "$GLOBAL_MANIFEST" ]]; then
#     echo -e "ts\tproject_id\tstatus\tbytes\trows\tmd5\ttsv_path\tmessage" > "$GLOBAL_MANIFEST"
#   fi
# }

# append_global_manifest(){
#   local ts_="$1" proj="$2" status="$3" bytes="$4" rows="$5" md5="$6" path="$7" msg="$8"
#   echo -e "${ts_}\t${proj}\t${status}\t${bytes}\t${rows}\t${md5}\t${path}\t${msg}" >> "$GLOBAL_MANIFEST"
# }

# init_run_log(){
#   local proj_dir="$1"
#   local logf="$proj_dir/_run_log.tsv"
#   if [[ ! -f "$logf" ]]; then
#     echo -e "ts\tstage\tstatus\tfield\tvalue\tmessage" > "$logf"
#   fi
# }

# append_run_log(){
#   local proj_dir="$1" stage="$2" status="$3" field="$4" value="$5" msg="$6"
#   echo -e "$(ts)\t${stage}\t${status}\t${field}\t${value}\t${msg}" >> "$proj_dir/_run_log.tsv"
# }

# # =========================
# # Utilities
# # =========================
# trim(){ awk '{$1=$1;print}' <<<"$*"; }

# csv_find_col_idx(){
#   local header="$1"; shift
#   local candidates=("$@")
#   local IFS=,
#   read -r -a cols <<<"$header"
#   for name in "${candidates[@]}"; do
#     for i in "${!cols[@]}"; do
#       if [[ "$(trim "${cols[$i]}")" == "$name" ]]; then
#         echo $((i+1)); return 0
#       fi
#     done
#   done
#   echo 0
# }

# bytes_of(){ [[ -f "$1" ]] && { stat -c '%s' "$1" 2>/dev/null || wc -c <"$1"; } || echo 0; }
# md5_of(){ [[ -f "$1" ]] && md5sum "$1" | awk '{print $1}' || echo ""; }

# tsv_has_required_fields(){
#   local header="$1"; shift
#   local missing=()
#   for f in "$@"; do
#     if ! awk -F'\t' -v want="$f" 'NR==1{for(i=1;i<=NF;i++) if($i==want) found=1} END{exit(found?0:1)}' <<<"$header"; then
#       missing+=("$f")
#     fi
#   done
#   if (( ${#missing[@]} )); then
#     printf '%s' "${missing[*]}"; return 1
#   fi
#   return 0
# }

# # =========================
# # ENA downloader (with retry)
# # =========================
# fetch_project_tsv(){
#   local project_id="$1"
#   local proj_dir="$META_ROOT/$project_id"
#   local final_tsv="$proj_dir/${project_id}.tsv"
#   local tmp_tsv="$proj_dir/${project_id}.tsv.part"
#   local url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${project_id}&result=read_run&fields=$(IFS=,; echo "${REQUIRED_FIELDS[*]}")&format=tsv&download=true&limit=0"

#   mkdir -p "$proj_dir"
#   init_run_log "$proj_dir"
#   append_run_log "$proj_dir" "prepare" "info" "url" "$url" "Built ENA request URL"

#   local attempt=1
#   while (( attempt <= MAX_RETRIES )); do
#     append_run_log "$proj_dir" "download" "start" "attempt" "$attempt" "Fetching TSV"
#     if curl -fL -A "$USER_AGENT" --retry "$MAX_RETRIES" --retry-delay "$RETRY_DELAY" \
#            --retry-connrefused -s "$url" -o "$tmp_tsv"; then
#       append_run_log "$proj_dir" "download" "ok" "attempt" "$attempt" "Downloaded"; break
#     else
#       append_run_log "$proj_dir" "download" "error" "attempt" "$attempt" "curl failed; will retry"
#       sleep "$RETRY_DELAY"
#     fi
#     ((attempt++))
#   done

#   if [[ ! -s "$tmp_tsv" ]]; then
#     append_run_log "$proj_dir" "download" "failed" "file" "$tmp_tsv" "Empty or missing after retries"
#     append_global_manifest "$(ts)" "$project_id" "failed" "0" "0" "" "" "Download failed"
#     return 1
#   fi

#   local header; header="$(head -n1 "$tmp_tsv")"
#   if [[ -z "$header" ]]; then
#     append_run_log "$proj_dir" "validate" "error" "header" "-" "Missing header"
#     append_global_manifest "$(ts)" "$project_id" "failed" "$(bytes_of "$tmp_tsv")" "0" "" "$tmp_tsv" "Missing header"
#     return 1
#   fi

#   if ! tsv_has_required_fields "$header" "${REQUIRED_FIELDS[@]}"; then
#     local missing; missing="$(tsv_has_required_fields "$header" "${REQUIRED_FIELDS[@]}"; true || true)"
#     append_run_log "$proj_dir" "validate" "error" "missing_fields" "$missing" "Header missing required fields"
#     append_global_manifest "$(ts)" "$project_id" "failed" "$(bytes_of "$tmp_tsv")" "0" "" "$tmp_tsv" "Missing fields: $missing"
#     return 1
#   fi

#   local rows; rows="$(($(wc -l <"$tmp_tsv") - 1))"; (( rows<0 )) && rows=0
#   append_run_log "$proj_dir" "validate" "ok" "rows" "$rows" "TSV looks valid"

#   mv -f "$tmp_tsv" "$final_tsv"
#   ln -sfn "$final_tsv" "$proj_dir/latest.tsv"   # (symlink is fine for the *source* tree)
#   printf '%s\n' "READY $(ts)" > "$proj_dir/READY"

#   local md5; md5="$(md5_of "$final_tsv")"
#   local bytes; bytes="$(bytes_of "$final_tsv")"
#   append_run_log "$proj_dir" "finalize" "ok" "md5" "$md5" "Saved TSV"
#   append_global_manifest "$(ts)" "$project_id" "success" "$bytes" "$rows" "$md5" "$final_tsv" "OK"
#   log "Project $project_id: rows=${rows} bytes=${bytes} md5=${md5}"
# }

# # =========================
# # Read CSV → unique project IDs
# # =========================
# read_projects_from_csv(){
#   local file="$1"
#   [[ -s "$file" ]] || { log "ERROR: CSV not found or empty: $file"; exit 1; }
#   local header; header="$(head -n1 "$file")"
#   local idx; idx=$(csv_find_col_idx "$header" "study_accession" "project_accession" "project_id" "study_id")
#   (( idx==0 )) && { log "WARNING: project column not found; fallback to column 2"; idx=2; }
#   awk -F',' -v i="$idx" 'NR>1 {gsub(/^[ \t"]+|[ \t"]+$/, "", $i); if($i!="") print $i}' "$file" | sort -u
# }

# # =========================
# # Main
# # =========================
# main(){
#   init_global_manifest
#   log "Reading projects from CSV: $CSV_FILE"
#   mapfile -t projects < <(read_projects_from_csv "$CSV_FILE")
#   (( ${#projects[@]}==0 )) && { log "No projects in CSV."; exit 0; }

#   log "Found ${#projects[@]} unique project IDs."
#   for pid in "${projects[@]}"; do
#     log "==> $pid"
#     fetch_project_tsv "$pid" || log "Project $pid failed."
#   done
#   log "Done. Global manifest: $GLOBAL_MANIFEST"
# }
# main "$@"

###################################
# After done step 1: Base on the file /mnt/12T/chibao/data/official_data/metadata/split_map.tsv to organize the split before step 2
###################################


#####################################################
# Section 2: Enhance metadata and saved to official/
#####################################################
#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

SRC_META_ROOT="/mnt/12T/chibao/data/stuff_data/special/official"
OFFICIAL_META_ROOT="/mnt/12T/chibao/data/official_data/metadata/official"
FASTQ_ROOT="/mnt/12T/chibao/data/official_data/fastq"
USER_AGENT="ena-meta-enhancer/1.0"

# Base columns (must exist)
BASE_COLS=(
  study_accession
  sample_accession
  experiment_accession
  run_accession
  scientific_name
  fastq_md5
  fastq_ftp
  sra_ftp
)

# Optional columns (fill empty if absent)
OPT_COLS=(
  library_name
  library_strategy
  library_layout
  instrument_platform
  instrument_model
  run_alias
  read_count
)

ts(){ date '+%Y-%m-%d %H:%M:%S%z'; }
log(){ printf '[%s] %s\n' "$(ts)" "$*" >&2; }
ensure_dir(){ mkdir -p "$1"; }
bytes_of(){ [[ -f "$1" ]] && { stat -c '%s' "$1" 2>/dev/null || wc -c <"$1"; } || echo 0; }
md5_of(){ [[ -f "$1" ]] && md5sum "$1" | awk '{print $1}' || echo ""; }

validate_header(){
  local header="$1"; shift
  local missing=()
  for f in "$@"; do
    if ! awk -F'\t' -v want="$f" 'NR==1{for(i=1;i<=NF;i++) if($i==want) ok=1} END{exit(ok?0:1)}' <<<"$header"; then
      missing+=("$f")
    fi
  done
  if (( ${#missing[@]} )); then printf '%s' "${missing[*]}"; return 1; fi
  return 0
}

split_from_path(){ basename "$(dirname "$1")"; }  # parent folder = split id

enhance_project(){
  local split_id="$1" project_id="$2" src_tsv="$3" dest_dir="$4"

  ensure_dir "$dest_dir"
  local runlog="$dest_dir/_run_log.tsv"
  [[ -f "$runlog" ]] || echo -e "ts\tstage\tstatus\tfield\tvalue\tmessage" > "$runlog"

  [[ -s "$src_tsv" ]] || { echo -e "$(ts)\tinit\terror\tsrc_tsv\t$src_tsv\tMissing"; return 1; }

  local header; header="$(head -n1 "$src_tsv")"
  if ! validate_header "$header" "${BASE_COLS[@]}"; then
    local miss; miss="$(validate_header "$header" "${BASE_COLS[@]}"; true || true)"
    echo -e "$(ts)\tvalidate\terror\tmissing_base\t$miss\tRequired columns missing" >> "$runlog"
    return 1
  fi

  # Build indices for all present columns
  # Any optional column not present will be emitted as empty string
  local latest_tsv="$dest_dir/latest.tsv"
  local tmp="$dest_dir/latest.tsv.part"

  awk -F'\t' -v OFS='\t' \
      -v split_id="$split_id" \
      -v fastq_root="$FASTQ_ROOT" \
      -v tmp_out="$tmp" '
    function idx(name,    i){
      # safe lookup: 0 if absent
      return (name in H ? H[name] : 0)
    }
    BEGIN{ }
    NR==1{
      for(i=1;i<=NF;i++) H[$i]=i
      # header
      print \
        "study_accession","sample_accession","experiment_accession","run_accession","scientific_name", \
        "fastq_ftp","fastq_md5","sra_ftp","library_name","library_strategy","library_layout", \
        "instrument_platform","instrument_model","run_alias","read_count", \
        "fastq_url_https","fastq_basename","md5_expected","split_id","target_path" > tmp_out
      next
    }
    NR>1{
      # required
      study = $(H["study_accession"])
      sample= $(H["sample_accession"])
      exper = $(H["experiment_accession"])
      run   = $(H["run_accession"])
      sci   = $(H["scientific_name"])
      fftp  = $(H["fastq_ftp"])
      fmd5  = $(H["fastq_md5"])
      sra   = $(H["sra_ftp"])

      # optional (guard missing -> "")
      lname = (idx("library_name") ? $(H["library_name"]) : "")
      lstr  = (idx("library_strategy") ? $(H["library_strategy"]) : "")
      llay  = (idx("library_layout") ? $(H["library_layout"]) : "")
      iplat = (idx("instrument_platform") ? $(H["instrument_platform"]) : "")
      imod  = (idx("instrument_model") ? $(H["instrument_model"]) : "")
      ralias= (idx("run_alias") ? $(H["run_alias"]) : "")
      rcnt  = (idx("read_count") ? $(H["read_count"]) : "")

      n1=split(fftp, U, /;/); n2=split(fmd5, M, /;/)
      if(n2<n1){ for(i=n2+1;i<=n1;i++) M[i]="" }

      for(i=1;i<=n1;i++){
        url_raw=U[i]; gsub(/^[ \t]+|[ \t]+$/, "", url_raw); if(url_raw=="") continue
        url_https=url_raw; sub(/^ftp:\/\//,"https://",url_https)
        n=split(url_https, P, "/"); fname=P[n]
        md5=M[i]; gsub(/^[ \t]+|[ \t]+$/, "", md5)

        target = fastq_root "/" split_id "/" study "/" sample "/" exper "/" run "/" fname

        print study,sample,exper,run,sci, \
              url_raw,md5,sra,lname,lstr,llay,iplat,imod,ralias,rcnt, \
              url_https,fname,md5,split_id,target >> tmp_out
      }
    }
  ' "$src_tsv"

  mv -f "$tmp" "$latest_tsv"
  echo -e "$(ts)\tfinalize\tok\tlatest\t$latest_tsv\tEnhanced" >> "$runlog"

  # Append to per-split _expected_files.tsv
  local split_root="$OFFICIAL_META_ROOT/$split_id"
  ensure_dir "$split_root"
  local expected="$split_root/_expected_files.tsv"
  if [[ ! -f "$expected" ]]; then
    echo -e "split_id\tproject_id\tstudy\tsample\texperiment\trun\tfile\turl\tmd5_expected\ttarget_path" > "$expected"
  fi
  awk -F'\t' -v OFS='\t' -v sid="$split_id" -v proj="$project_id" \
      'NR>1{print sid,proj,$1,$2,$3,$4,$17,$16,$18,$20}' "$latest_tsv" >> "$expected"
}

main(){
  log "Source:   $SRC_META_ROOT"
  log "Official: $OFFICIAL_META_ROOT"

  shopt -s nullglob
  mapfile -t project_dirs < <(find "$SRC_META_ROOT" -maxdepth 2 -mindepth 2 -type d -regex ".*/split[^/]+/PRJNA[0-9]+" | sort)
  shopt -u nullglob
  (( ${#project_dirs[@]}==0 )) && { log "No split*/PRJNA* found"; exit 0; }

  for src_dir in "${project_dirs[@]}"; do
    project_id="$(basename "$src_dir")"
    split_id="$(split_from_path "$src_dir")"

    # prefer real TSV inside split folder; else resolve symlink if present
    src_tsv=""
    if [[ -f "$src_dir/${project_id}.tsv" ]]; then
      src_tsv="$src_dir/${project_id}.tsv"
    elif [[ -L "$src_dir/latest.tsv" ]]; then
      target="$(readlink -f "$src_dir/latest.tsv" 2>/dev/null || true)"
      [[ -f "$target" ]] && src_tsv="$target"
    fi
    if [[ -z "$src_tsv" || ! -s "$src_tsv" ]]; then
      log "WARN: Missing TSV for $project_id under $src_dir; skipping."
      continue
    fi

    dest_dir="$OFFICIAL_META_ROOT/$split_id/$project_id"
    ensure_dir "$dest_dir"
    enhance_project "$split_id" "$project_id" "$src_tsv" "$dest_dir" || log "Enhance failed: $split_id/$project_id"
  done

  log "Done. Enhanced metadata under: $OFFICIAL_META_ROOT"
}
main "$@"
