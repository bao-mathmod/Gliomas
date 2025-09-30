#!/usr/bin/env bash
set -Euo pipefail
IFS=$'\n\t'

# ========= Config =========
SPLITS_PARENT="/mnt/12T/chibao/data/official_data/metadata"   # contains split1 ... split9
OUT_ROOT="/mnt/12T/chibao/data/official_data/metadata_enriched"
USER_AGENT="ena-enricher/1.1"
MAX_RETRIES=5
RETRY_SLEEP=2

# Extra ENA fields to fetch and merge (don’t remove join keys)
EXTRA_FIELDS=(
  study_accession
  sample_accession
  experiment_accession
  run_accession
  scientific_name
  fastq_md5
  fastq_ftp
  sra_ftp
  # Add the extra fields you want to fetch
  library_name
  library_strategy
  library_layout
  instrument_platform
  instrument_model
  run_alias
  read_count
  # ...and so on
)

JOIN_KEYS=(study_accession sample_accession experiment_accession run_accession)

# ========= Outputs =========
mkdir -p "$OUT_ROOT"
LOG="$OUT_ROOT/_enrich.log"
MANIFEST="$OUT_ROOT/_manifest_enriched.tsv"
ERRORS="$OUT_ROOT/_errors.tsv"
: > "$LOG"
: > "$ERRORS"
[[ -f "$MANIFEST" ]] || echo -e "ts\tproject_id\tstatus\trows_original\trows_added\trows_merged\torig_tsv\tadd_tsv\tmerged_tsv\tmessage" > "$MANIFEST"

ts(){ date '+%Y-%m-%d %H:%M:%S%z'; }
log(){ printf '[%s] %s\n' "$(ts)" "$*" | tee -a "$LOG" >&2; }

# ========= Helpers =========
url_fields(){
  local IFS=,
  echo "$*"
}
http_fetch(){
  # $1=url $2=outfile -> echoes HTTP code; returns curl exit code
  local url="$1" out="$2"
  local code
  code=$(curl -fL -A "$USER_AGENT" --retry "$MAX_RETRIES" --retry-delay "$RETRY_SLEEP" \
           --retry-connrefused -s -w '%{http_code}' -o "$out" "$url" || true)
  printf '%s' "$code"
}
best_original_tsv(){
  local projdir="$1" proj="$2"
  if [[ -f "$projdir/latest.tsv" ]]; then
    echo "$projdir/latest.tsv"
  elif [[ -f "$projdir/${proj}.tsv" ]]; then
    echo "$projdir/${proj}.tsv"
  else
    find "$projdir" -maxdepth 1 -type f -name '*.tsv' ! -name '_run_log.tsv' -printf '%T@ %p\n' \
      | sort -nr | head -n1 | awk '{print $2}'
  fi
}
rows_excl_header(){ local n=0; n=$(wc -l <"$1" 2>/dev/null || echo 0); echo $((n>0?n-1:0)); }

# Merge two TSVs on composite key (JOIN_KEYS). Keep original cols, add new ones.
merge_on_keys(){
  local orig="$1" extra="$2" out="$3"
  awk -F'\t' -v OFS='\t' '
    function idx(h, name,   i){ for(i=1;i<=length(h);i++) if(h[i]==name) return i; return 0 }
    function key_from(A, h,   a,b,c,d){
      a=A[idx(h,"study_accession")]; b=A[idx(h,"sample_accession")]
      c=A[idx(h,"experiment_accession")]; d=A[idx(h,"run_accession")]
      return a "|" b "|" c "|" d
    }
    FILENAME==ARGV[1] && FNR==1 { split($0,H1); next }
    FILENAME==ARGV[2] && FNR==1 { split($0,H2); next }
    # Store all extra rows by key/column
    FILENAME==ARGV[2] && FNR>1 {
      split("",A)
      for(i=1;i<=NF;i++) A[i]=$i
      K=key_from(A,H2); if(K=="|||") next
      for(i=1;i<=length(H2);i++) E[K SUBSEP H2[i]]=$i
      seenK[K]=1
      next
    }
    # Print merged header: original + any extra-only columns
    FILENAME==ARGV[1] && FNR==1 {
      split("",INMERGED); mN=0
      for(i=1;i<=length(H1);i++){ M[++mN]=H1[i]; INMERGED[H1[i]]=1 }
      for(i=1;i<=length(H2);i++){ if(!(H2[i] in INMERGED)){ M[++mN]=H2[i]; INMERGED[H2[i]]=1 } }
      line=M[1]; for(i=2;i<=mN;i++) line=line OFS M[i]
      print line > "'"$out"'"
      next
    }
    # Merge rows
    FILENAME==ARGV[1] && FNR>1 {
      split("",A)
      for(i=1;i<=NF;i++) A[i]=$i
      K=key_from(A,H1)
      line=""
      for(i=1;i<=mN;i++){
        col=M[i]; val=""
        # Prefer original value
        for(j=1;j<=length(H1);j++) if(H1[j]==col){ val=A[j]; break }
        if(val=="") val=E[K SUBSEP col]
        line = (i==1?val:line OFS val)
      }
      print line >> "'"$out"'"
      next
    }
  ' "$orig" "$extra"
}

concat_extras_same_header(){
  # $1 header_src, $2 dest, $3... parts (each TSV with same header)
  local header_src="$1"; shift
  local dest="$1"; shift
  : > "$dest"
  [[ -s "$header_src" ]] || return 0
  head -n1 "$header_src" > "$dest"
  local h; h="$(head -n1 "$header_src")"
  local f
  for f in "$@"; do
    [[ -s "$f" ]] || continue
    # skip if rows==0 (header only)
    if (( $(rows_excl_header "$f") == 0 )); then
      continue
    fi
    # ensure headers match
    if [[ "$(head -n1 "$f")" != "$h" ]]; then
      echo "[WARN] header mismatch while concatenating extras: $f" >&2
      continue
    fi
    tail -n +2 "$f" >> "$dest"
  done
}

# ========= Build project list from split1..split9 =========
log "Scanning splits under $SPLITS_PARENT"
mapfile -t PROJDIRS < <(find "$SPLITS_PARENT" -maxdepth 2 -type d -regex '.*/split[0-9]+/PRJ.*' | sort)
if (( ${#PROJDIRS[@]} == 0 )); then
  log "No project directories found (expected splitX/PRJNA*)."
  exit 0
fi
log "Found ${#PROJDIRS[@]} project directories."

# ========= Main loop =========
for projdir in "${PROJDIRS[@]}"; do
  proj="$(basename "$projdir")"
  outdir="$OUT_ROOT/$proj"
  mkdir -p "$outdir"

  orig="$(best_original_tsv "$projdir" "$proj")"
  if [[ -z "${orig:-}" || ! -s "$orig" ]]; then
    log "WARN: $proj has no original TSV. Skipping."
    echo -e "$(ts)\t$proj\tskipped\t0\t0\t0\t-\t-\t-\tNo original TSV" >> "$MANIFEST"
    continue
  fi
  rows_orig=$(rows_excl_header "$orig")

  # Determine which extra fields still missing
  declare -A present=()
  while IFS=$'\t' read -r -a cols; do
    for c in "${cols[@]}"; do present["$c"]=1; done
  done < <(head -n1 "$orig")
  to_fetch=()
  for f in "${EXTRA_FIELDS[@]}"; do
    [[ -z "${present[$f]+x}" ]] && to_fetch+=( "$f" )
  done

  if (( ${#to_fetch[@]} == 0 )); then
    merged="$outdir/${proj}_metadata.tsv"
    cp -f "$orig" "$merged"
    ln -sfn "$merged" "$outdir/latest.tsv"
    echo -e "$(ts)\t$proj\tsuccess\t$rows_orig\t0\t$rows_orig\t$orig\t-\t$merged\tNo extra fields needed" >> "$MANIFEST"
    log "INFO: $proj already has all requested fields; copied original."
    continue
  fi

  # Build requests
  fields_param="$(url_fields "${JOIN_KEYS[@]}" "${to_fetch[@]}")"
  url_proj="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${proj}&result=read_run&fields=${fields_param}&format=tsv&download=true&limit=0"
  add_proj="$outdir/${proj}.extra.by_project.tsv"
  log "FETCH[project]: $proj"
  code=$(http_fetch "$url_proj" "$add_proj")
  log "HTTP $code for $proj (project-level)"
  rows_added_proj=0
  if [[ -s "$add_proj" ]]; then
    rows_added_proj=$(rows_excl_header "$add_proj")
  fi

  # Fallback A: by each study_accession found in original
  studies=()
  mapfile -t studies < <(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="study_accession") {c=i; break}} NR>1 && c>0 && $c!=""{print $c}' "$orig" | sort -u)
  add_by_study_parts=()
  if (( rows_added_proj == 0 )) && (( ${#studies[@]} )); then
    log "FALLBACK[studies]: $proj studies=(${studies[*]})"
    for srp in "${studies[@]}"; do
      url_study="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${srp}&result=read_run&fields=${fields_param}&format=tsv&download=true&limit=0"
      part="$outdir/${proj}.extra.by_study.${srp}.tsv"
      code_s=$(http_fetch "$url_study" "$part")
      log "HTTP $code_s for $proj study=$srp"
      if [[ -s "$part" && $(rows_excl_header "$part") -gt 0 ]]; then
        add_by_study_parts+=( "$part" )
      else
        rm -f "$part" || true
      fi
    done
  fi

  # Fallback B: search endpoint by project_accession
  add_search="$outdir/${proj}.extra.by_search.tsv"
  if (( rows_added_proj == 0 )) && (( ${#add_by_study_parts[@]} == 0 )); then
    query=$(python3 - <<'PY'
import urllib.parse
print(urllib.parse.quote('project_accession="' + """'" + '""" + '"'))
PY
) # harmless placeholder – we’ll just hardcode below to avoid py dep
    # No python: do it manually
    url_search="https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=project_accession=%22${proj}%22&fields=${fields_param}&format=tsv&limit=0"
    log "FALLBACK[search]: $proj"
    code_srch=$(http_fetch "$url_search" "$add_search")
    log "HTTP $code_srch for $proj (search)"
    [[ -s "$add_search" && $(rows_excl_header "$add_search") -gt 0 ]] || rm -f "$add_search" || true
  fi

  # Consolidate extras (prefer project-level; else concat study parts; else search)
  add_all="$outdir/${proj}.extra.tsv"
  : > "$add_all"
  rows_added=0
  if (( rows_added_proj > 0 )); then
    # project-level has rows: use it directly
    cp -f "$add_proj" "$add_all"
    rows_added=$rows_added_proj
  elif (( ${#add_by_study_parts[@]} )); then
    # concat all study parts (same header)
    # use the first available header as template
    first_part="${add_by_study_parts[0]}"
    concat_extras_same_header "$first_part" "$add_all" "${add_by_study_parts[@]}"
    rows_added=$(rows_excl_header "$add_all")
  elif [[ -s "$add_search" ]]; then
    cp -f "$add_search" "$add_all"
    rows_added=$(rows_excl_header "$add_all")
  else
    # truly nothing fetched
    :
  fi

  # If still nothing, just copy original; otherwise merge.
  if (( rows_added == 0 )); then
    merged="$outdir/${proj}_metadata.tsv"
    cp -f "$orig" "$merged"
    ln -sfn "$merged" "$outdir/latest.tsv"
    echo -e "$(ts)\t$proj\tno_extra\t$rows_orig\t0\t$rows_orig\t$orig\t$add_all\t$merged\tNo extra rows from ENA (project/study/search all empty)" >> "$MANIFEST"
    echo -e "$proj\tEMPTY_EXTRA\tproject+study+search returned no rows" >> "$ERRORS"
    log "WARN: $proj extra rows=0; kept original only."
    continue
  fi

  merged="$outdir/${proj}_metadata.tsv"
  log "MERGE: $proj  orig_rows=$rows_orig  add_rows=$rows_added"
  merge_on_keys "$orig" "$add_all" "$merged" || {
    echo -e "$(ts)\t$proj\tmerge_failed\t$rows_orig\t$rows_added\t0\t$orig\t$add_all\t$merged\tMerge error" >> "$MANIFEST"
    echo -e "$proj\tMERGE_FAILED\t$orig\t$add_all" >> "$ERRORS"
    log "ERROR: merge failed for $proj"
    continue
  }
  rows_merged=$(rows_excl_header "$merged")
  ln -sfn "$merged" "$outdir/latest.tsv"
  echo -e "$(ts)\t$proj\tsuccess\t$rows_orig\t$rows_added\t$rows_merged\t$orig\t$add_all\t$merged\tOK" >> "$MANIFEST"
  log "DONE: $proj merged_rows=$rows_merged"
done

log "=== ENRICH COMPLETE ==="
log "Manifest: $MANIFEST"
log "Errors:   $ERRORS"
