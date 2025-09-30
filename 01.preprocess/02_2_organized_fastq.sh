#!/usr/bin/env bash
# ============================================================
# 04_organize_fastq_by_tech.sh
# Organize FASTQ files into technique-based folders using your
# enriched ENA metadata and the download index.
#
# USAGE:
#   bash 04_organize_fastq_by_tech.sh \
#     [--fastq-root /mnt/12T/chibao/data/official_data/fastq/split4] \
#     [--expected-index /mnt/12T/chibao/data/official_data/fastq/split4/_expected_files.tsv] \
#     [--enriched-root /mnt/12T/chibao/data/official_data/metadata/official] \
#     [--out-root /mnt/12T/chibao/data/official_data/fastq_by_tech] \
#     [--dry-run]
#
# You can safely re-run this script. It appends to TSV logs and
# only creates/updates missing symlinks. Nothing in your source
# folders is moved or deleted.
#
# OUTPUTS (under --out-root):
#   _organize.log                  # timestamped log
#   _tech_assign.tsv               # per-run classification
#   _symlink_manifest.tsv          # symlink actions taken
#   _pending_projects.tsv          # projects pending (no metadata or files missing)
#   _summary_by_tech.tsv           # counts by tech bucket
#
# Folder layout created (symlinked view):
#   <out-root>/
#     scRNA_GEX/<PRJNA>/... (mirrors source path)
#     scATAC/<PRJNA>/...
#     multiome_ARC/<PRJNA>/...
#     bulk_or_plate_RNA/<PRJNA>/...
#     unknown_no_metadata/<PRJNA>/...
#     unknown_unclassified/<PRJNA>/...
#
# Requirements: bash, awk, grep, sed, sort, ln, mkdir, readlink (usually present)
# ============================================================

# This script will organize FASTQ files into technique-based folders - Safe to re-run
set -Euo pipefail
IFS=$'\n\t'

# ---------- Defaults (override with CLI flags) ----------
FASTQ_ROOT="/mnt/12T/chibao/data/official_data/fastq/split5"
: "${SPLIT_ID:=$(basename "$FASTQ_ROOT")}"   # e.g., split4 if FASTQ_ROOT ends with .../fastq/split4
EXPECTED_INDEX="/mnt/12T/chibao/data/official_data/metadata/official/${SPLIT_ID}/_expected_files.tsv"
ENRICHED_ROOT="/mnt/12T/chibao/data/official_data/metadata/official"
OUT_ROOT="/mnt/12T/chibao/data/official_data/fastq_by_tech"
DRY_RUN="false"
QUIET_UNCHANGED="true"           # don’t log 'exists_ok' rows on re-runs
PRUNE_STALE="true"               # remove old symlink if bucket changes
STATE="$OUT_ROOT/_state_class.tsv"   # current classification state
CLEAN_DANGLING="true"   # remove existing symlink if source file is missing

# Projects known to have no extra metadata (keep in 'unknown_no_metadata')
NO_META_PROJECTS=("PRJNA577146" "PRJNA870065" "PRJNA941288" "PRJNA1009671" "PRJNA1134155")

# ---------- CLI ----------
while (( $# )); do
  case "$1" in
    --fastq-root)     FASTQ_ROOT="${2:-}"; shift 2;;
    --expected-index) EXPECTED_INDEX="${2:-}"; shift 2;;
    --enriched-root)  ENRICHED_ROOT="${2:-}"; shift 2;;
    --out-root)       OUT_ROOT="${2:-}"; shift 2;;
    --dry-run)        DRY_RUN="true"; shift 1;;
    *)                echo "WARN: unknown arg $1" >&2; shift 1;;
  esac
done


# ---------- Paths & logs ----------
mkdir -p "$OUT_ROOT"
LOG="$OUT_ROOT/_organize.log"
TECH_ASSIGN="$OUT_ROOT/_tech_assign.tsv"
SYMLINK_MANIFEST="$OUT_ROOT/_symlink_manifest.tsv"
PENDING_PROJ="$OUT_ROOT/_pending_projects.tsv"
SUMMARY_BY_TECH="$OUT_ROOT/_summary_by_tech.tsv"

ts(){ date '+%Y-%m-%d %H:%M:%S%z'; }
log(){ printf '[%s] %s\n' "$(ts)" "$*" | tee -a "$LOG" >&2; }

# Headers (append-safe)
[[ -f "$TECH_ASSIGN" ]] || echo -e "ts\tproject\tsample\trun\ttech_bucket\tsource\treason\tclass_fields" > "$TECH_ASSIGN"
[[ -f "$SYMLINK_MANIFEST" ]] || echo -e "ts\tproject\ttech_bucket\tsrc_path\tdst_path\taction\tnote" > "$SYMLINK_MANIFEST"
[[ -f "$PENDING_PROJ" ]] || echo -e "ts\tproject\treason\tnote" > "$PENDING_PROJ"

# ---------- Sanity ----------
[[ -s "$EXPECTED_INDEX" ]] || { log "ERROR: missing EXPECTED_INDEX: $EXPECTED_INDEX"; exit 1; }
[[ -d "$ENRICHED_ROOT"  ]] || { log "ERROR: enriched root not found: $ENRICHED_ROOT"; exit 1; }

# State file: project  run  last_bucket  last_dst
if [[ ! -f "$STATE" ]]; then
  echo -e "project\trun\tlast_bucket\tlast_dst" > "$STATE"
fi
declare -A LAST_BUCKET LAST_DST
while IFS=$'\t' read -r P R B D; do
  [[ "$P" == "project" ]] && continue
  LAST_BUCKET["$P|$R"]="$B"
  LAST_DST["$P|$R"]="$D"
done < "$STATE"

# ---------- Helpers ----------
in_array(){ # needle list...
  local n="$1"; shift
  for x in "$@"; do [[ "$x" == "$n" ]] && return 0; done
  return 1
}

# Prefer the extra.tsv (what you asked for). Fallback to merged; else none.
project_meta_source(){
  local proj="$1"
  local dir="$ENRICHED_ROOT/${SPLIT_ID}/${proj}"
  local latest="$dir/latest.tsv"
  if [[ -s "$latest" && $(wc -l <"$latest") -gt 1 ]]; then
    echo "$latest"; return 0
  fi
  echo ""   # no metadata for this project
}


# Build per-project run->tech map using enriched metadata (if present)
build_run2tech(){
  # Args: project_id
  # Emits: run_accession \t sample_accession \t tech_bucket \t source \t reason \t class_fields
  local proj="$1"
  local meta_tsv; meta_tsv="$(project_meta_source "$proj")"

  local no_meta="false"
  if [[ -z "$meta_tsv" || ! -s "$meta_tsv" ]]; then
    no_meta="true"
  fi
  if in_array "$proj" "${NO_META_PROJECTS[@]}"; then
    no_meta="true"
  fi
  [[ "$no_meta" == "true" ]] && return 0

  awk -F'\t' -v OFS='\t' '
    NR==1{
      for(i=1;i<=NF;i++) h[$i]=i
      rk=h["run_accession"]; sk=h["sample_accession"]; lk=h["library_strategy"]
      if(!rk || !lk){ exit 0 }   # need run + library_strategy to classify
      next
    }
    {
      run = (rk? $rk : "")
      samp = (sk? $sk : "")
      ls = (lk? $lk : "")
      l = tolower(ls)

      tech="unknown_unclassified"; reason="library_strategy=" ls

      if (l ~ /atac/) {
        tech="scATAC"
      } else if (l ~ /(single.?cell.*rna|scrna|snrna)/) {
        tech="scRNA_GEX"
      } else if (l ~ /rna/) {
        tech="bulk_or_plate_RNA"
      }

      print run, samp, tech, "library_strategy", reason, "ls=" ls
    }
  ' "$meta_tsv"
}

# Decide tech bucket when we DON’T have metadata (or run not found in metadata)
fallback_bucket_no_meta(){
  # Conservative: we don’t guess. Use "unknown_no_metadata".
  echo "unknown_no_metadata"
}

# Compute destination path under OUT_ROOT
dest_for(){
  # Args: tech_bucket, project_id, src_target_path
  local tech="$1" proj="$2" src="$3"
  local rel="$src"

  # Strip the FASTQ root
  if [[ "$rel" == "$FASTQ_ROOT/"* ]]; then
    rel="${rel#$FASTQ_ROOT/}"      # now starts: PRJNA.../SAMN.../SRX.../SRR.../file.gz
  fi

  # If rel already starts with this project, don't repeat it
  if [[ "$rel" == "$proj/"* ]]; then
    echo "$OUT_ROOT/$tech/$rel"
  else
    echo "$OUT_ROOT/$tech/$proj/$rel"
  fi
}

# Ensure directory, create/update symlink
make_link(){
  local proj="$1" tech="$2" src="$3" dst="$4"
  local action note

  # If the source file doesn't exist, optionally clean up a dangling link
  if [[ ! -e "$src" ]]; then
    if [[ -L "$dst" && "$CLEAN_DANGLING" == "true" ]]; then
      if [[ "$DRY_RUN" == "true" ]]; then
        action="would_remove_dangling"; note="src_missing"
      else
        rm -f "$dst"
        action="removed_dangling"; note="src_missing"
      fi
      echo -e "$(ts)\t$proj\t$tech\t$src\t$dst\t$action\t$note" >> "$SYMLINK_MANIFEST"
      return 0
    fi
    # No src and no dangling link to clean — just record the miss
    action="skip_missing_src"; note="source_not_found"
    echo -e "$(ts)\t$proj\t$tech\t$src\t$dst\t$action\t$note" >> "$SYMLINK_MANIFEST"
    return 0
  fi

  # Source exists -> ensure destination directory and link
  mkdir -p "$(dirname "$dst")"
  if [[ -L "$dst" ]]; then
    local cur; cur="$(readlink "$dst" || true)"
    if [[ "$cur" == "$src" ]]; then
      action="exists_ok"; note="-"
    else
      if [[ "$DRY_RUN" == "true" ]]; then
        action="would_update_link"; note="from:$cur"
      else
        rm -f "$dst"
        ln -s "$src" "$dst"
        action="updated_link"; note="from:$cur"
      fi
    fi
  elif [[ -e "$dst" ]]; then
    action="conflict_exists"; note="not_symlink"
  else
    if [[ "$DRY_RUN" == "true" ]]; then
      action="would_create"; note="-"
    else
      ln -s "$src" "$dst"
      action="created"; note="-"
    fi
  fi

  # Quiet unchanged logs if desired
  if [[ "$action" == "exists_ok" && "$QUIET_UNCHANGED" == "true" ]]; then
    return 0
  fi
  echo -e "$(ts)\t$proj\t$tech\t$src\t$dst\t$action\t$note" >> "$SYMLINK_MANIFEST"
}

save_state(){
  local tmp="$STATE.tmp.$$"
  echo -e "project\trun\tlast_bucket\tlast_dst" > "$tmp"
  for key in "${!LAST_BUCKET[@]}"; do
    P="${key%%|*}"; R="${key#*|}"
    echo -e "$P\t$R\t${LAST_BUCKET[$key]:-}\t${LAST_DST[$key]:-}" >> "$tmp"
  done
  mv -f "$tmp" "$STATE"
}

log "=== START organize ==="
log "FASTQ_ROOT=$FASTQ_ROOT"
log "EXPECTED_INDEX=$EXPECTED_INDEX"
log "ENRICHED_ROOT=$ENRICHED_ROOT"
log "OUT_ROOT=$OUT_ROOT"
log "DRY_RUN=$DRY_RUN"

# ---------- Load expected index rows once ----------
# Columns: 1 split_id 2 project 3 study 4 sample 5 experiment 6 run 7 file 8 url 9 md5 10 target_path
mapfile -t EXPECTED_ROWS < <(awk -F'\t' 'NR>1 && $10!="" {print $2"\t"$4"\t"$6"\t"$10}' "$EXPECTED_INDEX")
# Build project list from expected index
mapfile -t PROJECTS < <(printf '%s\n' "${EXPECTED_ROWS[@]}" | awk -F'\t' '{print $1}' | sort -u)

if (( ${#PROJECTS[@]} == 0 )); then
  log "No projects found in EXPECTED_INDEX (or target paths empty). Nothing to organize."
  exit 0
fi

# ---------- Main loop over projects ----------
for proj in "${PROJECTS[@]}"; do
# Build run->tech map from enriched metadata (if present)
declare -A RUN2TECH=()
declare -A RUN2SAMPLE=()
declare -A RUN2WHY=()
declare -A RUN2SOURCE=()
declare -A RUN2FIELDS=()

  while IFS=$'\t' read -r run_acc samp tech src reason fields; do
    [[ -z "${run_acc:-}" ]] && continue
    RUN2TECH["$run_acc"]="$tech"
    RUN2SAMPLE["$run_acc"]="$samp"
    RUN2SOURCE["$run_acc"]="$src"
    RUN2WHY["$run_acc"]="$reason"
    RUN2FIELDS["$run_acc"]="$fields"
  done < <(build_run2tech "$proj")

  # For multiome ARC: if same sample has scATAC AND scRNA_GEX runs, re-label both as multiome_ARC
  if (( ${#RUN2TECH[@]} )); then
    declare -A SAMP_HAS_ATAC SAMP_HAS_GEX
    for r in "${!RUN2TECH[@]}"; do
      s="${RUN2SAMPLE[$r]}"
      t="${RUN2TECH[$r]}"
      [[ "$t" == "scATAC" ]] && SAMP_HAS_ATAC["$s"]=1
      [[ "$t" == "scRNA_GEX" ]] && SAMP_HAS_GEX["$s"]=1
    done
    for r in "${!RUN2TECH[@]}"; do
      s="${RUN2SAMPLE[$r]}"
      if [[ -n "${SAMP_HAS_ATAC[$s]:-}" && -n "${SAMP_HAS_GEX[$s]:-}" ]]; then
        RUN2TECH["$r"]="multiome_ARC"
        RUN2WHY["$r"]="same_sample_has_ATAC_and_GEX"
      fi
    done
  fi

  # When enriched metadata is missing entirely
  local_meta_present="yes"
  if (( ${#RUN2TECH[@]} == 0 )); then
    local_meta_present="no"
    if ! in_array "$proj" "${NO_META_PROJECTS[@]}"; then
      echo -e "$(ts)\t$proj\tno_enriched_metadata\twill_use_unknown_bucket" >> "$PENDING_PROJ"
    fi
  fi

  # Walk expected rows for this project and symlink by bucket
  while IFS=$'\t' read -r p s r src; do
    [[ "$p" != "$proj" ]] && continue
    tech_bucket=""
    source="enriched_metadata"
    reason=""
    fields=""

    if [[ -n "${RUN2TECH[$r]:-}" ]]; then
      tech_bucket="${RUN2TECH[$r]}"
      reason="${RUN2WHY[$r]:-}"
      fields="${RUN2FIELDS[$r]:-}"
    else
      if [[ "$local_meta_present" == "no" ]]; then
        tech_bucket="$(fallback_bucket_no_meta)"
        source="no_metadata"
        reason="project_has_no_enriched_metadata"
      else
        tech_bucket="unknown_unclassified"
        source="partial_metadata"
        reason="run_not_found_in_enriched"
      fi
    fi

    # Emit classification (per run)
    echo -e "$(ts)\t$proj\t$s\t$r\t$tech_bucket\t$source\t$reason\t$fields" >> "$TECH_ASSIGN"

    # Destination path & symlink
    dst="$(dest_for "$tech_bucket" "$proj" "$src")"

    # If bucket changed since last run, optionally remove old symlink
    state_key="$proj|$r"
    prev_bucket="${LAST_BUCKET[$state_key]:-}"
    prev_dst="${LAST_DST[$state_key]:-}"
if [[ "$PRUNE_STALE" == "true" && -n "$prev_bucket" && "$prev_bucket" != "$tech_bucket" && -n "$prev_dst" ]]; then
  if [[ -L "$prev_dst" && "$(readlink "$prev_dst" || true)" == "$src" ]]; then
    if [[ "$DRY_RUN" == "true" ]]; then
      echo -e "$(ts)\t$proj\t$prev_bucket\t$src\t$prev_dst\twould_remove_stale\tbucket_change:${prev_bucket}->${tech_bucket}" >> "$SYMLINK_MANIFEST"
    else
      rm -f "$prev_dst"
      echo -e "$(ts)\t$proj\t$prev_bucket\t$src\t$prev_dst\tremoved_stale\tbucket_change:${prev_bucket}->${tech_bucket}" >> "$SYMLINK_MANIFEST"
    fi
  fi
fi

# Create/update the desired link
make_link "$proj" "$tech_bucket" "$src" "$dst"

# Update state in memory
LAST_BUCKET["$state_key"]="$tech_bucket"
LAST_DST["$state_key"]="$dst"

  done < <(printf '%s\n' "${EXPECTED_ROWS[@]}")

done

# ---------- Summary ----------
# Persist latest classification state before summarizing
save_state

# Count files per tech (based on symlink_manifest actions)
awk -F'\t' 'NR==1{next}
  $6=="created" || $6=="exists_ok" || $6=="updated_link" || $6=="would_create" || $6=="would_update_link" {
    cnt[$2 FS $3]++
  }
  END{
    print "project\ttech_bucket\tlinked_files"
    for (k in cnt) print k "\t" cnt[k]
  }' "$SYMLINK_MANIFEST" \
| sort -k1,1 -k2,2 > "$SUMMARY_BY_TECH"

log "=== DONE organize ==="
log "Log:               $LOG"
log "Tech assignments:  $TECH_ASSIGN"
log "Symlink manifest:  $SYMLINK_MANIFEST"
log "Pending projects:  $PENDING_PROJ"
log "Summary by tech:   $SUMMARY_BY_TECH"
