#!/usr/bin/env bash
set -Eeuo pipefail
# This script will move and organize into the logical structure the SRR IDs that is fetched from SRA and for output of cell ranger 

# Normalize SRR* directories based on source layout:
#  - If SRR* are under <PROJECT>/cellranger/count  (layered source)
#      -> DEST: <PROJECT>/cellranger/count/<SAMN>/<SRX>/<SRR>
#  - If SRR* are directly under <PROJECT>           (flat source)
#      -> DEST: <PROJECT>/<SAMN>/<SRX>/<SRR>
#
# Metadata columns (CSV or TSV): sample_accession, experiment_accession, run_accession

# bash /mnt/12T/chibao/code/01.preprocess/07_move.sh \
#     --project PRJNA683876 \
#     --root-base "/mnt/12T/chibao/data/official_data/fastq_by_tech/bulk_or_plate_RNA/official_cell" \
#     --metadata-root "/mnt/12T/chibao/data/official_data/metadata/official/split1" \
#     --log "/mnt/12T/chibao/data/official_data/logs/reorg_PRJNA683876.log" \
#     --dry-run true

ROOT_BASE="/mnt/12T/chibao/data/official_data/fastq_by_tech/bulk_or_plate_RNA/official_cell/PRJNA...."
METADATA_ROOT="/mnt/12T/chibao/data/official_data/metadata/official/split..."
DRY_RUN="true"
LOG="./reorg_cellranger_move_only.log"

PROJECT=""
META_PATH=""
COUNT_ROOT=""   # optional override to the parent holding SRR* (flat OR layered)

trap 'rc=$?; echo "[ERR] Line $LINENO failed: $BASH_COMMAND (exit $rc)" | tee -a "$LOG"; exit $rc' ERR

usage() {
  cat <<EOF
Usage:
  $(basename "$0") --project PRJNAxxxxxx
                   [--root-base DIR]
                   [--metadata-root DIR]
                   [--meta FILE]
                   [--count-root DIR]
                   [--dry-run true|false]
                   [--log FILE]

Notes:
  * Detects source layout:
      - flat source:  SRR* directly under <ROOT_BASE>/<PROJECT>
      - layered src:  SRR* directly under <ROOT_BASE>/<PROJECT>/cellranger/count
  * Destination depends on source layout:
      - from layered src -> <PROJECT>/cellranger/count/<SAMN>/<SRX>/<SRR>
      - from flat src    -> <PROJECT>/<SAMN>/<SRX>/<SRR>
  * Required metadata columns (case-insensitive): sample_accession, experiment_accession, run_accession
EOF
}

# ---------------- parse args ----------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --project) PROJECT="$2"; shift 2 ;;
    --root-base) ROOT_BASE="$2"; shift 2 ;;
    --metadata-root) METADATA_ROOT="$2"; shift 2 ;;
    --meta) META_PATH="$2"; shift 2 ;;
    --count-root) COUNT_ROOT="$2"; shift 2 ;;
    --dry-run) DRY_RUN="$2"; shift 2 ;;
    --log) LOG="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1"; usage; exit 1 ;;
  esac
done

[[ -z "$PROJECT" ]] && { echo "[ERR] --project required"; usage; exit 1; }

PROJECT_ROOT="${ROOT_BASE}/${PROJECT}"
mkdir -p "$(dirname "$LOG")"
: > "$LOG"

# ------------- find metadata ---------------
if [[ -z "${META_PATH:-}" ]]; then
  echo "[INFO] Auto-discovering metadata under $METADATA_ROOT" | tee -a "$LOG"
  META_PATH="$(find "$METADATA_ROOT" -type f -path "*/${PROJECT}/latest.tsv" | head -n1 || true)"
  [[ -z "$META_PATH" ]] && META_PATH="$(find "$METADATA_ROOT" -type f -path "*/${PROJECT}/*.tsv" | head -n1 || true)"
  [[ -z "$META_PATH" ]] && META_PATH="$(find "$METADATA_ROOT" -type f -path "*/${PROJECT}/*.csv" | head -n1 || true)"
fi
[[ -f "${META_PATH:-/dev/null}" ]] || { echo "[ERR] Metadata not found; pass --meta" | tee -a "$LOG"; exit 1; }

# ------------- decide source layout + COUNT_ROOT -------------
LAYOUT=""   # "flat" or "layered"
decide_source() {
  local explicit="${COUNT_ROOT:-}"
  if [[ -n "$explicit" ]]; then
    [[ -d "$explicit" ]] || { echo "[ERR] --count-root not found: $explicit" | tee -a "$LOG"; exit 1; }
    # Detect if explicit points to flat or layered by checking parent name
    if find "$explicit" -mindepth 1 -maxdepth 1 -type d -name 'SRR*' -print -quit >/dev/null; then
      # Heuristic: if explicit == <PROJECT_ROOT>/cellranger/count  => layered
      if [[ "$(basename "$explicit")" == "count" && "$(basename "$(dirname "$explicit")")" == "cellranger" ]]; then
        LAYOUT="layered"
      else
        # could still be flat; check if explicit == PROJECT_ROOT
        if [[ "$explicit" == "$PROJECT_ROOT" ]]; then
          LAYOUT="flat"
        else
          # fall back to content-based guess: if explicit ends with /count and has SRR*, treat layered; else flat
          if [[ "$(basename "$explicit")" == "count" ]]; then LAYOUT="layered"; else LAYOUT="flat"; fi
        fi
      fi
      echo "$explicit"
      return
    else
      echo "[ERR] No SRR* directories directly under --count-root: $explicit" | tee -a "$LOG"
      exit 1
    fi
  fi

  # No explicit override: probe common places
  local flat="${PROJECT_ROOT}"
  local layered="${PROJECT_ROOT}/cellranger/count"

  if [[ -d "$flat" ]] && find "$flat" -mindepth 1 -maxdepth 1 -type d -name 'SRR*' -print -quit >/dev/null; then
    LAYOUT="flat"
    echo "$flat"
    return
  fi
  if [[ -d "$layered" ]] && find "$layered" -mindepth 1 -maxdepth 1 -type d -name 'SRR*' -print -quit >/dev/null; then
    LAYOUT="layered"
    echo "$layered"
    return
  fi

  echo "[ERR] Could not find SRR* either in:
  - flat:    $flat
  - layered: $layered" | tee -a "$LOG"
  exit 1
}

COUNT_ROOT="$(decide_source)"

# ------------- pick destination base by layout -------------
# If source is layered -> keep layered destination base.
# If source is flat    -> use project root as destination base.
if [[ "$LAYOUT" == "layered" ]]; then
  DEST_BASE="${PROJECT_ROOT}/cellranger/count"
else
  DEST_BASE="${PROJECT_ROOT}"
fi

echo "[INFO] Project        : $PROJECT"       | tee -a "$LOG"
echo "[INFO] Source layout  : $LAYOUT"        | tee -a "$LOG"
echo "[INFO] Source root    : $COUNT_ROOT"    | tee -a "$LOG"
echo "[INFO] Dest base      : $DEST_BASE"     | tee -a "$LOG"
echo "[INFO] Metadata       : $META_PATH"     | tee -a "$LOG"
echo "[INFO] Dry-run        : $DRY_RUN"       | tee -a "$LOG"

# ------------- normalize metadata -> run \t sample \t experiment -------------
NORM="$(mktemp)"
awk '
  BEGIN{ OFS="\t"; FS="\t" }
  FNR==1{
    sub(/^\xef\xbb\xbf/,"",$0)
    if (index($0,",")>0 && index($0,"\t")==0) FS=","
    gsub(/\r/,"")
    for(i=1;i<=NF;i++){
      c=$i; gsub(/^[ \t]+|[ \t]+$/,"",c); tl=tolower(c)
      if(tl=="sample_accession") cs=i
      if(tl=="experiment_accession") ce=i
      if(tl=="run_accession") cr=i
    }
    next
  }
  {
    gsub(/\r/,"")
    if(cr && cs && ce){
      r=$cr; s=$cs; e=$ce
      gsub(/^[ \t]+|[ \t]+$/,"",r)
      gsub(/^[ \t]+|[ \t]+$/,"",s)
      gsub(/^[ \t]+|[ \t]+$/,"",e)
      if(r!="" && s!="" && e!="") print r, s, e
    }
  }' "$META_PATH" | sort -u > "$NORM"

[[ -s "$NORM" ]] || { echo "[ERR] Could not parse metadata headers/rows" | tee -a "$LOG"; exit 1; }
echo "[INFO] Unique runs in metadata: $(wc -l < "$NORM")" | tee -a "$LOG"

# ------------- discover SRR* in source -------------
mapfile -t RUN_DIRS < <(find "$COUNT_ROOT" -mindepth 1 -maxdepth 1 -type d -name 'SRR*' | sort)
if [[ ${#RUN_DIRS[@]} -eq 0 ]]; then
  echo "[WARN] No SRR* directories directly under $COUNT_ROOT" | tee -a "$LOG"
  exit 0
fi
echo "[INFO] SRR* directories found: ${#RUN_DIRS[@]}" | tee -a "$LOG"

# ------------- helpers -------------
lookup() { awk -v FS="\t" -v r="$1" '$1==r{print $2 "\t" $3; exit}' "$NORM"; }

do_mkdir() {
  if [[ "$DRY_RUN" == "true" ]]; then
    echo "[DRY] mkdir -p \"$1\"" | tee -a "$LOG"
  else
    mkdir -p "$1"
  fi
}

do_mv() {
  if [[ "$DRY_RUN" == "true" ]] ; then
    echo "[DRY] mv \"$1\" \"$2\"" | tee -a "$LOG"
  else
    mv "$1" "$2"
  fi
}

# ------------- plan -------------
echo "[PLAN] From -> To" | tee -a "$LOG"
MOVED=0; SKIPPED=0; MISSING=0; ALREADY=0

for rp in "${RUN_DIRS[@]}"; do
  r="$(basename "$rp")"
  se="$(lookup "$r")" || true
  if [[ -z "$se" ]]; then
    echo "[PLAN] $r -> (SKIP: no mapping in metadata)" | tee -a "$LOG"
    continue
  fi
  s="${se%%$'\t'*}"
  e="${se#*$'\t'}"
  echo "[PLAN] $r -> ${DEST_BASE}/${s}/${e}/${r}" | tee -a "$LOG"
done

# ------------- execute -------------
for rp in "${RUN_DIRS[@]}"; do
  r="$(basename "$rp")"
  se="$(lookup "$r")" || true
  if [[ -z "$se" ]]; then
    echo "[WARN] No mapping in metadata for ${r} — skipping" | tee -a "$LOG"
    ((MISSING++))
    continue
  fi
  s="${se%%$'\t'*}"
  e="${se#*$'\t'}"

  dst_dir="${DEST_BASE}/${s}/${e}"
  dst="${dst_dir}/${r}"

  if [[ "$rp" == "$dst" ]]; then
    echo "[INFO] Already placed: ${rp}" | tee -a "$LOG"
    ((ALREADY++))
    continue
  fi
  if [[ -e "$dst" ]]; then
    echo "[WARN] Destination exists, skip: ${dst}" | tee -a "$LOG"
    ((SKIPPED++))
    continue
  fi

  do_mkdir "$dst_dir"
  do_mv "$rp" "$dst"
  ((MOVED++)) || true
done

echo "[SUMMARY] moved=${MOVED} skipped=${SKIPPED} missing_meta=${MISSING} already_ok=${ALREADY}" | tee -a "$LOG"



# #!/usr/bin/env bash
# set -Eeuo pipefail

# # Move Cell Ranger runs:
# #   .../cellranger/count/SRR*  ->  .../cellranger/count/<SAMN>/<SRX>/<SRR>
# # Metadata columns (CSV or TSV): sample_accession, experiment_accession, run_accession

# ROOT_BASE="/mnt/12T/chibao/data/official_data/fastq_by_tech/bulk_or_plate_RNA/official_cell" # Edit as config
# METADATA_ROOT="/mnt/12T/chibao/data/official_data/metadata/official" 
# DRY_RUN="true"
# LOG="./reorg_cellranger_move_only.log"
# PROJECT=""
# META_PATH=""

# # ---- helpful trap so you SEE why it stops ----
# trap 'rc=$?; echo "[ERR] Line $LINENO failed: $BASH_COMMAND (exit $rc)" | tee -a "$LOG"; exit $rc' ERR

# usage() {
#   cat <<EOF
# Usage:
#   $(basename "$0") --project PRJNAxxxxxx [--root-base DIR] [--meta FILE]
#                    [--metadata-root DIR] [--dry-run true|false] [--log FILE]
# Notes:
#   * No symlinks; physical move only
#   * Handles CSV or TSV; strips BOM/CRLF; header names are case-insensitive
#   * Required columns: sample_accession, experiment_accession, run_accession
# EOF
# }

# # ---- parse args ----
# while [[ $# -gt 0 ]]; do
#   case "$1" in
#     --project) PROJECT="$2"; shift 2 ;;
#     --root-base) ROOT_BASE="$2"; shift 2 ;;
#     --meta) META_PATH="$2"; shift 2 ;;
#     --metadata-root) METADATA_ROOT="$2"; shift 2 ;;
#     --dry-run) DRY_RUN="$2"; shift 2 ;;
#     --log) LOG="$2"; shift 2 ;;
#     -h|--help) usage; exit 0 ;;
#     *) echo "Unknown arg: $1"; usage; exit 1 ;;
#   esac
# done

# [[ -z "$PROJECT" ]] && { echo "[ERR] --project required"; usage; exit 1; }

# PROJECT_ROOT="${ROOT_BASE}/${PROJECT}"
# COUNT_ROOT="${PROJECT_ROOT}/cellranger/count"
# [[ -d "$COUNT_ROOT" ]] || { echo "[ERR] Not found: $COUNT_ROOT"; exit 1; }

# mkdir -p "$(dirname "$LOG")"
# : > "$LOG"

# # ---- discover metadata if not provided ----
# if [[ -z "${META_PATH:-}" ]]; then
#   echo "[INFO] Auto-discovering metadata under $METADATA_ROOT" | tee -a "$LOG"
#   META_PATH="$(find "$METADATA_ROOT" -type f -path "*/${PROJECT}/latest.tsv" | head -n1 || true)"
#   [[ -z "$META_PATH" ]] && META_PATH="$(find "$METADATA_ROOT" -type f -path "*/${PROJECT}/*.tsv" | head -n1 || true)"
#   [[ -z "$META_PATH" ]] && META_PATH="$(find "$METADATA_ROOT" -type f -path "*/${PROJECT}/*.csv" | head -n1 || true)"
# fi
# [[ -f "${META_PATH:-/dev/null}" ]] || { echo "[ERR] Metadata not found; pass --meta" | tee -a "$LOG"; exit 1; }

# echo "[INFO] Project      : $PROJECT"    | tee -a "$LOG"
# echo "[INFO] Count root   : $COUNT_ROOT" | tee -a "$LOG"
# echo "[INFO] Metadata     : $META_PATH"  | tee -a "$LOG"
# echo "[INFO] Dry-run      : $DRY_RUN"    | tee -a "$LOG"

# # ---- normalize metadata to: run\t sample\texperiment ----
# NORM="$(mktemp)"
# awk '
#   BEGIN{ OFS="\t"; FS="\t" }
#   FNR==1{
#     sub(/^\xef\xbb\xbf/,"",$0)                   # BOM
#     if (index($0,",")>0 && index($0,"\t")==0) FS=",";  # sniff CSV
#     gsub(/\r/,"")                                # CRLF
#     for(i=1;i<=NF;i++){
#       c=$i; gsub(/^[ \t]+|[ \t]+$/,"",c); tl=tolower(c)
#       if(tl=="sample_accession") cs=i
#       if(tl=="experiment_accession") ce=i
#       if(tl=="run_accession") cr=i
#     }
#     next
#   }
#   {
#     gsub(/\r/,"")
#     if(cr && cs && ce){
#       r=$cr; s=$cs; e=$ce
#       gsub(/^[ \t]+|[ \t]+$/,"",r)
#       gsub(/^[ \t]+|[ \t]+$/,"",s)
#       gsub(/^[ \t]+|[ \t]+$/,"",e)
#       if(r!="" && s!="" && e!="") print r, s, e
#     }
#   }' "$META_PATH" | sort -u > "$NORM"

# [[ -s "$NORM" ]] || { echo "[ERR] Could not parse metadata headers/rows" | tee -a "$LOG"; exit 1; }
# echo "[INFO] Unique runs in metadata: $(wc -l < "$NORM")" | tee -a "$LOG"

# # ---- collect SRR* dirs one level under count ----
# mapfile -t RUN_DIRS < <(find "$COUNT_ROOT" -mindepth 1 -maxdepth 1 -type d -name 'SRR*' | sort)
# if [[ ${#RUN_DIRS[@]} -eq 0 ]]; then
#   echo "[WARN] No SRR* directories directly under $COUNT_ROOT" | tee -a "$LOG"
#   exit 0
# fi
# echo "[INFO] SRR* directories found: ${#RUN_DIRS[@]}" | tee -a "$LOG"

# # ---- helpers ----
# lookup() { awk -v FS="\t" -v r="$1" '$1==r{print $2 "\t" $3; exit}' "$NORM"; }

# do_mkdir() {
#   if [[ "$DRY_RUN" == "true" ]]; then
#     echo "[DRY] mkdir -p \"$1\"" | tee -a "$LOG"
#   else
#     mkdir -p "$1"
#   fi
# }

# do_mv() {
#   if [[ "$DRY_RUN" == "true" ]]; then
#     echo "[DRY] mv \"$1\" \"$2\"" | tee -a "$LOG"
#   else
#     mv "$1" "$2"
#   fi
# }

# # ---- preflight: show FULL plan so you can see all 8 mappings ----
# echo "[PLAN] From -> To" | tee -a "$LOG"
# MOVED=0
# SKIPPED=0
# MISSING=0
# ALREADY=0

# for rp in "${RUN_DIRS[@]}"; do
#   r="$(basename "$rp")"
#   se="$(lookup "$r")" || true
#   if [[ -z "$se" ]]; then
#     echo "[PLAN] $r -> (SKIP: no mapping)" | tee -a "$LOG"
#     continue
#   fi
#   s="${se%%$'\t'*}"
#   e="${se#*$'\t'}"
#   echo "[PLAN] $r -> ${s}/${e}/${r}" | tee -a "$LOG"
# done

# # ---- execute ----
# MOVED=0; SKIPPED=0; MISSING=0; ALREADY=0
# for rp in "${RUN_DIRS[@]}"; do
#   r="$(basename "$rp")"
#   se="$(lookup "$r")" || true
#   if [[ -z "$se" ]]; then
#     echo "[WARN] No mapping in metadata for ${r} — skipping" | tee -a "$LOG"
#     ((MISSING++))
#     continue
#   fi
#   s="${se%%$'\t'*}"
#   e="${se#*$'\t'}"
#   dst_dir="${COUNT_ROOT}/${s}/${e}"
#   dst="${dst_dir}/${r}"

#   if [[ "$rp" == "$dst" ]]; then
#     echo "[INFO] Already placed: ${rp}" | tee -a "$LOG"
#     ((ALREADY++))
#     continue
#   fi
#   if [[ -e "$dst" ]]; then
#     echo "[WARN] Destination exists, skip: ${dst}" | tee -a "$LOG"
#     ((SKIPPED++))
#     continue
#   fi

#   do_mkdir "$dst_dir"
#   do_mv "$rp" "$dst"
#   ((MOVED++)) || true
# done

# echo "[SUMMARY] moved=${MOVED} skipped=${SKIPPED} missing_meta=${MISSING} already_ok=${ALREADY}" | tee -a "$LOG"


# # # Paths
# # COUNT_ROOT="/mnt/12T/chibao/data/cellranger_data/PRJNA683876/cellranger/count"
# # TSV="/mnt/12T/chibao/data/official_data/metadata/official/split1/PRJNA683876/latest.tsv"

# # # A. Show first lines *with control chars* (reveals CRLF/BOM and commas vs tabs)
# # head -n 3 "$TSV" | sed -n 'l'

# # # B. Build a normalized (robust) run→(sample,exp) map from your metadata
# # awk '
# #   BEGIN{ OFS="\t"; FS="\t" }
# #   FNR==1{
# #     sub(/^\xef\xbb\xbf/,"",$0)           # strip BOM
# #     if (index($0, ",")>0 && index($0, "\t")==0) FS=",";  # sniff CSV vs TSV
# #     gsub(/\r/,"")                         # strip CR
# #     for(i=1;i<=NF;i++){
# #       t=$i; gsub(/^[ \t]+|[ \t]+$/,"",t); tl=tolower(t)
# #       if(tl=="run_accession") cr=i
# #       if(tl=="sample_accession") cs=i
# #       if(tl=="experiment_accession") ce=i
# #     }
# #     next
# #   }
# #   {
# #     gsub(/\r/,"")
# #     if(cr && cs && ce){
# #       r=$cr; s=$cs; e=$ce
# #       gsub(/^[ \t]+|[ \t]+$/,"",r)
# #       gsub(/^[ \t]+|[ \t]+$/,"",s)
# #       gsub(/^[ \t]+|[ \t]+$/,"",e)
# #       if(r!="" && s!="" && e!="") print r, s, e
# #     }
# #   }' "$TSV" | sort -u > /tmp/_meta_map.tsv

# # echo "Unique runs in metadata: $(wc -l < /tmp/_meta_map.tsv)"

# # # C. Compare filesystem SRRs vs metadata SRRs
# # cut -f1 /tmp/_meta_map.tsv | sort -u > /tmp/_meta_runs.txt
# # find "$COUNT_ROOT" -mindepth 1 -maxdepth 1 -type d -name 'SRR*' -printf '%f\n' | sort -u > /tmp/_fs_runs.txt

# # echo "== Present on disk but MISSING in metadata =="
# # comm -23 /tmp/_fs_runs.txt /tmp/_meta_runs.txt

# # echo "== Present in metadata but NOT on disk =="
# # comm -13 /tmp/_fs_runs.txt /tmp/_meta_runs.txt
