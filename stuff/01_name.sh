# This bash script is desined to prepare data for Cell Ranger 
# They will rename the FastQ files to Cell Ranger format and create symlinks in the suitable directory
# Input: 
# - Expected index file: _expected_files.tsv
# - Raw FastQ files directory: /mnt/rdisk/gliomas/data/official_data/fastq
# Output:
# - Cell Ranger FastQ directory: /mnt/rdisk/gliomas/data/official_data/fastq_cellranger
# - Sample manifest for Cell Ranger: _cr_samples.tsv
# - Symlink plan: _symlink_plan.tsv
# Usage: This will base on the path of expected index file to create symlinks for each FASTQ files => This link will lead to the root FASTQ Directory

# #!/usr/bin/env bash
# set -Euo pipefail
# IFS=$'\n\t'

# # ===== Config =====
# EXPECTED_INDEX="/mnt/12T/chibao/cell_ranger/_expected_files.tsv"
# RAW_FASTQ_ROOT="/mnt/12T/chibao/data/official_data/fastq"
# CR_FASTQ_ROOT="/mnt/12T/chibao/cell_ranger/fastq_cellranger"   # new workspace for symlinks
# LOG="$CR_FASTQ_ROOT/_prepare.log"
# MANIFEST="$CR_FASTQ_ROOT/_cr_samples.tsv"   # sample-level manifest for Cell Ranger
# PLAN="$CR_FASTQ_ROOT/_symlink_plan.tsv"     # per-read plan (debug/trace)

# mkdir -p "$CR_FASTQ_ROOT"
# : > "$LOG"
# echo -e "project_id\tsample\tfastq_dir" > "$MANIFEST"
# echo -e "project_id\tsample\trun\tlane\tread\tsrc\tdst" > "$PLAN"

# ts(){ date '+%Y-%m-%d %H:%M:%S%z'; }
# log(){ printf '[%s] %s\n' "$(ts)" "$*" | tee -a "$LOG" >&2; }

# [[ -s "$EXPECTED_INDEX" ]] || { log "ERROR: expected index not found: $EXPECTED_INDEX"; exit 1; }

# # Build symlink plan with lane assignment per (project,sample,run)
# # - lane is assigned sequentially per sample
# # - read determined by filename suffix (_1|_R1 -> R1, _2|_R2 -> R2)
# awk -F'\t' 'BEGIN{OFS="\t"}
#   NR==1{next}
#   {
#     proj=$1; study=$2; sample=$3; experiment_id=$4; run=$5; file=$6; src=$9
#     # Determine read (R1/R2)
#     r=0
#     if (file ~ /(_R?1)(\.fastq(\.gz)?)$/) r=1
#     else if (file ~ /(_R?2)(\.fastq(\.gz)?)$/) r=2
#     else if (file ~ /_1\.f(ast)?q(\.gz)?$/) r=1
#     else if (file ~ /_2\.f(ast)?q(\.gz)?$/) r=2
#     if (!r) next

#     key_samp = proj "|" sample
#     key_run  = proj "|" sample "|" run
#     if (!(key_run in lane_of_run)) {
#       lane_of_run[key_run] = ++lane_count[key_samp]
#     }
#     lane = lane_of_run[key_run]
#     # lane padded to 3 digits for Cell Ranger (e.g. 001, 002, ...)
#     printf "%s\t%s\t%s\t%03d\tR%d\t%s\n", proj, sample, run, lane, r, src
#   }
# ' "$EXPECTED_INDEX" \
# | sort -k1,1 -k2,2 -k3,3 -k5,5 \
# | while IFS=$'\t' read -r proj sample run lane read src; do
#     # Destination dirs by project/sample
#     ddir="$CR_FASTQ_ROOT/$proj/$sample"
#     mkdir -p "$ddir"

#     # Cell Ranger–style dest name
#     dst="$ddir/${sample}_S1_L${lane}_${read}_001.fastq.gz"

#     # Record plan
#     echo -e "${proj}\t${sample}\t${run}\t${lane}\t${read}\t${src}\t${dst}" >> "$PLAN"

#     # Create/update symlink (force replace if exists)
#     if [[ -e "$dst" || -L "$dst" ]]; then rm -f "$dst"; fi
#     ln -s "$src" "$dst"
# done

# # Build sample-level manifest (one row per sample)
# # Portable path dir extraction (no gawk-only gensub):
# # Take $7 (dst path), strip the last "/basename"
# awk -F'\t' 'NR>1{
#   p = $7
#   # remove trailing filename: keep up to last "/" (or whole if no slash)
#   lastslash = 0
#   for (i=1; i<=length(p); i++) { if (substr(p,i,1)=="/") lastslash=i }
#   if (lastslash>0) dir=substr(p,1,lastslash-1); else dir=p
#   print $1 "\t" $2 "\t" dir
# }' "$PLAN" | sort -u >> "$MANIFEST"

# log "Prepared Cell Ranger symlinks under: $CR_FASTQ_ROOT"
# log "Sample manifest:  $MANIFEST"
# log "Symlink plan:     $PLAN"

#!/usr/bin/env bash
set -Euo pipefail
IFS=$'\n\t'

# ===== Config =====
EXPECTED_INDEX="/mnt/12T/chibao/cell_ranger/_expected_files.tsv"
RAW_FASTQ_ROOT="/mnt/12T/chibao/data/official_data/fastq"
CR_FASTQ_ROOT="/mnt/12T/chibao/cell_ranger/fastq_cellranger"   # new workspace for symlinks
LOG="$CR_FASTQ_ROOT/_prepare.log"
MANIFEST="$CR_FASTQ_ROOT/_cr_samples.tsv"   # sample-level manifest for Cell Ranger
PLAN="$CR_FASTQ_ROOT/_symlink_plan.tsv"     # per-read plan (debug/trace)
ORPHANS="$CR_FASTQ_ROOT/_orphans.tsv"       # mates missing (R1 or R2 absent)

mkdir -p "$CR_FASTQ_ROOT"
: > "$LOG"
echo -e "project_id\tsample\tfastq_dir" > "$MANIFEST"
echo -e "project_id\tsample\trun\tlane\tread\tsrc\tdst" > "$PLAN"
echo -e "project_id\tsample\trun\tlane\tmissing\tpresent_src" > "$ORPHANS"

ts(){ date '+%Y-%m-%d %H:%M:%S%z'; }
log(){ printf '[%s] %s\n' "$(ts)" "$*" | tee -a "$LOG" >&2; }

[[ -s "$EXPECTED_INDEX" ]] || { log "ERROR: expected index not found: $EXPECTED_INDEX"; exit 1; }

# Build symlink plan with strong filtering & paired-only emission
# Your EXPECTED_INDEX columns:
# 1=split_id  2=project_id  3=study  4=sample  5=experiment  6=run
# 7=file      8=url         9=md5_expected  10=target_path
awk -F'\t' -v ORPHANS="$ORPHANS" 'BEGIN{OFS="\t"}
  NR==1{next}
  {
    proj=$2; study=$3; sample=$4; experiment_id=$5; run=$6; file=$7; src=$10

    # Require a real source path; skip if missing
    if (src=="" || src=="-") next

    # ---- Determine read (R1/R2) and reject non-reads (I1/I2, R3, unlabeled) ----
    r=0
    if (file ~ /(_I[12])(\.f(ast)?q(\.gz)?)$/) next   # skip index reads
    if (file ~ /(_R?3)(\.f(ast)?q(\.gz)?)$/)  next   # skip any R3

    if (file ~ /(_R?1)(\.f(ast)?q(\.gz)?)$/) r=1
    else if (file ~ /(_R?2)(\.f(ast)?q(\.gz)?)$/) r=2
    else if (file ~ /_1\.f(ast)?q(\.gz)?$/) r=1
    else if (file ~ /_2\.f(ast)?q(\.gz)?$/) r=2
    else next  # e.g. SRR10315838.fastq.gz -> ignored by design

    key_samp = proj "|" sample
    key_run  = proj "|" sample "|" run

    # Assign lane per (proj,sample,run)
    if (!(key_run in lane_of_run)) {
      lane_of_run[key_run] = ++lane_count[key_samp]
    }
    lane = lane_of_run[key_run]

    key_pair = key_run "|" sprintf("%03d", lane)
    srcs[key_pair "|R" r] = src
    files[key_pair "|R" r] = file      # keep file name to decide extension if needed
    meta[key_pair] = proj "\t" sample "\t" run "\t" sprintf("%03d", lane)
  }
  END{
    for (k in meta) {
      split(meta[k], m, "\t"); proj=m[1]; sample=m[2]; run=m[3]; lane=m[4]
      s1 = srcs[k "|R1"]; s2 = srcs[k "|R2"];

      if (s1!="" && s2!="") {
        printf "%s\t%s\t%s\t%s\t%s\t%s\n", proj, sample, run, lane, "R1", s1
        printf "%s\t%s\t%s\t%s\t%s\t%s\n", proj, sample, run, lane, "R2", s2
      } else if (s1!="" || s2!="") {
        missing = (s1=="" ? "R1" : "R2")
        present = (s1!="" ? s1 : s2)
        printf "%s\t%s\t%s\t%s\t%s\t%s\n", proj, sample, run, lane, missing, present >> ORPHANS
      }
    }
  }
' "$EXPECTED_INDEX" \
| sort -k1,1 -k2,2 -k3,3 -k5,5 \
| while IFS=$'\t' read -r proj sample run lane read src; do
    ddir="$CR_FASTQ_ROOT/$proj/$sample"
    mkdir -p "$ddir"

    # Extension: prefer the src path; if it lacks extension, fallback to read name
    if [[ "$src" =~ \.fastq\.gz$ ]]; then ext=".fastq.gz"
    elif [[ "$src" =~ \.fq\.gz$ ]]; then ext=".fq.gz"
    elif [[ "$src" =~ \.fastq$ ]]; then ext=".fastq"
    elif [[ "$src" =~ \.fq$ ]]; then ext=".fq"
    else ext=".fastq.gz"  # conservative default
    fi

    dst="$ddir/${sample}_S1_L${lane}_${read}_001${ext}"

    echo -e "${proj}\t${sample}\t${run}\t${lane}\t${read}\t${src}\t${dst}" >> "$PLAN"

    [[ -e "$dst" || -L "$dst" ]] && rm -f "$dst"
    ln -s "$src" "$dst"
done

# Build sample-level manifest (one row per sample)
awk -F'\t' 'NR>1{
  p = $7
  # remove trailing filename: keep up to last "/"
  lastslash = 0
  for (i=1; i<=length(p); i++) { if (substr(p,i,1)=="/") lastslash=i }
  if (lastslash>0) dir=substr(p,1,lastslash-1); else dir=p
  print $1 "\t" $2 "\t" dir
}' "$PLAN" | sort -u >> "$MANIFEST"

# Small summary
pairs=$(wc -l < "$PLAN"); pairs=$((pairs-1))   # subtract header count
orphs=$(wc -l < "$ORPHANS"); orphs=$((orphs-1))
log "Prepared Cell Ranger symlinks under: $CR_FASTQ_ROOT"
log "Sample manifest:  $MANIFEST"
log "Symlink plan:     $PLAN  (rows incl. R1+R2 = $pairs)"
log "Orphans (skipped): $ORPHANS  (rows = $orphs)"
