# #!/usr/bin/env bash
# set -Eeuo pipefail

# #############################################
# # CONFIG
# #############################################
# BASE_PATH="/mnt/12T/chibao/data/official_data/fastq_by_tech/bulk_or_plate_RNA/offical_run/PRJNA887805/SAMN31185512/SRX17821404/SRR21832248"
# OUTPUT_PATH="/mnt/12T/chibao/data/official_data/fastq_by_tech/bulk_or_plate_RNA/official_result"

# # Only treat these as “general” (not already split) FASTQ names:
# INCLUDE_EXT_REGEX='\.(fastq|fq)(\.gz)?$'
# # Exclude files that already look split into lanes/reads:
# EXCLUDE_NAME_REGEX='(_R[123]_|\bI[12]\b|_R1\.|_R2\.|_R3\.|_I1\.|_I2\.)'


# # Runtime options
# DRY_RUN="false"      # "true" = simulate
# MAX_SAMPLE_CHECK=100000   # sample pairs to check header match (speed/accuracy tradeoff)

# #############################################
# # LOGGING / OUTPUT
# #############################################
# mkdir -p "$OUTPUT_PATH"
# LOG_DIR="$OUTPUT_PATH/_logs"
# mkdir -p "$LOG_DIR"

# MAIN_LOG="$LOG_DIR/splitting_$(date +%Y%m%d_%H%M%S).log"
# SUMMARY_TSV="$OUTPUT_PATH/split_summary.tsv"

# log()   { echo "[$(date -Iseconds)] $*" | tee -a "$MAIN_LOG" >&2; }
# die()   { echo "[$(date -Iseconds)] ERROR: $*" | tee -a "$MAIN_LOG" >&2; exit 1; }

# # Header for summary once
# if [[ ! -f "$SUMMARY_TSV" ]]; then
#   printf "ts\trun_id\tinput_path\tinput_reads\tR1_reads\tR2_reads\tsingletons\tmismatch_sampled\tstatus\tnotes\n" > "$SUMMARY_TSV"
# fi

# #############################################
# # DEPENDENCY CHECKS
# #############################################
# need_cmd() {
#   command -v "$1" >/dev/null 2>&1 || die "Missing required command: $1"
# }

# need_cmd reformat.sh
# need_cmd repair.sh
# if command -v pigz >/dev/null 2>&1; then
#   GZIP="pigz -c"
# else
#   GZIP="gzip -c"
# fi

# #############################################
# # UTILS
# #############################################

# # Normalize FASTQ headers: strip trailing "/1" or "/2" from the FIRST token only;
# # keep everything after the first space unchanged.
# normalize_headers() {
#   local in="$1" out="$2"
#   if [[ "$in" =~ \.gz$ ]]; then
#     zcat "$in" | awk '
#       NR%4==1 {
#         header=$0
#         sp=index(header," ")
#         if (sp==0) { key=header; rest="" } else { key=substr(header,1,sp-1); rest=substr(header,sp+1) }
#         sub(/[\/][12]$/, "", key)   # /1 or /2
#         sub(/[.][12]$/,  "", key)   # .1 or .2
#         print key (rest==""? "":" "rest)
#         next
#       }
#       { print }
#     ' | $GZIP > "$out"
#   else
#     awk '
#       NR%4==1 {
#         header=$0
#         sp=index(header," ")
#         if (sp==0) { key=header; rest="" } else { key=substr(header,1,sp-1); rest=substr(header,sp+1) }
#         sub(/[\/][12]$/, "", key)
#         sub(/[.][12]$/,  "", key)
#         print key (rest==""? "":" "rest)
#         next
#       }
#       { print }
#     ' "$in" | $GZIP > "$out"
#   fi
# }

# # Count FASTQ records quickly
# fastq_count() {
#   local f="$1"
#   if [[ "$f" =~ \.gz$ ]]; then
#     echo $(( $(zcat "$f" | wc -l) / 4 ))
#   else
#     echo $(( $(wc -l < "$f") / 4 ))
#   fi
# }

# # Sample header mismatch count between R1/R2 after normalization (expect 0)
# sample_mismatch() {
#   local r1="$1" r2="$2" N="$3"
#   paste \
#     <(zcat "$r1" | sed -n '1~4p' | cut -d " " -f1 | head -"$N") \
#     <(zcat "$r2" | sed -n '1~4p' | cut -d " " -f1 | head -"$N") \
#   | awk '{if($1!=$2) c++} END{print (c+0)}'
# }

# # Symlink to 10x-style names for Cell Ranger friendliness
# make_10x_symlinks() {
#   local outdir="$1" sample="$2" r1="$3" r2="$4"
#   local r1_10x="${outdir}/${sample}_S1_L001_R1_001.fastq.gz"
#   local r2_10x="${outdir}/${sample}_S1_L001_R2_001.fastq.gz"
#   ln -sf "$(realpath --relative-to="$outdir" "$r1")" "$r1_10x"
#   ln -sf "$(realpath --relative-to="$outdir" "$r2")" "$r2_10x"
# }

# #############################################
# # CORE: PROCESS ONE GENERAL FASTQ
# #############################################
# process_one() {
#   local in_fastq="$1"
#   local base name run outdir logf
#   base="$(basename "$in_fastq")"
#   name="${base%.*}"; name="${name%.*}"
#   run="$name"
#   outdir="$OUTPUT_PATH/$run"
#   mkdir -p "$outdir"
#   logf="$LOG_DIR/${run}.log"

#   log "=== [$run] Start ==="
#   echo "[$(date -Iseconds)] in=$in_fastq" >> "$logf"

#   # Step 0: count input reads (records, not pairs)
#   local in_reads
#   if [[ "$in_fastq" =~ \.gz$ ]]; then
#     in_reads=$(( $(zcat "$in_fastq" | wc -l) / 4 ))
#   else
#     in_reads=$(( $(wc -l < "$in_fastq") / 4 ))
#   fi
#   log "[$run] Input reads: $in_reads"

#   # Step 1: de-interleave (very low RAM)
#   local r1_tmp="$outdir/${run}.R1.tmp.fastq.gz"
#   local r2_tmp="$outdir/${run}.R2.tmp.fastq.gz"
#   log "[$run] Deinterleave with reformat.sh ..."
#   reformat.sh -Xmx8g in="$in_fastq" out1="$r1_tmp" out2="$r2_tmp" int=t ow=t >> "$logf" 2>&1

#   # Step 2: normalize headers
#   local r1_norm="$outdir/${run}.R1.norm.fastq.gz"
#   local r2_norm="$outdir/${run}.R2.norm.fastq.gz"
#   log "[$run] Normalizing headers ..."
#   normalize_headers "$r1_tmp" "$r1_norm"
#   normalize_headers "$r2_tmp" "$r2_norm"

#   # Step 3: sample mismatches; if 0, skip repair
#   local mism cR1 cR2
#   mism=$(sample_mismatch "$r1_norm" "$r2_norm" "$MAX_SAMPLE_CHECK")
#   log "[$run] Sampled header mismatches: $mism"

#   local r1_pair="$outdir/${run}.R1.paired.norm.fastq.gz"
#   local r2_pair="$outdir/${run}.R2.paired.norm.fastq.gz"
#   local singles="$outdir/${run}.singletons.fastq.gz"

#   if [[ "$mism" -eq 0 ]]; then
#     # Just rename normalized files to final
#     mv -f "$r1_norm" "$r1_pair"
#     mv -f "$r2_norm" "$r2_pair"
#     : > "$singles"  # empty placeholder for summary consistency
#   else
#     log "[$run] Nonzero mismatches; attempting repair pass ..."
#     # Try a modest heap first; make it configurable
#     : "${BBMAP_JAVA_OPTS:="-Xmx90g -Xms90g"}"
#     export BBMAP_JAVA_OPTS

#     if repair.sh in1="$r1_norm" in2="$r2_norm" out1="$r1_pair" out2="$r2_pair" outs="$singles" ow=t >> "$logf" 2>&1; then
#       log "[$run] Repair completed."
#     else
#       log "[$run] Repair failed; keeping normalized outputs as-is."
#       mv -f "$r1_norm" "$r1_pair" || true
#       mv -f "$r2_norm" "$r2_pair" || true
#       : > "$singles" || true
#     fi
#   fi

#   # Step 4: counts & symlinks
#   cR1=$(fastq_count "$r1_pair"); cR2=$(fastq_count "$r2_pair")
#   log "[$run] Final counts R1=$cR1 R2=$cR2"
#   make_10x_symlinks "$outdir" "$run" "$r1_pair" "$r2_pair"

#   # Step 5: cleanup temps
#   rm -f "$r1_tmp" "$r2_tmp" "$r1_norm" "$r2_norm" 2>/dev/null || true

#   # Step 6: summary
#   local status="OK" notes=""
#   if [[ "$mism" -ne 0 ]]; then status="WARN"; notes="repair_attempted"; fi
#   printf "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\n" \
#     "$(date -Iseconds)" "$run" "$in_fastq" "$in_reads" "$cR1" "$cR2" "NA" "$mism" "$status" "$notes" \
#     >> "$SUMMARY_TSV"

#   log "=== [$run] Done ==="
# }

# #############################################
# # DISCOVER & RUN
# #############################################
# log "Scan input: $BASE_PATH"
# mapfile -t CANDIDATES < <(
#   find -L "$BASE_PATH" -type f \
#     | grep -E "$INCLUDE_EXT_REGEX" \
#     | grep -Ev "$EXCLUDE_NAME_REGEX" \
#     | sort
# )

# if [[ ${#CANDIDATES[@]} -eq 0 ]]; then
#   die "No general FASTQ files found under $BASE_PATH (extensions .fastq[.gz]/.fq[.gz], excluding pre-split R1/R2/I1)."
# fi

# log "Found ${#CANDIDATES[@]} general FASTQ files."
# for f in "${CANDIDATES[@]}"; do
#   # process each file in isolation; continue on error
#   if ! process_one "$f"; then
#     log "[FAIL] $f"
#     printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
#       "$(date -Iseconds)" "$(basename "$f")" "$f" "NA" "NA" "NA" "NA" "NA" "FAIL" "see log" \
#       >> "$SUMMARY_TSV"
#     continue
#   fi
# done

# log "All done. Summary at: $SUMMARY_TSV"


# RUN=SRR21832248
# OD=/mnt/12T/chibao/data/official_data/fastq_by_tech/bulk_or_plate_RNA/official_result/$RUN
# R1="$OD/${RUN}.R1.paired.norm.fastq.gz"
# R2="$OD/${RUN}.R2.paired.norm.fastq.gz"

# echo "R1 first-token headers:"
# zcat "$R1" | sed -n '1~4p' | cut -d" " -f1 | head -20

# echo "R2 first-token headers:"
# zcat "$R2" | sed -n '1~4p' | cut -d" " -f1 | head -20

# echo "First 10 mismatches after current norm:"
# paste \
#   <(zcat "$R1" | sed -n '1~4p' | cut -d" " -f1 | head -100000) \
#   <(zcat "$R2" | sed -n '1~4p' | cut -d" " -f1 | head -100000) \
# | awk '$1!=$2{print NR,$1,$2; if(++c==10) exit}'

RUN=SRR21832248
OD=/mnt/12T/chibao/data/official_data/fastq_by_tech/bulk_or_plate_RNA/official_result/$RUN
R1_IN="$OD/${RUN}.R1.paired.norm.fastq.gz"
R2_IN="$OD/${RUN}.R2.paired.norm.fastq.gz"
R1_OUT="$OD/${RUN}.R1.paired.norm.v2.fastq.gz"
R2_OUT="$OD/${RUN}.R2.paired.norm.v2.fastq.gz"

# Deterministic renaming by pair index: @<RUN>:<index>
reindex () {
  in="$1"; out="$2"
  zcat "$in" | awk -v RUN="$RUN" '
    NR%4==1 {
      i=(NR+3)/4
      header=$0
      sp=index(header," ")
      rest=(sp? substr(header,sp+1):"")
      print "@"RUN":"i (rest==""? "":" "rest)
      next
    }
    { print }
  ' | pigz -c > "$out"
}

reindex "$R1_IN" "$R1_OUT"
reindex "$R2_IN" "$R2_OUT"

# Verify first 100k pairs: expect mismatches=0
paste \
  <(zcat "$R1_OUT" | sed -n '1~4p' | cut -d" " -f1 | head -100000) \
  <(zcat "$R2_OUT" | sed -n '1~4p' | cut -d" " -f1 | head -100000) \
| awk '{if($1!=$2)c++} END{print "mismatches=" (c+0)}'

# (Optional) replace originals and refresh 10x-style symlinks
mv -f "$R1_OUT" "$OD/${RUN}.R1.paired.norm.fastq.gz"
mv -f "$R2_OUT" "$OD/${RUN}.R2.paired.norm.fastq.gz"
ln -sf "${RUN}.R1.paired.norm.fastq.gz" "$OD/${RUN}_S1_L001_R1_001.fastq.gz"
ln -sf "${RUN}.R2.paired.norm.fastq.gz" "$OD/${RUN}_S1_L001_R2_001.fastq.gz"
