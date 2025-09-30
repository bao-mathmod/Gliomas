#!/usr/bin/env bash
set -Eeuo pipefail
IFS=$'\n\t'

#############################################
# CONFIG (EDIT THESE)
#############################################
BASE_PATH="/mnt/12T/chibao/data/official_data/fastq_by_tech/bulk_or_plate_RNA/offical_run"
RAW_OUT="/mnt/12T/chibao/data/official_data/fastq_by_tech/bulk_or_plate_RNA/official_result"
CR_FASTQ_ROOT="/mnt/12T/chibao/data/official_data/fastq_by_tech/bulk_or_plate_RNA/official_cell"
CR_MANIFEST="${CR_FASTQ_ROOT}/_cr_samples.tsv"

# Performance knobs
MAX_JOBS=3          # process up to this many SRRs concurrently
THREADS=12          # threads per fasterq-dump job
PIGZ_THREADS=6      # threads per pigz in each job
TEMP_ROOT="/mnt/12T/chibao/data/stuff_data/single/scRNA/temp"  # fast temp (NVMe if available)

# Extra flags
DRY_RUN="false"
FASTERQ_EXTRA=( --progress )

#############################################
# INIT & UTILITIES
#############################################
mkdir -p "$RAW_OUT" "$CR_FASTQ_ROOT" "$TEMP_ROOT"
MASTER_LOG="${RAW_OUT}/sra_to_crfastq_$(date +%Y%m%d_%H%M%S).log"
SUMMARY="${RAW_OUT}/sra_to_crfastq.summary.tsv"

log(){ echo "[$(date -Iseconds)] $*" | tee -a "$MASTER_LOG" >&2; }
die(){ echo "[$(date -Iseconds)] ERROR: $*" | tee -a "$MASTER_LOG" >&2; exit 1; }

need(){ command -v "$1" >/dev/null || die "$1 not found in PATH"; }
need prefetch
need fasterq-dump
need pigz
if ! command -v vdb-dump >/dev/null; then log "WARN: vdb-dump not found; will skip that check"; fi

[[ -f "$CR_MANIFEST" ]] || printf "project\tsample\tfastq_dir\n" > "$CR_MANIFEST"
[[ -f "$SUMMARY"    ]] || printf "ts\taccession\tchosen_R1\tlen_R1\tchosen_R2\tlen_R2\tout_fastq_dir\tstatus\tnote\n" > "$SUMMARY"

# atomic append with flock
safe_append(){
  local file="$1"; shift
  (
    flock -w 60 9 || { echo "WARN: lock timeout for $file" >&2; exit 0; }
    printf "%s\n" "$*" >> "$file"
  ) 9>>"${file}.lock"
}

first_len(){
  local f="$1"
  if [[ "$f" =~ \.gz$ ]]; then zcat "$f" | sed -n '2p' | awk '{print length($0)}'
  else                         sed -n '2p' "$f" | awk '{print length($0)}'
  fi
}
mean_len(){
  local f="$1" N="${2:-10000}"
  if [[ "$f" =~ \.gz$ ]]; then zcat "$f" | sed -n '2~4p' | head -"$N" | awk '{s+=length($0)} END{if(NR) printf("%.1f\n", s/NR); else print 0}'
  else                         sed -n '2~4p' | head -"$N" | awk '{s+=length($0)} END{if(NR) printf("%.1f\n", s/NR); else print 0}'
  fi
}

# Decide R1/R2 among candidate files by first-read length
choose_r1_r2(){
  local -a files=("$@")
  declare -A L; local f l r1="" r2="" len1=0 len2=0
  for f in "${files[@]}"; do [[ -s "$f" ]] || continue; L["$f"]="$(first_len "$f" 2>/dev/null || echo 0)"; done
  for f in "${files[@]}"; do l="${L[$f]}"; [[ -z "$l" ]] && continue; if (( l>=26 && l<=36 )); then r1="$f"; len1="$l"; fi; done
  for f in "${files[@]}"; do l="${L[$f]}"; [[ -z "$l" ]] && continue; if (( l>=70 ));           then r2="$f"; len2="$l"; fi; done
  [[ -n "$r1" && -n "$r2" ]] && echo "$r1 $len1 $r2 $len2" || echo ""
}

gzip_glob(){
  local glob="$1"; shopt -s nullglob; local -a plains=( ${glob} )
  (( ${#plains[@]} )) && pigz -p "$PIGZ_THREADS" "${plains[@]}"
}

# Validate Cell Ranger-style files exist & look like 10x (R1≈28–30, R2≥70)
cr_files_ok(){
  local outdir="$1" acc="$2"
  local r1="${outdir}/${acc}_S1_L001_R1_001.fastq.gz"
  local r2="${outdir}/${acc}_S1_L001_R2_001.fastq.gz"
  [[ -s "$r1" && -s "$r2" ]] || return 1
  local l1 l2; l1="$(first_len "$r1" || echo 0)"; l2="$(first_len "$r2" || echo 0)"
  (( l1>=26 && l1<=36 && l2>=70 )) || return 1
  return 0
}

# NEW: Auto-detect project from the SRR path
# Priority: PRJNA/PRJEB/PRJDB -> SRP -> parent dir name
detect_project(){
  local srr_dir="$1"
  local path="$srr_dir"
  # scan path components from right to left
  IFS='/' read -r -a parts <<< "$path"
  local i
  for (( i=${#parts[@]}-1; i>=0; i-- )); do
    comp="${parts[$i]}"
    if [[ "$comp" =~ ^PRJN(A|EB|DB)[0-9]+$ ]]; then
      echo "$comp"; return 0
    fi
  done
  for (( i=${#parts[@]}-1; i>=0; i-- )); do
    comp="${parts[$i]}"
    if [[ "$comp" =~ ^SRP[0-9]+$ ]]; then
      echo "$comp"; return 0
    fi
  done
  # fallback to immediate parent directory (e.g., SRX.../SRR..., so parent is SRX... or project folder name)
  basename "$(dirname "$srr_dir")"
}

#############################################
# ONE-SAMPLE PIPE (with SKIP/REUSE)
#############################################
process_one(){
  local srr_path="$1"
  local accession; accession="$(basename "$srr_path")"
  local joblog="${RAW_OUT}/${accession}.log"
  local outdir="${CR_FASTQ_ROOT}/${accession}"
  local tempdir="${TEMP_ROOT}/${accession}_$$"
  mkdir -p "$outdir" "$tempdir"

  # auto project
  local project; project="$(detect_project "$srr_path")"

  {
    echo "[$(date -Iseconds)] === [$accession] START === (project=${project})"

    # ----- SKIP if CR-style outputs already valid -----
    if cr_files_ok "$outdir" "$accession"; then
      local r1="${outdir}/${accession}_S1_L001_R1_001.fastq.gz"
      local r2="${outdir}/${accession}_S1_L001_R2_001.fastq.gz"
      local l1="$(first_len "$r1")"; local l2="$(first_len "$r2")"
      echo "[$(date -Iseconds)] SKIP: valid CR files (R1≈$l1, R2≈$l2) at $outdir"
      if ! grep -Fq $'\t'"$accession"$'\t'"$outdir" "$CR_MANIFEST" 2>/dev/null; then
        safe_append "$CR_MANIFEST" "$project\t$accession\t$outdir"
      fi
      safe_append "$SUMMARY" "$(date -Iseconds)\t$accession\t$(basename "$r1")\t$l1\t$(basename "$r2")\t$l2\t$outdir\tSKIP_OK\talready_linked_and_valid"
      echo "[$(date -Iseconds)] === [$accession] END (SKIP) ==="
      return 0
    fi

    # ----- REUSE if RAW_OUT has valid _3/_4; just link -----
    shopt -s nullglob
    local cands=( "$RAW_OUT/${accession}"_*.fastq.gz )
    if (( ${#cands[@]} > 0 )); then
      local choice; choice="$(choose_r1_r2 "${cands[@]}")"
      if [[ -n "$choice" ]]; then
        local R1 R2 LEN1 LEN2; read -r R1 LEN1 R2 LEN2 <<< "$choice"
        ln -sf "$R1" "${outdir}/${accession}_S1_L001_R1_001.fastq.gz"
        ln -sf "$R2" "${outdir}/${accession}_S1_L001_R2_001.fastq.gz"
        local ML1 ML2; ML1="$(mean_len "$R1" 10000)"; ML2="$(mean_len "$R2" 10000)"
        echo "[$(date -Iseconds)] REUSE: linked existing dumps R1=$(basename "$R1") (~$LEN1), R2=$(basename "$R2") (~$LEN2) → $outdir"
        if ! grep -Fq $'\t'"$accession"$'\t'"$outdir" "$CR_MANIFEST" 2>/dev/null; then
          safe_append "$CR_MANIFEST" "$project\t$accession\t$outdir"
        fi
        safe_append "$SUMMARY" "$(date -Iseconds)\t$accession\t$(basename "$R1")\t$LEN1\t$(basename "$R2")\t$LEN2\t$outdir\tREUSE\tlinked_existing"
        echo "[$(date -Iseconds)] === [$accession] END (REUSE) ==="
        return 0
      fi
    fi

    # ----- Otherwise: prefetch → fasterq-dump → gzip → choose → link -----
    if command -v vdb-dump >/dev/null 2>&1; then
      echo "[$(date -Iseconds)] vdb-dump peek (first 20 spots)"
      vdb-dump -R 1 -C READ_LEN,READ_TYPE,READ_FILTER,READ_INDEX -T SEQUENCE "$accession" -n 20 || true
    fi

    if [[ "$DRY_RUN" == "true" ]]; then
      echo "DRY prefetch $accession"
      echo "DRY fasterq-dump --split-files --include-technical --threads $THREADS --temp $tempdir --outdir $RAW_OUT $accession"
    else
      echo "prefetch $accession"
      prefetch "$accession"
      echo "fasterq-dump $accession ..."
      fasterq-dump --split-files --include-technical \
        --threads "$THREADS" --temp "$tempdir" --outdir "$RAW_OUT" \
        "${FASTERQ_EXTRA[@]}" "$accession"
      gzip_glob "$RAW_OUT/${accession}_*.fastq"
    fi

    local files=( "$RAW_OUT/${accession}"_*.fastq.gz )
    if (( ${#files[@]} == 0 )); then
      echo "[$(date -Iseconds)] FAIL: no fastq.gz after fasterq-dump."
      safe_append "$SUMMARY" "$(date -Iseconds)\t$accession\t\t0\t\t0\t\tFAIL\tno_fastqs"
      rm -rf "$tempdir"
      return 0
    fi
    local choice2; choice2="$(choose_r1_r2 "${files[@]}")"
    if [[ -z "$choice2" ]]; then
      echo "[$(date -Iseconds)] FAIL: cannot identify 10x R1/R2 by length."
      for f in "${files[@]}"; do echo " - $(basename "$f") first_len=$(first_len "$f")"; done
      safe_append "$SUMMARY" "$(date -Iseconds)\t$accession\t\t0\t\t0\t\tFAIL\tcannot_select_R1R2"
      rm -rf "$tempdir"
      return 0
    fi

    local R1 R2 LEN1 LEN2; read -r R1 LEN1 R2 LEN2 <<< "$choice2"
    ln -sf "$R1" "${outdir}/${accession}_S1_L001_R1_001.fastq.gz"
    ln -sf "$R2" "${outdir}/${accession}_S1_L001_R2_001.fastq.gz"
    local ML1 ML2; ML1="$(mean_len "$R1" 10000)"; ML2="$(mean_len "$R2" 10000)"
    echo "[$(date -Iseconds)] OK: R1=$(basename "$R1") (~$LEN1), R2=$(basename "$R2") (~$LEN2) → $outdir"

    if ! grep -Fq $'\t'"$accession"$'\t'"$outdir" "$CR_MANIFEST" 2>/dev/null; then
      safe_append "$CR_MANIFEST" "$project\t$accession\t$outdir"
    fi
    safe_append "$SUMMARY" "$(date -Iseconds)\t$accession\t$(basename "$R1")\t$LEN1\t$(basename "$R2")\t$LEN2\t$outdir\tOK\tfresh_dump"

    rm -rf "$tempdir"
    echo "[$(date -Iseconds)] === [$accession] END (OK) ==="
  } >> "$joblog" 2>&1

  # short note to master log
  re=$'\t'"$accession"$'\t.*\t(OK|REUSE|SKIP_OK)'
  if grep -q -E "$re" "$SUMMARY" 2>/dev/null; then
    status="$(grep -m1 -E "$re" "$SUMMARY" | awk -F'\t' '{print $8}')"
    log "[$accession] DONE ($status) → $outdir  (project=${project})"
  else
    log "[$accession] FAIL (see $joblog)  (project=${project})"
  fi
}

export -f process_one log safe_append first_len mean_len choose_r1_r2 gzip_glob cr_files_ok detect_project
export RAW_OUT CR_FASTQ_ROOT CR_MANIFEST THREADS PIGZ_THREADS TEMP_ROOT DRY_RUN SUMMARY MASTER_LOG FASTERQ_EXTRA

#############################################
# DISCOVER & RUN IN PARALLEL
#############################################
log "START: BASE_PATH=$BASE_PATH  RAW_OUT=$RAW_OUT  CR_FASTQ_ROOT=$CR_FASTQ_ROOT  MAX_JOBS=$MAX_JOBS  THREADS=$THREADS"
mapfile -t SRR_DIRS < <(find "$BASE_PATH" -type d -name "SRR*" | sort)
(( ${#SRR_DIRS[@]} )) || die "No SRR directories found under $BASE_PATH"
log "Found ${#SRR_DIRS[@]} SRR directories."

active=0
for d in "${SRR_DIRS[@]}"; do
  bash -lc 'process_one "$@"' _ "$d" &
  (( active+=1 ))
  if (( active >= MAX_JOBS )); then
    wait -n
    (( active-=1 ))
  fi
done
wait
log "ALL DONE. Master log: $MASTER_LOG  Summary: $SUMMARY  Manifest: $CR_MANIFEST"
