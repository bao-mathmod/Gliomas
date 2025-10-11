#!/usr/bin/env bash
set -euo pipefail
IFS=$'\n\t'

# ─── Config ───────────────────────────────────────────────────────────
MANIFEST="/mnt/18T/chibao/gliomas/data/avail/manifest.tsv" # TSV: accession<tab>subdir<tab>url
DATA_DIR="/mnt/18T/chibao/gliomas/data/avail"           # root for all downloads
LOG_DIR="/mnt/18T/chibao/gliomas/data/avail/logs"     # per‐job logs
SUMMARY="$DATA_DIR/download_summary.csv"                   # final CSV

# ─── Init directories ─────────────────────────────────────────────────
mkdir -p "$DATA_DIR" "$LOG_DIR"

# ─── Logger ────────────────────────────────────────────────────────────
log(){ printf "[%s] %s\n" "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }

# ─── Download function (wget only) ────────────────────────────────────
download_file(){
  local url=$1 outdir=$2 fname=$3
  log "⇢ wget -c \"$url\" -> \"$outdir/$fname\""
  wget -c -O "$outdir/$fname" "$url"
}

# ─── Decompression ─────────────────────────────────────────────────────
decompress(){
  local file=$1 dir=$2
  case "$file" in
    *.tar.gz|*.tgz)    log "⇢ Untarring $file";       tar -xvzf "$file"    -C "$dir" ;;
    *.tar.bz2)         log "⇢ Untarring bz2 $file";   tar -xvjf "$file"    -C "$dir" ;;
    *.tar)             log "⇢ Untarring $file";       tar -xvf  "$file"    -C "$dir" ;;
    *.zip)             log "⇢ Unzipping $file";       unzip -o  "$file"    -d "$dir" ;;
    *.gz)              log "⇢ Gunzip $file";          gunzip -kf "$file"        ;;
    *)                 log "⇢ No decompression: $file" ;;
  esac
}

# ─── Prepare summary CSV ────────────────────────────────────────────────
echo "accession,subdir,filename,path,size_bytes,md5sum" > "$SUMMARY"

# ─── Main loop ─────────────────────────────────────────────────────────
log "Starting download & organize workflow"
# FIX: Use tail -n +2 to skip the header line of the manifest file
tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r ACCESSION SUBDIR URL; do
  log "→ [$ACCESSION | $SUBDIR]"
  SAMPLE_DIR="$DATA_DIR/$ACCESSION/$SUBDIR"
  mkdir -p "$SAMPLE_DIR"

  # per‐sample info.txt
  INFO="$DATA_DIR/$ACCESSION/info.txt"
  {
    echo "Accession: $ACCESSION"
    echo "Subdir:    $SUBDIR"
    echo "URL:       $URL"
    echo "Started:   $(date -u '+%Y-%m-%dT%H:%M:%SZ')"
  } > "$INFO"

  # derive filename: everything after “file=” or fallback to ACCESSION_RAW.tar
  if [[ "$URL" == *file=* ]]; then
    FNAME="${URL##*file=}"
  else
    FNAME="${ACCESSION}_RAW.tar"
  fi

  # download & log
  LOGF="$LOG_DIR/${ACCESSION}_${SUBDIR//\//_}.log"
  download_file "$URL" "$SAMPLE_DIR" "$FNAME" 2>&1 | tee "$LOGF"

  # decompress if needed
  DL="$SAMPLE_DIR/$FNAME"
  decompress "$DL" "$SAMPLE_DIR" 2>&1 | tee -a "$LOGF"

  # record files in info.txt & summary.csv
  log "Listing files in $SAMPLE_DIR"
  for f in "$SAMPLE_DIR"/*; do
    [ -e "$f" ] || continue
    echo " - $(basename "$f")" >> "$INFO"
    size=$(stat -c%s "$f")
    md5=$(md5sum "$f" | cut -d' ' -f1)
    echo "$ACCESSION,$SUBDIR,$(basename "$f"),$f,$size,$md5" >> "$SUMMARY"
  done

  echo "Completed: $(date -u '+%Y-%m-%dT%H:%M:%SZ')" >> "$INFO"
done

log "All done! Summary saved to $SUMMARY"