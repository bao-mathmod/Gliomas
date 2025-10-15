#!/usr/bin/env bash
set -Euo pipefail

# ===== CONFIGURATION (EDIT THESE PATHS) =====

# 1. The top-level directory containing your project's FASTQ files.
#    The script will search for sample folders (e.g., SAMN*) inside this path.
SEARCH_DIR="/mnt/18T/chibao/gliomas/data/fastq/official/PRJNA1081384"

# 2. A directory where the clean, organized symbolic links will be created.
#    This keeps your original data untouched. It will be created if it doesn't exist.
#    The final 'fastq_dir' column in your samplesheet will point to folders inside here.
LINK_ROOT="/mnt/18T/chibao/gliomas/data/fastq/cellranger_ready"

# 3. The name of the final output samplesheet file.
OUTPUT_SHEET="snRNA_samplesheet.tsv"

# ===== END OF CONFIGURATION =====

# --- Helper functions ---
log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
}

# --- Sanity Checks ---
if [[ ! -d "$SEARCH_DIR" ]]; then
  log "ERROR: Search directory not found: $SEARCH_DIR"
  exit 1
fi

# --- Main Logic ---
log "Starting samplesheet creation..."
log "Source FASTQ path: $SEARCH_DIR"
log "Organized links will be in: $LINK_ROOT"

# Create the root for our organized links
mkdir -p "$LINK_ROOT"
log "Created link directory at $LINK_ROOT"

# Get the project ID from the search directory name
PROJECT_ID=$(basename "$SEARCH_DIR")

# Prepare the output file with the correct header
printf "project_id\tsample\tfastq_dir\n" > "$OUTPUT_SHEET"

# Find all unique sample directories (assuming they start with 'SAMN')
# The `find ... -print -quit` trick gets the first match, and `dirname` gets its parent.
# This finds all directories at the level of 'SAMN40176204'.
find "$SEARCH_DIR" -mindepth 1 -maxdepth 1 -type d -name 'SAMN*' | sort -u | while IFS= read -r sample_path; do
  SAMPLE_ID=$(basename "$sample_path")
  
  log "Processing sample: $SAMPLE_ID"

  # Define the new, clean directory for this sample's FASTQ links
  clean_fastq_dir="$LINK_ROOT/$SAMPLE_ID"
  mkdir -p "$clean_fastq_dir"

  # Find ALL fastq.gz files for this sample, no matter how nested they are
  # and create a symbolic link in the clean directory.
  file_count=0
  while IFS= read -r fastq_file; do
    ln -sf "$fastq_file" "$clean_fastq_dir/"
    ((file_count++))
  done < <(find "$sample_path" -type f -name '*.fastq.gz')
  
  if (( file_count > 0 )); then
    log "-> Linked $file_count FASTQ files into $clean_fastq_dir"
    # Append the entry to the samplesheet
    printf "%s\t%s\t%s\n" "$PROJECT_ID" "$SAMPLE_ID" "$clean_fastq_dir" >> "$OUTPUT_SHEET"
  else
    log "-> WARNING: No FASTQ files found for sample $SAMPLE_ID. Skipping."
  fi
done

log "-----------------------------------------------------"
log "✅ All Done!"
log "Samplesheet successfully created: $OUTPUT_SHEET"
log "Your organized FASTQ directories are ready at: $LINK_ROOT"
log "You can now use '$OUTPUT_SHEET' as the MANIFEST input for your cellranger script."
