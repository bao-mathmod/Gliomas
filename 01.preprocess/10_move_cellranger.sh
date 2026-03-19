#!/bin/bash

# ==========================================
# TOGGLE DRY RUN HERE
# true  = Only print what would happen (Safe)
# false = Actually move the files
# ==========================================
DRY_RUN=false

# Define your source and destination paths
SOURCE_BASE="/mnt/18T/chibao/gliomas/data_official/00_raw_data_adult/1_raw_cell_ranger/new_cohort/nextflow/cellranger/count"
DEST_BASE="/mnt/18T/chibao/gliomas/data_official/00_raw_data_adult/1_raw_cell_ranger/cohort_official/PRJNA647809"

if [ "$DRY_RUN" = true ]; then
    echo "=========================================================="
    echo " DRY RUN MODE ENABLED: No files will be moved or created. "
    echo "=========================================================="
else
    # Create the destination base folder if not in dry-run mode
    mkdir -p "$DEST_BASE"
fi

# Loop through every SAMN folder in the source directory
for samn_path in "$SOURCE_BASE"/SAMN*; do
    
    # Check if it is a directory (ignores stray files)
    if [ -d "$samn_path" ]; then
        
        # Get just the folder name (e.g., "SAMN15590398")
        samn_id=$(basename "$samn_path")
        
        # Define the exact path to the filtered matrix folder
        target_matrix="$samn_path/outs/filtered_feature_bc_matrix"
        
        # Check if the filtered matrix folder actually exists
        if [ -d "$target_matrix" ]; then
            
            if [ "$DRY_RUN" = true ]; then
                # Just print what it would do
                echo "[DRY RUN] Would create folder: $DEST_BASE/$samn_id"
                echo "[DRY RUN] Would move: $target_matrix"
                echo "                 -->  $DEST_BASE/$samn_id/"
                echo "------------------------------------------------------"
            else
                # Actually do the work
                echo "Processing $samn_id..."
                mkdir -p "$DEST_BASE/$samn_id"
                mv "$target_matrix" "$DEST_BASE/$samn_id/"
                echo "Successfully moved to $DEST_BASE/$samn_id/"
            fi
            
        else
            echo "Warning: Filtered matrix not found for $samn_id. Skipping..."
        fi
    fi
done

if [ "$DRY_RUN" = true ]; then
    echo "Dry run complete. If everything looks correct, change DRY_RUN=false in the script and run it again."
else
    echo "All transfers complete!"
fi