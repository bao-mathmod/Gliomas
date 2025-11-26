#!/bin/bash

# ======================= CONFIG ==========================
DATA_DIR="/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/cNMF_myeloid"
RUN_NAME="Glioma_Adult_Myeloid"
K_RANGE="2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18" 
WORKERS=60 
# =========================================================

# Start global timer
SCRIPT_START=$SECONDS

# Helper function to format time (HH:MM:SS)
format_time() {
    local T=$1
    local H=$((T/3600))
    local M=$(( (T/60)%60 ))
    local S=$(( T%60 ))
    printf "%02d:%02d:%02d" $H $M $S
}

# Helper function for the progress bar
draw_progress_bar() {
    local current=$1
    local total=$2
    local width=50
    
    # Calculate percentage
    local percent=$(( 100 * current / total ))
    local num_chars=$(( width * current / total ))
    
    # Create bar string
    local bar=""
    for ((j=0; j<num_chars; j++)); do bar+="#"; done
    for ((j=num_chars; j<width; j++)); do bar+="."; done
    
    # Print bar with carriage return (\r) to overwrite the line
    printf "\r[%-50s] %3d%% (%d/%d Workers Finished) - Time Elapsed: %s" "$bar" "$percent" "$current" "$total" "$(format_time $((SECONDS - step_start_time)))"
}

cd "$DATA_DIR" || exit 1

echo "==================================================="
echo "Starting cNMF Pipeline | $(date)"
echo "==================================================="

# ---------------- STEP 1: PREPARE ----------------
echo ""
echo "[Step 1/3] Running Prepare..."
step_start_time=$SECONDS

cnmf prepare \
      --output-dir "$DATA_DIR" \
      --name "$RUN_NAME" \
      -c "$DATA_DIR/matrix.mtx" \
      --genes-file "$DATA_DIR/variable_genes.txt" \
      --n-iter 100 \
      --seed 123 \
      -k $K_RANGE

step_duration=$(( SECONDS - step_start_time ))
echo "✅ Prepare completed in $(format_time $step_duration)"

# ---------------- STEP 2: FACTORIZE ----------------
echo ""
echo "[Step 2/3] Running Factorize (Parallelized on $WORKERS workers)..."
step_start_time=$SECONDS

# Launch workers in background
for i in $(seq 0 $(($WORKERS - 1))); do
    cnmf factorize \
        --output-dir "$DATA_DIR" \
        --name "$RUN_NAME" \
        --worker-index "$i" \
        --total-workers "$WORKERS" &
done

# --- PROGRESS BAR LOOP ---
# Instead of a simple 'wait', we loop until all jobs are done
total_jobs=$WORKERS

while true; do
    # Count how many background jobs are currently running
    running_jobs=$(jobs -r | wc -l)
    
    # Calculate finished jobs
    finished_jobs=$(( total_jobs - running_jobs ))
    
    # Update the progress bar
    draw_progress_bar "$finished_jobs" "$total_jobs"
    
    # Break loop if no jobs are running
    if [ "$running_jobs" -eq 0 ]; then
        break
    fi
    
    # Wait a bit before updating again
    sleep 2
done
echo "" # Move to new line after progress bar finishes

step_duration=$(( SECONDS - step_start_time ))
echo "✅ Factorize completed in $(format_time $step_duration)"

# ---------------- STEP 3: COMBINE & PLOT ----------------
echo ""
echo "[Step 3/3] Combining results and plotting..."
step_start_time=$SECONDS

echo "   > Running Combine..."
cnmf combine --output-dir "$DATA_DIR" --name "$RUN_NAME"

echo "   > Running K-Selection Plot..."
cnmf k_selection_plot --output-dir "$DATA_DIR" --name "$RUN_NAME"

step_duration=$(( SECONDS - step_start_time ))
echo "✅ Combine & Plot completed in $(format_time $step_duration)"

# ---------------- FINAL SUMMARY ----------------
total_duration=$(( SECONDS - SCRIPT_START ))
echo ""
echo "==================================================="
echo "🎉 Pipeline Finished Successfully!"
echo "Total Runtime: $(format_time $total_duration)"
echo "Check output at: $DATA_DIR/$RUN_NAME"
echo "==================================================="