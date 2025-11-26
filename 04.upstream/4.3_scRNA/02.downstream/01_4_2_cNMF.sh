#!/bin/bash

DATA_DIR="/mnt/18T/chibao/gliomas/data/upstream/scRNA/official/integrated_v5_optimized/adult/cNMF_myeloid"
RUN_NAME="Glioma_Adult_Myeloid"
CHOSEN_K=14 

cd "$DATA_DIR" || exit 1

# Run 1: Standard Threshold
echo "Running Consensus (Standard 0.1)..."
cnmf consensus \
    --output-dir "$DATA_DIR" \
    --name "$RUN_NAME" \
    --components "$CHOSEN_K" \
    --local-density-threshold 0.1 \
    --show-clustering

# Run 2: Strict Threshold (Closer to Paper's 0.015)
echo "Running Consensus (Strict 0.02)..."
cnmf consensus \
    --output-dir "$DATA_DIR" \
    --name "$RUN_NAME" \
    --components "$CHOSEN_K" \
    --local-density-threshold 0.02 \
    --show-clustering

echo "Done! Check the two .png clustering plots to see which is cleaner."