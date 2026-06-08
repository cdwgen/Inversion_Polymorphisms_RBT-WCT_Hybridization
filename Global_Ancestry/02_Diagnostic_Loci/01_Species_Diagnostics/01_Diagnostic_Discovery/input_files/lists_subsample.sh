#!/bin/bash
set -euo pipefail

# === CONFIG ===
GROUP1="RBT.txt"
GROUP2="WCT.txt"
GROUP3="YCT.txt"
MIN_SIZE=51   # Smallest group size (YCT)
REPS=30

echo "Generating subsample lists..."

for i in $(seq 1 $REPS); do
    echo "Replicate $i"
    
    shuf "$GROUP1" | head -n "$MIN_SIZE" > subsample_RBT_${i}.txt
    shuf "$GROUP2" | head -n "$MIN_SIZE" > subsample_WCT_${i}.txt
    cp "$GROUP3" subsample_YCT_${i}.txt
done

echo "All subsample lists generated."
