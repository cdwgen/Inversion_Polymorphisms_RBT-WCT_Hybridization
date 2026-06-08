#!/bin/bash
set -euo pipefail

# === CONFIG ===
GROUP1="IRT.txt"
GROUP2="CRT.txt"
GROUP3="WCT.txt"
GROUP4="YCT.txt"   
MIN_SIZE=51        # Smallest group size (YCT)
REPS=30

echo "Generating subsample lists..."

for i in $(seq 1 $REPS); do
    echo "Replicate $i"
    
    shuf "$GROUP1" | head -n "$MIN_SIZE" > subsample_IRT_${i}.txt
    shuf "$GROUP2" | head -n "$MIN_SIZE" > subsample_CRT_${i}.txt
    shuf "$GROUP3" | head -n "$MIN_SIZE" > subsample_WCT_${i}.txt
    cp "$GROUP4" subsample_YCT_${i}.txt
done

echo "All subsample lists generated."
