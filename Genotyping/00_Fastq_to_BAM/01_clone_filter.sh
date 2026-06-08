#!/bin/bash

module load stacks/2.5 #using 2.5 to match Steve using 2.55

cd "$(dirname "$0")/.." || exit

RAW_DIR="01.RawData"
OUTPUT_BASE="01_clone_filter"

export OUTPUT_BASE

# Function for clone filter per lane
process_lane() {
    R1=$1
    PLATE_DIR=$(dirname "$R1")
    PLATE_NAME=$(basename "$PLATE_DIR")
    OUT_DIR="${OUTPUT_BASE}/${PLATE_NAME}"
    
    # Extract the lane base name (e.g., Lane1_1.fq.gz becomes Lane1)
    LANE_BASE=$(basename "$R1" _1.fq.gz)
    R2="${PLATE_DIR}/${LANE_BASE}_2.fq.gz"
    LOG_FILE="${OUT_DIR}/${LANE_BASE}.log"

    # Check if this specific lane is already successfully done
    if [ -f "$LOG_FILE" ] && grep -q "Calculating the distribution" "$LOG_FILE"; then
        echo "  -> $LANE_BASE in $PLATE_NAME already processed. Skipping..."
        return 0
    fi

    echo "  -> Running clone_filter on $LANE_BASE in $PLATE_NAME..."
    clone_filter -1 "$R1" -2 "$R2" -i gzfastq -o "$OUT_DIR" 2> "$LOG_FILE"
    
    # Capture the exit status in case of crashes
    status=$?
    if [ $status -ne 0 ]; then
        echo "ERROR: clone_filter failed on $LANE_BASE in $PLATE_NAME."
        return $status
    fi
}
export -f process_lane

# Loop through each plate directory
for PLATE_DIR in "$RAW_DIR"/IRT_Rap_p*; do
    [ -d "$PLATE_DIR" ] || continue 
    
    PLATE_NAME=$(basename "$PLATE_DIR")
    OUT_DIR="${OUTPUT_BASE}/${PLATE_NAME}"
    mkdir -p "$OUT_DIR"
    
    echo "========================================"
    echo "Processing Plate: $PLATE_NAME (3 lanes at a time)"
    echo "========================================"

    # Pass all R1 files for this plate to parallel, running 3 at once
    # --halt now,fail=1 ensures if the RAM maxes out and a lane crashes, the script stops immediately
    ls "$PLATE_DIR"/*_1.fq.gz | parallel --halt now,fail=1 -j 3 process_lane
    
    # If parallel halted due to an error, kill the whole pipeline
    if [ $? -ne 0 ]; then
        echo "Pipeline halted due to an error in $PLATE_NAME."
        exit 1
    fi

    echo "  -> Generating clone rate report for $PLATE_NAME..."
    grep "clone reads" "$OUT_DIR"/*.log > "${OUT_DIR}/${PLATE_NAME}_clone_rates.txt"
done

echo "========================================"
echo "All plates processed successfully."
echo "Generating global clone_rates.txt report..."
echo "========================================"

# Combine all the individual plate reports into one master report
cat "$OUTPUT_BASE"/*/*_clone_rates.txt > "${OUTPUT_BASE}/All_Plates_clone_rates.txt"

# Calculate and append global average clone rate
awk -F "," '{print $3}' "${OUTPUT_BASE}/All_Plates_clone_rates.txt" | \
    awk -F " " '{print $1}' | \
    awk '{sum+=$1} END {if (NR > 0) print "\nGlobal Average = ",sum/NR"%"}' \
    | tee -a "${OUTPUT_BASE}/All_Plates_clone_rates.txt"

echo "Done."