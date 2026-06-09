#!/bin/bash

cd /mnt/ceph/well0766/UMT_Redband/01_data/04_IRT11_RAW_DATA || exit

INPUT_DIR="01_clone_filter"
OUTPUT_DIR="02_merged_lanes"
mkdir -p "$OUTPUT_DIR"

# Loop through plates and merge
for plate_dir in "$INPUT_DIR"/IRT_Rap_p*; do
    [ -d "$plate_dir" ] || continue

    orig_name=$(basename "$plate_dir")
    proper_name="${orig_name/IRT_Rap_p/IRT_Rapture11_plate}"

    echo "Merging filtered lanes for $orig_name into $proper_name..."

    cat "$plate_dir"/*_1.1.fq.gz > "${OUTPUT_DIR}/${proper_name}_1.fq.gz"
    cat "$plate_dir"/*_2.2.fq.gz > "${OUTPUT_DIR}/${proper_name}_2.fq.gz"
done

echo "Lane merging complete."