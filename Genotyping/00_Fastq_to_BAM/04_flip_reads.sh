#!/bin/bash

set -euo pipefail

cd "$(dirname "$0")/.." || exit 1

# Make sure barcode files are present and the flip/trim script is accesible!
INPUT_DIR="03_trim_reads"
OUTPUT_DIR="04_flip_reads"
BARCODE_DIR="01.RawData/barcodes"
FLIP_SCRIPT="scripts/bRAD_flip_trim.pl"

mkdir -p "$OUTPUT_DIR"

export INPUT_DIR
export OUTPUT_DIR
export BARCODE_DIR
export FLIP_SCRIPT

# Flip reads for each pair for each lane
process_flip() {

    base_name=$1

    echo "[$(date)] Flipping ${base_name}"

    perl "$FLIP_SCRIPT" \
        "${BARCODE_DIR}/${base_name}_barcodes.txt" \
        <(pigz -dc "${INPUT_DIR}/${base_name}_paired.1.fq.gz") \
        <(pigz -dc "${INPUT_DIR}/${base_name}_paired.2.fq.gz") \
        "${OUTPUT_DIR}/${base_name}_flipped.1.fq" \
        "${OUTPUT_DIR}/${base_name}_flipped.2.fq"

    pigz -p 4 "${OUTPUT_DIR}/${base_name}_flipped.1.fq"
    pigz -p 4 "${OUTPUT_DIR}/${base_name}_flipped.2.fq"

    echo "[$(date)] Finished ${base_name}"

}

export -f process_flip

find "$INPUT_DIR" -name "*_paired.1.fq.gz" \
| sed 's/_paired\.1\.fq\.gz//' \
| xargs -n1 basename \
| parallel --bar -j 2 process_flip