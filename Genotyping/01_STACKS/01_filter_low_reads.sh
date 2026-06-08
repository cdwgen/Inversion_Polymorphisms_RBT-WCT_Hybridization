#!/bin/bash

# Create output dir
mkdir -p Less_Than_10K_Reads

# Function to count reads and move poor samples
process_bam() {
    bam="$1"
    count=$(samtools view -c "$bam")
    if [ "$count" -lt 10000 ]; then
        echo "$bam has $count reads — moving to Less_Than_10K_Reads"
        mv "$bam" Less_Than_10K_Reads/
    fi
}
export -f process_bam

# Run in parallel
parallel --bar --jobs 10 process_bam ::: *.bam