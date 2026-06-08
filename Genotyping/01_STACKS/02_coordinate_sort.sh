#!/bin/bash

# Make sure output directory exists
output_dir="/mnt/ceph/well0766/UMT_Redband/02_filtered_sorted_BAMS/01_sort"
mkdir -p "$output_dir"

export output_dir

# Function to sort BAMs
process_bam() {
    bam_file=$1
    base_name=$(basename "$bam_file" | sed 's/\..*//')
    sorted_bam_file="${output_dir}/${base_name}.pe.sorted.bam"

    echo "Sorting $bam_file by coordinate..."
    samtools sort "$bam_file" -o "$sorted_bam_file"
    echo "Coordinate-sorted $bam_file -> $sorted_bam_file"
}

export -f process_bam

# Find BAMs, sort in parallel
find . -name "*.bam" | parallel --bar -j 25 process_bam