#!/bin/bash
set -eo pipefail

# Log file for corrupted BAMs
corrupt_log="corrupted_bams.log"
> "$corrupt_log"

# Filter, sort, and write BAMs
filter_and_write() {
    input_file="$1"
    output_directory="$2"

    base_name=$(basename "$input_file" .bam)
    output_file="${output_directory}/${base_name}.filtered.bam"

    echo "Filtering, sorting $input_file -> $output_file"
    if ! samtools view -b -h -f2 -F2308 -q 20 "$input_file" \
        | samtools sort -o "$output_file"; then
        echo "Initial filter failed for $input_file"
        return 1
    fi

    # Index
    samtools index "$output_file"
}

# Validate BAM files
validate_bam() {
    bam_file="$1"
    if ! samtools view "$bam_file" > /dev/null 2>&1; then
        echo "$(date +%F_%T) $bam_file is corrupted" >> "$corrupt_log"
        return 1
    else
        return 0
    fi
}

# Retry failed filter once
retry_filter_if_corrupt() {
    filtered_bam="$1"
    input_directory="$2"
    base_name=$(basename "$filtered_bam" .filtered.bam)
    input_file="${input_directory}/${base_name}.bam"

    if ! validate_bam "$filtered_bam"; then
        echo "Retrying filter for $input_file..."
        if samtools view -b -h -f2 -F2308 -q 20 "$input_file" \
            | samtools sort -o "$filtered_bam"; then
            validate_bam "$filtered_bam" && samtools index "$filtered_bam" && echo "Retry succeeded for $filtered_bam"
        else
            echo "$(date +%F_%T) Retry failed for $filtered_bam" >> "$corrupt_log"
        fi
    fi
}

# Directories
input_directory="."  # Location of raw BAMs
output_directory="/mnt/ceph/well0766/UMT_Redband/02_filtered_sorted_BAMS/02_filter"
mkdir -p "$output_directory"

# Functions and variables for parallel
export -f filter_and_write
export -f validate_bam
export -f retry_filter_if_corrupt
export corrupt_log

# Filter + sort + index
echo "Filtering, sorting, and indexing BAM files..."
parallel --bar -j 30 filter_and_write {} "$output_directory" ::: "$input_directory"/*.bam

# Validate + Retry if corrupt
echo "Validating filtered BAM files..."
parallel --bar -j 30 'validate_bam {} || retry_filter_if_corrupt {} "'"$input_directory"'"' ::: "$output_directory"/*.bam

echo "Filtering, sorting, validation, and indexing complete."
echo "Any corrupted files have been logged in $corrupt_log"
