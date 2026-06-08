#!/bin/bash
set -eo pipefail

# Directories
bam_dir="/mnt/ceph/well0766/UMT_Redband/02_filtered_sorted_BAMS/02_filter"
output_dir="/mnt/ceph/well0766/UMT_Redband/02_filtered_sorted_BAMS/03_chromosomes"
mkdir -p "$output_dir"

# Chromosome whitelist based on norm_IRT_061722.fa (refseq Arlee assembly with renamed chromosomes)
chromosomes=("Omy01" "Omy02" "Omy03" "Omy04" "Omy05" "Omy06" "Omy07" "Omy08" "Omy09" "Omy10"
             "Omy11" "Omy12" "Omy13" "Omy14" "Omy15" "Omy16" "Omy17" "Omy18" "Omy19" "Omy20"
             "Omy21" "Omy22" "Omy23" "Omy24" "Omy25" "Omy26" "Omy27" "Omy28" "Omy30" "Omy31"
             "Omy32" "OmyY")

# Process one BAM and one chromosome
process_bam() {
    local bam_file="$1"
    local chr="$2"
    local bam_name
    bam_name=$(basename "$bam_file" .pe.sorted.filtered.bam)

    local chr_dir="$output_dir/$chr"
    mkdir -p "$chr_dir"
    local output_bam="$chr_dir/${bam_name}.bam"

    # Skip chromosome if no reads
    if [[ $(samtools view -c "$bam_file" "$chr") -gt 0 ]]; then
        # Extract, sort, and index in one step
        if ! samtools view -b "$bam_file" "$chr" | samtools sort -o "$output_bam"; then
            echo "Error processing $chr from $bam_file" >> "$output_dir/split_errors.log"
            return 1
        fi
        samtools index "$output_bam"
        echo "Processed $chr from $bam_file -> $output_bam"
    else
        echo "No reads for $chr in $bam_file" >> "$output_dir/empty_chroms.log"
    fi
}

export -f process_bam
export output_dir
export chromosomes

# Run in parallel
parallel --bar -j 32 process_bam ::: "$bam_dir"/*.pe.sorted.filtered.bam ::: "${chromosomes[@]}"

echo "All BAMs split by chromosome, sorted, and indexed."