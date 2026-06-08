#!/bin/bash

# Set input and output directories
input_dir="."     # <-- change this
output_dir="reformatted_freq"   # <-- change this

# Make sure output directory exists
mkdir -p "$output_dir"

# Loop over all .frq files in the input directory
for file in "$input_dir"/*.frq; do
    # Extract base filename
    base=$(basename "$file")
    outfile="$output_dir/$base"

    # Reformat with awk
    awk '
    BEGIN {
        OFS = "\t";
        print "CHROM", "POS", "N_ALLELES", "N_CHR", "REF_ALLELE", "REF_FREQ", "ALT_ALLELE", "ALT_FREQ"
    }
    NR > 1 {
        split($5, ref, ":");
        split($6, alt, ":");
        print $1, $2, $3, $4, ref[1], ref[2], alt[1], alt[2];
    }
    ' "$file" > "$outfile"

    echo "Reformatted: $base → $outfile"
done

