#!/bin/bash

indir="."  # Current directory
outdir="fst_averages"
mkdir -p "$outdir"

# Get all unique comparisons like RBTvsWCT, WCTvsYCT, etc.
comparisons=$(ls $indir/fst_*[0-9]*.weir.fst | sed -E 's|.*/fst_([A-Za-z]+vs[A-Za-z]+)_[0-9]+\.weir\.fst|\1|' | sort -u)

for comp in $comparisons; do
    echo "Processing $comp..."

    # List all replicate files for this comparison
    files=($(ls "$indir"/fst_${comp}_*.weir.fst | sort -V))
    if [ ${#files[@]} -eq 0 ]; then
        echo "  No FST files found for $comp"
        continue
    fi

    # Extract CHROM and POS (assumed same across replicates)
    awk 'NR>1 {print $1"\t"$2}' "${files[0]}" > "$outdir/${comp}_coords.tmp"

    # Extract FST values from each file
    for f in "${files[@]}"; do
        awk 'NR>1 {print $3}' "$f" > "$outdir/tmp_$(basename "$f").col"
    done

    paste "$outdir/${comp}_coords.tmp" "$outdir"/tmp_*.col > "$outdir/${comp}_all_fst.tsv"

    # Compute average and variance
    awk 'BEGIN {OFS="\t"} {
        chrom=$1; pos=$2;
        n=0; sum=0; sumsq=0;
        for (i=3; i<=NF; i++) {
            if ($i != "-nan") {
                val = $i + 0;
                sum += val;
                sumsq += val * val;
                n++;
            }
        }
        if (n == 0) {
            avg = "-nan"; var = "-nan";
        } else {
            avg = sum / n;
            var = (sumsq - sum*sum/n) / n;
        }
        print chrom, pos, avg >> "'$outdir/fst_${comp}.weir.fst'";
        print chrom, pos, var >> "'$outdir/fst_${comp}.weir.var'";
    }' "$outdir/${comp}_all_fst.tsv"

    # Clean up
    rm "$outdir"/tmp_*.col "$outdir/${comp}_coords.tmp" "$outdir/${comp}_all_fst.tsv"
done

echo "Done. See results in $outdir/"
