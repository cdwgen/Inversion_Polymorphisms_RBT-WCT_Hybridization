#!/bin/bash
LOG="vcftools_all_replicates.log"
OUT="fst_summary_table.txt"

echo -e "Replicate\tPair\tMean_Fst\tWeighted_Fst" > "$OUT"

awk '
# Capture replicate number
/^=== Replicate/ {rep=gensub(/[^0-9]/,"","g",$0)}

# Capture --out line to get pair name
/--out fst_/ {
    match($0, /--out fst_([^_]+)vs([^_]+)_[0-9]+/, arr)
    pair=arr[1]"vs"arr[2]
}

# Capture mean Fst
/Weir and Cockerham mean Fst estimate:/ {mean=$NF}

# Capture weighted Fst
/Weir and Cockerham weighted Fst estimate:/ {
    weighted=$NF
    # Output once we have both
    print rep "\t" pair "\t" mean "\t" weighted
}
' "$LOG" >> "$OUT"
