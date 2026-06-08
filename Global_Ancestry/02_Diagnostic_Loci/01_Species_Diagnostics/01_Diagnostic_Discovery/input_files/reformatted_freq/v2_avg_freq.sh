#!/bin/bash
set -euo pipefail
shopt -s nullglob

input_dir="."
output_dir="."
mkdir -p "$output_dir"

# collect groups from filenames (preserves order of discovery)
groups=()
for f in "$input_dir"/freq_*.frq; do
    [[ -e "$f" ]] || continue
    base=$(basename "$f")
    if [[ $base =~ ^freq_([A-Z]+)_[0-9]+\.frq$ ]]; then
        groups+=("${BASH_REMATCH[1]}")
    fi
done

# make unique while preserving order
declare -A _seen
uniq_groups=()
for g in "${groups[@]}"; do
    if [[ -z "${_seen[$g]:-}" ]]; then
        _seen[$g]=1
        uniq_groups+=("$g")
    fi
done

if [ ${#uniq_groups[@]} -eq 0 ]; then
    echo "No freq_*.frq files found in $input_dir"
    exit 0
fi

for group in "${uniq_groups[@]}"; do
    echo "Processing group: $group"

    files=( "$input_dir"/freq_"$group"_*.frq )
    # check for real files
    if [ ${#files[@]} -eq 0 ] || [ ! -e "${files[0]}" ]; then
        echo "  No files found for group $group"
        continue
    fi

    out="$output_dir/${group}_averaged.frq"

    # Use awk to:
    #  - capture key order from the first file
    #  - sum ref/alt frequencies and counts across all files
    #  - print means in original order
    awk -F '\t' -v first="${files[0]}" '
    BEGIN { OFS = "\t" }
    # record order from the first file (skip header)
    FILENAME == first && FNR > 1 {
        key = $1 "\t" $2 "\t" $5 "\t" $7
        if (!(key in seen)) { seen[key] = 1; order[++n] = key }
    }
    # accumulate from every file (skip header lines)
    FNR > 1 {
        key = $1 "\t" $2 "\t" $5 "\t" $7
        sum_ref[key] += ($6 + 0)
        sum_alt[key] += ($8 + 0)
        cnt[key] += 1
    }
    END {
        print "CHROM","POS","REF_ALLELE","ALT_ALLELE","MEAN_REF_FREQ","MEAN_ALT_FREQ"
        for (i = 1; i <= n; i++) {
            k = order[i]
            if (cnt[k] > 0) {
                mean_ref = sum_ref[k] / cnt[k]
                mean_alt = sum_alt[k] / cnt[k]
                # print key (which already contains 4 tab-separated fields) + the two means
                printf "%s\t%.6f\t%.6f\n", k, mean_ref, mean_alt
            } else {
                # should not normally happen, but guard anyway
                printf "%s\tNA\tNA\n", k
            }
        }
    }
    ' "${files[@]}" > "$out"

    echo "  Written $out"
done
