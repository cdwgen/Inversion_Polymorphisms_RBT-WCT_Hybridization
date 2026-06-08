#!/bin/bash
set -euo pipefail
shopt -s nullglob

input_dir="."
output_dir="."
mkdir -p "$output_dir"

# collect groups from filenames
groups=()
for f in "$input_dir"/freq_*.frq; do
    [[ -e "$f" ]] || continue
    base=$(basename "$f")
    if [[ $base =~ ^freq_([A-Z]+)_[0-9]+\.frq$ ]]; then
        groups+=("${BASH_REMATCH[1]}")
    fi
done

# make unique
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

summary_out="$output_dir/freq_summary_sd.tsv"
echo -e "GROUP\tMEAN_SD_REF_FREQ\tMEAN_SD_ALT_FREQ" > "$summary_out"

for group in "${uniq_groups[@]}"; do
    echo "Processing group: $group"

    files=( "$input_dir"/freq_"$group"_*.frq )
    if [ ${#files[@]} -eq 0 ]; then
        echo "  No files found for group $group"
        continue
    fi

    # Compute SD per locus and then summarize across loci
    awk -F '\t' -v first="${files[0]}" -v group="$group" '
    BEGIN { OFS="\t" }
    FILENAME == first && FNR > 1 {
        key = $1 "\t" $2 "\t" $5 "\t" $7
        if (!(key in seen)) { seen[key] = 1; order[++n] = key }
    }
    FNR > 1 {
        key = $1 "\t" $2 "\t" $5 "\t" $7
        ref = ($6 + 0)
        alt = ($8 + 0)
        sum_ref[key] += ref
        sum_alt[key] += alt
        sumsq_ref[key] += (ref * ref)
        sumsq_alt[key] += (alt * alt)
        cnt[key] += 1
    }
    END {
        total_sd_ref = 0
        total_sd_alt = 0
        locus_count = 0
        for (i = 1; i <= n; i++) {
            k = order[i]
            if (cnt[k] > 1) {
                mean_ref = sum_ref[k] / cnt[k]
                mean_alt = sum_alt[k] / cnt[k]
                var_ref = (sumsq_ref[k] / cnt[k]) - (mean_ref * mean_ref)
                var_alt = (sumsq_alt[k] / cnt[k]) - (mean_alt * mean_alt)
                sd_ref = (var_ref > 0 ? sqrt(var_ref) : 0)
                sd_alt = (var_alt > 0 ? sqrt(var_alt) : 0)
                total_sd_ref += sd_ref
                total_sd_alt += sd_alt
                locus_count++
            }
        }
        if (locus_count > 0) {
            mean_sd_ref = total_sd_ref / locus_count
            mean_sd_alt = total_sd_alt / locus_count
            printf "%s\t%.6f\t%.6f\n", group, mean_sd_ref, mean_sd_alt
        } else {
            printf "%s\tNA\tNA\n", group
        }
    }' "${files[@]}" >> "$summary_out"

    echo "  Added summary for $group"
done

echo "Summary written to $summary_out"
