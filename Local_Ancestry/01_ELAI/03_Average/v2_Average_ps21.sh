#!/bin/bash

# === USER VARIABLES ===
base_dir="/mnt/ceph/well0766/UMT_Redband/06_local_ancestry/01_ELAI/01_all_chrom_unphased/02_Results"
output_dir="${base_dir}/averaged_ps21"
generations=(10 20 30 40 50)
parallel_jobs=10

mkdir -p "$output_dir"

# === DETECT CHROMOSOMES ===
chrom_list=($(ls "${base_dir}/${generations[0]}_Generations/" | grep '^Omy'))
echo "Detected chromosomes: ${chrom_list[*]}"

# === FUNCTION TO AVERAGE ONE CHROMOSOME ===
average_chrom() {
  chrom="$1"
  base_dir="$2"
  output_dir="$3"
  generations=("${@:4}")

  # Collect files
  files=()
  for gen in "${generations[@]}"; do
    f="${base_dir}/${gen}_Generations/${chrom}/output/${chrom}.ps21.txt"
    if [[ ! -f "$f" ]]; then
      echo "⚠️ Missing file $f → skipping $chrom"
      return
    fi
    files+=("$f")
  done

  out="${output_dir}/${chrom}.ps21.txt"

  # Original working paste+awk logic
  paste "${files[@]}" | awk -v n="${#files[@]}" '
    BEGIN {OFS=" "}
    NR==1 {
      # Check for header
      is_header=0
      for (i=1;i<=NF;i++) { if ($i ~ /[^0-9\.\-eE]/) {is_header=1; break} }
      if (is_header) {print; next}
    }
    {
      for (i = 1; i <= NF/n; i++) {
        sum = 0
        for (j = 0; j < n; j++) sum += $(i + j*(NF/n))
        printf "%.6f ", sum/n
      }
      print ""
    }
  ' > "$out"

  echo "✅ Averaged: $chrom → $out"
}

export -f average_chrom
export base_dir output_dir

# === RUN IN PARALLEL WITH PROGRESS BAR ===
# GNU parallel --bar shows live progress for all chromosomes
parallel --bar -j "$parallel_jobs" average_chrom {} "$base_dir" "$output_dir" "${generations[*]}" ::: "${chrom_list[@]}"

echo "All done."
