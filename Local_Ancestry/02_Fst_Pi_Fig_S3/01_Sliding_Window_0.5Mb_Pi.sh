#!/bin/bash

# Exit on error, undefined variable, or pipe failure
set -euo pipefail
set -x

# Load modules
module load vcftools bcftools htslib

# ==== CONFIG ====
VCF="/mnt/ceph/well0766/UMT_Redband/04_vcf/01_STACKS_2025/01_global_ancestry/02_Concat_VCF/02_remove_CCT/all_RADCap_10x.filtered.noCCT.recode.vcf.gz"

# FIXED: removed duplicate IRT
POPLIST=(CRT_AA.tsv CRT_RR.tsv IRT.tsv WCT.tsv YCT.tsv)

# Parameters for window analysis
WINDOW_SIZE=500000
WINDOW_STEP=500000

# ==== STEP 1: Create sample list for VCF subset ====
echo "Gathering all sample names from population files..."
cat "${POPLIST[@]}" | sort | uniq > all_samples.txt

# ==== STEP 2: Subset VCF ====
echo "Subsetting VCF to relevant samples..."
vcftools --gzvcf "$VCF" \
         --keep all_samples.txt \
         --recode --recode-INFO-all \
         --out subset

# Optional but recommended: compress for faster downstream steps
bgzip -f subset.recode.vcf
tabix -p vcf subset.recode.vcf.gz

VCF_SUB="subset.recode.vcf.gz"

# ==== STEP 3: Compute nucleotide diversity (π) ====
echo "Computing nucleotide diversity (π) in ${WINDOW_SIZE}bp windows..."
for popfile in "${POPLIST[@]}"; do
    popname=$(basename "$popfile" .tsv)
    echo " - π for $popname"
    
    vcftools --gzvcf "$VCF_SUB" \
             --keep "$popfile" \
             --window-pi "$WINDOW_SIZE" \
             --window-pi-step "$WINDOW_STEP" \
             --out "pi_${popname}"
done

# ==== STEP 4: Compute pairwise FST ====
echo "Computing pairwise FST (windowed + genome-wide)..."

for ((i=0; i<${#POPLIST[@]}-1; i++)); do
  for ((j=i+1; j<${#POPLIST[@]}; j++)); do

    p1="${POPLIST[$i]}"
    p2="${POPLIST[$j]}"

    name1=$(basename "$p1" .tsv)
    name2=$(basename "$p2" .tsv)

    outname="fst_${name1}_vs_${name2}"

    echo " - $name1 vs $name2"

    # ---- (A) WINDOWED FST (for plots) ----
    vcftools --gzvcf "$VCF_SUB" \
             --weir-fst-pop "$p1" \
             --weir-fst-pop "$p2" \
             --fst-window-size "$WINDOW_SIZE" \
             --fst-window-step "$WINDOW_STEP" \
             --out "${outname}"

    # ---- (B) PER-SITE FST (for proper summaries) ----
    vcftools --gzvcf "$VCF_SUB" \
             --weir-fst-pop "$p1" \
             --weir-fst-pop "$p2" \
             --out "${outname}_site"

    # ---- (C) “GENOME-WIDE” APPROX (very large window) ----
    vcftools --gzvcf "$VCF_SUB" \
             --weir-fst-pop "$p1" \
             --weir-fst-pop "$p2" \
             --fst-window-size 1000000000 \
             --out "${outname}_global"

  done
done

echo "All analyses complete!"
