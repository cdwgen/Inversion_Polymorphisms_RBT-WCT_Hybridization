#!/bin/bash
set -euo pipefail

INPUT_VCF="../data/06_filtered_vcf/All_Koot_Ath_Norm_10X.filtered.snps.recode.vcf"
OUT_DIR="../data/06_filtered_vcf"
QC_SAMPLES_LIST="../data/metadata/CCT.txt" # N = 152 samples removed for QC

cd "$OUT_DIR" || exit

# Calculate missingness
echo "Calculating sample missingness..."
vcftools --vcf "$INPUT_VCF" --missing-indv --out sample_missingness

# Identify samples missing more than 30% of their genotype calls
awk '$5 > 0.3 {print $1}' sample_missingness.imiss > samples_over_30_missing.txt

# Remove high-miss samples
echo "Removing individuals with >30% missing genotypes..."
vcftools --vcf "$INPUT_VCF" \
    --remove samples_over_30_missing.txt \
    --recode --recode-INFO-all \
    --out All_Koot_Ath_Norm_10X.filtered.no_high_missing.snps

# Remove QC samples
echo "Removing QC/Outgroup samples..."
vcftools --vcf All_Koot_Ath_Norm_10X.filtered.no_high_missing.snps.recode.vcf \
    --remove "$QC_SAMPLES_LIST" \
    --recode --recode-INFO-all \
    --out all_RADCap_10x.filtered.final

echo "Final dataset generated: all_RADCap_10x.filtered.final.recode.vcf"