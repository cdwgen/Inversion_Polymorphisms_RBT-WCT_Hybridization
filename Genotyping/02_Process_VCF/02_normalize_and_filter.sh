#!/bin/bash
set -euo pipefail

INPUT_VCF="../data/05_merged_vcf/All_Koot_Ath_NoMO.snps.vcf.gz"
REFERENCE="../data/reference/norm_IRT_061722.fa"
OUTPUT_PREFIX="../data/06_filtered_vcf/All_Koot_Ath_Norm_10X"
TEMP_NORM="normalized.vcf"

mkdir -p "../data/06_filtered_vcf"

# Normalize VCF
echo "Normalizing VCF against reference genome..."
bcftools norm -c s -f "$REFERENCE" "$INPUT_VCF" -Ov -o "$TEMP_NORM"

# Apply VCFtools Filters
echo "Applying genotype and locus-level filters..."
vcftools --vcf "$TEMP_NORM" \
    --minGQ 30 \
    --minDP 10 \
    --mac 3 \
    --max-missing 0.8 \
    --min-alleles 2 \
    --max-alleles 2 \
    --remove-indels \
    --recode --recode-INFO-all \
    --out "${OUTPUT_PREFIX}.filtered.snps"

rm "$TEMP_NORM"
echo "Filtering complete."