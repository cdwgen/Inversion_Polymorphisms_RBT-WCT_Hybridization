#!/bin/bash

# Set input VCF and positions directory
VCF="/mnt/ceph/well0766/UMT_Redband/04_vcf/01_STACKS/01_global_ancestry/02_Concat_VCF/02_remove_CCT/all_RADCap_10x.filtered.noCCT.recode.vcf.gz"
POSITIONS_DIR="./sets_rapture"

# Load necessary modules
module load vcftools

# List of diagnostic SNP sets
SETS=("Inv1" "Inv2")

for SET in "${SETS[@]}"; do
  echo "Processing $SET..."

  # Extract only CHROM and POS to a temp file (no header assumed)
  cut -f1,2 "$POSITIONS_DIR/${SET}.txt" > "${SET}_positions.temp"

  # Filter VCF based on positions file
  vcftools --gzvcf "$VCF" \
           --positions "${SET}_positions.temp" \
           --out "${SET}_all_common_input.vcf" \
           --recode --recode-INFO-all

  # Extract genotype (GT) information
  vcftools --vcf "${SET}_all_common_input.vcf.recode.vcf" \
           --extract-FORMAT-info GT \
           --out "${SET}_all_common_input"

  # Remove temporary file
  rm "${SET}_positions.temp"
done

echo "Done!"
