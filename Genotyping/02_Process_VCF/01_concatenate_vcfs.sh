#!/bin/bash
set -euo pipefail

# Directories
STACKS_DIR="../data/04_stacks_output/populations"
VCF_DIR="../data/05_merged_vcf"
MT_DIR="${VCF_DIR}/OmyMT"

mkdir -p "$VCF_DIR"
mkdir -p "$MT_DIR"

# Copy and rename VCFs
echo "Extracting per-chromosome VCFs..."
for SUBDIR in "$STACKS_DIR"/*; do
    if [ -d "$SUBDIR" ]; then
        CHROM=$(basename "$SUBDIR")
        VCF_FILE="$SUBDIR/populations.snps.vcf"
        
        if [ -f "$VCF_FILE" ]; then
            cp "$VCF_FILE" "$VCF_DIR/${CHROM}.snps.vcf"
        fi
    fi
done

cd "$VCF_DIR" || exit

# Compress, sort, and index
echo "Compressing and sorting nuclear VCFs..."
for vcf in *.vcf; do
    bgzip "$vcf"
done

for vcf in *.vcf.gz; do
    sorted_name="${vcf%.vcf.gz}.sorted.vcf.gz"
    bcftools sort -Oz -o "$sorted_name" "$vcf"
    bcftools index -t "$sorted_name"
done

# Concatenate chromosomes
echo "Concatenating into global VCF..."
# Note that OmyY is the same as Omy29
bcftools concat \
    Omy01.snps.sorted.vcf.gz Omy02.snps.sorted.vcf.gz Omy03.snps.sorted.vcf.gz \
    Omy04.snps.sorted.vcf.gz Omy05.snps.sorted.vcf.gz Omy06.snps.sorted.vcf.gz \
    Omy07.snps.sorted.vcf.gz Omy08.snps.sorted.vcf.gz Omy09.snps.sorted.vcf.gz \
    Omy10.snps.sorted.vcf.gz Omy11.snps.sorted.vcf.gz Omy12.snps.sorted.vcf.gz \
    Omy13.snps.sorted.vcf.gz Omy14.snps.sorted.vcf.gz Omy15.snps.sorted.vcf.gz \
    Omy16.snps.sorted.vcf.gz Omy17.snps.sorted.vcf.gz Omy18.snps.sorted.vcf.gz \
    Omy19.snps.sorted.vcf.gz Omy20.snps.sorted.vcf.gz Omy21.snps.sorted.vcf.gz \
    Omy22.snps.sorted.vcf.gz Omy23.snps.sorted.vcf.gz Omy24.snps.sorted.vcf.gz \
    Omy25.snps.sorted.vcf.gz Omy26.snps.sorted.vcf.gz Omy27.snps.sorted.vcf.gz \
    Omy28.snps.sorted.vcf.gz Omy30.snps.sorted.vcf.gz Omy31.snps.sorted.vcf.gz \
    Omy32.snps.sorted.vcf.gz OmyY.snps.sorted.vcf.gz \
    -Oz -o All_Koot_Ath_NoMO.snps.vcf.gz

echo "Concatenation complete"