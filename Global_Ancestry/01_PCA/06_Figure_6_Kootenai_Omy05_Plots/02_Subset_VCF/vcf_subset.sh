#!/bin/bash

## These are the input files for plotting Omy05 hybrids in the Kootenai only
## samples.txt are Kootenai only samples (Gerrards retained) with no duplicates
## I subset with the following vcftools commands:

module load vcftools

vcftools --gzvcf ~/UMT_Redband/04_vcf/01_STACKS/01_global_ancestry/02_Concat_VCF/02_remove_CCT/all_RADCap_10x.filtered.noCCT.recode.vcf.gz \
         --keep samples.txt \
         --chr Omy05 \
         --from-bp 30102000 \
         --to-bp 60368999 \
         --recode --out Omy05_RBT_Inv1_subset

vcftools --gzvcf ~/UMT_Redband/04_vcf/01_STACKS/01_global_ancestry/02_Concat_VCF/02_remove_CCT/all_RADCap_10x.filtered.noCCT.recode.vcf.gz \
         --keep samples.txt \
         --chr Omy05 \
         --from-bp 60369000 \
         --to-bp 87502000 \
         --recode --out Omy05_RBT_Inv2_subset

vcftools --gzvcf ~/UMT_Redband/04_vcf/01_STACKS/01_global_ancestry/02_Concat_VCF/02_remove_CCT/all_RADCap_10x.filtered.noCCT.recode.vcf.gz \
         --keep samples.txt \
         --chr Omy05 \
         --from-bp 30102000 \
         --to-bp 87502000 \
         --recode --out Omy05_RBT_Full_subset
