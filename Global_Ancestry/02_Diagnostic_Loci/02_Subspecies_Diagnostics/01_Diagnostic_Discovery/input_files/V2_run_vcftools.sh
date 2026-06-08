#!/bin/bash
set -euo pipefail

VCF="/mnt/ceph/well0766/UMT_Redband/04_vcf/01_STACKS/01_global_ancestry/02_Concat_VCF/03_remove_dups/all_RADCap_10x.filtered.noCCT.nodups.recode.vcf.gz"
REPS=30
JOBS=30  # Number of parallel jobs
LOGFILE="vcftools_all_replicates.log"

module load vcftools

run_replicate() {
    i=$1
    echo "=== Replicate $i ==="

    # Pairwise FST
    vcftools --gzvcf "$VCF" --weir-fst-pop subsample_IRT_${i}.txt --weir-fst-pop subsample_CRT_${i}.txt --out fst_IRTvsCRT_${i}
    vcftools --gzvcf "$VCF" --weir-fst-pop subsample_IRT_${i}.txt --weir-fst-pop subsample_WCT_${i}.txt --out fst_IRTvsWCT_${i}
    vcftools --gzvcf "$VCF" --weir-fst-pop subsample_IRT_${i}.txt --weir-fst-pop subsample_YCT_${i}.txt --out fst_IRTvsYCT_${i}
    vcftools --gzvcf "$VCF" --weir-fst-pop subsample_CRT_${i}.txt --weir-fst-pop subsample_WCT_${i}.txt --out fst_CRTvsWCT_${i}
    vcftools --gzvcf "$VCF" --weir-fst-pop subsample_CRT_${i}.txt --weir-fst-pop subsample_YCT_${i}.txt --out fst_CRTvsYCT_${i}
    vcftools --gzvcf "$VCF" --weir-fst-pop subsample_WCT_${i}.txt --weir-fst-pop subsample_YCT_${i}.txt --out fst_WCTvsYCT_${i}

    # Allele frequencies
    vcftools --gzvcf "$VCF" --keep subsample_IRT_${i}.txt --freq --out freq_IRT_${i}
    vcftools --gzvcf "$VCF" --keep subsample_CRT_${i}.txt --freq --out freq_CRT_${i}
    vcftools --gzvcf "$VCF" --keep subsample_WCT_${i}.txt --freq --out freq_WCT_${i}
    vcftools --gzvcf "$VCF" --keep subsample_YCT_${i}.txt --freq --out freq_YCT_${i}

    echo "=== Replicate $i completed ==="
}

export -f run_replicate
export VCF

# Run all replicates in parallel, redirect all output to a single log file
seq 1 $REPS | parallel -j $JOBS run_replicate {} &> "$LOGFILE"
