### I followed these steps to create a BED file and run ADMIXTURE
## Now to run ADMIXTURE
## Need a bed file

## Rename chromosomes first

# a chromosome map is needed, use this one:

Omy01 1
Omy02 2
Omy03 3
Omy04 4
Omy05 5
Omy06 6
Omy07 7
Omy08 8
Omy09 9
Omy10 10
Omy11 11
Omy12 12
Omy13 13
Omy14 14
Omy15 15
Omy16 16
Omy17 17
Omy18 18
Omy19 19
Omy20 20
Omy21 21
Omy22 22
Omy23 23
Omy24 24
Omy25 25
Omy26 26
Omy27 27
Omy28 28
OmyY 29
Omy30 30
Omy31 31
Omy32 32

## use the map to rename the chromosomes in your vcf:

module load bcftools

bcftools annotate --rename-chrs chr_map.txt -o RADCap.reannotate.vcf All_Koot_Ath_Norm_10X.filtered.no_high_missing.snps.recode.vcf

## Now create the plink bed file:

module load plink

plink2 --vcf RADCap.reannotate.vcf --chr-set 32 --allow-extra-chr --make-bed --out All_KootAth_10X

# get sample Ids

grep "^#CHROM" RADCap.reannotate.vcf | cut -f10- | tr '\t' '\n' > sample_ids.txt