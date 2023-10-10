#!/bin/bash

#SBATCH --job-name=ADMIXTURE
#SBATCH --output=ADMIXTURE_2.out
##SBATCH --error=SNP_filtering.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=36:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu

##Remove Duplicate Variants

module load plink

#plink2 --bfile ALLCHR.phase3_v5.shapeit2_mvncall_integrated_int_subset --exclude dupIDs.txt  --make-bed --out ALLCHR.phase3_v5.shapeit2_mvncall_integrated_int_noDups_subset

##Filter out SNPs that are in LD with each other. scanning windows of 50kb with a 20kb offset, and pruning so that no pairs in each window had an R2>0.8. 
##Also filtering out variants with a MAF of <1% and variants which are missing from >25% of samples. 

#plink2 --bfile ALLCHR.phase3_v5.shapeit2_mvncall_integrated_int_noDups_subset --geno 0.25 --maf 0.01 --indep-pairwise 50 20 0.8

#plink2 --bfile ALLCHR.phase3_v5.shapeit2_mvncall_integrated_int_noDups_subset --exclude plink2.prune.out --make-bed --out ALLCHR.phase3_v5.shapeit2_mvncall_integrated_int_allFilters_subset

## SEND THESE OUTPUT FILES TO BE THE INPUT OF pca.sh

##Test for fine scale population structure using ADMIXTURE
for K in {6..6}
do 
/Genomics/grid/users/alea/programs/admixture_linux-1.3.0/admixture --cv ALLCHR.phase3_v5.shapeit2_mvncall_integrated_int_allFilters_subset.bed $K --seed=$RANDOM | tee ALLCHR.phase3_v5.shapeit2_mvncall_integrated_ADMIX.${K}.$RANDOM.out  ;
done


####NOTE TO SELF: NEED TO REPLACE CHR X NAME IN BIM FILE WITH "23" BECAUSE ADMIXTURE CAN ONLY READ CHR NAMES AS INTEGER NUMBERS. Can do this either by changing chr names in vcf file with
####	          bcf tools and including --output-chr 25 in all --make-bed plink commands OR by simple using sed to replace X with 23 in all newly produced bim files####
