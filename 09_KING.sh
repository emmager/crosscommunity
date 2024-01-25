#!/bin/bash
#SBATCH --job-name=KINGinput
#SBATCH --output=KING.out
#SBATCH --error=KING.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=1:00:00
##SBATCH --array=0-1
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu

module load conda
module load bcftools

conda activate plink2

snparcherdir=/Genomics/ayroleslab2/emma/snpArcher/past/results/hg38/
past_proj_dir=/Genomics/ayroleslab2/emma/pastoralist_project/
scratch_dir=/scratch/tmp/emmarg/PastGWAS/
 
chr=7

Impute2_output_dir=${scratch_dir}imputed

#filter chr plink file with plink1.9 to make it compatible with KING
plink --bfile ${Impute2_output_dir}/impute1panel.${chr}.ALL --make-bed --out ${Impute2_output_dir}/impute1panel.${chr}.ALL_filt1 --maf 0.01 --mind 0.05 --geno 0.05 --hwe 0.000001 --chr 1-22 --indep-pairwise 50 20 0.8

#Run KING
~/bin/king -b ${Impute2_output_dir}/impute1panel.${chr}.ALL_filt1.bed --related --degree 2
~/bin/king -b ${Impute2_output_dir}/impute1panel.${chr}.ALL_filt1.bed --unrelated --degree 2

