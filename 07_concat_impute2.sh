#!/bin/bash
#SBATCH --job-name=concatimpute2
#SBATCH --output=concatimpute2.out
#SBATCH --error=concatimpute2.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=1:00:00
##SBATCH --array=0-1
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu

module load conda
conda activate /Genomics/argo/users/emmarg/.conda/envs/plink2
module load bcftools

snparcherdir=/Genomics/ayroleslab2/emma/snpArcher/past/results/hg38/
past_proj_dir=/Genomics/ayroleslab2/emma/pastoralist_project/
scratch_dir=/scratch/tmp/emmarg/PastGWAS/
 
chr=7

Impute2_output_dir=${scratch_dir}imputed

#concatenate all the 5MB window Impute2 output files into 1 .gen file per chromosome
cat ${Impute2_output_dir}/impute1panel.${chr}.*.hg19_GQ10 > ${Impute2_output_dir}/impute1panel.${chr}.ALL.hg19_GQ10

#sort
sort -nk3 ${Impute2_output_dir}/impute1panel.${chr}.ALL.hg19_GQ10 > ${Impute2_output_dir}/impute1panel.${chr}.ALL.hg19_sorted

# convert imputed .gen file to plink binary format
plink2 --gen ${Impute2_output_dir}/impute1panel.${chr}.ALL.hg19_sorted ref-first --sample ${past_proj_dir}plate2_past.sample  --oxford-single-chr ${chr} --make-bed --out ${Impute2_output_dir}/imp$

head ${past_proj_dir}plate2_past.sample

