#!/bin/bash

#SBATCH --job-name=Plink
#SBATCH --output=pca_allchr_int.out
##SBATCH --error=plink.%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=4:00:00
##SBATCH --array=0-22
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu

module load plink


#vcfs=(/Genomics/ayroleslab2/emma/Turkana_Genotyping/1000G_data/*.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes_subset.vcf)

#this_vcf=${vcfs[${SLURM_ARRAY_TASK_ID}]}

#prefix=`echo ${this_vcf} | sed 's/.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes_subset.vcf//g'| sed 's/^\/Genomics\/ayroleslab2\/yushi\/ref\/public_datasets\/1000GP_VCF_P3\///'`

#echo ${this_vcf}

module load bcftools
bcftools index ALLCHR.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes_integernames.vcf.gz

plink2 --vcf ALLCHR.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes_integernames.vcf.gz --maf 0.05 --make-bed --chr 1-25 --out ALLCHR.phase3_v5.shapeit2_mvncall_integrated_int_subset

plink2 --bfile ALLCHR.phase3_v5.shapeit2_mvncall_integrated_int_subset --out ALLCHR.phase3_v5.shapeit2_mvncall_integrated_int_subset --pca 10

#These output eigenvalue and eigenvector files will be your input to the script pca.R. This script is for visualizing your principal components in ggplot
