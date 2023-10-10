#!/bin/bash

#SBATCH --job-name=Subset_vcfs
#SBATCH --output=Subset_vcfs.%a.out
#SBATCH --error=Subset_vcfs.%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=4:00:00
#SBATCH --array=0-22
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu


module load bcftools

vcfs=(/Genomics/ayroleslab2/yushi/ref/public_datasets/1000GP_VCF_P3/*.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes.vcf.gz)

this_vcf=${vcfs[${SLURM_ARRAY_TASK_ID}]}

prefix=`echo ${this_vcf} | sed 's/.vcf.gz//g'| sed 's/^\/Genomics\/ayroleslab2\/yushi\/ref\/public_datasets\/1000GP_VCF_P3\///'`

#echo ${this_vcf}
#echo ${prefix}


#bcftools view --samples-file 2023_08_16_1000G_P3_subset_IDs ${this_vcf} > ${prefix}_subset.vcf

bcftools view ${prefix}_subset.vcf -Oz -o ${prefix}_subset.vcf.gz

bcftools index -t ${prefix}_subset.vcf.gz
