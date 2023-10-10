#!/bin/bash

#SBATCH --job-name=Merge_Vcfs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=4:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu


module load gatk

gatk MergeVcfs \
	-I=input_variant_files.list \
	-O= ALLCHR.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes_RFmix_Query.vcf.gz
