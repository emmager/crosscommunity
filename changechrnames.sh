#!/bin/bash

#SBATCH --job-name=Plink
#SBATCH --output=changechrname
##SBATCH --error=plink.%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=4:00:00
##SBATCH --array=0-22
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu

module load bcftools
bcftools annotate --rename-chrs chr_name_conv.txt ALLCHR.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes_subset.vcf.gz -Oz -o ALLCHR.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes_integernames.vcf.gz
