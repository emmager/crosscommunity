#!/bin/bash

#SBATCH --job-name=RFmix
#SBATCH --output=RFmixChr.%a.out
#SBATCH --error=RFmixChr.%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=120G
#SBATCH --time=48:00:00
#SBATCH --array=1-6
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu


############
# RFMix Requirements
############

#	A phased VCF/BCF file containing "query" haplotypes which are to be analyzed.
#	A phased VCF/BCF file containing reference haplotypes (in any order)
#	A reference sample map file matching reference samples to their respective reference populations
#	A genetic map file
#	All files mapped or referenced to the same genome assembly
#	All diploid genotypes must be phase resolved
#	In addition to the above, you must at minimum also specify a basename (prefix) for output files, and the chromosome to analyze from the VCF/BCF inputs, even if they contain only one chromosome
#	If BCF files are used, bcftools must be installed and available in the PATH environment setting
#	VCF/BCF files may be gzip compressed, and should be indexed using bcftools

#	RFMIX upon completion will output two main files of interest: the most likely assignment of subpopulations per CRF point (<output basename>.msp.tsv), and the marginal probabilities of each subpopulation being the ancestral population of the corresponding CRF point (<output basename>.fb.tsv). 

#############
# RUN IN BASH
# NOTE: CHROM is replaced with the chromosome name
#############

module load bcftools 

CHROM=${SLURM_ARRAY_TASK_ID}

# run RFmix- running directory= /Genomics/ayroleslab2/emma/Turkana_Genotyping/1000G_data

map_path=/Genomics/ayroleslab2/alea/tsimane/int_data_files/temp.genetic_map_chr${CHROM}.txt

/Genomics/grid/users/alea/programs/rfmix-master/rfmix -f ALL.chr${CHROM}.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes_RFmix_Query.vcf.gz -r ALL.chr${CHROM}.phase3_v5.shapeit2_mvncall_integrated.noSingleton.genotypes_RFmix_REF.vcf.gz -m 1000G_RFmix_Ref_sample_map  -g $map_path -o RFmix_chr${CHROM} --chromosome=${CHROM} -G 10 

echo "done"
