#!/bin/bash
#SBATCH --job-name=filt_low_cov   # create a short name for your job
#SBATCH --output=filt_low_cov_%A_%a.out
#SBATCH --error=filt_low_cov_%A_%a.err
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=65G         # memory per cpu-core
#SBATCH --time=01:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=emmarg@princeton.edu
#SBATCH --array=1-22


module load bcftools/1.10.2
module load GATK/4.3.0.0
module load conda

conda activate plink2

snparcherdir=/Genomics/ayroleslab2/emma/snpArcher/past/results/hg38/
past_proj_dir=/Genomics/ayroleslab2/emma/pastoralist_project/
scratch_dir=/scratch/tmp/emmarg/PastGWAS/
 
chr=$(sed -n ${SLURM_ARRAY_TASK_ID}p ~/Chrs.txt)


# prep files for imputation 
# remove individuals not genotyped with GQ10 at 5% of loci
# extract good SNPs called from high cov

grep "PASS" ${scratch_dir}2024-01-22past_allCHR.SNP1.${chr}.vcf.gz.vcf | tail -n +3 | awk '{OFS=""; print $1,"\t",$2,"\t",$2+1,"\t",$1,"_",$2}' > temp_${chr}.range

#alea joint genotyping method
plink2 --vcf ${past_proj_dir}2024-01-18_past_filter1.vcf.gz --allow-extra-chr --chr 1-22 --geno 0.95 --recode vcf bgz --snps-only --vcf-min-gq 10 --out ${scratch_dir}mappingqual10/chr${chr}.GQ10 --extract range temp_${chr}.range --max-alleles 2


