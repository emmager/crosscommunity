#!/bin/bash
#SBATCH --job-name=vcf2gen   # create a short name for your job
#SBATCH --output=vcf2gen_%A_%a.out
#SBATCH --error=vcf2gen_%A_%a.err
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

# convert vcf to gen
perl ${past_proj_dir}scripts/vcf2impute_gen.pl \
 -vcf ${scratch_dir}mappingqual10/chr${chr}.hg19.GQ10.vcf \
 -gen ${scratch_dir}mappingqual10/chr${chr}.hg19.GQ10.gen \
 -chr ${chr}


