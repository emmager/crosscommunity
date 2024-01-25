#!/bin/bash
#SBATCH --job-name=filter_jointgenotyping    # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=64G         # memory per cpu-core
#SBATCH --time=24:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=emmarg@princeton.edu

module load GATK/4.3.0.0

snparcherdir=/Genomics/ayroleslab2/emma/snpArcher/past/results/hg38/
past_proj_dir=/Genomics/ayroleslab2/emma/pastoralist_project/
scratch_dir=/scratch/tmp/emmarg/PastGWAS/

#Filter to include only the reads which passed the filters that snpArcher instilled
gatk SelectVariants \
-V ${snparcherdir}past_raw.vcf.gz \
--exclude-filtered \
-O ${past_proj_dir}2024-01-18_past_filter1.vcf

bgzip ${past_proj_dir}2024-01-18_past_filter1.vcf
