#!/bin/bash
#SBATCH --job-name=G3_joint_call    # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=80G         # memory per cpu-core
#SBATCH --time=01:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=emmarg@princeton.edu
#SBATCH --array=2-587%40

# NOTE: INTERVAL was replaced with a genomic interval of 5000000bp

module load java
module load samtools

chrom=$(sed -n ${SLURM_ARRAY_TASK_ID}p /scratch/tmp/emmarg/TurkanaGWAS/Interval_list.txt)
in_gatk=/Genomics/grid/users/alea/programs/gatk-4.1.4.0
in_genome=/Genomics/ayroleslab2/yushi/ref/hg38_all_chr.fa
in_database=/scratch/tmp/emmarg/TurkanaGWAS/low_cov_all/$chrom
in_samples=/scratch/tmp/emmarg/TurkanaGWAS/all_low_cov.sample_map
out_vcf_gz=/scratch/tmp/emmarg/TurkanaGWAS/low_cov_all/$chrom.all_low_cov.vcf.gz

# merge GVCFs
#echo $chrom

#$in_gatk/gatk --java-options "-Xmx80g -Xms80g" GenomicsDBImport \
 #           --sample-name-map $in_samples \
  #          --genomicsdb-workspace-path /scratch/tmp/emmarg/TurkanaGWAS/low_cov_all/$chrom \
   #         --tmp-dir=/scratch/tmp/emmarg/TurkanaGWAS/ \
    #        -L $chrom \
     #     #  --batch-size 40 --reader-threads 16

# call genotypes

echo 'calling joint genotypes...'

echo $chrom > temp.${chrom}.list

$in_gatk/gatk --java-options "-Xmx80g" GenotypeGVCFs \
            -R $in_genome \
            -L temp.${chrom}.list \
            -V gendb://$in_database \
            -O $out_vcf_gz \
            --tmp-dir=/scratch/tmp/emmarg/TurkanaGWAS/ --max-alternate-alleles 2
