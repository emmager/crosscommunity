#!/bin/bash
#SBATCH --job-name=G2_bam_processing    # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=16G         # memory per cpu-core
#SBATCH --array=5-20
#SBATCH --time=03:00:00          # total run time limit (HH:MM:SS)

# define input
barcode=$(sed -n ${SLURM_ARRAY_TASK_ID}p TURKANA_SAMTOOLS_SAMPLE_IDs.txt)
in_picard=/Genomics/grid/users/alea/programs/picard-tools-1.141/picard.jar
in_gatk=/Genomics/grid/users/alea/programs/gatk-4.1.4.0
in_sambamba=/Genomics/grid/users/alea/programs/sambamba_v0.6.6
in_genome=/Genomics/ayroleslab2/alea/ref_genomes/hg38/hg38_all_chr.fa
in_variants=/Genomics/ayroleslab2/yushi/ref/public_datasets/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf

# define output
out_bam1=/scratch/tmp/ayroles/alea_files/turkana_ASU/$barcode.hg38.bam
out_bam2=/Genomics/ayroleslab2/alea/turkana_wgs/from_ASU/bams/$barcode.hg38.bam

out_bam_sorted1=/scratch/tmp/emmarg/TurkanaGWAS/$barcode.hg38.sorted1.bam
out_bam_sorted2=/scratch/tmp/emmarg/TurkanaGWAS/$barcode.hg38.sorted2.bam

out_bam_duplicates=/scratch/tmp/emmarg/TurkanaGWAS/$barcode.hg38.dup.bam
out_txt_dupmetrics=/scratch/tmp/emmarg/TurkanaGWAS/$barcode.hg38.dup.txt
out_bam_readgroups=/scratch/tmp/emmarg/TurkanaGWAS/$barcode.hg38.reg.bam
out_table_recalibr=/scratch/tmp/emmarg/TurkanaGWAS/$barcode.rbqs.table
out_bam_recalibrat=/scratch/tmp/emmarg/TurkanaGWAS/$barcode.rbqs.bam
out_g_vcf=/scratch/tmp/emmarg/TurkanaGWAS/$barcode.hg38.g.vcf

module load java
module load samtools

# Step   I: Identify duplicated reads with sorted bam or sam files
# Step  II: Recalibrate base quality scores
# Step III: Calling GATK variants and export VCF

# Step I
echo 'samtools sorting...' # sort bam file

samtools sort -m 3G -o $out_bam_sorted1 $out_bam1
samtools sort -m 3G -o $out_bam_sorted2 $out_bam2
