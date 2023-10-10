#!/bin/bash
#SBATCH --job-name=G2_bam_processing    # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=65G         # memory per cpu-core
#SBATCH --time=66:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=emmarg@princeton.edu
#SBATCH --array=1-17

# define input
bam1=$(sed -n ${SLURM_ARRAY_TASK_ID}p TURKANA_SAMTOOLS1_SAMPLE_IDs.txt)
bam2=$(sed -n ${SLURM_ARRAY_TASK_ID}p TURKANA_SAMTOOLS2_SAMPLE_IDs.txt)
barcode=$(sed -n ${SLURM_ARRAY_TASK_ID}p TURKANA_SAMPLE_IDs.txt)
in_picard=/Genomics/grid/users/alea/programs/picard-tools-1.141/picard.jar
in_gatk=/Genomics/grid/users/alea/programs/gatk-4.1.4.0
in_sambamba=/Genomics/grid/users/alea/programs/sambamba_v0.6.6
in_genome=/Genomics/ayroleslab2/alea/ref_genomes/hg38/hg38_all_chr.fa
in_variants=/Genomics/ayroleslab2/yushi/ref/public_datasets/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf

# define output
out_bam1=/scratch/tmp/ayroles/alea_files/turkana_ASU/$bam1.hg38.bam
out_bam2=/Genomics/ayroleslab2/alea/turkana_wgs/from_ASU/bams/$bam2.hg38.bam

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

#samtools sort -m 3G -o $out_bam_sorted1 $out_bam1
#samtools sort -m 3G -o $out_bam_sorted2 $out_bam2

echo 'picard duplications...' # run picard to mark duplicates

java -Xmx65g -jar $in_picard MarkDuplicates I=$out_bam_sorted1 I=$out_bam_sorted2 O=$out_bam_duplicates M=$out_txt_dupmetrics

# Step II
echo 'adding read groups...'

java -Xmx65g -jar $in_picard AddOrReplaceReadGroups I=$out_bam_duplicates O=$out_bam_readgroups SO=coordinate RGLB=$barcode RGPL=illumina RGPU=turkana RGSM=$barcode

rm -f $out_bam_duplicates
rm -f $out_txt_dupmetrics

echo 'pre-indexing...'

$in_sambamba index -t 4 $out_bam_readgroups

echo 'recalibrating base quality scores...' # sort bam file

$in_gatk/gatk BaseRecalibrator  -I $out_bam_readgroups -R $in_genome --known-sites $in_variants -O $out_table_recalibr

$in_gatk/gatk ApplyBQSR -R $in_genome -I $out_bam_readgroups -bqsr $out_table_recalibr -O $out_bam_recalibrat

rm -f $out_bam_readgroups

echo 'post-indexing...'

$in_sambamba index -t 4 $out_bam_recalibrat

# Step III
echo 'calling GATK variants as VCF and zip...' # sort bam file

$in_gatk/gatk --java-options "-Xmx65g" HaplotypeCaller -R $in_genome -I $out_bam_recalibrat -O $out_g_vcf -ERC GVCF --max-alternate-alleles 2

bgzip $out_g_vcf
tabix -p vcf $out_g_vcf.gz

rm $out_bam_recalibrat
