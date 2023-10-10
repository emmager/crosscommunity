#!/bin/bash
#SBATCH --job-name=G6_SNP_filtering    # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=65G         # memory per cpu-core
#SBATCH --time=01:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=emmarg@princeton.edu
#SBATCH --array=1-22

# NOTE: NUMBER was replaced with the chromosome number

module load bcftools
module load plink/2.00a3LM
module load java
 
chr=$(sed -n ${SLURM_ARRAY_TASK_ID}p ~/Chrs.txt)

# combine output of joint calling by chr
#cd /scratch/tmp/emmarg/TurkanaGWAS/low_cov_all/
$bcftools concat -Oz -a -o chr${chr}.low_cov.vcf.gz chr${chr}:*.vcf.gz 

# prep files for imputation 
# remove individuals not genotyped with GQ10 at 5% of loci
# extract good SNPs called from high cov

#grep "PASS" /scratch/tmp/ayroles/alea_files/turkana_high_cov_Mar2020/high_cov_allCHR.SNP1.${chr}.vcf.gz.vcf | tail -n +2 | awk '{OFS=""; print $1,"\t",$2,"\t",$2+1,"\t",$1,"_",$2}' > temp_${chr}.range

#$plink2 --vcf /scratch/tmp/emmarg/TurkanaGWAS/low_cov_all/chr${chr}.low_cov.vcf.gz --mind 0.95 --recode vcf bgz --snps-only --vcf-min-gq 10 --out /scratch/tmp/emmarg/TurkanaGWAS/low_cov_all/chr${chr}.low_cov.GQ10 --extract range temp_${chr}.range

# chr=1
# /Genomics/grid/users/alea/programs/plink_1.90 --vcf /scratch/tmp/ayroles/alea_files/low_cov_all/chr${chr}.low_cov.vcf.gz --mind 0.95 --snps-only --vcf-min-gq 10 --out /scratch/tmp/ayroles/alea_files/low_cov_all/chr${chr}.low_cov.GQ10_het --extract range temp_${chr}.range --het
# /Genomics/grid/users/alea/programs/plink_1.90 --vcf /scratch/tmp/ayroles/alea_files/low_cov_all/chr${chr}.low_cov.vcf.gz --mind 0.95 --snps-only --vcf-min-gq 10 --out /scratch/tmp/ayroles/alea_files/low_cov_all/chr${chr}.low_cov.GQ10_het --extract range temp_${chr}.range --missing
# /Genomics/grid/users/alea/programs/plink_1.90 --vcf /scratch/tmp/ayroles/alea_files/low_cov_all/chr${chr}.low_cov.vcf.gz --mind 0.95 --snps-only --vcf-min-gq 10 --out /scratch/tmp/ayroles/alea_files/low_cov_all/chr${chr}.low_cov.GQ10_het --extract range temp_${chr}.range --genome

# liftover low cov VCFs to hg19
#in_vcf=/scratch/tmp/emmarg/TurkanaGWAS/low_cov_all/chr${chr}.low_cov.GQ10.vcf.gz
#in_picard=/Genomics/grid/users/alea/programs/picard-tools-2.21.6/picard.jar
#in_chain=/Genomics/grid/users/alea/programs/hg38ToHg19.over.chain.gz
#in_genome=/Genomics/ayroleslab2/alea/ref_genomes/hg19/hg19_v0_Homo_sapiens_assembly19.fasta

#out_reject=/scratch/tmp/emmarg/TurkanaGWAS/low_cov_all/chr${chr}.low_cov.reject.GQ10.vcf.gz
#out_lifted_vcf=/scratch/tmp/emmarg/TurkanaGWAS/low_cov_all/chr${chr}.low_cov.hg19.GQ10.vcf.gz

#gzip -cd $in_vcf | grep -v '*' | sed -e s/^${chr}/chr${chr}/g > /scratch/tmp/emmarg/TurkanaGWAS/temp/temp1_${chr}.vcf

#java -jar $in_picard LiftoverVcf \
 #  I=/scratch/tmp/emmarg/TurkanaGWAS/temp/temp1_${chr}.vcf \
 #  O=$out_lifted_vcf \
 #  CHAIN=$in_chain \
 #  REJECT=$out_reject \
 #  R=$in_genome WARN_ON_MISSING_CONTIG=true

# get variants only lifted to the focal chromosome
#plink=/Genomics/grid/users/alea/programs/plink_1.90

#$plink --vcf $out_lifted_vcf --recode-vcf --chr chr${chr} --out /scratch/tmp/emmarg/TurkanaGWAS/temp/temp2_${chr}

#rm /scratch/tmp/emmarg/TurkanaGWAS/temp/temp1_${chr}.vcf
#rm $out_reject
#rm $out_lifted_vcf

#sed -e 's/chr//g' /scratch/tmp/emmarg/TurkanaGWAS/temp/temp2_${chr}.vcf > /scratch/tmp/emmarg/TurkanaGWAS/Turkana-low-cov/chr${chr}.hg19.low_cov.GQ10.vcf

# convert vcf to gen
perl /Genomics/grid/users/alea/programs/vcf2impute_gen.pl \
 -vcf /scratch/tmp/emmarg/TurkanaGWAS/Turkana-low-cov/chr${chr}.hg19.low_cov.GQ10.vcf \
 -gen /scratch/tmp/emmarg/TurkanaGWAS/Turkana-low-cov/chr${chr}.hg19.low_cov.GQ10.gen \
 -chr ${chr}

