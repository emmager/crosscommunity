#!/bin/bash
#SBATCH --job-name=highcovsnps    # create a short name for your job
#SBATCH --output=highcovsnps.out
#SBATCH --error=highcovsnps.err
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

# check for problematic samples (and only look in chr 1-22) on filtered vcf 
plink2  --vcf ${past_proj_dir}2024-01-18_past_filter1.vcf.gz --snps-only --missing --out ${past_proj_dir}missing --allow-extra-chr --chr 1-22


# check for problematic samples (and only look in chr 1-22) on raw vcf
plink2  --vcf ${snparcherdir}past_raw.vcf.gz --snps-only --missing --out ${past_proj_dir}rawpast_missing --allow-extra-chr --chr 1-22


##########
# USE HARD FILTERING FIRST TO REMOVE STUFF WE DEFINITELY WON'T ANALYZE
##########

# region filtering
cpgs=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/CpGs_hg38_v2.bed
mask=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/20160622.allChr.pilot_mask.bed # Check if there is a better version of this. 
super=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets/genomicSuperDups_hg38.bed
cat $cpgs $super > /Genomics/ayroleslab2/emma/refgenomes/temp.bed

awk '$4=(FNR FS)' /Genomics/ayroleslab2/emma/refgenomes/temp.bed > /Genomics/ayroleslab2/emma/refgenomes/temp_v2.bed

awk '$1 !~ /_/' /Genomics/ayroleslab2/emma/refgenomes/temp_v2.bed > /Genomics/ayroleslab2/emma/refgenomes/temp_v3.bed

# GATK reccomends filtering for excess hets before VQSR, HWE filtering should do something similar: https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR$
plink2 --vcf ${past_proj_dir}2024-01-18_past_filter1.vcf.gz --snps-only --hwe 0.000001 --make-just-bim --maf 0.01 --extract range $mask --exclude range /Genomics/ayroleslab2/emma/refgenomes/temp_v3.bed $

# turn into bed
awk '{OFS=""; print "chr",$1,"\t",$4-1,"\t",$4}' ${scratch_dir}2024-01-18_past_filter1SNPs.bim >  ${scratch_dir}2024-01-18_past_filter1.SNPs_loc1.bed

awk '{OFS=""; print "chr",$1,"\t",$4,"\t",$4+1}' ${scratch_dir}2024-01-18_past_filter1SNPs.bim >  ${scratch_dir}2024-01-18_past_filter1.SNPs_loc2.bed

tabix -p vcf ${past_proj_dir}2024-01-18_past_filter1.vcf.gz

# remove low MAF and problem regions, remove 2 problem samples
module load bcftools

regions=${scratch_dir}2024-01-18_past_filter1.SNPs_loc1.bed         
bcftools view -Oz -o ${past_proj_dir}2024-01-18_past_filter1.SNPs1.vcf.gz --samples ^2023-03-05_0093_08,2023-02-20_0009_21 --regions-file $regions ${past_proj_dir}2024-01-18_past_filter1.vcf.gz

regions=${scratch_dir}2024-01-18_past_filter1.SNPs_loc2.bed
bcftools view -Oz -o ${past_proj_dir}2024-01-18_past_filter1.SNPs2.vcf.gz --samples ^2023-03-05_0093_08,2023-02-20_0009_21 --regions-file $regions ${past_proj_dir}2024-01-18_past_filter1.vcf.gz

# filter and sort
bcftools sort -Oz -o ${past_proj_dir}2024-01-18_past_filter1.SNPs1_sort.vcf.gz ${past_proj_dir}2024-01-18_past_filter1.SNPs1.vcf.gz

gatk IndexFeatureFile -I ${past_proj_dir}2024-01-18_past_filter1.SNPs1_sort.vcf.gz
