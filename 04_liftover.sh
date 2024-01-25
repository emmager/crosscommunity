#!/bin/bash
#SBATCH --job-name=liftover    # create a short name for your job
#SBATCH --output=highcovsnps_%A_%a.out
#SBATCH --error=highcovsnps_%A_%a.err
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


# liftover low cov VCFs to hg19
in_vcf=${scratch_dir}mappingqual10/chr${chr}.GQ10.vcf.gz
#in_picard=/Genomics/grid/users/alea/programs/picard-tools-2.21.6/picard.jar -- I don't need this bc argo has GATK module built in
in_chain=~/hg38ToHg19.over.chain.gz
in_genome=/Genomics/ayroleslab2/emma/refgenomes/hg19/hg19_v0_Homo_sapiens_assembly19.fasta

out_reject=${scratch_dir}mappingqual10/chr${chr}.reject.GQ10.vcf.gz
out_lifted_vcf=${scratch_dir}mappingqual10/chr${chr}.hg19.GQ10.vcf.gz

gzip -cd $in_vcf | grep -v '*' | sed -e s/^${chr}/chr${chr}/g > ${scratch_dir}temp/temp1_${chr}.vcf

gatk LiftoverVcf \
   I=${scratch_dir}temp/temp1_${chr}.vcf \
   O=$out_lifted_vcf \
   CHAIN=$in_chain \
   REJECT=$out_reject \
   R=$in_genome WARN_ON_MISSING_CONTIG=true

# get variants only lifted to the focal chromosome
#plink=/Genomics/grid/users/alea/programs/plink_1.90 -- using the plink 2.0 I downloaded in the plink conda environment

plink2 --vcf $out_lifted_vcf --recode vcf --chr chr${chr} --out ${scratch_dir}temp/temp2_${chr}

rm ${scratch_dir}temp/temp1_${chr}.vcf
rm $out_reject
rm $out_lifted_vcf

#Remove "chr" from chromosome names
sed -e 's/chr//g' ${scratch_dir}temp/temp2_${chr}.vcf > ${scratch_dir}mappingqual10/chr${chr}.hg19.GQ10.vcf
