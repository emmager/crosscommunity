#!/bin/bash
#SBATCH --job-name=VQSR    # create a short name for your job
#SBATCH --output=VQSR.out
#SBATCH --error=VQSR.err
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=64G         # memory per cpu-core
#SBATCH --time=24:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=end          # send email when job ends
#SBATCH --mail-user=emmarg@princeton.edu

##########
# USE GATK MIXTURE MODEL FILTERING
##########

module load conda
module load GATK/4.3.0.0
module load R/4.3.1

conda activate plink2

snparcherdir=/Genomics/ayroleslab2/emma/snpArcher/past/results/hg38/
past_proj_dir=/Genomics/ayroleslab2/emma/pastoralist_project/
scratch_dir=/scratch/tmp/emmarg/PastGWAS/

echo 'calculating VQSLOD tranches for SNPs...'

in_vcf_gz=${past_proj_dir}2024-01-18_past_filter1.SNPs1_sort.vcf.gz

#in_gatk=/Genomics/grid/users/alea/programs/gatk-4.1.4.0
in_genome=/Genomics/ayroleslab2/emma/refgenomes/hg38/hg38.fa.gz
in_resourceDIR=/Genomics/ayroleslab2/alea/ref_genomes/public_datasets

out_r_plots=${scratch_dir}2024-01-22past_allCHR.SNP1.plots.R
out_snps_recal=${scratch_dir}2024-01-22past_allCHR.SNP1.recal
out_snps_tranc=${scratch_dir}2024-01-22past_allCHR.SNP1.tranches
out_vcf_gz=${scratch_dir}2024-01-22past_allCHR.SNP1.vqsr.vcf.gz

gatk VariantRecalibrator \
            -V $in_vcf_gz \
            --trust-all-polymorphic \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 \
            -an QD -an MQ -an ReadPosRankSum -an FS -an SOR -an DP \
            -mode SNP \
            --max-gaussians 6 \
            -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $in_resourceDIR/hapmap_3.3.hg38.vcf.gz \
            -resource:omni,known=false,training=true,truth=true,prior=12.0 $in_resourceDIR/1000G_omni2.5.hg38.vcf.gz \
            -resource:1000G,known=false,training=true,truth=false,prior=10.0 $in_resourceDIR/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
            -O $out_snps_recal \
            --tranches-file $out_snps_tranc \
            -R $in_genome \
            --rscript-file $out_r_plots

echo 'filtering SNPs with VQSLOD...'

gatk ApplyVQSR \
            -V $in_vcf_gz \
            --recal-file $out_snps_recal \
            --tranches-file $out_snps_tranc \
            --truth-sensitivity-filter-level 90.0 \
            --create-output-variant-index true \
            -mode SNP \
            -O $out_vcf_gz \
            -R $in_genome

# create one VCF per chromosome
for chr in {1..22}; do plink2 --vcf $out_vcf_gz --var-filter --recode vcf --out ${scratch_dir}2024-01-22past_allCHR.SNP1.${chr}.vcf.gz --chr ${chr}; echo ${chr}; done

