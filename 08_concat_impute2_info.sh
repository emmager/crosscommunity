#!/bin/bash
#SBATCH --job-name=concatimputeinfo
#SBATCH --output=concatimputeinfo.out
#SBATCH --error=concatimputeinfo.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G
#SBATCH --time=1:00:00
##SBATCH --array=0-1
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu

# concat imputation data for one chromosome
##########

touch /Genomics/ayroleslab2/emma/pastoralist_project/Jan2024_analyses/impute1panel_chr7.hg19_GQ10_info

for f in /scratch/tmp/emmarg/PastGWAS/imputed/impute1panel.7.*.hg19_GQ10_info; do awk '{OFS="\t"; print $0,FILENAME}' $f | tail -n+2 >> /Genomics/ayroleslab2/emma/pastoralist_project/Jan2024_analyses/impute1panel_chr7.hg19_GQ10_info; done

cat /scratch/tmp/emmarg/PastGWAS/imputed/impute1panel.7.*.hg19_GQ10_info_by_sample | grep -v 'concord_type'> /Genomics/ayroleslab2/emma/pastoralist_project/Jan2024_analyses/impute1panel_chr7.hg19_GQ10_info_by_sample

