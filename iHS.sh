#!/bin/bash
#SBATCH --job-name=iHS_test
#SBATCH --output=iHS.out
#SBATCH --error=iHS.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=80G
#SBATCH --time=06:00:00
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu

module load R/4.1.1

Rscript /Genomics/ayroleslab2/emma/Turkana_Genotyping/1000G_data/IHS.R



