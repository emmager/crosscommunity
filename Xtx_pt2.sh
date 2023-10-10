#!/bin/bash

#SBATCH --job-name=Xtx2
#SBATCH --output=Xtx2_test.out
#SBATCH --error=Xtx2_test.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=160G
#SBATCH --time=300:00:00
##SBATCH --array=1
#SBATCH --mail-type=end
#SBATCH --mail-user=emmarg@princeton.edu


# WGS order = ASW, CEU, JPT, YRI
bayenv2=/Genomics/grid/users/alea/programs/tguenther-bayenv2_public-2b2b7f20bb62/bayenv2
SNPSFILE=/Genomics/ayroleslab2/emma/Turkana_Genotyping/1000G_data/2023-09-06_BayEnv_chr6_w1000G_subset.txt
ENVFILE=/Genomics/ayroleslab2/emma/Turkana_Genotyping/1000G_data/2023-09-06_Bayenv2_ENVIRON.txt #The environmental file is required but not used for calculation so it can just include some dummy variables
POPNUM=4
ITNUM=100000
ENVNUM=1
rnd1=$(perl -e 'printf("%05d",rand(99999))')
rnd2=$(perl -e 'printf("%05d",rand(99999))')
rnd3=$(perl -e 'printf("%05d",rand(99999))')

#Make Covariance matrix
#$bayenv2 -i $SNPSFILE -p $POPNUM -k 100000 -r 63479 > Bayenv2_1000G_subset_MATFILE.out

MATFILE=/Genomics/ayroleslab2/emma/Turkana_Genotyping/1000G_data/Bayenv2_1000G_subset_MATFILE.txt


# rm 28May20_BayEnv_chrCHROMNUM_output* #Don't know what this line is for
#Make list of even numbers from 2:length of SNPSFILE. This is so you can run bayenv2 on 1 SNP at a time. This is the SNPFILE
f=$(< "$SNPSFILE" wc -l)
seq 2 2 $f > temp_chr6_loci.txt

#Calculate Bayes factor for 3 random seeds?
for f in `cat temp_chr6_loci.txt`; 
do head -$f $SNPSFILE | tail -2 > 2023-09-28_temp_chr6_SNP$f.txt 
$bayenv2 -i 2023-09-28_temp_chr6_SNP$f.txt -e $ENVFILE -m $MATFILE -k $ITNUM -r $RANDOM -p $POPNUM -n $ENVNUM -t -X -r $rnd1 -o 2023-09-28_BayEnv_chr6_output1; 
$bayenv2 -i 2023-09-28_temp_chr6_SNP$f.txt -e $ENVFILE -m $MATFILE -k $ITNUM -r $RANDOM -p $POPNUM -n $ENVNUM -t -X -r $rnd2 -o 2023-09-28_BayEnv_chr6_output2; 
$bayenv2 -i 2023-09-28_temp_chr6_SNP$f.txt -e $ENVFILE -m $MATFILE -k $ITNUM -r $RANDOM -p $POPNUM -n $ENVNUM -t -X -r $rnd3 -o 2023-09-28_BayEnv_chr6_output3; 
rm 2023-09-28_temp_chr6_SNP$f.txt; done
