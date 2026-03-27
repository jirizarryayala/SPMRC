#!/bin/bash
#SBATCH --job-name=tg_gwas_f    ### Job Name
#SBATCH --output=/lustre/project/crosslin/crosslin_team/ssalter/hem_results/tg_gwas_females_%A_%a.out       ### File in which to store job output
#SBATCH --error=/lustre/project/crosslin/crosslin_team/ssalter/hem_results/tg_gwas_females_%A_%a.err        ### File in which to store job error messages
#SBATCH --partition=centos7      ###Partition in which to start session
#SBATCH --qos=long          ### Quality of Service (like a queue in PBS)
#SBATCH --time=1-00:00:00     ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1             ### Node count required for the job
#SBATCH --ntasks-per-node=1   ### Number of tasks to be launched per Node
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ssalter@tulane.edu
#SBATCH --mem=128000
### Note: To run this code on Cypress, users need to download and configure plink2. Email jirizarry@cypress.tulane.edu so I can forward you the instructions.
###idev --partition=centos7 --mem=256000 -t 15

path="/lustre/project/crosslin/emerge/data/imputed_legacy/" #genotype data
path2="/lustre/project/crosslin/crosslin_team/ssalter" #where phenotype, keep, and covariate files are stored
path3="/lustre/project/crosslin/crosslin_team/ssalter/Biplob_Project/" #where plink writes gwas results
for chr in {1...22}
do
echo "Starting GWAS for Chromosome $chr at $(date)"
./plink2 \
--bfile ${path}emerge_chr${chr} \
--pheno ${path2}Phenotype_Generation_Script_GWAS.txt --1 \ #This is the file we need to change
--keep ${path2}hem_keep.txt \
--logistic hide-covar \
--maf 0.05 \
--covar ${path2}hem_covar2.txt \
--covar-col-nums 3-7 \
--out ${path3}chr${chr}.Biplob_GWAS.result \
--gwas-ssf
echo "Finished Chromosome $chr at $(date)"
done
date