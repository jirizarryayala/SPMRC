This folder contains files needed to run a GWAS. To run a GWAS, follow the below steps:
1. Download Plink 2.0 here: https://www.cog-genomics.org/plink/2.0/ 
2. Input ICD codes into your phenotype generation script.
3. Put phenotype generation script, GWAS script, covariate script, keep script in the same folder in your Crosslin Team home directory
4. Run covariate script in Cypress (type into cypress: sbatch Phenotype_Generation_script.sh) and confirm it ran (by checking the queue using squeue)
5. Run GWAS script (after going to the right folder with sh script, type into cypress: sbatch GWAS_Generation_script.sh) and confirm it ran (by checking the queue using squeue)


Possible errors and Fixes:

sbatch: error: instead of expected UNIX line breaks (\n).
Fix: dos2unix GWAS_Script.sh
