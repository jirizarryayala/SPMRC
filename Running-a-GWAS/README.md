This folder contains files needed to run a GWAS. To run a GWAS, follow the below steps:
1. Input ICD codes into your phenotype generation script.
2. Put phenotype generation script, GWAS script, covariate script, keep script in the same folder in your Crosslin Team home directory
3. Run covariate script in Cypress (type into cypress: sbatch Phenotype_Generation_script.sh) and confirm it ran (by checking the queue using squeue)
4. Run GWAS script (type into cypress: sbatch GWAS_Generation_script.sh) and confirm it ran (by checking the queue using squeue)


Possible errors and Fixes:
sbatch: error: instead of expected UNIX line breaks (\n).
Fix: dos2unix GWAS_Script.sh
