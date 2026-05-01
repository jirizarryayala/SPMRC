This folder contains files needed to run a GWAS. To run a GWAS, follow the below steps:

1) Download Plink 2.0 here: https://www.cog-genomics.org/plink/2.0/
2) Download all 3 files in this folder (excluding this one) and place in your eMERGE home directory via Cypress.
3) Adjust .sh file by inputting ICD codes, output file path, type of GWAS you'd like, etc.
4) Submit the job via Cypress (type sbatch GWAS_master_script.sh) and confirm it ran (by checking the queue using squeue)

List of possible Errors and Fixes:
1) sbatch: error: instead of expected UNIX line breaks (\n). Fix: dos2unix GWAS_Script.sh
