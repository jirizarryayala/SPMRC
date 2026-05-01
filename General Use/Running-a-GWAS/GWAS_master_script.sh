#!/bin/bash
#SBATCH --job-name=gwas_pipeline
#SBATCH --output=gwas_run_%j.out
#SBATCH --error=gwas_run_%j.err
#SBATCH --partition=centos7
#SBATCH --qos=long
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ssalter@tulane.edu
#SBATCH --mem=128000

# ==============================================================================
# CONFIGURABLE FLAGS (Update these for different diseases)
# ==============================================================================
OUT_DIR="/lustre/project/crosslin/crosslin_team/ssalter/hem_results_sex_stratified_reduced_controls/" #replace this with your directory
ICD_CODES="275.03,E83.110,E83.111,E83.119" # replace with your own comma-separated list of ICD codes of interest
COVAR_COLS="3-12"                   # Columns in the covariate file
DISEASE_NAME="Hemochromatosis"      #repalce with your own disease name of interest
MAF="0.05"                          #replace with your desired MAF

# STRAT_MODE options: "none" (Combined), "male", "female", or "both" (two GWAS)
STRAT_MODE="both"

# ==============================================================================
# PATHS & MODULES
# ==============================================================================
module load R/4.4.1
GENO_PATH="/lustre/project/crosslin/emerge/data/imputed_legacy/"
PLINK2_BIN="/lustre/project/crosslin/crosslin_team/ssalter/plink2"
R_GENERATE="generate_gwas_files.R"
R_PLOT="plot_gwas_results.R"

# ==============================================================================
# ORCHESTRATION LOGIC
# ==============================================================================
if [ "$STRAT_MODE" == "both" ]; then
    RUN_LIST="male female"
else
    RUN_LIST="$STRAT_MODE"
fi

for run in $RUN_LIST; do
    echo "========================================================================"
    echo "Starting GWAS Pipeline for $DISEASE_NAME | Mode: $run | $(date)"
    echo "========================================================================"

    # 1. FILE GENERATION
    # Passes: Output Directory, ICD Code String, and current Run Mode (male/female/none)
    echo "Step 1: Generating Keep, Phenotype, and Covariate files..."
    Rscript $R_GENERATE "$OUT_DIR" "$ICD_CODES" "$run"

    # 2. RUN PLINK 2 GWAS
    echo "Step 2: Running Plink 2 GWAS for Chromosomes 1-22..."
    for chr in {1..22}; do
        echo "Processing Chromosome $chr..."
        $PLINK2_BIN \
            --bfile ${GENO_PATH}emerge_chr${chr} \
            --pheno ${OUT_DIR}${run}_pheno.txt \
            --keep ${OUT_DIR}${run}_keep.txt \
            --keep-nosex \
            --logistic hide-covar \
            --maf $MAF \
            --covar ${OUT_DIR}${run}_covar.txt \
            --covar-col-nums $COVAR_COLS \
            --out ${OUT_DIR}chr${chr}.${DISEASE_NAME}_${run}.result \
            --gwas-ssf
    done

    # 3. VISUALIZATION
    # Passes: Output Directory, current Run Mode, and Disease Name
    echo "Step 3: Generating Manhattan and QQ plots..."
    Rscript $R_PLOT "$OUT_DIR" "$run" "$DISEASE_NAME"

    echo "Finished $DISEASE_NAME ($run) at $(date)"
done

echo "Workflow complete."
