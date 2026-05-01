# Load libraries
library(data.table)

# --- 1. CAPTURE COMMAND LINE ARGUMENTS ---
# Expected order: 1=OutputDir, 2=RunMode (male/female/none), 3=DiseaseName
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Missing arguments. Usage: Rscript plot_gwas.R <OutputDir> <RunMode> <DiseaseName>")
}

dir          <- args[1] # Results Directory
run_mode     <- args[2] # "male", "female", or "none"
disease_name <- args[3] # e.g., "Hemochromatosis"

print(paste("Working"))


# --- 2. LOAD PLINK RESULTS ---
# Standardize the file naming pattern to match your Bash/Plink script outputs
chroms <- 1:22
files  <- paste0(dir, "chr", chroms, ".", disease_name, "_", run_mode, ".result.PHENO1.glm.logistic.hybrid")

print(paste("Searching for results in:", dir))

list_of_tables <- lapply(seq_along(files), function(i) {
  f <- files[i]
  if (file.exists(f)) {
    return(fread(f))
  } else {
    warning(paste("File not found:", f))
    return(NULL)
  }
})
warnings()

# Merge all chromosome results into one table
p_all <- rbindlist(list_of_tables, fill = TRUE)

# --- 3. DATA CLEANING & STANDARDIZATION ---
# Ensure columns match the expected names for plotting functions
if ("#CHROM" %in% names(p_all)) { setnames(p_all, "#CHROM", "CHR") }
p_all[, CHR := as.numeric(CHR)]
p_all[, P := as.numeric(P)]
p_all <- p_all[!is.na(P)] # Remove rows with missing P-values

# --- 4. LOAD PLOTTING FUNCTIONS ---
source("/lustre/project/crosslin/crosslin_team/ssalter/manhattan_plots_package/manhattan.plot.drc.v1.R")
source("/lustre/project/crosslin/crosslin_team/ssalter/QQ_plots_package/qq.plot.v7.r")

# --- 5. GENERATE MANHATTAN PLOT ---
# Dynamically create the title and filename based on the flags provided
titl <- paste("Manhattan Plot of", run_mode, disease_name, "GWAS")
manhattan_filename <- paste0(disease_name, "_", run_mode, "_manhattan.png")

png(file=file.path(dir, manhattan_filename), width=1400, height=700, pointsize=18, type="cairo")
stripe(p=p_all$P, chromosome=p_all$CHR, chromcode=1:22, main=titl, win=FALSE)
dev.off()

# --- 6. GENERATE QQ PLOT ---
qpval  <- p_all$P
lambda <- median(-2 * log(qpval)) / 1.39
titl2  <- paste("QQ Plot of", run_mode, disease_name, "GWAS")
qq_filename <- paste0(disease_name, "_", run_mode, "_qq.png")

png(file=file.path(dir, qq_filename), width=384, height=384, type="cairo")
qq.plot(pval=qpval, trunc=FALSE, main=titl2, sub=paste("lambda =", format(lambda, digits=4)))
dev.off()

print(paste("Successfully generated plots for", disease_name, "under mode:", run_mode))