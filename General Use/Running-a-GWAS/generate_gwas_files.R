# Load libraries
library(data.table)
library(dplyr)
library(MatchIt)

# --- 1. CAPTURE COMMAND LINE ARGUMENTS ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) { 
  stop("Missing arguments. Usage: Rscript script.R <OutputDir> <ICD_Codes> <RunMode>") 
}

dir3      <- args[1] 
# Use trimws to prevent matching failures caused by spaces in the argument string
icd_list  <- trimws(unlist(strsplit(args[2], ","))) 
run_mode  <- args[3] 

# --- 2. LOAD DATA ---
if(!dir.exists(dir3)) dir.create(dir3, recursive = TRUE)

icdDir <- "/lustre/project/crosslin/emerge/data/phenotype_data/icd_codes/output/"
dir2   <- "/lustre/project/crosslin/emerge/data/accessory_files/"

icd <- fread(file=paste0(icdDir, "outputfile1.csv"))
man <- fread(file=paste0(dir2, "chr1-22.plink_pca.e123_imputation_sample_manifest.csv"))

# --- 3. DYNAMIC CASE IDENTIFICATION ---
# FIXED: Removed stray symbols ($) that were causing syntax errors
case_idx <- icd$ICD_CODE %in% icd_list
case_ids <- unique(icd$SUBJID[case_idx])

# --- 4. PREPARE AGE & MATCH POOL ---
age_df <- icd %>% 
  filter(!is.na(AGE_AT_EVENT), AGE_AT_EVENT > 0, AGE_AT_EVENT < 120) %>% 
  group_by(SUBJID) %>% 
  summarise(age_for_match = max(AGE_AT_EVENT), .groups = "drop")

match_pool <- man
if (run_mode == "male") { 
  match_pool <- match_pool %>% filter(sex2 == "male") 
} else if (run_mode == "female") { 
  match_pool <- match_pool %>% filter(sex2 == "female") 
}

match_pool <- match_pool %>% 
  select(IID, sex2, site, contains("PC")) %>% 
  inner_join(age_df, by = c("IID" = "SUBJID")) %>% 
  mutate(Status = ifelse(IID %in% case_ids, 1, 0)) %>% 
  filter(!is.na(site), !is.na(age_for_match))

# --- 5. EXECUTE MATCHING ---
num_cases <- sum(match_pool$Status == 1)
print(paste("Cases available for matching in mode", run_mode, ":", num_cases))

if(num_cases == 0) stop(paste("No cases found for mode:", run_mode))

match_obj <- matchit(Status ~ age_for_match, data = match_pool, 
                     method = "nearest", ratio = 4, caliper = 0.2, exact = ~ site)

matched_df <- match.data(match_obj)
fullData <- matched_df %>% 
  rename(iid = IID) %>% 
  mutate(fid = iid) %>% 
  mutate(CaseControl_PLINK = Status + 1)

# --- 6. ONE-HOT ENCODING FOR RACE & ETHNICITY ---
eth_raw <- man %>% select(IID, self_reported_race, ethnicity2)
fullData <- left_join(fullData, eth_raw, by = c("iid" = "IID"))

race_dummies <- model.matrix(~ 0 + self_reported_race, data = fullData)
eth_dummies  <- model.matrix(~ 0 + ethnicity2, data = fullData)
colnames(race_dummies) <- gsub("self_reported_race", "race_", colnames(race_dummies))
colnames(eth_dummies)  <- gsub("ethnicity2", "eth_", colnames(eth_dummies))

fullData <- cbind(fullData, race_dummies, eth_dummies)
new_covar_names <- c(colnames(race_dummies), colnames(eth_dummies))
reference_cols  <- c(new_covar_names[1], colnames(eth_dummies)[1])
final_one_hot_covars <- setdiff(new_covar_names, reference_cols)

# --- 7. DATA CLEANING FOR PLINK ---
# Convert to data.table to enable .() notation and optimize handling
setDT(fullData)

# Re-incorporate the missing value fix to ensure PLINK data integrity
cols_to_fix <- c("CaseControl_PLINK", paste0("PC", 1:10), final_one_hot_covars)
for (col in cols_to_fix) { 
  if (col %in% names(fullData)) { 
    fullData[is.na(get(col)), (col) := -9] 
  } 
}

# --- 8. FINAL FILE OUTPUTS ---
# Using data.table syntax .() which now works after setDT()
write.table(fullData[, .(fid, iid, CaseControl_PLINK)], 
            file = paste0(dir3, run_mode, "_pheno.txt"), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(fullData[, .(fid, iid)], 
            file = paste0(dir3, run_mode, "_keep.txt"), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

covar_cols <- c("fid", "iid", paste0("PC", 1:10), final_one_hot_covars)
covar_cols <- intersect(covar_cols, names(fullData))

write.table(fullData[, ..covar_cols], 
            file = paste0(dir3, run_mode, "_covar.txt"), 
            col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

print(paste("All files successfully generated for", run_mode, "in:", dir3))
