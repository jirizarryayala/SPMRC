library(data.table)
icd_data <- fread("/lustre/project/crosslin/emerge/data/phenotype_data/icd_codes/output/outputfile1.csv")
icd9_target <- "^496|^491\\.21|^492|^492\\.0|^492\\.8" #Replace with your ICD codes
icd10_target <- "^J44\\.9|^J44\\.1|^J43\\.9|^J44\\.0" #Replace with your ICD codes
cases_subset <- icd_data[
  (grepl(icd9_target, ICD_CODE) & ICD_FLAG == 9) | 
    (grepl(icd10_target, ICD_CODE) & ICD_FLAG == 10)
]
case_ids <- unique(cases_subset$SUBJID)
all_ids <- unique(icd_data$SUBJID)
final_cohort <- data.frame(ID = all_ids)
final_cohort$Status <- 0
final_cohort$Status[final_cohort$ID %in% case_ids] <- 1
write.csv(final_cohort, 
          "/lustre/project/crosslin/crosslin_team/ssalter/SophiaSalter", 
          row.names = FALSE)
# Q2: How many ICD codes for your disease were present?
nrow(cases_subset)
# Q3: How many unique individuals (cases)?
length(case_ids)
# Q4: Table of cases/controls
table(final_cohort$Status)