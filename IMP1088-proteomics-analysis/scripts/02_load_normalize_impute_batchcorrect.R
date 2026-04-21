##############################################################
# Script: 02_load_normalize_batchcorrect.R
# Purpose: Load, clean, normalize, impute, and batch-correct
#          proteomics data for IMP-1088 analysis
##############################################################
##########################################################
## Load the files
##########################################################
list.files()
df <- read_tsv("20230616-Saber_rerun_Report_BGS Factory Report (Normal).tsv")
##########################################################
## remove decoy, filter by Qvalue and remove contaminants
##########################################################
d1 <- df[!df$EG.IsDecoy, ]
d1 <- d1[d1$PG.Qvalue < 0.01, ]
d1_conrev_remove <- d1[!grepl("^(CON__|REV__)", d1$PG.ProteinAccessions), ]
data <- d1
##########################################################
# Log2-transform PG.Quantity & Rename Conditions
##########################################################
table(is.na(data$PG.Quantity))
data_na <- data %>%
  filter(is.na(PG.Quantity))
all_data <- data %>%
  mutate(log2_PG.Quantity = log2(PG.Quantity + 1)) %>%
  mutate(Sample_Name = paste0(R.Condition, "_", R.Replicate))  # Rename as Condition_Replicate
##########################################################
# Median normalization by sample (R.FileName)
##########################################################
all_data <- all_data %>%
  group_by(Sample_Name) %>%
  mutate(median_shift = median(log2_PG.Quantity, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(norm_log2_PG.Quantity = log2_PG.Quantity - median_shift)
# Check for missing values
table(is.na(all_data$norm_log2_PG.Quantity))
# Check if there's still NA
na_rows <- all_data %>%
  filter(is.na(norm_log2_PG.Quantity))
View(na_rows)
##########################################################
# Create metadata with renamed sample names
##########################################################
sample_info <- all_data %>%
  select(Sample_Name, R.Condition, R.Replicate) %>%
  distinct()
write_csv(sample_info, "metadata.csv")
##########################################################
# Create wide format matrix using new sample names
##########################################################
wide_input <- all_data %>%
  group_by(PG.ProteinAccessions, Sample_Name) %>%
  summarize(norm_log2_quantity = median(norm_log2_PG.Quantity, na.rm = TRUE)) %>%
  ungroup()
# Pivot to wide format (proteins as rows, samples as columns)
wide_matrix <- wide_input %>%
  pivot_wider(names_from = Sample_Name, values_from = norm_log2_quantity) %>%
  distinct(PG.ProteinAccessions, .keep_all = TRUE) %>%
  column_to_rownames("PG.ProteinAccessions")
wide_matrix <- wide_matrix %>%
  select(Ctrl_1, Ctrl_2, Ctrl_3, Ctrl_4, Ctrl_5, 
         IMP_1, IMP_2, IMP_3, IMP_4, IMP_5)
##########################################################
# keeping those that satisfy ≥ 80% in any one group
##########################################################
# Define your experimental groups based on sample names
groupings <- list(
  CTRL = grep("^Ctrl", colnames(wide_matrix), value = TRUE),
  IMP = grep("^IMP", colnames(wide_matrix), value = TRUE))

# Helper function to check if a row has data in ≥ 80% of group samples
keep_protein <- function(protein_row, group_samples) {
  observed <- sum(!is.na(protein_row[group_samples]))
  required <- ceiling(0.8 * length(group_samples))
  return(observed >= required)
}

# Apply the filter across all rows, keeping those that satisfy ≥ 80% in any one group
keep_rows <- apply(wide_matrix, 1, function(row) {
  any(sapply(groupings, function(samples) keep_protein(row, samples)))
})
# Filtered matrix
filtered_matrix <- wide_matrix[keep_rows, ]
# Turn into matrix 
exprs_mat <- as.matrix(filtered_matrix)

##########################################################
# impute 
##########################################################
library(imputeLCMD)
dim(exprs_mat)  # should be something like 3000 x 30
# Apply QRILC
set.seed(123)  # for reproducibility

imputed <- impute.QRILC(exprs_mat)
#Extract the imputed matrix
imputed_matrix <- as.data.frame(imputed[[1]])
imputed_matrix <- imputed_matrix %>%
  mutate(Protein = rownames(exprs_mat)) %>%
  select(Protein, everything())
library(openxlsx)
write.xlsx(imputed_matrix, "imputed_matrix.xlsx")