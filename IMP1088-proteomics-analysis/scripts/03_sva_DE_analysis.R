##############################################################
# Script: 03_DE_analysis.R
# Purpose: SVA correction and differential expression analysis
#          using limma for IMP-1088 proteomics data
##############################################################

# ========================
# Load required packages
# ========================
library(tidyverse)
library(readxl)
library(sva)
library(limma)
library(clusterProfiler)
library(org.Hs.eg.db)

# ========================
# 1. Load and prepare data
# ========================
sample_info <- read.csv("metadata.csv")
imputed_matrix <- read_xlsx("imputed_matrix.xlsx")

# Extract expression matrix
expr_matrix <- imputed_matrix[, -1]
rownames(expr_matrix) <- imputed_matrix$Protein

# Align metadata to expression matrix
sample_info_ordered <- sample_info %>%
  filter(Sample_Name %in% colnames(expr_matrix)) %>%
  arrange(match(Sample_Name, colnames(expr_matrix)))
expr_matrix <- expr_matrix[, sample_info_ordered$Sample_Name]

# ========================
# 2. SVA batch correction
# ========================
mod <- model.matrix(~ R.Condition, data = sample_info_ordered)
mod0 <- model.matrix(~ 1, data = sample_info_ordered)

sva_result <- sva(as.matrix(expr_matrix), mod, mod0)

# Append surrogate variables to metadata
sva_df <- as.data.frame(sva_result$sv)
colnames(sva_df) <- paste0("SV", seq_len(ncol(sva_df)))
sample_info_sva <- cbind(sample_info_ordered, sva_df)

# PCA plot of surrogate variables
ggplot(sample_info_sva, aes(x = SV1, y = SV2, color = R.Condition)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "Surrogate Variables from SVA")

# Remove SV effects from expression matrix
clean_expr <- removeBatchEffect(expr_matrix, covariates = sva_result$sv, design = mod)

# PCA plot after batch correction
pca_res <- prcomp(t(clean_expr), scale. = TRUE)
pca_df <- as.data.frame(pca_res$x[, 1:2]) %>%
  mutate(Sample_Name = rownames(.),
         Group = sample_info_sva$R.Condition)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, label = Sample_Name)) +
  geom_point(size = 3) +
  geom_text(hjust = 1.2, vjust = 0.5, size = 3) +
  theme_minimal() +
  labs(title = "PCA after SVA Batch Correction")

# ========================
# 3. Annotate proteins with UniProt → Gene info
# ========================
# Extract primary UniProt ID
clean_protein_ids <- sapply(strsplit(imputed_matrix$Protein, ";"), `[`, 1)
clean_expr_df <- as.data.frame(clean_expr) %>%
  mutate(Protein = clean_protein_ids) %>%
  select(Protein, everything())

# Map to ENTREZ and SYMBOL
hci_bitr <- bitr(clean_expr_df$Protein, 
                 fromType = "UNIPROT", 
                 toType = c("ENTREZID", "SYMBOL"), 
                 OrgDb = org.Hs.eg.db)

# Handle multiple mappings by selecting first match
hci_bitr_unique <- hci_bitr %>%
  group_by(UNIPROT) %>%
  slice(1) %>%
  ungroup()

# Join annotations to expression matrix
annotated_expr <- clean_expr_df %>%
  left_join(hci_bitr_unique, by = c("Protein" = "UNIPROT")) %>%
  select(Protein, ENTREZID, SYMBOL, everything())

# Save annotated batch-corrected expression matrix
write.csv(annotated_expr, "sva_corrected_expression_matrix_annotated.csv", row.names = FALSE)

# ========================
# 4. Differential Expression Analysis using limma
# ========================
expr_matrix <- annotated_expr %>%
  select(-Protein, -ENTREZID, -SYMBOL) %>%
  as.matrix()

# Re-align metadata
sample_info_ordered <- sample_info %>%
  filter(Sample_Name %in% colnames(expr_matrix)) %>%
  arrange(match(Sample_Name, colnames(expr_matrix)))

# Design matrix (no intercept)
design <- model.matrix(~ 0 + R.Condition, data = sample_info_ordered)
colnames(design) <- gsub("R.Condition", "", colnames(design))

# Contrast: IMP vs Ctrl
contrast_matrix <- makeContrasts(
  IMP_vs_Ctrl = IMP - Ctrl,
  levels = design
)

# Fit model and compute statistics
fit <- lmFit(expr_matrix, design) %>%
  contrasts.fit(contrast_matrix) %>%
  eBayes()

# Extract DE results
IMP_DE <- topTable(fit, coef = "IMP_vs_Ctrl", number = Inf)

# Add annotation back
results <- IMP_DE %>%
  mutate(Row_ID = as.integer(rownames(IMP_DE))) %>%
  bind_cols(annotated_expr[match(rownames(IMP_DE), rownames(annotated_expr)), c("Protein", "ENTREZID", "SYMBOL")]) %>%
  mutate(group = case_when(
    adj.P.Val <= 0.05 & logFC >= 1 ~ "Upregulated",
    adj.P.Val <= 0.05 & logFC <= -1 ~ "Downregulated",
    TRUE ~ "Non-significant"
  )) %>%
  select(Protein, ENTREZID, SYMBOL, everything(), -Row_ID)

# DE summary
cat("Upregulated:", nrow(results %>% filter(group == "Upregulated")), "\n")
cat("Downregulated:", nrow(results %>% filter(group == "Downregulated")), "\n")

# Check NMTs
results %>% filter(str_detect(SYMBOL, "NMT"))

# Save DE results
write.csv(results, "DE_results_limma_annotated.csv", row.names = FALSE)

# ========================
# Citation
# ========================
# Tyanova S, Temu T, Cox J. The MaxQuant computational platform for mass spectrometry-based 
# shotgun proteomics. Nature Protocols (2016).
