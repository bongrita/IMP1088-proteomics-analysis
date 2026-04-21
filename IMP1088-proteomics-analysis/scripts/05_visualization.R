


##############################################################
# Script: 05_visualization.R
# Purpose: Generate heatmap
##############################################################

# ========================
# Load libraries
# ========================
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(Cairo)

conflicts_prefer(dplyr::slice)
df <- read.csv("EnrichR_pathways_DE.csv")
react <- df %>%
  filter(Database %in% "Reactome_2022")  %>%
  separate(Term, sep = " R-", into = c("Term", "id"))

# ========================
# 1. Define relevant pathways
# ========================
relevant_pathways <- c(
  "Unfolded Protein Response (UPR)",
  "COPI-dependent Golgi-to-ER Retrograde Traffic",
  "Golgi Cisternae Pericentriolar Stack Reorganization",
  "Diseases Of Glycosylation",
  "Post-translational Protein Modification",
  "Lipophagy",
  "Extracellular Matrix Organization"
)

# Extract gene-to-pathway mapping (one term per gene)
er <- react %>%
  filter(Term %in% relevant_pathways) %>%
  separate_rows(Genes, sep = ";") %>%
  group_by(Genes) %>%
  arrange(Adjusted.P.value) %>%
  slice(1) %>%
  ungroup() %>%
  select(Term, Genes)
gene <- er$Genes
pathway <- er$Term

# Term                                                id         
# <chr>                                               <chr>      
#   1 Extracellular Matrix Organization                   HSA-1474244
# 2 Diseases Of Glycosylation                           HSA-3781865
# 3 Post-translational Protein Modification             HSA-597592 
# 4 COPI-dependent Golgi-to-ER Retrograde Traffic       HSA-6811434
# 5 Unfolded Protein Response (UPR)                     HSA-381119 
# 6 Golgi Cisternae Pericentriolar Stack Reorganization HSA-162658 
# 7 Lipophagy                                           HSA-9613354


# ========================
# 2. Prepare expression matrix and z-score transform
# ========================
mat <- read.csv("sva_corrected_expression_matrix_annotated.csv")
expr_mat <- mat %>%
  filter(SYMBOL %in% gene) %>%
  distinct(SYMBOL, .keep_all = TRUE)

expr_mat_named <- expr_mat %>%
  select(SYMBOL, starts_with("Ctrl"), starts_with("IMP")) %>%
  column_to_rownames("SYMBOL")

# Z-score by gene
zscore_mat <- t(scale(t(as.matrix(expr_mat_named))))
zscore_mat <- t(zscore_mat)  # rows = samples, columns = genes

# # Collapse replicates by condition
# ctrl_avg <- rowMeans(zscore_mat[, grep("^Ctrl", colnames(zscore_mat))])
# imp_avg <- rowMeans(zscore_mat[, grep("^IMP", colnames(zscore_mat))])
# collapsed_mat <- rbind(CTRL = ctrl_avg, IMP = imp_avg)

# ========================
# 3. Column annotation (Pathway per gene)
# ========================
column_annot <- data.frame(
  SYMBOL = colnames(zscore_mat),
  Pathway = pathway[match(colnames(zscore_mat), gene)]
)

# Recode pathway names for plotting
column_annot$Pathway <- recode(column_annot$Pathway,
                               "COPI-dependent Golgi-to-ER Retrograde Traffic" = "COPI-dependent\nGolgi-to-ER\nRetrograde Traffic",
                               "Golgi Cisternae Pericentriolar Stack Reorganization" = "Golgi Cisternae\nPericentriolar \nStack\nReorganization",
                               "Diseases Of Glycosylation" = "Diseases of\nGlycosylation",
                               "Post-translational Protein Modification" = "Post-translational\nProtein Modification",
                               "Extracellular Matrix Organization" = "Extracellular Matrix\nOrganization",
                               "Unfolded Protein Response (UPR)" = "Unfolded Protein\nResponse (UPR)"
)

# Define pathway order
desired_order <- c(
  "COPI-dependent\nGolgi-to-ER\nRetrograde Traffic",
  "Golgi Cisternae\nPericentriolar \nStack\nReorganization",
  "Lipophagy",
  "Diseases of\nGlycosylation",
  "Post-translational\nProtein Modification",
  "Extracellular Matrix\nOrganization",
  "Unfolded Protein\nResponse (UPR)"
)

column_annot$Pathway <- factor(column_annot$Pathway, levels = desired_order, ordered = TRUE)

# Reorder genes based on pathway
gene_order <- order(column_annot$Pathway)
zscore_mat <- zscore_mat[, gene_order]
column_annot <- column_annot[gene_order, , drop = FALSE]

# ========================
# 4. Top annotation (Pathway)
# ========================
color_palette <- brewer.pal(length(levels(column_annot$Pathway)), "Set3")
top_ha <- HeatmapAnnotation(
  Pathway = column_annot$Pathway,
  gp = gpar(fontsize = rel(26), fontfamily = "Times", col = "black"),
  col = list(Pathway = setNames(color_palette, levels(column_annot$Pathway))),
  show_annotation_name = TRUE,
  annotation_name_side = "left",
  show_legend = FALSE,
  annotation_name_gp = gpar(fontsize = rel(20), fontface = "bold", fontfamily = "Times")
)

# ========================
# 5. Row annotation (Sample condition)
# ========================
# Define sample group with descriptive labels
sample_group <- ifelse(grepl("^Ctrl", rownames(zscore_mat)), "DMSO", "IMP-1088")
sample_group <- factor(sample_group, levels = c("DMSO", "IMP-1088"))



# ========================
# 6. Plot heatmap
# ========================
col_fun <- colorRamp2(seq(-2, 2, length.out = 11), rev(brewer.pal(11, "RdYlBu")))

# ===============================
# 🎨 Font Setup for Publication
# ===============================
font_add("Times", regular = "C:/Windows/Fonts/times.ttf")
showtext_auto()

CairoPNG("Heatmap_pathways.png", width = 2600, height = 500)

# CairoPDF("Heatmap_pathways.pdf", width = 26, height = 8)
Heatmap(
  zscore_mat,
  name = "Z-score",
  row_split = sample_group,
  col = col_fun,
  top_annotation = top_ha,
  # left_annotation = row_ha,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  column_split = column_annot$Pathway,
  column_gap = unit(10, "mm"),
  show_row_names = FALSE,
  show_column_names = TRUE,
  # width = unit(60, "cm"),
  # height = unit(8, "cm"),
  row_title_gp = gpar(fontsize = rel(20), fontfamily = "Times"), # this is DMSO and IMP
  column_names_gp = gpar(fontsize = rel(20), fontfamily = "Times"), # this is the gene name at the bottom
  column_title_gp = gpar(fontsize = rel(20), fontfamily = "Times"), # this is the pathway names on the top
  heatmap_legend_param = list(
    legend_height = unit(6, "cm"),
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 10)
  )
)
dev.off()


#
# ========================
# 1. Define pathway and extract genes
# ========================
# Load proteomics expression matrix
matrix <- read.csv("sva_corrected_expression_matrix_annotated.csv")

# Define databases to use for enrichment
dbs_to_use <- c("KEGG_2021_Human", "Reactome_2022", "WikiPathways_2024_Human")

# Pull all gene symbols from expression matrix
genes <- matrix %>% dplyr::pull(SYMBOL)

# Exit if no genes found
if (length(genes) == 0) stop("No gene symbols found in matrix")

# Run enrichment
results <- enrichr(genes, dbs_to_use)

# Combine enrichment results into one dataframe
combined_results <- bind_rows(
  lapply(names(results), function(db_name) {
    if (!is.null(results[[db_name]]) && nrow(results[[db_name]]) > 0) {
      results[[db_name]] %>% mutate(Database = db_name)
    } else {
      NULL
    }
  }),
  .id = NULL
)

# ========================
# Extract SARS-CoV-2–related terms
# ========================
sars <- combined_results %>%
  filter(Database == "WikiPathways_2024_Human", grepl("SARS", Term)) %>%
  separate_rows(Genes, sep = ";") %>%
  group_by(Genes) %>%
  arrange(Adjusted.P.value) %>%
  slice(1) %>%
  ungroup()

# ========================
# Isolate genes from specific SARS-CoV-2 autophagy pathway
# ========================
virus_pathway <- "Perturbations Host Cell Autophagy Induced By SARS CoV 2 Prots WP4936"

View(virus)
virus <- sars %>%
  filter(Term == virus_pathway) %>%
  select(Term, Genes) %>%
  separate_rows(Genes, sep = ";") %>%
  distinct()

# Output clean gene and pathway lists
gene <- virus$Genes
pathway <- virus$Term
# ========================
# 2. Prepare expression matrix and Z-score transform
# ========================
mat <- read.csv("sva_corrected_expression_matrix_annotated.csv")

expr_mat <- mat %>%
  filter(SYMBOL %in% gene) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  select(SYMBOL, starts_with("Ctrl"), starts_with("IMP")) %>%
  column_to_rownames("SYMBOL")

# Remove genes with zero variance to avoid NaN in z-score
expr_mat_clean <- expr_mat[rowSds(as.matrix(expr_mat)) != 0, ]

# Z-score transform
zscore_mat <- t(scale(t(as.matrix(expr_mat_clean))))
zscore_mat <- t(zscore_mat)

# ========================
# 3. Column annotation (Pathway per gene)
# ========================
column_annot <- data.frame(
  SYMBOL = colnames(zscore_mat),
  Pathway = rep("Perturbation of Host Cell Autophagy (WP4936)", ncol(zscore_mat))
)
rownames(column_annot) <- column_annot$SYMBOL
column_annot$Pathway <- factor(column_annot$Pathway)

# ========================
# 4. Top annotation
# ========================
top_ha <- HeatmapAnnotation(
  Pathway = column_annot$Pathway,
  col = list(Pathway = setNames(brewer.pal(3, "Set3")[1], levels(column_annot$Pathway))),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = rel(20), fontface = "bold", fontfamily = "Times"),
  show_annotation_name = TRUE,
  show_legend = FALSE
)

# ========================
# 5. Row annotation (sample group)
# ========================
sample_group <- ifelse(grepl("^Ctrl", rownames(zscore_mat)), "DMSO", "IMP-1088")
sample_group <- factor(sample_group, levels = c("DMSO", "IMP-1088"))

# ========================
# 6. Plot heatmap
# ========================
col_fun <- colorRamp2(seq(-2, 2, length.out = 11), rev(brewer.pal(11, "RdYlBu")))

font_add("Times", regular = "C:/Windows/Fonts/times.ttf")
showtext_auto()

CairoPDF("Heatmap_SARS_perturbation_final.pdf", width = 26, height = 8)

Heatmap(
  zscore_mat,
  name = "Z-score",
  row_split = sample_group,
  col = col_fun,
  top_annotation = top_ha,
  cluster_rows = FALSE,
  cluster_columns = TRUE,
  column_gap = unit(10, "mm"),
  show_row_names = FALSE,
  show_column_names = TRUE,
  row_title_gp = gpar(fontsize = rel(20), fontfamily = "Times"),
  column_names_gp = gpar(fontsize = rel(20), fontfamily = "Times"),
  column_title_gp = gpar(fontsize = rel(20), fontfamily = "Times"),
  heatmap_legend_param = list(
    legend_height = unit(6, "cm"),
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 10)
  )
)

dev.off()
sessionInfo()
