##############################################################
# Script: 04_enrichment_analysis.R
# Purpose: Perform GO (enrichGO) and pathway (Enrichr) enrichment
#          analyses for differentially expressed proteins
##############################################################

# ========================
# Load required packages
# ========================
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(rrvgo)
library(enrichR)
library(openxlsx)

# ========================
# 1. GO Enrichment using enrichGO + rrvgo
# ========================
results <- read.csv("DE_results_limma_annotated.csv")

ontologies <- c("BP", "MF", "CC")
sig_groups <- c("Upregulated", "Downregulated")
all_results <- data.frame()

for (significance in sig_groups) {
  for (ontology in ontologies) {
    
    df <- results %>% filter(group == significance)
    
    go_result <- enrichGO(
      gene = unique(df$ENTREZID),
      OrgDb = org.Hs.eg.db,
      ont = ontology,
      pvalueCutoff = 0.05
    )
    
    if (is.null(go_result) || nrow(as.data.frame(go_result)) == 0) next
    
    go_df <- as.data.frame(go_result@result)
    go_ids <- go_df$ID
    go_scores <- setNames(-log10(go_df$p.adjust), go_ids)
    
    sim_matrix <- calculateSimMatrix(
      go_ids, orgdb = org.Hs.eg.db,
      ont = ontology, method = "Rel"
    )
    
    reduced_go <- reduceSimMatrix(
      sim_matrix, scores = go_scores,
      threshold = 0.7, orgdb = org.Hs.eg.db
    )
    
    go_reduced_df <- as.data.frame(reduced_go)
    
    go_merged <- go_df %>%
      left_join(go_reduced_df, by = c("ID" = "go")) %>%
      mutate(Ont = ontology, group = significance)
    
    all_results <- bind_rows(all_results, go_merged)
  }
}

# Fill in missing parent terms with Description
all_results_final <- all_results %>%
  mutate(parentTerm = ifelse(is.na(parentTerm), Description, parentTerm)) %>%
  mutate(
    GeneRatio_val = as.numeric(sapply(strsplit(GeneRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))),
    BgRatio_val = as.numeric(sapply(strsplit(BgRatio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))),
    ratio = GeneRatio_val / BgRatio_val
  )

# Save GO enrichment results
write.xlsx(all_results_final, file = "IMP_differential_GO.xlsx")

# ========================
# 2. Pathway Enrichment using Enrichr
# ========================
data <- read.csv("DE_results_limma_annotated.csv")
sig_groups <- c("Upregulated", "Downregulated")
dbs_to_use <- c("KEGG_2021_Human", "Reactome_2022", "WikiPathways_2024_Human")

all_enrichr_results <- list()

for (grp in sig_groups) {
  genes <- data %>% filter(group == grp) %>% pull(SYMBOL)
  
  if (length(genes) == 0) next
  
  enr_results <- enrichr(genes, dbs_to_use)
  
  combined <- bind_rows(
    lapply(names(enr_results), function(db_name) {
      db_df <- enr_results[[db_name]]
      if (!is.null(db_df) && nrow(db_df) > 0) {
        db_df %>%
          mutate(Database = db_name, Group = grp)
      } else {
        NULL
      }
    }),
    .id = NULL
  )
  
  all_enrichr_results[[grp]] <- combined
}

combined_results_all <- bind_rows(all_enrichr_results)

# Filter for significant terms
sig_combined_results <- combined_results_all %>%
  filter(Adjusted.P.value <= 0.05)

# Save significant pathway enrichment results
write.csv(sig_combined_results, "EnrichR_pathways_DE.csv", row.names = FALSE)
