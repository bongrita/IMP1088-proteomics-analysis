# ================================
# Load or Install All Libraries
# ================================
packages <- c(
  "readxl","conflicted", "readr", "tidyr", "tibble", "pheatmap", "ggbreak",
  "imputeLCMD", "openxlsx", "limma", "writexl", "showtext", "jsonlite", "curl",
  "ggplot2", "scales", "ggrepel", "sva", "biomaRt", "org.Hs.eg.db", 
  "clusterProfiler", "GOSemSim", "rrvgo", "enrichR", "Cairo",
  "ComplexHeatmap", "circlize", "RColorBrewer", "forcats", "dplyr"
)

install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

invisible(lapply(packages, install_if_missing))
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
#######################################################################################