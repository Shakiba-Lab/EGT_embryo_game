# Standalone script to calculate DEGs from hex_data
# Run this ONCE, then load the cached results in MAC_et_al_2026_scRNAseq_Analysis.Rmd
# Params: only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25

library(Seurat)
library(tidyverse)

# Set directories
data_dir   <- "../data/input_rna_data/"
output_dir <- "../output/"

cat("Loading and preparing data...\n")

# Load hex_data
hex_data_raw <- readRDS(file = file.path(data_dir, "D0toD5_WTonly.rds")) %>% UpdateSeuratObject()

# Rename clusters with biological annotations
hex_data <- hex_data_raw
hex_data$seurat_clusters <- dplyr::case_when(
  as.character(hex_data$seurat_clusters) == "3"  ~ "D3-D5 Epiblast",
  as.character(hex_data$seurat_clusters) == "9"  ~ "D0-D2 Epiblast",
  as.character(hex_data$seurat_clusters) == "10" ~ "Low Quality",
  as.character(hex_data$seurat_clusters) == "15" ~ "D3-D5 Amnion-like",
  as.character(hex_data$seurat_clusters) == "16" ~ "D3-D5 Posterior-like",
  TRUE ~ as.character(hex_data$seurat_clusters)
)
hex_data <- subset(hex_data, subset = seurat_clusters != "Low Quality")

# Calculate DEGs
cat("=== Calculating DEGs ===\n\n")

options(future.globals.maxSize = 5e9)  # 5 GiB for parallel workers

cat("[1/1] hex_data...\n")
degs_hex_data <- suppressWarnings(
  FindAllMarkers(
    hex_data,
    assay            = "SCT",
    group.by         = "seurat_clusters",
    only.pos         = TRUE,
    min.pct          = 0.25,
    logfc.threshold  = 0.25,
    verbose          = FALSE
  )
)
cat("✓", nrow(degs_hex_data), "genes across", length(unique(degs_hex_data$cluster)), "clusters\n\n")

# Compile DEGs into list
DEGs <- list(
  hex_data = degs_hex_data
)

# Summary
cat("=== DEG Summary ===\n")
for (name in names(DEGs)) {
  cat(sprintf("%-20s %d genes across %d clusters\n",
    name,
    nrow(DEGs[[name]]),
    length(unique(DEGs[[name]]$cluster))
  ))
}

# Save to RDS
data_output_dir <- file.path(output_dir, "data")
if (!dir.exists(data_output_dir)) {
  dir.create(data_output_dir, recursive = TRUE)
}

degs_file <- file.path(data_output_dir, "degs_results.rds")
saveRDS(DEGs, file = degs_file)
cat("\n✓ DEGs saved to:", degs_file, "\n")
