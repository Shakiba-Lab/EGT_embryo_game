# ============================================================================
# Abou Chakra et al. 2026 scRNA-seq Analysis
# Complete pipeline: data prep, visualization, enrichment, coexpression
# ============================================================================

# === Setup & Data Preparation ===

# Load required libraries
suppressMessages({
  library(tidyverse)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(ggridges)
  library(ggtext)
  library(parallel)
})

# Set directories
data_dir <- "../data/input_rna_data/"
output_dir <- "../output/"

# Initialize parameters list
params <- list()

cat("Loading hex_data Seurat object...\n")
hex_data_raw <- readRDS(file = file.path(data_dir, "D0toD5_WTonly.rds")) %>% UpdateSeuratObject()
cat("  ✓ Loaded hex_data:", ncol(hex_data_raw), "cells ×", nrow(hex_data_raw), "genes\n")

cat("Loading Lima et al supplementary data...\n")
winner_epi_lima <- readxl::read_excel("../data/lima_et_al_data/EMS127075-supplement-Supplementary_Tables.xlsx", sheet = 1)
loser_epi_lima <- readxl::read_excel("../data/lima_et_al_data/EMS127075-supplement-Supplementary_Tables.xlsx", sheet = 2)
cat("  ✓ Loaded winner_epi:", nrow(winner_epi_lima), "genes\n")
cat("  ✓ Loaded loser_epi:", nrow(loser_epi_lima), "genes\n")

# Assign datasets for analysis
hex_data <- hex_data_raw

# Rename hex_data clusters with biological annotations
hex_data$seurat_clusters <- dplyr::case_when(
  as.character(hex_data$seurat_clusters) == "3" ~ "D3-D5 Epiblast",
  as.character(hex_data$seurat_clusters) == "9" ~ "D0-D2 Epiblast",
  as.character(hex_data$seurat_clusters) == "10" ~ "Low Quality",
  as.character(hex_data$seurat_clusters) == "15" ~ "D3-D5 Amnion-like",
  as.character(hex_data$seurat_clusters) == "16" ~ "D3-D5 Posterior-like",
  TRUE ~ as.character(hex_data$seurat_clusters)
)

# Filter out Low Quality cluster
cat("Filtering out Low Quality cells...\n")
cells_before <- ncol(hex_data)
hex_data <- subset(hex_data, subset = seurat_clusters != "Low Quality")
cells_after <- ncol(hex_data)
cat("  ✓ Removed", cells_before - cells_after, "Low Quality cells\n")
cat("  ✓ Remaining:", cells_after, "cells\n")
cat("  ✓ Clusters:", paste(unique(hex_data$seurat_clusters), collapse = ", "), "\n")

# Organize datasets into params list
params$data$hex_data   <- hex_data

# Datasets used for visualization (gene expression, coexpression plots)
params$viz_datasets <- c("hex_data")

# Store Lima et al supplementary data (winner/loser epiblast only)
params$lima_et_al$win_epi <- winner_epi_lima
params$lima_et_al$lose_epi <- loser_epi_lima

# Define genes to analyze
params$genes$to_plot <- c("NANOG", "POU5F1", "ISL1", "TFAP2A", "MIXL1", "TBXT", "MYC")

# Build gene mapping and validate presence in visualization datasets
params$genes$mapping <- list()
params$genes$final <- params$genes$to_plot

for (dataset_name in params$viz_datasets) {
  dataset <- params$data[[dataset_name]]
  available_genes <- rownames(dataset)
  mapped_genes <- params$genes$to_plot

  # Handle TBXT ↔ T conversion
  if ("TBXT" %in% params$genes$to_plot && !("TBXT" %in% available_genes) && "T" %in% available_genes) {
    mapped_genes[mapped_genes == "TBXT"] <- "T"
  }

  params$genes$mapping[[dataset_name]] <- mapped_genes
  params$genes$final <- params$genes$final[mapped_genes %in% available_genes]
}

# Keep only essential objects (everything is organized in params)
rm(list = setdiff(ls(), c("params", "output_dir")))
gc()  # Force garbage collection to free memory

# === Initialize Figures Object ===

# Create figures object to store all publication-ready plots
# Structure: figures$FIGURE_NUMBER$PANEL_LETTER
# Example: figures$1$A, figures$1$B, figures$2$C, etc.

figures <- list()

cat("✓ Figures object initialized\n")
cat("  Store plots as: figures$1$A, figures$1$B, figures$2$A, etc.\n")

# ============================================================================
# VISUALIZATION: HEX DATA UMAP
# ============================================================================

# === Hex Clusters UMAP (S4a) ===

hex_data_obj <- params$data$hex_data

# Extract UMAP coordinates and metadata
umap_coords_hex <- Embeddings(hex_data_obj, reduction = "umap") %>%
  as.data.frame() %>%
  rownames_to_column("cell_id") %>%
  rename(umap_1 = UMAP_1, umap_2 = UMAP_2)

metadata_hex <- hex_data_obj@meta.data %>%
  rownames_to_column("cell_id") %>%
  select(cell_id, seurat_clusters)

# Define hex clusters to highlight
hex_named_clusters <- c("D0-D2 Epiblast", "D3-D5 Epiblast", "D3-D5 Amnion-like", "D3-D5 Posterior-like")

# Create plotting data frame
plot_data_hex <- umap_coords_hex %>%
  left_join(metadata_hex, by = "cell_id") %>%
  mutate(
    cell_group = case_when(
      seurat_clusters %in% hex_named_clusters ~ seurat_clusters,
      TRUE ~ "Other"
    )
  ) %>%
  # Sort so "Other" is plotted first (background)
  arrange(desc(cell_group))

# Create plot with layered approach: Other cells first (background), then colored clusters
p_hex_umap <- ggplot(plot_data_hex, aes(x = umap_1, y = umap_2)) +
  # Colored layer: clusters of interest
  geom_point(
    data = plot_data_hex %>% filter(cell_group != "Other"),
    aes(color = cell_group),
    size = 0.5,
    alpha = 0.8
  ) +
  scale_color_manual(
    values = c(
      "D0-D2 Epiblast" = "#6BB1E4",
      "D3-D5 Epiblast" = "#E5D98B",
      "D3-D5 Amnion-like" = "#AA70AF",
      "D3-D5 Posterior-like" = "#C4789A"
    )
  ) +
  labs(
    x = "UMAP 1",
    y = "UMAP 2",
    title = "Hex-Embryoid Data",
    color = "Cluster"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    aspect.ratio = 1
  )

# Store in figures object (Supplementary Figure 4, Panel A)
figures$S4$a <- p_hex_umap

cat("✓ Stored in figures$S4$a (Hex Data UMAP)\n")

# === Gene Expression UMAPs - Panel B (NANOG, MYC, POU5F1) ===

dataset    <- params$data$hex_data
gene_names <- c("NANOG", "MYC", "POU5F1")

# Handle TBXT <-> T mapping for hex_data
available_genes <- rownames(dataset)
mapped_genes <- setNames(gene_names, gene_names)

figures$S4$b <- list()

for (gene in gene_names) {
  gene_plot <- mapped_genes[[gene]]

  p <- FeaturePlot(dataset, features = gene_plot, reduction = "umap",
                   pt.size = 0.8, order = TRUE, label = FALSE) &
    labs(x = "UMAP 1", y = "UMAP 2", title = gene) &
    theme_classic() +
    theme(
      plot.title  = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title  = element_text(face = "bold"),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks  = element_blank(),
      aspect.ratio = 1
    )

  figures$S4$b[[gene]] <- p
  cat("✓ Stored in figures$S4$b$", gene, "\n", sep = "")
}

# === Gene Expression UMAPs - Panel C (ISL1, TFAP2A) ===

dataset    <- params$data$hex_data
gene_names <- c("ISL1", "TFAP2A")

figures$S4$c <- list()

for (gene in gene_names) {
  p <- FeaturePlot(dataset, features = gene, reduction = "umap",
                   pt.size = 0.8, order = TRUE, label = FALSE) &
    labs(x = "UMAP 1", y = "UMAP 2", title = gene) &
    theme_classic() +
    theme(
      plot.title  = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title  = element_text(face = "bold"),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks  = element_blank(),
      aspect.ratio = 1
    )

  figures$S4$c[[gene]] <- p
  cat("✓ Stored in figures$S4$c$", gene, "\n", sep = "")
}

# === Gene Expression UMAPs - Panel D (MIXL1, TBXT) ===

dataset    <- params$data$hex_data
gene_names <- c("MIXL1", "TBXT")

# Handle TBXT <-> T: hex_data may store it as "T"
available_genes <- rownames(dataset)
gene_map <- setNames(gene_names, gene_names)
if (!"TBXT" %in% available_genes && "T" %in% available_genes) {
  gene_map[["TBXT"]] <- "T"
}

figures$S4$d <- list()

for (gene in gene_names) {
  gene_plot <- gene_map[[gene]]

  p <- FeaturePlot(dataset, features = gene_plot, reduction = "umap",
                   pt.size = 0.8, order = TRUE, label = FALSE) &
    labs(x = "UMAP 1", y = "UMAP 2", title = gene) &
    theme_classic() +
    theme(
      plot.title  = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title  = element_text(face = "bold"),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks  = element_blank(),
      aspect.ratio = 1
    )

  figures$S4$d[[gene]] <- p
  cat("✓ Stored in figures$S4$d$", gene, "\n", sep = "")
}

# ============================================================================
# ENRICHMENT ANALYSIS
# ============================================================================

# === DEG Calculation ===

# Run calculate_degs.R to compute DEGs and save to degs_results.rds
# Then load result into params$DEGs for downstream enrichment analysis

cat("\n--- Calculating Differentially Expressed Genes ---\n")
source("calculate_degs.R")
params$DEGs <- DEGs

cat("\n✓ DEG Summary:\n")
for (name in names(DEGs)) {
  n_genes <- nrow(DEGs[[name]])
  n_clusters <- length(unique(DEGs[[name]]$cluster))
  cat(sprintf("  %s: %d genes across %d clusters\n", name, n_genes, n_clusters))
  cat("  Sample DEGs (first 5):\n")
  print(head(DEGs[[name]][, c("cluster", "gene", "avg_log2FC", "p_val_adj")], 5))
}

# === Lima Ortholog Conversion ===

# Convert Lima et al mouse genes to human orthologs (winner/loser epiblast only)
params$lima_orthologs <- list()

params$lima_orthologs$win_epi <- orthogene::convert_orthologs(
  gene_df          = params$lima_et_al$win_epi,
  gene_input       = "Gene ID",
  gene_output      = "column",
  input_species    = "mouse",
  output_species   = "human",
  non121_strategy  = "drop_both_species"
)

params$lima_orthologs$lose_epi <- orthogene::convert_orthologs(
  gene_df          = params$lima_et_al$lose_epi,
  gene_input       = "Gene ID",
  gene_output      = "column",
  input_species    = "mouse",
  output_species   = "human",
  non121_strategy  = "drop_both_species"
)

cat("✓ Human orthologs generated for Lima et al win_epi and lose_epi\n")

# === Hypergeometric Enrichment Testing ===

# hex_data only — background genes = genes expressed in each cluster's cells

n_cores <- max(1L, detectCores() - 1L)

data_obj    <- params$data$hex_data
cluster_col <- "seurat_clusters"
assay_use   <- "SCT"
degs_data   <- DEGs[["hex_data"]]

valid_clusters <- unique(data_obj[[cluster_col, drop = TRUE]])
degs_data      <- degs_data %>% filter(cluster %in% valid_clusters)
cluster_ids    <- unique(degs_data$cluster)

cat("\n--- Running Hypergeometric Enrichment ---\n")
cat(sprintf("Testing %d clusters against 2 Lima signatures (win_epi, lose_epi)\n", length(cluster_ids)))
cat("Clusters to test:", paste(cluster_ids, collapse = ", "), "\n")

cluster_results <- mclapply(cluster_ids, function(cluster_id,
                                                   data_obj, cluster_col, assay_use,
                                                   degs_data, params) {
  cluster_cells    <- Cells(data_obj)[data_obj[[cluster_col, drop = TRUE]] == cluster_id]
  counts_mat       <- GetAssayData(data_obj, assay = assay_use, layer = "counts")[, cluster_cells, drop = FALSE]
  background_genes <- rownames(counts_mat)[rowSums(counts_mat) > 0]

  cluster_degs <- degs_data %>% filter(cluster == cluster_id)
  K <- nrow(cluster_degs)
  N <- length(background_genes)

  lima_rows <- lapply(names(params$lima_orthologs), function(lima_set) {
    ortholog_gene_df <- params$lima_orthologs[[lima_set]]
    M <- sum(ortholog_gene_df$ortholog_gene %in% background_genes)
    overlap_genes <- ortholog_gene_df$ortholog_gene[
      ortholog_gene_df$ortholog_gene %in% cluster_degs$gene
    ]
    k <- length(overlap_genes)
    p_value <- phyper(k - 1, M, N - M, K, lower.tail = FALSE)

    tibble(
      dataset_base  = "hex_data",
      cluster       = as.character(cluster_id),
      lima_set      = lima_set,
      k = k, M = M, K = K, N = N,
      p_val         = p_value,
      overlap_genes = list(overlap_genes)
    )
  })

  bind_rows(lima_rows)
}, mc.cores = n_cores,
   data_obj = data_obj, cluster_col = cluster_col, assay_use = assay_use,
   degs_data = degs_data, params = params)

# Combine and apply Bonferroni correction
enrichment_results <- bind_rows(cluster_results) %>%
  mutate(
    n_tests         = n_distinct(cluster),
    p_val_corrected = pmin(p_val * n_tests, 1)
  )

params$enrichment_results <- enrichment_results

cat("✓ Enrichment complete:", nrow(enrichment_results), "tests\n")

# Show sample enrichment results
cat("\nEnrichment results summary:\n")
cat("  Significant hits (p_val_corrected < 0.05):", sum(enrichment_results$p_val_corrected < 0.05), "\n")
cat("  Sample results (top 6 by p-value):\n")
top_enrich <- enrichment_results %>% arrange(p_val_corrected) %>% slice_head(n = 6)
print(top_enrich[, c("cluster", "lima_set", "k", "M", "K", "p_val_corrected")])

# Create data output directory if it doesn't exist
data_output_dir <- file.path(output_dir, "data")
if (!dir.exists(data_output_dir)) {
  dir.create(data_output_dir, recursive = TRUE)
}

enrich_file <- file.path(data_output_dir, "enrichment_results.rds")
saveRDS(enrichment_results, enrich_file)
cat("\n✓ Enrichment results saved to:", enrich_file, "\n")

# === Enrichment Heatmap - Hex Data (S4e) ===

p_display_cap <- 50

# Cluster ordering: named clusters first, any others sorted alphabetically after
hex_named_clusters <- c("D0-D2 Epiblast", "D3-D5 Epiblast", "D3-D5 Amnion-like", "D3-D5 Posterior-like")
hex_other_clusters <- params$enrichment_results %>%
  filter(!cluster %in% hex_named_clusters) %>%
  pull(cluster) %>% unique() %>% sort()
hex_cluster_levels <- c(hex_named_clusters, hex_other_clusters)

# Prepare display data
enrich_data <- params$enrichment_results %>%
  mutate(
    p_val_display = pmin(pmax(0, -log2(pmin(p_val_corrected, 1))), p_display_cap),
    lima_set = recode(lima_set,
      win_epi  = "Winner\nEpiblast",
      lose_epi = "Loser\nEpiblast"
    )
  )

p_enrich_hex <- enrich_data %>%
  mutate(cluster = factor(cluster, levels = hex_cluster_levels)) %>%
  ggplot(aes(x = cluster, y = lima_set, fill = p_val_display)) +
  geom_tile(linewidth = 0.6) +
  scale_fill_gradientn(
    colors = c("white", "white", "red"),
    values = c(0, 4.32 / p_display_cap, 1),
    limits = c(0, p_display_cap),
    breaks = c(0, 25, 50),
    labels = c("0", "25", "50+")
  ) +
  coord_fixed(ratio = 1) +
  theme_classic() +
  labs(
    title = "Hypergeometric Enrichment - hex-Embryoid Data",
    x     = "",
    y     = "Gene signature\nfrom Lima et al. 2021\nNat. Metabolism",
    fill  = expression(-log[2](italic(p[value])))
  ) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.text.x  = element_text(angle = 45, hjust = 1),
    axis.title.y = element_text(size = 11),
    axis.title.x = element_blank(),
    axis.line    = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
    legend.title = element_text(face = "bold")
  )

figures$S4$e <- p_enrich_hex

cat("✓ Stored in figures$S4$e\n")

# ============================================================================
# COEXPRESSION ANALYSIS
# ============================================================================

# === Coexpression Setup ===

# 4-corner color scheme (gene1 = blue, gene2 = red, both = magenta)
color_ll <- "gray85"
color_hl <- "blue"
color_lh <- "red"
color_hh <- "magenta"

create_4corner_color <- function(val1, val2, color_ll, color_hl, color_lh, color_hh) {
  rgb_ll <- col2rgb(color_ll) / 255
  rgb_hl <- col2rgb(color_hl) / 255
  rgb_lh <- col2rgb(color_lh) / 255
  rgb_hh <- col2rgb(color_hh) / 255

  r <- (1-val1)*(1-val2)*rgb_ll[1] + val1*(1-val2)*rgb_hl[1] +
       (1-val1)*val2*rgb_lh[1]     + val1*val2*rgb_hh[1]
  g <- (1-val1)*(1-val2)*rgb_ll[2] + val1*(1-val2)*rgb_hl[2] +
       (1-val1)*val2*rgb_lh[2]     + val1*val2*rgb_hh[2]
  b <- (1-val1)*(1-val2)*rgb_ll[3] + val1*(1-val2)*rgb_hl[3] +
       (1-val1)*val2*rgb_lh[3]     + val1*val2*rgb_hh[3]

  rgb(pmin(1, pmax(0, r)), pmin(1, pmax(0, g)), pmin(1, pmax(0, b)))
}

# Helper: generate one coexpression UMAP + legend composite for a gene1 × gene2 pair
make_coexp_plot <- function(g1, g2, umap_coords, data_obj,
                             color_ll, color_hl, color_lh, color_hh) {
  available <- rownames(data_obj)
  g1_mapped <- if (g1 == "TBXT" && !"TBXT" %in% available && "T" %in% available) "T" else g1
  g2_mapped <- if (g2 == "TBXT" && !"TBXT" %in% available && "T" %in% available) "T" else g2

  if (!g1_mapped %in% available || !g2_mapped %in% available) return(NULL)

  expr_data <- FetchData(data_obj, vars = c(g1_mapped, g2_mapped)) %>%
    rownames_to_column("cell_id")

  plot_data <- umap_coords %>%
    left_join(expr_data, by = "cell_id") %>%
    mutate(
      gene1_norm      = (!!sym(g1_mapped) - min(!!sym(g1_mapped), na.rm = TRUE)) /
                        (max(!!sym(g1_mapped), na.rm = TRUE) - min(!!sym(g1_mapped), na.rm = TRUE)),
      gene2_norm      = (!!sym(g2_mapped) - min(!!sym(g2_mapped), na.rm = TRUE)) /
                        (max(!!sym(g2_mapped), na.rm = TRUE) - min(!!sym(g2_mapped), na.rm = TRUE)),
      blend_color     = create_4corner_color(gene1_norm, gene2_norm,
                                             color_ll, color_hl, color_lh, color_hh),
      blend_intensity = gene1_norm + gene2_norm
    ) %>%
    arrange(blend_intensity)

  title_text <- paste0(
    "<span style='color:", color_hl, "'><b>", g1, "</b></span> \u00d7 ",
    "<span style='color:", color_lh, "'><b>", g2, "</b></span>"
  )

  p_umap <- ggplot(plot_data, aes(x = umap_1, y = umap_2, color = blend_color)) +
    geom_point(size = 0.8, alpha = 0.8) +
    scale_color_identity() +
    theme_classic() +
    labs(title = title_text, x = "UMAP 1", y = "UMAP 2") +
    theme(
      plot.title   = element_markdown(face = "bold", size = 12, hjust = 0.5),
      axis.title   = element_text(face = "bold"),
      axis.text.x  = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks   = element_blank(),
      aspect.ratio = 1
    )

  legend_grid <- expand.grid(g1v = seq(0, 1, length.out = 12),
                              g2v = seq(0, 1, length.out = 12)) %>%
    mutate(blend_color = create_4corner_color(g1v, g2v,
                                              color_ll, color_hl, color_lh, color_hh))

  p_legend <- ggplot(legend_grid, aes(x = g2v, y = g1v, fill = blend_color)) +
    geom_tile(linewidth = 0, color = NA) +
    scale_fill_identity() +
    scale_x_continuous(name = g2, breaks = c(0, 1),
                       labels = c("Low", "High"), expand = c(0, 0)) +
    scale_y_continuous(name = g1, breaks = c(0, 1),
                       labels = c("Low", "High"), expand = c(0, 0)) +
    theme_minimal() +
    theme(
      axis.title   = element_text(face = "bold", size = 8),
      axis.text    = element_text(size = 7),
      aspect.ratio = 1,
      plot.margin  = margin(2, 2, 2, 2, "pt")
    )

  p_umap | p_legend + plot_layout(widths = c(3, 1))
}

# Shared UMAP coordinates for all coexpression plots
hex_data_coexp <- params$data$hex_data

umap_coords_coexp <- Embeddings(hex_data_coexp, reduction = "umap") %>%
  as.data.frame() %>%
  rownames_to_column("cell_id")
colnames(umap_coords_coexp)[2:3] <- c("umap_1", "umap_2")

coexp_gene1 <- c("ISL1", "MIXL1")

cat("\n--- Coexpression Analysis Setup ---\n")
cat("✓ UMAP coordinates loaded:", nrow(umap_coords_coexp), "cells\n")
cat("✓ Gene1 (blue channel):", paste(coexp_gene1, collapse = ", "), "\n")

# === Gene2 Lists ===

coexp_winner_genes <- c("FGF8", "CDA", "CNTNAP2", "IFITM1", "TUBB2B")
coexp_loser_genes  <- c("KRT19", "IER5L", "BAMBI", "DLL1", "EZR")

cat("✓ Gene2 lists defined:\n")
cat("  Winners (red channel):", paste(coexp_winner_genes, collapse = ", "), "\n")
cat("  Losers  (red channel):", paste(coexp_loser_genes, collapse = ", "), "\n")
cat("  Total coexpression plots to generate:", length(coexp_gene1) * (length(coexp_winner_genes) + length(coexp_loser_genes)), "\n")

# === Coexpression UMAPs - S5A (Winner Genes) ===

figures$S5$A <- list()

for (g1 in coexp_gene1) {
  figures$S5$A[[g1]] <- list()

  for (g2 in coexp_winner_genes) {
    p <- make_coexp_plot(g1, g2, umap_coords_coexp, hex_data_coexp,
                         color_ll, color_hl, color_lh, color_hh)

    if (!is.null(p)) {
      figures$S5$A[[g1]][[g2]] <- p
      cat("✓ figures$S5$A$", g1, "$", g2, "\n", sep = "")
    } else {
      cat("⚠ Skipped (gene not found):", g1, "×", g2, "\n")
    }
  }
}

cat("\n✓ S5A complete:", sum(sapply(figures$S5$A, length)), "plots stored\n")

# === Coexpression UMAPs - S5B (Loser Genes) ===

figures$S5$B <- list()

for (g1 in coexp_gene1) {
  figures$S5$B[[g1]] <- list()

  for (g2 in coexp_loser_genes) {
    p <- make_coexp_plot(g1, g2, umap_coords_coexp, hex_data_coexp,
                         color_ll, color_hl, color_lh, color_hh)

    if (!is.null(p)) {
      figures$S5$B[[g1]][[g2]] <- p
      cat("✓ figures$S5$B$", g1, "$", g2, "\n", sep = "")
    } else {
      cat("⚠ Skipped (gene not found):", g1, "×", g2, "\n")
    }
  }
}

cat("\n✓ S5B complete:", sum(sapply(figures$S5$B, length)), "plots stored\n")

# ============================================================================
# PLOT SAVING (COMMENTED OUT)
# ============================================================================
#
# === SAVE ALL FIGURES TO DISK ===
# Saves all plots from figures object to ../output/plots/ (600 DPI, Google Slides size)
# Uncomment the code below to activate plot saving
#
# cat("\n--- Saving Publication Plots ---\n")
# cat("Format: PNG @ 600 DPI, Google Slides size (10\" × 5.625\")\n\n")
#
# # Create output directory if it doesn't exist
# plots_output_dir <- file.path(output_dir, "plots")
# if (!dir.exists(plots_output_dir)) {
#   dir.create(plots_output_dir, recursive = TRUE)
# }
#
# # Google Slides 16:9 dimensions (inches)
# gs_width  <- 10
# gs_height <- 5.625
# gs_dpi    <- 600
#
# # Helper: save a single plot to PNG
# save_figure <- function(plot_obj, figure_name, output_dir, width, height, dpi) {
#   png_path <- file.path(output_dir, paste0(figure_name, ".png"))
#
#   if (!is.null(plot_obj)) {
#     # Save PNG
#     ggsave(png_path, plot_obj, width = width, height = height, dpi = dpi, bg = "white")
#     cat("✓ Saved", figure_name, ".png\n")
#   }
# }
#
# # Iterate through figures object and save all plots
# # S4a, S4e: single plots
# save_figure(figures$S4$a, "Figure_S4_a_hex_clusters_umap", plots_output_dir, gs_width, gs_height, gs_dpi)
# save_figure(figures$S4$e, "Figure_S4_e_enrichment_heatmap", plots_output_dir, gs_width, gs_height, gs_dpi)
#
# # S4b, S4c, S4d: nested lists by gene
# for (gene in names(figures$S4$b)) {
#   save_figure(figures$S4$b[[gene]], paste0("Figure_S4_b_", gene, "_umap"), plots_output_dir, gs_width, gs_height, gs_dpi)
# }
# for (gene in names(figures$S4$c)) {
#   save_figure(figures$S4$c[[gene]], paste0("Figure_S4_c_", gene, "_umap"), plots_output_dir, gs_width, gs_height, gs_dpi)
# }
# for (gene in names(figures$S4$d)) {
#   save_figure(figures$S4$d[[gene]], paste0("Figure_S4_d_", gene, "_umap"), plots_output_dir, gs_width, gs_height, gs_dpi)
# }
#
# # S5A, S5B: nested lists by gene1 and gene2
# for (gene1 in names(figures$S5$A)) {
#   for (gene2 in names(figures$S5$A[[gene1]])) {
#     save_figure(
#       figures$S5$A[[gene1]][[gene2]],
#       paste0("Figure_S5_A_", gene1, "_", gene2, "_coexpression"),
#       plots_output_dir, gs_width, gs_height, gs_dpi
#     )
#   }
# }
# for (gene1 in names(figures$S5$B)) {
#   for (gene2 in names(figures$S5$B[[gene1]])) {
#     save_figure(
#       figures$S5$B[[gene1]][[gene2]],
#       paste0("Figure_S5_B_", gene1, "_", gene2, "_coexpression"),
#       plots_output_dir, gs_width, gs_height, gs_dpi
#     )
#   }
# }
#
# cat("\n✓ All figures saved to:", plots_output_dir, "\n")

# === SUMMARY ===
cat("\n", strrep("=", 70), "\n", sep = "")
cat("ANALYSIS PIPELINE COMPLETE\n")
cat(strrep("=", 70), "\n\n", sep = "")

cat("Output Summary:\n")
cat("  ✓ Data files saved to: ../output/data/\n")
cat("    - degs_results.rds\n")
cat("    - enrichment_results.rds\n\n")

cat("  ℹ Plot saving is currently disabled (commented out)\n")
cat("    - Uncomment the PLOT SAVING section above to generate PNG files\n\n")

cat("Figures object structure:\n")
cat("  $S4$a  : Hex Data UMAP\n")
cat("  $S4$b  : Gene expression UMAPs (NANOG, MYC, POU5F1)\n")
cat("  $S4$c  : Gene expression UMAPs (ISL1, TFAP2A)\n")
cat("  $S4$d  : Gene expression UMAPs (MIXL1, TBXT)\n")
cat("  $S4$e  : Enrichment heatmap\n")
cat("  $S5$A  : Coexpression UMAPs (Winners) -", sum(sapply(figures$S5$A, length)), "plots\n")
cat("  $S5$B  : Coexpression UMAPs (Losers) -", sum(sapply(figures$S5$B, length)), "plots\n")
cat("Code running is complete. Live long and prosper. 🖖\n\n")
