# Author: Victoria Flanary
# Date: 250208
# Objective: Run gene ontology analysis on neuroendocrine clusters.

# Set seed for reproducibility
set.seed(42)

# Load required packages
library(topGO)
library(org.Hs.eg.db)
library(tidyverse)
library(here)
library(ggplot2)
library(gridExtra)

# Set filepaths
data_dir <- here("NBAtlas", "data", "alldata")
results_dir <- here("NBAtlas", "05_cluster_analyses")

# Load DEA results (ensure this includes all tested genes, not just significant ones)
dea_results <- read_csv(
  here(results_dir, "01_annotated_dea", "all_marker_genes.csv")
)

# Extract the background genes (all genes used in the DEA)
background <- unique(dea_results$gene)

# Subset clusters to only include those starting with "NE_"
ne_clusters <- paste0("NE_", 1:15)
dea_results <- dea_results |> filter(cluster %in% ne_clusters)

# Function to perform GO analysis for a given cluster
run_GO_analysis <- function(cluster, dea_results, background) {
  # Extract **only significant** genes for the current cluster
  significant_genes <- dea_results |> 
    filter(cluster == !!cluster, p_val_adj <= 0.05) |>  # Adjust significance threshold as needed
    pull(gene) |> 
    unique()
  
  # Skip GO analysis if there are no significant genes for the cluster
  if (length(significant_genes) == 0) {
    return(NULL)
  }

  # Create gene list (foreground = significant genes, background = all genes tested)
  genelist <- factor(as.integer(background %in% significant_genes))
  names(genelist) <- background
  
  # Initialize GOdata
  GOdata <- new("topGOdata",
                ontology = "BP",
                allGenes = genelist,
                annot = annFUN.org,
                mapping = "org.Hs.eg.db",
                ID = "symbol",
                nodeSize = 10)
  
  # Run classic Fisher exact test
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  
  # Extract top GO terms
  gobp <- GenTable(GOdata, classicFisher = resultFisher, topNodes = 15)
  
  # Add cluster information
  gobp$Cluster <- cluster
  
  return(gobp)
}

# Run GO analysis for each cluster and combine results
go_results <- bind_rows(lapply(ne_clusters, run_GO_analysis, dea_results = dea_results, background = background))

# Remove clusters with no significant terms
go_results <- go_results |> filter(!is.na(classicFisher))

# Save results
write_csv(
  go_results,
  here(results_dir, "04_gene_ontology", "topgo_results.csv")
)

# Convert p-values to -log10 for better visualization and handle conversion warnings
go_results <- go_results |> 
  mutate(classicFisher = suppressWarnings(as.numeric(classicFisher))) |> # Suppress conversion warnings
  filter(!is.na(classicFisher)) |>  # Remove rows where conversion failed
  filter(classicFisher <= 0.05) |> # Only keep significant results
  mutate(log_pval = -log10(classicFisher)) |> 
  arrange(Cluster, desc(log_pval)) |>
  mutate(Cluster = factor(Cluster, levels = paste0("NE_", 1:15)))

# Keep only the **top 5** GO terms per cluster
top_go_per_cluster <- go_results |> 
  group_by(Cluster) |> 
  slice_max(order_by = log_pval, n = 5) |> 
  ungroup()

# Define colors dynamically based on number of clusters
num_clusters <- length(unique(top_go_per_cluster$Cluster))
if (num_clusters > 12) {
  cluster_colors <- scales::hue_pal()(num_clusters)  # Generate distinct colors if >12 clusters
} else {
  cluster_colors <- RColorBrewer::brewer.pal(num_clusters, "Set3")
}
names(cluster_colors) <- unique(top_go_per_cluster$Cluster)

# Create a multi-panel PDF plot with 5 columns and 3 rows
pdf(
  here(results_dir, "04_gene_ontology", "topgo_bar_plots.pdf"),
  height = 12, width = 30
)

# Generate plots for each cluster
plots <- list()
for (i in unique(top_go_per_cluster$Cluster)) {
  top_go <- top_go_per_cluster |> filter(Cluster == i)
  
  p <- ggplot(top_go, aes(x = reorder(Term, -log_pval), y = log_pval, fill = Cluster)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_minimal(base_size = 12) +
    labs(title = paste0(i, " - Top 5 GO Terms"),
         x = "GO Biological Process",
         y = "-log10(p-value)") +
    theme(legend.position = "none") +
    scale_fill_manual(values = cluster_colors)  # Use dynamic color assignment
  
  plots[[i]] <- p
}

# Arrange the plots in a 5-column, 3-row layout
gridExtra::grid.arrange(grobs = plots, ncol = 5, nrow = 3)

dev.off()

# End of Script

