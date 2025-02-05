# Author: Victoria Flanary
# Date: 250129
# Objective: Name TME cell clusters and run gene ontology analysis on the
# neuroendocrine clusters.

# Set.seed
set.seed(42)

# Load packages
library(Seurat)
library(topGO)
library(org.Hs.eg.db)
library(tidyverse)
library(here)

# Set filepaths
data_dir <- here("NB_ITH", "NBAtlas", "data", "alldata")
results_dir <- here("NB_ITH", "NBAtlas", "04_cluster_annotation")

# Load data
seurat_obj <- readRDS(here(data_dir, "07_init_anno.rds"))

# Subset data for non-annotated clusters
seurat_obj <- subset(
  seurat_obj,
  idents = grep("^Cluster_", levels(Idents(seurat_obj)), value = TRUE)
  )

# Load DEA results
marker_df <- read_csv(
  here(results_dir, "marker_genes_wilcox.csv")
)

# Run TopGO on the unknown clusters
# Function to perform GO analysis for a given cluster
run_GO_analysis <- function(cluster, marker_df, bg_genes) {
  # Extract genes for the cluster
  genes <- marker_df %>%
    filter(cluster == !!cluster) %>%
    pull(gene) %>%
    unique()
  
  # Create gene list
  genelist <- factor(as.integer(bg_genes %in% genes))
  names(genelist) <- bg_genes
  
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

# Get background gene list (all genes in the dataset)
bg_genes <- rownames(GetAssayData(seurat_obj, layer = "counts"))

# Get unique clusters from marker results
clusters <- unique(marker_df$cluster)

# Run GO analysis for each cluster and combine results
go_results <- bind_rows(lapply(clusters, run_GO_analysis, marker_df = marker_df, bg_genes = bg_genes))

# Save results
write_csv(
  go_results,
  here(results_dir, "topgo_results.csv")
)

# Visualize results
# Convert p-values to -log10 for better visualization
go_results <- go_results |> 
  mutate(classicFisher = as.numeric(classicFisher)) |> # Convert properly
  filter(!is.na(classicFisher)) |>  # Remove rows where p-value conversion failed
  mutate(log_pval = -log10(classicFisher)) |> 
  arrange(Cluster, classicFisher)

# Keep only the **top 10** GO terms per cluster
top_go_per_cluster <- go_results |>
  group_by(Cluster) |>
  slice_min(order_by = classicFisher, n = 10) |>
  ungroup()

# Plot the results for each cluster
pdf(
  here(results_dir, "topgo_bar_plots.pdf"),
  height = 6, width = 8
)
for (i in unique(top_go_per_cluster$Cluster)) {
  # Filter top 10 GO terms for the current cluster
  top_go <- go_results |> 
    filter(Cluster == i) |> 
    slice_min(order_by = classicFisher, n = 10)
  
  # Generate plot
  p <- ggplot(top_go, aes(x = reorder(Term, -log_pval), y = log_pval, fill = Term)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_minimal(base_size = 12) +
    labs(title = paste0("Cluster ", i, " - Top 10 GO Terms"),
         x = "GO Biological Process",
         y = "-log10(p-value)") +
    theme(legend.position = "none")
  
  print(p)
}

dev.off()

# End of script
