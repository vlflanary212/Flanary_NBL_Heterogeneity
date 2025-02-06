# Author: Victoria Flanary
# Date: 250206
# Objective: Run DEA to identify cluster-specific marker genes for annotated clusters.
# Re-running to get correct cluster labels in DEA results and save both
# significant and non-significant results for gene ontology analysis.

# Set.seed
set.seed(42)

# Load packages
library(Seurat)
library(tidyverse)
library(here)
library(future)

# Future parameters
options(future.globals.maxSize = 256 * 1024^3)  # 256 GB maximum RAM
plan("multisession", workers = 4)  # parallelize across 4 cores

# Set filepaths
data_dir <- here("NBAtlas", "data", "alldata")
results_dir <- here("NBAtlas", "04_cluster_annotation", "04_annotated_dea")

# Load data
seurat_obj <- readRDS(here(data_dir, "08_final_anno.rds"))

# Check cluster sizes
table(seurat_obj@active.ident)

# Clusters drastically range in size from several thousand to a few 100
# Want to balance groups to reduce bias toward larger clusters while 
# Preserving biological signals from these large clusters

# Subsample clusters so that no cluster has more than 1000 cells
seurat_subset <- subset(
  seurat_obj,
  cells = WhichCells(seurat_obj, downsample = 1000)
)
table(seurat_subset@active.ident)

# DEA
## Wilcoxon (default)
marker_genes <- FindAllMarkers(
  seurat_subset,
  test.use = "wilcox",
  min.pct = 0.1,
  min.diff.pct = 0.1,
  only.pos = TRUE,
  assay = "RNA"
)

## Save all results
write_csv(
marker_genes,
here(results_dir, "all_marker_genes.csv")
)

## Filter and arrange DEGs
marker_genes_filt <- marker_genes |>
  filter(p_val_adj <= 0.05) |>
  filter(avg_log2FC >= 2) |>
  arrange(cluster, desc(avg_log2FC)) |>
  distinct()

## Save the DEGs
write_csv(
  marker_genes_filt,
  here(results_dir, "filt_marker_genes.csv")
)

## Plot the top 25 cluster-specific DEGs
top25 <- marker_genes_filt |>
  group_by(cluster) |>
  arrange(p_val_adj) |>
  slice_head(n = 25)

pdf(
  here(results_dir, "marker_genes_bar_plots_labeled.pdf"),
  height = 7, width = 9
)

par(mfrow = c(2, 4), mar = c(4, 6, 3, 1))
for (i in unique(top25$cluster)) {
  barplot(
    sort(setNames(top25$avg_log2FC, top25$gene)[top25$cluster == i], FALSE),
    horiz = TRUE,
    las = 1,
    main = paste0(i, " vs. rest"),
    border = "white",
    yaxs = "i"
  )
  abline(v = c(0, 0.25), lty = c(1, 2))
}

dev.off()

# End of script
