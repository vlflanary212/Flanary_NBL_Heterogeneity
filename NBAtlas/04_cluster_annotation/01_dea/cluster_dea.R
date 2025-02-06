# Author: Victoria Flanary
# Date: 250129
# Objective: Run DEA to identify cluster-specific marker genes.

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
results_dir <- here("NBAtlas", "04_cluster_annotation", "01_dea")

# Load data
seurat_obj <- readRDS(here(data_dir, "05_harmony_clust.rds"))

# JoinLayers for DEA
seurat_obj <- JoinLayers(seurat_obj)

# Set Ident to desired cluster size
# 0.3 is max resolution before cells look overclustered on UMAP
# Also the lowest cluster that distinguishes between pDCs and plasma cells
seurat_obj <- SetIdent(seurat_obj, value = "RNA_snn_res.0.3")
table(seurat_obj@active.ident)

# Save current object
saveRDS(
  seurat_obj,
  here(data_dir, "06_dea_obj.rds")
)

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
marker_genes_wilcox <- FindAllMarkers(
  seurat_subset,
  test.use = "wilcox",
  min.pct = 0.1,
  min.diff.pct = 0.1,
  only.pos = TRUE,
  assay = "RNA"
)

## Filter and arrange DEGs
wilcox_filt <- marker_genes_wilcox |>
  filter(p_val_adj < 0.05) |>
  filter(avg_log2FC > 2) |>
  arrange(cluster, desc(avg_log2FC)) |>
  distinct()

## Save the DEGs
write_csv(
  wilcox_filt,
  here(results_dir, "marker_genes_wilcox.csv")
)

## Plot the top 25 cluster-specific DEGs
top25 <- wilcox_filt |>
  group_by(cluster) |>
  arrange(p_val_adj) |>
  slice_head(n = 25)

pdf(
  here(results_dir, "marker_genes_wilcox_bar_plots.pdf"),
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
