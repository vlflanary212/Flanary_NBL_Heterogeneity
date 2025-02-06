# Author: Victoria Flanary
# Date: 250129
# Objective: Plot feature plots for select genes.

# Set.seed
set.seed(42)

# Load packages
library(Seurat)
library(tidyverse)
library(patchwork)
library(here)
library(future)

# Future parameters
options(future.globals.maxSize = 256 * 1024^3)  # 256 GB maximum RAM
plan("multisession", workers = 4)  # parallelize across 4 cores

# Set filepaths
data_dir <- here("NBAtlas", "data", "alldata")
results_dir <- here("NBAtlas", "04_cluster_annotation", "02_feature_plots")

# Load data
seurat_obj <- readRDS(here(data_dir, "06_dea_obj.rds"))

# Plot marker genes on feature plots
## Side-by-side
pdf(
  here(results_dir, "mycn_twist1_feature_plots.pdf"),
  height = 4, width = 10
)

p <- FeaturePlot(seurat_obj, reduction = "umap_harmony", raster = TRUE, ncol = 2,
              features = c("MYCN", "TWIST1"))
p <- p & labs(x = "UMAP_Harmony_1", y = "UMAP_Harmony_2")  # edit axis labels
print(p)

dev.off()

## Same plot
pdf(
  here(results_dir, "mycn_twist1_feature_plot_combined.pdf"),
  height = 4, width = 16
)

print(
  FeaturePlot(seurat_obj, reduction = "umap_harmony", raster = TRUE, blend = TRUE,
              features = c("MYCN", "TWIST1"))
)

dev.off()

# End of script
