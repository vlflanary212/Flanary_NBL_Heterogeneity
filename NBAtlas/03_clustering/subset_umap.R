# Author: Victoria Flanary
# Date: 250205
# Objective: Plot clusters with similar colors on UMAP to determine 
# whether separate clusters are being grouped together.

# Set.seed
set.seed(42)

# Load packages
library(Seurat)
library(harmony)
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

# Set the correct resolution ident
seurat_obj <- SetIdent(seurat_obj, value = "RNA_snn_res.0.2")

# Subset the data for select clusters
subset <- subset(seurat_obj, idents = 18:19)

# Plot subsetted data on UMAP
pdf(
  here(results_dir, "subset_umap.pdf"),
  height = 4, width = 6
)

p <- DimPlot(subset, reduction = "umap_harmony", raster = TRUE, shuffle = TRUE) +
    labs(
      title = "Unclear Clusters",
      x = "UMAP_Harmony_1", y = "UMAP_Harmony_2"
    ) +
    theme(plot.title = element_text(hjust = 0.5))
  
  print(p)

dev.off()
