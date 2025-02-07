# Author: Victoria Flanary
# Date: 250206
# Objective: Plot dot plots of cluster-specific marker genes.

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
results_dir <- here("NBAtlas", "05_cluster_analyses", "02_dot_plots")

# Load data
seurat_obj <- readRDS(here(data_dir, "08_final_anno.rds"))
ne_obj <- readRDS(here(data_dir, "09_ne_subset.rds"))
tme_obj <- readRDS(here(data_dir, "09_tme_subset.rds"))

# Define markers 
tme_markers <- c("CD3D", "CD3E",  # T cells
                 "KLRF1", "XCL1", # NK cells
                 "MS4A1", "CD19", # B cells
                 "JCHAIN", "MZB1", # Plasma cells
                 "CLEC4C", "LILRA4", # pDCs
                 "CXCL8", "LYZ",  # Macrophages/Myeloid cells
                 "CD34", "VWF",   # Endothelial cells
                 "COL1A1", "COL1A2", # Fibroblasts
                 "CDH19", "PLP1",    # Schwann cell precursors
                 "STAR", "CYP11B1")  # Steroidogenic adrenocortical cells

ne_markers <- c("MYCN", "ALK",      # Commonly mutated genes
                "PHOX2B", "GATA3",  # Adrenergic markers
                "CHGB", "PRPH",     # Neuroendocrine markers
                "SNAI2", "ZEB1", "TWIST1",  # Mesenchymal/EMT markers
                "MKI67", "UBE2C")   # Proliferation markers

# Plot marker genes on dot plots
## All clusters
pdf(
  here(results_dir, "all_clusters_dot_plot.pdf"),
  height = 12, width = 16
)

print(
  DotPlot(seurat_obj, features = c(tme_markers, ne_markers)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)

dev.off()

## NE clusters
pdf(
  here(results_dir, "ne_clusters_dot_plot.pdf"),
  height = 6, width = 8
)

print(
  DotPlot(ne_obj, features = ne_markers) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)

dev.off()

## TME clusters
pdf(
  here(results_dir, "tme_clusters_dot_plot.pdf"),
  height = 6, width = 8
)

print(
  DotPlot(tme_obj, features = tme_markers) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)

dev.off()

# End of script
