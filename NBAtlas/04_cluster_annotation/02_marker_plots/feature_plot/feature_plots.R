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
results_dir <- here("NBAtlas", "04_cluster_annotation", "02_marker_plots", "feature_plot")

# Load data
seurat_obj <- readRDS(here(data_dir, "06_dea_obj.rds"))

# Plot marker genes on feature plots
pdf(
  here(results_dir, "marker_gene_feature_plots.pdf"),
  height = 40, width = 15
)

print(
  FeaturePlot(seurat_obj, reduction = "umap_harmony", raster = TRUE, ncol = 3,
              features = c("CD2", "CD3D", "LCK",    # T cells
                           "KLRF1", "XCL1", "GNLY",  # NK cells
                           "CD19", "MS4A1", "CD22", # B cells
                           "JCHAIN", "TNFRSF17", "MZB1", # Plasma cells
                           "CD163", "CXCL8", "LYZ", # Macrophages
                           "CD34", "KDR", "CDH5", # Endothelial cells
                           "COL1A1", "COL1A2", "DCN",  # Fibroblasts
                           "CDH19", "ERBB3", "SOX10",  # Schwann cell precursors
                           "PRPH", "ELAVL4", "STMN2", # Sympathoblasts
                           "CHGB", "DBH", "PNMT"  # Chrommafin cells
                           ))
)

dev.off()

# Plot ADR vs MES TFs on feature plots
## Extract MES and ADR TFs
nb_tfs <- read.table(
  here("docs", "nb_phenotype_tfs.txt"),
  header = TRUE
)

mes_tfs <- nb_tfs$gene[nb_tfs$group == "MES"]
adr_tfs <- nb_tfs$gene[nb_tfs$group == "ADRN"]

## ADR TFs
pdf(
  here(results_dir, "adr_tf_feature_plots.pdf"),
  height = 20, width = 20
)

print(
  FeaturePlot(seurat_obj, reduction = "umap_harmony", raster = TRUE, ncol = 4,
              features = adr_tfs)
)

dev.off()

## MES TFs
pdf(
  here(results_dir, "mes_tf_feature_plots.pdf"),
  height = 20, width = 20
)

print(
  FeaturePlot(seurat_obj, reduction = "umap_harmony", raster = TRUE, ncol = 4,
              features = mes_tfs)
)

dev.off()

## Other genes of interest
pdf(
  here(results_dir, "misc_feature_plots.pdf"),
  height = 4, width = 15
)

print(
  FeaturePlot(seurat_obj, reduction = "umap_harmony", raster = TRUE, ncol = 3,
              features = c("MYCN", "VIM", "FOXJ3"))
)

dev.off()

# Plot QC metrics on feature plot to ensure these don't define any cluster
features <- c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo")

pdf(
  here(results_dir, "qc_feature_plots.pdf"),
  height = 8, width = 10
)

print(
  FeaturePlot(seurat_obj, reduction = "umap_harmony", raster = TRUE,
              ncol = 2, features = features)
)

dev.off()

# End of script

