# Plot relevant UMAP figures

# Set seed
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
data_dir <- here("NB_ITH", "NBAtlas", "data", "alldata")
results_dir <- here("NB_ITH", "NBAtlas", "misc_figures")

# Load data
seurat_unint <- readRDS(here(data_dir, "03_unint_umap_obj.rds"))
seurat_int <- readRDS(here(data_dir, "05_harmony_clust.rds"))

# Plot the original UMAP and cell type annotations
## Representative image
pdf(
  here(results_dir, "umap_init.pdf"),
  height = 4, width = 6
)

print(
  DimPlot(seurat_int, reduction = "scvi_umap", group.by = "Cell_type",
          raster = TRUE, shuffle = TRUE) + 
    labs(
      title = "NBAtlas Cell Types",
      x = "UMAP_scVI_1", y = "UMAP_scVI_2"
    ) +
    theme(plot.title = element_text(hjust = 0.5))
)

dev.off()

## For integration comparison
pdf(
  here(results_dir, "umap_scvi.pdf"),
  height = 4, width = 6
)

print(
  DimPlot(seurat_int, reduction = "scvi_umap", group.by = "Cell_type",
          raster = TRUE, shuffle = TRUE) + 
    labs(
      title = "scVI Integration",
      x = "UMAP_scVI_1", y = "UMAP_scVI_2"
    ) +
    theme(plot.title = element_text(hjust = 0.5))
)

dev.off()

# Plot the same annotations on the Harmony-integrated embedding
pdf(
  here(results_dir, "umap_harmony.pdf"),
  height = 4, width = 6
)

print(
  DimPlot(seurat_int, reduction = "umap_harmony", group.by = "Cell_type",
          raster = TRUE, shuffle = TRUE) + 
    labs(
      title = "Harmony Integration",
      x = "UMAP_Harmony_1", y = "UMAP_Harmony_2"
    ) +
    theme(plot.title = element_text(hjust = 0.5))
)

dev.off()

# Compare integrated vs unintegrated umaps by batches
## unintegrated
pdf(
  here(results_dir, "umap_batch_unint.pdf"),
  height = 4.5, width = 9
)

print(
  DimPlot(seurat_unint, reduction = "umap", group.by = "Batch",
          raster = TRUE, shuffle = TRUE) + 
    labs(
      title = "Unintegrated UMAP",
      x = "UMAP_1", y = "UMAP_2"
    ) +
    theme(plot.title = element_text(hjust = 0.5))
)

dev.off()

## harmony-integrated
pdf(
  here(results_dir, "umap_batch_int_harmony.pdf"),
  height = 4.5, width = 9
)

print(
  DimPlot(seurat_int, reduction = "umap_harmony", group.by = "Batch",
          raster = TRUE, shuffle = TRUE) + 
    labs(
      title = "Harmony-Integrated UMAP",
      x = "UMAP_Harmony_1", y = "UMAP_Harmony_2"
    ) +
    theme(plot.title = element_text(hjust = 0.5))
)

dev.off()

## scvi-integrated
pdf(
  here(results_dir, "umap_batch_int_scvi.pdf"),
  height = 4.5, width = 9
)

print(
  DimPlot(seurat_int, reduction = "scvi_umap", group.by = "Batch",
          raster = TRUE, shuffle = TRUE) + 
    labs(
      title = "scVI-Integrated UMAP",
      x = "UMAP_scVI_1", y = "UMAP_scVI_2"
    ) +
    theme(plot.title = element_text(hjust = 0.5))
)

dev.off()

# End of script
