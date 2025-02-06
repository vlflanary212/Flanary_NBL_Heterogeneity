# Author: Victoria Flanary
# Date: 250206
# Objective: Add cluster labels to annotation comparison UMAP.

# Set.seed
set.seed(42)

# Load packages
library(Seurat)
library(tidyverse)
library(patchwork)
library(here)

# Set filepaths
data_dir <- here("NBAtlas", "data", "alldata")
results_dir <- here("NBAtlas", "04_cluster_annotation", "03_annotation")

# Load data
seurat_obj <- readRDS(here(data_dir, "08_final_anno.rds"))

# Set current idents to Cell_type
seurat_obj$Cell_type <- seurat_obj@active.ident

# Set idents to original annotation
seurat_obj <- SetIdent(seurat_obj, value = "orig.ident")

# Remove cells originally annotated as RBCs from the object
cell_types <- as.character(setdiff(seurat_obj$orig.ident, "RBCs"))
seurat_obj <- subset(seurat_obj, idents = cell_types)

# Compare known annotations to the original NBAtlas annotations
pdf(
  here(results_dir, "cluster_annotation_comparison_labeled.pdf"),
  height = 4, width = 13 
)

p1 <- DimPlot(seurat_obj, reduction = "umap_harmony", raster = TRUE, 
              shuffle = TRUE, group.by = "orig.ident", label = TRUE) +
  labs(title = "Original Annotations", x = "UMAP_Harmony_1", y = "UMAP_Harmony_2") +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(seurat_obj, reduction = "umap_harmony", raster = TRUE, 
              shuffle = TRUE, group.by = "Cell_type", label = TRUE) +
  labs(title = "New Annotations", x = "UMAP_Harmony_1", y = "UMAP_Harmony_2") +
  theme(plot.title = element_text(hjust = 0.5))

print(p1 + p2)

dev.off()

# Remove the legend, since directly labelling clusters
pdf(
  here(results_dir, "cluster_annotation_comparison_labeled_no_legend.pdf"),
  height = 4, width = 9
)

p1 <- DimPlot(seurat_obj, reduction = "umap_harmony", raster = TRUE,
              shuffle = TRUE, group.by = "orig.ident", label = TRUE) +
  labs(title = "Original Annotations", x = "UMAP_Harmony_1", y = "UMAP_Harmony_2") +
  theme(plot.title = element_text(hjust = 0.5)) + NoLegend()

p2 <- DimPlot(seurat_obj, reduction = "umap_harmony", raster = TRUE,
              shuffle = TRUE, group.by = "Cell_type", label = TRUE) +
  labs(title = "New Annotations", x = "UMAP_Harmony_1", y = "UMAP_Harmony_2") +
  theme(plot.title = element_text(hjust = 0.5)) + NoLegend()

print(p1 + p2)

dev.off()


# End of script
