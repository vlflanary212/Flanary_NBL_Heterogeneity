# Author: Victoria Flanary
# Date: 250129
# Objective: Name TME cell clusters and run gene ontology analysis on the
# neuroendocrine clusters.

# Set.seed
set.seed(42)

# Load packages
library(Seurat)
library(tidyverse)
library(patchwork)
library(here)

# Set filepaths
data_dir <- here("NB_ITH", "NBAtlas", "data", "alldata")
results_dir <- here("NB_ITH", "NBAtlas", "04_cluster_annotation")

# Load data
seurat_obj <- readRDS(here(data_dir, "06_dea_obj.rds"))

# Rename non-tumor cell clusters
seurat_anno <- RenameIdents(
  seurat_obj,
  '0' = "T cells",
  '12' = "NK cells",
  '10' = "B cells",
  '19' = "Plasma cells",
  '5' = "Myeloid cells",
  '11' = "Endothelial cells",
  '7' = "Fibroblasts",
  '20' = "Schwann Cell Precursors"
)

# Split the object by annotated and non-annotated clusters
known <- c("T cells", "NK cells", "B cells", "Plasma cells", "Myeloid cells",
           "Endothelial cells", "Fibroblasts", "Schwann Cell Precursors")

# Get current cluster identities
current_clusters <- levels(Idents(seurat_anno))

# Add "Neuroendocrine_" prefix to clusters not in `known`
new_cluster_names <- ifelse(current_clusters %in% known, 
                            current_clusters, 
                            paste0("Cluster_", current_clusters))

names(new_cluster_names) <- current_clusters
seurat_anno <- RenameIdents(seurat_anno, new_cluster_names)

# Compare known annotations to the original NBAtlas annotations
pdf(
  here(results_dir, "init_anno_comparison.pdf"),
  height = 4, width = 12 
)

p1 <- DimPlot(seurat_anno, reduction = "umap_harmony", raster = TRUE, 
              shuffle = TRUE, group.by = "orig.ident") +
  labs(title = "Original Annotations", x = "UMAP_1", y = "UMAP_2") +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(seurat_anno, reduction = "umap_harmony", raster = TRUE, 
              shuffle = TRUE) +
  labs(title = "New Annotations", x = "UMAP_1", y = "UMAP_2") +
  theme(plot.title = element_text(hjust = 0.5))

print(p1 + p2)

dev.off()

# Save the initially annotated object
saveRDS(
  seurat_anno,
  here(data_dir, "07_init_anno.rds")
)

# End of script