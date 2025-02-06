# Author: Victoria Flanary
# Date: 250129
# Objective: Name clusters with clear markers by DEA results and
# Compare with original NBAtlas annotations.

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
seurat_obj <- readRDS(here(data_dir, "06_dea_obj.rds"))

# Rename non-tumor cell clusters
seurat_anno <- RenameIdents(
  seurat_obj,
  '1' = "T cells",
  '23' = "T cells",
  '13' = "NK cells",
  '10' = "B cells",
  '22' = "Plasma cells",
  '25' = "pDCs",
  '5' = "Myeloid cells",
  '12' = "Endothelial cells",
  '6' = "Fibroblasts",
  '19' = "Schwann",
  '20' = "Adrenocortical cells"
)

# Save the initial annotations
saveRDS(
seurat_anno,
here(data_dir, "07_init_anno.rds")
)

# Create vector of known cell types
known <- c("T cells", "NK cells", "B cells", "Plasma cells", "pDCs", "Myeloid cells",
           "Endothelial cells", "Fibroblasts", "Schwann", "Adrenocortical cells")

# Get current cluster identities
current_clusters <- levels(Idents(seurat_anno))

# Identify numeric clusters
numeric_clusters <- setdiff(current_clusters, known)

# Create a mapping for numeric clusters to new sequential numbers
renumbered_clusters <- setNames(seq_along(numeric_clusters), numeric_clusters)

# Generate new cluster names
new_cluster_names <- ifelse(current_clusters %in% known,
                            current_clusters,
                            paste0("NE_", renumbered_clusters[current_clusters]))

# Rename the cluster identities in the Seurat object
Idents(seurat_anno) <- factor(Idents(seurat_anno), levels = current_clusters, labels = new_cluster_names)

# Compare known annotations to the original NBAtlas annotations
pdf(
  here(results_dir, "cluster_annotation_comparison.pdf"),
  height = 4, width = 13 
)

p1 <- DimPlot(seurat_anno, reduction = "umap_harmony", raster = TRUE, 
              shuffle = TRUE, group.by = "orig.ident") +
  labs(title = "Original Annotations", x = "UMAP_Harmony_1", y = "UMAP_Harmony_2") +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- DimPlot(seurat_anno, reduction = "umap_harmony", raster = TRUE, 
              shuffle = TRUE) +
  labs(title = "New Annotations", x = "UMAP_Harmony_1", y = "UMAP_Harmony_2") +
  theme(plot.title = element_text(hjust = 0.5))

print(p1 + p2)

dev.off()

# Split the object by neuroendocrine (ne) vs tumor microenvironment (tme) clusters
tme_obj <- subset(seurat_anno, idents = known)
ne_obj <- subset(seurat_anno, idents = setdiff(new_cluster_names, known))

# Plot TME and NE cell types seperately
## TME cell types
pdf(
  here(results_dir, "tumor_microenvironment_umap.pdf"),
  height = 5, width = 7.5
)

p <- DimPlot(tme_obj, reduction = "umap_harmony", raster = TRUE, shuffle = TRUE) +
  labs(title = "NBAtlas Tumor Microenvironment", x = "UMAP_Harmony_1", y = "UMAP_Harmony_2") +
  theme(plot.title = element_text(hjust = 0.5))

print(p)

dev.off()

## Tumor cell clusters
pdf(
  here(results_dir, "tumor_subclusters_umap.pdf"),
  height = 5, width = 7.5
)

p <- DimPlot(ne_obj, reduction = "umap_harmony", raster = TRUE, shuffle = TRUE) +
  labs(title = "NBAtlas Tumor Populations", x = "UMAP_Harmony_1", y = "UMAP_Harmony_2") +
  theme(plot.title = element_text(hjust = 0.5))

print(p)

dev.off()


# Save the initially annotated object
saveRDS(
  seurat_anno,
  here(data_dir, "08_final_anno.rds")
)

# Save subsetted objects
saveRDS(
  ne_obj,
  here(data_dir, "09_ne_subset.rds")
)

saveRDS(
  tme_obj,
  here(data_dir, "09_tme_subset.rds")
)

# End of script
