# Author: Victoria Flanary
# Date: 250129
# Objective: Cluster cells in the NBAtlas for cell type annotation.

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
data_dir <- here("NB_ITH", "NBAtlas", "data", "alldata")
results_dir <- here("NB_ITH", "NBAtlas", "03_clustering")

# Load data
seurat_obj <- readRDS(here(data_dir, "01_filt_obj.rds"))

## Split the Seurat object into layers by batch
seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj$Batch)

# Normalize the data
seurat_obj <- NormalizeData(seurat_obj)

# Find highly variable features
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = TRUE
)

# Plot highly variable features
hvf_info <- seurat_obj |>
  HVFInfo(method = "vst") |>
  arrange(desc(variance.standardized))

top20 <- hvf_info |> 
  slice_head(n = 20) |>
  rownames()

pdf(
  here(results_dir, "top_hvfs.pdf"),
  height = 4, width = 6
)

hvf_plot <- VariableFeaturePlot(seurat_obj)
print(
  LabelPoints(hvf_plot, points = top20, repel = TRUE)
)

dev.off()

# Scale the data
seurat_obj <- ScaleData(seurat_obj)

# Run PCA
seurat_obj <- RunPCA(seurat_obj, npcs = 50, approx = FALSE)

# Plot PCs on elbow plot
pdf(
  here(results_dir, "elbow_plot.pdf"),
  height = 4, width = 6
)

print(ElbowPlot(seurat_obj, reduction = "pca", ndims = 50))

dev.off()

# Determine the number of calculated PCs that capture at least 95% of the variance
## Extract the standard dev for PCs used to run RunPCA
pca_stdev <- seurat_obj@reductions[["pca"]]@stdev[1:50]

## Calculate the percent variance explained by each calculated PC
variance_explained <- (pca_stdev^2) / sum(pca_stdev^2) * 100

## Calculate the cumulative percent variance explained by each calculated PC
cumulative_variance <- cumsum(variance_explained)

## Identify the minimum number of PCs that capture 95% of the variance
npcs <- which(cumulative_variance >= 95)[1]  # 40 PCs

# Save the object
saveRDS(
  seurat_obj,
  here(data_dir, "02_pca_obj.rds")
)

# Compute the KNN graph
seurat_clust <- FindNeighbors(seurat_obj, dims = 1:npcs, verbose = TRUE)

# Cluster cells at multiple resolutions
resolutions <- seq(from = 0.1, to = 0.9, by = 0.1)

for (i in resolutions) {
  print(paste0("Clustering at Louvain resolution ", i))
  seurat_clust <- FindClusters(seurat_clust, resolution = i)
  print(paste0("Finished clustering at Louvain resolution ", i))
}

# Run UMAP on the clustered object
seurat_clust <- RunUMAP(
  seurat_clust, reduction = "pca", dims = 1:npcs, reduction_name = "umap",
  min.dist = 0.5, spread = 1.5
)

# Plot unintegrated UMAPs by different clustering resolutions
res_names <- paste0("RNA_snn_res.", resolutions)
plot_list <- list()

pdf(
  here(results_dir, "unintegrated_umaps_by_resolution.pdf"),
  height = 4, width = 6
)

for (i in res_names){
  p <- DimPlot(seurat_clust, reduction = "umap", group.by = i, 
               raster = TRUE, shuffle = TRUE) +
    labs(
      x = "UMAP_1", y = "UMAP_2"
    ) +
    theme(plot.title = element_text(hjust = 0.5))
  
  plot_list[[i]] <- p
  print(plot_list[[i]])
}

dev.off()

# Plot unintegrated UMAPs by metadata
groups <- c("Cell_type", "Batch", "Study", "Assay", "Platform", 
            "Cell_condition", "Patient_No", "Gender", "MYCN_amplification", 
            "Risk", "INSS_stage", "Timepoint")
plot_list <- list()

pdf(
  here(results_dir, "unintegrated_umaps_by_metadata.pdf"),
  height = 5, width = 7.5
)

for (i in groups){
  p <- DimPlot(seurat_clust, reduction = "umap", group.by = i, 
               raster = TRUE, shuffle = TRUE) +
    labs(
      title = paste0("Unintegrated UMAP by ", i),
      x = "UMAP_1", y = "UMAP_2"
    ) +
    theme(plot.title = element_text(hjust = 0.5))
  
  plot_list[[i]] <- p
  print(plot_list[[i]])
}

dev.off()

saveRDS(
  seurat_clust,
  here(data_dir, "03_unint_umap_obj.rds")
)

# Integration
seurat_harmony <- IntegrateLayers(
  seurat_obj, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  dims = 1:npcs, verbose = TRUE
)

saveRDS(
  seurat_harmony,
  here(data_dir, "04_harmony_obj.rds")
)

# Compute the KNN graph
seurat_harmony_clust <- FindNeighbors(
  seurat_harmony, reduction = "harmony", dims = 1:npcs, verbose = TRUE
)

# Cluster cells at multiple resolutions
resolutions <- seq(from = 0.1, to = 0.9, by = 0.1)

for (i in resolutions) {
  print(paste0("Clustering at Louvain resolution ", i))
  seurat_harmony_clust <- FindClusters(seurat_harmony_clust, resolution = i)
  print(paste0("Finished clustering at Louvain resolution ", i))
}

# Run UMAP on the clustered object
seurat_harmony_clust <- RunUMAP(
  seurat_harmony_clust, reduction = "harmony", dims = 1:npcs, 
  reduction.name = "umap_harmony", min.dist = 0.5, spread = 1.5
)

# Plot integrated UMAPs by different clustering resolutions
res_names <- paste0("RNA_snn_res.", resolutions)
plot_list <- list()

pdf(
  here(results_dir, "integrated_umaps_by_resolution.pdf"),
  height = 4, width = 6
)

for (i in res_names){
  p <- DimPlot(seurat_harmony_clust, reduction = "umap_harmony", group.by = i, 
               raster = TRUE, shuffle = TRUE) +
    labs(
      x = "UMAP_Harmony_1", y = "UMAP_Harmony_2"
    ) +
    theme(plot.title = element_text(hjust = 0.5))
  
  plot_list[[i]] <- p
  print(plot_list[[i]])
}

dev.off()

# Plot integrated UMAPs by metadata
groups <- c("Cell_type", "Batch", "Study", "Assay", "Platform", 
            "Cell_condition", "Patient_No", "Gender", "MYCN_amplification", 
            "Risk", "INSS_stage", "Timepoint")
plot_list <- list()

pdf(
  here(results_dir, "integrated_umaps_by_metadata.pdf"),
  height = 5, width = 7.5
)

for (i in groups){
  p <- DimPlot(seurat_harmony_clust, reduction = "umap_harmony", group.by = i, 
               raster = TRUE, shuffle = TRUE) +
    labs(
      title = paste0("Integrated UMAP by ", i),
      x = "UMAP_Harmony_1", y = "UMAP_Harmony_2"
    ) +
    theme(plot.title = element_text(hjust = 0.5))
  
  plot_list[[i]] <- p
  print(plot_list[[i]])
}

dev.off()

# Save the integrated object
saveRDS(
  seurat_harmony_clust.rds,
  here(data_dir, "05_harmony_clust.rds")
)

# End of script  