# Author: Victoria Flanary
# Date: 240122
# Objective: Visualize the number of tumors belonging to each metadata
# variable in the NBAtlas

# Set-up
## Load  packages
library(Seurat)
library(tidyverse)
library(here)

## Load the seurat_obj
seurat_obj <- readRDS(here("NBAtlas", "data", "alldata", "00_formatted_obj.rds"))

## Load metadata
metadata <- read_csv(
  here("NBAtlas", "01_dataset_eda", "metadata_filt.csv")
) |> as.data.frame()

## Set results directory
results_dir <- here("NBAtlas", "01_dataset_eda", "umaps")

# Visualize degree of integration by Batch
pdf(here(results_dir, "UMAP_by_Batch.pdf"), 
height = 4, width = 9)
  
  print(
    DimPlot(seurat_obj, reduction = "scvi_umap", group.by = "Batch",
            shuffle = TRUE, raster = TRUE) +
      labs(
        title = "NBAtlas UMAP by Batch",
        x = "scVI_UMAP_1",
        y = "scVI_UMAP_2"
      ) +
      theme(plot.title = element_text(hjust = 0.5))
  )
  
  dev.off()

# End of script
