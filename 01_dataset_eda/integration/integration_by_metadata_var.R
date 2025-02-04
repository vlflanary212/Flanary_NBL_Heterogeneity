# Author: Victoria Flanary
# Date: 240122
# Objective: Visualize the number of tumors belonging to each metadata
# variable in the NBAtlas

# Set-up
## Load  packages
library(tidyverse)
library(here)

## Load the seurat_obj
seurat_obj <- readRDS(
  here("NB_ITH", "NBAtlas", "data", "alldata", "00_init_obj.rds")
)

## Load metadata
metadata <- read_csv(
  here("NB_ITH", "NBAtlas", "metadata_filt.csv")
) |> as.data.frame()

## Set results directory
results_dir <- here(
  "NB_ITH", "NBAtlas", "01_dataset_eda", "integration"
)

# Visualize degree of integration by different metadata vars on UMAP
## Plot in a single pdf
pdf(
    here(results_dir, "Integration_by_Metadata.pdf"),
    height = 4, width = 6
  )

for (i in colnames(metadata)) {
   print(
    DimPlot(seurat_obj, reduction = "scvi_umap", group.by = i,
            shuffle = TRUE, raster = TRUE) +
      labs(
        title = paste0("NBAtlas UMAP by ", i),
        x = "scVI_UMAP_1",
        y = "scVI_UMAP_2"
      ) +
      theme(plot.title = element_text(hjust = 0.5))
  )
}

dev.off()

## Plot as individual files
for (i in colnames(metadata)) {
  pdf(
    here(results_dir, paste0("UMAP_by_", i, ".pdf")),
    height = 4, width = 6
  )
  
  print(
    DimPlot(seurat_obj, reduction = "scvi_umap", group.by = i,
            shuffle = TRUE, raster = TRUE) +
      labs(
        title = paste0("NBAtlas UMAP by ", i),
        x = "scVI_UMAP_1",
        y = "scVI_UMAP_2"
      ) +
      theme(plot.title = element_text(hjust = 0.5))
  )
  
  dev.off()
}

# End of script