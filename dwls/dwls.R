# Author: Victoria Flanary
# Date: 250130
# Objective: Deconvolute PCGP RNA-seq data of patient neuroblastoma tumors
# using the DWLS algorithm and clustered NBAtlas as a reference.


## Set a seed
set.seed(42)

## Load packages
library(DWLS)
library(circlize)
library(ComplexHeatmap)
library(tidyverse)
library(here)
library(future)

# Future parameters
options(future.globals.maxSize = 256 * 1024^3)  # 256 GB maximum RAM
plan("multisession", workers = 4)

# Set results directory
results_dir <- here("NB_ITH", "PCGP_RNA", "dwls", "outputs")

# Load data as matrices 
signature_matrix <- read.table(
  here("NB_ITH", "PCGP_RNA", "dwls", "inputs", "single_cell_ref.txt")
) |> as.matrix()

bulk_data <- read.table(
  here("NB_ITH", "PCGP_RNA", "dwls", "inputs", "bulk_matrix.txt"),
) |> as.matrix()

## Ensure both objects are matrices with numeric values
storage.mode(signature_matrix) <- "numeric"
storage.mode(bulk_data) <- "numeric"

# Run DWLS deconvolution using solveDampenedWLS
## Initialize an empty matrix to store results
deconvolution_results <- matrix(
  NA, nrow = ncol(bulk_data), ncol = ncol(signature_matrix)
)
rownames(deconvolution_results) <- colnames(bulk_data)
colnames(deconvolution_results) <- colnames(signature_matrix)

# Loop through each bulk sample and perform deconvolution
for (i in 1:ncol(bulk_data)) {
  deconvolution_results[i, ] <- DWLS::solveDampenedWLS(signature_matrix, bulk_data[, i])
}

# Convert results to a data frame for easier inspection
deconvolution_results <- as.data.frame(deconvolution_results)

# Save the results
write.table(
  deconvolution_results,
  here(results_dir, "dwls_results.txt"),
  sep="\t",
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE
)

# Plot results as a heatmap
pdf(
  here(results_dir, "dwls_heatmap.pdf"), 
  height = 6, width = 9
)
deconvolution_matrix <- t(deconvolution_results)

heatmap <- Heatmap(deconvolution_matrix,
                   name = "Cell Type Proportion",
                   col = colorRampPalette(c("white", "aquamarine4"))(100),
                   cluster_rows = FALSE,            
                   cluster_columns = FALSE,        
                   show_row_names = TRUE,
                   show_column_names = TRUE,
                   row_names_gp = gpar(fontsize = 8),
                   row_names_side = "left",
                   column_names_gp = gpar(fontsize = 8),
                   column_names_side = "top",
                   heatmap_legend_param = list(title = "Proportion", 
                                               color_bar = "continuous", 
                                               legend_direction = "horizontal",
                                               legend_width = unit(5, "cm"))
)
draw(heatmap, heatmap_legend_side = "bottom")
dev.off()

# End of script