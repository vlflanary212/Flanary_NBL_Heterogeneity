# Author: Victoria Flanary
# Date: 250129
# Objective: Generate counts matrix for PCGP neuroblastoma RNA-seq samples.

# Load packages
library(tidyverse)
library(data.table)
library(here)

# Define file paths
results_dir <- here("PCGP_RNA", "01_counts_matrices")

## Read in the raw counts
merged_counts <- read.table(here(results_dir, "raw_counts.txt"))

# Save rownames for later
gene_ids <- rownames(merged_counts)

# Normalize the Counts
## Use the same normalization method as the Seurat default.

## Convert all columns in the counts matrix to numeric
merged_counts <- apply(merged_counts, 2, as.numeric)

## Calculate library size for each sample
library_sizes <- colSums(merged_counts)

## Normalize counts by library size and scale factor (10,000)
scale_factor <- 10000
normalized_counts <- sweep(merged_counts, 2, library_sizes, "/") * scale_factor

## Log-transform the normalized counts
log_normalized_counts <- log1p(normalized_counts)

## Add gene ids as row names
log_normalized_counts <- as.data.frame(log_normalized_counts)
rownames(log_normalized_counts) <- gene_ids

## Save the normalized counts
write.table(
  log_normalized_counts,
  here(results_dir, "log_normalized_counts.txt"),
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE
)

# End of script
