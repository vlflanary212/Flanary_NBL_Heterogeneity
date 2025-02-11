# Author: Victoria Flanary
# Date: 250129
# Objective: Generate counts matrix for PCGP neuroblastoma RNA-seq samples.

# Load packages
library(tidyverse)
library(data.table)
library(here)

# Raw Counts
## Define file paths
data_dir <- "/data/project/sen-lab/external/processed/St_Jude/RNAseq"
sample_file <- file.path(data_dir, "sample.txt")
rnaseq_dir <- file.path(data_dir, "realigned_bam")
results_dir <- here("PCGP_RNA", "01_counts_matrices")

## Read in the sample list
samples <- readLines(sample_file)

## Initialize an empty object to store the merged counts
merged_counts <- NULL

## Get counts matrix for each sample in sample_list
for (sample in samples) {
  # Define the directory and count file for the sample
  sample_dir <- file.path(rnaseq_dir, sample)
  count_file <- file.path(sample_dir, "counts.txt")
  
  # Check if the count file exists
  if (file.exists(count_file)) {
    # Read in the count file
    sample_counts <- read_tsv(
      count_file, 
      col_names = FALSE, 
      show_col_types = FALSE
    )
    sample_counts <- sample_counts[-1, ] # Rm old colnames
    colnames(sample_counts) <- c("gene", sample)
    # First column is 'gene', second column named after the sample
    
    # Merge the data into the main matrix
    if (is.null(merged_counts)) {
      merged_counts <- sample_counts
    } else {
      merged_counts <- full_join(merged_counts, sample_counts, by = "gene")
    }
  } else {
    cat("Warning: Counts file not found for sample", sample, "\n")
  }
  
  # Print progress
  print(paste0(sample, " added to merged_counts"))
}

## Set rownames to formatted gene ids
# Remove version numbers from ensembl ids
merged_counts$gene <- sub("\\.\\d+$", "", merged_counts$gene)

# Save the formatted gene ids
gene_ids <- merged_counts$gene

# Remove the gene column
merged_counts <- merged_counts[, colnames(merged_counts) != "gene"]

# Add gene ids as rownames
merged_counts <- as.data.frame(merged_counts)
rownames(merged_counts) <- gene_ids

## Save the merged counts
write.table(
  merged_counts,
  file = here(results_dir, "raw_counts.txt"),
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE
)

# End of script
