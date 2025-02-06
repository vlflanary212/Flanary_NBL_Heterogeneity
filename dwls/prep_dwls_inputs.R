# Author: Victoria Flanary
# Date: 250130
# Objective: Prepare inputs for DWLS deconvolution

# Load packages
library(tidyverse)
library(here)

# Load data
bulk_counts <- read.table(
  here("NB_ITH", "PCGP_RNA", "counts_matrices", "log_normalized_counts.txt")
)

sc_counts <- read_csv(
  here("NB_ITH", "NBAtlas", "dwls_signature_matrix", "signature_matrix.csv")
) |> as.data.frame()

# Format Single-Cell Signatures
## Filter to keep genes in both matrices
overlapping_genes <- intersect(rownames(bulk_counts), sc_counts$gene_id)
sc_counts_filt <- subset(sc_counts, gene_id %in% overlapping_genes)

## Move gene ids to rownames
rownames(sc_counts_filt) <- sc_counts_filt$gene_id

## Rm gene annotation columns
sc_counts_filt <- sc_counts_filt[, setdiff(colnames(sc_counts_filt), c("gene_id", "gene_name"))]

## Save the filtered dataframe
write.table(
  sc_counts_filt,
  here("NB_ITH", "PCGP_RNA", "dwls", "inputs", "single_cell_ref.txt"),
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE
)

# Format the Bulk Counts
## Filter to keep genes in both matrices
bulk_counts_filt <- subset(bulk_counts, rownames(bulk_counts) %in% overlapping_genes)

## Match row order to genes in sc_count
bulk_counts_ordered <- bulk_counts_filt[match(rownames(sc_counts_filt), rownames(bulk_counts_filt)), ]
identical(rownames(sc_counts_filt), rownames(bulk_counts_ordered)) # TRUE

## Save the ordered dataframe
write.table(
  bulk_counts_ordered,
  here("NB_ITH", "PCGP_RNA", "dwls", "inputs", "bulk_matrix.txt"),
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE,
  quote = FALSE
)

# End of script