# Author: Victoria Flanary
# Date: 250207
# Objective: Subset neuroendocrine clusters and generate a pseudobulked
# raw counts matrix for comparative analyses with other bulk RNA-seq data.

# Load packages
library(Seurat)
library(Matrix)
library(tidyverse)
library(here)

# Set filepaths
data_dir <- here("NBAtlas", "data", "alldata")
results_dir <- here("NBAtlas", "06_pseudobulking")

# Load data
seurat_obj <- readRDS(here(data_dir, "08_final_anno.rds"))

# Subset neuroendocrine cells
seurat_obj <- subset(seurat_obj, idents = grep("^NE_", Idents(seurat_obj), value = TRUE))

# Pseudo-bulking
## Aggregate gene expression per cell type
agg_exp <- AggregateExpression(
  seurat_obj,
  assay = "RNA",
  slot = "counts",  # retrieve raw counts
  group.by = "ident",
  return.seurat = FALSE
)

# Convert result to a matrix
agg_exp <- as.matrix(agg_exp[["RNA"]])

# Normalize by number of cells per cell type annotation
cell_counts <- table(seurat_obj@active.ident)
agg_exp_normalized <- sweep(agg_exp, 2, cell_counts, FUN = "/")

# Move gene names to a separate column
genes <- rownames(agg_exp_normalized)

agg_exp_df <- data.frame(
  gene_name = genes,
  agg_exp_normalized,
  row.names = NULL
)

# Save the aggregated counts
write_csv(
  agg_exp_df,
  here(results_dir, "ne_agg_norm_counts.csv")
)

# Convert gene names to gene ids
## Read in gtf annotations
gtf_file <- "/data/project/sen-lab/genome/hg38/gencode.v22.annotation.gtf"
gtf <- as.data.frame(rtracklayer::import(gtf_file))

## Filter gtf file for gene entries in agg_exp
gene_mappings <- gtf |>
  filter(type == "gene") |>
  dplyr::select(gene_id, gene_name) |>
  distinct() |>
  filter(gene_name %in% agg_exp_df$gene_name)

## Rm duplicates
gene_mappings <- gene_mappings[!duplicated(gene_mappings$gene_name), ]

## Subset agg_exp to only include mappable genes
signature_matrix <- agg_exp_df[agg_exp_df$gene_name %in% gene_mappings$gene_name, ]

## Match row order in the agg counts and gene_mapping
signature_matrix <- signature_matrix[match(gene_mappings$gene_name, signature_matrix$gene_name), ]

## Update rownames to gene ids without version numbers
signature_matrix$gene_id <- sub("\\..*", "", gene_mappings$gene_id)
signature_matrix <- relocate(signature_matrix, gene_id, .before = gene_name)
rownames(signature_matrix) <- NULL

## Save the signature matrix
write_csv(
  signature_matrix,
  here(results_dir, "single_cell_tumor_pseudobulk_matrix.csv")
)

# End of script