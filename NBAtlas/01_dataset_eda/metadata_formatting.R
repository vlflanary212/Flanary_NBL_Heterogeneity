# Author: Victoria Flanary
# Date: 250122
# Objective: Extract metadata from the NBAtlas for plotting

# Set-up
## Load  packages
library(Seurat)
library(tidyverse)
library(here)

# Set filepaths
data_dir <- here("NBAtlas", "data", "alldata")
results_dir <- here("NBAtlas", "01_dataset_eda")

## Load data
seurat_obj <- readRDS(
  file.path(data_dir, "00_init_obj.rds")
)

# Split the object by batch
seurat_obj@meta.data <- mutate(
  seurat_obj@meta.data,
  Batch = paste(Study, Assay, Platform, Cell_condition, sep = "_")
)

## Save batches in a text file
write.table(
  unique(seurat_obj$Batch),
  here("NBAtlas", "batches.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# Extract metadata
metadata <- seurat_obj@meta.data

# Note: would save and track all the metadata at this point, but the file size
# Exceeds GitHub's limit of 50 MB

#  Note the types of metadata available
colnames(metadata)
# [1] "orig.ident"         "nCount_RNA"         "nFeature_RNA"       "percent_mito"      
# [5] "Study"              "Assay"              "Platform"           "Sample"            
# [9] "Patient_No"         "Timepoint"          "INSS_stage"         "MYCN_amplification"
# [13] "Gender"             "Risk"               "Cell_condition"     "Cell_type" 
# [17] "Batch"

# Continuous technical: nCount_RNA, nFeature_RNA, percent_mito
# Discrete technical: Study, Assay, Platform, Cell_condition (included in Batch)
# Discrete biological: orig.ident, Sample, Patient_No, Timepoint, INSS_stage,
# MYCN_amplification, Gender, Risk, Cell_type

# orig.ident seems to refer to original cell type annotations given by the 
# NBAtlas authors instead of patient or tumor sample
identical(metadata$orig.ident, metadata$Cell_type)  # TRUE

# Edit the metadata in the object to your desired formatting
## Combine Slyper2020_cell and Slyper2020_nucleus - split by batch later
metadata <- mutate(
  metadata, 
  Study = case_when(grepl("Slyper2020", Study) ~ "Slyper2020", TRUE ~ Study)
  )

## Edit Patient_No
## More easily lets you see which patients have more than one associated sample
metadata$Patient_No <- gsub("^([0-9]+).*", "\\1", metadata$Patient_No)
metadata$Patient_No <- as.factor(metadata$Patient_No)

## Factorize and re-order levels for Timepoint
metadata$Timepoint <- factor(metadata$Timepoint, levels = c("pre-treatment", "post-treatment", "relapse", "unknown"))

# Save the metadata in the Seurat object
seurat_obj@meta.data <- metadata

# Save the edited Seurat object
saveRDS(
  seurat_obj,
  file.path(data_dir, "00_formatted_obj.rds")
)

# Isolate metadata of interest for plotting
# Keeping all but orig.ident and the continuous technical variables
vars <- c("Batch", "Study", "Assay", "Platform", "Cell_condition", "Sample", "Patient_No", 
          "Gender", "MYCN_amplification", "Risk", "INSS_stage", "Timepoint", "Cell_type")

metadata_filt <- select(metadata, all_of(vars))

## Save the filtered metadata
write_csv(
  metadata_filt, 
  file.path(results_dir, "metadata_filt.csv")
)

# End of Script
