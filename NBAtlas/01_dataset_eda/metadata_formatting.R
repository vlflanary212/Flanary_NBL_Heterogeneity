# Author: Victoria Flanary
# Date: 250122
# Objective: Extract metadata from the NBAtlas for plotting

# Set-up
## Load  packages
library(Seurat)
library(tidyverse)
library(here)

## Load data
seurat_obj <- readRDS(
  here("NB_ITH", "NBAtlas", "data", "subset", "00_init_obj.rds")
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

# Continuous technical: nCount_RNA, nFeature_RNA, percent_mito
# Discrete technical: Study, Assay, Platform, Cell_condition
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

## Split Patient_No by Number and Timepoint, and remove Timepoint (already recorded)
## More easily lets you see which patients have more than one associated sample
metadata <- metadata |>
  separate(col = "Patient_No", into = c("Patient_No", "T")) |>
  dplyr::select(all_of(setdiff(colnames(metadata), "T")))

# Save the metadata in the Seurat object
seurat_obj@meta.data <- metadata

# Save the edited Seurat object
saveRDS(
  seurat_obj,
  here("NB_ITH", "NBAtlas", "data", "subset", "00_formatted_obj.rds")
)

# Isolate metadata of interest for plotting
# Keeping all but orig.ident and the continuous technical variables
vars <- c("Study", "Assay", "Platform", "Cell_condition", "Sample", "Patient_No", 
          "Gender", "MYCN_amplification", "Risk", "INSS_stage", "Timepoint", "Cell_type")

metadata_filt <- select(metadata, all_of(vars))

## Save the filtered metadata
write_csv(
  metadata_filt, 
  here("NB_ITH", "NBAtlas", "metadata_filt.csv")
)

# End of Script
