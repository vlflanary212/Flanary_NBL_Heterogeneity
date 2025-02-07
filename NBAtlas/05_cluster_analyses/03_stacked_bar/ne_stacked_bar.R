# Author: Victoria Flanary
# Date: 250205
# Objective: Visualize distribution of tumor cell clusters across patient
# risk group and treatment timepoint.

# Set-up
## Load  packages
library(Seurat)
library(tidyverse)
library(here)

## Set filepaths
data_dir <- here("NBAtlas", "data", "alldata")
results_dir <- here("NBAtlas", "05_cluster_analyses", "03_stacked_bar")

## Load data
seurat_obj <- readRDS(file.path(data_dir, "09_ne_subset.rds"))

# Set current ident to Cell_type
seurat_obj$Cell_type <- seurat_obj@active.ident

## Extract metadata
metadata <- seurat_obj@meta.data

# Format metadata to add Westermann subgroups (combines Risk and MYCN status)
### Edit Risk column 
unique(metadata$Risk)  # "Low/Intermediate" "High"             "_unknown"   
metadata$Risk <- ifelse(metadata$Risk == "Low/Intermediate", "Low", metadata$Risk)
metadata$Risk <- factor(metadata$Risk, levels = c("High", "Low", "_unknown"))
unique(metadata$Risk) 

### Plot MYCN status for different Risk levels
#### Group and summarize data
proportion_data <- metadata |>
  group_by(Risk, MYCN_amplification) |>
  summarize(count = n(), .groups = "drop") |>
  group_by(Risk) |>
  mutate(percentage = count / sum(count)) |>
  ungroup()

#### Create the stacked bar plot
pdf(
  file.path(results_dir, "mycn_status_by_risk.pdf")
)

p <- ggplot(
  proportion_data,
  aes(x = Risk, y = percentage, fill = MYCN_amplification)
) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "MYCN Amplification by Risk",
    x = "Risk",
    y = "Percentage",
    fill = "MYCN Amplification"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1) 
  )

print(p)

dev.off()

### Add Westermann Groups
metadata <- metadata |> 
  mutate(Group = ifelse(
    Risk != "unknown" & MYCN_amplification != "unknown",
    paste0(Risk, "-Risk", "_", "MYCN-", MYCN_amplification), "Unknown")) |>
  relocate(Group, .after = Risk)

metadata$Group <- ifelse(metadata$Group == "_unknown-Risk_MYCN-_unknown", "Unknown", metadata$Group)

### Plot the number of tumors in each Westermann Group on bar plot
#### Select the relevant columns and summarize data
metadata_summary <- metadata[, c("Sample", "Group")] |>  # Select Sample and Group column
  distinct() |>
  group_by(Group) |>   # Group by Group
  summarise(n = n(), .groups = "drop")  # Count the number of Samples per Group value

#### Assign consistent column names
colnames(metadata_summary) <- c("metadata_value", "n")

#### Adjust y-axis limit
max_y <- max(metadata_summary$n, na.rm = TRUE)
adjusted_max_y <- max_y * 1.1

#### Create and print the plot
pdf(
  file.path(results_dir, "num_tumors_per_westermann_group.pdf"),
  height = 4, width = 6
)

print(
  ggplot(metadata_summary, aes(x = metadata_value, y = n, fill = metadata_value)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = n), vjust = -0.5, size = 4) +
    ylim(0, adjusted_max_y) +
    labs(
      title = "Number of Samples by Westermann Group",
      x = NULL,
      y = "Number of Samples"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),  # Center the plot title
      legend.position = "none",  # Remove the plot legend
      axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis text
    )
)

dev.off()

### Save the tumor seurat object with updated metadata
seurat_obj@meta.data <- metadata
saveRDS(seurat_obj, file.path(data_dir, "10_ne_westermann.rds"))

# Determine the distribution of tumor cell clusters across Westermann groups
## Calculate proportion data
proportion_data <- metadata |>
  group_by(Group, Cell_type) |>
  summarize(count = n(), .groups = "drop") |>
  group_by(Group) |>
  mutate(percentage = count / sum(count)) |>
  ungroup()

## Plot the results on a stacked bar plot
pdf(
  file.path(results_dir, "ne_clusters_by_westermann_group.pdf")
)

p <- ggplot(
  proportion_data,
  aes(x = Group, y = percentage, fill = Cell_type)
) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Tumor Cell Clusters by Westermann Group",
    x = "Westermann Group",
    y = "Percentage",
    fill = "Tumor Cluster"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1) 
  )

print(p)

dev.off()

# Determine tumor cell clusters across treatment timepoints
## Edit Timepoint
metadata$Timepoint <- factor(metadata$Timepoint, levels = c("pre-treatment", "post-treatment", "relapse", "_unknown"))

## Calculate proportion data
proportion_data <- metadata |>
  group_by(Timepoint, Cell_type) |>
  summarize(count = n(), .groups = "drop") |>
  group_by(Timepoint) |>
  mutate(percentage = count / sum(count)) |>
  ungroup()

## Plot the results on a stacked bar plot
pdf(
  file.path(results_dir, "ne_clusters_by_treatment_timepoint.pdf")
)

p <- ggplot(
  proportion_data,
  aes(x = Timepoint, y = percentage, fill = Cell_type)
) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Tumor Cell Clusters by Treatment Timepoint",
    x = "Timepoint",
    y = "Percentage",
    fill = "Tumor Cluster"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1) 
  )

print(p)

dev.off()

# End of script