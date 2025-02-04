# Author: Victoria Flanary
# Date: 250124
# Objective: Visualize the proportion of each cell type representated across
# different metadata variables

# Set-up
## Load  packages
library(tidyverse)
library(here)

## Load metadata
metadata <- read_csv(
  here("NB_ITH", "NBAtlas", "metadata_filt.csv")
) |> as.data.frame()

## Set results directory
results_dir <- here(
  "NB_ITH", "NBAtlas", "01_dataset_eda", "cell_type_proportions"
)

# Plot the cell type distributions across metadata variables
## Plot in a single pdf
pdf(
  here(results_dir, "all_cells", paste0("Cell_Types_by_Metadata.pdf")),
  height = 4, width = 6
)

# Loop through metadata columns
for (i in setdiff(colnames(metadata), "Cell_type")) {
  # Select relevant columns and calculate proportions
  proportion_data <- metadata[, c("Cell_type", i)] %>%
    group_by(.data[[i]], Cell_type) %>% # Group by metadata variable and cell type
    summarize(count = n(), .groups = "drop") %>% # Count cells
    group_by(.data[[i]]) %>% # Regroup by metadata variable
    mutate(percentage = count / sum(count)) %>% # Calculate percentage
    ungroup()
  
  # Create the stacked bar plot
  p <- ggplot(
    proportion_data,
    aes(x = .data[[i]], y = percentage, fill = Cell_type)
  ) +
    geom_bar(stat = "identity") +
    scale_y_continuous(labels = scales::percent) + # Format y-axis as percents
    labs(
      title = paste("Cell Type Distributions by", i),
      x = NULL,
      y = "Percentage",
      fill = "Cell Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels
    )
  
  # Print the plot to the PDF device
  print(p)
}

dev.off()

## Plot as individual files
for (i in setdiff(colnames(metadata), "Cell_type")) {
  pdf(
    here(results_dir, "all_cells", paste0("Cell_Types_by_", i, ".pdf")),
    height = 4, width = 6
  )
  
   # Select relevant columns and calculate proportions
  proportion_data <- metadata[, c("Cell_type", i)] %>%
    group_by(.data[[i]], Cell_type) %>% # Group by metadata variable and cell type
    summarize(count = n(), .groups = "drop") %>% # Count cells
    group_by(.data[[i]]) %>% # Regroup by metadata variable
    mutate(percentage = count / sum(count)) %>% # Calculate percentage
    ungroup()
  
  # Create the stacked bar plot
  p <- ggplot(
    proportion_data,
    aes(x = .data[[i]], y = percentage, fill = Cell_type)
  ) +
    geom_bar(stat = "identity") +
    scale_y_continuous(labels = scales::percent) + # Format y-axis as percents
    labs(
      title = paste("Cell Type Distributions by", i),
      x = NULL,
      y = "Percentage",
      fill = "Cell Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels
    )
  
  # Print the plot to the PDF device
  print(p)
  
  dev.off()
}

# Immune cells only
## Hard to see changes in immune and strome cell populations with how many more
## neuroendocrine cells there are - plot separately
tme_metadata <- metadata |>
  filter(Cell_type != "Neuroendocrine")

## Plot in a single pdf
pdf(
  here(results_dir, "tme_cells", paste0("Cell_Types_by_Metadata.pdf")),
  height = 4, width = 6
)

# Loop through metadata columns
for (i in setdiff(colnames(tme_metadata), "Cell_type")) {
  # Select relevant columns and calculate proportions
  proportion_data <- tme_metadata[, c("Cell_type", i)] %>%
    group_by(.data[[i]], Cell_type) %>% # Group by metadata variable and cell type
    summarize(count = n(), .groups = "drop") %>% # Count cells
    group_by(.data[[i]]) %>% # Regroup by metadata variable
    mutate(percentage = count / sum(count)) %>% # Calculate percentage
    ungroup()
  
  # Create the stacked bar plot
  p <- ggplot(
    proportion_data,
    aes(x = .data[[i]], y = percentage, fill = Cell_type)
  ) +
    geom_bar(stat = "identity") +
    scale_y_continuous(labels = scales::percent) + # Format y-axis as percents
    labs(
      title = paste("TME Cell Type Distributions by", i),
      x = NULL,
      y = "Percentage",
      fill = "Cell Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels
    )
  
  # Print the plot to the PDF device
  print(p)
}

dev.off()

## Plot as individual files
for (i in setdiff(colnames(tme_metadata), "Cell_type")) {
  pdf(
    here(results_dir, "tme_cells", paste0("Cell_Types_by_", i, ".pdf")),
    height = 4, width = 6
  )
  
   # Select relevant columns and calculate proportions
  proportion_data <- tme_metadata[, c("Cell_type", i)] %>%
    group_by(.data[[i]], Cell_type) %>% # Group by metadata variable and cell type
    summarize(count = n(), .groups = "drop") %>% # Count cells
    group_by(.data[[i]]) %>% # Regroup by metadata variable
    mutate(percentage = count / sum(count)) %>% # Calculate percentage
    ungroup()
  
  # Create the stacked bar plot
  p <- ggplot(
    proportion_data,
    aes(x = .data[[i]], y = percentage, fill = Cell_type)
  ) +
    geom_bar(stat = "identity") +
    scale_y_continuous(labels = scales::percent) + # Format y-axis as percents
    labs(
      title = paste("TME Cell Type Distributions by", i),
      x = NULL,
      y = "Percentage",
      fill = "Cell Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1) # Rotate x-axis labels
    )
  
  # Print the plot to the PDF device
  print(p)
  
  dev.off()
}

# End of script