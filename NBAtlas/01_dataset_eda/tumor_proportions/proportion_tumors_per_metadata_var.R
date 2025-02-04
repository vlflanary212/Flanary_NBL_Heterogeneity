# Author: Victoria Flanary
# Date: 240122
# Objective: Visualize the number of tumors belonging to each metadata
# variable in the NBAtlas

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
  "NB_ITH", "NBAtlas", "01_dataset_eda", "dataset_overview", "tumor_proportions"
)

# Plot the number of tumors per metadata variable as bar plots
## Plot in a single pdf
pdf(
  here(results_dir, "Proportion_Tumors_by_Metadata.pdf"),
  height = 4, width = 6
)

for (i in colnames(metadata)) {
  # Select the relevant columns and summarize data
  metadata_summary <- metadata[, c("Sample", i)] |>  # Select Sample and current metadata column
    distinct() |>
    group_by(.data[[i]]) |>   # Group by current metadata variable
    summarise(n = n(), .groups = "drop")  # Count the number of Samples per metadata value
  
  # Assign consistent column names
  colnames(metadata_summary) <- c("metadata_value", "n")
  
  # Calculate percentages for annotation
  metadata_summary <- metadata_summary %>%
    mutate(
      percentage = n / sum(n) * 100, # Calculate percentage
      label = paste0(n, " (", round(percentage, 1), "%)")
    ) # Combine count and percentage
  
  print(
    # Create the pie chart
    ggplot(metadata_summary, aes(x = "", y = n, fill = metadata_value)) +
      geom_bar(stat = "identity", width = 1) + # Bar chart for polar transformation
      coord_polar(theta = "y") + # Transform to pie chart
      geom_text(aes(label = label), # Add annotations
        position = position_stack(vjust = 0.5), size = 4
      ) + # Center labels
      labs(
        title = paste("Proportion of Samples by", i),
        x = NULL, # Remove x-axis label
        y = NULL # Remove y-axis label
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5), # Center the plot title
        axis.text.x = element_blank(), # Remove x-axis text
        axis.ticks = element_blank(), # Remove axis ticks
        panel.grid = element_blank(), # Remove grid lines
        legend.title = element_blank() # Remove legend title
      )
  )
}

dev.off()

## Plot as individual files
for (i in colnames(metadata)) {
  pdf(
  here(results_dir, paste0("Proportion_of_Tumors_by_", i, ".pdf")),
  height = 4, width = 6
  )
  
  # Select the relevant columns and summarize data
  metadata_summary <- metadata[, c("Sample", i)] |>  # Select Sample and current metadata column
    distinct() |>
    group_by(.data[[i]]) |>   # Group by current metadata variable
    summarise(n = n(), .groups = "drop")  # Count the number of Samples per metadata value
  
  # Assign consistent column names
  colnames(metadata_summary) <- c("metadata_value", "n")
  
  # Calculate percentages for annotation
  metadata_summary <- metadata_summary %>%
    mutate(
      percentage = n / sum(n) * 100, # Calculate percentage
      label = paste0(n, " (", round(percentage, 1), "%)")
    ) # Combine count and percentage
  
  print(
    # Create the pie chart
    ggplot(metadata_summary, aes(x = "", y = n, fill = metadata_value)) +
      geom_bar(stat = "identity", width = 1) + # Bar chart for polar transformation
      coord_polar(theta = "y") + # Transform to pie chart
      geom_text(aes(label = label), # Add annotations
        position = position_stack(vjust = 0.5), size = 4
      ) + # Center labels
      labs(
        title = paste("Proportion of Samples by", i),
        x = NULL, # Remove x-axis label
        y = NULL # Remove y-axis label
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5), # Center the plot title
        axis.text.x = element_blank(), # Remove x-axis text
        axis.ticks = element_blank(), # Remove axis ticks
        panel.grid = element_blank(), # Remove grid lines
        legend.title = element_blank() # Remove legend title
      )
  )
  
  dev.off()
}

# ENd of script