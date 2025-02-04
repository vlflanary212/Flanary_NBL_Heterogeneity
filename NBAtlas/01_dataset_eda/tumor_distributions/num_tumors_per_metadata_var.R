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
  "NB_ITH", "NBAtlas", "01_dataset_eda", "dataset_overview", "tumor_distributions"
)

# Plot the number of tumors per metadata variable as bar plots
## Plot in a single pdf
pdf(
    here(results_dir, paste0("Num_Tumors_by_Metadata.pdf")),
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
  
  # Adjust y-axis limit
  max_y <- max(metadata_summary$n, na.rm = TRUE)
  adjusted_max_y <- max_y * 1.1
  
  # Create and print the plot
  print(
    ggplot(metadata_summary, aes(x = metadata_value, y = n, fill = metadata_value)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = n), vjust = -0.5, size = 4) +
      ylim(0, adjusted_max_y) +
      labs(
        title = paste("Number of Samples by", i),
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
}

dev.off()

## Plot as individual files
for (i in colnames(metadata)) {
  pdf(
    here(results_dir, paste0("Num_Tumors_by_", i, ".pdf")),
    height = 4, width = 6
  )
  
  # Select the relevant columns and summarize data
  metadata_summary <- metadata[, c("Sample", i)] |>  # Select Sample and current metadata column
    distinct() |>
    group_by(.data[[i]]) |>   # Group by current metadata variable
    summarise(n = n(), .groups = "drop")  # Count the number of Samples per metadata value
  
  # Assign consistent column names
  colnames(metadata_summary) <- c("metadata_value", "n")
  
  # Adjust y-axis limit
  max_y <- max(metadata_summary$n, na.rm = TRUE)
  adjusted_max_y <- max_y * 1.1
  
  # Create and print the plot
  print(
    ggplot(metadata_summary, aes(x = metadata_value, y = n, fill = metadata_value)) +
      geom_bar(stat = "identity") +
      geom_text(aes(label = n), vjust = -0.5, size = 4) +
      ylim(0, adjusted_max_y) +
      labs(
        title = paste("Number of Samples by", i),
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
}

# End of script