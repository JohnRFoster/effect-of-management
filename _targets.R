# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed.

# Set target options:
tar_option_set(
  packages = c(
    "tibble",
    "dplyr",
    "tidyr",
    "readr",
    "ggplot2",
    "ggpubr",
    "stringr",
    "lubridate"
  ) # Packages that your targets need for their tasks.

  # Set other options as needed.
)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()

# Replace the target list below with your own:
list(
  tar_target(
    name = model_data_last,
    command = read_rds(file.path(data_dir, last_iter, "modelData.rds")),
    format = "file"
    # format = "qs" # Efficient storage for general data objects.
  ),

  tar_target(
    name = density_last,
    command = read_rds(file.path(data_dir, last_iter, "densitySummaries.rds")),
    format = "file"
  )
)
