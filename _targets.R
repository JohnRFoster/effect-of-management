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
    name = model_data_last_file,
    command = file.path(data_dir, last_iter, "modelData.rds"),
    format = "file"
    # format = "qs" # Efficient storage for general data objects.
  ),

  tar_target(
    name = model_data_last,
    command = get_rds(model_data_last_file)
  ),

  tar_target(
    name = density_last_file,
    command = file.path(data_dir, last_iter, "densitySummaries.rds"),
    format = "file"
  ),

  tar_target(
    name = density_last,
    command = get_rds(density_last_file)
  ),

  tar_target(
    name = parameter_last_file,
    command = file.path(data_dir, last_iter, "posteriorSamples.rds"),
    format = "file"
  ),

  tar_target(
    name = parameter_last,
    command = get_rds(parameter_last_file)
  ),

  tar_target(
    name = data_mis,
    command = data_prep(model_data_last, density_last, cutoff_date)
  ),

  tar_target(
    name = realized_growth,
    command = calculate_realized_r(data_mis)
  ),

  tar_target(
    name = mu_lambda,
    command = mean_growth(realized_growth)
  ),

  tar_target(
    name = plot_lambda_hat_by_state,
    command = plot_mean_lambdas(mu_lambda)
  ),

  tar_target(
    name = actual_percent_take,
    command = prop_take_by_method(model_data_last, mu_lambda, density_last)
  ),

  tar_target(
    name = theoretical_percent_take,
    command = prop_to_remove(realized_growth)
  ),

  tar_target(
    name = plot_prop_take,
    command = plot_percent_take(actual_percent_take, theoretical_percent_take)
  )
)
