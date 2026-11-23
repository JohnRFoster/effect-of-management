# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline

# Load packages required to define the pipeline:
source("packages.R") |>
  suppressPackageStartupMessages()

# Run the R scripts in the R/ folder with your custom functions:
tar_source()

# project setup
setup_targets <- tar_plan(
  # for the bootstrap analysis
  bootstrap_parameters = list(
    n_ens = 50, # the number of ensemble members to use in each bootstrap iteration
    sim_days = 365 * 5, # the number of days to simulate in each bootstrap iteration
    n_sim = 50 # the number of bootstrap iterations to run for each state
  ),

  # directories
  d = list(
    plot_dir = "plots",
    data_dir = "data",
    data_store = "../data-store",
    state_fits = file.path("stateFit", "states2-2026-03-25")
  )
)

# data files
data_targets <- tar_plan(
  tar_file_read(
    name = model_data,
    command = file.path(d$data_store, d$state_fits, "model_data.rds"),
    read = get_model_data(file = !!.x) |> fix_state_names()
  ),

  tar_file_read(
    name = densities,
    command = file.path(d$data_store, d$state_fits, "all_density.rds"),
    read = get_density(file = !!.x) |> fix_state_names()
  ),

  tar_file_read(
    name = parameter_samples,
    command = file.path(d$data_store, d$state_fits, "all_params_samples.rds"),
    read = read_rds(file = !!.x)
  ),

  tar_file_read(
    name = df_params,
    command = file.path(d$data_store, d$state_fits, "all_params.rds"),
    read = read_rds(file = !!.x) |> fix_state_names2()
  ),

  tar_file_read(
    name = df_priors,
    command = file.path(d$data_store, d$state_fits, "priors_df.rds"),
    read = read_rds(file = !!.x)
  ),

  tar_file_read(
    name = abundance_samples,
    command = file.path(
      d$data_store,
      d$state_fits,
      "all_abundance_samples.rds"
    ),
    read = read_rds(file = !!.x) |> fix_state_names()
  )
)

analysis_targets <- tar_plan(
  data_mis = data_prep(model_data, densities),
  realized_growth = calculate_realized_r(data_mis),
  mu_lambda = mean_growth(realized_growth),
  actual_percent_take = prop_take_by_method(model_data, densities),
  bootstrap_mgmt = run_bootstrap(
    parameter_samples,
    abundance_samples,
    data_mis,
    bootstrap_parameters
  ),
  epr_data = make_epr_data(realized_growth, actual_percent_take),
  density_cis = density_summaries(densities, model_data)
)

figure_targets <- tar_plan(
  fig_lambda_hat_by_state = plot_mean_lambdas(mu_lambda),
  fig_prop_take = plot_percent_take(actual_percent_take),
  fig_density_prevented = plot_bootstrap(bootstrap_mgmt),
  fig_epr = plot_effort(epr_data),
  fig_state_removal_rates_by_method = plot_state_intercepts_rem_rate(
    df_params,
    df_priors
  ),
  fig_states_road_density = plot_states_road_density(
    df_params,
    df_priors
  ),
  fig_states_ruggedness = plot_states_ruggedness(
    df_params,
    df_priors
  ),
  fig_states_canopy_cover = plot_states_canopy_cover(
    df_params,
    df_priors
  ),
  fig_effort_saturation = plot_effort_saturation(
    df_params,
    df_priors
  ),
  fig_effort_unique_area = plot_effort_unique_area(
    df_params,
    df_priors
  ),
  fig_effort_scaling = plot_effort_scaling(
    df_params,
    df_priors
  ),
  fig_demographic_parameters = plot_demographic_parameters(
    df_params,
    df_priors
  )
)

# List targets ----
# comment/uncomment to include in tar_make() or tar_visnetwork()
list(
  setup_targets,
  data_targets,
  analysis_targets,
  figure_targets
)
