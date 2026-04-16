#---------
#
# Workflow for checking property-level Bayes model fits of the MIS data
# - combines mcmc chunks
# - checks for convergence
#
#---------

library(dplyr)
library(tidyr)
library(readr)
library(parallel)
library(coda)
library(ggplot2)

source("nimble/functions_nimble.R")
source("nimble/functions_prep_nimble.R")

config_name <- "hpc_dev"
config <- config::get(config = config_name)

out_dir <- config$out_dir
files_in_out_dir <- list.files(out_dir)
n_props_fit <- as.numeric(stringr::str_extract(files_in_out_dir, "\\d*(?=\\D)"))
last_fit <- max(n_props_fit)
last_fit <- 1
dest_mcmc <- file.path(out_dir, paste0(last_fit, "_mcmc"))
dest_posterior <- file.path(out_dir, paste0(last_fit, "_posterior"))

# get samples
mcmc_list <- collate_mcmc_chunks(dest_mcmc, start = 1)
params_mcmc_list <- mcmc_list$params

# calculate psrf (convergence stat) and effective sample size
message("==== Pre-burnin diagnostics ====")
diagnostic <- continue_mcmc(
  params_mcmc_list,
  effective_size = 1000,
  max_psrf = 15,
  verbose = TRUE
)


total_iter <- nrow(params_mcmc_list[[1]])
n_chains <- length(params_mcmc_list)
GBR <- gelman.plot(params_mcmc_list)
burnin <- GBR$last.iter[
  tail(which(apply(GBR$shrink[,, 2] > 1.1, 1, any)), 1) + 1
]
message("Burnin: ", burnin)
if (is.na(burnin)) {
  burnin <- round(total_iter / 2)
}
params_burnin <- window(params_mcmc_list, start = burnin)

# calculate psrf (convergence stat) and effective sample size
message("==== Post-burnin diagnostics ====")
diagnostic <- continue_mcmc(
  params_burnin,
  effective_size = 1000,
  max_psrf = 15,
  verbose = TRUE
)

# X11 is not on the HPC, so need to use ggplot2 to create traceplots
# which means we need to create a tibble of samples and identify chains and iterations
posterior <- tibble()
message("Create posterior tibble...")
pb <- txtProgressBar(min = 1, max = n_chains, style = 1)
for (i in seq_len(n_chains)) {
  chain_i <- as.matrix(params_mcmc_list[[i]]) |>
    as_tibble() |>
    mutate(chain = i)
  posterior <- bind_rows(posterior, chain_i)

  setTxtProgressBar(pb, i)
}
close(pb)
message("  done")

# get a random sample of the posterior and save
posterior_burnin <- params_burnin |>
  as.matrix()

draws <- sample.int(nrow(posterior_burnin), 5000, replace = TRUE)
posterior_samples <- posterior_burnin |>
  as_tibble() |>
  slice(draws) |>
  mutate(np = config$np)

write_rds(posterior_samples, file.path(dest_posterior, "posteriorSamples.rds"))

# function to create traceplots from thinned posterior
trace_plot <- function(post, nodes_2_plot, thin = 5000) {
  df <- post |>
    select(chain, all_of(nodes_2_plot)) |>
    group_by(chain) |>
    mutate(iteration = 1:n()) |>
    ungroup()

  posterior_mat <- df |>
    select(-chain, -iteration) |>
    as.matrix()

  print(apply(posterior_mat, 2, quantile, c(0.025, 0.5, 0.975)))

  total_iterations <- max(df$iteration)
  thin_interval <- floor(seq(1, total_iterations, length.out = thin))

  gg <- df |>
    filter(iteration %in% thin_interval) |>
    pivot_longer(cols = -c(iteration, chain), names_to = "node") |>
    mutate(chain = as.character(chain)) |>
    ggplot() +
    aes(x = iteration, y = value, color = chain) +
    geom_line() +
    facet_wrap(~node, scales = "free_y") +
    labs(x = "Iteration", y = "Value") +
    theme_bw()
  return(gg)
}

nodes <- setdiff(colnames(posterior), "chain")
n_plots_per_page <- 4 # want to put 4 nodes on a single plot
idx <- rep(
  seq(1, ceiling(length(nodes) / n_plots_per_page), by = 1),
  each = n_plots_per_page
)[1:length(nodes)]

plots <- tibble(
  nodes = nodes,
  idx = idx
)

message("Creating traceplots...")
# options(bitmapType = 'Xlib')
pb <- txtProgressBar(min = 1, max = max(plots$idx), style = 1)
for (i in seq_along(unique(plots$idx))) {
  n2p <- plots |>
    filter(idx == i) |>
    pull(nodes)

  gg <- trace_plot(posterior, n2p)

  filename <- file.path(
    dest_posterior,
    paste0("mcmcTimeseries_", sprintf("%03d", i), ".pdf")
  )
  ggsave(filename, gg)

  setTxtProgressBar(pb, i)
}
close(pb)

message("Parameters Done")

states_mcmc_list <- mcmc_list$states

states_burnin <- window(states_mcmc_list, start = burnin)

posterior_burnin2 <- states_burnin |>
  as.matrix()

state_samples <- posterior_burnin2 |>
  as_tibble() |>
  slice(draws) |>
  mutate(np = config$np)

write_rds(state_samples, file.path(dest_posterior, "stateSamples.rds"))

data_path <- file.path(dest_posterior, "modelData.rds")
model_data <- read_rds(data_path)

density_stats <- function(mcmc_list, data) {
  mcmc_matrix <- as.matrix(mcmc_list)
  draws <- sample.int(nrow(mcmc_matrix), 5000, replace = TRUE)
  abundance_sample <- mcmc_matrix[draws, ] |>
    as_tibble() |>
    pivot_longer(cols = everything(), names_to = "node", values_to = "value") |>
    mutate(n_id = as.numeric(stringr::str_extract(node, "(?<=\\[)\\d*(?=\\])")))

  all_pp <- create_all_primary_periods(data) |>
    select(-timestep)

  property_info <- data |>
    select(agrp_prp_id, property, primary_period, property_area_km2) |>
    left_join(all_pp) |>
    distinct()

  property_match <- left_join(abundance_sample, property_info) |>
    mutate(density = value / property_area_km2) |>
    group_by(
      node,
      n_id,
      agrp_prp_id,
      property,
      primary_period,
      property_area_km2
    ) |>
    summarise(
      mean = mean(density),
      variance = var(density),
      `0.025` = quantile(density, 0.025),
      `0.05` = quantile(density, 0.05),
      `0.1` = quantile(density, 0.1),
      `0.25` = quantile(density, 0.25),
      `0.5` = quantile(density, 0.5),
      `0.75` = quantile(density, 0.75),
      `0.9` = quantile(density, 0.9),
      `0.95` = quantile(density, 0.95),
      `0.975` = quantile(density, 0.975)
    ) |>
    ungroup()
  return(property_match)
}

density <- density_stats(states_mcmc_list, model_data)
write_rds(density, file.path(dest_posterior, "densitySummaries.rds"))

message("Density Done")
