###############################################################################
# Fit a state-specific removal model
#
# Purpose:
# Fit the model for one state and one ecoregion. The script prepares the input
# data, sets model options, runs MCMC, and saves diagnostic results.
#
# Run requirements:
# Run this script from the project directory in the configured HPC environment.
# The environment must provide the `prod` config profile and the helper
# functions used below. This script does not source all of those functions.
#
# Input files:
# - data_store/MIS/2026-03-25/processed/
#   dev_MIS.Effort.Take.all_methods.Daily.Events.csv
# - data_store/county_ecoregion.csv
# - data/original_ids.rds
#
# Environment variables:
# - data_store: path to the data store.
# - STATENAME: uppercase state name. The default is OKLAHOMA.
# - SLURM_ARRAY_TASK_ID: ecoregion index from 1 to 5. The default is 1.
# - SLURM_CPUS_PER_TASK: number of MCMC chains.
# - pbStyle: progress-bar style used by boaR.
#
# Outputs:
# The script writes the method lookup table, MCMC output, and diagnostics under
# out/states2/states2-2026-03-25/<STATE>_<ECOREGION>/.
###############################################################################

library(tidyr)
library(dplyr)
library(readr)
library(nimble)
library(coda)
library(lubridate)
library(boaR)
library(parallel)

source("R/state_inits.R")

set_boaR_options(pbStyle = as.numeric(Sys.getenv("pbStyle")))
data_store <- Sys.getenv("data_store")

# Set the project name and output directory.
project <- "states2"

write_dir <- file.path("out", project)

n_chains <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
message("n_chains: ", n_chains)

pull_date <- "2026-03-25"
post_round <- "first"

project_pull <- paste0(project, "-", pull_date)

# Read MCMC settings from the external production config profile.
config_name <- "prod"
config <- config::get(config = config_name)
n_iter <- config$n_iter

# Read the primary-period length and data preparation option.
interval <- config$interval

create_new <- config$create_new

# Read the processed MIS event data.
mis <- "MIS"
mis_processed <- "processed"
file_name <- "dev_MIS.Effort.Take.all_methods.Daily.Events.csv"

fname <- file.path(data_store, mis, pull_date, mis_processed, file_name)
df <- readr::read_csv(fname, show_col_types = FALSE)

data_mis <- get_data(df, interval, create_new)

original_ids <- read_rds("data/original_ids.rds")

# Keep records with complete landscape covariates. Exclude original-fit properties.
data_complete <- data_mis |>
	filter(
		!is.na(c_road_den),
		!is.na(c_rugged),
		!is.na(c_canopy),
		!propertyID %in% original_ids
	)

# write_rds(
#   data_complete,
#   file.path(data_store, "stateFit", project_pull, "model_data.rds")
# )

# Select the state from the job environment.
st <- Sys.getenv("STATENAME")
st <- if_else(st == "", "OKLAHOMA", st) # for testing
message("\n")
message("STATENAME: ", st)

county_ecoregion <- readr::read_csv(
	file.path(data_store, "county_ecoregion.csv")
)

# Add ecoregion labels and reset primary-period indexes for this state.
data_for_nimble <- data_complete |>
	filter(st_name == st) |>
	select(-p) |>
	mutate(primary_period = primary_period - min(primary_period) + 1) |>
	left_join(county_ecoregion, by = "county_code")


ecoregion_vec <- c(
	"SOUTH CENTRAL SEMIARID PRAIRIES",
	"TEXAS-LOUISIANA COASTAL PLAIN",
	"TAMAULIPAS-TEXAS SEMIARID PLAIN",
	"WARM DESERTS",
	"SOUTHEASTERN USA PLAINS"
)

array_num <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
eco_sub <- ecoregion_vec[array_num]

message("\nSelected ecoregion: ", eco_sub, "\n")

data_for_nimble <- data_for_nimble |>
	mutate(
		county_code = sprintf("%05d", as.integer(county_code))
	) |>
	filter(ecoregion_name == eco_sub)

eco_path <- gsub("\\s+", "_", trimws(eco_sub))
st_path <- paste0(st, "_", eco_path)
path <- file.path(write_dir, project_pull, st_path)

if (!dir.exists(path)) {
	dir.create(path, showWarnings = FALSE, recursive = TRUE)
}

n <- nrow(data_for_nimble)
np <- length(unique(data_for_nimble$propertyID))
nm <- length(unique(data_for_nimble$method))
method_names <- data_for_nimble$method
methods <- unique(method_names)

message("n events: ", n)
message("n properties: ", np)
message("n methods: ", nm, " (", paste(methods, collapse = ", "), ")")

constants <- nimble_constants(
	df = data_for_nimble,
	interval = 28,
	post_round = "first"
)
data <- nimble_data(data_for_nimble)

# Create a method lookup table for the methods in this state and ecoregion.
method_lookup_table <- tibble(
	method_names = method_names,
	m_vec = boaR:::method_factors(data_for_nimble),
	ts_id = constants$ts_id
) |>
	mutate(method_idx = as.numeric(as.factor(m_vec))) |>
	distinct() |>
	arrange(method_idx)

write_rds(method_lookup_table, file.path(path, "method_lookup_table.rds"))

# Set model options and the MCMC output directory.
model_flags <- get_model_flags(data_for_nimble)
params_check <- config$params_check
dest <- file.path(path, "mcmc")

# Set initial values for the selected state.
init_list <- state_inits(st)
# Use NULL initial values for ecoregion array task 3.
if (array_num == 3) {
	init_list <- NULL
}

# Set monitored abundance nodes and the MCMC buffer size.
if (st == "TEXAS" || st == "OKLAHOMA") {
	nH_p <- constants$nH_p
	monitors_add <- c("N", paste0("N[", nH_p, "]"))
	buffer <- 5000
} else {
	# otherwise, monitor all primary periods
	monitors_add <- "N"
	buffer <- 600
}

if (st == "GEORGIA") {
	buffer <- 800
}

beta_p_samplers <- list()

if (model_flags$use_beta_p) {
	beta_p_samplers <- lapply(seq_len(constants$n_betaP), function(k) {
		list(
			node = sprintf(
				"beta_p[%d, %d]",
				constants$beta_p_row[k],
				constants$beta_p_col[k]
			),
			type = "ess",
			control = list(
				maxContractions = 1000,
				maxContractionsWarning = TRUE
			)
		)
	})

	message(
		"Adding individual ESS samplers for ",
		length(beta_p_samplers),
		" beta_p nodes"
	)
}

area_block_samplers <- list()

trap_snare_methods <- method_lookup_table |>
	filter(ts_id > 0) |>
	distinct(method_idx, ts_id) |>
	arrange(method_idx)

if (nrow(trap_snare_methods) > 0) {
	for (k in seq_len(nrow(trap_snare_methods))) {
		method_idx <- trap_snare_methods$method_idx[k]
		ts_idx <- trap_snare_methods$ts_id[k]

		if (model_flags$use_traps_and_snares) {
			gamma_node <- sprintf("log_gamma[%d]", ts_idx)
			p_mu_node <- sprintf("p_mu[%d]", ts_idx)

			gamma_var <- 1 / constants$log_gamma_tau[ts_idx]
			p_mu_var <- 1 / constants$p_mu_tau[ts_idx]
		} else {
			# Only one of traps or snares occurs, so these are scalar nodes.
			gamma_node <- "log_gamma"
			p_mu_node <- "p_mu"

			gamma_var <- 1 / constants$log_gamma_tau[1]
			p_mu_var <- 1 / constants$p_mu_tau[1]
		}

		nodes <- c(
			sprintf("log_rho[%d]", method_idx),
			gamma_node,
			p_mu_node
		)

		prior_var <- c(
			1 / constants$log_rho_tau[method_idx],
			gamma_var,
			p_mu_var
		)

		area_block_samplers[[length(area_block_samplers) + 1L]] <- list(
			node = nodes,
			type = "RW_block",
			control = list(
				adaptive = TRUE,
				adaptScaleOnly = FALSE,
				adaptInterval = 100,
				adaptFactorExponent = 0.4,
				scale = 0.1,
				propCov = diag(prior_var),
				tries = 2
			)
		)

		message(
			"Adding trap/snare area block: ",
			paste(nodes, collapse = ", ")
		)
	}
}

custom_samplers <- c(beta_p_samplers, area_block_samplers)

# Run Markov chain Monte Carlo (MCMC) and save the samples in chunks.
mcmc_parallel(
	n_chains = n_chains,
	model_constants = constants,
	model_data = data,
	model_flags = model_flags,
	params_check = params_check,
	n_iters = n_iter,
	dest = dest,
	monitors_add = monitors_add,
	custom_samplers = custom_samplers,
	export = "calc_log_area",
	buffer = buffer,
	beta1 = init_list$beta1,
	beta_p = init_list$beta_p,
	p_mu = init_list$p_mu,
	log_gamma = init_list$log_gamma,
	log_rho = init_list$log_rho,
	psi_phi = init_list$psi_phi,
	phi_mu = init_list$phi_mu,
	log_nu = init_list$log_nu
)

# Save MCMC diagnostics.
analysis_dir <- file.path(path, "analysis")
if (!dir.exists(analysis_dir)) {
	dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)
}

mcmc_diagnostics(
	mcmc_dir = dest,
	dest = analysis_dir,
	data = data_for_nimble,
	params_check = params_check
)
