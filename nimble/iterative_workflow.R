#---------
#
# Workflow for fitting property-level Bayes model to MIS data
#
#---------

library(dplyr)
library(tidyr)
library(readr)
library(parallel)

# config_name <- "hpc_dev"
config_name <- "default"
config <- config::get(config = config_name)

source("nimble/functions_data.R")
source("nimble/prep_and_run_mcmc.R")
source("nimble/iterative_fitting.R")

# ===================================================
# Data ingest ----
# ===================================================

data_repo <- "../data-store"

## MIS data ----
file <- file.path(data_repo, config$file_mis)
interval <- config$interval
dev <- config$dev

data_mis <- get_data(file, interval, data_repo)

## observation covariates ----
file <- file.path(data_repo, config$file_land)
data_obs <- get_obs_covars(file)

## join MIS with observation covariates ----
data_join <- left_join(data_mis, data_obs, by = join_by(county_code))

## filter missing states ----
data_join2 <- data_join |>
	filter(!st_name %in% c("CALIFORNIA", "ALABAMA", "ARIZONA", "ARKANSAS"))

targets::tar_assert_true(!any(is.na(data_join2$c_road_den)))
targets::tar_assert_true(!any(is.na(data_join2$c_rugged)))
targets::tar_assert_true(!any(is.na(data_join2$c_canopy)))

## join with farm bill properties ----
data_farm_bill <- read_csv(file.path(
	data_repo,
	"All_FB_Agreements_long_2024-05-30.csv"
))
farm_bill_properties <- data_farm_bill |>
	rename(alws_agrprop_id = propertyID) |>
	select(-agreement_name, -property_name) |>
	mutate(farm_bill = 1)

data_final <- left_join(data_join2, farm_bill_properties) |>
	mutate(
		property = as.numeric(as.factor(propertyID)),
		county = as.numeric(as.factor(county_code))
	)

data_final |>
	filter(farm_bill == 1) |>
	pull(propertyID) |>
	unique() |>
	length()

print_info <- function(df) {
	fit_properties <- length(unique(df$propertyID))
	nfb <- df |>
		filter(farm_bill == 1) |>
		pull(propertyID) |>
		unique() |>
		length()
	message("=======================================")
	message("Total properties in data: ", fit_properties)
	message("Total Farm Bill properties in data: ", nfb)
	message("Total counties in data: ", length(unique(df$county_code)))
	message("=======================================")
}

print_info(data_final)

dataIDs <- data_final |>
	mutate(
		ppID = paste0(propertyID, "-", primary_period),
		method_fac = as.numeric(as.factor(method)),
		eventID = paste0(ppID, "-", order, "-", method_fac)
	)

write_csv(dataIDs, "data/originalFitDataWithIDs.csv")

params_check <- config$params_check
out_dir <- config$out_dir
files_in_out_dir <- list.files(out_dir)
first_fit <- length(files_in_out_dir) == 0

monitors_add <- "N"
custom_samplers <- NULL

if (first_fit) {
	# run first fit

	informed <- FALSE

	data_for_nimble <- subset_data_for_development(data_final) |>
		mutate(primary_period = primary_period - min(primary_period) + 1)

	print_info(data_for_nimble)
	data_for_nimble |>
		group_by(method) |>
		summarise(n = n(), take = sum(take)) |>
		print()

	first_fit_properties <- unique(data_for_nimble$propertyID)
	all_properties <- unique(data_final$propertyID)
	not_fit_properties <- setdiff(all_properties, first_fit_properties)

	np <- length(first_fit_properties)
	dest_mcmc <- file.path(out_dir, paste0(1, "_mcmc"))
	dest_posterior <- file.path(out_dir, paste0(1, "_posterior"))

	fit_successfully <- tibble(
		property = first_fit_properties,
		round = 1,
		fit = TRUE,
		n_properties = np,
		dir = dest_posterior
	)

	last_fit_properties <- unique(data_for_nimble$propertyID)
	all_properties <- unique(data_final$propertyID)
	n_total_properties <- length(all_properties)
	n_properties_to_fit <- n_total_properties -
		length(unique(data_for_nimble$propertyID))

	to_fit_properties <- setdiff(all_properties, last_fit_properties)

	sample_size <- np

	set.seed(9)
	property_order <- tibble(
		property = to_fit_properties,
		round = sample.int(
			ceiling(n_properties_to_fit / sample_size) + 1,
			length(to_fit_properties),
			replace = TRUE
		) +
			1,
		fit = NA,
		n_properties = NA,
		dir = NA
	)

	fit_successfully <- bind_rows(fit_successfully, property_order)

	write_rds(fit_successfully, file.path(data_repo, "iterativeFitting.rds"))

	finished <- prep_and_run_mcmc(
		informed,
		NULL,
		data_repo,
		dest_mcmc,
		dest_posterior,
		data_for_nimble,
		monitors_add,
		custom_samplers
	)

	source("nimble/check_mcmc.R")
} else {
	# run iterative fitting

	informed <- TRUE

	fit_successfully <- read_rds(file.path(data_repo, "iterativeFitting.rds"))
	good_fits <- fit_successfully |> filter(fit)

	start <- good_fits |>
		pull(round) |>
		max() +
		1

	end <- max(fit_successfully$round)

	for (i in start:end) {
		# need to add a check so that if a property failed, we don't try to read the posterior from that fit
		# need to keep track of posterior file paths and order of properties
		# look back to last fit group and use
		fit_successfully <- read_rds(file.path(data_repo, "iterativeFitting.rds"))
		good_fits <- fit_successfully |> filter(fit)

		post_path <- good_fits |>
			filter(round == max(round)) |>
			pull(dir) |>
			unique()
		last_good_fit <- basename(post_path)

		data_last_fit <- read_rds(file.path(
			out_dir,
			last_good_fit,
			"modelData.rds"
		))
		old_properies <- unique(data_last_fit$propertyID)

		current_properties <- fit_successfully |>
			filter(round == i) |>
			pull(property)

		data_for_nimble <- data_final |>
			filter(propertyID %in% c(old_properies, current_properties)) |>
			mutate(primary_period = primary_period - min(primary_period) + 1)

		message("\n==========================================================")
		print_info(data_for_nimble)
		message("==========================================================\n")

		np <- length(unique(data_for_nimble$propertyID))
		dest_mcmc <- file.path(out_dir, paste0(i, "_mcmc"))
		dest_posterior <- file.path(out_dir, paste0(i, "_posterior"))

		finished <- prep_and_run_mcmc(
			informed,
			post_path,
			data_repo,
			dest_mcmc,
			dest_posterior,
			data_for_nimble,
			monitors_add,
			custom_samplers,
			post_file = NULL
		)
		# post_file = file.path(post_path, "posteriorSamples.rds"))

		fit_successfully <- fit_successfully |>
			mutate(
				fit = if_else(property %in% current_properties, finished, fit),
				n_properties = if_else(
					property %in% current_properties,
					np,
					n_properties
				),
				dir = if_else(property %in% current_properties, dest_posterior, dir)
			)

		write_rds(fit_successfully, file.path(data_repo, "iterativeFitting.rds"))

		if (finished) source("nimble/check_mcmc.R")
	}
}
