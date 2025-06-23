library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(knitr)
library(lubridate)

source("R/functions_misc.R")

# project directories
plot_dir <- "plots"
data_dir <- "data"

# shared data across projects live here
data_store <- "../data-store"

source("R/functions_data.R")

cutoff_date <- ymd("2023-12-31")

# last batch iteration was number 11
# from this batch we need:
# - model data
# - density summaries
if_dir <- "11_posterior"
posterior_path_dir <- file.path(data_dir, if_dir)

fname <- "modelData.rds"
model_data <- read_rds(file.path(posterior_path_dir, fname))
data <- fix_method_names(model_data)

fname <- "densitySummaries.rds"
density <- read_rds(file.path(posterior_path_dir, fname))

# will need to summarize by primary period, method used, and property
make_methods_used <- function(df) {
	df |>
		select(propertyID, primary_period, method) |>
		distinct() |>
		pivot_wider(names_from = method, values_from = method) |>
		unite(method, -c(propertyID, primary_period), sep = ", ", na.rm = TRUE) |>
		rename(methods_used = method)
}

methods_used <- make_methods_used(data)

# for analysis by primary period, we need to calculate totals
# by property, primary period, and method for take and effort
# then join to density summaries and methods_used
make_model_data <- function(df) {
	df |>
		group_by(
			propertyID,
			agrp_prp_id,
			start_dates,
			end_dates,
			st_name,
			cnty_name,
			farm_bill,
			alws_agrprop_id,
			property,
			primary_period,
			property_area_km2,
			county_code
		) |>
		summarise(
			total_take = sum(take),
			total_effort_per = sum(effort_per),
			n_events = n()
		) |>
		ungroup() |>
		mutate(take_density = total_take / property_area_km2) |>
		mutate(farm_bill = if_else(is.na(farm_bill), 0, farm_bill))
}

model_data <- make_model_data(data)

# last bit of data prep
data_prep <- function(df, df_density, df_methods, cutoff_date) {
	df |>
		left_join(df_density) |>
		left_join(df_methods) |>
		mutate(year = year(end_dates)) |>
		filter(end_dates <= cutoff_date) |>
		mutate(abundance_estimate = round(`0.5` * property_area_km2)) |>
		rename(
			density_estimate = `0.5`,
			density_quantile_0.05 = `0.05`,
			density_quantile_0.95 = `0.95`
		) |>
		select(
			propertyID,
			end_dates,
			year,
			st_name,
			cnty_name,
			county_code,
			property_area_km2,
			total_take,
			take_density,
			density_estimate,
			abundance_estimate,
			n_events,
			methods_used
		)
}

data_mis <- data_prep(model_data, density, methods_used, cutoff_date)

# need to create a data frame with all primary periods for each property
# then join to data_mis to get pp for each end_dates
make_all_timesteps <- function(df) {
	prop_vec <- df |>
		pull(propertyID) |>
		unique()

	all_timesteps <- tibble()
	for (i in seq_along(prop_vec)) {
		tmp <- df |>
			filter(propertyID == prop_vec[i])

		min_time <- min(tmp$end_dates)
		max_time <- max(tmp$end_dates)

		prop_info <- tmp |>
			dplyr::slice(1) |>
			select(propertyID, st_name, county_code, property_area_km2)

		tmpp <- tibble(
			propertyID = prop_vec[i],
			end_dates = seq.Date(min_time, max_time, 28)
		) |>
			left_join(prop_info, by = "propertyID")

		all_timesteps <- bind_rows(all_timesteps, tmpp)
	}
	all_timesteps |>
		mutate(pp = as.numeric(as.factor(end_dates)))
}

all_pp <- make_all_timesteps(data_mis)

# only keep properties with more than one estimate
p2 <- data_mis |>
	count(propertyID) |>
	filter(n > 1) |>
	pull(propertyID)

# calculate realized growth rate between primary periods
# need to get abundance from previous primary period
# and the number of primary periods between estimates
# and the catch from the previous primary period
create_data_growth <- function(df, df_pp, p2) {
	left_join(df, df_pp) |>
		filter(propertyID %in% p2) |>
		select(
			propertyID,
			end_dates,
			pp,
			st_name,
			total_take,
			abundance_estimate
		) |>
		arrange(propertyID, pp) |>
		group_by(propertyID) |>
		mutate(
			take_m1 = c(NA, total_take[1:(n() - 1)]),
			N_m1 = c(NA, abundance_estimate[1:(n() - 1)]),
			delta_pp = c(NA, diff(pp))
		) |>
		ungroup() |>
		filter(!is.na(take_m1))
}

data_for_growth <- create_data_growth(data_mis, all_pp, p2)

# lambda = N[t] / (N[t-1] - Catch[t-1])
calculate_lambda <- function(df) {
	df |>
		mutate(
			Z = N_m1 - take_m1,
			lambda_star = abundance_estimate / Z,
			lambda = lambda_star^(1 / delta_pp)
		) |>
		filter(!is.nan(lambda))
}

realized_growth <- calculate_lambda(data_for_growth)

# look at proportion at population taken within the last x timeframe
# add astricks where sample size is small

glimpse(realized_growth)

## the mean population growth rate for each property ----
## separated by state
mean_growth <- function(df) {
	df |>
		group_by(propertyID, st_name) |>
		summarise(
			n = n(),
			mu = mean(lambda)
		) |>
		ungroup()
}

mu_lambda <- mean_growth(realized_growth)

plot_mean_lambdas <- function(df) {
	y_order <- mu_lambda |>
		group_by(st_name) |>
		summarise(y = median(mu), n = n()) |>
		arrange(y) |>
		mutate(st_name2 = paste0(st_name, " (n=", n, ")")) |>
		select(-n)

	y_fac <- y_order$st_name2

	mu_lambda |>
		left_join(y_order) |>
		mutate(st_name2 = factor(st_name2, levels = y_fac)) |>
		ggplot() +
		aes(y = mu, x = st_name2) +
		geom_jitter(
			size = 1,
			# alpha = 0.2,
			position = position_jitter(width = 0.25),
			color = "grey",
			fill = "grey"
		) +
		geom_boxplot(outliers = FALSE, alpha = 0) +
		geom_hline(yintercept = 1, linetype = "dashed") +
		coord_flip() +
		labs(
			title = "Rate of change (with removals)",
			y = "Population growth rate / 28-days",
			x = ""
		) +
		theme_bw()
}


gname <- rateOfChange.jpeg
ggsave(
	gname,
	path = plot_dir,
	units = "cm",
	width = 18,
	height = 12
)


realized_growth |>
	filter(delta_pp == 13) |>
	mutate(x = "x") |>
	ggplot(aes(x = x, y = lambda_star)) +
	geom_jitter(size = 1, color = "grey", fill = "grey") +
	geom_boxplot(outliers = FALSE, alpha = 0) +
	labs(
		x = "Annual population growth rate",
		y = ""
	) +
	theme_bw() +
	theme(axis.text.x = element_blank())

ggsave(
	"plots/effectOfMgmt/rateOfChangeAnnual.jpeg",
	units = "cm",
	width = 6,
	height = 12
)

realized_growth |>
	filter(lambda >= 1) |>
	mutate(
		`Maintenance` = -1 * ((1 - lambda) / lambda),
		`50% reduction` = -1 * ((0.5 - lambda) / lambda),
		`Elimination` = -1 * ((0 - lambda) / lambda)
	) |>
	pivot_longer(
		cols = c(`Elimination`, `50% reduction`, `Maintenance`),
		names_to = "Management goal",
		values_to = "y"
	) |>
	mutate(
		`Management goal` = factor(
			`Management goal`,
			levels = c("Elimination", "50% reduction", "Maintenance")
		)
	) |>
	ggplot() +
	aes(x = lambda, y = y, color = `Management goal`) +
	geom_line(linewidth = 2) +
	labs(
		title = "Theoretical expectation",
		x = "Population growth rate / 28-days",
		y = "Proportion needed to reduce population"
	) +
	# coord_cartesian(ylim = c(0, 1)) +
	theme_bw()

ggsave(
	"plots/effectOfMgmt/mgmtObjectiveTheory.jpeg",
	units = "cm",
	width = 12,
	height = 12
)

# population growth
# need to get Z for each primary period
# need the number of PPs between each estimate
# lambda <- N[t+1] / (N[t] - C[t])

lambda <- realized_growth$lambda
round(quantile(lambda, c(0.05, 0.5, 0.95)), 3)

lambda_year <- lambda^13
round(quantile(lambda_year, c(0.05, 0.5, 0.95)), 3)

xpop <- (1 - (1 / lambda[lambda >= 1])) * 100
round(quantile(xpop, c(0.05, 0.5, 0.95)), 3)

xpop_year <- (1 - (1 / lambda_year[lambda_year >= 1])) * 100
round(quantile(xpop_year, c(0.05, 0.5, 0.95)), 3)


realized_growth |>
	count(st_name)


percent_take <- data_mis |>
	filter(abundance_estimate > 0) |>
	mutate(
		percent_take = total_take / abundance_estimate * 100,
		methods_used = stringr::str_replace(
			methods_used,
			"Firearms",
			"Sharpshooting"
		)
	) |>
	group_by(methods_used, st_name) |>
	summarise(
		n = n(),
		low = quantile(percent_take, 0.05),
		q1 = quantile(percent_take, 0.25),
		med = quantile(percent_take, 0.5),
		q3 = quantile(percent_take, 0.75),
		high = quantile(percent_take, 0.95)
	) |>
	ungroup() # |>
# filter(n >= 10)

percent_take <- left_join(data, mu_lambda) |>
	filter(mu >= 1) |>
	group_by(
		propertyID,
		agrp_prp_id,
		start_dates,
		end_dates,
		st_name,
		cnty_name,
		farm_bill,
		alws_agrprop_id,
		property,
		primary_period,
		property_area_km2,
		county_code,
		method
	) |>
	summarise(
		total_take = sum(take),
		total_effort_per = sum(effort_per),
		total_trap_count = sum(trap_count),
		n_events = n()
	) |>
	left_join(density) |>
	ungroup() |>
	mutate(abundance_estimate = round(`0.5` * property_area_km2)) |>
	filter(abundance_estimate > 0) |>
	mutate(
		mu = total_take / abundance_estimate,
		method = stringr::str_replace(method, "Firearms", "Sharpshooting")
	) |>
	select(
		method,
		st_name,
		propertyID,
		end_dates,
		total_take,
		abundance_estimate,
		mu
	) |>

	group_by(method, propertyID, st_name) |>
	summarise(
		n = n(),
		mu = mean(mu)
	) |>
	ungroup() |>

	mutate(group = "Actual")

# what they should be doing
mu_x <- realized_growth |>
	filter(lambda >= 1) |>
	mutate(
		mu_maintain = -1 * ((1 - lambda) / lambda),
		mu_reduction_50 = -1 * ((0.5 - lambda) / lambda)
	) |>

	group_by(propertyID, st_name) |>
	summarise(
		n = n(),
		mu_maintain = mean(mu_maintain),
		mu_reduction_50 = mean(mu_reduction_50),
	) |>
	ungroup() |>

	select(mu_maintain, mu_reduction_50, propertyID, st_name) |>
	pivot_longer(
		cols = c(mu_maintain, mu_reduction_50),
		names_to = "method",
		values_to = "mu"
	) |>
	mutate(group = "Theoretical")

bind_rows(mu_x, percent_take) |>
	filter(!st_name %in% c("NORTH CAROLINA", "VIRGINIA", "WEST VIRGINIA")) |>
	mutate(
		method = if_else(grepl("maintain", method), "Maintenance", method),
		method = if_else(grepl("reduction", method), "50% reduction", method),
		method = factor(
			method,
			levels = c(
				"50% reduction",
				"Maintenance",
				"Fixed wing",
				"Helicopter",
				"Sharpshooting",
				"Snare",
				"Trap"
			)
		)
	) |>
	ggplot() +
	aes(x = mu, y = method, color = group) +
	geom_boxplot() +
	facet_wrap(~st_name) +
	coord_cartesian(xlim = c(0, 1)) +
	labs(
		x = "Proportion of population removed in a 28-day period",
		y = "Method"
	) +
	theme_bw()

ggsave(
	"plots/effectOfMgmt/proportionRemovedMethodState_growingPops_property.jpeg",
	units = "cm",
	width = 24,
	height = 18
)

facs <- percent_take |>
	select(methods_used, med) |>
	arrange(med) |>
	pull(methods_used)

effort <- data |>
	select(propertyID, end_dates, method, st_name) |>
	distinct() |>
	mutate(effort = 1)

property_effort <- left_join(all_timesteps, effort) |>
	mutate(
		method = if_else(is.na(method), "No effort", method),
		effort = if_else(is.na(effort), 0, effort)
	)


property_effort2 <- property_effort |>
	group_by(propertyID, st_name) |>
	summarise(n_pp_with_effort = sum(effort), n_pp = n())

# make proportion of PPs figure for farm bill only and all mis

percent_effort <- property_effort |>
	mutate(effort = 1) |>
	group_by(method, propertyID) |>
	summarise(n_method = sum(effort)) |>
	left_join(property_effort2) |>
	ungroup() |>
	mutate(
		percent_method_used = n_method / n_pp,
		method = stringr::str_replace(method, "Firearms", "Sharpshooting")
	) |>
	group_by(method, st_name) |>
	summarise(
		n = n(),
		low = quantile(percent_method_used, 0.05),
		q1 = quantile(percent_method_used, 0.25),
		med = quantile(percent_method_used, 0.5),
		q3 = quantile(percent_method_used, 0.75),
		high = quantile(percent_method_used, 0.95)
	) |>
	ungroup() |>
	# filter(n >= 10) |>
	mutate(
		method = factor(
			method,
			levels = c(
				"No effort",
				"Fixed wing",
				"Helicopter",
				"Sharpshooting",
				"Snare",
				"Trap"
			)
		)
	)

percent_effort |>
	ggplot() +
	aes(x = med, y = method) +
	geom_linerange(aes(xmin = q1, xmax = q3), linewidth = 2) +
	geom_linerange(aes(xmin = low, xmax = high), linewidth = 0.5) +
	geom_point(size = 3, color = "black") +
	geom_point(size = 2, color = "white") +
	facet_wrap(~st_name) +
	coord_cartesian(xlim = c(0, 1)) +
	labs(
		x = "Proportion of PPs with effort",
		y = "Method"
	) +
	theme_bw()

ggsave(
	"plots/effectOfMgmt/effortState.jpeg",
	units = "cm",
	width = 24,
	height = 10
)
