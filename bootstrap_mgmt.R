library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(lubridate)

set.seed(803)

plot_dir <- "plots"

source("R/functions_data.R")

run_date <- today()
cutoff_date <- ymd("2023-12-31")


if_dir <- "11_posterior"
posterior_path <- file.path(top_dir, if_dir, "posteriorSamples.rds")
params <- read_rds(posterior_path)

posterior_path <- file.path(top_dir, if_dir, "modelData.rds")
data <- read_rds(posterior_path) |>
	mutate(
		method = if_else(method == "FIREARMS", "Firearms", method),
		method = if_else(method == "FIXED WING", "Fixed wing", method),
		method = if_else(method == "HELICOPTER", "Helicopter", method),
		method = if_else(method == "SNARE", "Snare", method),
		method = if_else(method == "TRAPS", "Trap", method)
	)

posterior_path <- file.path(top_dir, if_dir, "densitySummaries.rds")
density <- read_rds(posterior_path)

n_return <- data |>
	select(propertyID, primary_period, method) |>
	distinct() |>
	pivot_wider(names_from = method, values_from = method) |>
	unite(method, -c(propertyID, primary_period), sep = ", ", na.rm = TRUE) |>
	group_by(propertyID, method) |>
	mutate(return_interval = c(0, diff(primary_period))) |>
	ungroup() |>
	rename(methods_used = method)

model_data <- data |>
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
	mutate(
		take_density = total_take / property_area_km2,
		farm_bill = if_else(is.na(farm_bill), 0, farm_bill)
	) |>
	left_join(density) |>
	left_join(n_return) |>
	mutate(year = year(end_dates))

postal_codes_df <- read_csv(
	"../pigs-statistical/data/counties/statePostalCodes.csv"
)
postal_codes <- postal_codes_df |>
	mutate(st_name = toupper(State))

data_mis <- left_join(model_data, postal_codes) |>
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
		Postal,
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

n_ens <- 1000

dem <- params[sample.int(nrow(params), n_ens, replace = TRUE), ]
phi <- dem$phi_mu
psi <- dem$psi_phi

zeta <- 28 * exp(dem$log_nu) / 365

all_out <- tibble()
counter <- 1
n_sim <- 10000
while (counter < n_sim) {
	i <- sample.int(nrow(data_mis), 1)

	tmp_start <- data_mis |>
		slice(i) |>
		select(
			propertyID,
			property_area_km2,
			abundance_estimate,
			end_dates,
			year,
			st_name,
			property_area_km2
		)

	start_date <- tmp_start$end_dates
	property <- tmp_start$propertyID
	area <- tmp_start$property_area_km2

	date_seq <- seq.Date(start_date, to = start_date + 365, by = "4 weeks")
	N1 <- tmp_start$abundance_estimate

	Kd <- 30 # carrying capacity as a density
	K <- round(Kd * area) # carrying capacity as abundance

	if (N1 == 0 | N1 > K) {
		next
	}

	tmp_compare <- data_mis |>
		filter(propertyID %in% property, end_dates %in% date_seq) |>
		select(
			propertyID,
			property_area_km2,
			abundance_estimate,
			end_dates,
			year,
			st_name
		)

	timestep_df <- tibble(
		start_date = start_date,
		end_dates = date_seq
	) |>
		filter(end_dates <= max(tmp_compare$end_dates)) |>
		mutate(
			delta_time = as.numeric(end_dates - start_date),
			year = year(end_dates)
		)

	n_time <- nrow(timestep_df)

	if (n_time == 1) {
		next
	}

	N <- matrix(NA, n_ens, n_time)
	N[, 1] <- pmin(K, rpois(n_ens, N1))

	for (j in 2:n_time) {
		# phi_t <- rbeta(n_ens, a_phi, b_phi)
		r <- phi + zeta / 2
		nm1 <- N[, j - 1]

		lambda <- nm1 + (r - 1) * nm1 * (1 - nm1 / K)
		lambda <- pmax(0, lambda)

		if (any(is.na(lambda))) {
			print(i)
			stop()
		}

		N[, j] <- rpois(n_ens, lambda)
	}

	timestep_df$abundanceNoMgmt <- apply(N, 2, median)

	out <- left_join(tmp_compare, timestep_df, by = c("end_dates", "year")) |>
		filter(delta_time > 0) |>
		mutate(
			density_estimate = abundance_estimate / area,
			densityNoMgmt_estimate = abundanceNoMgmt / area,
			area = area,
			K = K,
			D1 = N1 / area,
			density_cat = if_else(D1 < 2, "Low", "Medium"),
			density_cat = if_else(D1 > 6, "High", density_cat),
			density_cat = factor(density_cat, levels = c("Low", "Medium", "High"))
		)
	all_out <- bind_rows(all_out, out)

	counter <- counter + 1

	if (counter %% 1000 == 0) message(counter / n_sim * 100, "%")
}

glimpse(all_out)

keep <- seq(28, 28 * 30, by = 28 * 4)
weeks_f <- sort(keep / 7)

with_mgmt <- all_out |>
	select(delta_time, st_name, density_estimate, D1, density_cat) |>
	mutate(
		mgmt = "Yes",
		delta_density = density_estimate - D1,
		percent_change = delta_density / D1 * 100
	)

without_mgmt <- all_out |>
	select(delta_time, st_name, densityNoMgmt_estimate, D1, density_cat) |>
	rename(density_estimate = densityNoMgmt_estimate) |>
	mutate(
		mgmt = "No",
		delta_density = density_estimate - D1,
		percent_change = delta_density / D1 * 100
	)

df <- bind_rows(with_mgmt, without_mgmt)

good_enough <- df |>
	group_by(st_name, delta_time, density_cat) |>
	count() |>
	ungroup() |>
	filter(n >= 10) |>
	select(-n)

all_3_cats <- good_enough |>
	select(st_name, density_cat) |>
	distinct() |>
	count(st_name) |>
	filter(n == 3) |>
	pull(st_name)

df_plot <- left_join(good_enough, df) |>
	filter(st_name %in% all_3_cats)

df_plot |>
	filter(delta_time %in% keep, !st_name %in% c("NORTH CAROLINA", "VIRGINIA")) |>
	mutate(weeks = delta_time / 7, weeks = factor(weeks, levels = weeks_f)) |>
	ggplot() +
	aes(
		x = weeks,
		y = delta_density,
		fill = mgmt
	) +
	geom_boxplot(outliers = FALSE) +
	geom_hline(yintercept = 0, linetype = "dashed") +
	facet_grid(density_cat ~ st_name, scales = "free_y") +
	labs(x = "Weeks", y = "Change in density", fill = "Management") +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave(
	"plots/effectOfMgmt/bootstapAll.jpeg",
	units = "cm",
	width = 24,
	height = 18
)
