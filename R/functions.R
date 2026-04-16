library(readr)
library(dplyr)
library(tidyr)


# read rds files
get_rds <- function(file) {
	read_rds(file)
}


# rename methods for better plotting
fix_method_names <- function(df) {
	df |>
		mutate(
			method = if_else(method == "FIREARMS", "Ground-shooting", method),
			method = if_else(method == "FIXED WING", "Fixed wing", method),
			method = if_else(method == "HELICOPTER", "Helicopter", method),
			method = if_else(method == "SNARE", "Snare", method),
			method = if_else(method == "TRAPS", "Trap", method)
		) |>
		left_join(state_names, by = "st_name") |>
		select(-st_name) |>
		rename(st_name = state)
}


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

# will need to summarize by primary period, method used, and property
make_methods_used <- function(df) {
	df |>
		select(propertyID, primary_period, method) |>
		distinct() |>
		pivot_wider(names_from = method, values_from = method) |>
		unite(method, -c(propertyID, primary_period), sep = ", ", na.rm = TRUE) |>
		rename(methods_used = method)
}

# prepare data for analysis and plotting
data_prep <- function(df, df_density) {
	tmp <- df |>
		fix_method_names() |>
		make_model_data()

	df_methods <- make_methods_used(df)

	tmp |>
		left_join(df_density) |>
		left_join(df_methods) |>
		mutate(year = year(end_dates)) |>
		mutate(abundance_estimate = round(`0.5` * property_area_km2)) |>
		rename(
			density_estimate = `0.5`,
			density_quantile_0.05 = `0.05`,
			density_quantile_0.95 = `0.95`
		) |>
		select(
			propertyID,
			farm_bill,
			end_dates,
			`year`,
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
		filter(!is.na(delta_pp))
}

calculate_lambda <- function(df) {
	df |>
		mutate(
			Z = N_m1 - take_m1,
			lambda_star = abundance_estimate / Z,
			lambda = lambda_star^(1 / delta_pp)
		) |>
		filter(!is.nan(lambda))
}

calculate_realized_r <- function(df) {
	all_pp <- make_all_timesteps(df)

	# only keep properties with more than one estimate
	p2 <- df |>
		count(propertyID) |>
		filter(n > 1) |>
		pull(propertyID)

	data_for_growth <- create_data_growth(df, all_pp, p2)

	calculate_lambda(data_for_growth)
}

# the mean population growth rate for each property
# separated by state
mean_growth <- function(df) {
	df |>
		group_by(propertyID, st_name) |>
		summarise(
			n = n(),
			mu = mean(lambda)
		) |>
		ungroup()
}

plot_mean_lambdas <- function(df) {
	y_order <- df |>
		group_by(st_name) |>
		summarise(y = median(mu), n = n()) |>
		arrange(y) |>
		select(-n)

	y_fac <- y_order$st_name

	g <- df |>
		left_join(y_order) |>
		mutate(st_name = factor(st_name, levels = y_fac)) |>
		ggplot() +
		aes(y = mu, x = st_name) +
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
		theme_bw() +
		theme(
			axis.text = element_text(size = 14),
			axis.title = element_text(size = 16),
			legend.title = element_text(size = 16),
			strip.text = element_text(size = 14),
			legend.text = element_text(size = 14)
		)

	ggsave(
		filename = "plots/lambda_hat_by_state.jpeg",
		plot = g,
		width = 10,
		height = 6,
		units = "in",
		dpi = 300
	)

	g
}

prop_take_by_method <- function(df, df_density) {
	df |>
		fix_method_names() |>
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
			total_effort = sum(effort),
			total_log_effort = sum(log(effort)),
			total_effort_per = sum(effort_per),
			total_log_effort_per = sum(log(effort_per)),
			total_trap_count = sum(trap_count),
			total_log_trap_count = sum(log(trap_count)),
			total_events = n()
		) |>
		left_join(df_density) |>
		ungroup() |>
		mutate(abundance_estimate = round(`0.5` * property_area_km2)) |>
		filter(abundance_estimate > 0) |>
		mutate(
			prop_take = total_take / abundance_estimate,
			density = `0.5`,
			method = stringr::str_replace(method, "Firearms", "Ground-shooting")
		) |>
		group_by(method, propertyID, st_name) |>
		summarise(
			mu_events = mean(total_events),
			mu_prop_take = mean(prop_take),
			mu_total_effort = mean(total_effort),
			mu_total_log_effort = mean(total_log_effort),
			mu_total_effort_per = mean(total_effort_per),
			mu_total_log_effort_per = mean(total_log_effort_per),
			mu_total_trap_count = mean(total_trap_count),
			mu_total_log_trap_count = mean(total_log_trap_count)
		) |>
		ungroup()
}

plot_percent_take <- function(df_actual) {
	tmp <- bind_rows(df_actual) |>
		# removing these states because their properties are all declining
		filter(!st_name %in% c("North Carolina", "Virginia", "West Virginia"))

	plot_states <- function(df) {
		ta <- df |>
			mutate(
				method = factor(
					method,
					levels = c(
						"Fixed wing",
						"Helicopter",
						"Ground-shooting",
						"Snare",
						"Trap"
					)
				)
			)

		ta |>
			ggplot() +
			aes(x = mu_prop_take, y = method) +
			geom_boxplot() +
			facet_wrap(~st_name) +
			coord_cartesian(xlim = c(0, 1)) +
			labs(
				x = "Mean proportion of individuals removed in a 28-day period",
				y = "Method",
				color = "Actual"
			) +
			scale_y_discrete(drop = FALSE) +
			scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
			facet_wrap(~st_name) +
			theme_bw() +
			theme(
				axis.text = element_text(size = 14),
				axis.title = element_text(size = 16),
				legend.title = element_text(size = 16),
				strip.text = element_text(size = 14),
				legend.text = element_text(size = 14)
			)
	}

	g <- plot_states(tmp)

	ggsave(
		filename = "plots/prop_take.jpeg",
		plot = g,
		width = 10,
		height = 6,
		units = "in",
		dpi = 300
	)

	g
}

run_bootstrap <- function(
	df_params,
	df_data,
	n_ens,
	sim_days = 365 * 5,
	n_sim = 10000
) {
	# get random sample of parameters
	set.seed(803)
	dem <- df_params[sample.int(nrow(df_params), n_ens, replace = TRUE), ]
	phi <- dem$phi_mu # survival
	zeta <- 28 * exp(dem$log_nu) / 365 # per capita growth rate
	r <- phi + zeta / 2 # finite rate of increase

	all_out <- tibble()
	counter <- 1
	while (counter < n_sim) {
		i <- sample.int(nrow(df_data), 1)

		tmp_start <- df_data |>
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

		date_seq <- seq.Date(start_date, to = start_date + sim_days, by = "4 weeks")
		n_1 <- tmp_start$abundance_estimate

		k_dens <- 30 # carrying capacity as a density
		K <- round(k_dens * area) # carrying capacity as abundance

		if (n_1 == 0 || n_1 > K) {
			next
		}

		tmp_compare <- df_data |>
			filter(propertyID %in% property, end_dates %in% date_seq) |>
			select(
				propertyID,
				property_area_km2,
				abundance_estimate,
				methods_used,
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

		process_model <- function(n_init, r, K, n_time, n_ens) {
			N <- matrix(NA, n_ens, n_time)
			N[, 1] <- n_init

			for (j in 2:n_time) {
				nm1 <- N[, j - 1]

				lambda <- nm1 + (r - 1) * nm1 * (1 - nm1 / K)
				lambda <- pmax(0, lambda)

				if (any(is.na(lambda))) {
					print(i)
					stop()
				}

				N[, j] <- rpois(n_ens, lambda)
			}
			apply(N, 2, median)
		}

		n_init <- pmin(K, rpois(n_ens, n_1))
		N <- process_model(n_init, r, K, n_time, n_ens)

		timestep_df$abundanceNoMgmt <- N

		out <- left_join(tmp_compare, timestep_df, by = c("end_dates", "year")) |>
			filter(delta_time > 0) |>
			mutate(
				simID = counter,
				density_estimate = abundance_estimate / area,
				densityNoMgmt_estimate = abundanceNoMgmt / area,
				area = area,
				K = K,
				D1 = n_1 / area,
				density_cat = if_else(D1 < 1, "Low", "Very high"),
				density_cat = if_else(D1 >= 1 & D1 < 3, "Medium", density_cat),
				density_cat = if_else(D1 >= 3 & D1 < 5, "High", density_cat),
				density_cat = factor(
					density_cat,
					levels = c("Low", "Medium", "High", "Very high")
				)
			)
		all_out <- bind_rows(all_out, out)

		counter <- counter + 1

		if (counter %% 1000 == 0) message(counter / n_sim * 100, "%")
	}

	cols_needed <- c(
		"simID",
		"delta_time",
		"st_name",
		"D1",
		"density_cat",
		"propertyID",
		"property_area_km2",
		"methods_used"
	)

	with_mgmt <- all_out |>
		select(
			all_of(cols_needed),
			density_estimate
		) |>
		mutate(
			mgmt = "With",
			delta_density = density_estimate - D1,
			percent_change = delta_density / D1 * 100
		)

	without_mgmt <- all_out |>
		select(
			all_of(cols_needed),
			densityNoMgmt_estimate
		) |>
		rename(density_estimate = densityNoMgmt_estimate) |>
		mutate(
			mgmt = "Without",
			delta_density = density_estimate - D1,
			percent_change = delta_density / D1 * 100
		)

	bind_rows(with_mgmt, without_mgmt)
}

plot_bootstrap <- function(df) {
	keep <- c(364, 728, 1092, 1456, 1820)
	years_f <- sort(keep / 364)

	my_plot <- function(df) {
		df |>
			filter(
				delta_time %in% keep,
				!st_name %in% c("North Carolina", "Virginia", "West Virginia")
			) |>
			mutate(
				weeks = delta_time / 364,
				weeks = factor(weeks, levels = years_f)
			) |>
			ggplot() +
			aes(
				x = weeks,
				y = delta_density,
				fill = mgmt
			) +
			geom_boxplot(outliers = FALSE) +
			geom_hline(yintercept = 0, linetype = "dashed") +
			labs(x = "Years", y = "Change in density", fill = "Management") +
			scale_fill_manual(values = c("With" = "#998ec3", "Without" = "#f1a340")) +
			theme_bw() +
			theme(
				legend.position = "right",
				axis.text = element_text(size = 14),
				axis.title = element_text(size = 16),
				legend.title = element_text(size = 16),
				strip.text = element_text(size = 14),
				legend.text = element_text(size = 14)
			)
	}

	g1 <- my_plot(df) + facet_grid(~density_cat)
	g2 <- my_plot(df) + facet_grid(density_cat ~ st_name)

	ggsave(
		filename = "plots/management_trajectory.jpeg",
		plot = g1,
		width = 10,
		height = 6,
		units = "in",
		dpi = 300
	)

	ggsave(
		filename = "plots/management_trajectory_by_state.jpeg",
		plot = g2,
		width = 10,
		height = 6,
		units = "in",
		dpi = 300
	)

	list(density = g1, state = g2)
}

make_epr_data <- function(mu_lambda, df_actual, model_data_last) {
	df <- mu_lambda |>
		rename(lambda = mu) |>
		select(-n)

	df_epr <- df_actual |>
		left_join(df) |>
		filter(!is.na(lambda), lambda > 0) |>
		mutate(
			EPR = lambda * (1 - mu_prop_take) - 1,
			epr_cat = if_else(EPR > 0, "Grow", "Shrink"),
			epr_cat = if_else(EPR <= -0.5, "Slash", epr_cat),
			epr_cat = if_else(EPR == -1, "Eradicate", epr_cat),
			epr_cat = factor(
				epr_cat,
				levels = c("Grow", "Shrink", "Slash", "Eradicate")
			)
		) |>
		filter(epr_cat != "Eradicate")

	dat <- model_data_last |>
		select(propertyID, property_area_km2) |>
		# mutate(property_area_km2 = log(property_area_km2)) |>
		distinct()

	df_plot <- df_epr |>
		left_join(dat, by = "propertyID") |>
		mutate(effort_per_unit_per_area = mu_total_effort_per / property_area_km2)
}

plot_effort <- function(epr_data) {
	ept_cols <- c(
		"Grow" = "#d73027",
		"Shrink" = "#4393c3",
		"Slash" = "#2166ac",
		"Eradicate" = "#053061"
	)

	my_box <- function(x) {
		epr_data |>
			ggplot() +
			aes(x = .data[[x]], y = method, color = epr_cat) +
			labs(color = "Population trajectory") +
			geom_boxplot(outliers = FALSE) +
			scale_color_manual(
				values = ept_cols
			) +
			theme_bw() +
			theme(
				axis.title.y = element_blank(),
				axis.text = element_text(size = 10),
				axis.title = element_text(size = 14),
				legend.title = element_text(size = 12),
				strip.text = element_text(size = 12),
				legend.text = element_text(size = 10)
			)
	}

	g1 <- my_box("mu_prop_take") +
		labs(x = "Average APR") +
		scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))
	g2 <- my_box("property_area_km2") +
		scale_x_log10() +
		labs(x = expression("Property area " ~ (km^2)))
	g3 <- my_box("effort_per_unit_per_area") +
		scale_x_log10() +
		# fmt: skip
		labs(x = expression("Avg. TEPU/"~km^2))

	gg <- ggarrange(
		g1,
		g2,
		g3,
		nrow = 1,
		common.legend = TRUE,
		legend = "right",
		labels = c("A", "B", "C"),
		align = "hv"
	)

	d1 <- expand_grid(
		epr = c(-1, -0.5, 0),
		lambda = seq(0, 2, by = 0.01)
	) |>
		mutate(apr = 1 - ((1 + epr) / lambda)) |>
		filter(apr >= 0, apr <= 1) |>
		mutate(
			epr_cat = if_else(epr == 0, "Maintenance", "50% reduction"),
			epr_cat = if_else(epr == -1, "Eradication", epr_cat),
			epr_cat = factor(
				epr_cat,
				levels = c("Maintenance", "50% reduction", "Eradication")
			)
		)

	g4 <- expand_grid(
		apr = seq(0, 1, by = 0.05),
		lambda = seq(0, 2, by = 0.05)
	) |>
		mutate(
			epr = lambda * (1 - apr) - 1
		) |>
		ggplot() +
		aes(x = lambda, y = apr, fill = epr) +
		geom_tile() +
		geom_line(
			data = d1,
			aes(x = lambda, y = apr, linetype = epr_cat),
			linewidth = 1,
			color = "#f33"
		) +
		scale_fill_gradientn(
			breaks = c(1, 0.5, 0, -0.5, -1),
			labels = c("Double", "", "Maintain", "", "Eradicate"),
			colors = c("#053061", "#2166ac", "#f7f7f7", "#f4a582", "#d73027")
		) +
		scale_linetype_manual(
			values = c(
				"Maintenance" = "dotted",
				"50% reduction" = "dashed",
				"Eradication" = "solid"
			)
		) +
		labs(
			y = "APR",
			x = "Populuation growth rate",
			fill = "EPT",
			linetype = element_blank()
		) +
		theme_bw() +
		theme(
			axis.text = element_text(size = 12),
			axis.title = element_text(size = 14),
			legend.title = element_text(size = 12),
			strip.text = element_text(size = 12),
			legend.text = element_text(size = 10)
		)

	g <- ggarrange(
		gg,
		g4,
		nrow = 2,
		labels = c("", "D"),
		heights = c(2, 1)
	)

	ggsave(
		filename = "plots/epr_effort_compAdditive.jpeg",
		plot = g,
		width = 10,
		height = 8,
		units = "in",
		dpi = 300
	)
	g
}

extract_parameters <- function(post_1, post_2) {
	post_1 <- read_rds(post_1) |>
		mutate(round = "First") |>
		pivot_longer(
			cols = -round,
			names_to = "parameter",
			values_to = "value"
		)

	post_11 <- read_rds(post_2) |>
		mutate(round = "Last") |>
		pivot_longer(
			cols = -round,
			names_to = "parameter",
			values_to = "value"
		)

	post_all <- bind_rows(post_1, post_11) |>
		mutate(
			method_idx = as.numeric(stringr::str_extract(parameter, "(?<=\\[)\\d")),
			method_idx = if_else(
				grepl("log_gamma", parameter) | grepl("p_mu", parameter),
				method_idx + 3,
				method_idx
			),
			land_idx = as.numeric(stringr::str_extract(parameter, "(?<=\\,\\s)\\d"))
		)

	method_vector <- c(
		"Ground-shooting",
		"Fixed wing",
		"Helicopter",
		"Snare",
		"Trap"
	)

	method_lookup <- tibble(
		method_idx = 1:5,
		method = method_vector
	)

	land_lookup <- tibble(
		land_idx = 1:3,
		land = c("Road density", "Ruggedness", "Canopy cover")
	)

	param_names <- function(df) {
		df |>
			mutate(
				name = case_when(
					grepl("gamma", parameter) ~ "Gamma",
					grepl("rho", parameter) ~ "Rho",
					grepl("p_mu", parameter) ~ "Omega",
					grepl("beta", parameter) ~ "Beta",
					grepl("phi_mu", parameter) ~ "Phi",
					grepl("psi_phi", parameter) ~ "Psi",
					grepl("log_nu", parameter) ~ "Nu"
				)
			)
	}

	post_all |>
		left_join(method_lookup, by = "method_idx") |>
		left_join(land_lookup, by = "land_idx") |>
		mutate(
			land = if_else(grepl("beta1", parameter), "Intercept", land),
			value = if_else(grepl("log", parameter), exp(value), value),
			value = if_else(grepl("p_mu", parameter), nimble::ilogit(value), value)
		) |>
		param_names() |>
		group_by(name, round, parameter, method, land) |>
		summarise(
			med = median(value),
			low = quantile(value, 0.025),
			high = quantile(value, 0.975),
			.groups = "drop"
		) |>
		ungroup()
}

create_parameter_table <- function(post_all) {
	post_all |>
		mutate(
			est = paste0(
				round(med, 2),
				" (",
				round(low, 2),
				", ",
				round(high, 2),
				")"
			)
		) |>
		select(parameter, round, method, land, name, est) |>
		pivot_wider(
			names_from = round,
			values_from = est
		)
}

plot_params_compare <- function(post_all) {
	round_cols <- c(
		"First" = "#1b9e77",
		"Last" = "#7570b3"
	)

	plot_lines <- function(df, m) {
		df |>
			filter(method == m) |>
			ggplot() +
			aes(x = round, ymin = low, ymax = high, y = med) +
			geom_pointrange() +
			facet_wrap(~name, scales = "free") +
			labs(
				x = "Round",
				y = "Value",
				title = m
			) +
			scale_color_manual(values = round_cols) +
			theme_bw() +
			theme(
				axis.text = element_text(size = 14),
				axis.title = element_text(size = 16),
				legend.title = element_text(size = 16),
				title = element_text(size = 12),
				legend.text = element_text(size = 14),
				strip.text = element_text(size = 12)
			)
	}

	g1 <- post_all |> filter(is.na(land)) |> plot_lines("Ground-shooting")
	g2 <- post_all |> filter(is.na(land)) |> plot_lines("Fixed wing")
	g3 <- post_all |> filter(is.na(land)) |> plot_lines("Helicopter")

	ggarrange(
		g1,
		g2,
		g3,
		ncol = 3,
		nrow = 1,
		legend = "bottom",
		labels = "AUTO"
	)

	ggsave(
		filename = "plots/posterior_comparison_active_methods.jpeg",
		width = 7,
		height = 5,
		units = "in",
		dpi = 300
	)

	g1 <- post_all |> filter(is.na(land)) |> plot_lines("Trap")
	g2 <- post_all |> filter(is.na(land)) |> plot_lines("Snare")

	ggarrange(
		g1,
		g2,
		nrow = 2,
		legend = "bottom",
		labels = "AUTO"
	)

	ggsave(
		filename = "plots/posterior_comparison_passive_methods.jpeg",
		width = 7,
		height = 5,
		units = "in",
		dpi = 300
	)

	g1 <- post_all |>
		filter(!is.na(land)) |>
		plot_lines("Ground-shooting") +
		facet_wrap(~land, scales = "free")
	g2 <- post_all |>
		filter(!is.na(land)) |>
		plot_lines("Snare") +
		facet_wrap(~land, scales = "free")
	g3 <- post_all |>
		filter(!is.na(land)) |>
		plot_lines("Helicopter") +
		facet_wrap(~land, scales = "free")
	g4 <- post_all |>
		filter(!is.na(land)) |>
		plot_lines("Fixed wing") +
		facet_wrap(~land, scales = "free")
	g5 <- post_all |>
		filter(!is.na(land)) |>
		plot_lines("Trap") +
		facet_wrap(~land, scales = "free")

	ggarrange(
		g1,
		g2,
		g3,
		g4,
		g5,
		ncol = 3,
		nrow = 2,
		legend = "bottom",
		labels = "AUTO"
	)

	ggsave(
		filename = "plots/posterior_comparison_betas.jpeg",
		width = 12,
		height = 10,
		units = "in",
		dpi = 300
	)
}

density_summaries <- function(data_density, data_model) {
	properties <- data_model |>
		select(
			propertyID,
			agrp_prp_id,
			property_area_km2,
			end_dates,
			STATE_NAME,
			cnty_name,
			primary_period,
			property
		) |>
		distinct() |>
		mutate(year = year(end_dates))

	data <- left_join(data_density, properties) |>
		rename(density = `0.5`) |>
		mutate(abundance = round(density * property_area_km2))

	get_ci <- function(df) {
		df |>
			summarise(
				n = n(),
				med_dens = round(median(density), 2),
				low_dens = round(quantile(density, 0.05), 2),
				high_dens = round(quantile(density, 0.95), 2),
				q1_dens = round(quantile(density, 0.25), 2),
				q3_dens = round(quantile(density, 0.75), 2),
				total_abun = sum(abundance)
			) |>
			ungroup()
	}

	# all
	overall <- data |>
		get_ci()

	# by year
	dens_year <- data |>
		group_by(year) |>
		get_ci()

	list(overall = overall, by_year = dens_year)
}

plot_density_year <- function(dens_ls) {
	dens_ls$by_year |>
		ggplot(aes(x = year, y = med_dens, ymin = low_dens, ymax = high_dens)) +
		geom_linerange(aes(ymin = q1_dens, ymax = q3_dens), linewidth = 2) +
		geom_linerange() +
		geom_point(size = 5) +
		geom_point(size = 3, color = "white") +
		labs(
			x = "Year",
			y = "Median Density (pigs/km2)",
			title = "Median Density of Pigs by Year"
		) +
		theme_bw() +
		theme(
			axis.text = element_text(size = 14),
			axis.title = element_text(size = 18),
			title = element_text(size = 18)
		)
}
