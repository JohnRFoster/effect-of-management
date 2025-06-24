# read rds files
get_rds <- function(file) {
	read_rds(file)
}


# rename methods for better plotting
fix_method_names <- function(df) {
	df |>
		mutate(
			method = if_else(method == "FIREARMS", "Sharpshooting", method),
			method = if_else(method == "FIXED WING", "Fixed wing", method),
			method = if_else(method == "HELICOPTER", "Helicopter", method),
			method = if_else(method == "SNARE", "Snare", method),
			method = if_else(method == "TRAPS", "Trap", method)
		)
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
data_prep <- function(df, df_density, cutoff_date) {
	tmp <- df |>
		fix_method_names() |>
		make_model_data()

	df_methods <- make_methods_used(df)

	tmp |>
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
		filter(!is.na(take_m1))
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

plot_mean_lambdas <- function(df) {
	y_order <- df |>
		group_by(st_name) |>
		summarise(y = median(mu), n = n()) |>
		arrange(y) |>
		mutate(st_name2 = paste0(st_name, " (n=", n, ")")) |>
		select(-n)

	y_fac <- y_order$st_name2

	df |>
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

prop_take_by_method <- function(df, df_mu_lambda, df_density) {
	tmp <- df |>
		fix_method_names()

	left_join(tmp, df_mu_lambda) |>
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
		left_join(df_density) |>
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
}

prop_to_remove <- function(df_realized_growth) {
	df_realized_growth |>
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
}

plot_percent_take <- function(df_actual, df_theoretical) {
	tmp <- bind_rows(df_actual, df_theoretical) |>
		# removing these states because their properties are all declining
		filter(!st_name %in% c("NORTH CAROLINA", "VIRGINIA", "WEST VIRGINIA"))

	plot_states <- function(df, states) {
		ta <- df |>
			filter(st_name %in% states, group == "Actual") |>
			mutate(
				method = factor(
					method,
					levels = c(
						"Fixed wing",
						"Helicopter",
						"Sharpshooting",
						"Snare",
						"Trap"
					)
				)
			)

		th <- df |>
			filter(st_name %in% states, group == "Theoretical") |>
			mutate(
				method = if_else(grepl("maintain", method), "Maintenance", method),
				method = if_else(grepl("reduction", method), "50% reduction", method),
				method = factor(
					method,
					levels = c(
						"50% reduction",
						"Maintenance"
					)
				)
			)

		g1 <- ta |>
			ggplot() +
			aes(x = mu, y = method) +
			geom_boxplot() +
			facet_wrap(~st_name) +
			coord_cartesian(xlim = c(0, 1)) +
			labs(
				x = "Proportion of population removed in a 28-day period",
				y = "Method",
				color = "Actual"
			) +
			scale_y_discrete(drop = FALSE) +
			theme_bw() +
			theme(axis.title.x = element_blank(), axis.text.x = element_blank())

		g2 <- th |>
			ggplot() +
			aes(x = mu, y = method) +
			geom_boxplot() +
			facet_wrap(~st_name) +
			coord_cartesian(xlim = c(0, 1)) +
			labs(
				x = "Proportion of population removed in a 28-day period",
				y = "Goal",
				color = "Theoretical"
			) +
			theme_bw() +
			theme(strip.text = element_blank())

		ggarrange(
			g1,
			g2,
			ncol = 1,
			common.legend = TRUE,
			legend = "right",
			heights = c(5, 2)
		)
	}

	g1 <- plot_states(tmp, c("FLORIDA", "GEORGIA", "LOUISIANA"))
	g2 <- plot_states(tmp, c("MISSISSIPPI", "MISSOURI", "NEW MEXICO"))
	g3 <- plot_states(tmp, c("OKLAHOMA", "SOUTH CAROLINA", "TEXAS"))

	ggarrange(
		g1,
		g2,
		g3,
		ncol = 1,
		common.legend = TRUE
	)
}
