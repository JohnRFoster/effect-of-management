#' Calculate finite population growth rates
#'
#' Calculates realized population growth between consecutive primary periods after
#' accounting for removals in the previous period.
#'
#' @param df A data frame containing `abundance_estimate`, `N_m1`, `take_m1`, and
#'   `delta_pp` columns.
#'
#' @return A data frame with intermediate `Z`, `lambda_star`, and `lambda` columns,
#'   excluding rows with undefined lambda values.
calculate_lambda <- function(df) {
	df |>
		mutate(
			Z = N_m1 - take_m1,
			lambda_star = abundance_estimate / Z,
			lambda = lambda_star^(1 / delta_pp)
		) |>
		filter(!is.nan(lambda))
}

#' Calculate realized growth rates
#'
#' Creates primary-period intervals for each property and estimates realized
#' population growth rates for properties with more than one abundance estimate.
#'
#' @param df A prepared model data frame with property identifiers, end dates,
#'   total take, and abundance estimates.
#'
#' @return A data frame of realized growth-rate records with lagged abundance,
#'   lagged take, primary-period spacing, and lambda estimates.
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

#' Create data for realized growth calculations
#'
#' Joins prepared model data to primary-period indices, keeps properties with
#' repeated estimates, and adds lagged abundance, lagged take, and elapsed primary
#' periods.
#'
#' @param df A prepared model data frame containing property-level abundance and
#'   take summaries.
#' @param df_pp A data frame of all property primary-period dates, usually from
#'   `make_all_timesteps()`.
#' @param p2 A vector of property IDs with more than one abundance estimate.
#'
#' @return A data frame filtered to intervals with non-missing elapsed primary
#'   periods.
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

#' Create posterior parameter table
#'
#' Formats posterior summaries as median and interval strings and widens the first
#' and last model rounds into separate columns.
#'
#' @param post_all A posterior summary data frame from `extract_parameters()`.
#'
#' @return A wide data frame with parameter metadata and formatted first/last
#'   posterior summaries.
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

#' Prepare analysis data
#'
#' Standardizes method and state names, aggregates removals and effort to the
#' property-primary-period scale, joins density summaries and methods used, and
#' derives abundance estimates.
#'
#' @param df A raw model data frame containing removal events, effort, methods,
#'   property metadata, and primary periods.
#' @param df_density A density summary data frame with quantile columns named
#'   `0.05`, `0.5`, and `0.95`.
#'
#' @return A property-primary-period data frame ready for downstream analysis and
#'   plotting.
data_prep <- function(df, df_density) {
	tmp <- df |>
		make_model_data()

	df_methods <- make_methods_used(df)
	df_density <- df_density

	tmp |>
		left_join(df_density) |>
		left_join(df_methods) |>
		mutate(
			year = year(end_dates),
			abundance_estimate = round(`0.5` * property_area_km2)
		) |>
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
			methods_used,
			n_id
		)
}

#' Summarize density estimates
#'
#' Joins density summaries to property metadata, converts median densities to
#' abundance estimates, and summarizes density intervals overall and by year.
#'
#' @param data_density A density summary data frame with a median density column
#'   named `0.5`.
#' @param data_model A model data frame containing property metadata and dates.
#'
#' @return A list with `overall` and `by_year` density summary data frames.
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

#' Extract posterior parameter summaries
#'
#' Reads first- and last-round posterior samples, reshapes them to long format,
#' maps model parameter indices to method and landscape labels, transforms log and
#' logit-scale parameters, and calculates posterior medians and intervals.
#'
#' @param post_1 File path to the first-round posterior samples RDS file.
#' @param post_2 File path to the last-round posterior samples RDS file.
#'
#' @return A data frame of posterior parameter summaries by parameter, method,
#'   landscape covariate, and fitting round.
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

#' Standardize method names
#'
#' Converts all-caps method labels into publication-friendly names and maps state
#' abbreviations/names using the shared `state_names` lookup.
#'
#' @param df A data frame containing `method` and `st_name` columns.
#'
#' @return The input data frame with standardized method labels and state names.
fix_method_names <- function(df) {
	df |>
		mutate(
			method = case_when(
				method == "FIREARMS" ~ "Ground-shooting",
				method == "FIXED WING" ~ "Fixed wing",
				method == "HELICOPTER" ~ "Helicopter",
				method == "SNARE" ~ "Snare",
				method == "TRAPS" ~ "Trap",
				TRUE ~ method
			)
		)
}

#' Standardize Texas state names
#'
#' Converts Texas state names to a standardized format.
#'
#' @param df A data frame containing `st_name` column.
#'
#' @return The input data frame with standardized Texas state names.
fix_state_names <- function(df) {
	mutate(df, st_name = if_else(grepl("Texas", st_name), "Texas", st_name))
}

fix_state_names2 <- function(df) {
	df |>
		mutate(
			st_name = case_when(
				st_name == "Texas_south_central_semiarid_prairies" ~ "Texas SCSP",
				st_name == "Texas_texas-Louisiana_coastal_plain" ~ "Texas TLCP",
				st_name == "Texas_tamaulipas-Texas_semiarid_plain" ~ "Texas TTSP",
				st_name == "Texas_warm_deserts" ~ "Texas WD",
				st_name == "Texas_southeastern_usa_plains" ~ "Texas SEUP",
				TRUE ~ st_name
			)
		)
}

#' Read density summaries
#'
#' @param file Path to the RDS file containing density summaries.
#'
#' @return A data frame of density summaries with standardized state names.
get_density <- function(file) {
	read_rds(file) |>
		mutate(st_name = str_to_title(st_name)) |>
		select(-primary_period)
}

#' Read model data
#'
#' @param file Path to the RDS file containing model data.
#'
#' @return A data frame of model data with standardized state names and method labels.
get_model_data <- function(file) {
	read_rds(file) |>
		mutate(st_name = str_to_title(st_name)) |>
		fix_method_names()
}

#' Create complete property timesteps
#'
#' Creates a 28-day date sequence from the first to last observed primary-period
#' end date for each property and assigns sequential primary-period indices.
#'
#' @param df A prepared model data frame containing `propertyID`, `end_dates`,
#'   `st_name`, `county_code`, and `property_area_km2` columns.
#'
#' @return A data frame of property-date combinations with metadata and `pp`
#'   primary-period indices.
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

#' Create effort-population-rate data
#'
#' Combines mean realized growth rates with method-specific removal summaries to
#' calculate effort-population trajectory categories for plotting.
#'
#' @param mu_lambda A data frame of mean growth rates from `mean_growth()`.
#' @param df_actual A data frame of observed method-level removal proportions from
#'   `prop_take_by_method()`.
#' @return A data frame with estimated population trajectory (`EPR`), trajectory
#'   category, property area, and effort per unit area.
make_epr_data <- function(df, df_actual) {
	df_epr <- df_actual |>
		left_join(df) |>
		filter(!is.na(lambda), lambda > 0) |>
		mutate(
			EPR = lambda * (1 - mu_prop_take) - 1,
			epr_cat = if_else(EPR > 0, "Increasing", "Moderately declining"),
			epr_cat = if_else(EPR <= -0.5, "Strongly declining", epr_cat),
			epr_cat = if_else(EPR == -1, "Eradicate", epr_cat),
			epr_cat = factor(
				epr_cat,
				levels = c(
					"Increasing",
					"Moderately declining",
					"Strongly declining",
					"Eradicate"
				)
			)
		) |>
		filter(epr_cat != "Eradicate")

	df_plot <- df_epr |>
		mutate(effort_per_unit_per_area = mu_total_effort_per / property_area_km2)
	df_plot
}

#' Summarize methods used by primary period
#'
#' Collapses distinct management methods used at each property and primary period
#' into a comma-separated methods-used field.
#'
#' @param df A model data frame containing `propertyID`, `primary_period`, and
#'   `method` columns.
#'
#' @return A data frame with one row per property-primary-period combination and a
#'   `methods_used` column.
make_methods_used <- function(df) {
	df |>
		select(propertyID, primary_period, method) |>
		distinct() |>
		pivot_wider(names_from = method, values_from = method) |>
		unite(method, -c(propertyID, primary_period), sep = ", ", na.rm = TRUE) |>
		rename(methods_used = method)
}

#' Aggregate model data
#'
#' Aggregates event-level removals and effort to the property-primary-period scale
#' and derives take density.
#'
#' @param df A model data frame containing property metadata, primary periods,
#'   removal counts, effort-per-unit values, and property area.
#'
#' @return A property-primary-period data frame with total take, total effort per
#'   unit, event counts, take density, and missing farm bill values set to zero.
make_model_data <- function(df) {
	df |>
		group_by(
			propertyID,
			agrp_prp_id,
			start_dates,
			end_dates,
			st_name,
			cnty_name,
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
		mutate(take_density = total_take / property_area_km2)
}

#' Calculate mean growth by property
#'
#' Averages realized lambda values for each property and state.
#'
#' @param df A realized growth data frame containing `propertyID`, `st_name`, and
#'   `lambda` columns.
#'
#' @return A data frame with sample size `n` and mean growth rate `mu` for each
#'   property-state combination.
mean_growth <- function(df) {
	df |>
		group_by(propertyID, st_name) |>
		summarise(
			n = n(),
			mu = mean(lambda)
		) |>
		ungroup()
}

#' Plot land-cover parameter intervals
#'
#' Filters parameter summaries to a named landscape covariate and plots matching
#' prior and posterior intervals.
#'
#' @param df_params A posterior parameter summary data frame.
#' @param df_priors A prior parameter summary data frame.
#' @param x A character pattern matched against `land_type`.
#' @param w Width passed to `position_dodge()` for points and lineranges.
#' @param lw Linerange width passed to `geom_linerange()`.
#'
#' @return A ggplot object for the matching landscape parameter intervals.
my_linerange <- function(df_params, df_priors, x, w = 0.5, lw = 1.25) {
	df <- prepare_prior_posterior_params(df_params, df_priors) |>
		filter(grepl(x, land_type))

	plot_prior_posterior_linerange(df, w, lw)
}

#' Plot node-matched parameter intervals
#'
#' Filters parameter summaries by node name or node pattern and plots matching
#' prior and posterior intervals.
#'
#' @param df_params A posterior parameter summary data frame.
#' @param df_priors A prior parameter summary data frame.
#' @param node_pattern Optional character pattern matched against `node` when
#'   `nodes` is `NULL`.
#' @param nodes Optional character vector of exact node names to retain.
#' @param w Width passed to `position_dodge()` for points and lineranges.
#' @param lw Linerange width passed to `geom_linerange()`.
#'
#' @return A ggplot object for the matching node parameter intervals.
my_node_linerange <- function(
	df_params,
	df_priors,
	node_pattern = NULL,
	nodes = NULL,
	w = 0.5,
	lw = 1.25
) {
	df <- prepare_prior_posterior_params(df_params, df_priors)

	if (is.null(nodes)) {
		df <- df |>
			filter(grepl(node_pattern, node))
	} else {
		df <- df |>
			filter(node %in% nodes)
	}

	plot_prior_posterior_linerange(df, w, lw)
}

#' Plot management bootstrap trajectories
#'
#' Creates density-change boxplots comparing simulated trajectories with and
#' without management, both by starting-density class and by state.
#'
#' @param df A bootstrap output data frame from `run_bootstrap()`.
#'
#' @return A list containing the density-class plot and the state-by-density-class
#'   plot. Also writes PDF files to the `plots/` directory.
plot_bootstrap <- function(df) {
	keep <- c(364, 728, 1092, 1456, 1820)
	years_f <- sort(keep / 364)

	my_plot <- function(df) {
		df |>
			filter(
				delta_time %in% keep,
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
		filename = "plots/management_trajectory.tiff",
		plot = g1,
		width = 10,
		height = 6,
		units = "in",
		dpi = 300
	)

	ggsave(
		filename = "plots/management_trajectory_by_state.tiff",
		plot = g2,
		width = 10,
		height = 6,
		units = "in",
		dpi = 300
	)

	list(density = g1, state = g2)
}

#' Plot demographic parameter intervals by state
#'
#' Creates a faceted prior/posterior interval plot for selected demographic nodes
#' and relabels node names for publication.
#'
#' @param df_params A posterior parameter summary data frame.
#' @param df_priors A prior parameter summary data frame.
#'
#' @return A ggplot object. Also writes `plots/demographic_parameters.tiff`.
plot_demographic_parameters <- function(df_params, df_priors) {
	nodes <- c("phi_mu", "psi_phi", "zeta")

	g <- my_node_linerange(df_params, df_priors, nodes = nodes) +
		facet_wrap(~node, scales = "free_x") +
		labs(
			x = "Parameter estimate",
			y = "State",
			title = "Demographic parameters by state",
		) +
		theme(legend.position = "none")

	g$data <- g$data |>
		mutate(
			node = case_when(
				node == "phi_mu" ~ "Mean 28-day survival rate",
				node == "psi_phi" ~ "Shrinkage",
				node == "zeta" ~ "28-day per-capita recruitment rate",
				node == "nu" ~ "Litter size",
				node == "lambda" ~ "28-day growth rate",
				node == "lambda_annual" ~ "Annual growth rate",
				.default = node
			)
		)

	ggsave(
		filename = "plots/demographic_parameters.tiff",
		plot = g,
		width = 10,
		height = 6,
		units = "in",
		dpi = 300
	)

	g
}

#' Plot annual density summaries
#'
#' Creates a plot of median density by year with interquartile and 90 percent
#' interval ranges.
#'
#' @param dens_ls A list of density summaries from `density_summaries()`.
#'
#' @return A ggplot object showing annual median density and uncertainty intervals.
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

#' Plot effort and population trajectory relationships
#'
#' Creates a multi-panel figure comparing annual proportion removed, property
#' area, and effort per unit area by trajectory category, along with an APR and
#' growth-rate response surface.
#'
#' @param epr_data A data frame from `make_epr_data()`.
#'
#' @return A ggplot/ggarrange object. Also writes a PDF file to the `plots/`
#'   directory.
plot_effort <- function(epr_data) {
	ept_cols <- c(
		"Increasing" = "#d73027",
		"Moderately declining" = "#4393c3",
		"Strongly declining" = "#2166ac",
		"Eradicate" = "#053061"
	)

	my_box <- function(x) {
		w <- 0.75
		epr_data |>
			group_by(method, epr_cat) |>
			summarise(
				median = median(.data[[x]], na.rm = TRUE),
				l = quantile(.data[[x]], 0.05, na.rm = TRUE),
				q1 = quantile(.data[[x]], 0.25, na.rm = TRUE),
				q3 = quantile(.data[[x]], 0.75, na.rm = TRUE),
				u = quantile(.data[[x]], 0.95, na.rm = TRUE)
			) |>
			ggplot() +
			aes(x = median, y = method, color = epr_cat, group = epr_cat) +
			labs(color = "Population trajectory") +
			geom_linerange(
				aes(xmin = l, xmax = u),
				linewidth = 1,
				position = position_dodge(width = w)
			) +
			geom_linerange(
				aes(xmin = q1, xmax = q3),
				linewidth = 2,
				position = position_dodge(width = w)
			) +
			geom_point(
				aes(x = median),
				size = 3,
				position = position_dodge(width = w)
			) +
			scale_color_manual(
				values = ept_cols
			) +
			theme_bw() +
			theme(
				axis.title.y = element_blank(),
				axis.text = element_text(size = 10),
				axis.title = element_text(size = 12),
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
		filename = "plots/epr_effort_compAdditive.tiff",
		plot = g,
		width = 10,
		height = 8,
		units = "in",
		dpi = 300
	)
	g
}

#' Plot effort-saturation effects on removal rate
#'
#' Creates a faceted prior/posterior interval plot for `gamma` nodes, interpreted
#' as saturation effects of effort by method and state.
#'
#' @param df_params A posterior parameter summary data frame.
#' @param df_priors A prior parameter summary data frame.
#'
#' @return A ggplot object. Also writes `plots/effort_saturation.tiff`.
plot_effort_saturation <- function(df_params, df_priors) {
	g <- my_node_linerange(df_params, df_priors, "gamma") +
		facet_wrap(~method_names, scales = "free_x") +
		labs(
			x = "Effect",
			y = "State",
			title = "Saturation effect of effort on removal rate",
		) +
		theme(legend.position = "none")

	ggsave(
		filename = "plots/effort_saturation.tiff",
		plot = g,
		width = 10,
		height = 6,
		units = "in",
		dpi = 300
	)

	g
}

#' Plot effort-scaling effects on removal rate
#'
#' Creates a faceted prior/posterior interval plot for `rho` nodes, interpreted as
#' scaling effects of effort by method and state.
#'
#' @param df_params A posterior parameter summary data frame.
#' @param df_priors A prior parameter summary data frame.
#'
#' @return A ggplot object. Also writes `plots/effort_scaling.tiff`.
plot_effort_scaling <- function(df_params, df_priors) {
	g <- my_node_linerange(df_params, df_priors, "rho") +
		facet_wrap(~method_names, scales = "free_x") +
		labs(
			x = "Effect",
			y = "State",
			title = "Scaling effect of effort on removal rate",
		) +
		theme(legend.position = "none")

	ggsave(
		filename = "plots/effort_scaling.tiff",
		plot = g,
		width = 10,
		height = 6,
		units = "in",
		dpi = 300
	)

	g
}

#' Plot unique-area effort effects on removal rate
#'
#' Creates a faceted prior/posterior interval plot for `p_mu` nodes, interpreted
#' as the proportion of unique area surveyed by additional deployed units.
#'
#' @param df_params A posterior parameter summary data frame.
#' @param df_priors A prior parameter summary data frame.
#'
#' @return A ggplot object. Also writes `plots/effort_unique_area.tiff`.
plot_effort_unique_area <- function(df_params, df_priors) {
	g <- my_node_linerange(df_params, df_priors, "p_mu") +
		facet_wrap(~method_names) +
		labs(
			x = "Effect",
			y = "State",
			title = "Proportion of unique area surveyed by additional units deployed",
		) +
		theme(legend.position = "none")

	ggsave(
		filename = "plots/effort_unique_area.tiff",
		plot = g,
		width = 10,
		height = 6,
		units = "in",
		dpi = 300
	)

	g
}

#' Plot mean lambda by state
#'
#' Creates a state-level boxplot of property mean growth rates and saves it to the
#' plots directory.
#'
#' @param df A mean-growth data frame from `mean_growth()`.
#'
#' @return A ggplot object. Also writes `plots/lambda_hat_by_state.tiff`.
plot_mean_lambdas <- function(df) {
	y_order <- df |>
		group_by(st_name) |>
		summarise(y = median(mu, na.rm = TRUE), n = n()) |>
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
		filename = "plots/lambda_hat_by_state.tiff",
		plot = g,
		width = 10,
		height = 6,
		units = "in",
		dpi = 300
	)

	g
}

#' Plot posterior parameter comparisons
#'
#' Creates posterior comparison plots for active methods, passive methods, and
#' landscape coefficients across the first and last fitting rounds.
#'
#' @param post_all A posterior summary data frame from `extract_parameters()`.
#'
#' @return Invisibly returns `NULL`. Writes posterior comparison PDF files to the
#'   `plots/` directory.
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
		filename = "plots/posterior_comparison_active_methods.tiff",
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
		filename = "plots/posterior_comparison_passive_methods.tiff",
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
		filename = "plots/posterior_comparison_betas.tiff",
		width = 12,
		height = 10,
		units = "in",
		dpi = 300
	)
}

#' Plot annual proportion removed by method
#'
#' Creates boxplots of mean proportion removed by management method for each
#' state, excluding states where all properties are declining.
#'
#' @param df_actual A data frame of observed method-level removal proportions from
#'   `prop_take_by_method()`.
#'
#' @return A ggplot object. Also writes `plots/prop_take.tiff`.
plot_percent_take <- function(df_actual) {
	tmp <- df_actual

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
			) |>
			group_by(st_name, method) |>
			reframe(
				y = median(mu_prop_take),
				n = n(),
				interval_type = if_else(n < 10, "Range", "90% interval"),
				interval_type = if_else(n == 1, "Point", interval_type),
				ymin = if_else(
					n < 10,
					min(mu_prop_take),
					quantile(mu_prop_take, 0.05)
				),
				ymax = if_else(
					n < 10,
					max(mu_prop_take),
					quantile(mu_prop_take, 0.95)
				)
			)

		sts <- rev(sort(unique(ta$st_name)))

		ta |>
			mutate(st_name = factor(st_name, levels = sts)) |>
			ggplot() +
			aes(
				x = y,
				xmin = ymin,
				xmax = ymax,
				y = method,
				color = interval_type
			) +
			geom_point(size = 3) +
			geom_linerange(
				linewidth = 1,
				stat = "identity"
			) +
			facet_wrap(~st_name) +
			labs(
				x = "Mean proportion of individuals removed in a 28-day period",
				y = "Method",
				color = "Statistic"
			) +
			scale_color_manual(
				values = c(
					"Range" = "#1b9e77",
					"90% interval" = "#7570b3",
					"Point" = "#d95f02"
				)
			) +
			theme_bw() +
			theme(
				axis.text.y = element_text(size = 12),
				axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
				axis.title = element_text(size = 16),
				legend.title = element_text(size = 16),
				strip.text = element_text(size = 12),
				legend.text = element_text(size = 14)
			)
	}

	g <- plot_states(tmp)

	ggsave(
		filename = "plots/prop_take.tiff",
		plot = g,
		width = 10,
		height = 6,
		units = "in",
		dpi = 300
	)

	g
}

#' Plot prior and posterior parameter intervals
#'
#' Creates a common linerange plot comparing prior and posterior median estimates
#' with 90 percent intervals by state.
#'
#' @param df A prepared prior/posterior parameter data frame from
#'   `prepare_prior_posterior_params()`.
#' @param w Width passed to `position_dodge()` for points and lineranges.
#' @param lw Linerange width passed to `geom_linerange()`.
#'
#' @return A ggplot object with prior and posterior parameter intervals.
plot_prior_posterior_linerange <- function(df, w = 0.5, lw = 1.25) {
	ggplot(df) +
		aes(
			x = median,
			xmin = quantile_0.05,
			xmax = quantile_0.95,
			y = st_name,
			colour = group,
			group = st_name
		) +
		geom_point(position = position_dodge(width = w), size = 3) +
		geom_linerange(position = position_dodge(width = w), linewidth = lw) +
		scale_color_manual(
			values = c("Prior" = "#7570b3", "Posterior" = "#1b9e77")
		) +
		theme_bw() +
		labs_pubr()
}

#' Plot state intercept removal rates by method
#'
#' Creates a faceted prior/posterior interval plot of method-specific intercepts,
#' interpreted as mean removal rates by state.
#'
#' @param df_params A posterior parameter summary data frame.
#' @param df_priors A prior parameter summary data frame.
#'
#' @return A ggplot object. Also writes `plots/state_intercepts_rem_rate.tiff`.
plot_state_intercepts_rem_rate <- function(df_params, df_priors) {
	g <- my_linerange(df_params, df_priors, "Intercept") +
		facet_wrap(~method_names, scales = "free_x") +
		labs(
			x = "Mean removal rate",
			y = "State",
			title = "Mean removal rate by method",
		) +
		theme(legend.position = "none")

	ggsave(
		filename = "plots/state_intercepts_rem_rate.tiff",
		plot = g,
		width = 10,
		height = 6,
		units = "in",
		dpi = 300
	)

	g
}

#' Plot canopy-cover effects on removal rate
#'
#' Creates a faceted prior/posterior interval plot for canopy-cover effects by
#' method and state.
#'
#' @param df_params A posterior parameter summary data frame.
#' @param df_priors A prior parameter summary data frame.
#'
#' @return A ggplot object. Also writes `plots/states_canopy_cover.tiff`.
plot_states_canopy_cover <- function(df_params, df_priors) {
	g <- my_linerange(df_params, df_priors, "Canopy cover") +
		facet_wrap(~method_names) +
		labs(
			x = "Effect",
			y = "State",
			title = "Effect of canopy cover on removal rate",
		) +
		geom_vline(xintercept = 0, linetype = "dashed") +
		theme(legend.position = "none")

	ggsave(
		filename = "plots/states_canopy_cover.tiff",
		plot = g,
		width = 10,
		height = 6,
		units = "in",
		dpi = 300
	)

	g
}

#' Plot road-density effects on removal rate
#'
#' Creates a faceted prior/posterior interval plot for road-density effects by
#' method and state.
#'
#' @param df_params A posterior parameter summary data frame.
#' @param df_priors A prior parameter summary data frame.
#'
#' @return A ggplot object. Also writes `plots/states_road_density.tiff`.
plot_states_road_density <- function(df_params, df_priors) {
	g <- my_linerange(df_params, df_priors, "Road density") +
		facet_wrap(~method_names) +
		labs(
			x = "Effect",
			y = "State",
			title = "Effect of road density on removal rate",
		) +
		geom_vline(xintercept = 0, linetype = "dashed") +
		theme(legend.position = "none")

	ggsave(
		filename = "plots/states_road_density.tiff",
		plot = g,
		width = 10,
		height = 6,
		units = "in",
		dpi = 300
	)

	g
}

#' Plot ruggedness effects on removal rate
#'
#' Creates a faceted prior/posterior interval plot for terrain-ruggedness effects
#' by method and state.
#'
#' @param df_params A posterior parameter summary data frame.
#' @param df_priors A prior parameter summary data frame.
#'
#' @return A ggplot object. Also writes `plots/states_ruggedness.tiff`.
plot_states_ruggedness <- function(df_params, df_priors) {
	g <- my_linerange(df_params, df_priors, "Ruggedness") +
		facet_wrap(~method_names) +
		labs(
			x = "Effect",
			y = "State",
			title = "Effect of terrain ruggedness on removal rate",
		) +
		geom_vline(xintercept = 0, linetype = "dashed") +
		theme(legend.position = "none")

	ggsave(
		filename = "plots/states_ruggedness.tiff",
		plot = g,
		width = 10,
		height = 6,
		units = "in",
		dpi = 300
	)

	g
}

#' Prepare prior and posterior parameter summaries for plotting
#'
#' Combines posterior state-level parameter summaries with prior summaries,
#' assigns plotting groups, and orders states with the prior shown last.
#'
#' @param df_params A posterior parameter summary data frame containing `st_name`,
#'   `position`, and parameter summary columns.
#' @param df_priors A prior parameter summary data frame containing `st_name`,
#'   `position`, `method_idx`, and parameter summary columns.
#'
#' @return A data frame combining prior and posterior summaries with `group` and
#'   ordered `st_name` columns for plotting.
prepare_prior_posterior_params <- function(df_params, df_priors) {
	states <- c(sort(unique(df_params$st_name)), "Prior")

	bind_rows(
		select(df_params, -position),
		select(df_priors, -position, -method_idx),
	) |>
		mutate(
			group = case_when(
				st_name == "Prior" ~ "Prior",
				.default = "Posterior"
			),
			st_name = factor(st_name, levels = states),
		)
}

#' Summarize proportion removed by method
#'
#' Aggregates removals and effort by method, property, and primary period, joins
#' density estimates, converts density to abundance, and averages observed
#' removal proportions and effort metrics by method, property, and state.
#'
#' @param df A raw model data frame containing removal events, effort metrics,
#'   method names, property metadata, and primary periods.
#' @param df_density A density summary data frame with a median density column
#'   named `0.5`.
#'
#' @return A data frame of method-property-state averages for removal proportion,
#'   effort, effort-per-unit, trap count, and event counts.
prop_take_by_method <- function(df, df_density) {
	z_take <- df |>
		fix_method_names() |>
		select(
			propertyID,
			method,
			take,
			effort,
			effort_per,
			trap_count,
			start_dates,
			end_dates,
			st_name,
			order
		) |>
		left_join(df_density) |>
		filter(!is.na(`0.5`)) |>
		mutate(abundance_estimate = round(`0.5` * property_area_km2)) |>
		filter(abundance_estimate > 0) |>
		group_by(
			propertyID,
			st_name,
			start_dates,
			end_dates,
			method,
			property_area_km2,
			abundance_estimate
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
		ungroup() |>
		mutate(
			mu_events = total_events,
			mu_prop_take = total_take / abundance_estimate,
			mu_total_effort = total_effort,
			mu_total_log_effort = total_log_effort,
			mu_total_effort_per = total_effort_per,
			mu_total_log_effort_per = total_log_effort_per,
			mu_total_trap_count = total_trap_count,
			mu_total_log_trap_count = total_log_trap_count
		)

	z_take
}

#' Run management bootstrap simulations
#'
#' Simulates property-level abundance trajectories with and without management by
#' resampling demographic parameters and observed property starting conditions.
#'
#' @param df_params A posterior parameter data frame containing `phi_mu` and
#'   `log_nu` columns.
#' @param df_abundance A posterior abundance-sample data frame with state names
#'   and abundance nodes indexed by `n_id`.
#' @param df_data A prepared model data frame containing abundance estimates,
#'   dates, property metadata, and methods used.
#' @param bootstrap_parameters A list containing `n_ens`, `sim_days`, and `n_sim` elements.
#'
#' @return A data frame comparing simulated density trajectories with and without
#'   management, including density change and percent change from starting
#'   density.
run_bootstrap <- function(
	df_params,
	df_abundance,
	df_data,
	bootstrap_parameters
) {
	process_model <- function(n_init, r, K, n_time, n_ens) {
		N <- matrix(NA, n_ens, n_time)
		N[, 1] <- n_init

		for (j in 2:n_time) {
			nm1 <- N[, j - 1]

			lambda <- nm1 + (nm1 * (r - 1) * (1 - (nm1 / K)))
			lambda <- pmax(0, lambda)

			if (any(is.na(lambda))) {
				print(i)
				stop()
			}

			N[, j] <- round(lambda)
		}
		N
	}

	set.seed(803)

	n_ens <- bootstrap_parameters$n_ens
	sim_days <- bootstrap_parameters$sim_days
	n_sim <- bootstrap_parameters$n_sim

	all_out <- tibble()

	df_data <- df_data |> filter(!is.na(n_id))
	df_params <- df_params |> fix_state_names()

	state_vec <- df_data |>
		pull(st_name) |>
		unique()
	pb <- progress::progress_bar$new(
		format = "  running bootstrap [:bar] :percent eta: :eta",
		total = length(state_vec),
		clear = FALSE,
		width = 60
	)
	for (s in seq_along(state_vec)) {
		state <- state_vec[s]

		state_data <- df_data |>
			filter(st_name == state)

		state_params <- df_params |>
			filter(st_name == state)

		state_abundance <- df_abundance |>
			filter(st_name == state)

		counter <- 1
		while (counter < n_sim) {
			draws <- sample.int(nrow(state_params), n_ens, replace = TRUE)
			phi <- state_params$phi_mu[draws] # survival
			log_nu <- state_params$log_nu[draws] # log per capita growth rate
			zeta <- 28 * exp(log_nu) / 365 # per capita growth rate
			r <- phi + zeta / 2 # finite rate of increase

			i <- sample.int(nrow(state_data), 1)

			tmp_start <- state_data |>
				slice(i) |>
				select(
					propertyID,
					property_area_km2,
					abundance_estimate,
					end_dates,
					year,
					st_name,
					n_id
				)

			start_date <- tmp_start$end_dates
			property <- tmp_start$propertyID
			area <- tmp_start$property_area_km2

			date_seq <- seq.Date(
				start_date,
				to = start_date + sim_days,
				by = "4 weeks"
			)
			node <- paste0("N[", tmp_start$n_id, "]")
			k_dens <- 30 # carrying capacity as a density
			K <- round(k_dens * area) # carrying capacity as abundance

			n_1 <- state_abundance |> slice(draws) |> pull(node)

			n_idx <- which(n_1 > 0 & n_1 <= K)
			if (length(n_idx) == 0) {
				next
			}

			# make sure we use joint samples
			n_1 <- n_1[n_idx]
			r <- r[n_idx]

			tmp_compare <- state_data |>
				filter(propertyID %in% property, end_dates %in% date_seq) |>
				select(
					propertyID,
					property_area_km2,
					abundance_estimate,
					methods_used,
					end_dates,
					year,
					st_name,
					n_id
				)

			n_ids <- tmp_compare$n_id

			timestep_df <- tibble(
				start_date = start_date,
				end_dates = date_seq
			) |>
				filter(end_dates <= max(tmp_compare$end_dates)) |>
				mutate(
					delta_time = as.numeric(end_dates - start_date),
					year = year(end_dates),
					n_id = n_ids[1]:n_ids[length(n_ids)]
				)

			n_time <- nrow(timestep_df)

			if (n_time == 1) {
				next
			}

			N <- process_model(n_1, r, K, n_time, length(r))
			D <- N / area
			colnames(D) <- timestep_df$n_id
			D1 <- n_1 / area

			density_long <- D |>
				as_tibble() |>
				mutate(ens = seq_along(r), D1 = n_1 / area) |>
				pivot_longer(
					cols = -c(ens, D1),
					names_to = "n_id",
					values_to = "density_no_mgmt"
				) |>
				filter(n_id %in% n_ids) |>
				mutate(n_id = as.numeric(n_id))

			out_lut <- left_join(
				tmp_compare,
				timestep_df,
				by = c("end_dates", "year", "n_id")
			)

			out <- left_join(density_long, out_lut, by = "n_id") |>
				filter(delta_time %in% c(364, 728, 1092, 1456, 1820)) |>
				mutate(
					simID = counter,
					density_mgmt = abundance_estimate / area,
					area = area,
					K = K,
					D1 = D1,
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
		}
		pb$tick()
	}

	all_out |>
		pivot_longer(
			cols = c(density_mgmt, density_no_mgmt),
			names_to = "management",
			values_to = "density_estimate"
		) |>
		mutate(
			mgmt = if_else(management == "density_mgmt", "With", "Without"),
			delta_density = density_estimate - D1,
			percent_change = delta_density / D1 * 100,
			unique_sim = paste0(st_name, "-", simID, "-", ens)
		)
}
