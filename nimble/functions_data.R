## functions called throughout the target pipeline

## ---------------------------- Data ingest -----------------------------------

create_primary_periods <- function(df, interval, data_repo) {
	require(lubridate)

	end_dates <- unique(sort(df$end.date))
	min_date <- min(end_dates)
	max_date <- max(end_dates)

	start_dates <- seq(min_date, max_date, by = paste(interval, "week"))
	end_dates <- c(start_dates[-1] - 1, max_date)

	targets::tar_assert_identical(length(start_dates), length(end_dates))
	targets::tar_assert_true(min(df$start.date) >= min_date)
	targets::tar_assert_true(max(df$start.date) <= max_date)

	timestep_df <- tibble(start_dates, end_dates) |>
		mutate(timestep = 1:n())
	timestep_df$month <- month(timestep_df$end_dates)
	timestep_df$year <- year(timestep_df$end_dates)

	# for each row in the merged data, insert the integer primary period timestep
	df$timestep <- NA
	message("Assign timesteps...")
	pb <- txtProgressBar(max = nrow(df), style = 1)
	for (i in 1:nrow(df)) {
		after_start <- which(timestep_df$start_dates <= df$start.date[i]) |> max()
		before_end <- which(timestep_df$end_dates >= df$end.date[i]) |> min()
		if (after_start == before_end) {
			# then the start and end date is contained within a primary period
			df$timestep[i] <- timestep_df$timestep[before_end]
		} # otherwise, timestep[i] will be left as NA and filtered out later
		setTxtProgressBar(pb, i)
	}
	close(pb)

	write_csv(timestep_df, file.path(data_repo, "calibration_timestep_df.csv"))

	df |>
		filter(!is.na(timestep)) |>
		left_join(timestep_df) |>
		arrange(propertyID, timestep)
}

# resolve duplicate property values - when there are multiple values, take max
resolve_duplicate <- function(insitu_data) {
	property_areas <- insitu_data |>
		distinct(propertyID, property_area_km2) |>
		group_by(propertyID) |>
		summarize(
			n_areas = length(unique(property_area_km2)),
			property_area_max = max(property_area_km2, na.rm = TRUE),
			# properties with all NA areas get -Inf
			# but the following line changes -Inf to NA
			property_area_km2 = ifelse(
				is.infinite(property_area_max),
				NA,
				property_area_km2
			)
		) |>
		ungroup()

	insitu_data |>
		left_join(property_areas, by = join_by(propertyID, property_area_km2)) |>
		filter(!is.na(property_area_km2), property_area_km2 >= 1.8, effort > 0)
}

# need to filter events so that there are at least two events in a timestep
# and there are at least 2 timesteps for each property
dynamic_filter <- function(df) {
	good_events <- df |>
		select(propertyID, timestep) |>
		group_by(propertyID, timestep) |>
		mutate(two_plus_takes = n() >= 2) |>
		filter(two_plus_takes) |>
		group_by(propertyID) |>
		arrange(propertyID, timestep) |>
		distinct() |>
		mutate(n_timesteps = length(unique(timestep))) |>
		filter(n_timesteps >= 2) |>
		ungroup() |>
		mutate(event_id = paste0(propertyID, "-", timestep)) |>
		pull(event_id)

	df |>
		mutate(event_id = paste0(propertyID, "-", timestep)) |>
		filter(event_id %in% good_events) |>
		arrange(propertyID, timestep)
}

take_filter <- function(df) {
	zero_take_prp <- df |>
		group_by(propertyID) |>
		summarise(sum_take = sum(take)) |>
		ungroup() |>
		filter(sum_take == 0) |>
		pull(propertyID)

	df |> filter(!propertyID %in% zero_take_prp)
}

# compute ordering based on time interval midpoints
order_interval <- function(df) {
	df |>
		distinct() |>
		mutate(midpoint = as.numeric(end.date)) |>
		ungroup()
}

# impose stochastic ordering of events by adding jitter
# we are assuming the following order of events when the events have the same midpoint
# e.g., are on the same day:
# 1. (trap or snare), with order random
# 2. (heli or plane), with order random
# 3. hunting
order_stochastic <- function(order.df) {
	n_methods_mid <- order.df |>
		select(propertyID, midpoint, method) |>
		group_by(propertyID, midpoint) |>
		count() |>
		arrange(propertyID, midpoint)

	set.seed(8)

	order_df <- order.df
	order_df$jittered_midpoint <- NA
	message("Stochastic ordering...")
	pb <- txtProgressBar(max = nrow(order_df), style = 1)
	for (i in 1:nrow(order_df)) {
		if (order_df$method[i] %in% c('Trap', 'Snare')) {
			order_df$jittered_midpoint[i] <- order_df$midpoint[i] +
				runif(1, min = 0, max = .01)
		} else if (order_df$method[i] %in% c('Fixed Wing', 'Helicopter')) {
			order_df$jittered_midpoint[i] <- order_df$midpoint[i] +
				runif(1, min = .02, max = .03)
		} else {
			order_df$jittered_midpoint[i] <- order_df$midpoint[i] +
				runif(1, min = .04, max = .05)
		}
		setTxtProgressBar(pb, i)
	}
	return(order_df)
}

# now compute orders of events based on jittered midpoints
order_of_events <- function(order_df) {
	df <- order_df |>
		ungroup() |>
		group_by(propertyID, timestep) |>
		mutate(
			order = order(jittered_midpoint),
			has_multi = any(order > 1),
			any_ties = any(duplicated(jittered_midpoint)),
			n_survey = n()
		) |>
		ungroup() |>
		arrange(propertyID, timestep, order) |>
		mutate(p = 1:n())

	targets::tar_assert_true(all(df$has_multi))
	targets::tar_assert_true(!all(df$any_ties))
	targets::tar_assert_true(all(df$n_survey > 1))

	return(df)
}

county_codes <- function(df) {
	df |>
		rename(statefp = st_gsa_state_cd, countyfp = cnty_gsa_cnty_cd) |>
		mutate(
			countyfp = sprintf("%03d", countyfp),
			countyfp = ifelse(cnty_name == "HUMBOLDT (E)", "013", countyfp),
			county_code = as.numeric(paste0(statefp, countyfp)),
			county_code = sprintf("%05d", county_code)
		) |> # need this for joining downstream
		select(
			propertyID,
			agrp_prp_id,
			alws_agrprop_id,
			start_dates,
			end_dates,
			st_name,
			cnty_name,
			county_code,
			method,
			trap_count,
			take,
			property_area_km2,
			effort,
			effort_per,
			timestep,
			order,
			n_survey,
			p
		) |>
		rename(primary_period = timestep)
}

get_fips <- function(file = "data/fips/national_county.txt") {
	fips <- read_csv(
		file,
		show_col_types = FALSE,
		col_names = c("state", "statefp", "countyfp", "countyname", "classfp"),
		comment = '#'
	)
}

get_ts_length <- function(df) {
	df |>
		filter(!st_name %in% c("CALIFORNIA", "ALABAMA", "ARIZONA", "ARKANSAS")) |>
		select(propertyID, primary_period) |>
		distinct() |>
		group_by(propertyID) |>
		filter(
			primary_period == min(primary_period) |
				primary_period == max(primary_period)
		) |>
		mutate(delta = c(0, diff(primary_period) + 1)) |>
		ungroup() |>
		filter(delta != 0)
}

get_n_observations <- function(df, good_props) {
	df |>
		filter(propertyID %in% good_props) |>
		group_by(propertyID, primary_period) |>
		summarise(take = sum(take)) |>
		ungroup() |>
		# filter(take > 0) |>
		group_by(propertyID) |>
		count()
}

subset_data_for_development <- function(df) {
	require(rsample)

	observed_pp <- df |>
		select(propertyID, primary_period) |>
		distinct() |>
		count(propertyID) |>
		rename(n_observed_pp = n)

	ts_length <- get_ts_length(df) |>
		rename(ts_length = delta) |>
		select(-primary_period)

	best_properties <- df |>
		group_by(propertyID, property_area_km2) |>
		summarise(
			take = sum(take),
			effort = sum(effort_per),
			unit_count = sum(trap_count),
			n_total_events = n()
		) |>
		ungroup() |>
		left_join(ts_length) |>
		left_join(observed_pp) |>
		mutate(prop_observed = n_observed_pp / ts_length) |>
		filter(
			property_area_km2 >= 4.05,
			property_area_km2 <= 405,
			take >= 13,
			take <= 858,
			n_total_events >= 7,
			n_total_events <= 117,
			n_observed_pp >= 2,
			n_observed_pp <= 18,
			ts_length >= 3,
			ts_length <= 34,
			prop_observed >= 0.24,
			effort >= 13.9,
			effort <= 263,
			unit_count >= 10,
			unit_count <= 653
		)

	property_order <- best_properties$propertyID

	# create strata by decile
	# each property will belong to a decile of each land cover variable
	df_strata <- df |>
		mutate(
			canopy_strata = make_strata(c_canopy, breaks = 10),
			rugged_strata = make_strata(c_rugged, breaks = 10),
			road_den_strata = make_strata(c_road_den, breaks = 10)
		) |>
		select(propertyID, contains("strata")) |>
		distinct()

	df_sample <- df_strata |>
		filter(propertyID %in% property_order)

	message("Number of properties in each canopy strata:")
	print(table(df_sample$canopy_strata))
	message("Number of properties in each rugged strata:")
	print(table(df_sample$rugged_strata))
	message("Number of properties in each road density strata:")
	print(table(df_sample$road_den_strata))

	new_data <- df |>
		filter(propertyID %in% property_order)
	return(new_data)
}

condition_first_capture <- function(df) {
	good_events <- df |>
		group_by(propertyID, timestep) |>
		summarise(take = sum(take)) |>
		ungroup() |>
		group_by(propertyID) |>
		mutate(cumulative_take = cumsum(take)) |>
		ungroup() |>
		filter(cumulative_take > 0) |>
		mutate(event_id = paste0(propertyID, "-", timestep)) |>
		pull(event_id)

	df |>
		mutate(event_id = paste0(propertyID, "-", timestep)) |>
		filter(event_id %in% good_events) |>
		arrange(propertyID, timestep)
}

create_timestep_df <- function(df) {
	df |>
		select(propertyID, primary_period) |>
		unique() |>
		arrange(propertyID, primary_period) |>
		group_by(propertyID) |>
		mutate(observed_timestep = 1:n()) |> # timestep is the sequence of primary periods within a property
		ungroup() |>
		mutate(primary_period = primary_period - min(primary_period) + 1)
}


get_data <- function(file, interval, data_repo) {
	all_take <- read_csv(file, show_col_types = FALSE) |>
		filter(start.date >= lubridate::ymd("2014-01-01")) |>
		mutate(
			cnty_name = if_else(
				grepl("ST ", cnty_name),
				gsub("ST ", "ST. ", cnty_name),
				cnty_name
			),
			cnty_name = if_else(grepl("KERN", cnty_name), "KERN", cnty_name)
		)

	data_mis <- all_take |>
		mutate(property_area_km2 = round(property.size * 0.00404686, 2)) |>
		filter(property_area_km2 >= 1.8, st_name != "HAWAII") |>
		mutate(
			effort = if_else(
				cmp_name %in% c("TRAPS, CAGE", "SNARE"),
				cmp.days,
				cmp.hours
			),
			effort_per = effort / cmp.qty,
			cmp_name = if_else(cmp_name == "TRAPS, CAGE", "TRAPS", cmp_name)
		) |>
		rename(method = cmp_name, trap_count = cmp.qty) |>
		select(-wt_work_date, -hours, -cmp.hours, -cmp.days) |>
		distinct() |>
		mutate(propertyID = paste0(agrp_prp_id, "-", alws_agrprop_id)) |>
		arrange(propertyID, start.date, end.date)

	# create PP of length [interval]
	data_timestep <- create_primary_periods(data_mis, interval, data_repo) |>
		resolve_duplicate() |> # resolve duplicate property areas
		take_filter() |> # remove properties with zero pigs taken
		dynamic_filter() |> # filter out bad events & properties
		condition_first_capture() |> # condition on first positive removal event for each property
		dynamic_filter() |> # filter out bad events & properties
		order_interval() |> # determine midpoints from start/end dates
		order_stochastic() |> # randomly order events
		order_of_events() |> # assign order number, check
		county_codes() # county codes and renaming

	# now we have two columns for time
	# primary_period is how [interval] sequences are aligned across the data set
	# timestep is the sequence of primary periods within a property
	timestep_df <- create_timestep_df(data_timestep)

	data_pp <- left_join(
		data_timestep,
		timestep_df,
		by = join_by(propertyID, primary_period)
	) |>
		mutate(primary_period = primary_period - min(primary_period) + 1)

	return(data_pp)
}

# for missing values, impute with state mean
mean_impute <- function(x) {
	ifelse(is.na(x), mean(x, na.rm = TRUE), x)
}

# generate centered and scaled versions of these numeric variables
center_scale <- function(x) {
	(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Deal with observation covariates, a small percentage of which are missing
get_obs_covars <- function(file) {
	obs_covs <- file |>
		read_csv(show_col_types = FALSE) |>
		mutate(county_code = sprintf("%05d", FIPS)) |>
		dplyr::select(-starts_with('sd'), -NAME, -FIPS)

	obs_covs <- obs_covs |>
		group_by(STATE_NAME) |>
		mutate(
			rural.road.density = mean_impute(rural.road.density),
			prop.pub.land = mean_impute(prop.pub.land),
			mean.ruggedness = mean_impute(mean.ruggedness),
			mean.canopy.density = mean_impute(mean.canopy.density)
		) |>
		ungroup() |>
		mutate(
			c_road_den = center_scale(rural.road.density),
			c_rugged = center_scale(mean.ruggedness),
			c_canopy = center_scale(mean.canopy.density)
		)

	targets::tar_assert_true(!any(is.na(obs_covs$c_road_den)))
	targets::tar_assert_true(!any(is.na(obs_covs$c_rugged)))
	targets::tar_assert_true(!any(is.na(obs_covs$c_canopy)))

	return(obs_covs)
}
