

get_next_property <- function(all_data, data_last_fit){
  require(rsample)

  # find properties not in first fit
  last_fit_properties <- unique(data_last_fit$propertyID)
  all_properties <- unique(all_data$propertyID)
  not_fit_properties <- setdiff(all_properties, last_fit_properties)
  targets::tar_assert_true(
    length(not_fit_properties) + length(last_fit_properties) == length(all_properties)
  )

  ts_length <- get_ts_length(all_data)
  n_obs <- get_n_observations(all_data, ts_length$propertyID)

  ts_sorted <- all_data |>
    group_by(propertyID, property_area_km2) |>
    summarise(take = sum(take),
              effort = sum(effort_per),
              unit_count = sum(trap_count),
              n_total_events = n()) |>
    ungroup() |>
    left_join(ts_length) |>
    left_join(n_obs)

  group1 <- ts_sorted |>
    filter(property_area_km2 <= 2339,
           take > 1, take <= 16443,
           n_total_events >= 4, n_total_events <= 269,
           n >= 2, n <= 35,
           delta >= 2, delta <= 40,
           effort >= 1.3, effort <= 614,
           unit_count >= 4, unit_count <= 1947) |>
    pull(propertyID)

  get_strata <- function(df){
    df_strata <- df |>
      mutate(canopy_strata = as.numeric(as.factor(make_strata(c_canopy, breaks = 10))),
             rugged_strata = as.numeric(as.factor(make_strata(c_rugged, breaks = 10))),
             road_den_strata = as.numeric(as.factor(make_strata(c_road_den, breaks = 10)))) |>
      select(propertyID, st_name, cnty_name, contains("strata")) |>
      distinct()

  }

  strata1 <- all_data |>
    filter(propertyID %in% group1) |>
    get_strata()
  data1 <- ts_sorted |>
    filter(propertyID %in% group1)

  strata2 <- all_data |>
    filter(!propertyID %in% group1) |>
    get_strata()
  data2 <- ts_sorted |>
    filter(!propertyID %in% group1)

  get_order <- function(df, strat){
    df_join <- left_join(df, strat)

    varmax <- df_join |>
      group_by(rugged_strata, road_den_strata) |>
      summarise(gmax = max(canopy_strata)) |>
      ungroup()

    # this is the order of properties if we didn't fit the initial batch
    property_order <- left_join(df_join, varmax) |>
      arrange(desc(gmax))
  }

  # rank remaining properties
  order1 <- get_order(data1, strata1)
  order2 <- get_order(data2, strata2)

  # remaining properties
  remaining_property_order <- bind_rows(order1, order2) |>
    filter(!propertyID %in% last_fit_properties)

  if(nrow(remaining_property_order) == 0){
    stop("All properties have been fit!")
  } else {
    next_property <- remaining_property_order |>
      slice(1) |>
      pull(propertyID)

    new_data <- all_data |> filter(propertyID == next_property)

    to_fit_properties <- unique(new_data$propertyID)
    not_fit_properties <- setdiff(all_properties, to_fit_properties)
    targets::tar_assert_true(
      length(not_fit_properties) + length(to_fit_properties) == length(all_properties)
    )

    return(bind_rows(data_last_fit, new_data))
  }

}
