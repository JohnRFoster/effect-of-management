# download and process county and ecoregion shapefiles

library(dplyr)

counties <- tigris::counties(
	cb = TRUE,
	year = 2024,
	class = "sf"
) |>
	select(
		county_code = GEOID,
		county_name = NAME
	) |>
	sf::st_make_valid() |>
	sf::st_transform(5070) |>
	mutate(
		county_area_km2 = as.numeric(sf::st_area(geometry)) / 1e6
	)

ecoregion_data <- file.path(
	data_store,
	"NorthAmericaEcoregionsLevel2",
	"NA_CEC_Eco_Level2.shp"
)
my_shapefile <- sf::st_read(ecoregion_data)

ecoregions <- my_shapefile |>
	select(
		ecoregion_code = NA_L2CODE,
		ecoregion_name = NA_L2NAME
	) |>
	sf::st_make_valid() |>
	sf::st_transform(5070)

county_ecoregion <- sf::st_intersection(
	counties,
	ecoregions
) |>
	mutate(
		overlap_km2 = as.numeric(sf::st_area(geometry)) / 1e6,
		overlap_proportion = overlap_km2 / county_area_km2
	) |>
	sf::st_drop_geometry() |>
	group_by(county_code) |>
	slice_max(overlap_km2, n = 1, with_ties = FALSE) |>
	ungroup() |>
	select(
		county_code,
		county_name,
		ecoregion_code,
		ecoregion_name,
		overlap_proportion
	)

readr::write_csv(county_ecoregion, "../data-store/county_ecoregion.csv")
