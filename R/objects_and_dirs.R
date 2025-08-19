# project directories
plot_dir <- "plots"
data_dir <- "data"
last_iter <- "11_posterior"
first_iter <- "1_posterior"

# shared data across projects live here
data_store <- "../data-store"

# cuttoff date for analysis
# 2023 is the last full year with complete data
cutoff_date <- lubridate::ymd("2023-12-31")

# the number of bootstrap samples to draw
n_ens <- 1000

# fmt: skip
state_names <- dplyr::tribble(
	~st_name,         ~state,
	"OKLAHOMA",       "Oklahoma",
	"TEXAS",          "Texas",
	"NEW MEXICO",     "New Mexico",
	"LOUISIANA",      "Louisiana",
	"VIRGINIA",       "Virginia",
	"GEORGIA",        "Georgia",
	"MISSISSIPPI",    "Mississippi",
	"MISSOURI",       "Missouri",
	"NORTH CAROLINA", "North Carolina",
	"SOUTH CAROLINA", "South Carolina",
	"FLORIDA",        "Florida",
	"WEST VIRGINIA",  "West Virginia"
)
