# miscelanueous functions

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
