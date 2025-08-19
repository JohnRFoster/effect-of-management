library(targets)
library(ggplot2)

g <- tar_read("plot_management")

ggsave(
	filename = "plots/management_trajectory.jpeg",
	plot = g,
	width = 10,
	height = 6,
	units = "in",
	dpi = 300
)

g <- tar_read("plot_lambda_hat_by_state")

ggsave(
	filename = "plots/lambda_hat_by_state.jpeg",
	plot = g,
	width = 10,
	height = 6,
	units = "in",
	dpi = 300
)

g <- tar_read("plot_prop_take")

ggsave(
	filename = "plots/prop_take.jpeg",
	plot = g,
	width = 10,
	height = 6,
	units = "in",
	dpi = 300
)

g <- tar_read("plot_epr")

ggsave(
	filename = "plots/epr_effort_compAdditive.jpeg",
	plot = g,
	width = 12,
	height = 10,
	units = "in",
	dpi = 300
)
