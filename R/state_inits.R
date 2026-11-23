#' Function to get MCMC initial values depending on the state being modeled
#'
#' @title state_inits
#'
#' @param STATE_NAME The state being modeled, uppercase (i.e. "FLORIDA")
#' @details if no `STATE_NAME` is given, return NULL for each parameter, initial conditions based
#' on posterior distributions from the original fit in Foster et al. 2026
#'
#' @author John Foster

state_inits <- function(state_name) {
	out <- list()
	out$beta1 <- NULL
	out$beta_p <- NULL
	out$p_mu <- NULL
	out$log_gamma <- NULL
	out$log_rho <- NULL
	out$psi_phi <- NULL
	out$phi_mu <- NULL
	out$log_nu <- NULL

	# Florida ----
	if (state_name == "FLORIDA") {
		out$beta1 <- tribble(
			~min   , ~max   ,
			-2     , -1     ,
			-0.182 ,  0.418 ,
			-7.58  , -5.1   ,
			-0.01  ,  0.01
		)
		out$beta_p <- tribble(
			~min , ~max ,
			-0.1 ,  0.1 , # [1, 1]
			 0.5 ,  2   , # [1, 2]
			-0.1 ,  0.1 , # [1, 3]
			-0.1 ,  0.1 , # [2, 1]
			-0.1 ,  0.1 , # [2, 2]
			-0.1 ,  0.1 , # [2, 3]
			 0.5 ,  1   , # [3, 1]
			-1.5 , -0.5 , # [3, 2]
			 0.5 ,  1.5 , # [3, 3]
			-0.1 ,  0.1 , # [4, 1]
			-0.1 ,  0.1 , # [4, 2]
			-0.1 ,  0.1 # [4, 3]
		)

		out$p_mu <- data.frame(
			min = c(-3, -0.3),
			max = c(-1, 1.93)
		)
		out$log_gamma <- data.frame(
			min = c(-2, -4),
			max = c(-1, -3)
		)
		out$log_rho <- tribble(
			~min  , ~max   ,
			-1    , -0.5   ,
			 1.5  ,  1.7   ,
			-3    , -2     ,
			 0.09 ,  0.341
		)
		out$psi_phi <- c(0.817, 0.976)
		out$phi_mu <- c(0.653, 0.69)
		out$log_nu <- c(2.31, 2.4)
	}

	# Georgia ----
	if (state_name == "GEORGIA") {
		out$beta1 <- tribble(
			~min , ~max  ,
			-2.5 , -1.5  ,
			-1   , -0.75 ,
			-5.5 , -4.5  ,
			-3.1 , -2.9
		)
		out$beta_p <- tribble(
			~min  , ~max  ,
			 0.75 ,  1    , # [1, 1]
			 0.5  ,  1    , # [1, 2]
			-0.75 , -0.25 , # [1, 3]
			-0.2  , -0.1  , # [2, 1]
			-0.75 , -0.25 , # [2, 2]
			-0.75 , -0.25 , # [2, 3]
			-1.25 , -0.5  , # [3, 1]
			-2.5  , -1.5  , # [3, 2]
			 0.25 ,  1    , # [3, 3]
			 0.1  ,  0.4  , # [4, 1]
			 0.1  ,  0.2  , # [4, 2]
			-0.5  , -0.3 # [4, 3]
		)

		out$p_mu <- data.frame(
			min = c(-2, 0.25),
			max = c(-1, 0.75)
		)
		out$log_gamma <- data.frame(
			min = c(-2, -3.75),
			max = c(-1, -3.25)
		)
		out$log_rho <- tribble(
			~min , ~max ,
			 0.1 ,  0.5 ,
			 1.6 ,  1.7 ,
			-3   , -2.5 ,
			 0.7 ,  0.9
		)
		out$psi_phi <- c(0.9, 1.0)
		out$phi_mu <- c(0.65, 0.675)
		out$log_nu <- c(2, 3)
	}

	# Mississippi ----
	if (state_name == "MISSISSIPPI") {
		out$beta1 <- tribble(
			~min , ~max ,
			-3.5 , -2.5 ,
			 1.2 ,  2   ,
			-7   , -5   ,
			-3   , -2.8
		)
		out$beta_p <- tribble(
			~min , ~max ,
			-1   , -0.5 , # [1, 1]
			 0.5 ,  1.5 , # [1, 2]
			 0   ,  0.5 , # [1, 3]
			 0.1 ,  0.3 , # [2, 1]
			 2   ,  3   , # [2, 2]
			-0.5 , -0.1 , # [2, 3]
			-2   , -1   , # [3, 1]
			-3   , -1   , # [3, 2]
			 0   ,  1   , # [3, 3]
			 0.1 ,  0.3 , # [4, 1]
			 1.5 ,  2   , # [4, 2]
			 0   ,  0.2 # [4, 3]
		)

		out$p_mu <- data.frame(
			min = c(-1, -0.1),
			max = c(0.5, 0.1)
		)
		out$log_gamma <- data.frame(
			min = c(-2, -2.5),
			max = c(-1, -2)
		)
		out$log_rho <- tribble(
			~min , ~max ,
			0.5  , 0.8  ,
			0.7  , 0.9  ,
			1.9  , 2    ,
			0.1  , 0.3
		)
		out$psi_phi <- c(0.7, 0.75)
		out$phi_mu <- c(0.625, 0.65)
		out$log_nu <- c(2.58, 2.62)
	}

	# Oklahoma ----
	if (state_name == "OKLAHOMA") {
		out$beta1 <- tribble(
			~min  , ~max  ,
			-2    , -1    ,
			-5.9  , -5.6  ,
			-1    , -0.5  ,
			-3    , -2    ,
			-2.26 , -2.25
		)
		out$beta_p <- tribble(
			~min  , ~max  ,
			-0.3  , -0.1  , # [1, 1]
			-0.01 ,  0.01 , # [1, 2]
			 0.2  ,  0.5  , # [1, 3]
			 0.3  ,  0.5  , # [2, 1]
			-0.01 ,  0.01 , # [2, 2]
			-1.2  , -0.8  , # [2, 3]
			 0.1  ,  0.2  , # [3, 1]
			-0.01 ,  0.01 , # [3, 2]
			-0.2  , -0.1  , # [3, 3]
			-0.01 ,  0.01 , # [4, 1]
			-2.25 , -1.5  , # [4, 2]
			-0.01 ,  0.01 , # [4, 3]
			 0.1  ,  0.2  , # [5, 1]
			 0.1  ,  0.3  , # [5, 2]
			-0.2  , -0.1 # [5, 3]
		)
		out$p_mu <- data.frame(
			min = c(-0.5, -0.01),
			max = c(0.5, 0.01)
		)
		out$log_gamma <- data.frame(
			min = c(-2, -3.75),
			max = c(-1, -3.3)
		)
		out$log_rho <- tribble(
			~min  , ~max  ,
			-0.01 ,  0.01 ,
			 2.6  ,  2.8  ,
			 0.98 ,  2    ,
			-2    , -1.5  ,
			 0.35 ,  0.4
		)

		out$psi_phi <- c(0.7, 0.75)
		out$phi_mu <- c(0.5, 0.6)
		out$log_nu <- c(2.5, 2.60)
	}

	# Texas ----
	if (state_name == "TEXAS") {
		out$beta1 <- tribble(
			~min , ~max ,
			-3   , -2.5 ,
			-2.9 , -2.6 ,
			-0.3 , -0.1 ,
			-3.1 , -2.9 ,
			-2.6 , -2.4
		)
		out$beta_p <- tribble(
			~min  , ~max  ,
			-0.01 ,  0.01 , # [1, 1]
			-1.1  , -0.9  , # [1, 2]
			-0.01 ,  0.01 , # [1, 3]
			-0.01 ,  0.01 , # [2, 1]
			 2.9  ,  3.1  , # [2, 2]
			-0.7  , -0.4  , # [2, 3]
			 0.9  ,  1.1  , # [3, 1]
			 0.9  ,  1.1  , # [3, 2]
			-0.01 ,  0.01 , # [3, 3]
			 0.1  ,  0.3  , # [4, 1]
			-1.75 , -1.4  , # [4, 2]
			 0.4  ,  0.6  , # [4, 3]
			-0.01 ,  0.01 , # [5, 1]
			 1.9  ,  2.1  , # [5, 2]
			-0.6  , -0.4 # [5, 3]
		)
		out$p_mu <- data.frame(
			min = c(0.8, -3),
			max = c(1.1, -2.6)
		)
		out$log_gamma <- data.frame(
			min = c(-2.5, -3.5),
			max = c(-2, -3)
		)
		out$log_rho <- tribble(
			~min , ~max ,
			-0.6 , -0.4 ,
			 2.7 ,  2.9 ,
			 1.8 ,  1.9 ,
			-2.1 , -1.9 ,
			 0.3 ,  0.4
		)

		out$psi_phi <- c(0.7, 0.9)
		out$phi_mu <- c(0.66, 0.68)
		out$log_nu <- c(2.1, 2.3)
	}

	# Georgia ----
	if (state_name == "GEORGIA") {
		out$beta1 <- tribble(
			~min  , ~max   ,
			-2.25 , -1.750 ,
			-1.2  , -0.8   ,
			-5.5  , -4.5   ,
			-3.1  , -2.9
		)
		out$beta_p <- tribble(
			~min  , ~max  ,
			 0.5  ,  1    , # [1, 1]
			 0.5  ,  1    , # [1, 2]
			-0.8  , -0.2  , # [1, 3]
			-0.25 , -0.1  , # [2, 1]
			-0.75 , -0.5  , # [2, 2]
			-0.75 , -0.25 , # [2, 3]
			-1.5  , -0.5  , # [3, 1]
			-2.25 , -1.75 , # [3, 2]
			 0.25 ,  0.75 , # [3, 3]
			 0.1  ,  0.2  , # [4, 1]
			 0.39 ,  0.41 , # [4, 2]
			-0.5  , -0.25 # [4, 3]
		)

		out$p_mu <- data.frame(
			min = c(-2, 0.5),
			max = c(-0.5, 1.5)
		)
		out$log_gamma <- data.frame(
			min = c(-2, -3.6),
			max = c(-1, -3.25)
		)
		out$log_rho <- tribble(
			~min  , ~max  ,
			-0.01 ,  0.01 ,
			 1.7  ,  1.8  ,
			-3    , -2.5  ,
			 0.8  ,  0.9
		)
		out$psi_phi <- c(0.8, 0.9)
		out$phi_mu <- c(0.65, 0.7)
		out$log_nu <- c(2, 2.4)
	}

	# Ohio ----
	if (state_name == "OHIO") {
		out$beta1 <- tribble(
			~min , ~max ,
			-1   ,  1   ,
			 0.5 ,  1.5 ,
			-5   , -3   ,
			-2.2 , -2
		)
		out$beta_p <- tribble(
			~min  , ~max  ,
			-1    , -1    , # [1, 1]
			 0    ,  2    , # [1, 2]
			 0    ,  2    , # [1, 3]
			-0.75 , -0.25 , # [2, 1]
			 1    ,  2    , # [2, 2]
			-0.4  ,  0.4  , # [2, 3]
			-2    , -1    , # [3, 1]
			-3    , -1    , # [3, 2]
			 0.5  ,  2    , # [3, 3]
			 0.8  ,  1.2  , # [4, 1]
			 3.5  ,  3.7  , # [4, 2]
			-0.25 , -0.1 # [4, 3]
		)
		out$p_mu <- data.frame(
			min = c(-1, 1),
			max = c(1, 3)
		)
		out$log_gamma <- data.frame(
			min = c(-2, 0.5),
			max = c(-1, 1)
		)
		out$log_rho <- tribble(
			~min , ~max ,
			-1   ,  1   ,
			 1.6 ,  2   ,
			-2.5 , -1   ,
			-0.2 , -0.1
		)

		out$psi_phi <- c(0.8, 1)
		out$phi_mu <- c(0.66, 0.72)
		out$log_nu <- c(1.9, 2)
	}

	# Louisiana ----
	if (state_name == "LOUISIANA") {
		out$beta1 <- tribble(
			~min , ~max ,
			-1   ,  0   ,
			 0.1 ,  1   ,
			-4   , -3   ,
			-2   , -0.5
		)
		out$beta_p <- tribble(
			~min  , ~max  ,
			 0.5  ,  0.75 , # [1, 1]
			 1.75 ,  2.25 , # [1, 2]
			-0.1  ,  0.1  , # [1, 3]
			-0.01 ,  0.01 , # [2, 1]
			-0.01 ,  0.01 , # [2, 2]
			-0.4  , -0.2  , # [2, 3]
			-0.5  , -0.1  , # [3, 1]
			-2    , -1    , # [3, 2]
			 1.5  ,  2.5  , # [3, 3]
			 0.6  ,  0.7  , # [4, 1]
			 1.45 ,  1.55 , # [4, 2]
			-0.01 ,  0.01 # [4, 3]
		)
		out$p_mu <- data.frame(
			min = c(-1, -2.25),
			max = c(0.5, -1.75)
		)
		out$log_gamma <- data.frame(
			min = c(-2, -2.2),
			max = c(-1, -2)
		)
		out$log_rho <- tribble(
			~min , ~max ,
			-0.5 , -0.1 ,
			 0.8 ,  0.1 ,
			 0.2 ,  0.5 ,
			-0.1 ,  0.1
		)

		out$psi_phi <- c(0.8, 1)
		out$phi_mu <- c(0.65, 0.675)
		out$log_nu <- c(2.35, 2.4)
	}
	# Tennessee ----
	if (state_name == "TENNESSEE") {
		out$beta1 <- tribble(
			~min , ~max  ,
			-2   ,  0    ,
			-1   , -0.5  ,
			-4   , -2    ,
			-1.5 , -1.25
		)
		out$beta_p <- tribble(
			~min  , ~max ,
			-2    ,  0   , # [1, 1]
			 0.5  ,  1.5 , # [1, 2]
			-2    , -1   , # [1, 3]
			-1    , -0.5 , # [2, 1]
			-0.3  , -0.2 , # [2, 2]
			-0.4  ,  0.2 , # [2, 3]
			-2    , -1   , # [3, 1]
			-3    , -1   , # [3, 2]
			 1    ,  2   , # [3, 3]
			-0.75 , -0.5 , # [4, 1]
			 3    ,  3.5 , # [4, 2]
			 0.1  ,  0.5 # [4, 3]
		)
		out$p_mu <- data.frame(
			min = c(-1, -3.5),
			max = c(0.5, -3)
		)
		out$log_gamma <- data.frame(
			min = c(-2, -3.5),
			max = c(-1, -3)
		)
		out$log_rho <- tribble(
			~min  , ~max  ,
			-1    ,  0    ,
			 1.75 ,  2    ,
			-2.5  , -1.5  ,
			-1.3  , -1.25
		)

		out$psi_phi <- c(0.8, 0.9)
		out$phi_mu <- c(0.65, 0.675)
		out$log_nu <- c(2.35, 2.4)
	}

	# Indiana ----
	if (state_name == "INDIANA") {
		out$beta1 <- tribble(
			~min , ~max ,
			   0 ,  1   ,
			  -3 , -1   ,
			  -3 , -2.5
		)
		out$beta_p <- tribble(
			~min , ~max ,
			-0.8 , 0    , # [1, 1]
			 0   , 2    , # [1, 2]
			-0.2 , 0.2  , # [1, 3]
			-1.4 , 0    , # [2, 1]
			-2   , 0    , # [2, 2]
			 1   , 2    , # [2, 3]
			-0.5 , 0    , # [3, 1]
			 0.5 , 1.5  , # [3, 2]
			-0.5 , 0 # [3, 3]
		)
		out$p_mu <- data.frame(
			min = c(-1, 1),
			max = c(0.5, 3)
		)
		out$log_gamma <- data.frame(
			min = c(-2, -2.4),
			max = c(-1, -2)
		)
		out$log_rho <- tribble(
			~min , ~max ,
			 1.7 ,  1.9 ,
			-2   , -1   ,
			 0.1 ,  0.3
		)

		out$psi_phi <- c(0.8, 1)
		out$phi_mu <- c(0.65, 0.7)
		out$log_nu <- c(1.8, 2)
	}
	# Kentucky ----
	if (state_name == "KENTUCKY") {
		out$beta1 <- tribble(
			~min , ~max ,
			-0.5 ,  0.5 ,
			-3   , -1   ,
			-2.5 , -2
		)
		out$beta_p <- tribble(
			~min , ~max ,
			-0.8 ,  0   , # [1, 1]
			 0   ,  1   , # [1, 2]
			-0.2 ,  0.2 , # [1, 3]
			-2   , -1   , # [2, 1]
			-3   , -1   , # [2, 2]
			 1   ,  2   , # [2, 3]
			 0.5 ,  1   , # [3, 1]
			 3   ,  4   , # [3, 2]
			-0.4 ,  0 # [3, 3]
		)
		out$p_mu <- data.frame(
			min = c(-1, -1),
			max = c(1, 1)
		)
		out$log_gamma <- data.frame(
			min = c(-2, -1.4),
			max = c(-1, -0.6)
		)
		out$log_rho <- tribble(
			~min , ~max ,
			 1.6 ,  1.8 ,
			-2.5 , -1.5 ,
			-0.7 , -0.5
		)

		out$psi_phi <- c(0.8, 1)
		out$phi_mu <- c(0.65, 0.7)
		out$log_nu <- c(1.8, 2)
	}
	# West Virginia ----
	if (state_name == "WEST VIRGINIA") {
		out$beta1 <- tribble(
			~min , ~max ,
			-0.5 ,  0.5 ,
			-4   , -2   ,
			-2.5 , -2
		)
		out$beta_p <- tribble(
			~min , ~max  ,
			-2   , -1    , # [1, 1]
			-0.5 ,  0    , # [1, 2]
			 0.5 ,  2    , # [1, 3]
			-1.5 , -0.5  , # [2, 1]
			-2   , -0.5  , # [2, 2]
			 1   ,  2    , # [2, 3]
			 0.5 ,  0.75 , # [3, 1]
			-0.1 ,  0.1  , # [3, 2]
			-0.4 , -0.1 # [3, 3]
		)
		out$p_mu <- data.frame(
			min = c(-1, 1),
			max = c(1, 3)
		)
		out$log_gamma <- data.frame(
			min = c(-2, -3),
			max = c(-1, -2)
		)
		out$log_rho <- tribble(
			~min , ~max  ,
			-0.5 ,  0.5  ,
			-2.5 , -1.5  ,
			-0.2 ,  0.25
		)

		out$psi_phi <- c(0.8, 1)
		out$phi_mu <- c(0.65, 0.7)
		out$log_nu <- c(1.9, 2)
	}
	out
}

# Template if all methods used
# out$beta1 <- tribble(
# 	~min , ~max ,
# 	  NA ,   NA ,
# 	  NA ,   NA ,
# 	  NA ,   NA ,
# 	  NA ,   NA
# )
# out$beta_p <- tribble(
# 	~min , ~max ,
# 	NA   , NA   , # [1, 1]
# 	NA   , NA   , # [1, 2]
# 	NA   , NA   , # [1, 3]
# 	NA   , NA   , # [2, 1]
# 	NA   , NA   , # [2, 2]
# 	NA   , NA   , # [2, 3]
# 	NA   , NA   , # [3, 1]
# 	NA   , NA   , # [3, 2]
# 	NA   , NA   , # [3, 3]
# 	NA   , NA   , # [4, 1]
# 	NA   , NA   , # [4, 2]
# 	NA   , NA   , # [4, 3]
# 	NA   , NA   , # [5, 1]
# 	NA   , NA   , # [5, 2]
# 	NA   , NA     # [5, 3]
# )

# out$p_mu <- data.frame(
# 	min = c(NA, NA),
# 	max = c(NA, NA)
# )
# out$log_gamma <- data.frame(
# 	min = c(NA, NA),
# 	max = c(NA, NA)
# )
# out$log_rho <- tribble(
# 	~min , ~max ,
# 	NA   ,    NA ,
# 	NA   ,    NA ,
# 	NA   ,    NA ,
# 	NA   ,    NA
# )
# out$psi_phi <- c(NA, NA)
# out$phi_mu <- c(NA, NA)
# out$log_nu <- c(NA, NA)
