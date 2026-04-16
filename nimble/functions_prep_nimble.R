## functions called throughout the target pipeline

# ==========================================
# Helper functions ----
# ==========================================

# need all timesteps whether there are observations or not
create_all_primary_periods <- function(df) {
  pp_min_max <- df |>
    select(property, primary_period) |>
    distinct() |>
    group_by(property) |>
    filter(
      primary_period == min(primary_period) |
        primary_period == max(primary_period)
    ) |>
    ungroup()

  properties <- unique(pp_min_max$property)

  all_pp <- tibble()
  message("\nInclude all primary periods")
  pb <- txtProgressBar(max = length(properties), style = 1)
  for (i in seq_along(properties)) {
    pid <- pp_min_max |> filter(property == properties[i])
    p_min <- min(pid$primary_period)
    p_max <- max(pid$primary_period)
    pp <- tibble(
      property = properties[i],
      primary_period = p_min:p_max,
      timestep = 1:length(p_min:p_max)
    )
    all_pp <- bind_rows(all_pp, pp)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  all_pp |> mutate(n_id = 1:n())
}

# need to know the total number of timesteps in each property (sampled or not) for indexing
n_timesteps <- function(df) {
  df |>
    group_by(property) |>
    filter(timestep == max(timestep)) |>
    pull(timestep)
}

# index (as a matrix) for tracking abundance in long format (converting wide to long)
# used in the process model
N_lookup_table <- function(df) {
  df |>
    mutate(n_id = 1:n()) |>
    select(-primary_period) |>
    pivot_wider(names_from = timestep, values_from = n_id) |>
    select(-property) |>
    as.matrix()
}

# calculate the cumulative number of pigs taken as a primary period progresses
removed_in_pp_cumsum <- function(df) {
  df |>
    group_by(property, primary_period) |>
    mutate(ysum = cumsum(take) - take) |>
    ungroup() |>
    pull(ysum)
}

# the total number of pigs taken in a primary period across all methods
# including periods without removals effort (equal to 0)
# wide format
total_take <- function(df_take, df_pp) {
  sum_take <- df_take |>
    group_by(property, primary_period) |>
    summarise(sum_take = sum(take)) |>
    ungroup()

  left_join(df_pp, sum_take, by = join_by(property, primary_period)) |>
    mutate(sum_take = if_else(is.na(sum_take), 0, sum_take)) |>
    select(-primary_period, -n_id) |>
    pivot_wider(names_from = timestep, values_from = sum_take) |>
    select(-property) |>
    as.matrix()
}

# index (long format) for which primary periods have removal effort
# and therefore included in the data model
N_lookup_data <- function(df_take, df_pp) {
  tH <- df_take |>
    select(property, primary_period)

  df_pp |>
    select(property, primary_period, n_id) |>
    right_join(tH, by = join_by(property, primary_period)) |>
    pull(n_id)
}

# need start and end indicies for data model
# used for estimiating p
create_start_end <- function(df_take, df_pp) {
  start <- end <- numeric(nrow(df_take))

  df <- left_join(df_take, df_pp, by = join_by(primary_period, property))

  message("Creating start/end indicies")
  pb <- txtProgressBar(max = nrow(df), style = 1)
  for (i in seq_len(nrow(df))) {
    if (df$order[i] > 1) {
      idx <- which(
        df$county == df$county[i] &
          df$property == df$property[i] &
          df$timestep == df$timestep[i] &
          df$order < df$order[i]
      )
      start[i] <- idx[1]
      end[i] <- idx[length(idx)]
      assertthat::are_equal(idx, start[i]:end[i])
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  tibble(start = start, end = end)
}

# informed hyper parameters for beta distribution on global pig survival
create_surv_prior <- function(interval, data_repo) {
  require(lubridate)
  require(readr)
  require(dplyr)
  require(tidyr)

  data <- read_csv(
    file.path(data_repo, "insitu/Vital_Rate_Data.csv"),
    show_col_types = FALSE
  )

  data_usa <- data |>
    filter(
      country == "USA",
      time.period.end != "null",
      time.period.start != "null",
      !paper.ID %in% c(128, 1007, 130, 136)
    ) |> # these papers don't have specified date ranges or are meta-analysis
    mutate(
      time.period.end = mdy(time.period.end),
      time.period.start = mdy(time.period.start)
    )

  surv_data <- data_usa |>
    filter(!is.na(survival.prop)) |>
    select(
      unique.ID,
      paper.ID,
      N.hogs.in.study,
      contains("survival"),
      contains("hunting"),
      state,
      contains("time"),
      method.for.data
    )

  surv_mu <- surv_data |>
    mutate(
      weeks = as.numeric(time.period.end - time.period.start) / 7,
      weeks4 = weeks / interval,
      survival.per.4week = survival.prop^(1 / weeks4),
      logit.survival.per.4week = boot::logit(survival.per.4week)
    ) |>
    filter(survival.per.4week > 0) |>
    mutate(scale_factor = survival.per.4week / survival.prop)

  surv_mu_summary <- surv_mu |>
    summarise(
      mu = mean(survival.per.4week),
      mu.logit = mean(logit.survival.per.4week)
    )

  surv_var <- surv_data |>
    filter(survival.var.type %in% c("SD", "95% CI"))

  surv_sd <- surv_var |>
    filter(survival.var.type == "SD") |>
    mutate(sd = as.numeric(survival.var))

  surv_sd_calc <- surv_var |>
    filter(survival.var.type == "95% CI") |>
    mutate(
      low.CI = as.numeric(stringr::str_extract(
        survival.var,
        "[[:graph:]]*(?=\\-)"
      )),
      high.CI = as.numeric(stringr::str_extract(
        survival.var,
        "(?<=\\-)[[:graph:]]*"
      )),
      sd_low = (low.CI - survival.prop) / -1.96,
      sd_high = (high.CI - survival.prop) / 1.96
    ) |>
    group_by(unique.ID) |>
    summarise(sd = max(sd_high, sd_low))

  surv_var_join <- left_join(surv_var, surv_sd_calc, by = join_by(unique.ID)) |>
    filter(survival.var.type != "SD")

  scale_ids <- surv_mu |>
    select(unique.ID, scale_factor)

  surv_variance <- bind_rows(surv_var_join, surv_sd) |>
    left_join(scale_ids, by = join_by(unique.ID)) |>
    mutate(
      variance = sd^2,
      variance.4week = variance * scale_factor^2,
      sd.4week = sqrt(variance.4week)
    )

  surv_sd_summary <- surv_variance |>
    pull(sd.4week) |>
    mean()

  mu <- surv_mu_summary$mu
  psi <- 1 / mean(surv_variance$variance.4week)
  alpha <- mu * psi
  beta <- (1 - mu) * psi

  return(list(
    alpha = alpha,
    beta = beta
  ))
}

# matrix of landscape covariates for data model
create_X <- function(df, cols = c("c_road_den", "c_rugged", "c_canopy")) {
  df |>
    select(all_of(cols)) |>
    as.matrix()
}

get_prior_hyperparams <- function(
  informed,
  posterior_path = NULL,
  data_repo = NULL
) {
  if (informed) {
    rds <- file.path(posterior_path, "posteriorSamples.rds")
    post <- read_rds(rds)

    get_vec <- function(df, node) {
      df |>
        select(contains(node)) |>
        pivot_longer(cols = everything(), names_to = "node") |>
        group_by(node) |>
        summarise(mu = mean(value), tau = 1 / var(value))
    }

    log_rho <- get_vec(post, "log_rho")
    p_mu <- get_vec(post, "p_mu")
    log_gamma <- get_vec(post, "log_gamma")
    beta1 <- get_vec(post, "beta1")
    beta_p <- get_vec(post, "beta_p")
    log_nu <- get_vec(post, "log_nu")

    phi_mu <- get_vec(post, "phi_mu")
    mu <- phi_mu$mu
    v <- 1 / phi_mu$tau
    w <- ((mu * (1 - mu)) / v) - 1
    alpha <- mu * w
    beta <- (1 - mu) * w

    psi_phi <- get_vec(post, "psi_phi")
    mu <- psi_phi$mu
    v <- 1 / psi_phi$tau
    shape <- (mu^2) / v
    rate <- mu / v

    hyperparams <- list(
      log_rho_mu = log_rho$mu,
      log_rho_tau = log_rho$tau,
      p_mu_mu = p_mu$mu,
      p_mu_tau = p_mu$tau,
      log_gamma_mu = log_gamma$mu,
      log_gamma_tau = log_gamma$tau,
      beta1_mu = beta1$mu,
      beta1_tau = beta1$tau,
      beta_p_mu = beta_p$mu,
      beta_p_tau = beta_p$tau,
      phi_mu_a = alpha,
      phi_mu_b = beta,
      psi_shape = shape,
      psi_rate = rate,
      log_nu_mu = log_nu$mu,
      log_nu_tau = log_nu$tau
    )
  } else {
    survival_prior <- create_surv_prior(4, data_repo)

    hyperparams <- list(
      log_rho_mu = rep(0, 5),
      log_rho_tau = c(2, 1, 1, 3, 3),
      p_mu_mu = rep(0, 2),
      p_mu_tau = rep(1, 2),
      log_gamma_mu = rep(0, 2),
      log_gamma_tau = rep(3, 2),
      beta1_mu = rep(0, 5),
      beta1_tau = rep(1, 5),
      beta_p_mu = rep(0, 15),
      beta_p_tau = rep(1, 15),
      phi_mu_a = survival_prior$alpha,
      phi_mu_b = survival_prior$beta,
      psi_shape = 1,
      psi_rate = 0.1,
      log_nu_mu = 2,
      log_nu_tau = 1
    )
  }
  return(hyperparams)
}


# ==========================================
# Create lists for nimble ----
# ==========================================

nimble_constants <- function(
  df,
  interval,
  data_repo,
  informed,
  posterior_path = NULL
) {
  # mean litter size year from VerCauteren et al. 2019 pg 63
  data_ls <- round(
    c(
      5.6,
      6.1,
      5.6,
      6.1,
      4.2,
      5.0,
      5.0,
      6.5,
      5.5,
      6.8,
      5.6,
      5.9,
      4.9,
      5.1,
      4.5,
      4.7,
      5.3,
      5.7,
      7.4,
      8.4,
      4.7,
      4.9,
      3.0,
      3.0,
      4.8,
      4.8,
      4.2,
      5.4,
      4.7,
      5.2,
      5.4
    )
  )

  all_primary_periods <- create_all_primary_periods(df)
  n_time_prop <- n_timesteps(all_primary_periods)
  n_method <- length(unique(df$method))
  nH <- N_lookup_table(all_primary_periods)
  nH_p <- N_lookup_data(df, all_primary_periods)
  N_full_unique <- nH_p |> unique()
  N_quant_unique <- setdiff(seq(1, max(N_full_unique)), N_full_unique)
  rem <- total_take(df, all_primary_periods)
  targets::tar_assert_true(all(rem[, 1] != 0))
  X <- create_X(df)
  start_end <- create_start_end(df, all_primary_periods)

  constants <- list(
    n_survey = nrow(df),
    n_ls = length(data_ls),
    n_property = length(unique(df$property)),
    n_first_survey = length(which(df$order == 1)),
    n_not_first_survey = length(which(df$order != 1)),
    n_method = n_method,
    n_betaP = n_method * ncol(X),
    beta_p_row = rep(1:n_method, each = ncol(X)),
    beta_p_col = rep(1:ncol(X), n_method),
    n_time_prop = n_time_prop,
    nH = nH,
    nH_p = nH_p,
    N_full_unique = N_full_unique,
    N_quant_unique = N_quant_unique,
    rem = rem,
    log_pi = log(pi),
    first_survey = which(df$order == 1),
    not_first_survey = which(df$order != 1),
    m_p = ncol(X),
    start = start_end$start,
    end = start_end$end,
    method = as.numeric(as.factor(df$method)),
    pp_len = interval * 7
  )

  prior_hyperparams <- get_prior_hyperparams(
    informed = informed,
    posterior_path = posterior_path,
    data_repo = data_repo
  )

  append(constants, prior_hyperparams)
}

nimble_data <- function(df) {
  # mean litter size year from VerCauteren et al. 2019 pg 63
  data_ls <- round(
    c(
      5.6,
      6.1,
      5.6,
      6.1,
      4.2,
      5.0,
      5.0,
      6.5,
      5.5,
      6.8,
      5.6,
      5.9,
      4.9,
      5.1,
      4.5,
      4.7,
      5.3,
      5.7,
      7.4,
      8.4,
      4.7,
      4.9,
      3.0,
      3.0,
      4.8,
      4.8,
      4.2,
      5.4,
      4.7,
      5.2,
      5.4
    )
  )

  X <- create_X(df)
  y_sum <- removed_in_pp_cumsum(df)

  # need these for inits
  unique_log_areas <- df |>
    select(propertyID, property_area_km2) |>
    distinct() |>
    pull(property_area_km2) |>
    log()

  list(
    y = df$take,
    y_sum = y_sum,
    J = data_ls,
    X_p = X,
    effort_per = df$effort_per,
    log_effort_per = log(df$effort_per),
    n_trap_m1 = df$trap_count - 1,
    log_survey_area_km2 = log(df$property_area_km2),
    unique_log_areas = unique_log_areas
  )
}


# ==========================================
# Create inits for nimble ----
# ==========================================

nimble_inits <- function(constants_nimble, data_nimble, buffer = 1000) {
  with(append(constants_nimble, data_nimble), {
    beta1 <- jitter(c(-1, -3.75, -0.25, 0.3, -1.5))
    beta_p <- matrix(
      jitter(
        c(1.75, 1.5, 0, 0, -0.5, -1.2, 0.15, 0, 0.1, -1.75, -1, 0, -0.75, 0, 0)
      ),
      n_method,
      m_p
    )
    p_mu <- jitter(c(-4, -3))
    log_gamma <- jitter(c(-3, -2.1))
    log_rho <- jitter(c(0.8, 2.25, 2.15, -1.35, -0.55))
    psi_phi <- runif(1, 0.65, 0.7)
    phi_mu <- runif(1, 0.57, 0.59)
    mean_ls <- jitter(12.5)

    a <- phi_mu * psi_phi
    b <- (1 - phi_mu) * psi_phi
    mean_lpy <- 1
    zeta <- mean_lpy / 365 * pp_len * mean_ls
    N <- phi <- lambda <- rep(NA, max(nH, na.rm = TRUE))
    n_init <- rep(NA, n_property)
    for (i in 1:n_property) {
      n_init[i] <- round(exp(unique_log_areas[i])) + sum(rem[i, ], na.rm = TRUE)
      # if(n_init[i] > 5000) n_init[i] <- rpois(1, 500)
      N[nH[i, 1]] <- n_init[i]
      for (j in 2:n_time_prop[i]) {
        phi[nH[i, j - 1]] <- max(0.05, min(rbeta(1, a, b), 0.95))
        z <- N[nH[i, j - 1]] - rem[i, j - 1]
        z <- max(1, z)
        lambda[nH[i, j - 1]] <- z * zeta / 2 + z * phi[nH[i, j - 1]]

        N[nH[i, j]] <- rpois(1, lambda[nH[i, j - 1]])
      }
    }

    list(
      log_lambda_1 = log(n_init + buffer),
      beta_p = beta_p,
      beta1 = beta1,
      p_mu = p_mu,
      p_unique = boot::inv.logit(p_mu),
      phi_mu = phi_mu,
      psi_phi = psi_phi,
      a_phi = a,
      b_phi = b,
      N = N + buffer,
      # lambda = lambda + buffer,
      log_nu = log(mean_ls),
      nu = mean_ls,
      log_gamma = log_gamma,
      log_rho = log_rho,
      phi = phi,
      zeta = zeta,
      log_zeta = log(zeta)
    )
  })
}

nimble_inits_sample <- function(
  posterior_file,
  constants_nimble,
  data_nimble,
  buffer = 1000
) {
  params <- read_rds(posterior_file)
  draw <- sample.int(nrow(params), 1)
  init <- params |> slice(draw)

  get_vec <- function(df, node) {
    df |>
      select(contains(node)) |>
      pivot_longer(cols = everything()) |>
      pull(value)
  }

  beta_p_vec <- get_vec(init, "beta_p")
  beta_p <- matrix(beta_p_vec, 5, 3)

  beta1 <- get_vec(init, "beta1")
  p_mu <- get_vec(init, "p_mu")
  phi_mu <- get_vec(init, "phi_mu")
  psi_phi <- get_vec(init, "psi_phi") + 2
  log_nu <- get_vec(init, "log_nu")
  log_gamma <- get_vec(init, "log_gamma")
  log_rho <- get_vec(init, "log_rho")

  unique_log_areas <- data_nimble$unique_log_areas
  rem <- constants_nimble$rem
  n_property <- constants_nimble$n_property
  nH <- constants_nimble$nH
  n_time_prop <- constants_nimble$n_time_prop
  pp_len <- constants_nimble$pp_len

  a <- phi_mu * psi_phi
  b <- (1 - phi_mu) * psi_phi
  mean_lpy <- 1
  mean_ls <- exp(log_nu)
  zeta <- mean_lpy / 365 * pp_len * mean_ls
  N <- phi <- lambda <- rep(NA, max(nH, na.rm = TRUE))
  n_init <- rep(NA, n_property)
  for (i in 1:n_property) {
    n_init[i] <- round(exp(unique_log_areas[i])) + sum(rem[i, ], na.rm = TRUE)
    # if(n_init[i] > 5000) n_init[i] <- rpois(1, 500)
    N[nH[i, 1]] <- n_init[i]
    for (j in 2:n_time_prop[i]) {
      phi[nH[i, j - 1]] <- max(0.05, min(rbeta(1, a, b), 0.95))
      z <- N[nH[i, j - 1]] - rem[i, j - 1]
      z <- max(1, z)
      lambda[nH[i, j - 1]] <- z * zeta / 2 + z * phi[nH[i, j - 1]]

      N[nH[i, j]] <- rpois(1, lambda[nH[i, j - 1]])
    }
  }

  list(
    log_lambda_1 = log(n_init + buffer),
    beta_p = beta_p,
    beta1 = beta1,
    p_mu = p_mu,
    p_unique = boot::inv.logit(p_mu),
    phi_mu = phi_mu,
    psi_phi = psi_phi,
    a_phi = a,
    b_phi = b,
    N = N + buffer,
    # lambda = lambda + buffer,
    log_nu = log_nu,
    nu = exp(log_nu),
    log_gamma = log_gamma,
    log_rho = log_rho,
    phi = phi,
    zeta = zeta,
    log_zeta = log(zeta)
  )
}
