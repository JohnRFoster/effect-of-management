## model and helper functions for nimble
library(nimble)


# calculate the potential area (on the log scale) for a survey given method
calc_log_potential_area <- nimbleFunction(
  run = function(
    log_rho = double(1),
    log_gamma = double(1),
    p_unique = double(1),
    log_effort_per = double(0),
    effort_per = double(0),
    n_trap_m1 = double(0),
    log_pi = double(0),
    method = double(0)
  ){
    m <- method

    if(m <= 3){ # firearms, fixed wing, and helicopter
      log_potential_area <- log_rho[m] + log_effort_per
    } else { # traps and snares
      log_potential_area <- log_pi +
        (2 * (log_rho[m] + log_effort_per -
                log(exp(log_gamma[m-3]) + effort_per))) +
        log(1 + (p_unique[m-3] * n_trap_m1))
    }
    return(log_potential_area)
    returnType(double(0))
  },
  buildDerivs = TRUE
)


# figure out if we need to keep sampling
continue_mcmc <- function(mcmc, effective_size, max_psrf, verbose){
  require(coda)
  require(purrr)

  message("Checking convergence and sample size")
  psrf <- gelman.diag(mcmc, multivariate = FALSE)$psrf
  effective_samples <- effectiveSize(mcmc)

  converged <- all(psrf[, 2] < 1.1)
  enough_samples <- all(effective_samples >= effective_size)
  funky <- any(is.nan(psrf)) | max(psrf) > max_psrf

  message('Convergence [', converged, ']')
  message('Enough effective samples [', enough_samples, ']')
  message('Mixing [', !funky, ']')

  done <- converged & enough_samples
  if(done) message("MCMC complete!")

  # TODO determine burnin and effective sample size POST burnin

  if(funky){
    message("\n*** Something is wrong with the mcmc! ***")
    done <- FALSE
  }

  if(verbose){
    print(psrf)
    print(effective_samples)
  }

  return(list(
    done = done,
    psrf = psrf
  ))

}

test_build <- function(code, constants, data, inits){

  message("\n\n=== Test build ===")

  Rmodel <- nimbleModel(
    code = code,
    constants = constants,
    data = data,
    inits = inits
  )

  Rmodel$initializeInfo()

  N <- Rmodel$N
  nH_p <- constants$nH_p
  n_survey <- constants$n_survey
  y_sum <- data$y_sum

  for(i in 1:n_survey){
    N_model <- N[nH_p[i]]
    n <- round(N_model - y_sum[i])
    if(n <= 0){
      print(i)
      Rmodel$N[nH_p[i]] <- N_model + (abs(n) + 200)
    }
  }

  # Rmodel$simulate()

  calc <- Rmodel$calculate()
  if(is.infinite(calc) | is.nan(calc) | is.na(calc)){
    stop(paste0("Model log probability is ", calc))
  }

  message("==================\n")
}

collate_mcmc_chunks <- function(dest, start = 1){
  mcmc_dirs <- list.files(dest)
  mcmc_dirs <- setdiff(mcmc_dirs, "modelData.rds")
  mcmc_dirs <- mcmc_dirs[start:length(mcmc_dirs)]
  param_file_name <- "paramSamples.rds"
  state_file_name <- "observedAbundanceSamples.rds"

  # use the first mcmc chunk to initialize storage for each chain
  mcmc_rds <- file.path(dest, mcmc_dirs[1], param_file_name)
  rds <- read_rds(mcmc_rds)
  mcmc <- rds$params
  n_chains <- length(mcmc)
  store_mcmc <- list()
  store_mcmc_state <- list()
  state_count <- 0

  # store each chain, will append below
  for(j in seq_len(n_chains)){
    store_mcmc[[j]] <- as.matrix(mcmc[[j]])
  }

  state_rds <- file.path(dest, mcmc_dirs[1], state_file_name)
  if(file.exists(state_rds)){
    state_count <- 1
    mcmc <- read_rds(state_rds)
    for(j in seq_len(n_chains)){
      store_mcmc_state[[j]] <- as.matrix(mcmc[[j]])
    }
  }

  if(length(mcmc_dirs) >= 2){

    use_pb <- if_else(length(mcmc_dirs) == 2, FALSE, TRUE)
    if(use_pb){
      pb <- txtProgressBar(min = 2, max = length(mcmc_dirs), style = 1)
    }

    # read each mcmc chunk, store each chain from the chunk as a matrix
    for(i in 2:length(mcmc_dirs)){
      mcmc_rds <- file.path(dest, mcmc_dirs[i], param_file_name)
      rds <- read_rds(mcmc_rds)
      mcmc <- rds$params
      for(j in seq_len(n_chains)){
        store_mcmc[[j]] <- rbind(store_mcmc[[j]], as.matrix(mcmc[[j]]))
      }

      state_rds <- file.path(dest, mcmc_dirs[i], state_file_name)
      if(file.exists(state_rds)){

        mcmc <- read_rds(state_rds)

        if(state_count == 0){
          for(j in seq_len(n_chains)){
            store_mcmc_state[[j]] <- as.matrix(mcmc[[j]])
          }
        } else {
          for(j in seq_len(n_chains)){
            store_mcmc_state[[j]] <- rbind(store_mcmc_state[[j]], as.matrix(mcmc[[j]]))
          }
        }
        state_count <- 1
      }
      if(use_pb) setTxtProgressBar(pb, i)
    }
    if(use_pb) close(pb)
  }


  # need to create an mcmc list object to check for convergence
  params <- as.mcmc.list(lapply(store_mcmc, as.mcmc))

  if(state_count == 1){
    states <- as.mcmc.list(lapply(purrr::compact(store_mcmc_state), as.mcmc))
  } else {
    states <- list()
  }


  return(list(params = params, states = states))

}
