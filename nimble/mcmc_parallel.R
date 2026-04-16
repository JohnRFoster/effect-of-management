#---------
#
# Workhorse function for running nimble in parallel
#
#---------

#'@description Run nimble on specified cluster, check for convergence and save samples periodically
#'@param cl Cluster made using the parallel package
#'@param model_code the nimble model code
#'@param model_constants model constants, as a list
#'@param model_data model data, as a list
#'@param model_inits initial model values, as a list, one list element for each chain
#'@param params_check vector of parameters (nodes) to assess convergence
#'@param n_iters the number of iterations to run in each chunk after the model has compiled
#'@param dest output filepath. within this directory each 'chunk' will be saved into it's own directory
#'@param monitors_add vector of nodes to monitor in addition to default nodes
#'@param custom_samplers data frame of mcmc config. Nodes in one column, sampler types in another

mcmc_parallel <- function(
  cl,
  model_code,
  model_constants,
  model_data,
  model_inits,
  params_check,
  n_iters,
  dest,
  monitors_add = NULL,
  custom_samplers = NULL
) {
  require(foreach)
  require(coda)
  require(doParallel)

  # add function arguments to the local environment, needed for exporting to 'cl'
  as.list(environment(), all = TRUE)

  # compile model
  # this function gets passed as an argument and runs on each cluster of 'cl'
  # runs through typical nimble set-up precidure
  # can add nodes to monitors and change samplers from the default configuration
  single_mcmc_chain <- function() {
    require(nimble)
    source("nimble/functions_nimble.R")

    Rmodel <- nimbleModel(
      code = model_code,
      constants = model_constants,
      data = model_data,
      inits = init,
      calculate = TRUE
    )

    Rmodel$initializeInfo()

    N <- Rmodel$N
    nH_p <- model_constants$nH_p
    n_survey <- model_constants$n_survey
    y_sum <- model_data$y_sum

    for (i in 1:n_survey) {
      N_model <- N[nH_p[i]]
      n <- round(N_model - y_sum[i])
      if (n <= 0) {
        print(i)
        Rmodel$N[nH_p[i]] <- N_model + (abs(n)^2)
      }
    }

    # calc <- Rmodel$calculate()
    # if(is.infinite(calc) | is.nan(calc) | is.na(calc)){
    #   stop(paste0("Model log probability is ", calc))
    # }

    # default MCMC configuration
    mcmcConf <- configureMCMC(Rmodel, useConjugacy = TRUE)

    control_rw <- list(
      adaptInterval = 100,
      adaptFactorExponent = 0.6
    )

    mcmcConf$removeSamplers("beta1")
    mcmcConf$addSampler(
      target = "beta1",
      type = "RW_block",
      control = control_rw
    )

    mcmcConf$removeSamplers("log_nu")
    mcmcConf$addSampler(
      target = "log_nu",
      type = "slice"
    )

    # if specified, change nodes to specified parameters
    if (!is.null(custom_samplers)) {
      for (i in seq_len(nrow(custom_samplers))) {
        node <- custom_samplers$node[i]
        type <- custom_samplers$type[i]
        mcmcConf$removeSampler(node)
        mcmcConf$addSampler(node, type)
      }
    }

    if (!is.null(monitors_add)) {
      mcmcConf$addMonitors(monitors_add)
    }

    Rmcmc <- buildMCMC(mcmcConf)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <<- compileNimble(Rmcmc) # need to define this in the global environment to continue mcmc
    Cmcmc$run(niter = 1000, nburnin = 0, thin = 10)
    samples <- as.matrix(Cmcmc$mvSamples)
    return(samples)
  }

  # sample for the specified duration 'n_iters'
  # drop previous samples (we save them after each time this function is called, so we don't need them)
  continue_sampling <- function() {
    require(nimble)
    require(coda)
    Cmcmc$run(niter = n_iters, reset = FALSE, resetMV = TRUE, thin = 10)
    samples <<- as.matrix(Cmcmc$mvSamples) # need to define this in the global environment to use in subset_mcmc
    return(samples)
  }

  # for taking nodes and splitting them from an mcmc list
  # as is, does one node at a time but can be passed to lapply to use on a vector of nodes
  subset_mcmc <- function(node) {
    require(dplyr)
    s <- samples[, grep(node, colnames(samples), value = TRUE, fixed = TRUE)] |>
      as.matrix() |>
      as_tibble()
    if (ncol(s) == 1) {
      colnames(s) <- node
    }
    return(s)
  }

  # get parameter nodes using subset_mcmc
  subset_params <- function() {
    require(dplyr)
    require(purrr)
    map_dfc(lapply(params_check, subset_mcmc), as_tibble) |> as.matrix()
  }

  # get observed abundance nodes using subset_mcmc
  subset_N_observed <- function() {
    nodes <- paste0("N[", model_constants$N_full_unique, "]")
    map_dfc(lapply(nodes, subset_mcmc), as_tibble) |> as.matrix()
  }

  # get unobserved abundance nodes using subset_mcmc
  subset_N_unobserved <- function() {
    nodes <- paste0("N[", model_constants$N_quant_unique, "]")
    N_unobserved <- map_dfc(lapply(nodes, subset_mcmc), as_tibble) |>
      as.matrix()
    t(apply(N_unobserved, 2, quantile, c(0.025, 0.5, 0.975)))
  }

  export <- c(
    "single_mcmc_chain",
    "continue_sampling",
    "subset_mcmc",
    "subset_params",
    "subset_N_observed",
    "subset_N_unobserved",
    "model_code",
    "model_data",
    "model_constants",
    "n_iters",
    "custom_samplers",
    "monitors_add",
    "params_check"
  )

  clusterExport(cl, export, envir = environment())

  for (i in seq_along(cl)) {
    init <- model_inits[[i]]
    clusterExport(cl[i], "init", envir = environment())
  }

  # initialize model and first samples
  c <- 1
  start <- Sys.time()
  out <- clusterEvalQ(cl, single_mcmc_chain())
  message("Model compile and initial 1000 iterations completed in:")
  print(round(Sys.time() - start, 2))

  start2 <- Sys.time()
  out2 <- clusterEvalQ(cl, continue_sampling())
  message("Additional ", n_iters, " iterations completed in:")
  print(round(Sys.time() - start2, 2))

  message("\n", n_iters * c, " total iterations completed in:")
  print(round(Sys.time() - start, 2))

  # use mcmc on clusters to subset parameters, observed states, and unobserved states
  params <- clusterEvalQ(cl, subset_params())
  params <- as.mcmc.list(lapply(params, as.mcmc))

  N_observed <- clusterEvalQ(cl, subset_N_observed())
  N_observed <- as.mcmc.list(lapply(N_observed, as.mcmc))

  N_unobserved <- clusterEvalQ(cl, subset_N_unobserved()) |> as.matrix()

  params_mcmc_list <- as.mcmc.list(lapply(params, as.mcmc))
  efsize <- 1000
  diagnostic <- continue_mcmc(
    params_mcmc_list,
    effective_size = efsize,
    max_psrf = 15,
    verbose = TRUE
  )

  c_dir <- sprintf("%04d", c)
  path <- file.path(dest, c_dir)
  if (!dir.exists(path)) {
    dir.create(path, showWarnings = FALSE, recursive = TRUE)
  }

  write_out_p <- function(p, d, dest) {
    f <- "paramSamples.rds"
    write_rds(list(params = p, diagnostic = d), file.path(dest, f))
  }

  write_abundnace <- function(no, nu, dest) {
    f <- "observedAbundanceSamples.rds"
    write_rds(no, file.path(dest, f))

    f <- "unobservedAbundanceQuantiles.rds"
    write_rds(nu, file.path(dest, f))
  }

  write_out_p(params, diagnostic, path)

  converged <- all(diagnostic$psrf[, 2] < 1.1)
  if (converged) {
    write_abundnace(N_observed, N_unobserved, path)
  }

  continue <- !diagnostic$done
  write_N <- FALSE
  while (continue) {
    c <- c + 1
    resetMV <- TRUE
    clusterExport(cl, "resetMV", envir = environment())

    start2 <- Sys.time()
    out2 <- clusterEvalQ(cl, continue_sampling())
    message("Additional ", n_iters, " iterations completed in:")
    print(round(Sys.time() - start2, 2))

    total_iters <- n_iters * c
    message("\n", total_iters, " total iterations completed in:")
    print(round(Sys.time() - start, 2))

    # use mcmc on clusters to subset parameters, observed states, and unobserved states
    params <- clusterEvalQ(cl, subset_params())
    params <- as.mcmc.list(lapply(params, as.mcmc))

    c_dir <- sprintf("%04d", c)
    path <- file.path(dest, c_dir)
    if (!dir.exists(path)) {
      dir.create(path, showWarnings = FALSE, recursive = TRUE)
    }

    diagnostic <- NULL
    write_out_p(params, diagnostic, path)

    N_observed <- clusterEvalQ(cl, subset_N_observed())
    N_observed <- as.mcmc.list(lapply(N_observed, as.mcmc))
    N_unobserved <- clusterEvalQ(cl, subset_N_unobserved()) |> as.matrix()

    params_mcmc_list <- collate_mcmc_chunks(dest)$params
    diagnostic <- continue_mcmc(
      params_mcmc_list,
      effective_size = efsize,
      max_psrf = 15,
      verbose = TRUE
    )
    converged <- all(diagnostic$psrf[, 2] < 1.1)
    write_N <- if_else(converged, TRUE, FALSE)

    if (converged | write_N) {
      write_abundnace(N_observed, N_unobserved, path)
    }

    continue <- !diagnostic$done
    message("=================================================")

    if (c == 100) continue <- FALSE
  }

  if (diagnostic$done) {
    return(TRUE)
  }
  if (!diagnostic$done) return(FALSE)
}
