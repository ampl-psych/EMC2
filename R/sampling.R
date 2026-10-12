pmwgs <- function(dadm, type, pars = NULL, prior = NULL, ...) {
  if(is.data.frame(dadm)) dadm <- list(dadm)
  if(is.null(pars)) pars <- names(sampled_pars(attr(dadm[[1]], "prior")))
  if(is.null(prior)) prior <- attr(dadm[[1]], "prior")
  dadm <- extractDadms(dadm)

  dadm_list <-dadm$dadm_list
  # Storage for the samples.
  subjects <- sort(as.numeric(unique(dadm$subjects)))
  samples <- sample_store(dadm, pars, type, ...)
  sampler <- list(
    data = dadm_list,
    par_names = pars,
    subjects = subjects,
    n_pars = length(pars),
    n_subjects = length(subjects),
    samples = samples,
    type = type,
    init = FALSE
  )
  class(sampler) <- "pmwgs"
  sampler <- add_info(sampler, prior, type, ...)
  return(sampler)
}

#' Initialize Chains
#'
#' Adds a set of start points to each chain. These start points are sampled from a user-defined multivariate
#' normal across subjects.
#'
#' @param emc An emc object made by `make_emc()`
#' @param start_mu A vector. Mean of multivariate normal used in proposal distribution
#' @param start_var A matrix. Variance covariance matrix of multivariate normal used in proposal distribution.
#' Smaller values will lead to less deviation around the mean.
#' @param cores_per_chain An integer. How many cores to use per chain. Parallelizes across participant calculations.
#' @param cores_for_chains An integer. How many cores to use to parallelize across chains. Default is the number of chains.
#' @param particles An integer. Number of starting values
#' @param ... optional additional arguments
#'
#' @return An emc object
#' @examples \donttest{
#' # Make a design and an emc object
#' design_DDMaE <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1),
#'                            constants=c(s=log(1)))
#'
#' DDMaE <- make_emc(forstmann, design_DDMaE, compress = FALSE)
#' # set up our mean starting points (same used across subjects).
#' mu <- c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
#'        t0=log(.2))
#' # Small variances to simulate start points from a tight range
#' var <- diag(0.05, length(mu))
#' # Initialize chains, 4 cores per chain, and parallelizing across our 3 chains as well
#' # so 4*3 cores used.
#' DDMaE <- init_chains(DDMaE, start_mu = mu, start_var = var,
#'                      cores_per_chain = 1, cores_for_chains = 1, particles = 3)
#' # Afterwards we can just use fit
#' # DDMaE <- fit(DDMaE, cores_per_chain = 4)
#' }
#' @export
init_chains <- function(emc, start_mu = NULL, start_var = NULL, particles = 1000,
                        cores_per_chain=1,cores_for_chains = length(emc),
                        ...)
{
  emc <- restore_duplicates(emc)
  emc <- lapply(emc, function(x) if (x$type == "SEM") normalise_SEM_prior(x) else x)
  context <- sampler_contexts(emc, cores_per_chain * min(length(emc), cores_for_chains))
  emc <- native_init_chains(emc, context, start_mu, start_var, particles, TRUE)
  class(emc) <- "emc"
  return(emc)
}

check_tune_settings <- function(tune, n_pars, stage, particles){
  if(is.null(tune$search_width)) tune$search_width <- 1
  # Acceptance ratio tuning
  tune$alphaStar <- ifelse(stage == "sample", 2, 3)
  tune$p_accept <- set_p_accept(stage, tune$search_width)
  # Potential blocking settings
  if(is.null(tune$components)) tune$components <- rep(1, n_pars)
  if(is.null(tune$shared_ll_idx)) tune$shared_ll_idx <- tune$components
  # Tuning of number of particles, might be a bit arbitrary
  if(is.null(tune$target_ESS)) tune$target_ESS <- 2.5*sqrt(n_pars)
  if(is.null(tune$ESS_scale)) tune$ESS_scale <- .05
  if(is.null(tune$max_particles)) tune$max_particles <- particles*1.2
  # Mix tuning settings
  if(is.null(tune$mix_adapt)) tune$mix_adapt <- .05
  # After n0 all the tuning kicks in
  tune$n0 <- 25
  return(tune)
}

check_sampling_settings <- function(pm_settings, stage, n_pars, particles){
  for(i in 1:length(pm_settings)){
    # Mix settings
    pm_settings[[i]]$mix <- check_mix(pm_settings[[i]]$mix, stage)
    # For p_accept
    pm_settings[[i]]$epsilon <- check_epsilon(pm_settings[[i]]$epsilon, n_pars, pm_settings[[i]]$mix)
    # Only local components have an adapted epsilon.
    pm_settings[[i]]$epsilon[!local_components(stage)[-1]] <- 1
    # For mix and p_accept tuning
    pm_settings[[i]]$proposal_counts <- check_prop_performance(pm_settings[[i]]$proposal_counts, stage)
    pm_settings[[i]]$acc_counts <- check_prop_performance(pm_settings[[i]]$acc_counts, stage)
    # Setting particles
    if(is.null(pm_settings[[i]]$n_particles)) pm_settings[[i]]$n_particles <- particles
    if(is.null(pm_settings[[i]]$iter)) pm_settings[[i]]$iter <- 1
  }
  return(pm_settings)
}

# Retry numerically singular group draws.
resolve_on_singular <- function(on_singular) {
  defaults <- list(max_retries = 3)
  if (is.null(on_singular)) return(defaults)
  if (!is.list(on_singular)) stop("`on_singular` must be a list or NULL")
  bad <- setdiff(names(on_singular), names(defaults))
  if (length(bad)) stop("unknown `on_singular` field(s): ", paste(bad, collapse = ", "),
                        ". Valid fields: ", paste(names(defaults), collapse = ", "))
  modifyList(defaults, on_singular)
}

# Components centred on the current value.
local_components <- function(stage){
  switch(stage,
         preburn = c(FALSE, TRUE),
         burn = c(FALSE, TRUE, TRUE),
         adapt = c(FALSE, TRUE, FALSE),
         c(FALSE, TRUE, FALSE, FALSE))
}

extend_sampler <- function(sampler, n_samples, stage) {
  # This function takes the sampler and extends it along the intended number of
  # iterations, to ensure that we're not constantly increasing our sampled object
  # by 1. Big shout out to the rapply function
  sampler$samples$stage <- c(sampler$samples$stage, rep(stage, n_samples))
  last_theta_var_inv <- sampler$samples$last_theta_var_inv
  sampler$samples$last_theta_var_inv <- NULL
  sampler$samples <- rapply(sampler$samples, f = function(x) extend_obj(x, n_samples), how = "replace")
  sampler$samples$last_theta_var_inv <- last_theta_var_inv
  return(sampler)
}

extend_obj <- function(obj, n_extend){
  old_dim <- dim(obj)
  n_dimensions <- length(old_dim)
  if(is.null(old_dim) | n_dimensions == 1) return(obj)
  new_dim <- c(rep(0, (n_dimensions -1)), n_extend)
  extended <- array(NA_real_, dim = old_dim +  new_dim, dimnames = dimnames(obj))
  extended[slice.index(extended,n_dimensions) <= old_dim[n_dimensions]] <- obj
  return(extended)
}

sample_store_base <- function(data, par_names, iters = 1, stage = "init", ...) {
  subject_ids <- unique(data$subjects)
  n_pars <- length(par_names)
  n_subjects <- length(subject_ids)
  samples <- list(
    alpha = array(NA_real_,dim = c(n_pars, n_subjects, iters),dimnames = list(par_names, subject_ids, NULL)),
    stage = array(stage, iters),
    subj_ll = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL))
  )
}

set_p_accept <- function(stage, search_width){
  # Proposals distributions:
  # 1. Prior - unscaled: all stages
  # 2. Prev particle - scaled chain variance: all stages in preburn scaled by prior variance
  # 3. Chain mean - scaled chain variance: burn onwards
  # 4. Eff mean - scaled eff variance: sample onwards
  # adapt/sample: the local component targets ~3% of particles beating the
  # current value (where the ensemble move's ESS per iteration peaks);
  # independent components are not scaled.
  if(stage == "preburn") return(0.02 * (1/search_width))
  if(stage == "burn") return(c(0.02, 0.25)* (1/search_width))
  if(stage == "adapt") return(c(0.03, 0.25)* (1/search_width))
  if(stage == "sample") return(c(0.03, 0.3, 0.3)* (1/search_width))
}

get_default_mix <- function(stage){
  if (stage == "burn") {
    default_mix <- c(0.15, 0.35, 0.5)
  } else if(stage == "adapt"){
    default_mix <- c(0.1, 0.4, 0.4)
  } else if(stage == "sample"){
    default_mix <- c(0.05, 0.3, 0.3, 0.35)
  }  else{
    default_mix <- c(0.5, 0.5)
  }
  return(default_mix)
}

check_mix <- function(mix = NULL, stage) {
  default_mix <- get_default_mix(stage)
  if(is.null(mix)) mix <- default_mix
  if(length(mix) < length(default_mix)){
    if(length(default_mix) - length(mix) == 1){
      mix <- mix*(1-default_mix[length(default_mix)])
      mix <- c(mix, default_mix[length(default_mix)])
    } else{
      stop("mix settings are not compatible with stage")
    }
  }
  return(mix)
}

check_epsilon <- function(epsilon, n_pars, mix) {
  if (is.null(epsilon)) { # In preburn case, there's only one epsilon here
    if (n_pars > 15) {
      epsilon <- .5
    } else if (n_pars > 10) {
      epsilon <- .6
    } else {
      epsilon <- .7
    }
    # The first proposal is always unscaled
    # Every subsequent phase at most adds one epsilon
  } else if(length(epsilon) < (length(mix) -1)){
    epsilon <- c(epsilon, epsilon[length(epsilon)])
  }
  return(epsilon)
}

check_prop_performance <- function(prop_performance, stage){
  default_mix <- get_default_mix(stage)
  if(is.null(prop_performance) || stage == "adapt") prop_performance <- rep(0, length(default_mix))
  if(length(prop_performance) < length(default_mix)){
    if(length(default_mix) - length(prop_performance) == 1){
      n_total <- sum(prop_performance)
      prop_performance <- prop_performance*(1-default_mix[length(default_mix)])
      prop_performance <- c(prop_performance, default_mix[length(default_mix)]*n_total)
    } else{
      stop("prop_performance settings are not compatible with stage")
    }
  }
  return(round(prop_performance))
}

calc_ll_manager <- function(proposals, dadm, model, component = NULL, r_cores = 1, return_trialwise=FALSE){
  if(!is.data.frame(dadm)){
    lls <- log_likelihood_joint(proposals, dadm, model, component)
  } else{
    model <- model()
    if(is.null(model$c_name)){ # use the R implementation
      lls <- apply(proposals,1, calc_ll_R, model, dadm = dadm)
    } else {
      p_types <- names(model$p_types)
      designs <- get_designs_expanded(dadm, model)
      constants <- attr(dadm, "constants")
      if(is.null(constants)) constants <- NA

      backend   <- getOption("emc.ll_backend", default = "multiprocess")
      n_threads <- getOption("emc.n_threads", default = 1)

      if (backend == "multithreaded") {
        lls <- calc_ll_multithreaded(proposals, dadm, constants = constants, designs = designs,
                                     type = model$c_name, model$bound, model$transform,
                                     model$pre_transform, p_types = p_types,
                                     min_ll = log(1e-10), model$trend, n_threads = n_threads, return_trialwise=return_trialwise)
      } else {
        lls <- calc_ll(proposals, dadm, constants = constants, designs = designs,
                       type = model$c_name, model$bound, model$transform,
                       model$pre_transform, p_types = p_types, min_ll = log(1e-10),
                       model$trend, return_trialwise=return_trialwise)
      }
    }
  }
  return(lls)
}

#' Run a Group-level Model.
#'
#' Separate function for running only the group-level model. This can be useful in a
#' two-step analysis. Works similar in functionality to make_emc,
#' except also does the fitting and returns an emc object that works with
#' most posterior checking tests (but not the data generation/posterior predictives).
#'
#' @param prior an emc.prior object.
#' @param iter Number of MCMC samples to collect.
#' @inheritParams make_emc
#'
#' @returns an emc object with only group-level samples
#' @export
run_hyper <- function(type = "standard", data, prior = NULL, iter = 1000, n_chains =3, ...){
  args <- list(...)
  if(length(dim(data)) == 3){
    data_input <- data
    data <- as.data.frame(t(data_input[,,1]))
    data$subjects <- 1:nrow(data)
    iter <- dim(data_input)[3]
    is_mcmc <- T
    pars <- rownames(data_input)
  } else{
    data_input <- data[,colnames(data)!= "subjects"]
    is_mcmc <- F
    pars <- colnames(data_input)
  }
  emc <- list()
  for(j in 1:n_chains){
    samples <- sample_store(data = data ,par_names = pars, integrate = F, type = type, ...)
    subjects <- unique(data$subjects)
    sampler <- list(
      data = split(data, data$subjects),
      par_names = pars,
      subjects = subjects,
      n_pars = length(pars),
      n_subjects = length(subjects),
      samples = samples,
      init = TRUE
    )
    class(sampler) <- "pmwgs"
    sampler <- add_info(sampler, prior, type = type, ...)
    sampler$type <- type
    sampler <- extend_sampler(sampler, iter - 1, "sample")
    draws <- if (is_mcmc) data_input else array(t(as.matrix(data_input)), c(length(pars), length(subjects), 1L))
    emc[[j]] <- native_run_hyper(sampler, draws, iter)
  }
  emc[[1]]$type <- type
  class(emc) <- "emc"
  emc <- subset(emc, filter = 1)
  return(emc)
}

check_CR <- function(emc, p_vector, range = .2, N = 500){
  covs <- diag(length(p_vector)) * range
  props <- mvtnorm::rmvnorm(N, mean = p_vector, sigma = covs)
  model <- emc[[1]]$model
  if(is.null(model()$c_name)) stop("C not implemented yet for this model")
  dat <- emc[[1]]$data[[1]]
  modelRlist <- model()
  modelRlist$c_name <- NULL
  modelR <- function()return(modelRlist)
  t1 <- system.time(
    R <- calc_ll_manager(props, dat, modelR)
  )
  t2 <- system.time(
    C <- calc_ll_manager(props, dat, model)
  )
  print(paste0("C ", t1$elapsed/t2$elapsed, " times faster"))
  if(!identical(C, R)){
    warning("C and R results differ")
  }
  return(list(C = C, R = R))
}
