# One potential issue, some data will have lapse responses that will never get
# an acceptance for the latent states.

add_cov_noise_ll <- function(pars, dadm, model, min_ll = log(1e-10)){
  # This step adds the likelihood of the noisy covariate model.
  ll <- 0
  for(cov in model$noisy_cov){
    # Add the latent states to the dadm
    dadm <- fill_dadm(dadm, cov)
    # dadm[,cov] is actually the latent now!
    # This will return a vector of names, the first the mean, 2nd the latent sd, 3rd the error sd
    latent_pars <- get_noisy_cov_pnames(cov)
    # True latent mean with latent variability
    latent_ll <- dnorm(dadm[,cov], mean = pars[,latent_pars[1]], sd = pars[,latent_pars[2]], log = TRUE)
    # Measurement error going from the latent mean to the measurement noise
    noise_ll <- dnorm(dadm$obs, mean = dadm[,cov], sd = pars[,latent_pars[3]], log = TRUE)
    ll <- ll + pmax(min_ll, latent_ll + noise_ll)
  }
  return(ll)
}

fill_dadm <- function(dadm, cov){
  noisy_cov <- attr(dadm, cov)
  # Here we're replacing the observed variable with the latent variable!
  dadm[,cov] <- noisy_cov$latent
  # Add the actually observed noisy covariate
  dadm$obs <- noisy_cov$obs
  return(dadm)
}

update_latent_cov <- function(pars, dadm, model){
  # Overall strategy: update the true ''parameters'', not data-augmented
  # latent states, and then update the latent states.

  # Here we have previously updated our true parameters,
  # now we need to update the latent states.
  # We'll update the latent states using multiple draws per trial that are all evaluated
  # Ideally we'd do the generation of the states in C++ so that we can loop over
  # all proposals in C++.

  # This function should update the latent parameter attribute
  # of dadm, as well as the mean and variance of the latent parameter.

  # To construct the weights of each proposal we need to first evalue the likelihood
  # for the different proposals, holding the true parameters constant, and only
  # updating each latent state. We'll make a proposal for each latent state
  # and evaluate the likelihood for all these simultaneously and add density of proposals
  # to account for detailed balance.
  # We can use the proposal density to weight the likelihood of each proposal
  # and then sample from the weighted likelihoods to get the new latent states.
  n_particles <- 5 # Number of particles per proposal distribution
  n_dists <- 3 # Number of proposal distributions
  weights <- matrix(NA, nrow = nrow(dadm), ncol = n_particles*n_dists + 1)
  latents <- matrix(NA, nrow = nrow(dadm), ncol = n_particles*n_dists +1)
  prev_latent <- dadm$latent
  # Now for k proposals we evaluate the likelihood given the new latent states
  # For now we'll sample from 3 densities.
  # 1. Density centered on the observed value with some adaptive covariance
  # 2. Density centered on previous z_i, with same adaptive covariance
  # 3. Density centered on running z_i mean, with some adaptive covariance
  cur_sd <- dadm$eps *dadm$running_sd
  for(i in 1:(n_particles*n_dists + 1)){
    if(i != 1){ # The old state is particle 1
      if(i < (1 + n_particles)){
        dadm$latent <- rnorm(nrow(dadm), mean = dadm$obs, sd = cur_sd)
      } else if(i < (1 + 2*n_particles)){
        dadm$latent <- rnorm(nrow(dadm), mean = prev_latent, sd = cur_sd)
      } else{
        dadm$latent <- rnorm(nrow(dadm), mean = dadm$running_mu, sd = cur_sd)
      }
    }
    # Calculate the likelihood
    ll <- calc_ll_R(pars, model, dadm)
    # For now assign each distribution a 1/3 weight
    lm <- 1/3*dnorm(dadm$latent, mean = dadm$obs, sd = cur_sd) +
      1/3*dnorm(dadm$latent, mean = prev_latent, sd = cur_sd) +
      1/3*dnorm(dadm$latent, mean = dadm$running_mu, sd = cur_sd)
    # Importance correction of all proposal densities
    l <- ll - log(lm)
    latents[,i] <- dadm$latent
    # For now these weights are uncorrected for
    weights[,i] <- l
  }
  # First update acceptance
  dadm$acceptance <- calc_acc(weights, dadm$acceptance, idx)
  update_eps <- update_epsilon_continuous(dadm$epsilon, dadm$acceptance, target = .2,
                                          iter =iter, d = 1, alphaStar = 3,
                                          clamp = c(.5, 2))
  # Now we need to sample from the weighted likelihoods
  # using n metropolis steps, with n being the number of observations
  mh_idx <- apply(weights, 1, function(x){
    weights <- exp(x - max(x))
    sample(x = n_particles*n_dists + 1, size = 1, prob = weights)
  })
  # Sample each row of latents based on weight and replace
  dadm$latent <- latents[cbind(seq_len(nrow(dadm)), mh_idx)]
  # For race models we also need to keep a tracker that only the trials (not the lR)
  # are updated (so 1 per R accumulators, not per row)
  return(dadm)
}

calc_acceptance <- function(weights, prev_acc, idx){
  # Subtract weight from previous particle
  n <- ncol(weights)
  weights <- weights[,2:n] - weights[,1]
  curr_acc <- rowSums(weights > 0)
  # Calculate total accept counts, add the current one
  updated_acc <- ((prev_acc * idx * (n-1)) + curr_acc) / ((idx + 1) * (n-1))
  return(updated_acc)
}

get_noisy_cov_pnames <- function(noisy_cov){
  return(paste0(c("mu.t", "sigma.t", "sigma.e"), "_", noisy_cov))
}

update_model_noisy_cov <- function(noisy_cov, model) {
  # Get model list to modify
  model_list <- model()

  # For each noisy_covariate
  for (cov in noisy_cov) {
    p_names <- get_noisy_cov_pnames(cov)
    transforms <- setNames(c("identity", "exp", "exp"), p_names)
    model_list$transform$func <- c(model_list$transform$func, transforms)
    model_list$p_types <- c(model_list$p_types, setNames(numeric(length(p_names)), p_names))
  }
  model_list$noisy_cov <- noisy_cov
  # Return updated model function
  model <- function() { return(model_list) }
  return(model)
}

check_noisy_cov <- function(noisy_cov, covariates, formula = NULL) {
  if (!all(noisy_cov %in% names(covariates))){
    stop("noisy covariate not in defined covariates")
  }
  for(cov in noisy_cov){
    noisy_cov_pnames <- get_noisy_cov_pnames(cov)
    if (!is.null(formula)) {
      isin <-  noisy_cov_pnames %in% unlist(lapply(formula,function(x)all.vars(x)[1]))
      if(any(!isin)){
        # Add missing noise parameters to formula with intercept-only model
        formula <- c(formula, lapply(noisy_cov_pnames[!isin], function(x) as.formula(paste(x, "~ 1"))))
        message("Intercept formula added for noise_pars: ", paste(noisy_cov_pnames[!isin],collapse=", "))
      }
    }
  }
  return(formula)
}

latent_manager <- function(proposals, dadm, model, component = NULL){
  if(!is.data.frame(dadm)){
    lls <- log_likelihood_joint(proposals, dadm, model, component)
  } else{
    model <- model()
    if(is.null(model$c_name)){ # use the R implementation
      lls <- apply(proposals,1, calc_ll_R, model, dadm = dadm)
    } else{
      p_types <- names(model$p_types)
      designs <- list()
      for(p in p_types){
        designs[[p]] <- attr(dadm,"designs")[[p]][attr(attr(dadm,"designs")[[p]],"expand"),,drop=FALSE]
      }
      constants <- attr(dadm, "constants")
      if(is.null(constants)) constants <- NA
      lls <- calc_ll(proposals, dadm, constants = constants, designs = designs, type = model$c_name,
                     model$bound, model$transform, model$pre_transform, p_types = p_types, min_ll = log(1e-10),
                     model$trend)
    }
  }
  return(dadm)
}



