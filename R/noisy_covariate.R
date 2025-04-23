# Major issue:
# calc_ll_manager only returns lls, but I need to update running mean, variance
# and scaling. Current idea: Attach it as an attribute to lls, and pass this
# attribute to  the dadm in calc_ll_manger, which is then passed back to the next ll calculation.


# Sampling ----------------------------------------------------------------

# One potential minor issue, some data will have lapse responses that will never get
# an acceptance for the latent states.
cov_noise_ll <- function(pars, dadm, model, min_ll = log(1e-10)){
  # This step adds the likelihood of the noisy covariate model.
  ll <- 0
  for(cov in model$noisy_cov){
    # Add the latent states to the dadm
    dadm <- fill_dadm(dadm, attr(dadm, "noisy_cov")[[cov]], cov)
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

fill_dadm <- function(dadm, noisy_cov, cov_name){
  # Here we could perform some indexing needed for race models!
  ### SOME INDEXING ###
  # Here we're replacing the observed variable with the latent variable!
  dadm[,cov_name] <- noisy_cov$latent
  des <- attr(dadm, "designs")
  des$v[,2] <- noisy_cov$latent
  attr(dadm, "designs") <- des
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
  noisy_covs <- list()
  for(cov in model$noisy_cov){
    latent_pars <- get_noisy_cov_pnames(cov)

    noisy_cov <- attr(dadm, "noisy_cov")[[cov]]
    prev_latent <- noisy_cov$latent
    # Now for k proposals we evaluate the likelihood given the new latent states
    # For now we'll sample from 3 densities.
    # 1. Density centered on the observed value with some adaptive covariance
    # 2. Density centered on previous z_i, with same adaptive covariance
    # 3. Density centered on running z_i mean, with some adaptive covariance
    cur_sd <- noisy_cov$epsilon * noisy_cov$running_sd
    # For now we'll just assume that each distributions carries 1/3 of the particles
    for(i in 1:(n_particles*n_dists + 1)){
      if(i != 1){ # The old state is particle 1
        if(i < (1 + n_particles)){
          noisy_cov$latent <- rnorm(nrow(noisy_cov), mean = noisy_cov$obs, sd = exp(pars[latent_pars[3]]))
        } else if(i < (1 + 2*n_particles)){
          noisy_cov$latent <- rnorm(nrow(noisy_cov), mean = prev_latent, sd = exp(pars[latent_pars[3]]))
        } else{
          noisy_cov$latent <- rnorm(nrow(noisy_cov), mean = noisy_cov$running_mu, sd = cur_sd)
        }
      }
      # Calculate the likelihood
      ll <- calc_ll_R(pars, model, fill_dadm(dadm, noisy_cov, cov), return_sum = FALSE)
      # For now assign each distribution a 1/3 weight
      lm <- 1/3*dnorm(noisy_cov$latent, mean = noisy_cov$obs, sd = exp(pars[latent_pars[3]])) +
        1/3*dnorm(noisy_cov$latent, mean = prev_latent, sd = exp(pars[latent_pars[3]])) +
        1/3*dnorm(noisy_cov$latent, mean = noisy_cov$running_mu, sd = cur_sd)
      # Importance correction of all proposal densities
      l <- ll - log(lm)
      latents[,i] <- noisy_cov$latent
      # Unnormalized weights
      weights[,i] <- l
    }

    # Now we need to sample from the weights
    # using n metropolis steps, with n being the number of observations
    mh_idx <- apply(weights, 1, function(x){
      weights <- exp(x - max(x))
      sample(x = n_particles*n_dists + 1, size = 1, prob = weights)
    })
    # Sample each row of latents based on weight and replace
    noisy_cov$latent <- latents[cbind(seq_len(nrow(noisy_cov)), mh_idx)]
    # Finally let's update some sampling tuning parameters
    # First update acceptance
    noisy_cov$acceptance <- calc_acceptance(weights, noisy_cov$acceptance, noisy_cov$idx)
    # Scaling parameter
    noisy_cov$epsilon <- update_epsilon_continuous(noisy_cov$epsilon, noisy_cov$acceptance, target = .2,
                                              iter =noisy_cov$idx, d = 1, alphaStar = 3,
                                              clamp = c(.5, 2))
    # Running mean:
    new_moments <- running_moments(noisy_cov$latent, noisy_cov$running_mu, noisy_cov$running_sd, noisy_cov$idx)
    noisy_cov$running_mu <- new_moments$mu
    noisy_cov$running_sd <- new_moments$sd
    noisy_cov$idx <- noisy_cov$idx + 1
    noisy_covs[[cov]] <- noisy_cov
    dadm <- fill_dadm(dadm, noisy_cov, cov)
  }
  attr(dadm, "noisy_cov") <- noisy_covs
  return(dadm)
}

calc_acceptance <- function(weights, prev_acc, idx){
  # Subtract weight from previous particle
  n <- ncol(weights)
  new <- weights[,2:n] - rep(weights[,1], times = n-1)
  curr_acc <- rowSums(new > 0)
  # Calculate total accept counts, add the current one
  updated_acc <- ((prev_acc * idx * (n-1)) + curr_acc) / ((idx + 1) * (n-1))
  return(updated_acc)
}

running_moments <- function(new_latent, old_mu, old_sd, idx){
  n_new   <- idx + 1
  delta   <- new_latent - old_mu
  meanNew <- old_mu + delta / n_new
  new_var <- ((idx - 1)*(old_sd^2) + (delta^2 * idx / n_new)) / (n_new - 1)
  list(
    mu = meanNew,
    sd   = sqrt(new_var)
  )
}


add_noisy_cov_dadm <- function(dadm, model, noisy_cov){
  if(is.data.frame(dadm)){
    if(!is.null(model()$noisy_cov)){
      noisy_covs <- model()$noisy_cov
      attr(dadm, "noisy_cov") <- noisy_cov[noisy_covs]
    }
  } else{
    for(j in 1:length(dadm)){
      if(!is.null(model[[j]]()$noisy_cov)){
        noisy_covs <- model[[j]]()$noisy_cov
        attr(dadm[[j]], "noisy_cov") <- noisy_cov[noisy_covs]
      }
    }
  }
  return(dadm)
}

find_noisy_cov <- function(pm_settings){
  return(!is.null(pm_settings$noisy_cov))
}



# Design ------------------------------------------------------------------


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



