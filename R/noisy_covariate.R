add_cov_noise_ll <- function(pars, dadm, min_ll = log(1e-10)){
  # This step adds the likelihood of the noisy covariate model.
  mu_latent <- pars[,'mu_latent']
  sigma_latent <- pars[,'sigma_latent']
  sigma_noise <- pars[,'sigma_noise']
  # True latent mean with latent variability
  latent_ll <- dnorm(dadm$latent, mean = mu_latent, sd = sigma_latent, log = TRUE)
  # Measurement error going from the latent mean to the measurement noise
  noise_ll <- dnorm(dadm$obs, mean = dadm$latent, sd = sigma_noise, log = TRUE)
  return(pmax(min_ll, latent_ll + noise_ll))
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
    ll <- calc_ll_R(pars, model, dadm)
    lm <- 1/3*dnorm(dadm$latent, mean = dadm$obs, sd = cur_sd) +
      1/3*dnorm(dadm$latent, mean = prev_latent, sd = cur_sd) +
      1/3*dnorm(dadm$latent, mean = dadm$running_mu, sd = cur_sd)
    # Some importance correction of all proposal densities
    l <- ll - log(lm)
    latents[,i] <- dadm$latent
    # For now these weights are uncorrected for
    weights[,i] <- l
  }
  # Now we need to sample from the weighted likelihoods
  # using n metropolis steps, with n being the number of observations
  mh_idx <- apply(weights, 1, function(x){
    weights <- exp(x - max(x))
    sample(x = n_particles*n_dists + 1, size = 1, prob = weights)
  })
  # Sample each row based on weight
  dadm$latent <- latents[cbind(seq_len(nrow(dadm)), mh_idx)]
  # For race models we also need to keep a tracker that only the trials (not the lR)
  # are updated (so 1 per x accumulators, not per row)
  return(dadm)
}

sample_latent_cov <- function(dadm){
  # This function samples the latent states given the current parameters
  # and the data. It returns a vector of length nrow(dadm) with the
  # sampled latent states.

  latent <- cbind(dadm$latent, latent_new)
  return(latent)
}


