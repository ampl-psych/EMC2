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
  n_particles <- 10
  lls <- matrix(NA, nrow = nrow(dadm), ncol = n_particles + 1)
  latents <- matrix(NA, nrow = nrow(dadm), ncol = n_particles +1)
  lls[,1] <- calc_ll_R(pars, model, dadm)
  latents[,1] <- dadm$latent
  # Now for k proposals we evaluate the likelihood given the new latent states
  for(i in 1:n_particles){
    dadm$latent <- rnorm(nrow(dadm), mean = dadm$obs, sd = dadm$running_sd)
    lls[,i+1] <- calc_ll_R(pars, model, dadm)
    # Some importance correction of all proposal densities
    #<here>
    latents[,i+1] <- dadm$latent
  }

  # Then do n metropolis steps, with n being the number of observations
  # Sample each row based on weight
  dadm$latent <- sampled_latent
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


