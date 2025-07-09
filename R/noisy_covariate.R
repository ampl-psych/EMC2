# Steps:
# At each iteration, for each subject:
# 1. Update latent states (update_latent_cov)
#   - For i in 1:n_particles
#     a. Draw a particle for each row (data point)
#     b. Update design with that particle
#     c. Calculate ll for each row
#   - Importance correction then to weight per row for each column
#   - MH step and take 1 particle (columns) per row
#   - Update running moments, epsilon and acceptance
#   - Update dadm design with winning particle per row
# 2. Update other parameters
#   - Keep winning state of particles as stable design
#   - Do the standard steps
# 3. Return noisy covariate states as attribute to ll and update in new_particle

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
    # This will return a vector of names: mean, total_sd, and error_proportion
    latent_pars <- get_noisy_cov_pnames(cov)

    # Transform reparameterized parameters back to original scale
    sigma_total <- pars[,latent_pars[2]]  # sigma_total
    rho <- pars[,latent_pars[3]]          # rho (measurement error proportion)

    # Recover original parameters
    sigma_e <- sigma_total * sqrt(rho)           # measurement error SD
    tau <- sigma_total * sqrt(1 - rho)           # trial-to-trial variability SD

    # True latent mean with latent variability
    latent_ll <- dnorm(dadm[,cov], mean = pars[,latent_pars[1]], sd = tau, log = TRUE)
    # Measurement error going from the latent mean to the measurement noise
    noise_ll <- dnorm(dadm$obs, mean = dadm[,cov], sd = sigma_e, log = TRUE)
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
  n_particles <- 25 # Number of particles per proposal distribution
  n_dists <- 3 # Number of proposal distributions

  weights <- matrix(NA, nrow = nrow(dadm), ncol = n_particles*n_dists + 1)
  latents <- matrix(NA, nrow = nrow(dadm), ncol = n_particles*n_dists +1)
  noisy_covs <- list()
  all_pars <- c(pars, attr(dadm, "constants"))
  for(cov in model$noisy_cov){
    latent_pars <- get_noisy_cov_pnames(cov)

    noisy_cov <- attr(dadm, "noisy_cov")[[cov]]
    prev_latent <- noisy_cov$latent

    # Transform reparameterized parameters back to original scale
    # mu_z       <- pars[latent_pars[1]]                    # μ_z
    # sigma_total <- pars[latent_pars[2]]                   # total SD
    # rho        <- pars[latent_pars[3]]                    # error proportion
    # sigma_e    <- sigma_total * sqrt(rho)                 # σ_error  > 0
    # tau        <- sigma_total * sqrt(1 - rho)             # σ_true   > 0
    #
    # prec_true <- 1 / tau^2
    # prec_err  <- 1 / sigma_e^2
    # post_prec <- prec_true + prec_err
    # post_sd   <- 1 / sqrt(post_prec)
    # post_mu   <- (prec_true * mu_z + prec_err * noisy_cov$obs) / post_prec
    #

    # Now for k proposals we evaluate the likelihood given the new latent states
    # For now we'll sample from 3 densities.
    # 1. Density centered on the observed value with some adaptive covariance
    # 2. Density centered on previous z_i, with same adaptive covariance
    # 3. Density centered on running z_i mean, with some adaptive covariance
    # 4. Density centered on approximate posterior mean giving measurement and error normal-modeled
    cur_sd <- noisy_cov$epsilon * noisy_cov$running_sd
    # For now we'll just assume that each distributions carries 1/3 of the particles
    for(i in 1:(n_particles*n_dists + 1)){
      if(i != 1){ # The old state is particle 1
        if(i < (1 + n_particles)){
          noisy_cov$latent <- rnorm(nrow(noisy_cov), mean = noisy_cov$obs, sd = cur_sd)
        } else if(i < (1 + 2*n_particles)){
          noisy_cov$latent <- rnorm(nrow(noisy_cov), mean = prev_latent, sd = cur_sd)
        } else {#if(i < (1 + 3*n_particles)){
          noisy_cov$latent <- rnorm(nrow(noisy_cov), mean = noisy_cov$running_mu, cur_sd)
        }
        # } else{
        #   noisy_cov$latent <- rnorm(nrow(noisy_cov), mean = post_mu, post_sd)
        # }
      }
      # Calculate the likelihood
      ll <- calc_ll_R(pars, model, fill_dadm(dadm, noisy_cov, cov), return_sum = FALSE)
      # For now assign each distribution a 1/3 weight
      lm <- 1/n_dists*dnorm(noisy_cov$latent, mean = noisy_cov$obs, sd = cur_sd) +
        1/n_dists*dnorm(noisy_cov$latent, mean = prev_latent, sd = cur_sd) +
        1/n_dists*dnorm(noisy_cov$latent, mean = noisy_cov$running_mu, sd = cur_sd) #+
        # 1/n_dists*dnorm(noisy_cov$latent, mean = post_mu, sd = post_sd)
        # 1/n_dists*dnorm(noisy_cov$latent, mean = noisy_cov$obs, sd = exp(pars[latent_pars[3]]))
      # Importance correction of all proposal densities
      l <- ll - log(lm)
      latents[,i] <- noisy_cov$latent
      # Unnormalized weights
      weights[,i] <- l
    }

    # Now we need to sample from the weights
    # using n metropolis steps, with n being the number of observations
    mh_idx <- apply(weights, 1, function(logw) {        # loop over rows
      p <- exp(logw - max(logw))                     # stabilise & exponentiate
      p <- p / sum(p)                                # normalise to probabilities
      sample.int(length(p), 1, prob = p)             # draw a single column index
    })
    # Sample each row of latents based on weight and replace
    noisy_cov$latent <- latents[cbind(seq_len(nrow(noisy_cov)), mh_idx)]
    # Finally let's update some sampling tuning parameters
    # First update acceptance
    noisy_cov$acceptance <- calc_acceptance(weights, noisy_cov$acceptance, noisy_cov$idx)
    # Scaling parameter
    noisy_cov$epsilon <- update_epsilon_continuous(noisy_cov$epsilon, noisy_cov$acceptance, target = .1,
                                              iter =noisy_cov$idx, d = 1, alphaStar = 3,
                                              clamp = c(.5, 15))

    noisy_cov$is1[mh_idx ==1] <- noisy_cov$is1[mh_idx ==1] + 1
    noisy_cov$is2[mh_idx %in% (2:(1+n_particles))] <- noisy_cov$is2[mh_idx %in% (2:(1+n_particles))] + 1
    noisy_cov$is3[mh_idx %in% (n_particles+2):(2*n_particles+1)] <- noisy_cov$is3[mh_idx %in% (n_particles+2):(2*n_particles+1)] + 1
    noisy_cov$is4[mh_idx %in% (2*n_particles+2):(3*n_particles+1)] <- noisy_cov$is4[mh_idx %in% (2*n_particles+2):(3*n_particles+1)] + 1
    # noisy_cov$is5[mh_idx %in% (3*n_particles+2):(4*n_particles+1)] <- noisy_cov$is5[mh_idx %in% (3*n_particles+2):(4*n_particles+1)] + 1
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

# update_latent_cov <- function(pars, dadm, model){
#   # This function updates the latent states using a particle Metropolis-Hastings
#   # approach with 3 proposal distributions.
#   n_particles <- 10 # Number of particles per proposal distribution
#   n_dists <- 3
#
#   noisy_covs <- attr(dadm, "noisy_cov")
#
#   for(cov in model$noisy_cov){
#     ## ---------- global (per‑covariate) parameters ----------
#     latent_pars <- get_noisy_cov_pnames(cov)
#     mu_z    <-           pars[ latent_pars[1] ]       # μ_z
#     tau     <-  exp(     pars[ latent_pars[2] ] )     # σ_true   > 0
#     sigma_e <-  exp(     pars[ latent_pars[3] ] )     # σ_error  > 0
#
#     ## ---------- current latent states & data ----------
#     noisy_cov   <- noisy_covs[[cov]]
#     prev_latent <- noisy_cov$latent
#     c_i         <- noisy_cov$obs
#     n_obs       <- length(prev_latent)
#
#     ## ---------- Set up storage for particles and weights ----------
#     n_total_particles <- n_particles * n_dists + 1
#     weights <- matrix(NA, nrow = n_obs, ncol = n_total_particles)
#     latents <- matrix(NA, nrow = n_obs, ncol = n_total_particles)
#     latents[, 1] <- prev_latent
#
#     ## ---------- p(z | c, θ) for independence proposal ----------
#     prec_true <- 1 / tau^2
#     prec_err  <- 1 / sigma_e^2
#     post_prec <- prec_true + prec_err
#     post_sd   <- noisy_cov$epsilon * noisy_cov$running_sd #1 / sqrt(post_prec)
#     post_mu   <- (prec_true * mu_z + prec_err * c_i) / post_prec
#
#     ## ---------- Generate particles from 3 proposal distributions ----------
#     scaled_running_sd <- noisy_cov$epsilon * noisy_cov$running_sd
#
#     for(p in 1:n_particles) {
#       # Proposal 1: Adaptive Mean
#       latents[, 1 + p] <- rnorm(n_obs, mean = noisy_cov$running_mu, sd = scaled_running_sd)
#       # Proposal 2: Adaptive Random Walk
#       latents[, 1 + n_particles + p] <- rnorm(n_obs, mean = prev_latent, sd = scaled_running_sd)
#       # Proposal 3: Conjugate Posterior
#       latents[, 1 + 2*n_particles + p] <- rnorm(n_obs, mean = noisy_cov$obs, sd = post_sd)
#     }
#
#     ## ---------- Calculate importance weights for all particles ----------
#     for(i in 1:n_total_particles) {
#         current_latents <- latents[,i]
#
#         # 1. Calculate Likelihood (EAM + Measurement)
#         ncov_tmp <- noisy_cov
#         ncov_tmp$latent <- current_latents
#         ll <- calc_ll_R(pars, model, fill_dadm(dadm, ncov_tmp, cov), return_sum = FALSE)
#
#         # 2. Calculate proposal mixture density
#         d1 <- dnorm(current_latents, mean = noisy_cov$running_mu, sd = scaled_running_sd)
#         d2 <- dnorm(current_latents, mean = prev_latent, sd = scaled_running_sd)
#         d3 <- dnorm(current_latents, mean = post_mu, sd = post_sd)
#         log_q <- log( (d1 + d2 + d3) / n_dists )
#
#         # 3. Calculate log prior of z
#         log_prior <- dnorm(current_latents, mu_z, tau, log = TRUE)
#
#         # 4. Full unnormalized log weight
#         weights[,i] <- ll + log_prior - log_q
#     }
#
#     ## ---------- Sample one particle per trial based on weights ----------
#     mh_idx <- apply(weights, 1, function(logw) {
#       logw[is.na(logw)] <- -Inf
#       p <- exp(logw - max(logw))
#       p <- p / sum(p)
#       if(any(is.na(p)) || all(p == 0)) return(1)
#       sample.int(length(p), 1, prob = p)
#     })
#     z_new <- latents[cbind(seq_len(n_obs), mh_idx)]
#
#     ## ---------- Bookkeeping and tuning ----------
#     noisy_cov$acceptance <- calc_acceptance(weights, noisy_cov$acceptance, noisy_cov$idx)
#
#     # Update scaling parameter epsilon
#     noisy_cov$epsilon <- update_epsilon_continuous(noisy_cov$epsilon, noisy_cov$acceptance, target = .25,
#                                                 iter = noisy_cov$idx, d = 1, alphaStar = 0.25,
#                                                 clamp = c(.1, 10))
#
#     # Update running moments
#     rmom <- running_moments(z_new, noisy_cov$running_mu,
#                             noisy_cov$running_sd, noisy_cov$idx)
#     noisy_cov$running_mu <- rmom$mu
#     noisy_cov$running_sd <- rmom$sd
#
#     noisy_cov$latent <- z_new
#     noisy_cov$idx <- noisy_cov$idx + 1
#     noisy_covs[[cov]] <- noisy_cov
#     dadm <- fill_dadm(dadm, noisy_cov, cov)
#   }
#   attr(dadm, "noisy_cov") <- noisy_covs
#   return(dadm)
# }

# prev_acc : running average acceptance rate (vector length = nrow(weights))
# iter     : how many *previous* batches of proposals have already been folded in
calc_acceptance <- function(weights, prev_acc, iter) {
  k <- ncol(weights) - 1                        # proposals per row this batch

  diff      <- weights[ , -1, drop = FALSE] - weights[ , 1]
  curr_acc  <- rowSums(diff > 0)               # accepted in this batch (0 … k)

  # running average over (iter + 1) batches
  updated_acc <- (prev_acc * iter * k + curr_acc) / ((iter + 1) * k)
  updated_acc
}



running_moments <- function(new_latent, old_mu, old_sd, idx){
  n_new   <- idx + 1
  delta   <- new_latent - old_mu
  meanNew <- old_mu + delta / n_new
  new_var <- ((idx - 1)*(old_sd^2) + (delta^2 * idx / n_new)) / (n_new - 1)
  list(
    mu = meanNew,
    sd   = pmax(sqrt(new_var), 0.01)
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

# Design ------------------------------------------------------------------


get_noisy_cov_pnames <- function(noisy_cov){
  return(paste0(c("mu.t", "sigma.total", "rho"), "_", noisy_cov))
}

update_model_noisy_cov <- function(noisy_cov, model) {
  # Get model list to modify
  model_list <- model()

  # For each noisy_covariate
  for (cov in noisy_cov) {
    p_names <- get_noisy_cov_pnames(cov)
    # Transforms: identity for mean, exp for total_sd, pnorm for proportion
    transforms <- setNames(c("identity", "exp", "pnorm"), p_names)
    # Bounds: mean unrestricted, total_sd > 0, proportion in (0,1)
    minmax <- setNames(list(mu = c(-Inf,Inf), st = c(0.01,Inf), rho = c(0.01,0.99)), p_names)
    model_list$transform$func <- c(model_list$transform$func, transforms)
    model_list$bound$minmax <- cbind(model_list$bound$minmax, do.call(cbind, minmax))

    model_list$p_types <- c(model_list$p_types, setNames(numeric(length(p_names)), p_names))
  }
  model_list$noisy_cov <- noisy_cov
  # Return updated model function
  model <- function() { return(model_list) }
  return(model)
}


check_noisy_cov <- function(noisy_cov, covariates, formula = NULL) {
  if (!all(noisy_cov %in% covariates)){
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



