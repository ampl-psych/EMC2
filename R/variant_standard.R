sample_store_standard <- function(data, par_names, iters = 1, stage = "init", integrate = T, is_nuisance, ...) {
  subject_ids <- unique(data$subjects)
  n_subjects <- length(subject_ids)
  base_samples <- sample_store_base(data, par_names, iters, stage)
  par_names <- par_names[!is_nuisance]
  n_pars <- length(par_names)

  # Get total parameters including regressors
  group_designs <- list(...)$group_design
  total_par_names <- add_group_par_names(par_names, group_designs)

  # Check for random effects
  n_rand <- get_n_random(par_names, group_designs)
  n_s <- get_n_random_variance(par_names, group_designs) # Number of variance components

  # Create samples structure
  samples <- list(
    theta_mu = array(NA_real_,
      dim = c(length(total_par_names), iters),
      dimnames = list(total_par_names, NULL)
    ),
    theta_var = array(NA_real_,
      dim = c(n_pars, n_pars, iters),
      dimnames = list(par_names, par_names, NULL)
    ),
    a_half = array(NA_real_,
      dim = c(n_pars, iters),
      dimnames = list(par_names, NULL)
    )
  )
  if (n_rand > 0) {
    rand_names <- get_random_names(par_names, group_designs)
    samples$theta_u <- array(NA_real_, dim = c(n_rand, iters), dimnames = list(rand_names$u_names, NULL)) # Random effects
    samples$theta_s <- array(NA_real_, dim = c(n_s, iters), dimnames = list(rand_names$s_names, NULL)) # Variance components
    samples$a_half_s <- array(NA_real_, dim = c(n_s, iters), dimnames = list(rand_names$s_names, NULL)) # Auxiliary for s
  }

  if (integrate) samples <- c(samples, base_samples)
  return(samples)
}

add_info_standard <- function(sampler, prior = NULL, ...) {
  n_pars <- sum(!sampler$nuisance)
  group_design <- list(...)$group_design

  # Get all parameter names including regressors
  sampler$par_names_all <- add_group_par_names(sampler$par_names, group_design)

  sampler$par_group <- list(...)$par_groups
  if (is.null(sampler$par_group)) sampler$par_group <- rep(1, n_pars)
  sampler$is_blocked <- sampler$par_group %in% which(table(sampler$par_group) > 1)
  sampler$prior <- get_prior_standard(prior, n_pars,
    sample = F,
    group_design = group_design
  )
  sampler$group_designs <- group_design
  return(sampler)
}

calculate_mean_design <- function(group_designs, n_pars) {
  # Initialize the mean design matrix for each parameter
  mean_designs <- list()

  for (k in 1:n_pars) {
    if (is.null(group_designs) || is.null(group_designs[[k]])) {
      # If no design matrix, use a simple intercept
      mean_designs[[k]] <- matrix(1, nrow = 1, ncol = 1)
    } else {
      # Calculate the mean of the design matrix across subjects
      mean_designs[[k]] <- colMeans(group_designs[[k]], dims = 1)
    }
  }

  return(mean_designs)
}

calculate_implied_means <- function(mean_designs, beta_params, n_pars) {
  mu_implied <- numeric(n_pars)

  par_idx <- 0
  for (k in 1:n_pars) {
    x_k_mean <- mean_designs[[k]]
    mu_implied[k] <- x_k_mean %*% beta_params[par_idx + 1:length(x_k_mean)]
    par_idx <- par_idx + length(x_k_mean)
  }

  return(mu_implied)
}


get_prior_standard <- function(prior = NULL, n_pars = NULL, sample = TRUE, N = 1e5, selection = "mu", design = NULL, group_design = NULL,
                               par_groups = NULL) {
  # Checking and default priors
  if (is.null(prior)) {
    prior <- list()
  }
  if (!is.null(design)) {
    n_pars <- length(sampled_pars(design, doMap = F))
  }
  if (!is.null(prior$A) & is.null(n_pars)) {
    n_pars <- length(prior$A)
  }

  # Number of additional parameters from design matrices
  if (!is.null(group_design)) {
    n_additional <- n_additional_group_pars(group_design)
  } else {
    n_additional <- 0
  }

  # Set up combined theta_mu_mean to include both intercepts and slopes
  if (is.null(prior$theta_mu_mean)) {
    prior$theta_mu_mean <- rep(0, n_pars + n_additional)
  }
  if (is.null(prior$theta_mu_var)) {
    prior$theta_mu_var <- diag(rep(1, n_pars + n_additional))
  }
  if (is.null(prior$v)) {
    prior$v <- 2
  }
  if (is.null(prior$A)) {
    prior$A <- rep(.3, n_pars)
  }
  n_s <- 0
  if (!is.null(group_design)) {
    # Use parameter names from design or group_design to count random effects
    pnames <- if (!is.null(design)) names(sampled_pars(design, doMap = FALSE)) else names(group_design)
    n_s <- get_n_random_variance(pnames, group_design)
    n_rand <- get_n_random(pnames, group_design)

    if (n_s > 0) {
      if (is.null(prior$v_Z)) {
        prior$v_Z <- 2
      }
      if (is.null(prior$A_z)) {
        prior$A_z <- rep(0.15, n_s)
      }
    }
  }
  # Things I save rather than re-compute inside the loops.
  prior$theta_mu_invar <- ginv(prior$theta_mu_var) # Inverse of the matrix
  attr(prior, "type") <- "standard"
  out <- prior
  if (sample) {
    par_names <- names(sampled_pars(design, doMap = F))
    samples <- list()
    if (selection %in% c("mu", "beta", "alpha")) {
      # Sample beta (all parameters including regressors)
      beta <- t(mvtnorm::rmvnorm(N,
        mean = prior$theta_mu_mean,
        sigma = prior$theta_mu_var
      ))
      rownames(beta) <- add_group_par_names(par_names, group_design)
      if (selection %in% c("beta", "mu", "alpha")) {
        samples$theta_mu <- beta
        if (selection %in% c("mu", "alpha")) {
          samples$par_names <- par_names
          samples$group_designs <- group_design
        }
      }
    }
    if (selection %in% c("sigma2", "covariance", "correlation", "Sigma", "alpha")) {
      vars <- array(NA_real_, dim = c(n_pars, n_pars, N))
      colnames(vars) <- rownames(vars) <- par_names
      for (i in 1:N) {
        a_half <- 1 / rgamma(
          n = n_pars, shape = 1 / 2,
          rate = 1 / (prior$A^2)
        )
        attempt <- tryCatch(
          {
            vars[, , i] <- riwish(prior$v + n_pars - 1, 2 * prior$v * diag(1 / a_half))
          },
          error = function(e) e,
          warning = function(w) w
        )
        if (any(class(attempt) %in% c("warning", "error", "try-error"))) {
          sample_idx <- sample(1:(i - 1), 1)
          vars[, , i] <- vars[, , sample_idx]
        }
      }
      if (is.null(par_groups)) par_groups <- rep(1, n_pars)
      constraintMat <- matrix(0, n_pars, n_pars)
      for (i in 1:length(unique(par_groups))) {
        idx <- par_groups == i
        constraintMat[idx, idx] <- Inf
      }
      vars <- constrain_lambda(vars, constraintMat)
      if (selection != "alpha") samples$theta_var <- vars
    }

    # Sample random effects if present
    u <- NULL
    random_designs <- NULL
    random_maps <- NULL
    if (n_s > 0) {
      if (selection %in% c("s", "u", "alpha")) {
        # Sample s (variances)
        # s2 ~ IG(v_Z, A_z)
        # Actually s2_j ~ IG((v_Z+1)/2, A_z_j^2 * v_Z / 2)? No wait.
        # Gibbs step: s2 ~ IG((v+n)/2, (v*A^2 + sum(u^2))/2) ?
        # Prior for s2 implies:
        # We generally use similar to Sigma prior?
        # But here it's univariate variance per component.
        # Let's match the bridge sampler: prior_s <- dinvgamma(s2, shape = (v_Z + 1)/2, scale = (v_Z * A_z^2)/2)

        s2 <- matrix(NA, nrow = n_s, ncol = N)
        a_half_s <- matrix(NA, nrow = n_s, ncol = N)
        rand_names <- get_random_names(par_names, group_design)
        rownames(s2) <- rand_names$s_names
        rownames(a_half_s) <- rand_names$s_names

        for (j in 1:n_s) {
          # a_half_s ~ IG(1/2, 1/A_z^2)
          # shape = 1/2, rate = 1/(A_z^2)
          # R's rgamma uses rate, so IG sample is 1/rgamma(shape, rate)
          a_half_s[j, ] <- 1 / rgamma(N, shape = 1 / 2, rate = 1 / (prior$A_z[j]^2))

          # s2 | a_half_s ~ IG(v_Z/2, v_Z/a_half_s)
          # shape = v_Z/2, rate = v_Z/a_half_s
          # Note: rate for 1/rgamma is (v_Z / a_half_s)
          # s2 = 1 / rgamma(..., rate = v_Z / a_half_s)
          s2[j, ] <- 1 / rgamma(N, shape = prior$v_Z / 2, rate = prior$v_Z / a_half_s[j, ])
        }
        if (selection == "s") {
          samples$theta_s <- s2
          samples$a_half_s <- a_half_s
        }

        # Sample u
        # u ~ N(0, s2)
        if (selection %in% c("u", "alpha")) {
          # u is n_rand x N
          u <- matrix(NA, nrow = n_rand, ncol = N)
          rownames(u) <- rand_names$u_names

          # Map s2 to u components
          # Need to know which s corresponds to which u
          # get_random_names iterates same way
          idx_u <- 1
          idx_s <- 1

          # Re-iterate structure to check sizes
          for (p in names(group_design)) {
            if (is.list(group_design[[p]]) && !is.null(group_design[[p]]$random)) {
              for (gname in names(group_design[[p]]$random)) {
                z <- group_design[[p]]$random[[gname]]
                k_z <- ncol(z)
                # s2[idx_s, ] is the variance for this block
                sd_vec <- sqrt(s2[idx_s, ])

                # For each of k_z levels, sample N times
                # u_block: k_z x N
                for (k in 1:k_z) {
                  u[idx_u + k - 1, ] <- rnorm(N, mean = 0, sd = sd_vec)
                }

                idx_u <- idx_u + k_z
                idx_s <- idx_s + 1
              }
            }
          }
          if (selection == "u") samples$theta_u <- u

          # Prepare for alpha
          if (selection == "alpha") {
            random_designs <- add_group_design_random(names(group_design), group_design)
            random_maps <- add_group_design_map(names(group_design), group_design)
          }
        }
      }
    }

    if (selection %in% "alpha") {
      samples <- list()
      samples$alpha <- get_alphas(beta, vars, design[[1]]$Ffactors$subjects,
        group_design = group_design, u = u, random_designs = random_designs, random_maps = random_maps
      )
    }
    out <- samples
  }
  return(out)
}

get_startpoints_standard <- function(pmwgs, start_mu, start_var) {
  n_pars <- sum(!pmwgs$nuisance)
  n_total_pars <- length(pmwgs$prior$theta_mu_mean) # Includes regressor parameters
  if (is.null(start_mu)) start_mu <- rmvnorm(1, mean = pmwgs$prior$theta_mu_mean, sigma = pmwgs$prior$theta_mu_var)[1, ]
  # If no starting point for group var just sample some
  if (is.null(start_var)) start_var <- riwish(n_pars * 3, diag(n_pars))
  start_a_half <- 1 / rgamma(n = n_pars, shape = 2, rate = 1)

  # Calculate subject-specific means using design matrices
  group_designs <- add_group_design(pmwgs$par_names[!pmwgs$nuisance], pmwgs$group_designs, pmwgs$n_subjects)
  # Needs to include random effects if present
  # For startpoints, if u is not provided, we can either sample it or set to 0.
  # Code below needs subj_mu to be computed including u.

  n_rand <- get_n_random(pmwgs$par_names, pmwgs$group_designs)
  u <- NULL
  s <- NULL
  a_half_s <- NULL
  if (n_rand > 0) {
    n_s <- get_n_random_variance(pmwgs$par_names, pmwgs$group_designs)

    # Initialize a_half_s from prior: IG(1/2, 1/A_z^2)
    # A_z might be scalar or vector in prior.
    A_z <- if (is.null(pmwgs$prior$A_z)) rep(0.3, n_s) else rep(pmwgs$prior$A_z, length.out = n_s)
    a_half_s <- 1 / rgamma(n_s, shape = 1 / 2, rate = 1 / (A_z^2))

    # Initialize s (variance) from prior: IG(v_Z/2, v_Z/a_half_s)
    v_Z <- if (is.null(pmwgs$prior$v_Z)) 2 else pmwgs$prior$v_Z
    s <- 1 / rgamma(n_s, shape = v_Z / 2, rate = v_Z / a_half_s)

    # Re-using logic from get_prior_standard for u generation
    u <- numeric(n_rand)
    idx_u <- 1
    idx_s <- 1
    # Iterate group designs to match s to u
    for (p in names(pmwgs$group_designs)) {
      if (is.list(pmwgs$group_designs[[p]]) && !is.null(pmwgs$group_designs[[p]]$random)) {
        for (gname in names(pmwgs$group_designs[[p]]$random)) {
          z <- pmwgs$group_designs[[p]]$random[[gname]]
          k_z <- ncol(z)
          sd_val <- sqrt(s[idx_s])
          u[idx_u:(idx_u + k_z - 1)] <- rnorm(k_z, 0, sd_val)
          idx_u <- idx_u + k_z
          idx_s <- idx_s + 1
        }
      }
    }
  }

  random_designs <- add_group_design_random(pmwgs$par_names[!pmwgs$nuisance], pmwgs$group_designs)
  random_maps <- add_group_design_map(pmwgs$par_names[!pmwgs$nuisance], pmwgs$group_designs)
  subj_mu <- calculate_subject_means_RE(group_designs, start_mu, random_designs, u, random_maps)

  out <- list(
    tmu = start_mu, tvar = start_var, tvinv = ginv(start_var),
    a_half = start_a_half, subj_mu = subj_mu
  )

  if (!is.null(u)) {
    out$u <- u
    out$s <- s
    out$a_half_s <- a_half_s
  }
  return(out)
}

get_group_level_standard <- function(parameters, s) {
  # This function is modified for other versions
  mu <- parameters$subj_mu[, s]
  var <- parameters$tvar
  out <- list(mu = mu, var = var)
  if (!is.null(parameters$u)) out$u <- parameters$u
  if (!is.null(parameters$s)) {
    out$s <- parameters$s
    out$a_half_s <- parameters$a_half_s
  }
  return(out)
}

calculate_subject_means_RE <- function(group_designs, tmu, random_designs = NULL, u = NULL, random_maps = NULL) {
  # Calculate fixed effects
  mu <- calculate_subject_means(group_designs, tmu)

  # Add random effects if present
  if (!is.null(random_designs) && !is.null(u)) {
    # u is a concatenated vector. We need to split it based on random_maps.
    u_idx <- 1
    par_names <- names(group_designs)

    for (p in par_names) {
      if (!is.null(random_designs[[p]])) {
        # Loop over random components for this parameter
        z_list <- random_designs[[p]]
        maps <- random_maps[[p]]
        par_idx <- which(par_names == p) # Assuming row order matches list name order

        for (j in seq_along(z_list)) {
          z <- z_list[[j]]
          n_u <- ncol(z)
          u_curr <- u[u_idx:(u_idx + n_u - 1)]

          # Add to the appropriate row (parameter p)
          # mu is (n_pars x n_subjects)
          # z is (n_subjects x n_u)
          # u_curr is (n_u)
          # z %*% u_curr -> (n_subjects x 1)

          offset <- as.vector(z %*% u_curr)
          mu[par_idx, ] <- mu[par_idx, ] + offset

          u_idx <- u_idx + n_u
        }
      }
    }
  }
  return(mu)
}

fill_samples_standard <- function(samples, group_level, proposals, j = 1, n_pars) {
  samples$a_half[, j] <- group_level$a_half
  samples$last_theta_var_inv <- group_level$tvinv
  if (!is.null(group_level$u) && !is.null(samples$theta_u)) samples$theta_u[, j] <- group_level$u
  if (!is.null(group_level$s) && !is.null(samples$theta_s)) {
    samples$theta_s[, j] <- group_level$s
    samples$a_half_s[, j] <- group_level$a_half_s
  }
  samples <- fill_samples_base(samples, group_level, proposals, j = j, n_pars)
  return(samples)
}

gibbs_step_standard <- function(sampler, alpha) {
  #
  # alpha:  (p x n) subject-level parameter matrix
  #
  # sampler$group_designs:   list of length p, each (n x m_k)
  # sampler$is_blocked:    length p logical; is_blocked[k] = TRUE => param k is in some covariance block
  # sampler$par_group:     length p integer group labels (only used if is_blocked=TRUE)
  #
  # sampler$prior contains at least:
  #   - theta_mu_mean (length total_params) and theta_mu_invar (total_params x total_params)
  #   - v, A => hyperparams for the IW/Gamma updates
  #
  # last_sample_standard(...) must provide:
  #   - last$tmu (total_params-vector) of group means (includes both intercepts and slopes)
  #   - last$tvar, last$tvinv (p x p)
  #   - last$a_half (p)
  #
  # The function returns an updated (tmu, tvar, tvinv, a_half, alpha),
  # doing partial block updates for the covariance tvar, and a joint
  # update for all mu parameters (including former beta parameters)


  group_designs <- sampler$group_designs
  is_blocked <- sampler$is_blocked[!sampler$nuisance]
  par_group <- sampler$par_group[!sampler$nuisance]
  prior <- sampler$prior

  last <- last_sample_standard(sampler$samples)
  tmu <- last$tmu # group means, dimension total_params (includes former beta)
  tvar <- last$tvar # (p x p)
  tvinv <- last$tvinv
  a_half <- last$a_half # length p

  # Extract Random Effects state
  u <- last$u
  s <- last$s
  a_half_s <- last$a_half_s

  p <- nrow(alpha) # number of subject-level parameters
  n <- ncol(alpha) # number of subjects

  # Some backwards compatibility
  if (is.null(is_blocked) || length(is_blocked) == 0) {
    is_blocked <- rep(T, p)
  }
  if (is.null(par_group)) {
    par_group <- rep(1, p)
  }

  # Prepare designs
  fixed_designs <- add_group_design(sampler$par_names[!sampler$nuisance], group_designs, n)
  random_designs <- add_group_design_random(sampler$par_names[!sampler$nuisance], group_designs)
  random_maps <- add_group_design_map(sampler$par_names[!sampler$nuisance], group_designs)



  n_rand <- get_n_random(sampler$par_names[!sampler$nuisance], group_designs)
  has_random <- n_rand > 0

  M <- length(tmu) # Total parameters (former mu + beta)

  ## --------------------------------------------------
  ## 1) Update Beta (Fixed Effects)
  ## --------------------------------------------------
  # Target: alpha - Zu
  alpha_star <- alpha
  if (has_random && !is.null(u)) {
    # Compute Zu offset
    # Compute Zu offset manually relative to alpha to get residuals
    u_idx <- 1
    for (k in seq_len(p)) {
      pname <- names(fixed_designs)[k]
      if (!is.null(random_designs[[pname]])) {
        for (z in random_designs[[pname]]) {
          n_u <- ncol(z)
          offset <- as.vector(z %*% u[u_idx:(u_idx + n_u - 1)])
          alpha_star[k, ] <- alpha_star[k, ] - offset
          u_idx <- u_idx + n_u
        }
      }
    }
  }

  prec_data <- matrix(0, M, M)
  mean_data <- numeric(M)

  # For each subject i, we build an (p x M) matrix that maps tmu -> predicted alpha_i
  for (i in seq_len(n)) {
    par_idx <- 0
    M_i <- matrix(0, nrow = p, ncol = M)
    for (k in seq_len(p)) {
      x_ik <- fixed_designs[[k]][i, , drop = FALSE]
      # place them in row k:
      M_i[k, par_idx + 1:ncol(fixed_designs[[k]])] <- x_ik
      par_idx <- par_idx + ncol(fixed_designs[[k]])
    }
    # Accumulate
    prec_data <- prec_data + crossprod(M_i, tvinv %*% M_i)
    mean_data <- mean_data + crossprod(M_i, tvinv %*% alpha_star[, i, drop = FALSE])
  }

  prec_post <- prior$theta_mu_invar + prec_data
  cov_post <- solve(prec_post)
  mean_post <- cov_post %*% (prior$theta_mu_invar %*% prior$theta_mu_mean + mean_data)

  # Draw new parameter vector
  L <- t(chol(cov_post))
  z <- rnorm(M)
  tmu_new <- as.vector(mean_post + L %*% z)

  ## --------------------------------------------------
  ## 2) Update Random Effects u
  ## --------------------------------------------------
  u_new <- u
  s_new <- s
  a_half_s_new <- a_half_s

  if (has_random) {
    # Calculate mu_fixed = X * beta
    mu_fixed <- calculate_subject_means(fixed_designs, tmu_new)

    # We iterate over each random component
    u_ptr <- 1
    s_ptr <- 1

    for (k in seq_len(p)) { # For each parameter
      pname <- names(fixed_designs)[k]
      if (!is.null(random_designs[[pname]])) {
        # Prior precision in parameter space
        prec_noise <- tvinv[k, k]

        # Calculate total random effect contribution for this parameter
        total_rand_offset <- numeric(n)
        u_temp_ptr <- u_ptr
        for (z in random_designs[[pname]]) {
          n_u <- ncol(z)
          total_rand_offset <- total_rand_offset + as.vector(z %*% u_new[u_temp_ptr:(u_temp_ptr + n_u - 1)])
          u_temp_ptr <- u_temp_ptr + n_u
        }

        # Calculate cross-correlation adjustment offset
        # This accounts for correlations with other parameters in the blocked covariance structure
        mu_full <- calculate_subject_means_RE(fixed_designs, tmu_new, random_designs, u_new, random_maps)
        devs <- alpha - mu_full

        weighted_devs <- (tvinv %*% devs)[k, ]
        other_terms <- weighted_devs - tvinv[k, k] * devs[k, ]
        adjustment <- -other_terms / tvinv[k, k]

        # Base target for this parameter is the residual after removing fixed effects,
        # total random effects, and the cross-correlation adjustment.
        base_target <- alpha[k, ] - mu_fixed[k, ] - total_rand_offset - adjustment

        # Now iterate over components for THIS parameter k
        for (z in random_designs[[pname]]) {
          n_u <- ncol(z)
          curr_s2 <- s_new[s_ptr]

          current_u_segment <- u_new[u_ptr:(u_ptr + n_u - 1)]

          # Isolate the current random effect component (u_curr) from the residuals.
          # base_target excludes all random effects, so we add back the current
          # component (Z * u_segment) to get the target for this update step.
          target <- base_target + as.vector(z %*% current_u_segment)

          # Bayesian Linear Regression Update for u
          # Prior u ~ N(0, s^2), Error Precision = prec_noise
          ztz <- crossprod(z)
          # If z is indicator, ztz is diagonal.

          post_prec <- prec_noise * ztz
          diag(post_prec) <- diag(post_prec) + 1 / curr_s2

          post_cov <- solve(post_prec)

          post_mean <- post_cov %*% (prec_noise * crossprod(z, target))

          # Sample new u
          L_u <- t(chol(post_cov))
          u_drawn <- as.vector(post_mean + L_u %*% rnorm(n_u))

          # Update u_new
          u_new[u_ptr:(u_ptr + n_u - 1)] <- u_drawn

          # Update s for this component
          # s^2 | u ~ IG( v/2 + J/2, v/a_half_s + sum(u^2)/2 )
          # We sample the Variance (s^2) from the Inverse Gamma posterior.

          J_dim <- n_u
          shape_s <- (prior$v_Z + J_dim) / 2
          rate_s <- prior$v_Z / a_half_s_new[s_ptr] + 0.5 * sum(u_drawn^2)

          s2_new <- 1 / rgamma(1, shape = shape_s, rate = rate_s)

          # Store as Variance (s^2)
          s_new[s_ptr] <- s2_new

          # Update a_half_s (auxiliary variable for s variance)
          # Uses Inverse Gamma sampler (via 1/rgamma)
          # Selects A_z based on component index s_ptr
          val_A <- if (length(prior$A_z) >= s_ptr) prior$A_z[s_ptr] else prior$A_z[1]
          shape_a <- (prior$v_Z + 1) / 2

          a_half_s_new[s_ptr] <- 1 / rgamma(1, shape = shape_a, rate = prior$v_Z * (1 / s2_new) + 1 / (val_A^2))


          u_ptr <- u_ptr + n_u
          s_ptr <- s_ptr + 1
        }
      }
    }

    # Final update of subj_mu with NEW u for residuals
    subj_mu <- calculate_subject_means_RE(fixed_designs, tmu_new, random_designs, u_new, random_maps)
  } else {
    subj_mu <- calculate_subject_means(fixed_designs, tmu_new)
  }

  ## --------------------------------------------------
  ## 3) Update Covariance (tvar)
  ## --------------------------------------------------

  resid <- alpha - subj_mu

  tvar_new <- matrix(0, p, p)
  blocked_idx <- which(is_blocked)
  unblocked_idx <- which(!is_blocked)

  if (length(blocked_idx) > 0) {
    block_groups <- unique(par_group[blocked_idx])
    for (g in block_groups) {
      group_idx <- blocked_idx[par_group[blocked_idx] == g]
      d <- length(group_idx)
      resid_block <- resid[group_idx, , drop = FALSE]
      cov_block <- resid_block %*% t(resid_block)

      B_half_block <- 2 * prior$v * diag(1 / a_half[group_idx], d) + cov_block
      df_block <- prior$v + d - 1 + n

      Sigma_block <- riwish(df_block, B_half_block)
      tvar_new[group_idx, group_idx] <- Sigma_block
    }
  }

  if (length(unblocked_idx) > 0) {
    sum_sq_vec <- rowSums(resid[unblocked_idx, , drop = FALSE]^2)
    shape_vec <- prior$v / 2 + n / 2
    rate_vec <- prior$v / a_half[unblocked_idx] + 0.5 * sum_sq_vec

    tvinv_diag <- rgamma(length(unblocked_idx), shape = shape_vec, rate = rate_vec)
    tvar_diag <- 1 / tvinv_diag

    tvar_new[cbind(unblocked_idx, unblocked_idx)] <- tvar_diag
  }

  tvinv_new <- solve(tvar_new)

  block_dim <- integer(p)
  block_dim[unblocked_idx] <- 1
  if (length(blocked_idx) > 0) {
    block_groups <- unique(par_group[blocked_idx])
    for (g in block_groups) {
      group_idx <- blocked_idx[par_group[blocked_idx] == g]
      block_dim[group_idx] <- length(group_idx)
    }
  }

  shape_vec <- (prior$v + block_dim) / 2
  rate_vec <- prior$v * diag(tvinv_new) + 1 / (prior$A^2)
  a_half_new <- 1 / rgamma(p, shape = shape_vec, rate = rate_vec)

  out <- list(
    tmu = tmu_new,
    tvar = tvar_new,
    tvinv = tvinv_new,
    a_half = a_half_new,
    subj_mu = subj_mu,
    alpha = alpha
  )
  if (has_random) {
    out$u <- u_new
    out$s <- s_new
    out$a_half_s <- a_half_s_new
  }
  return(out)
}


last_sample_standard <- function(store) {
  out <- list(
    tmu = store$theta_mu[, store$idx], # Now includes beta
    tvar = store$theta_var[, , store$idx],
    tvinv = store$last_theta_var_inv,
    a_half = store$a_half[, store$idx]
  )
  if (!is.null(store$theta_u)) out$u <- store$theta_u[, store$idx]
  if (!is.null(store$theta_s)) {
    out$s <- store$theta_s[, store$idx]
    out$a_half_s <- store$a_half_s[, store$idx]
  }
  return(out)
}

get_conditionals_standard <- function(s, samples, n_pars, iteration = NULL, idx = NULL) {
  iteration <- ifelse(is.null(iteration), samples$iteration, iteration)
  if (is.null(idx)) idx <- 1:n_pars
  pts2_unwound <- apply(samples$theta_var[idx, idx, ], 3, unwind)
  all_samples <- rbind(samples$alpha[idx, s, ], samples$theta_mu[idx, ], pts2_unwound)
  mu_tilde <- rowMeans(all_samples)
  var_tilde <- stats::cov(t(all_samples))
  condmvn <- condMVN(
    mean = mu_tilde, sigma = var_tilde,
    dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
    X.given = c(samples$theta_mu[idx, iteration], unwind(samples$theta_var[idx, idx, iteration]))
  )
  return(list(eff_mu = condmvn$condMean, eff_var = condmvn$condVar))
}

unwind <- function(var_matrix, ...) {
  y <- t(chol(var_matrix))
  diag(y) <- log(diag(y))
  y[lower.tri(y, diag = TRUE)]
}

filtered_samples_standard <- function(sampler, filter, ...) {
  out <- list(
    theta_mu = sampler$samples$theta_mu[, filter, drop = F],
    theta_var = sampler$samples$theta_var[, , filter, drop = F],
    alpha = sampler$samples$alpha[, , filter, drop = F],
    iteration = length(filter)
  )
}


calc_log_jac_chol <- function(x) {
  ## work out dimension -------------------------------------------------
  n <- floor(sqrt(2 * length(x) + 0.25) - 0.5) # solves n(n+1)/2 = length(x)

  ## rebuild the raw lower‑triangular matrix ----------------------------
  mat <- matrix(NA_real_, n, n)
  mat[lower.tri(mat, diag = TRUE)] <- x # raw params (diag on log‑scale)

  ## step 1: create L ---------------------------------------------------
  L <- mat
  diag(L) <- exp(diag(mat)) # positive diagonal
  # (off‑diagonals stay on natural scale)

  ## log‑Jacobian -------------------------------------------------------
  logL <- diag(mat) # log of L_ii
  log_jac <- n * log(2) + sum((n - seq_len(n) + 2) * logL)

  return(log_jac)
}

unwind_chol <- function(x, reverse = FALSE) {
  if (reverse) {
    n <- sqrt(2 * length(x) + 0.25) - 0.5 ## Dim of matrix.
    out <- array(0, dim = c(n, n))
    out[lower.tri(out, diag = TRUE)] <- x
    diag(out) <- exp(diag(out))
    out <- out %*% t(out)
  } else {
    y <- t(base::chol(x))
    diag(y) <- log(diag(y))
    out <- y[lower.tri(y, diag = TRUE)]
  }
  return(out)
}

# bridge_sampling ---------------------------------------------------------
bridge_group_and_prior_and_jac_standard <- function(
    proposals_group, # (n_iter x [theta_mu, theta_a, var1, var2, u, s, a_s])
    proposals_list, # list of length n_subj, each item is (n_iter x p) alphas
    info) {
  # 'info' should contain:
  #   info$n_pars        = p        (# of rows in alpha)
  #   info$n_subjects    = n_subj   (# of subjects)
  #   info$is_blocked    = logical p  (which rows are blocked)
  #   info$par_group     = integer p  (group IDs for blocked rows)
  #   info$group_designs = list of length p, each (n_subj x m_k) [Full structure]
  #   info$prior:
  #       - $theta_mu_mean (p+B), $theta_mu_var ((p+B)x(p+B))
  #       - $v, $A (IW & gamma hyperparams)

  par_group <- info$par_group
  has_cov <- info$is_blocked
  block_groups <- unique(par_group[has_cov])
  prior <- info$prior

  p <- info$n_pars # dimension of each alpha_i
  n_subj <- info$n_subjects

  # Prepare designs
  # Fixed designs for mean calculation
  group_designs_fixed <- add_group_design(info$par_names, info$group_designs, n_subj)
  # Random designs for offset calculation
  random_designs <- add_group_design_random(info$par_names, info$group_designs)

  n_random <- get_n_random(info$par_names, info$group_designs)
  n_random_var <- get_n_random_variance(info$par_names, info$group_designs)

  # Total fixed parameters (including regressors)
  total_pars <- length(prior$theta_mu_mean)

  # Extract columns from proposals_group:
  # 1. theta_mu
  theta_mu <- proposals_group[, seq_len(total_pars), drop = FALSE]

  # 2. theta_a (log a_half)
  theta_a <- proposals_group[, total_pars + seq_len(p), drop = FALSE]
  min_idx <- total_pars + p

  # 3. theta_var1 (unblocked diagonals)
  if (any(!has_cov)) {
    theta_var1 <- proposals_group[, (min_idx + 1):(min_idx + sum(!has_cov)), drop = FALSE]
    min_idx <- min_idx + sum(!has_cov)
  }

  # 4. theta_var2 (blocked cholesky)
  if (any(has_cov)) {
    theta_var2_list <- list()
    for (block in block_groups) {
      cur_idx <- has_cov[par_group == block]
      n_cur <- sum(cur_idx)
      len_block <- (n_cur * (n_cur + 1)) / 2
      theta_var2_list[[block]] <- proposals_group[, (min_idx + 1):(min_idx + len_block), drop = FALSE]
      min_idx <- min_idx + len_block
    }
  }

  # 5. Random Effects (u, s, a_s)
  theta_u <- NULL
  theta_s <- NULL
  theta_a_s <- NULL

  if (n_random > 0) {
    # u (normal)
    theta_u <- proposals_group[, (min_idx + 1):(min_idx + n_random), drop = FALSE]
    min_idx <- min_idx + n_random

    # s (log SD)
    theta_s <- proposals_group[, (min_idx + 1):(min_idx + n_random_var), drop = FALSE]
    min_idx <- min_idx + n_random_var

    # a_half_s (log aux)
    theta_a_s <- proposals_group[, (min_idx + 1):(min_idx + n_random_var), drop = FALSE]
    min_idx <- min_idx + n_random_var
  }

  n_iter <- nrow(theta_mu)
  sum_out <- numeric(n_iter)

  # Precompute the vectorized prior on theta_mu
  prior_mu_log <- dmvnorm(theta_mu,
    mean  = prior$theta_mu_mean,
    sigma = prior$theta_mu_var,
    log   = TRUE
  )

  # Jacobians for var1, a
  if (any(!has_cov)) {
    jac_var1 <- rowSums(theta_var1)
  } else {
    jac_var1 <- 0
  }
  jac_a_vec <- rowSums(theta_a)

  # Jacobians for s and a_s (log transform)
  jac_s <- 0
  jac_a_s <- 0
  if (n_random > 0) {
    jac_s <- rowSums(theta_s)
    jac_a_s <- rowSums(theta_a_s)
  }

  #--- Main loop over MCMC draws ---
  for (i in seq_len(n_iter)) {
    # 1) Reconstruct partial-block var_curr
    var_curr <- matrix(0, p, p)

    # Unblocked diagonal => exp(theta_var1[i,])
    if (any(!has_cov)) {
      var_curr[!has_cov, !has_cov] <- diag(exp(theta_var1[i, ]), sum(!has_cov))
    }

    # Jacobian for var2
    jac_var2 <- 0
    prior_var2 <- 0
    # Blocked => unwind_chol
    if (any(has_cov)) {
      # Fill in theta_var2
      for (block in block_groups) {
        cur_idx <- par_group == block
        var_block <- unwind_chol(theta_var2_list[[block]][i, ], reverse = TRUE)
        var_curr[cur_idx, cur_idx] <- var_block
        # Prior influence
        d_block <- sum(cur_idx)
        prior_var2 <- prior_var2 + log(robust_diwish(
          W = var_block,
          v = prior$v + d_block - 1,
          S = 2 * prior$v * diag(1 / exp(theta_a[i, cur_idx]))
        ))
        # Jacobian influence
        jac_var2 <- jac_var2 + calc_log_jac_chol(theta_var2_list[[block]][i, ])
      }
    }

    # 2) Calculate Full Subject Means (incorporating u if present)
    u_vec <- NULL
    if (!is.null(theta_u)) u_vec <- theta_u[i, ]

    # Needs to be efficient: vectorized over subjects
    mu_full <- calculate_subject_means_RE(group_designs_fixed, theta_mu[i, ], random_designs, u_vec)

    # 3) Compute group likelihood = sum_{s=1..n_subj} dmvnorm(alpha_s, mu_s, var_curr)
    group_ll <- 0
    for (s in seq_len(n_subj)) {
      alpha_s <- proposals_list[[s]][i, ] # length p
      # mu_full is (p x n_subj)
      mu_s <- mu_full[, s]
      group_ll <- group_ll + dmvnorm(alpha_s, mu_s, var_curr, log = TRUE)
    }

    # 4) Prior on var1, var2, and a
    prior_var1 <- 0
    if (any(!has_cov)) {
      prior_var1 <- sum(logdinvGamma(
        x     = exp(theta_var1[i, ]),
        shape = prior$v / 2,
        rate  = prior$v / exp(theta_a[i, !has_cov, drop = FALSE])
      ))
    }

    prior_a <- sum(logdinvGamma(
      x     = exp(theta_a[i, ]),
      shape = 1 / 2,
      rate  = 1 / (prior$A^2)
    ))

    # 5) Prior on u, s, a_s
    prior_u <- 0
    prior_s <- 0
    prior_a_s <- 0

    if (n_random > 0) {
      # Current values
      s2_val <- exp(theta_s[i, ]) # vector of s^2 (Variance)
      a_s_val <- exp(theta_a_s[i, ]) # vector of a_half_s
      u_val <- theta_u[i, ]

      # prior_a_s: a_s ~ InverseGamma(1/2, 1/A_z^2)
      # calculated using logdinvGamma
      prior_a_s <- sum(logdinvGamma(
        x = a_s_val,
        shape = 1 / 2,
        rate = 1 / (prior$A_z^2)
      ))

      # prior_s: s^2 | a_s ~ InverseGamma(v_Z/2, v_Z/a_s)
      # s^2 is the parameter. We work with log(s^2).
      # The Jacobian for log(s^2) -> s^2 change of variable is handled by sum(theta_s).
      # So we calculate density of s^2 directly.
      prior_s <- sum(logdinvGamma(
        x = s2_val,
        shape = prior$v_Z / 2,
        rate = prior$v_Z / a_s_val
      ))

      # prior_u: u_k ~ N(0, s_k^2)
      # We need to map u indices to s indices.
      # Loop through parameters to handle mapping
      u_ptr <- 1
      s_ptr <- 1

      for (k in seq_along(info$par_names)) {
        pname <- info$par_names[k]
        if (!is.null(random_designs[[pname]])) {
          Z_parts <- random_designs[[pname]]
          for (j in seq_along(Z_parts)) {
            n_cols <- ncol(Z_parts[[j]])
            # Current s (SD) for this component
            curr_s <- sqrt(s2_val[s_ptr])

            # Current u chunk
            u_chunk <- u_val[u_ptr:(u_ptr + n_cols - 1)]

            # Density sum(dnorm(u, 0, s, log=T))
            prior_u <- prior_u + sum(dnorm(u_chunk, 0, curr_s, log = TRUE))

            u_ptr <- u_ptr + n_cols
            s_ptr <- s_ptr + 1
          }
        }
      }
    }

    # Combine
    sum_out[i] <- group_ll +
      prior_var1 + prior_var2 + jac_var2 + prior_a +
      prior_u + prior_s + prior_a_s
  } # end loop i

  # Add vectorized priors/Jacobians
  sum_out <- sum_out + prior_mu_log + jac_var1 + jac_a_vec + jac_s + jac_a_s

  return(sum_out)
}


bridge_add_info_standard <- function(info, samples) {
  if (is.null(samples$is_blocked)) {
    has_cov <- rep(T, samples$n_pars)
  } else {
    has_cov <- samples$is_blocked
  }
  info$par_group <- samples$par_group
  if (is.null(info$par_group)) {
    info$par_group <- rep(1, samples$n_pars)
  }
  info$is_blocked <- has_cov
  info$group_designs <- samples$group_designs
  info$par_names <- samples$par_names

  # Calculate total parameters (including former beta)
  n_total_pars <- nrow(samples$samples$theta_mu)

  # Base counts
  n_mu <- n_total_pars
  n_a <- samples$n_pars
  n_var <- sum(!has_cov) + (sum(has_cov) * (sum(has_cov) + 1)) / 2

  # Random effects counts
  n_u <- get_n_random(samples$par_names, samples$group_designs)
  n_s <- get_n_random_variance(samples$par_names, samples$group_designs)
  # a_half_s has same dim as s

  total_group_pars <- n_mu + n_a + n_var + n_u + n_s + n_s

  # Calculate group index including all parameters (theta_mu now includes beta, plus random effects)
  start_idx <- samples$n_pars * samples$n_subjects + 1
  info$group_idx <- start_idx:(start_idx + total_group_pars - 1)

  return(info)
}

bridge_add_group_standard <- function(all_samples, samples, idx) {
  # Add theta_mu (now includes all parameters)
  all_samples <- cbind(all_samples, t(samples$samples$theta_mu[, idx]))
  all_samples <- cbind(all_samples, t(log(samples$samples$a_half[, idx])))

  # Handle variance matrices
  par_group <- samples$par_group
  has_cov <- samples$is_blocked
  if (is.null(has_cov)) {
    has_cov <- rep(T, samples$n_pars)
  }
  if (is.null(par_group)) {
    par_group <- rep(1, samples$n_pars)
  }

  if (any(!has_cov)) {
    all_samples <- cbind(all_samples, t(log(matrix(apply(samples$samples$theta_var[!has_cov, !has_cov, idx, drop = F], 3, diag), ncol = nrow(all_samples)))))
  }
  if (any(has_cov)) {
    for (block in unique(par_group[has_cov])) {
      cur_idx <- par_group == block
      all_samples <- cbind(all_samples, t(matrix(apply(samples$samples$theta_var[cur_idx, cur_idx, idx, drop = F], 3, unwind_chol), ncol = nrow(all_samples))))
    }
  }

  # Add Random Effects u, s, a_half_s
  if (!is.null(samples$samples$theta_u)) {
    all_samples <- cbind(all_samples, t(samples$samples$theta_u[, idx, drop = FALSE]))
    # s (stored as SD)
    all_samples <- cbind(all_samples, t(log(samples$samples$theta_s[, idx, drop = FALSE])))
    # a_half_s
    all_samples <- cbind(all_samples, t(log(samples$samples$a_half_s[, idx, drop = FALSE])))
  }

  return(all_samples)
}

# for IC ------------------------------------------------------------------

group__IC_standard <- function(emc, stage = "sample", filter = NULL) {
  # 1) Retrieve draws
  alpha <- get_pars(emc,
    selection = "alpha", stage = stage, filter = filter,
    return_mcmc = FALSE, merge_chains = TRUE
  )
  theta_mu <- get_pars(emc,
    selection = "mu", stage = stage, filter = filter,
    return_mcmc = FALSE, merge_chains = TRUE
  )
  theta_var <- get_pars(emc,
    selection = "Sigma", stage = stage, filter = filter,
    return_mcmc = FALSE, merge_chains = TRUE, remove_constants = F
  )

  p <- dim(alpha)[1] # number of subject-level parameters
  n_subj <- dim(alpha)[2] # number of subjects
  N <- dim(alpha)[3] # number of samples

  # 2) Averages
  mean_alpha <- apply(alpha, c(1, 2), mean)
  mean_mu <- rowMeans(theta_mu)
  mean_var <- apply(theta_var, c(1, 2), mean)

  if (is.null(emc[[1]]$group_designs)) {
    lls <- numeric(N)
    for (i in 1:N) {
      lls[i] <- sum(dmvnorm(t(alpha[, , i]), theta_mu[, i], theta_var[, , i], log = T))
    }
    mean_pars_ll <- sum(dmvnorm(t(mean_alpha), mean_mu, mean_var, log = TRUE))
  } else {
    # 3) Build design matrices
    group_designs <- add_group_design(emc[[1]]$par_names, emc[[1]]$group_designs, n_subj)
    # 4) Per-draw log-likelihood
    lls <- standard_subj_ll(group_designs, theta_var, theta_mu, alpha, n_subj)
    # 5) Likelihood at posterior mean
    # Calculate subject-level means using the mean parameters
    mean_subj_means <- calculate_subject_means(group_designs, mean_mu)
    mean_pars_ll <- 0
    for (s in seq_len(n_subj)) {
      mean_pars_ll <- mean_pars_ll +
        dmvnorm(mean_alpha[, s], mean_subj_means[, s], mean_var, log = TRUE)
    }
  }

  minD <- -2 * max(lls)
  mean_ll <- mean(lls)
  Dmean <- -2 * mean_pars_ll
  list(mean_ll = mean_ll, Dmean = Dmean, minD = minD)
}

standard_subj_ll <- function(theta_var, theta_mu, alpha, n_subj, group_designs) {
  N <- dim(theta_var)[3]
  p <- nrow(theta_mu)
  log2pi <- log(2 * pi)

  # pre‑allocate
  ll <- numeric(N)

  for (i in seq_len(N)) {
    Sigma <- theta_var[, , i]
    U <- chol(Sigma) # will error if non‑PD
    rooti <- backsolve(U, diag(p))
    log_const <- sum(log(diag(rooti))) - 0.5 * p * log2pi

    mu <- calculate_subject_means(group_designs, theta_mu[, i])
    z <- t(alpha[, , i] - mu) %*% rooti # n_subj × p
    ll[i] <- -0.5 * sum(z^2) + n_subj * log_const
  }
  ll
}
