
# PMwG --------------------------------------------------------------------


add_info_lm <- function(sampler, prior = NULL, ...){
  args <- list(...)
  # Checking and default priors
  if (is.null(prior)) { # deliberately high
    prior <- list(theta_mu_mean = rep(0, sampler$n_pars), theta_mu_var = rep(100000, sampler$n_pars))
  }
  # Things I save rather than re-compute inside the loops.
  prior$theta_mu_invar <- 1/prior$theta_mu_var #Inverse of the prior

  #Hyper parameters
  attr(sampler, "a_g") <- 1/2
  attr(sampler, "b_g") <- 1/2
  attr(sampler, "a_0") <- 0.001
  attr(sampler, "b_0") <- 0.001
  attr(sampler, "formula") <- args$formula
  sampler$prior <- prior
  sampler$aggr_data <- args$aggr_data
  effects <- get_effects(args$aggr_data, args$formula)
  sampler$n_effects <- length(effects$effect_mapping)
  attr(sampler, "effect_mapping") <- effects$effect_mapping
  attr(sampler, "effect_grouping") <- effects$effect_grouping
  attr(sampler, "effect_types") <- effects$effect_types
  attr(sampler, "dms") <- effects$dms
  return(sampler)
}

get_effects <- function(aggr_data, formula){
  n_effects <- 0
  effect_mapping <- numeric(0)
  effect_grouping <- c(0)
  effect_types <- numeric(0)
  dms <- list()
  for(form in formula){
    vars <- split_form(form)
    # Find out which are factors and set contrasts
    factor_vars <- names(Filter(is.factor, aggr_data[,vars$dep, drop = F]))
    contrasts <- replicate(length(factor_vars), contr.bayes)
    names(contrasts) <- factor_vars
    # drop independent variable
    newform <- update(form, NULL ~ .)
    if(!is.null(factor_vars)){
      m_matrix <- model.matrix(newform, aggr_data, contrasts.arg = contrasts)
    } else{
      m_matrix <- model.matrix(newform, aggr_data)
    }
    dms[[vars$ind]] <- m_matrix[,-1, drop =F]
    groups <- attr(m_matrix, "assign")[-1]
    effect_grouping <- c(effect_grouping, max(effect_grouping) + groups)
    effect_mapping <- c(effect_mapping, rep(vars$ind, length(groups)))
    effect_types <- c(effect_types, get_effect_types(form, m_matrix, vars$dep %in% factor_vars))
  }
  return(list(effect_grouping = effect_grouping[-1], effect_mapping = effect_mapping,
              effect_types = effect_types, dms = dms))
}

get_effect_types <- function(form, m_matrix, is_factor){
  all_vars <- split_form(form)$dep # remove dependent variable
  effect_names <- colnames(m_matrix)[-1] # remove intercept
  matches <- matrix(NA_integer_, nrow = length(all_vars), ncol = length(effect_names))
  for(i in 1:length(all_vars)){
    matches[i,grepl(all_vars[i], effect_names)] <- is_factor[i]
  }
  effect_types <- numeric(length(effect_names))
  for(j in 1:ncol(matches)){
    idx <- !is.na(matches[,j])
    if(sum(matches[idx, j]) < sum(idx)){
      effect_types[j] <- 0
    } else{
      effect_types[j] <- 1
    }
  }
  return(effect_types)
}

sample_store_lm <- function(data, par_names, iters = 1, stage = "init", integrate = T, ...) {
  args <- list(...)
  effects <- get_effects(args$aggr_data, args$formula)
  n_effects <- length(effects$effect_mapping)
  subject_ids <- unique(data$subjects)
  n_pars <- length(par_names)
  n_subjects <- length(subject_ids)
  base_samples <- sample_store_base(data, par_names, iters, stage)
  samples <- list(
    theta_mu = array(NA_real_,dim = c(n_pars, iters), dimnames = list(par_names, NULL)),
    theta_var = array(NA_real_,dim = c(n_pars, n_pars, iters),dimnames = list(par_names, par_names, NULL)),
    theta_theta = array(NA_real_, dim = c(n_effects, iters), dimnames = list(NULL, NULL)),
    theta_g = array(NA_real_, dim = c(n_effects, iters), dimnames = list(NULL, NULL)),
    theta_g_untr = array(NA_real_, dim = c(n_effects, iters), dimnames = list(NULL, NULL))
  )
  if(integrate) samples <- c(samples, base_samples)
  return(samples)
}

get_startpoints_lm <- function(pmwgs, start_mu, start_var){
  if (is.null(start_mu)) start_mu <- rnorm(pmwgs$n_pars, mean = pmwgs$prior$theta_mu_mean, sd = 1)
  # If no starting point for group var just sample some
  if (is.null(start_var)) start_var <- diag(1/rgamma(pmwgs$n_pars, 10, 5)) #Bit stupid maybe as startpoint
  start_theta <- rnorm(pmwgs$n_effects)
  start_g <- rep(1, pmwgs$n_effects)
  sub_mu <- matrix(rep(start_mu, pmwgs$n_subjects), ncol = pmwgs$n_subjects)
  return(list(tmu = start_mu, tvar = start_var, theta = start_theta,
              g = start_g, g_untr = start_g, sub_mu = sub_mu))
}

fill_samples_lm <- function(samples, group_level, proposals, epsilon, j = 1, n_pars){
  samples$theta_theta[,j] <- group_level$theta
  samples$theta_g[,j] <- group_level$g
  samples$theta_g_untr[,j] <- group_level$g_untr
  samples <- fill_samples_base(samples, group_level, proposals, epsilon, j = j, n_pars)
  return(samples)
}

gibbs_step_lm <- function(sampler, alpha){
  last <- last_sample_lm(sampler$samples)
  hyper <- attributes(sampler)
  prior <- sampler$prior

  mu <- last$mu
  sigma2 <- diag(last$var)
  theta <- last$theta
  g <- last$g
  g_untr <- last$g
  sub_mu <- matrix(0, nrow = nrow(alpha), ncol = ncol(alpha))
  for(i in 1:sampler$n_pars){
    y <- alpha[i,]
    X <- hyper$dms[[rownames(alpha)[i]]]
    idx <- attr(sampler, "effect_mapping") == rownames(alpha)[i]

    if(is.null(X)){
      Xtheta <- 0
      theta_g_theta <- 0
    } else{
      Xtheta <- X %*% theta[idx]
      theta_g_theta <- t(theta[idx]) %*% diag(1/g[idx], sum(idx)) %*% theta[idx]
    }
    # mu
    mu_sigma2 <- 1/(sampler$n_subjects/sigma2[i] + 1/prior$theta_mu_var[i])
    mu_mu <- mu_sigma2*((sum(y) - Xtheta)/sigma2[i] + prior$theta_mu_mean[i]/prior$theta_mu_var[i])
    mu[i] <- rnorm(1, mu_mu, sqrt(mu_sigma2))

    # sigma
    sigma_shape <- hyper$a_0 + sampler$n_subjects/2 + sum(idx)/2 #sum(idx) = number of effects for this parameter
    sigma_rate <- hyper$b_0 + (sum((y - mu[i] - Xtheta)^2) +theta_g_theta)/2
    sigma2[i] <- 1/rgamma(1, sigma_shape, sigma_rate)

    if(!is.null(X)){
      # Theta
      theta_var <- solve(t(X) %*% X + diag(1/g[idx], nrow = sum(idx)))
      theta_mu <- theta_var %*% t(t(y) %*% X - t(rep(mu[i], sampler$n_subjects)) %*% X)
      theta[idx] <- as.vector(mvtnorm::rmvnorm(1, theta_mu, theta_var))

      # G
      effect_groups <- attr(sampler, "effect_grouping")
      effect_types <- attr(sampler, "effect_types")
      for(group in unique(effect_groups[idx])){
        group_idx <- effect_groups == group
        factor_type <- effect_types[group_idx][1]
        if(factor_type){
          # G
          g_shape <- hyper$a_g + 1/2
          g_rate <- hyper$b_g + (theta[group_idx]^2)/(2*sigma2[i])
          g[group_idx] <- rep(1/rgamma(1, g_shape, g_rate), sum(group_idx))
          g_untr[group_idx] <- g[group_idx]
        } else{
          # G
          g_shape <- hyper$a_g + 1/2
          g_rate <- sampler$n_subjects * hyper$b_g + (theta[group_idx]^2)/(2*sigma2[i])
          g_untr[group_idx] <-  rep(1/rgamma(1, g_shape, g_rate), sum(group_idx))
          g[group_idx] <- g_untr[group_idx] %*% solve(t(X[,which(group_idx[idx])]) %*% X[,which(group_idx[idx])])
        }
      }
    }
    if(!is.null(X)){
      sub_mu[i,] <- mu[i] + X %*% theta[idx]
    } else{
      sub_mu[i,] <- mu[i]
    }
  }
  return(list(tmu = mu, tvar = diag(sigma2), theta = theta, g = g, g_untr = g_untr, sub_mu = sub_mu, alpha = alpha))
}

last_sample_lm <- function(store) {
  list(
    mu = store$theta_mu[, store$idx],
    var = store$theta_var[,,store$idx],
    theta = store$theta_theta[,store$idx],
    g = store$theta_g[,store$idx]
  )
}

get_group_level_lm <- function(parameters, s){
  mu <- parameters$sub_mu[,s]
  var <- parameters$tvar
  return(list(mu = mu, var = var))
}

get_conditionals_lm <- function(s, samples, n_pars, iteration = idx = NULL){
  iteration <- ifelse(is.null(iteration), samples$iteration, iteration)
  if(is.null(idx)) idx <- 1:n_pars
  theta_var <-log(apply(samples$theta_var[idx,idx,],3,diag))
  theta_mu <- samples$theta_mu[idx]
  theta_theta <- samples$theta_theta
  all_samples <- rbind(samples$alpha[, s,],theta_mu, theta_var, theta_theta)
  mu_tilde <- rowMeans(all_samples)
  var_tilde <- cov(t(all_samples))
  condmvn <- condMVN(mean = mu_tilde, sigma = var_tilde,
                     dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                     X.given = c(theta_mu[idx,iteration],
                                 log(diag(samples$theta_var[idx,idx,iteration])),
                                 theta_theta[,iteration]))
  return(list(eff_mu = condmvn$condMean, eff_var = condmvn$condVar))
}

filtered_samples_lm <- function(sampler, filter){
  out <- list(
    theta_mu = sampler$samples$theta_mu[, filter],
    theta_var = sampler$samples$theta_var[, , filter],
    theta_theta = sampler$samples$theta_theta[, filter, drop = F],
    alpha = sampler$samples$alpha[, , filter],
    iteration = length(filter)
  )
}

contr.bayes <- function(n, contrasts = TRUE) {
  if (length(n) <= 1L) {
    if (is.numeric(n) && length(n) == 1L && n > 1L)
      TRUE
    else stop("not enough degrees of freedom to define contrasts")
  } else n <- length(n)
  cont <- diag(n)
  if (contrasts) {
    a <- n
    I_a <- diag(a)
    J_a <- matrix(1, nrow = a, ncol = a)
    Sigma_a <- I_a - J_a/a
    cont <- eigen(Sigma_a)$vectors[,seq_len(a-1), drop = FALSE]
  }
  cont
}

# IS2 ---------------------------------------------------------------------



get_all_pars_lm <- function(samples, idx, info){
  n_subjects <- samples$n_subjects
  n_iter <- length(samples$samples$stage[idx])
  info$n_effects <- samples$n_effects
  info$effect_mapping <- attr(samples, "effect_mapping")
  info$effect_grouping <- attr(samples, "effect_grouping")
  info$effect_types <- attr(samples, "effect_types")
  info$dms <- attr(samples, "dms")

  # Extract relevant objects
  alpha <- samples$samples$alpha[,,idx]
  theta_mu <- samples$samples$theta_mu[,idx]
  theta_theta <- samples$samples$theta_theta[,idx, drop = F]
  theta_var <- samples$samples$theta_var[,,idx, drop = F]
  theta_g <- log(samples$samples$theta_g[,idx, drop = F])
  theta_g_untr <- log(samples$samples$theta_g_untr[,idx, drop = F])
  theta_g[!info$effect_types,] <- theta_g_untr[!info$effect_types,]
  theta_g <- theta_g[!duplicated(info$effect_grouping),]
  theta_var.unwound <- log(apply(theta_var,3,diag))
  # Set up
  n_params<- nrow(alpha) + nrow(theta_mu) + nrow(theta_theta) + nrow(theta_var.unwound)
  all_samples=array(dim=c(n_subjects,n_params,n_iter))
  mu_tilde=array(dim = c(n_subjects,n_params))
  var_tilde=array(dim = c(n_subjects,n_params,n_params))

  for (j in 1:n_subjects){
    all_samples[j,,] = rbind(alpha[,j,],theta_mu[,],theta_theta[,],theta_var.unwound[,])
    # calculate the mean for re, mu and sigma
    mu_tilde[j,] =apply(all_samples[j,,],1,mean)
    # calculate the covariance matrix for random effects, mu and sigma
    var_tilde[j,,] = cov(t(all_samples[j,,]))
  }

  for(i in 1:n_subjects){ #RJI_change: this bit makes sure that the sigma tilde is pos def
    if(!corpcor::is.positive.definite(var_tilde[i,,], tol=1e-8)){
      var_tilde[i,,]<-corpcor::make.positive.definite(var_tilde[i,,], tol=1e-6)
    }
  }

  X <- cbind(t(theta_mu),t(theta_theta),t(theta_g), t(theta_var.unwound))
  info$n_params <- n_params
  info$given.ind <- (info$n_randeffect+1):n_params
  info$X.given_ind <- 1:length(info$given.ind)
  return(list(X = X, mu_tilde = mu_tilde, var_tilde = var_tilde, info = info))
}

group_dist_lm = function(random_effect = NULL, parameters, sample = FALSE, n_samples = NULL, info){
  n_randeffect <- info$n_randeffect
  par_names <- info$par_names
  param.theta_mu <- parameters[1:n_randeffect]
  param.theta_theta <- parameters[(n_randeffect + 1):(n_randeffect + info$n_effects)]
  for(i in 1:n_randeffect){
    X <- info$dms[[par_names[i]]]
    if(!is.null(X)){
      idx <- info$effect_mapping == par_names[i]
      param.theta_mu[i] <- param.theta_mu[i] + sum(X[info$sub,] * param.theta_theta[idx])
    }
  }
  param.var <- parameters[(length(parameters) - n_randeffect + 1):length(parameters)]
  if (sample){
    return(matrix(rnorm(n_samples*length(param.theta_mu), param.theta_mu, sqrt(exp(param.var))),
                  ncol = length(param.theta_mu),byrow = T))
  }
  else{
    logw_second<-max(-5000*info$n_randeffect, sum(dnorm(t(random_effect), param.theta_mu, sqrt(exp(param.var)), log = T)))
    return(logw_second)
  }
}

prior_dist_lm = function(parameters, info){
  n_randeffect <- info$n_randeffect
  n_effects <- info$n_effects
  prior <- info$prior
  hyper <- info$hyper
  #Extract and when necessary transform back
  param.theta_mu <- parameters[1:n_randeffect]
  if(info$n_effects > 0){
    param.theta_theta <- exp(parameters[(n_randeffect+1):(n_randeffect + n_effects)])
    param.theta_g <- exp(parameters[(n_randeffect+n_effects+1):(n_randeffect + n_effects + max(info$effect_grouping))])
  }

  param.theta_var <- exp(parameters[(length(parameters) - n_randeffect + 1):length(parameters)])

  names(param.theta_var) <- info$par_names
  log_prior_mu <- sum(dnorm(param.theta_mu, mean = prior$theta_mu_mean, sd = sqrt(prior$theta_mu_var), log =TRUE))
  log_prior_var <- sum(logdinvGamma(param.theta_var, shape = hyper$a_0, rate = hyper$b_0))
  if(info$n_effects > 0){
    log_prior_theta <- sum(dnorm(param.theta_theta, mean = 0, sd = sqrt(param.theta_g[info$effect_grouping]
                                                                        *param.theta_var[info$effect_mapping])))
  } else{
    log_prior_theta <- 0
  }
  log_prior_g_factor <- 0
  log_prior_g_regr <- 0
  jac_g_factor <- 0
  jac_g_regr <- 0
  if(info$n_effects > 0){
    for(is_factor in unique(info$effect_types)){
      if(is_factor){
        log_prior_g_factor <- sum(logdinvGamma(param.theta_g[info$effect_types], shape = 1/2, rate = hyper$b_g/2))
        jac_g_factor <- -sum(log(param.theta_g[info$effect_types[!duplicated(info$effect_grouping)]]))
      } else{
        log_prior_g_regr <- sum(logdinvGamma(param.theta_g[!info$effect_types[!duplicated(info$effect_grouping)]], shape = 1/2, rate = (info$n_subjects*hyper$b_g)/2))
        jac_g_factor <- -sum(log(param.theta_g[!info$effect_types[!duplicated(info$effect_grouping)]]))
      }
    }
  }
  log_prior_g <- log_prior_g_factor + log_prior_g_regr
  jac_g <- jac_g_factor + jac_g_regr
  jac_var <- -sum(log(param.theta_var)) # Jacobian determinant of transformation of log of the var
  # Jacobians are actually part of the denominator (dnorm(prop_theta)) since transformations of the data (rather than parameters),
  # warrant a jacobian added. But we add the jacobians here for ease of calculations.
  return(log_prior_mu + log_prior_var + log_prior_theta + log_prior_g - jac_var - jac_g)
}

split_form <- function(formula) {
  tt <- terms(formula)
  vars <- as.character(attr(tt, "variables"))[-1] ## [1] is the list call
  response <- attr(tt, "response") # index of response var
  return(list(dep = vars[-response], ind = vars[response]))
}

