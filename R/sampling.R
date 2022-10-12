# type=c("standard","diagonal","blocked","factor","factorRegression","single")[1]
# n_chains=3; rt_resolution=0.02
# prior_list = NULL;par_groups=NULL;n_factors=NULL;constraintMat = NULL;covariates=NULL
# data_list=miletic1_rdm_simdat; design_list=design;model_list=NULL; rt_resolution=.001
make_samplers <- function(data_list,design_list,model_list=NULL,
                          type=c("standard","diagonal","blocked","factor","factorRegression","single")[1],
                          n_chains=3,rt_resolution=0.02,
                          prior_list = NULL,
                          par_groups=NULL,
                          subject_covariates = NULL,
                          n_factors=NULL,constraintMat = NULL,covariates=NULL)

{
  if (!(type %in% c("standard","diagonal","blocked","factor","factorRegression","single")))
    stop("type must be one of: standard,diagonal,blocked,factor,factorRegression,single")
  if (class(data_list)=="data.frame") data_list <- list(data_list)
  # Sort subject together and add unique trial within subject integer
  # create overarching data list with one list per subject
  data_list <- lapply(data_list,function(d){
    if (!is.factor(d$subjects)) d$subjects <- factor(d$subjects)
    d <- d[order(d$subjects),]
    add_trials(d)
  })
  if (!is.null(names(design_list)[1]) && names(design_list)[1]=="Flist")
    design_list <- list(design_list)
  if (length(design_list)!=length(data_list))
    design_list <- rep(design_list,length(data_list))
  if (is.null(model_list)) model_list <- lapply(design_list,function(x){x$model})
  if (any(unlist(lapply(model_list,is.null))))
    stop("Must supply model_list if model is not in all design_list components")
  if (!is.null(names(model_list)[1]) && names(model_list)[1]=="type")
    model_list <- list(model_list)
  if (length(model_list)!=length(data_list))
    model_list <- rep(model_list,length(data_list))
  if (!is.null(names(prior_list)) && any(names(prior_list)=="theta_mu_mean"))
    prior_list <- list(prior_list)
  if (length(prior_list)!=length(data_list))
    prior_list <- rep(prior_list,length(data_list))
  dadm_list <- vector(mode="list",length=length(data_list))
  rt_resolution <- rep(rt_resolution,length.out=length(data_list))
  for (i in 1:length(dadm_list)) {
    message("Processing data set ",i)
    # if (!is.null(design_list[[i]]$Ffunctions)) {
    #   pars <- attr(data_list[[i]],"pars")
    #   data_list[[i]] <- cbind.data.frame(data_list[[i]],data.frame(lapply(
    #     design_list[[i]]$Ffunctions,function(f){f(data_list[[i]])})))
    #   if (!is.null(pars)) attr(data_list[[i]],"pars") <- pars
    # }
    # create a design model
    dadm_list[[i]] <- design_model(data=data_list[[i]],design=design_list[[i]],
                                   model=model_list[[i]],rt_resolution=rt_resolution[i],prior=prior_list[[i]])
  }
  if(!is.null(subject_covariates)) attr(dadm_list, "subject_covariates") <- subject_covariates
  if (type == "standard") {
    variant_funs <- get_variant_funs(type = "standard")
    out <- pmwgs(dadm_list, variant_funs)
  } else if(type == "diagonal"){
    source("samplers/pmwg/variants/diag.R")
    out <- pmwgs(dadm_list)
    out$source <- "samplers/pmwg/variants/diag.R"
  } else if (type == "blocked") {
    if (is.null(par_groups)) stop("Must specify par_groups for blocked type")
    source("samplers/pmwg/variants/blocked.R")
    out <- pmwgs(dadm_list,par_groups=par_groups)
    out$source <- "samplers/pmwg/variants/blocked.R"
  } else if (type=="factor") {
    if (is.null(n_factors)) stop("Must specify n_factors for factor type")
    source("samplers/pmwg/variants/factor.R")
    out <- pmwgs(dadm_list,n_factors=n_factors,constraintMat=constraintMat)
    out$source <- "samplers/pmwg/variants/factor.R"
  } else if (type=="factorRegression") {
    if (is.null(n_factors)) stop("Must specify n_factors for factorRegression type")
    if (is.null(covariates)) stop("Must specify covariates for factorRegression type")
    source("samplers/pmwg/variants/factorRegr.R")
    out <- pmwgs(dadm_list,n_factors=n_factors,constraintMat=constraintMat,
                 covariates=covariates)
    out$source <- "samplers/pmwg/variants/factorRegr.R"
  } else if (type == "single") {
    source("samplers/pmwg/variants/single.R")
    out <- pmwgs(dadm_list)
    out$source <- "samplers/pmwg/variants/single.R"
  }
  # replicate chains
  dadm_lists <- rep(list(out),n_chains)
  # For post predict
  attr(dadm_lists,"data_list") <- data_list
  attr(dadm_lists,"design_list") <- design_list
  attr(dadm_lists,"model_list") <- model_list
  return(dadm_lists)
}

# data generation

# Used in make_data and make_samplers
add_trials <- function(dat)
  # Add trials column, 1:n for each subject
{
  n <- table(dat$subjects)
  if (!any(names(dat)=="trials")) dat <- cbind.data.frame(dat,trials=NA)
  for (i in names(n)) dat$trials[dat$subjects==i] <- 1:n[i]
  dat
}

pmwgs <- function(dadm, variant_funs, pars = NULL, ll_func = NULL, prior = NULL, ...) {
  if(is.data.frame(dadm)) dadm <- list(dadm)
  dadm <- extractDadms(dadm)
  if(is.null(pars)) pars <- dadm$pars
  if(is.null(ll_func)) ll_func <- dadm$ll_func
  if(is.null(prior)) prior <- dadm$prior
  dadm_list <-dadm$dadm_list

  # Storage for the samples.
  subjects <- sort(as.numeric(unique(dadm$subjects)))
  samples <- variant_funs$sample_store(dadm, pars)
  sampler <- list(
    data = dadm_list,
    par_names = c(pars, names(dadm$subject_covariates)),
    subjects = subjects,
    n_pars = length(pars) + length(dadm$subject_covariates),
    subject_covariates = dadm$subject_covariates,
    n_subjects = length(subjects),
    ll_func = ll_func,
    samples = samples,
    init = FALSE
  )
  class(sampler) <- "pmwgs"
  sampler <- variant_funs$add_info(sampler, prior)
  attr(sampler, "variant_funs") <- variant_funs
  return(sampler)
}

extractDadms <- function(dadms, names = 1:length(dadms)){
  N_models <- length(dadms)
  subject_covariates <- attr(dadms, "subject_covariates")
  pars <- attr(dadms[[1]], "sampled_p_names")
  prior <- attr(dadms[[1]], "prior")
  ll_func <- attr(dadms[[1]], "model")$log_likelihood
  subjects <- unique(factor(sapply(dadms, FUN = function(x) levels(x$subjects))))
  dadm_list <- dm_list(dadms[[1]])
  components <- length(pars)
  if(N_models > 1){
    total_dadm_list <- vector("list", length = N_models)
    k <- 1
    pars <- paste(names[1], pars, sep = "|")
    dadm_list[as.character(which(!subjects %in% unique(dadms[[1]]$subjects)))] <- NA
    total_dadm_list[[1]] <- dadm_list
    for(dadm in dadms[-1]){
      k <- k + 1
      tmp_list <- vector("list", length = length(subjects))
      tmp_list[as.numeric(unique(dadm$subjects))] <- dm_list(dadm)
      total_dadm_list[[k]] <- tmp_list
      curr_pars <- attr(dadm, "sampled_p_names")
      components <- c(components, max(components) + length(curr_pars))
      pars <- c(pars, paste(names[k], curr_pars, sep = "|"))
      prior$theta_mu_mean <- c(prior$theta_mu_mean, attr(dadm, "prior")$theta_mu_mean)
      if(is.matrix(prior$theta_mu_var)){
        prior$theta_mu_var <- adiag(prior$theta_mu_var, attr(dadm, "prior")$theta_mu_var)
      } else{
        prior$theta_mu_var <- c(prior$theta_mu_var, attr(dadm, "prior")$theta_mu_var)
      }
    }
    ll_func <- jointLL
    dadm_list <- do.call(mapply, c(list, total_dadm_list, SIMPLIFY = F))
  }
  subject_covariates_ok <- unlist(lapply(subject_covariates, FUN = function(x) length(x) == length(subjects)))
  if(!is.null(subject_covariates_ok)) if(any(!subject_covariates_ok)) stop("subject_covariates must be as long as the number of subjects")
  attr(dadm_list, "components") <- components
  return(list(ll_func = ll_func, pars = pars, prior = prior,
              dadm_list = dadm_list, subjects = subjects, subject_covariates = subject_covariates))
}


# some sort of variants, used to be called via a list called variant_funs

get_variant_funs <- function(type = "standard") {
  if(type == "standard") {
    list_fun <- list(# store functions
      sample_store = sample_store_standard,
      add_info = add_info_standard,
      get_startpoints = get_startpoints_standard,
      get_group_level = get_group_level_standard,
      fill_samples = fill_samples_standard,
      gibbs_step = gibbs_step_standard,
      filtered_samples = filtered_samples_standard,
      get_conditionals = get_conditionals_standard
    )
  }
  return(list_fun)
}


sample_store_standard <- function(data, par_names, iters = 1, stage = "init", integrate = T, ...) {
  subject_ids <- unique(data$subjects)
  n_pars <- length(par_names)
  n_subjects <- length(subject_ids)
  base_samples <- sample_store_base(data, par_names, iters, stage)
  samples <- list(
    theta_mu = array(NA_real_,dim = c(n_pars, iters), dimnames = list(par_names, NULL)),
    theta_var = array(NA_real_,dim = c(n_pars, n_pars, iters),dimnames = list(par_names, par_names, NULL)),
    a_half = array(NA_real_,dim = c(n_pars, iters),dimnames = list(par_names, NULL))
  )
  if(integrate) samples <- c(samples, base_samples)
  return(samples)
}

add_info_standard <- function(sampler, prior = NULL, ...){
  # Checking and default priors
  if (is.null(prior)) {
    prior <- list(theta_mu_mean = rep(0, sampler$n_pars), theta_mu_var = diag(rep(1, sampler$n_pars)))
  }
  # Things I save rather than re-compute inside the loops.
  prior$theta_mu_invar <- MASS::ginv(prior$theta_mu_var) #Inverse of the matrix

  #Hyper parameters
  attr(sampler, "v_half") <- 2
  attr(sampler, "A_half") <- 1
  sampler$prior <- prior
  return(sampler)
}

get_startpoints_standard <- function(pmwgs, start_mu, start_var){
  if (is.null(start_mu)) start_mu <- mvtnorm::rmvnorm(1, mean = pmwgs$prior$theta_mu_mean, sigma = pmwgs$prior$theta_mu_var)
  # If no starting point for group var just sample some
  if (is.null(start_var)) start_var <- MCMCpack::riwish(pmwgs$n_pars * 3,diag(pmwgs$n_pars))
  start_a_half <- 1 / rgamma(n = pmwgs$n_pars, shape = 2, rate = 1)
  return(list(tmu = start_mu, tvar = start_var, tvinv = MASS::ginv(start_var), a_half = start_a_half))
}

get_group_level_standard <- function(parameters, s){
  # This function is modified for other versions
  mu <- parameters$tmu
  var <- parameters$tvar
  return(list(mu = mu, var = var))
}

fill_samples_standard <- function(samples, group_level, proposals, epsilon, j = 1, n_pars){
  samples$a_half[, j] <- group_level$a_half
  samples$last_theta_var_inv <- group_level$tvinv
  samples <- fill_samples_base(samples, group_level, proposals, epsilon, j = j, n_pars)
  return(samples)
}



gibbs_step_standard <- function(sampler, alpha){
  # Gibbs step for group means, with full covariance matrix estimation
  # tmu = theta_mu, tvar = theta_var
  last <- last_sample_standard(sampler$samples)
  hyper <- attributes(sampler)
  prior <- sampler$prior

  # Here mu is group mean, so we are getting mean and variance
  var_mu <- MASS::ginv(sampler$n_subjects * last$tvinv + prior$theta_mu_invar)
  mean_mu <- as.vector(var_mu %*% (last$tvinv %*% apply(alpha, 1, sum) +
                                     prior$theta_mu_invar %*% prior$theta_mu_mean))
  chol_var_mu <- t(chol(var_mu)) # t() because I want lower triangle.
  tmu <- mvtnorm::rmvnorm(1, mean_mu, chol_var_mu %*% t(chol_var_mu))[1, ]
  names(tmu) <- sampler$par_names

  # New values for group var
  theta_temp <- alpha - tmu
  cov_temp <- (theta_temp) %*% (t(theta_temp))
  if(!is.null(hyper$std_df)){
    B_half <- hyper$std_scale * diag(1, nrow = sampler$n_pars) + cov_temp # nolint
    tvar <- MCMCpack::riwish(hyper$std_df + sampler$n_subjects, B_half) # New sample for group variance
    tvinv <- MASS::ginv(tvar)
    # Sample new mixing weights.
    a_half <- NULL
  } else{
    B_half <- 2 * hyper$v_half * diag(1 / last$a_half) + cov_temp # nolint
    tvar <- MCMCpack::riwish(hyper$v_half + sampler$n_pars - 1 + sampler$n_subjects, B_half) # New sample for group variance
    tvinv <- MASS::ginv(tvar)

    # Sample new mixing weights.
    a_half <- 1 / rgamma(n = sampler$n_pars,shape = (hyper$v_half + sampler$n_pars) / 2,
                         rate = hyper$v_half * diag(tvinv) + 1/(hyper$A_half^2))
  }
  return(list(tmu = tmu,tvar = tvar,tvinv = tvinv,a_half = a_half,alpha = alpha))
}

get_conditionals_standard <- function(s, samples, n_pars){
  iteration <- samples$iteration
  pts2_unwound <- apply(samples$theta_var,3,unwind)
  all_samples <- rbind(samples$alpha[, s,],samples$theta_mu,pts2_unwound)
  mu_tilde <- rowMeans(all_samples)
  var_tilde <- cov(t(all_samples))
  condmvn <- condMVN(mean = mu_tilde, sigma = var_tilde,
                     dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                     X.given = c(samples$theta_mu[,iteration], unwind(samples$theta_var[,,iteration])))
  return(list(eff_mu = condmvn$condMean, eff_var = condmvn$condVar))
}

unwind <- function(var_matrix, ...) {
  y <- t(chol(var_matrix))
  diag(y) <- log(diag(y))
  y[lower.tri(y, diag = TRUE)]
}

last_sample_standard <- function(store) {
  list(
    tmu = store$theta_mu[, store$idx],
    tvar = store$theta_var[, , store$idx],
    tvinv = store$last_theta_var_inv,
    a_half = store$a_half[, store$idx]
  )
}

filtered_samples_standard <- function(sampler, filter){
  out <- list(
    theta_mu = sampler$samples$theta_mu[, filter],
    theta_var = sampler$samples$theta_var[, , filter],
    alpha = sampler$samples$alpha[, , filter],
    iteration = length(filter)
  )
}


dm_list <- function(dadm)
  # Makes data model into subjects list for use by likelihood
  # Assumes each subject has the same design.
{

  sub_design <- function(designs,isin)
    lapply(designs,function(x) {
      attr(x,"expand") <- attr(x,"expand")[isin]
      x
    })


  model <- attr(dadm,"model")
  p_names <- attr(dadm,"p_names")
  sampled_p_names <- attr(dadm,"sampled_p_names")
  designs <- attr(dadm,"designs")
  expand <- attr(dadm,"expand")
  s_expand <- attr(dadm,"s_expand")
  unique_nort <- attr(dadm,"unique_nort")
  expand_nort <- attr(dadm,"expand_nort")
  unique_nortR <- attr(dadm,"unique_nortR")
  expand_nortR <- attr(dadm,"expand_nortR")
  ok_trials <- attr(dadm,"ok_trials")
  ok_dadm_winner <- attr(dadm,"ok_dadm_winner")
  ok_dadm_looser <- attr(dadm,"ok_dadm_looser")
  ok_da_winner <- attr(dadm,"ok_da_winner")
  ok_da_looser <- attr(dadm,"ok_da_looser")
  expand_uc <- attr(dadm,"expand_uc")
  expand_lc <- attr(dadm,"expand_lc")
  adapt <- attr(dadm,"adapt")

  # winner on expanded dadm
  expand_winner <- attr(dadm,"expand_winner")
  # subjects for first level of lR in expanded dadm
  slR1=dadm$subjects[expand][dadm$lR[expand]==levels(dadm$lR)[[1]]]

  dl <- setNames(vector(mode="list",length=length(levels(dadm$subjects))),
                 levels(dadm$subjects))
  for (i in levels(dadm$subjects)) {
    isin <- dadm$subjects==i         # dadm
    isin1 <- s_expand==i             # da
    isin2 <- attr(dadm,"s_data")==i  # data
    dl[[i]] <- dadm[isin,]
    dl[[i]]$subjects <- factor(as.character(dl[[i]]$subjects))

    attr(dl[[i]],"model") <- model
    attr(dl[[i]],"p_names") <- p_names
    attr(dl[[i]],"sampled_p_names") <- sampled_p_names
    attr(dl[[i]],"designs") <- sub_design(designs,isin)
    attr(dl[[i]],"expand") <- expand[isin1]-min(expand[isin1]) + 1
    attr(dl[[i]],"s_expand") <- NULL

    attr(dl[[i]],"ok_dadm_winner") <- ok_dadm_winner[isin]
    attr(dl[[i]],"ok_dadm_looser") <- ok_dadm_looser[isin]

    attr(dl[[i]],"ok_da_winner") <- ok_da_winner[isin1]
    attr(dl[[i]],"ok_da_looser") <- ok_da_looser[isin1]

    attr(dl[[i]],"unique_nort") <- unique_nort[isin]
    attr(dl[[i]],"unique_nortR") <- unique_nortR[isin]

    isinlR1 <- slR1==i
    if (!is.null(expand_nort))
      attr(dl[[i]],"expand_nort") <-  expand_nort[isinlR1]- min( expand_nort[isinlR1]) + 1
    if (!is.null(expand_nortR))
      attr(dl[[i]],"expand_nortR") <- expand_nortR[isinlR1]-min(expand_nortR[isinlR1]) + 1

    attr(dl[[i]],"ok_trials") <- ok_trials[isin2]
    if (!is.null(expand_winner))
      attr(dl[[i]],"expand_winner") <- expand_winner[isin2]-min(expand_winner[isin2]) + 1

    if (!is.null(attr(dadm,"expand_uc")))
      attr(dl[[i]],"expand_uc") <- as.numeric(factor(expand_uc[isin2]))
    if (!is.null(attr(dadm,"expand_lc")))
      attr(dl[[i]],"expand_lc") <- as.numeric(factor(expand_lc[isin2]))

    if (!is.null(attr(dadm,"LT")))
      attr(dl[[i]],"LT") <- attr(dadm,"LT")[names(attr(dadm,"LT"))==i]
    if (!is.null(attr(dadm,"UT")))
      attr(dl[[i]],"UT") <- attr(dadm,"UT")[names(attr(dadm,"UT"))==i]
    if (!is.null(attr(dadm,"LC")))
      attr(dl[[i]],"LC") <- attr(dadm,"LC")[names(attr(dadm,"LC"))==i]
    if (!is.null(attr(dadm,"UC")))
      attr(dl[[i]],"UC") <- attr(dadm,"UC")[names(attr(dadm,"UC"))==i]

    # adapt models
    if (!is.null(adapt))
      attr(dl[[i]],"adapt") <- setNames(list(adapt[[i]],adapt$design),c(i,"design"))

  }
  dl
}


sample_store_base <- function(data, par_names, iters = 1, stage = "init") {
  subject_ids <- unique(data$subjects)
  n_pars <- length(par_names)
  n_subjects <- length(subject_ids)
  samples <- list(
    epsilon = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL)),
    origin = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL)),
    alpha = array(NA_real_,dim = c(n_pars, n_subjects, iters),dimnames = list(par_names, subject_ids, NULL)),
    stage = array(stage, iters),
    subj_ll = array(NA_real_,dim = c(n_subjects, iters),dimnames = list(subject_ids, NULL))
  )
}

fill_samples_base <- function(samples, group_level, proposals, epsilon, j = 1, n_pars){
  # Fill samples both group level and random effects
  samples$theta_mu[, j] <- group_level$tmu
  samples$theta_var[, , j] <- group_level$tvar
  if(!is.null(proposals)) samples <- fill_samples_RE(samples, proposals, epsilon,j, n_pars)
  return(samples)
}

fill_samples_RE <- function(samples, proposals, epsilon, j = 1, n_pars){
  # Only for random effects, separated because group level sometimes differs.
  samples$alpha[, , j] <- proposals[1:n_pars,]
  samples$subj_ll[, j] <- proposals[n_pars + 1,]
  samples$origin[,j] <- proposals[n_pars + 2,]
  samples$idx <- j
  samples$epsilon[,j] <- epsilon
  return(samples)
}

init <- function(pmwgs, start_mu = NULL, start_var = NULL,
                 verbose = FALSE, particles = 1000, n_cores = 1, epsilon = NULL) {
  # Gets starting points for the mcmc process
  # If no starting point for group mean just use zeros
  variant_funs <- attr(pmwgs, "variant_funs")
  startpoints <- variant_funs$get_startpoints(pmwgs, start_mu, start_var)
  if(n_cores > 1){
    proposals <- parallel::mclapply(X=1:pmwgs$n_subjects,FUN=start_proposals,
                          parameters = startpoints, n_particles = particles,
                          pmwgs = pmwgs, variant_funs = variant_funs,mc.cores = n_cores)
  } else{
    proposals <- lapply(X=1:pmwgs$n_subjects,FUN=start_proposals,
                        parameters = startpoints, n_particles = particles,
                        pmwgs = pmwgs, variant_funs = variant_funs)
  }
  proposals <- array(unlist(proposals), dim = c(pmwgs$n_pars + 2, pmwgs$n_subjects))

  # Sample the mixture variables' initial values.
  if(is.null(epsilon)) epsilon <- rep(set_epsilon(pmwgs$n_pars, verbose), pmwgs$n_subjects)
  if(length(epsilon) == 1) epsilon <- rep(epsilon, pmwgs$n_subjects)
  pmwgs$samples <- variant_funs$fill_samples(samples = pmwgs$samples, group_level = startpoints, proposals = proposals,
                                             epsilon = epsilon, j = 1, n_pars = pmwgs$n_pars)
  pmwgs$init <- TRUE
  return(pmwgs)
}


set_mix <- function(stage, mix, verbose) {
  if (stage == "sample") {
    mix <- c(0.1, 0.2, 0.7)
  } else {
    mix <- c(0.5, 0.5, 0.0)
  }
  if(verbose) message(sprintf("mix has been set to c(%s) based on the stage being run",  paste(mix, collapse = ", ")))
  return(mix)
}

start_proposals <- function(s, parameters, n_particles, pmwgs, variant_funs){
  #Draw the first start point
  group_pars <- variant_funs$get_group_level(parameters, s)
  proposals <- particle_draws(n_particles, group_pars$mu, group_pars$var)
  colnames(proposals) <- rownames(pmwgs$samples$theta_mu) # preserve par names
  lw <- apply(proposals,1,pmwgs$ll_func,dadm = pmwgs$data[[which(pmwgs$subjects == s)]])
  weight <- exp(lw - max(lw))
  idx <- sample(x = n_particles, size = 1, prob = weight)
  return(list(proposal = proposals[idx,], ll = lw[idx], origin = 2))
}


particle_draws <- function(n, mu, covar) {
  if (n <= 0) {
    return(NULL)
  }
  return(mvtnorm::rmvnorm(n, mu, covar))
}

set_epsilon <- function(n_pars, verbose = T) {
  if (n_pars > 15) {
    epsilon <- 0.1
  } else if (n_pars > 10) {
    epsilon <- 0.3
  } else {
    epsilon <- 0.5
  }
  if(verbose) message(sprintf("Epsilon has been set to %.1f based on number of parameters",epsilon))
  return(epsilon)
}


extend_sampler <- function(sampler, n_samples, stage) {
  # This function takes the sampler and extends it along the intended number of
  # iterations, to ensure that we're not constantly increasing our sampled object
  # by 1. Big shout out to the rapply function
  sampler$samples$stage <- c(sampler$samples$stage, rep(stage, n_samples))
  sampler$samples <- rapply(sampler$samples, f = function(x) extend_obj(x, n_samples), how = "replace")
  return(sampler)
}

extend_obj <- function(obj, n_extend){
  old_dim <- dim(obj)
  n_dimensions <- length(old_dim)
  if(is.null(old_dim) | n_dimensions == 1) return(obj)
  if(n_dimensions == 2){
    if(isSymmetric(round(obj, 5))) return(obj) #Don't extend priors and theta_mu_var_inv
  }
  new_dim <- c(rep(0, (n_dimensions -1)), n_extend)
  extended <- array(NA_real_, dim = old_dim +  new_dim, dimnames = dimnames(obj))
  extended[slice.index(extended,n_dimensions) <= old_dim[n_dimensions]] <- obj
  return(extended)
}



run_stage <- function(pmwgs,
                      stage,
                      iter = 1000,
                      particles = 100,
                      verbose = TRUE,
                      force_prev_epsilon = TRUE,
                      n_cores = 1,
                      epsilon = NULL,
                      pstar = NULL,
                      mix = NULL,
                      pdist_update_n = ifelse(stage == "sample", 50, NA),
                      epsilon_upper_bound = 15,
                      n_cores_conditional = 1,
                      eff_mu = NULL, eff_var = NULL,
                      thin = NULL,
                      thin_eff_only = FALSE) {
  # Set defaults for NULL values
  mix <- set_mix(stage, mix, verbose)
  # Set necessary local variables
  # Set stable (fixed) new_sample argument for this run
  n_pars <- pmwgs$n_pars
  components <- attr(pmwgs$data, "components")
  # Display stage to screen
  if(verbose){
    msgs <- list(
      burn = "Phase 1: Burn in\n",
      adapt = "Phase 2: Adaptation\n",
      sample = "Phase 3: Sampling\n"
    )
    cat(msgs[[stage]])
  }

  alphaStar=-qnorm(pstar/2) #Idk about this one, Fan 2014
  n0=round(5/(pstar*(1-pstar))) #Also not questioning this math for now
  if(is.null(epsilon) | force_prev_epsilon){
    epsilon <- pmwgs$samples$epsilon[,ncol(pmwgs$samples$epsilon)]
  }
  if(length(epsilon) == 1){
    epsilon <- rep(epsilon, pmwgs$n_subjects)
  }
  if(length(particles == 1)){
    particles <- rep(particles, pmwgs$n_subjects)
  }
  epsilon <- replicate(length(components), epsilon)
  if(stage!= "burn") components <- pmwgs$n_pars - length(pmwgs$subject_covariates)
  # Build new sample storage
  pmwgs <- extend_sampler(pmwgs, iter, stage)
  # create progress bar
  if (verbose) {
    pb <- accept_progress_bar(min = 0, max = iter)
  }
  start_iter <- pmwgs$samples$idx

  data <- pmwgs$data
  subjects <- pmwgs$subjects
  variant_funs <- attr(pmwgs, "variant_funs")
  # Main iteration loop
  for (i in 1:iter) {
    if (verbose) {
      accRate <- mean(accept_rate(pmwgs))
      update_progress_bar(pb, i, extra = accRate)
    }

    j <- start_iter + i
    # Gibbs step
    pars <- variant_funs$gibbs_step(pmwgs, rbind(pmwgs$samples$alpha[,,j-1], pmwgs$subject_covariates))
    # Particle step
    proposals_raw <- parallel::mclapply(X=1:pmwgs$n_subjects,FUN = new_particle, data, particles, pars, eff_mu,
                       eff_var, mix, pmwgs$ll_func, epsilon, subjects, components, pmwgs$subject_covariates, variant_funs, mc.cores =n_cores)
    proposals <- array(unlist(proposals_raw), dim = c(pmwgs$n_pars + 2, pmwgs$n_subjects))

    #Fill samples
    pmwgs$samples <- variant_funs$fill_samples(samples = pmwgs$samples, group_level = pars,
                                               proposals = proposals, epsilon = rowMeans(epsilon), j = j, n_pars = pmwgs$n_pars)

    # Update epsilon
    if(!is.null(pstar)){
      if(j > n0){
        prev_comp <- 0
        comp_n <- 0
        for(component in components){
          comp_n <- comp_n + 1
          acc <-  pmwgs$samples$alpha[component,,j] != pmwgs$samples$alpha[component,,(j-1)]
          epsilon[,comp_n] <-pmin(update.epsilon(epsilon[,comp_n]^2, acc, pstar, j, component-prev_comp, alphaStar), epsilon_upper_bound)

          prev_comp <- component
        }
      }
    }

  }
  if (verbose) close(pb)
  if(!is.null(thin) & !thin_eff_only){
    pmwgs <- thin_objects(pmwgs, thin, i)
  }
  return(pmwgs)
}

new_particle <- function (s, data, num_particles, parameters, eff_mu = NULL,
                          eff_var = NULL, mix_proportion = c(0.5, 0.5, 0),
                          likelihood_func = NULL, epsilon = NULL, subjects,
                          components, subject_covariates = 0, variant_funs)
{
  num_particles <- num_particles[s]
  start_par <- length(subject_covariates) + 1
  group_pars <- variant_funs$get_group_level(parameters, s)
  n_components <- length(components)
  proposal_out <- numeric(length(group_pars$mu))
  ll <- 0
  for(i in 1:n_components){
    idx <- start_par:components[i]
    eff_mu <- eff_mu[idx, s]
    eff_var <- eff_var[idx, idx , s]
    mu <- group_pars$mu[idx]
    var <- group_pars$var[idx, idx]
    subj_mu <- parameters$alpha[idx, s]
    particle_numbers <- numbers_from_proportion(mix_proportion, num_particles)
    cumuNumbers <- cumsum(particle_numbers)
    pop_particles <- particle_draws(particle_numbers[1], mu,
                                    var)
    ind_particles <- particle_draws(particle_numbers[2], subj_mu,
                                    var * epsilon[s,i]^2)
    if(mix_proportion[3] == 0){
      eff_particles <- NULL
    } else{
      eff_particles <- particle_draws(particle_numbers[3], eff_mu, eff_var)
    }
    proposals <- rbind(pop_particles, ind_particles, eff_particles)
    colnames(proposals) <- names(subj_mu)
    proposals[1, ] <- subj_mu
    if(n_components > 1){
      lw <- apply(proposals, 1, likelihood_func, dadm = data[[which(subjects == s)]], component = i)
    } else{
      lw <- apply(proposals, 1, likelihood_func, dadm = data[[which(subjects == s)]])
    }
    lp <- mvtnorm::dmvnorm(x = proposals, mean = mu, sigma = var,
                           log = TRUE)
    prop_density <- mvtnorm::dmvnorm(x = proposals, mean = subj_mu,
                                     sigma = var * (epsilon[s,i]^2))
    if (mix_proportion[3] == 0) {
      eff_density <- 0
    }
    else {
      eff_density <- mvtnorm::dmvnorm(x = proposals, mean = eff_mu, sigma = eff_var)
    }
    lm <- log(mix_proportion[1] * exp(lp) + (mix_proportion[2] *
                                               prop_density) + (mix_proportion[3] * eff_density))
    l <- lw + lp - lm
    weights <- exp(l - max(l))
    idx_ll <- sample(x = num_particles, size = 1, prob = weights)
    origin <- min(which(idx_ll <= cumuNumbers))
    ll <- ll + lw[idx_ll]
    proposal_out[idx] <- proposals[idx_ll,]
    start_par <- components[i] + 1
  }
  proposal_out <- c(proposal_out, subject_covariates)
  return(list(proposal = proposal_out, ll = ll, origin = origin))
}

update.epsilon<- function(epsilon2, acc, p, i, d, alpha) {
  c <- ((1-1/d)*sqrt(2*pi)*exp(alpha^2/2)/(2*alpha) + 1/(d*p*(1-p)))
  Theta <- log(sqrt(epsilon2))
  Theta <- Theta+c*(acc-p)/max(200, i/d)
  return(exp(Theta))
}


numbers_from_proportion <- function(mix_proportion, num_particles = 1000) {
  numbers <- stats::rbinom(n = 2, size = num_particles, prob = mix_proportion)
  if (mix_proportion[3] == 0) {
    numbers[3] <- 0
    numbers[2] <- num_particles - numbers[1]
  } else {
    numbers[3] <- num_particles - sum(numbers)
  }
  return(numbers)
}
