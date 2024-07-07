## standard IS2 written by Reilly Innes and David Gunawan, from Tran et al. 2021
## set up environment and packages


IS2 <- function(samples, stage = "sample", filter = 0, IS_samples = 1000, stepsize_particles = 500, max_particles = 5000, n_cores = 1, df = 5){
  ###### set up variables #####
  info <- add_info_base(samples)
  idx <- which(samples$samples$stage == stage)
  if (length(filter)==1) {
    if(filter > 0) idx <- idx[-c(1:filter)]
  } else idx <- idx[filter]
  all_pars <- info$variant_funs$get_all_pars_IS2(samples, idx, info)
  if(!is.null(all_pars$X)){
    muX<-apply(all_pars$X,2,mean)
    varX<-stats::cov(all_pars$X)
  } else{
    muX<-info$prior$theta_mu_mean
    varX<-info$prior$theta_mu_var
  }
  prop_theta <- mvtnorm::rmvt(IS_samples,sigma = varX, df=df, delta=muX)



  logw_num <- auto_mclapply(X=1:IS_samples,
                            FUN = compute_lw_num,
                            prop_theta = prop_theta,
                            stepsize_particles = stepsize_particles,
                            max_particles = max_particles,
                            mu_tilde=all_pars$mu_tilde,
                            var_tilde = all_pars$var_tilde,
                            info = all_pars$info,
                            mc.cores = n_cores)
  # sub_and_group <- simplify2array(lapply(logw_num, FUN = function(x) return(x$sub_and_group)))
  # prior_and_jac <- sapply(logw_num, FUN = function(x) return(x$prior_and_jac))
  #
  logw_den <- mvtnorm::dmvt(prop_theta, delta=muX, sigma=varX,df=df, log = TRUE)
  finished <- unlist(logw_num) - logw_den
  max.lw <- max(finished)
  centred.lw <- exp(finished-max.lw) #takes off the max and gets mean (avoids infs)
  lw <- log(centred.lw)+max.lw #puts max back on to get the lw
  # return(list(sub_and_group = sub_and_group, prior_and_jac = sub_and_group, logw_den = logw_den))
  # logw_den <- mvtnorm::dmvt(prop_theta, delta=muX, sigma=varX,df=df, log = TRUE)
  #
  # finished <- unlist(logw_num) - logw_den
  return(finished)

}

get_sub_weights <- function(stepsize_particles, condMean, condVar, prop_theta, info, sub){
  wmix <- .95
  n1 <- stats::rbinom(n=1,size=stepsize_particles,prob=wmix)
  if (n1<2) n1 <- 2
  if (n1>(stepsize_particles-2)) n1 <- stepsize_particles-2 ## These just avoid degenerate arrays.
  n2 <- stepsize_particles-n1
  data <- info$data
  particles1 <- mvtnorm::rmvnorm(n1, condMean,condVar)
  # Group level
  info$sub <- sub # Needed for variant lm
  particles2 <- info$variant_funs$group_dist_IS2(n_samples=n2, parameters = prop_theta,
                                                 sample=TRUE, info = info)
  particles <- rbind(particles1,particles2)
  # names for ll function to work
  colnames(particles) <- info$par_names
  # do lba log likelihood with given parameters for each subject, gets density of particle from ll func
  lw_first <- calc_ll_manager(particles, dadm = data[[sub]], info$ll_func)
  # below gets second part of equation 5 numerator ie density under prop_theta
  lw_second <- apply(particles, 1, info$variant_funs$group_dist_IS2, prop_theta, FALSE, NULL, info)
  # below is the denominator - ie mix of density under conditional and density under pro_theta
  lw_third <- log(wmix*pmax(1e-25 * info$n_randeffect, mvtnorm::dmvnorm(particles, condMean, condVar)) + (1-wmix) * exp(lw_second))
  # does equation 5
  lw <- lw_first+lw_second-lw_third
  return(lw)
}

get_logp=function(prop_theta,stepsize_particles, max_particles, mu_tilde,var_tilde, info){
  # Unload for legibility
  n_subjects <- info$n_subjects
  var_opt_sub <- 1/n_subjects

  # for each subject, get N IS samples (particles) and find log weight of each
  lw_subs <- numeric(n_subjects)
  for (j in 1:n_subjects){
    var_z_sub <- 999 # just a place holder > var_opt_sub
    n_total <- 0
    lw <- numeric()
    # generate the particles from the conditional MVnorm AND mix of group level proposals
    conditional = condMVN(mean=mu_tilde[j,],sigma=var_tilde[j,,],
                          dependent.ind=c(1:info$n_randeffect),
                          given.ind=info$given.ind,
                          X.given=prop_theta[info$X.given_ind],
                          check.sigma = F)
    while((var_z_sub > var_opt_sub) & n_total < max_particles){
      n_total <- n_total + stepsize_particles
      lw_tmp <- get_sub_weights(stepsize_particles = stepsize_particles, condMean = conditional$condMean,
                                condVar = conditional$condVar, prop_theta = prop_theta,
                                info = info, sub = j)

      # lw <- -weights(psis(-lw, log = F)) # default args are log=TRUE, normalize=TRUE
      lw <- c(lw, lw_tmp)
      max_lw <- max(lw)
      var_z_sub = sum(exp(2*(lw-max_lw)))/(sum(exp(lw-max_lw)))^2-1/n_total
    }
    weight <- exp(lw-max_lw)
    # lw_subs[j, 1:length(weight)] <- max_lw+log((weight))
    lw_subs[j] <- max_lw+log(mean(weight))
  }
  # sum the logp and return
  return(sum(lw_subs))
}

compute_lw_num=function(i, prop_theta,stepsize_particles, max_particles, mu_tilde,var_tilde,info){
  # sub_and_group <- get_logp(prop_theta[i,], stepsize_particles, max_particles, mu_tilde, var_tilde, info)
  # prior_and_jac <- info$variant_funs$prior_dist_IS2(parameters = prop_theta[i,], info)
  # return(list(sub_and_group = sub_and_group, prior_and_jac = prior_and_jac))
  logp.out <- get_logp(prop_theta[i,], stepsize_particles, max_particles, mu_tilde, var_tilde, info)
  logw_num <- logp.out+info$variant_funs$prior_dist_IS2(parameters = prop_theta[i,], info)
  return(logw_num)
}

add_info_base <- function(samples){
  info <- list(
    n_randeffect = samples$n_pars,
    n_subjects = samples$n_subjects,
    variant_funs = attr(samples, "variant_funs"),
    par_names = samples$par_names,
    data = samples$data,
    ll_func = samples$ll_func,
    prior = samples$prior,
    hyper = attributes(samples)
  )
  return(info)
}
