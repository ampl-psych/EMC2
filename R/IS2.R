## standard IS2 written by Reilly Innes and David Gunawan, from Tran et al. 2021
## set up environment and packages


IS2 <- function(samples, filter = "sample", subfilter = 0, IS_samples = 1000, stepsize_particles = 500, max_particles = 5000, n_cores = 1, df = 5){
  ###### set up variables #####
  info <- add_info_base(samples)
  idx <- which(samples$samples$stage == filter)
  if (length(subfilter)==1) {
    if(subfilter > 0) idx <- idx[-c(1:subfilter)]
  } else idx <- idx[subfilter]
  all_pars <- variant_funs$get_all_pars(samples, idx, info)
  muX<-apply(all_pars$X,2,mean)
  varX<-cov(all_pars$X)

  prop_theta=mvtnorm::rmvt(IS_samples,sigma = varX, df=df, delta=muX)
  # out <- numeric(nrow(prop_theta))
  # for(i in 1:nrow(prop_theta)){
  #   out[i] <- variant_funs$prior_dist(parameters = prop_theta[i,], all_pars$info)
  # }
  #do the sampling
  logw_num <- mclapply(X=1:IS_samples,
                       FUN = compute_lw_num,
                       prop_theta = prop_theta,
                       stepsize_particles = stepsize_particles,
                       max_particles = max_particles,
                       mu_tilde=all_pars$mu_tilde,
                       var_tilde = all_pars$var_tilde,
                       info = all_pars$info,
                       mc.cores = n_cores)
  logw_den <- mvtnorm::dmvt(prop_theta, delta=muX, sigma=varX,df=df, log = TRUE)

  finished <- unlist(logw_num) - logw_den
  max.lw <- max(finished)
  mean.centred.lw <- mean(exp(finished-max.lw)) #takes off the max and gets mean (avoids infs)
  lw <- log(mean.centred.lw)+max.lw #puts max back on to get the lw
  return(list(lw = lw, finished = finished))
}

get_sub_weights <- function(stepsize_particles, condMean, condVar, prop_theta, info, sub){
  wmix = c(.05, .95)
  n1=rbinom(n=1,size=stepsize_particles,prob=wmix)
  if (n1<2) n1=2
  if (n1>(stepsize_particles-2)) n1=stepsize_particles-2 ## These just avoid degenerate arrays.
  n2=stepsize_particles-n1
  data <- info$data
  particles1 <- mvtnorm::rmvnorm(n1, condMean,condVar)
  # Group level
  particles2 <- variant_funs$group_dist(n_samples=n2, parameters = prop_theta,
                                        sample=TRUE, info = info)
  particles <- rbind(particles1,particles2)
  # names for ll function to work
  colnames(particles) <- info$par_names
  # do lba log likelihood with given parameters for each subject, gets density of particle from ll func
  lw_first <- apply(particles, 1, info$ll_func, dadm = data[[sub]])
  # below gets second part of equation 5 numerator ie density under prop_theta
  lw_second <- apply(particles, 1, variant_funs$group_dist, prop_theta, FALSE, NULL, info)
  # below is the denominator - ie mix of density under conditional and density under pro_theta
  lw_third <- log(wmix*pmax(1e-25 * info$n_randeffect, dmvnorm(particles, condMean, condVar)) + (1-wmix) * exp(lw_second))
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
    lw_subs[j] <- max_lw+log(mean(weight))
  }
  # sum the logp and return
  return(sum(lw_subs))
}

compute_lw_num=function(i, prop_theta,stepsize_particles, max_particles, mu_tilde,var_tilde,info){
  logp.out <- get_logp(prop_theta[i,], stepsize_particles, max_particles, mu_tilde, var_tilde, info)
  logw_num <- logp.out+variant_funs$prior_dist(parameters = prop_theta[i,], info)
  return(logw_num)
}

add_info_base <- function(samples){
  info <- list(
    n_randeffect = samples$n_pars,
    n_subjects = samples$n_subjects,
    par_names = samples$par_names,
    data = samples$data,
    ll_func = samples$ll_func,
    prior = samples$prior,
    hyper = attributes(samples)
  )
  return(info)
}

condMVN <- function (mean, sigma, dependent.ind, given.ind, X.given, check.sigma = TRUE)
{
  if (missing(dependent.ind))
    return("You must specify the indices of dependent random variables in `dependent.ind'")
  if (missing(given.ind) & missing(X.given))
    return(list(condMean = mean[dependent.ind], condVar = as.matrix(sigma[dependent.ind,
                                                                          dependent.ind])))
  if (length(given.ind) == 0)
    return(list(condMean = mean[dependent.ind], condVar = as.matrix(sigma[dependent.ind,
                                                                          dependent.ind])))
  if (length(X.given) != length(given.ind))
    stop("lengths of `X.given' and `given.ind' must be same")
  if (check.sigma) {
    if (!isSymmetric(sigma))
      stop("sigma is not a symmetric matrix")
    eigenvalues <- eigen(sigma, only.values = TRUE)$values
    if (any(eigenvalues < 1e-08))
      stop("sigma is not positive-definite")
  }
  B <- sigma[dependent.ind, dependent.ind]
  C <- sigma[dependent.ind, given.ind, drop = FALSE]
  D <- sigma[given.ind, given.ind]
  CDinv <- C %*% chol2inv(chol(D))
  cMu <- c(mean[dependent.ind] + CDinv %*% (X.given - mean[given.ind]))
  cVar <- B - CDinv %*% t(C)
  list(condMean = cMu, condVar = cVar)
}

