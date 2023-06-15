prior_samples_alpha <- function(theta_mu,theta_var,n=1e3)
  # samples from prior for alpha implied by hyper model
  # (mixture over hyper posterior samples)
{
  if(is.null(dim(theta_mu))){
    out <- mvtnorm::rmvnorm(n,theta_mu,theta_var)
  } else{
    out <- t(sapply(sample.int(dim(theta_mu)[2],n,TRUE),function(x){
      mvtnorm::rmvnorm(1,theta_mu[x,],theta_var[,,x])
    }))
  }
  return(out)
}


# prior_samples <- function(samps,type=c("mu","variance","covariance","correlation","sigma")[1],n=1e4)
#   # Samples from prior for standard at hyper levels for three types,
#   # or for single at subject level (in which case type="mu")
# {
#   if (type=="mu") {
#     out <- mvtnorm::rmvnorm(n, mean = samps$prior$theta_mu_mean,
#                             sigma = samps$prior$theta_mu_var)
#     dimnames(out)[[2]] <- samps$par_names
#     return(out)
#   } else {
#     var <- array(NA_real_, dim = c(samps$n_pars, samps$n_pars, n),
#                  dimnames = list(samps$par_names, samps$par_names, NULL))
#     # Can't specify n in riwish?
#     hyper <- attributes(samps)
#     for(i in 1:n){
#       a_half <- 1 / rgamma(n = samps$n_pars,shape = 1/2,
#                            rate = 1/(hyper$A_half^2))
#       var[,,i] <- riwish(hyper$v_half + samps$n_pars - 1, 2 * hyper$v_half * diag(1 / a_half))
#     }
#     if (type=="sigma") return(var)
#     if (type=="variance") return(t(apply(var,3,diag)))
#     if (type=="correlation")
#       var <- array(apply(var,3,cov2cor),dim=dim(var),dimnames=dimnames(var))
#     lt <- lower.tri(var[,,1])
#     return(t(apply(var,3,function(x){x[lt]})))
#   }
# }

#' Convenience function to plot from the prior
#'
#' @param prior A list of prior samples
#' @param type Optional. Otherwised inferred from the prior samples
#'
#' @return NULL. Makes a plot of the prior samples
#' @export
#'
#' @examples
plot_prior <- function(prior, type = NULL,add_density=FALSE,adjust=1,breaks=50,
                       layout=c(3,3)){
  if(is.null(type)) type <- names(prior)
  for(typ in type){
    samples <- prior[["alpha"]]
    par(mfrow = layout)
    par_names <- colnames(samples)
    for(i in 1:ncol(samples)){
      if(!any(samples[,i] < 0) || !any(samples[,i] > 0)){
        quants <- quantile(abs(samples[,i]), probs = 0.95)
      } else{
        quants <- quantile(abs(samples[,i]), probs = 0.995)
      }
      filtered <- samples[,i][abs(samples[,i]) < quants]
      hist(filtered, breaks = breaks, main = par_names[i], prob = TRUE,
           xlab = type, cex.lab = 1.25, cex.main = 1.5)
      if (add_density) lines(density(filtered,adjust=adjust), col = "red")
    }
  }
}

get_prior_samples <- function(samples,selection,filter,thin,subfilter,n_prior)
  # get matrix of prior samples for different parameter types
{
  if (inherits(samples, "pmwgs")) samps <- samples else samps <- samples[[1]]
  if (selection=="alpha") {
    if(!is.null(samples[[1]]$samples$theta_mu)){
      if (inherits(samples, "pmwgs")) {
        theta_mu <- as_Mcmc(samples,selection="mu",filter=filter,
                            thin=thin,subfilter=subfilter)
        theta_var <- get_sigma(samples,filter=filter,
                               thin=thin,subfilter=subfilter)
      } else {
        theta_mu <- do.call(rbind,as_mcmc.list(samples,selection="mu",
                                               filter=filter,thin=thin,subfilter=subfilter))
        theta_var <-   abind(lapply(samples,get_sigma,filter=filter,thin=thin,subfilter=subfilter))
      }
      psamples <- prior_samples_alpha(theta_mu,theta_var,n_prior)
      colnames(psamples) <- colnames(theta_mu)
    } else{
      theta_mu <- samples[[1]]$prior$theta_mu_mean
      theta_var <- samples[[1]]$prior$theta_mu_var
      psamples <- prior_samples_alpha(theta_mu,theta_var,n_prior)
      colnames(psamples) <- samples[[1]]$par_names
    }
    return(psamples)
  } else {
    variant_funs <- attr(samples[[1]], "variant_funs")
    design <- list(model = attr(samples[[1]]$data[[1]], "model"))
    attr(design, "p_vector") <- samplers[[1]]$samples$alpha[,1,1] # just need the names in the right format
    psamples <- variant_funs$get_prior(design = design, N = n_prior, prior = samples[[1]]$prior, type = selection, map = FALSE)[[1]]
    return(psamples)
  }
}
