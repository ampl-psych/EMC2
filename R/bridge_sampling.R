make_info <- function(samples, variant_funs){
  info <- list(
    subjects = samples$subjects,
    par_names = samples$par_names,
    n_pars = samples$n_pars,
    n_subjects = samples$n_subjects,
    prior = samples$prior,
    variant_funs = variant_funs
  )
  info <- variant_funs$bridge_add_info(info, samples)

}



eval.unnormalized.posterior <- function(samples_iter, gen_samples, data, m, L, info, cores_for_props = 1, cores_per_prop = 1) {

  ### evaluate unnormalized posterior for posterior and generated samples
  n_post <- nrow(samples_iter)
  e <- Brobdingnag::as.brob( exp(1) )
  twom_min_s <- matrix(2*m, nrow = n_post, ncol = length(m), byrow = TRUE) - samples_iter
  m_min_gen <- matrix(m, nrow = n_post, ncol = length(m), byrow = TRUE) - gen_samples %*% t(L)
  m_plus_gen <- matrix(m, nrow = n_post, ncol = length(m), byrow = TRUE) + gen_samples %*% t(L)
  samples_list <- list(
    q11.a = samples_iter,
    q11.b = twom_min_s,
    q21.a = m_min_gen,
    q21.b = m_plus_gen
  )
  mls <- mclapply(samples_list, h.unnormalized.posterior, data = data, info = info, n_cores = cores_per_prop,
                  mc.cores = cores_for_props)
  q11 <- log(e^(mls$q11.a) + e^(mls$q11.b))
  q21 <- log(e^(mls$q21.a) + e^(mls$q21.b))
  return(list(q11 = q11, q21 = q21))

}

# Q11 = samples + 2*m - samples
# Q21 = samples + cov


h.unnormalized.posterior <- function(proposals, data, info, n_cores) {
  n_pars <- info$n_pars
  n_subjects <- info$n_subjects
  proposals_list <- vector("list", length = n_subjects)
  for(i in 1:n_subjects){
    tmp <- proposals[,((i-1)*n_pars + 1):(i*n_pars)]
    colnames(tmp) <- info$par_names
    proposals_list[[i]] <- tmp
  }
  lws <- parallel::mcmapply(calc_ll_for_group, proposals_list, data, MoreArgs = list(ll = info$ll_func), mc.cores = n_cores)
  lw <- rowSums(lws)
  proposals_group <- proposals[,info$group_idx]
  gr_pr_jac <- info$variant_funs$bridge_group_and_prior_and_jac(proposals_group, proposals_list, info)
  return(lw + gr_pr_jac)
}

run.iterative.scheme <- function(q11, q12, q21, q22, r0, tol,
                                 criterion, L, silent, maxiter,
                                 neff) {

  ### run iterative updating scheme (using "optimal" bridge function,
  ### see Meng & Wong, 1996)

  l1 <- -log(2) + determinant(L)$modulus + q11 - q12 # log(l)
  l2 <- -log(2) + determinant(L)$modulus + q21 - q22 # log(ltilde)

  lstar <- median(l1)
  n.1 <- length(l1)
  n.2 <- length(l2)

  if (is.null(neff)) {
    neff <- n.1
  }

  s1 <- neff/(neff + n.2)
  s2 <- n.2/(neff + n.2)
  r <- r0
  logml <- log(r) + lstar
  i <- 1

  r_vals <- r
  logml_vals <-  logml
  e <- as.brob( exp(1) )

  criterion_val <- 1 + tol

  while (criterion_val > tol && i <= maxiter) {

    if (! silent) {
      cat(paste0("Iteration: ", i, "\n"))
    }

    rold <- r
    logmlold <- logml
    numi <- as.numeric( e^(l2 - lstar)/(s1 * e^(l2 - lstar) + s2 *  r) )
    deni <- as.numeric( 1/(s1 * e^(l1 - lstar) + s2 * r) )
    r <- (n.1/n.2) * sum(numi)/sum(deni)
    logml <- log(r) + lstar
    r_vals <- c(r_vals, r)
    logml_vals <- c(logml_vals, logml)
    i <- i + 1

    criterion_val <- switch(criterion, "r" = abs((r - rold)/r),
                            "logml" = abs((logml - logmlold)/logml))

  }

  if (i > maxiter) {
    logml <- NA
  }

  return(list(logml = logml, niter = i - 1, r_vals = r_vals,
              logml_vals = logml_vals, tol = tol, neff = neff,
              criterion = criterion, maxiter = maxiter))

}


bridge_sampling <- function(samples, n_eff, split_idx, cores_for_props = 1, cores_per_prop = 1, maxiter = 500,
                            r0 = 1e-5, tol1 = 1e-10, tol2 = 1e-6){
  variant_funs <- attr(samples, "variant_funs")
  data <- samples$data
  info <- make_info(samples, variant_funs)
  n_pars <- samples$n_pars
  n_subjects <- samples$n_subjects
  idx <- samples$samples$stage == "sample"
  all_samples <- matrix(NA_real_, nrow = sum(idx), ncol = n_pars * n_subjects)
  for(i in 1:n_subjects){
    all_samples[,((i-1)*n_pars + 1):(i*n_pars)] <- t(samples$samples$alpha[,i,idx])
  }
  all_samples <- variant_funs$bridge_add_group(all_samples, samples, idx)
  samples_fit <- all_samples[split_idx,]
  samples_iter <- all_samples[-split_idx,]
  if(nrow(samples_fit) != nrow(samples_iter)){
    samples_fit <- samples_fit[1:nrow(samples_iter),]
  }

  m <- colMeans(samples_fit)
  V <- as.matrix(Matrix::nearPD(cov(samples_fit))$mat)
  L <- t(chol(V))
  gen_samples <- mvtnorm::rmvnorm(nrow(samples_fit), sigma = diag(ncol(samples_fit)))

  q12 <- mvtnorm::dmvnorm((samples_iter - matrix(m, nrow = nrow(samples_iter),
                                                 ncol = length(m),
                                                 byrow = TRUE)) %*%
                            t(solve(L)), sigma = diag(ncol(samples_iter)), log = TRUE)
  q22 <- mvtnorm::dmvnorm(gen_samples, sigma = diag(ncol(samples_fit)), log = TRUE)
  info <- make_info(samples, variant_funs)
  qList <- eval.unnormalized.posterior(samples_iter = samples_iter, gen_samples = gen_samples,
                                       data = data, m = m, L =L, info = info, cores_for_props = cores_for_props, cores_per_prop = cores_per_prop)

  q11 <- qList$q11
  q21 <- qList$q21
  tmp <- run.iterative.scheme(q11 = q11, q12 = q12, q21 = q21,
                              q22 = q22, r0 = r0, tol = tol1,
                              L = L, silent = T,
                              maxiter = maxiter, neff = n_eff,
                              criterion = "r")

  if (is.na(tmp$logml)) {
    lr <- length(tmp$r_vals)
    # use geometric mean as starting value
    r0_2 <- sqrt(tmp$r_vals[lr - 1]*tmp$r_vals[lr])
    tmp <- run.iterative.scheme(q11 = q11, q12 = q12, q21 = q21,
                                q22 = q22, r0 = r0_2, tol = tol2,
                                L = L, silent = T,
                                maxiter = maxiter, neff = n_eff,
                                criterion = "logml")
  }

  tmp$logml

}


#' Estimating Marginal likelihoods using WARP-III bridge sampling
#'
#' @param samplers A list with a set of converged samplers
#' @param filter A character indicating which stage to use, default is sample
#' @param subfilter An integer or vector indicating whether to exclude up to (integer case), or which to exclude
#' @param repetitions An integer. How many times to repeat the bridge sampling scheme. Can help get an estimate of stability of the estimate.
#' @param cores_for_props Integer. Warp-III evaluates the posterior over 4 different proposal densities. If you have the power, 4 cores will do this in parallel, 2 also helps, 3 not really.
#' @param cores_per_prop Integer. Per density we can also parallelize across subjects. Eventual cores will be cores_for_props * cores_per_prop. Prioritize cores_for_props being 4.
#' @param both_splits Boolean. Bridge sampling uses a proposal density and a target density. We can estimate the stability of our samples and therefore MLL estimate, by running 2 bridge sampling iteratations.
#' The first one uses the first half of the samples as the proposal and the second half as the target, the second run uses the opposite. If this is is FALSE, it will only run bridge samplign once and
#'  it will instead do a odd-even iterations split to get a more reasonable estimate for just one run.
#' @param maxiter Integer, how many iterations to run in a post-processing scheme, best left at default
#' @param r0 Numeric. Hyperparameter in a post-processing scheme, best left at default
#' @param tol1 Numeric. Hyperparameter in a post-processing scheme, best left at default
#' @param tol2 Numeric. Hyperparameter in a post-processing scheme, best left at default
#'
#' @return A vector of length repetitions which contains the Marginal Log Likelihood estimates per repetition
#' @export
#'
run_bridge_sampling <- function(samplers, filter = "sample", subfilter = 0, repetitions = 1, cores_for_props = 4,  cores_per_prop = 1, both_splits = T, maxiter = 500,
                                r0 = 1e-5, tol1 = 1e-10, tol2 = 1e-6){

  if(subfilter != 0){
    samplers <- lapply(samplers, remove_iterations, select = subfilter, filter = filter)
  }
  n_eff <- round(median(es_pmwg(as_mcmc.list(samplers, selection = "alpha",
                                             filter = "sample"), print_summary = F))/2)
  samples <- merge_samples(samplers)
  idx <- samples$samples$stage == filter
  mls <- numeric(repetitions)
  for(i in 1:repetitions){
    if(both_splits){
      split1 <- seq(1, round(sum(idx)/2))
      s1 <- bridge_sampling(samples, n_eff, split1, cores_for_props = cores_for_props, cores_per_prop = cores_per_prop,
                            maxiter = maxiter,
                            r0 = r0, tol1 = tol1, tol2 = tol2)
      split2 <- seq(round(sum(idx)/2 + 1) : sum(idx))
      s2 <- bridge_sampling(samples, n_eff, split2, cores_for_props = cores_for_props, cores_per_prop = cores_per_prop,
                            maxiter = maxiter,
                            r0 = r0, tol1 = tol1, tol2 = tol2)
      if(abs(s1 - s2) > 1) warning("First split and second split marginal likelihood estimates are off by 1 log point. \n This usually means that your MCMC chains aren't completely stable yet. \n Consider running the MCMC chain longer if you need more precise estimates (e.g. when comparing different priors)")
      mls[i] <- mean(c(s1, s2))
    } else{
      split_idx <- seq(1, sum(idx), by = 2)
      mls[i] <- bridge_sampling(samples, n_eff, split_idx, cores_for_props = cores_for_props, cores_per_prop = cores_per_prop,
                            maxiter = maxiter,
                            r0 = r0, tol1 = tol1, tol2 = tol2)
    }
  }
  return(mls)
}
