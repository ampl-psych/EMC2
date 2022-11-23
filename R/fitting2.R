run_samplers <- function(samplers, stage, iter = NULL, max_gd = NULL, min_es = NULL, min_unique,
                         p_accept = .8, step_size = 100, verbose = F,
                         particles = NULL, particle_factor = 50, cores_per_chain = 1,
                         cores_for_chains = length(samplers), epsilon = NULL, max_trys = NULL){
  attributes <- get_attributes(samplers)
  progress <- check_progress(samplers, stage, iter, max_gd, min_es, min_unique, max_trys, step_size)
  samplers <- progress$samplers
  while(!progress$done){
    samplers <- add_proposals(samplers, stage, cores_per_chain)
    samplers <- parallel::mclapply(samplers,run_stage, stage = stage, iter= progress$step_size,
                                   verbose=verbose,  particles=particles,particle_factor=particle_factor,
                                   p_accept=p_accept, epsilon = epsilon,  mix=mix,
                                   n_cores=cores_per_chain, mc.cores = cores_for_chains)
    progress <- check_progress(stage, max_gd, min_es, min_unique, stage, max_trys, progress)
    samplers <- progress$samplers
  }
  samplers <- get_attributes(samplers, give_back = TRUE)
  return(samplers)
}

add_proposals <- function(samplers, stage){
  if(stage != "preburn"){
    samplers <- create_cov_proposals(samplers_new, verbose)
  }
  if(stage == "sample"){

  }
}

check_progress <- function(samplers, stage, iter, max_gd, min_es, min_unique, max_trys, step_size, progress = NULL){
  if(is.null(progress)){
    progress$n_remove <- step_size
    iters_total <- 0
    trys <- 0
  } else{
    iters_total <- progress$iters_total + step_size
    trys <- progress$trys + 1
  }
  gd <- check_gd(samplers, stage, progress$n_remove, max_gd)
  iter_done <- iters_total >= iter
  es_done <- min(es_pmwg(as_mcmc.list(samplers,selection="alpha",filter=stage))) > min_es
  step_size <- min(iter - iters_total, step_size)
  trys_done <- trys > max_trys
  done <- (es_done & iter_done & gd$gd_done) | max_trys
  return(list(samplers = gd$samplers, done = done, step_size = step_size, n_remove = n_remove))
}

check_gd <- function(samplers, stage, n_remove, max_gd){
  if(is.null(max_gd)){
    return(list(gd_done = TRUE, samplers = samplers))
  }
  gd <- gd_pmwg(as_mcmc.list(samplers,filter=stage), return_summary = FALSE,print_summary = FALSE,filter=stage,mapped=FALSE)
  samplers_short <- lapply(samplers,remove_iterations,select=n_remove,filter=stage)
  gd_short <- gd_pmwg(as_mcmc.list(samplers_short,filter=stage), filter=stage,mapped=FALSE)
  if (mean(gd_short) < mean(gd)) {
    gd <- gd_short
    samplers <- samplers_short
    n_remove <- step_size
  } else {
    n_remove <- round(n_remove + step_size/2)
  }
  if(all(is.finite(gd))) {
    ok_gd <- all(gd < max_gd)
    shorten <- !ok_gd
  } else {
    ok_gd <- FALSE
  })
  return(list(gd = gd, gd_done = ok_gd, n_remove = n_remove, samplers = samplers))
}

run_stage <- function(sampler, stage, iter=0, verbose,
                     particles=NA,particle_factor=50, p_accept= NULL,
                     epsilon = NULL, mix=NULL, n_cores=1)
{

  if (is.na(particles))
    particles <- round(particle_factor*sqrt(length(sampler$par_names)))
  if (!sampler$init) {
    sampler <- init(sampler, n_cores = n_cores, epsilon = epsilon)
  }
  if (iter == 0) return(sampler)
  if (verbose) message(paste0("Running ", stage, " stage"))
  sampler <- run_stage(sampler, stage = stage,iter = iter, particles = particles,
                       n_cores = n_cores, pstar = p_accept, epsilon = epsilon,
                       verbose = verbose, mix=mix)
  return(sampler)
}

create_eff_proposals <- function(samplers, n_cores){
  samples_merged <- merge_samples(samplers)
  test_samples <- extract_samples(samples_merged, stage = c("adapt", "sample"), max_sample_n)
  for(i in 1:length(samplers)){
    iteration = test_samples$iteration * i/length(samplers)
    conditionals <- parallel::mclapply(X = 1:samplers[[1]]$n_subjects,
                                       FUN = variant_funs$get_conditionals,samples = test_samples,
                                       samplers[[1]]$n_pars, iteration =  iteration,
                                       mc.cores = n_cores_conditional)
    conditionals <- array(unlist(conditionals), dim = c(samplers[[1]]$n_pars,
                                                        samplers[[1]]$n_pars + 1, samplers[[1]]$n_subjects))
    attr(samplers[[i]], "eff_mu") <- conditionals[,1,] #First column is the means
    attr(samplers[[i]], "eff_var") <- conditionals[,2:(samplers[[1]]$n_pars+1),] #Other columns are the variances
  }
  return(samplers)
}

create_cov_proposals <- function(samplers, samples_idx = NULL){
  get_covs <- function(sampler, samples_idx, sub){
    return(var(t(sampler$samples$alpha[,sub, samples_idx])))
  }
  n_pars <- samplers[[1]]$n_pars
  n_subjects <- samplers[[1]]$n_subjects
  n_chains <- length(samplers)
  if(is.null(samples_idx)){
    idx_subtract <- min(250, samplers[[1]]$samples$idx/2)
    samples_idx <- round(samplers[[1]]$samples$idx - idx_subtract):samplers[[1]]$samples$idx
  }
  preburned <- !is.null(attr(samplers[[1]], "chains_cov"))
  for(j in 1:n_chains){
    chains_cov <- array(NA_real_, dim = c(n_pars, n_pars, n_subjects))
    for(sub in 1:n_subjects){
      mean_covs <- get_covs(samplers[[j]], samples_idx, sub)
      # curr_sum <- sum(abs(mean_covs))
      # if(preburned){
      #   prev_sum <- sum(abs(attr(samplers[[j]], "chains_cov")[,,sub]))
      # } else{
      #   prev_sum <- curr_sum
      # }
      # mean_covs <- mean_covs/(curr_sum/prev_sum)
      if(is.negative.semi.definite(mean_covs)){
        chains_cov[,,sub] <- attr(samplers[[j]], "chains_cov")[,,sub]
      } else{
        chains_cov[,,sub] <-  as.matrix(nearPD(mean_covs)$mat)
      }
    }
    attr(samplers[[j]], "chains_cov") <- chains_cov
  }
  return(samplers)
}


