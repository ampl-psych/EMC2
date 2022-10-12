#### Fitting automation

run_stages <- function(sampler,iter=c(300,0,0),
                       verbose=FALSE,verbose_run_stage=FALSE,
                       particles=NA,particle_factor=50, p_accept= NULL, n_cores=1,
                       epsilon = NULL, start_mu = NULL, start_var = NULL, mix=NULL,
                       epsilon_upper_bound=15, eff_var = NULL, eff_mu = NULL)
  # Initializes (if needed) then runs burn, adapt and sample depending on iter
  # Verbose and verbose_run_stage are whether to keep you updated on EMC and PMwG output
  # Number of particles is set at particle_factor*sqrt(number of parameters) if particles is NA
  # P_accept is the target acceptance ratio for pmwg
  # Epsilon is a start value for the p_star tuning, best left at default
  # start_mu and start_var are startpoints for the sampler. If let NULL we'll sample startpoints
  # Epsilon_upper_bound is a safety net
  # Eff_var and Eff_mu are the efficient distibutions as calculated in run_adapt and run_sample

{

  if (is.na(particles))
    particles <- round(particle_factor*sqrt(length(sampler$par_names)))
  if (!sampler$init) {
    sampler <- init(sampler, n_cores = n_cores,
                    epsilon = epsilon, start_mu = start_mu, start_var = start_var)
  }
  if (all(iter==0)) return(sampler)
  if ( iter[1] != 0 ) {
    if (verbose) message("Running burn stage")
    sampler <- run_stage(sampler, stage = "burn",iter = iter[1], particles = particles,
                         n_cores = n_cores, pstar = p_accept, epsilon = epsilon, verbose = verbose_run_stage,
                         epsilon_upper_bound=epsilon_upper_bound, mix=mix)
    if (all(iter[2:3]==0)) return(sampler)
  }
  if (iter[2] != 0) {
    if (verbose) message("Running adapt stage")
    sampler <- run_stage(sampler, stage = "adapt",iter = iter[2], particles = particles,
                         n_cores = n_cores, pstar = p_accept, epsilon = epsilon, verbose = verbose_run_stage,
                         epsilon_upper_bound=epsilon_upper_bound, mix=mix)
    if (iter[3]==0) return(sampler)
  }
  if (verbose) message("Running sample stage")
  sampler <- run_stage(sampler, stage = "sample", iter = iter[3], epsilon = epsilon,
                       particles = particles, n_cores = n_cores, # pdist_update_n=pdist_update_n,
                       pstar = p_accept, verbose = verbose_run_stage, mix=mix, eff_mu = eff_mu,
                       eff_var = eff_var,  epsilon_upper_bound=epsilon_upper_bound)
  sampler
}



run_burn <- function(samplers,iter=300,
                     verbose=TRUE,verbose_run_stage=FALSE,mix=NULL,
                     particles=NA,particle_factor=50, p_accept= 0.7,
                     cores_per_chain=1,cores_for_chains=NULL,epsilon_upper_bound=15,
                     epsilon = NULL, start_mu = NULL, start_var = NULL)
  # This is the multi-chain version of run_chains but only for Burn
  # cores_per_chain is how many cores to give each chain
  # if cores_for chain > 1 will parallelize across chains
  # Initializes (if needed)
  # For the other arguments see run_stages
{
  if (is.null(cores_for_chains)) cores_for_chains <- length(samplers)
  data_list <- attr(samplers,"data_list")
  design_list <- attr(samplers,"design_list")
  model_list <- attr(samplers,"model_list")
  run_try <- 0
  repeat {
    samplers_new <- parallel::mclapply(samplers,run_stages,iter=c(iter, 0,0),
                             verbose=verbose,verbose_run_stage=verbose_run_stage,
                             particles=particles,particle_factor=particle_factor,
                             p_accept=p_accept, mix=mix, epsilon = epsilon, start_mu = start_mu, start_var = start_var,
                             epsilon_upper_bound=epsilon_upper_bound,
                             n_cores=cores_per_chain, mc.cores = cores_for_chains)
    if (class(samplers_new)=="try-error" ||
        any(lapply(samplers_new,class)=="try-error")) {
      save(samplers,iter,particles,particle_factor,p_accept, # pdist_update_n,
           epsilon_upper_bound,file="fail_run_stage.RData")
      run_try <- run_try + 1
      if (verbose) message("run_stage try error", run_try)
      if (run_try > 10)
        stop("Fail after 10 run_stage try errors, see saved samplers in fail_run_stage.RData")
    } else break
  }
  samplers <- samplers_new
  attr(samplers,"data_list") <- data_list
  attr(samplers,"design_list") <- design_list
  attr(samplers,"model_list") <- model_list
  samplers
}

# New version fail handling
run_gd <- function(samplers,iter=NA,max_trys=100,verbose=FALSE,burn=TRUE,
                   max_gd=1.1,thorough=TRUE, mapped=FALSE, shorten = TRUE,
                   epsilon = NULL, verbose_run_stage = FALSE,
                   particles=NA,particle_factor=50, p_accept=NULL,
                   cores_per_chain=1,cores_for_chains=NA,mix=NULL,
                   min_es=NULL,min_iter=NULL,max_iter=NULL)
  # Repeatedly runs burn or sample to get subject average multivariate
  # gelman.diag of alpha samples < max_gd (if !through) or if all univariate
  # psrf for every subject and parameter and multivariate < max_gd and
  # if min_es specified at least that effective size in worst case and if
  # min_iter specified at least that many iterations. If max_iter specified
  # pulls out after max_iter.
  # Trys adding iter (if NA 1/3 initial length) of the length and keeps all or
  # extra minus first n_remove (initially iter, then grows by iter/2 on each
  # expansion based on mean gd of alphas (mpsrf alone or with psrf depending on thorough)
  # until success or max_trys. Verbose prints out progress after each try
  # Cores used = cores_per_chain*cores_for_chains, latter set to number of chains by default
{

  enough_samples <- function(samplers,min_es,min_iter,max_iter,filter="burn") {
    if (!is.null(max_iter) && (sum(samplers[[1]]$samples$stage==filter) >= max_iter)) return(TRUE)
    if (is.null(min_iter)) ok <- TRUE else
      ok <- (sum(samplers[[1]]$samples$stage==filter)) > min_iter
    if (!is.null(min_es)) {
      es_min <- min(es_pmwg(as_mcmc.list(samplers,selection="alpha",filter=filter)))
      ok <- ok & (es_min > min_es)
      attr(ok,"es") <- es_min
    }
    ok
  }

  if (thorough) {
    max_name <- ", Max alpha mpsrf/psrf = "
    message("Exit on max alpha psrf/mpsrf < ",max_gd)
  } else {
    max_name <- ", Max alpha msrf = "
    message("Exit on max alpha mpsrf < ",max_gd)
  }
  n_chains <- length(samplers)
  n <- chain_n(samplers)
  if (is.na(cores_for_chains)) cores_for_chains <- n_chains
  if (is.na(iter)) iter <- round(n[,1][1]/3)
  n_remove <- iter
  # Iterate until criterion
  trys <- 0
  data_list <- attr(samplers,"data_list")
  design_list <- attr(samplers,"design_list")
  model_list <- attr(samplers,"model_list")
  gd <- gd_pmwg(as_mcmc.list(samplers,filter="burn"),!thorough,FALSE,
                filter="burn",mapped=mapped)
  if (all(is.finite(gd))) gd_diff <- apply(gd, 1, max) - 1.5*max_gd else gd_diff <- NA
  repeat {
    new_particle_n <- adaptive_particles(gd_diff, max_gd, particle_factor, particles)
    particle_factor <- new_particle_n$particle_factor
    particles <- new_particle_n$particles
    if(!any(is.na(gd_diff)) & any(gd_diff > .5)) samplers <- check_stuck(samplers) # Maybe a dumb heuristic
    samplers_new <- parallel::mclapply(samplers,run_stages,iter=c(iter,0,0),
                             n_cores=cores_per_chain,p_accept = p_accept, mix=mix,
                             particles=particles,particle_factor=particle_factor,epsilon=epsilon,
                             verbose=FALSE,verbose_run_stage=verbose_run_stage,
                             mc.cores=cores_for_chains)
    trys <- trys+1
    bad_new <- !class(samplers_new)=="list" || !all(unlist(lapply(samplers_new,class))=="pmwgs")
    if (bad_new)  gd <- matrix(Inf) else
      gd <- gd_pmwg(as_mcmc.list(samplers_new,filter="burn"),!thorough,FALSE,
                    filter="burn",mapped=mapped)
    if (all(is.finite(gd))) {
      samplers <- samplers_new
      if (shorten) {
        samplers_short <- lapply(samplers,remove_iterations,select=n_remove,filter="burn")
        if (!class(samplers_short)=="list" || !all(unlist(lapply(samplers_short,class))=="pmwgs"))
          gd_short <- matrix(Inf) else
            gd_short <- gd_pmwg(as_mcmc.list(samplers_short,filter="burn"),!thorough,FALSE,
                                filter="burn",mapped=mapped)
          if (mean(gd_short) < mean(gd)) {
            gd <- gd_short
            samplers <- samplers_short
            n_remove <- iter
          } else {
            n_remove <- round(n_remove + iter/2)
          }
      }
    } else {
      if (bad_new) {
        message("Failed to get new samples.")
        # bad_samplers <<- samplers
        # bad_samplers_new <<- samplers_new
      } else {
        message("Gelman diag try error.")
        samplers <- lapply(samplers_new,remove_iterations,select=n_remove,filter="burn")
      }
    }
    enough <- enough_samples(samplers,min_es,min_iter,max_iter,filter=filter)
    if (is.null(attr(enough,"es"))) es_message <- NULL else
      es_message <- paste(", Effective samples =",round(attr(enough,"es")))
    if (all(is.finite(gd))) {
      gd_diff <- (gd[,ncol(gd)] - 2*max_gd)
      ok_gd <- all(gd < max_gd)
      shorten <- !ok_gd
    } else ok_gd <- FALSE
    if (trys == max_trys || (ok_gd & enough)) {
      if (verbose) {
        message("\nFinal multivariate gelman.diag per participant")
        message("Iterations = ",samplers[[1]]$samples$idx,", Mean mpsrf= ",
                round(mean(gd),3),max_name,round(max(gd),3))
      }
      attr(samplers,"data_list") <- data_list
      attr(samplers,"design_list") <- design_list
      attr(samplers,"model_list") <- model_list
      if (ok_gd) attr(samplers,"burnt") <- max_gd else attr(samplers,"burnt") <- NA
      return(samplers)
    }
    if (verbose) {
      chain_n(samplers)[,"burn"][1]
      message(trys,": Iterations burn = ",chain_n(samplers)[,"burn"][1],", Mean mpsrf= ",
              round(mean(gd),3),max_name,round(max(gd),3),es_message)
    }
  }
}

# ndiscard=100;nstart=300;
#                       particles=NA; particle_factor = 50; start_mu = NULL; start_var = NULL;
#                       mix = NULL; verbose=TRUE;verbose_run_stage=FALSE;
#                       epsilon = NULL; epsilon_upper_bound=15; p_accept=0.7;
#                       cores_per_chain=1;cores_for_chains=NULL;
#                       max_gd_trys=100;max_gd=1.2;
#                       thorough=TRUE;mapped=FALSE; step_size = NA;
#                       min_es=NULL;min_iter=NULL;max_iter=NULL
auto_burn <- function(samplers,ndiscard=100,nstart=300,
                      particles=NA, particle_factor = 50, start_mu = NULL, start_var = NULL,
                      mix = NULL, verbose=TRUE,verbose_run_stage=FALSE,
                      epsilon = NULL, epsilon_upper_bound=15, p_accept=0.7,
                      cores_per_chain=1,cores_for_chains=NULL,
                      max_gd_trys=100,max_gd=1.2,
                      thorough=TRUE,mapped=FALSE, step_size = NA,
                      min_es=NULL,min_iter=NULL,max_iter=NULL
)
  # Will run burn until convergence through run_stages and run_gd
  # ndiscard and nstart together form the first batch of samples ran, of which ndiscard then discarded
  # For the first sets of arguments see run_burn and run_stages
  # Max_gd trys is how many times it will run gelman diags in fixed stepsizes
  # If stepsize not provided will use a heuristic in run_gd
  # If thorough gd should be lower than criterion for all parameters
  # If mapped parameters are transformed back to the 'interpretable' scale before gd is run.
  # Min_iter, min_es (effective sample size), max_iter and max_gd determine when run_gd is done

{
  if (is.null(cores_for_chains)) cores_for_chains <- length(samplers)
  data_list <- attr(samplers,"data_list")
  design_list <- attr(samplers,"design_list")
  model_list <- attr(samplers,"model_list")
  if (is.null(samplers[[1]]$samples$idx) || samplers[[1]]$samples$idx==1) {
    discard_message <- paste("(first ",ndiscard," then removed)")
    message("Getting initial ",ndiscard + nstart," samples ",discard_message)
    run_try <- 0
    samplers_new <- parallel::mclapply(samplers,run_stages,iter=c(nstart + ndiscard,0,0),
                             n_cores=cores_per_chain,p_accept = p_accept, mix=mix,
                             particles=particles,particle_factor=particle_factor,epsilon=epsilon,
                             verbose=FALSE,verbose_run_stage=verbose_run_stage,
                             mc.cores=cores_for_chains)
    samplers <- samplers_new
    if(ndiscard!= 0) samplers <- lapply(samplers,remove_iterations,select=ndiscard+1)
    message("Finished initial run")
    attr(samplers,"data_list") <- data_list
    attr(samplers,"design_list") <- design_list
    attr(samplers,"model_list") <- model_list
  }
  if (max_gd_trys==0) return(samplers)
  message("Beginning iterations to achieve Rhat < ",max_gd)
  samplers <- run_gd(samplers, max_trys=max_gd_trys,verbose=verbose,
         max_gd=max_gd,thorough=thorough, mapped=mapped, p_accept = p_accept,
         epsilon=epsilon, particles=particles,particle_factor=particle_factor, min_es=min_es,
         min_iter=min_iter,
         max_iter=max_iter, mix=mix, iter = step_size, verbose_run_stage = verbose_run_stage,
         cores_per_chain=cores_per_chain,cores_for_chains=cores_for_chains)
  return(samplers)
}


run_adapt <- function(samplers,max_trys=25,epsilon = NULL,
                      particles_adapt=NA,particle_factor_adapt=40, p_accept=.7,
                      cores_per_chain=1,cores_for_chains=NULL,mix=NULL,
                      n_cores_conditional = 1, min_unique = 200,
                      step_size_adapt = 25, thin = NULL,
                      verbose=TRUE,verbose_run_stage = FALSE){
  # Uses all the same arguments as run_stages and run_burn, but with:
  # max_trys, the amount of times to run adapt in step_size
  # After every step_size iterations it will check whether min_unique particles are there and adapt based on that
  # If the conditional distribution can be created it will finish adaptation and you can start sampling.
  # Thin whether to ONLY thin the samples passed to the creation of the conditional. Probably not necessary
  # Thin doesn't thin the actual ouptut
  if(verbose){
    message("Running adapt stage")
  }
  source(samplers[[1]]$source)
  if (is.null(cores_for_chains)) {
    cores_for_chains <- length(samplers)
  }
  data_list <- attr(samplers,"data_list")
  design_list <- attr(samplers,"design_list")
  model_list <- attr(samplers,"model_list")
  burnt <- attr(samplers,"burnt")
  trys <- 0
  samplers_new <- parallel::mclapply(samplers,run_stages,iter=c(0,min_unique/(length(samplers)*p_accept) - step_size_adapt,0),
                           n_cores=cores_per_chain,p_accept = p_accept, mix=mix,
                           particles=particles_adapt,particle_factor=particle_factor_adapt,epsilon=epsilon,
                           verbose=FALSE,verbose_run_stage=verbose_run_stage,
                           mc.cores=cores_for_chains)
  repeat {
    trys <- trys + 1
    samplers_new <- parallel::mclapply(samplers_new,run_stages,iter=c(0,step_size_adapt,0),
                             n_cores=cores_per_chain,p_accept = p_accept, mix=mix,
                             particles=particles_adapt,particle_factor=particle_factor_adapt,epsilon=epsilon,
                             verbose=FALSE,verbose_run_stage=verbose_run_stage,
                             mc.cores=cores_for_chains)
    samples_merged <- merge_samples(samplers_new)
    test_samples <- extract_samples(samples_merged, stage = "adapt", thin = thin, samples_merged$samples$idx, thin_eff_only = F)
    # Only need information like n_pars & n_subjects, so only need to pass the first chain
    adapted <- test_adapted(samplers_new[[1]], test_samples, min_unique, n_cores_conditional, verbose_run_stage)
    if(trys > max_trys | adapted) break
  }
  samplers <- samplers_new
  attr(samplers,"data_list") <- data_list
  attr(samplers,"design_list") <- design_list
  attr(samplers,"model_list") <- model_list
  attr(samplers,"burnt") <- burnt
  attr(samplers, "adapted") <- adapted
  return(samplers)
}

# verbose=TRUE;epsilon = NULL; particles_sample=NA;particle_factor_sample=25; p_accept=.7;
#                         cores_per_chain=1;cores_for_chains=NULL;mix=NULL;
#                         n_cores_conditional = 1; step_size_sample = 50; thin = NULL;
#                         verbose_run_stage = FALSE
# iter=1000
run_sample <- function(samplers,iter=NA,verbose=TRUE,
                       epsilon = NULL, particles_sample=NA,particle_factor_sample=25, p_accept=.7,
                       cores_per_chain=1,cores_for_chains=NULL,mix=NULL,
                       n_cores_conditional = 1, step_size_sample = 50, thin = NULL,
                       verbose_run_stage = FALSE)
  # Uses all the same arguments as run_stages and run_burn, but with:
  # iter: The amount of samples desired at the end from the sample stage
  # After step_size iterations it will update the conditional efficient distribution
  # Thin whether to ONLY thin the samples passed to the creation of the conditional. Probably not necessary
  # Thin doesn't thin the actual ouptut
{
  if(!attr(samplers, "adapted")){
    warning("samplers should be adapted before you can run sample stage")
    return(samplers)
  }
  if(verbose) message("Running sample stage")
  source(samplers[[1]]$source)
  if (is.null(cores_for_chains)) cores_for_chains <- length(samplers)
  data_list <- attr(samplers,"data_list")
  design_list <- attr(samplers,"design_list")
  model_list <- attr(samplers,"model_list")
  variant_funs <- attr(samplers[[1]], "variant_funs")
  burnt <- attr(samplers,"burnt")
  adapted <- attr(samplers,"adapted")
  samplers_new <- samplers
  n_steps <- ceiling(iter/step_size_sample)
  for(step in 1:n_steps){
    if(step == n_steps){
      step_size <- ifelse(iter %% step_size_sample == 0, step_size_sample, iter %% step_size_sample)
    }
    samples_merged <- merge_samples(samplers_new)
    test_samples <- extract_samples(samples_merged, stage = c("adapt", "sample"), thin = thin, samples_merged$samples$idx, thin_eff_only = F)

    conditionals <- parallel::mclapply(X = 1:samplers_new[[1]]$n_subjects,
                          FUN = variant_funs$get_conditionals,samples = test_samples,
                          samplers_new[[1]]$n_pars, mc.cores = n_cores_conditional)
    conditionals <- array(unlist(conditionals), dim = c(samplers_new[[1]]$n_pars,
                                                        samplers_new[[1]]$n_pars + 1, samplers_new[[1]]$n_subjects))
    eff_mu <- conditionals[,1,] #First column is the means
    eff_var <- conditionals[,2:(samplers_new[[1]]$n_pars+1),] #Other columns are the variances
    samplers_new <- parallel::mclapply(samplers_new,run_stages,iter=c(0,0,step_size_sample),
                             n_cores=cores_per_chain,p_accept = p_accept, mix=mix,
                             particles=particles_sample,particle_factor=particle_factor_sample,epsilon=epsilon,
                             verbose=FALSE,verbose_run_stage=verbose_run_stage, eff_mu = eff_mu,
                             eff_var = eff_var,
                             mc.cores=cores_for_chains)
    if(verbose) {
      cat(paste(step*step_size_sample," "))
      if (step %% 15 == 0) cat("\n")
    }

    # tmp <- run_stage(samplers_new[[2]],stage="sample",iter=50,epsilon = epsilon,
    #   particles = particles, n_cores = cores_per_chain,
    #    pstar = p_accept, verbose = TRUE, mix=mix, eff_mu = eff_mu,
    #    eff_var = eff_var,  epsilon_upper_bound=15)

  }
  if(verbose) cat("\n")
  samplers <- samplers_new
  attr(samplers,"data_list") <- data_list
  attr(samplers,"design_list") <- design_list
  attr(samplers,"model_list") <- model_list
  attr(samplers,"burnt") <- burnt
  attr(samplers, "adapted") <- adapted
  return(samplers)
}

run_IS2 <- function(samples, filter = "sample", subfilter = 0, IS_samples = 1000,
                    stepsize_particles = 500, max_particles = 5000, n_cores = 1, df = 5){
  variant <- basename(samples[[1]]$source)
  source(paste0("samplers/IS2/variants/", variant))
  samples_merged <- merge_samples(samples)
  IS2(samples_merged, filter, subfilter = subfilter, IS_samples, stepsize_particles, max_particles, n_cores, df)
}


test_adapted <- function(sampler, test_samples, min_unique, n_cores_conditional = 1,
                         verbose = FALSE)
{
  # Function used by run_adapt to check whether we can create the conditional.

  # Only need to check uniqueness for one parameter
  first_par <- test_samples$alpha[1, , ]
  # Split the matrix into a list of vectors by subject
  # Needed for the case where every sample is unique for all subjects
  first_par_list <- split(first_par, seq(NROW(first_par)))
  # Get unique pars (new accepted particles) and check length for
  # all subjects is greater than unq_vals
  n_unique_sub <- lapply(lapply(first_par_list, unique), length)
  n_pars <- sampler$n_pars
  variant_funs <- attr(sampler, "variant_funs")
  if (all(n_unique_sub > min_unique)) {
    if(verbose){
      message("Enough unique values detected: ", min_unique)
      message("Testing proposal distribution creation")
    }
    attempt <- tryCatch({
      parallel::mclapply(X = 1:sampler$n_subjects,FUN = variant_funs$get_conditionals,samples = test_samples,
               n_pars, mc.cores = n_cores_conditional)
    },error=function(e) e, warning=function(w) w)
    if (any(class(attempt) %in% c("warning", "error", "try-error"))) {
      if(verbose){
        message("A problem was encountered creating proposal distribution")
        message("Increasing required unique values and continuing adaptation")
      }
      return(FALSE)
    }
    else {
      if(verbose) message("Successfully adapted after ", test_samples$iteration, "iterations - stopping adaptation")
      return(TRUE)
    }
  } else{
    return(FALSE) # Not enough unique particles found
  }
}

# min_particles = 50; max_particles = 500;min_factor = 25; max_factor = 100; percent_up = 10; percent_down = 5
adaptive_particles <- function(gd_diff, max_gd, particle_factor = NA, particles = NA,
                               min_particles = 50, max_particles = 500,
                               min_factor = 25, max_factor = 100,
                               percent_up = 10, percent_down = 5){
  # This function adaptively tunes the particles per participant using broad,
  # so that we can lower the number of particles if we're closer to gd_criterion,
  # percent_up and down are relative to the max. Percent up is scaled by sqrt(gd_diff)
  if (any(is.na(gd_diff))) {
    if (is.na(particles))
      return(list(particles=100,particle_factor = particle_factor)) else
        return(list(particles=particles,particle_factor = particle_factor))
  }
  if(is.na(particles)){
    gd_diff[gd_diff > 0] <- sqrt(gd_diff[gd_diff > 0])*(percent_up/100)*max_factor
    gd_diff[gd_diff < 0] <- -(percent_down/100)*max_factor
    particle_factor <- pmin(pmax(min_factor, particle_factor + gd_diff), max_factor)
  } else{
    gd_diff[gd_diff > 0] <- sqrt(gd_diff[gd_diff > 0])*(percent_up/100)*max_particles
    gd_diff[gd_diff < 0] <- -(percent_down/100)*min_particles
    particles <- pmin(pmax(min_particles, particles + gd_diff), max_particles)
  }
  return(list(particles = round(particles), particle_factor = particle_factor))
}

check_stuck <- function(samples,filter=c("burn","sample")[1], # dont use adapt
                        start=1,last=TRUE,n=90, # n, from start or last n
                        flat=0.1,dfac=3)
  # flat: criterion on % change in median first to last 1/3 relative to last 1/3 Let me know what you think.

  # dfac: bad if  median(best)-median(chain) > dfac*IQR(best)
  # returns numeric(0) if none stuck, otherwise idices of stuck chains

{
  samplell <- function(sampler,filter,subfilter)
    sampler$samples$subj_ll[,sampler$samples$stage==filter][,subfilter]

  ns <- unlist(lapply(samples,function(x){sum(x$samples$stage==filter)}))
  if (!all(ns[1]==ns[-1])) stop("Filtered chains not of same length")
  if (n>ns[1]) stop("n is larger than available samples")
  if (last) start <- ns[1]-n + 1
  iter <- 1:n
  lls <- do.call(abind, list(lapply(samples,samplell,filter=filter,subfilter=iter+start-1), along = 3))
  first <- 1:(round(n/3))
  last <- round(2*n/3):n
  first <- apply(lls[,first,],c(1,3),median)
  last <- apply(lls[,last,],c(1,3),median)
  isFlat <- 100*abs((last-first)/last) < flat
  # if(any(isFlat)){
  diff <- apply(last, 1, FUN = function(x) return(max(x) - min(x)))
  IQRs <- apply(lls, c(1,3), FUN = IQR)
  best <- max.col(last)
  worst <- max.col(-last)
  n_subs <- samples[[1]]$n_subjects
  bad_subs <- which(diff > dfac*IQRs[matrix(c(1:n_subs, best), nrow = n_subs)]) # Yay R tricks
  if(any(bad_subs)) samples <- fix_stuck(samples, bad_subs, best, worst)
  # }

  return(samples)
}

std_error_IS2 <- function(IS_samples, n_bootstrap = 50000){
  log_marglik_boot= array(dim = n_bootstrap)
  for (i in 1:n_bootstrap){
    log_weight_boot = sample(IS_samples, length(IS_samples), replace = TRUE) #resample with replacement from the lw
    log_marglik_boot[i] <- median(log_weight_boot)
  }
  return(sd(log_marglik_boot))
}


fix_stuck <- function(samples, bad_subs, best, worst){
  # This function takes the bad subjects and replaces the last entry of the
  # worst chain with the best chain.
  idx <- samples[[1]]$samples$idx
  for(i in 1:length(bad_subs)){
    samples[[worst[i]]]$samples$alpha[,i,idx] <- samples[[best[i]]]$samples$alpha[,i,idx]
  }
  return(samples)
}



# cores_for_chains=1; verbose_run_stage=TRUE; cores_per_chain=28


run_emc <- function(file_name,nsample=1000, ...)
  # Combined auto fitting functions, getting and saving samples to disk.
  #    NB: samples must be first object in file_name file on disk
  # OR if file_name is not a character vector file_name is treated as a samplers
  #    object and results of fitting returned by the function.
{
  if (is.character(file_name)) {
    sname <- load(file_name)
    if (is.null(attr(get(sname[1]),"burnt")) || is.na(attr(get(sname[1]),"burnt"))) {
      assign(sname[1],auto_burn(get(sname[1]), ...))
      save(list=sname,file=file_name)
    }
    if (is.null(attr(get(sname[1]),"adapted")) && !is.null(attr(get(sname[1]),"burnt")) &&
        !is.na(attr(get(sname[1]),"burnt"))) {
      assign(sname[1],run_adapt(get(sname[1]), ...))
      save(list=sname,file=file_name)
    }
    if (!is.null(attr(get(sname[1]),"adapted")) && !is.null(attr(get(sname[1]),"adapted")) &&
        attr(get(sname[1]),"adapted") && nsample > 0) {
      assign(sname[1],run_sample(get(sname[1]),iter=nsample, ...))
      save(list=sname,file=file_name)
    }
  } else { # file_name is actually a samplers object
    if (is.null(attr(file_name,"burnt")) || is.na(attr(file_name,"burnt")))
      file_name <- auto_burn(file_name, ...)
    if (is.null(attr(file_name,"adapted")) && !is.null(attr(file_name,"burnt")) &&
        !is.na(attr(file_name,"burnt")))
      file_name <- run_adapt(file_name, ...)
    if (!is.null(attr(file_name,"adapted")) && !is.null(attr(file_name,"adapted")) &&
        attr(file_name,"adapted") && nsample > 0)
      file_name <- run_sample(file_name,iter=nsample, ...)
    return(file_name)
  }
}
