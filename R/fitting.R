#### Fitting automation
#' Generic function to run samplers with default settings.
#'
#' Calls `auto_burn`, `run_adapt` and `run_sample` in order with default settings.
#'
#' @param samplers A list of samplers, or a fileName of where the samplers are stored.
#' @param stage A string. Indicates which stage is to be run, either preburn, burn, adapt or sample. If unspecified will assume the next unrun stage.
#' @param iter An integer. Indicates how many iterations to run sampling stage.
#' @param max_gd A double. The maximum gelman diagnostic convergence allowed. Will stop sample stage if below this number.
#' @param mean_gd A double. The mean gelman diagnostic convergence allowed. Will stop burn-in if below this number.
#' @param min_es An integer. The minimal effective size required to stop sampling.
#' @param min_unique An integer. The minimal number of samples required. Only works in adaptation.
#' @param preburn An integer. Specifies how many iterations to run preburn stage.
#' @param p_accept A double. The target acceptance probability of the MCMC process. This will fine tune the width of the search space. Default = .8
#' @param step_size An integer. After each of these steps, the requirements will be checked if they are met and proposal distributions will be updated. Default = 100.
#' @param verbose Logical. Whether to print emc related messages
#' @param verboseProgress Logical. Whether to print sampling related messages
#' @param fileName A string. If specified will autosave samplers at this location.
#' @param particles An integer. How many particles to use, default is NULL and particle_factor is used. If specified will override particle_factor
#' @param particle_factor An integer. Particle factor multiplied by the square root of the number of sampled parameters will determine the number of particles used.
#' @param cores_per_chain An integer. How many cores to use per chain. Parallelizes across participant calculations.
#' @param cores_for_chains An integer. How many cores to use across chains. Default is the number of chains.
#' @param max_trys An integer. How many times it will try to meet the finish conditions. Default is 50.
#'
#' @return A list of samplers
#' @export
#'
#' @examples
run_emc <- function(samplers, stage = NULL, iter = 1000, max_gd = 1.1, mean_gd = 1.1, min_es = 0, min_unique = 600, preburn = 150,
                    p_accept = .8, step_size = 100, verbose = FALSE, verboseProgress = FALSE, fileName = NULL,
                    particles = NULL, particle_factor = 50, cores_per_chain = 1,
                    cores_for_chains = length(samplers), max_trys = 50, useC = FALSE, n_blocks = NULL){
  if (is.character(samplers)) {
    samplers <- fix_fileName(samplers)
    if(is.null(fileName)) fileName <- samplers
    samplers <- loadRData(samplers)
  }
  if(is.null(stage)){
    nstage <- colSums(chain_n(samplers))
    if (nstage["preburn"]==0) stage <- "preburn" else
      if (nstage["sample"]>0) stage <- "sample" else
       stage <- names(nstage)[which(nstage==0)[1]-1]
    }
  if(stage == "preburn"){
    samplers <- run_samplers(samplers, stage = "preburn", iter = preburn, cores_for_chains = cores_for_chains, p_accept = p_accept,
                             step_size = step_size,  verbose = verbose, verboseProgress = verboseProgress,
                             fileName = fileName,
                             particles = particles, particle_factor =  particle_factor,
                             cores_per_chain = cores_per_chain, max_trys = max_trys, useC = useC, n_blocks = n_blocks)
  }
  if(any(stage %in% c("preburn", "burn"))){
    samplers <-  run_samplers(samplers, stage = "burn", mean_gd = mean_gd, cores_for_chains = cores_for_chains, p_accept = p_accept,
                              step_size = step_size,  verbose = verbose, verboseProgress = verboseProgress,
                              fileName = fileName,
                              particles = particles, particle_factor =  particle_factor,
                              cores_per_chain = cores_per_chain, max_trys = max_trys, useC = useC, n_blocks = n_blocks)
  }
  if(any(stage %in% c("preburn", "burn", "adapt"))){
    samplers <-  run_samplers(samplers, stage = "adapt", min_unique = min_unique, cores_for_chains = cores_for_chains, p_accept = p_accept,
                              step_size = step_size,  verbose = verbose, verboseProgress = verboseProgress,
                              fileName = fileName,
                              particles = particles, particle_factor =  particle_factor,
                              cores_per_chain = cores_per_chain, max_trys = max_trys, useC = useC, n_blocks = n_blocks)
  }
  if(any(stage %in% c("preburn", "burn", "adapt", "sample")) ){
    samplers <-  run_samplers(samplers, stage = "sample", iter = iter, max_gd = max_gd, cores_for_chains = cores_for_chains, p_accept = p_accept,
                              step_size = step_size,  verbose = verbose, verboseProgress = verboseProgress,
                              fileName = fileName,
                              particles = particles, particle_factor = particle_factor,
                              cores_per_chain = cores_per_chain, max_trys = max_trys, useC = useC, n_blocks = n_blocks)
  }
  return(samplers)
}

#' Generic function to run samplers for any stage and any requirements.
#'
#' Used by `run_emc`, `auto_burn`, `run_adapt` and `run_sample`.
#' Will break if you skip a stage, the stages have to be run in order (preburn, burn, adapt, sample).
#' Either iter, max_gd, min_es or min_unique has to be specified. Multiple conditions for finishing can be specified. Will finish if all conditions are met.
#' @param samplers A list of samplers, could be in any stage, as long as they've been initialized with make_samplers
#' @param stage A string. Indicates which stage is to be run, either preburn, burn, adapt or sample
#' @param iter An integer. Indicates how many iterations to run,
#' @param max_gd A double. The maximum gelman diagnostic convergence allowed. Will stop if below this number.
#' @param mean_gd A double. The mean gelman diagnostic convergence allowed. Will stop if below this number.
#' @param min_es An integer. The minimal effective size required.
#' @param min_unique An integer. The minimal number of samples required. Only works in adaptation.
#' @param p_accept A double. The target acceptance probability of the MCMC process. This will fine tune the width of the search space. Default = .8
#' @param step_size An integer. After each of these steps, the requirements will be checked if they are met and proposal distributions will be updated. Default = 100.
#' @param verbose Logical. Whether to print emc related messages
#' @param verboseProgress Logical. Whether to print sampling related messages
#' @param fileName A string. If specified will autosave samplers at this location.
#' @param particles An integer. How many particles to use, default is NULL and particle_factor is used. If specified will override particle_factor
#' @param particle_factor An integer. Particle factor multiplied by the square root of the number of sampled parameters will determine the number of particles used.
#' @param cores_per_chain An integer. How many cores to use per chain. Parallelizes across participant calculations.
#' @param cores_for_chains An integer. How many cores to use across chains. Default is the number of chains.
#' @param max_trys An integer. How many times it will try to meet the finish conditions. Default is 50.
#'
#' @return A list of samplers
#' @export
#'
#' @examples
run_samplers <- function(samplers, stage, iter = NULL, max_gd = NULL, mean_gd = NULL, min_es = 0, min_unique = 600,
                         p_accept = .8, step_size = 100, verbose = FALSE, verboseProgress = FALSE,
                         fileName = NULL,
                         particles = NULL, particle_factor = 50, cores_per_chain = 1,
                         cores_for_chains = length(samplers), max_trys = 50, useC = FALSE, n_blocks = NULL){
  if (verbose) message(paste0("Running ", stage, " stage"))
  attributes <- get_attributes(samplers)
  total_iters_stage <- chain_n(samplers)[,stage][1]
  iter <- iter + total_iters_stage
  progress <- check_progress(samplers, stage, iter, max_gd, mean_gd, min_es, min_unique, max_trys, step_size, cores_per_chain, verbose)
  samplers <- progress$samplers
  while(!progress$done){
    if(!is.numeric(progress$step_size) | progress$step_size < 1) warning("Something wrong with the stepsize again, Niek's to blame")
    samplers <- add_proposals(samplers, stage, cores_per_chain, n_blocks)
    samplers <- parallel::mclapply(samplers,run_stages, stage = stage, iter= progress$step_size,
                                   verbose=verbose,  verboseProgress = verboseProgress,
                                   particles=particles,particle_factor=particle_factor,
                                   p_accept=p_accept, n_cores=cores_per_chain, useC = useC, mc.cores = cores_for_chains)
    progress <- check_progress(samplers, stage, iter, max_gd, mean_gd, min_es, min_unique, max_trys, step_size, cores_per_chain,
                               verbose, progress)
    samplers <- progress$samplers
    if(!is.null(fileName)){
      fileName <- fix_fileName(fileName)
      save(samplers, file = fileName)
    }
  }
  samplers <- get_attributes(samplers, attributes)
  return(samplers)
}

run_stages <- function(sampler, stage = "preburn", iter=0, verbose = TRUE, verboseProgress = TRUE,
                       particles=NULL,particle_factor=50, p_accept= NULL, n_cores=1, useC = FALSE)
{

  if (is.null(particles))
    particles <- round(particle_factor*sqrt(length(sampler$par_names)))
  if (!sampler$init) {
    sampler <- init(sampler, n_cores = n_cores, useC = useC)
  }
  if (iter == 0) return(sampler)
  sampler <- run_stage(sampler, stage = stage,iter = iter, particles = particles,
                       n_cores = n_cores, p_accept = p_accept, verbose = verbose, verboseProgress = verboseProgress,
                       useC = useC)
  return(sampler)
}

add_proposals <- function(samplers, stage, n_cores, n_blocks){
  if(stage != "preburn"){
    samplers <- create_cov_proposals(samplers)
    if(!is.null(n_blocks)){
      components <- sub_blocking(samplers, n_blocks)
      for(i in 1:length(samplers)){
        attr(samplers[[i]]$data, "components") <- components
      }
    }
  }
  if(stage == "sample"){
    samplers <- create_eff_proposals(samplers, n_cores)
  }
  return(samplers)
}

check_progress <- function (samplers, stage, iter, max_gd, mean_gd, min_es, min_unique,
                            max_trys, step_size, n_cores, verbose, progress = NULL)
{
  total_iters_stage <- chain_n(samplers)[, stage][1]
  if (is.null(progress)) {
    iters_total <- 0
    trys <- 0
  }
  else {
    iters_total <- progress$iters_total + step_size
    trys <- progress$trys + 1
    if (verbose)
      message(trys, ": Iterations ", stage, " = ", total_iters_stage)
  }
  gd <- check_gd(samplers, stage, max_gd, mean_gd, trys, verbose)
  iter_done <- ifelse(is.null(iter) || length(iter) == 0, TRUE, total_iters_stage >= iter)
  if (min_es == 0) {
    es_done <- TRUE
  }
  else if (iters_total != 0) {
    curr_min_es <- min(es_pmwg(as_mcmc.list(samplers, selection = "alpha",
                                            filter = stage), print_summary = F))
    if (verbose)
      message("Smallest effective size = ", round(curr_min_es))
    es_done <- ifelse(!samplers[[1]]$init, FALSE, curr_min_es >
                        min_es)
  }
  else {
    es_done <- FALSE
  }
  trys_done <- ifelse(is.null(max_trys), FALSE, trys >= max_trys)
  if (trys_done) {
    warning("Max trys reached. If this happens in burn-in while trying to get gelman diagnostics < 1.2, you might have a particularly hard model. Make sure your model is well specified. If so, you can run adapt and sample, if run for long enough, sample usually converges eventually.")
  }
  if (stage == "adapt") {
    samples_merged <- merge_samples(samplers)
    test_samples <- extract_samples(samples_merged, stage = "adapt",
                                    samples_merged$samples$idx)
    adapted <- test_adapted(samplers[[1]], test_samples,
                            min_unique, n_cores, verbose)
  }
  else {
    adapted <- TRUE
  }
  done <- (es_done & iter_done & gd$gd_done & adapted) | trys_done
  if(es_done & gd$gd_done & adapted & !iter_done){
    step_size <- min(step_size, abs(iter - total_iters_stage))[1]
  }
  return(list(samplers = gd$samplers, done = done, step_size = step_size,
              trys = trys, iters_total = iters_total))
}

check_gd <- function(samplers, stage, max_gd, mean_gd, trys, verbose){
  if(is.null(max_gd) & is.null(mean_gd)) return(list(gd_done = TRUE, samplers = samplers))
  if(!samplers[[1]]$init | !stage %in% samplers[[1]]$samples$stage) return(list(gd_done = FALSE, samplers = samplers))
  gd <- gd_pmwg(as_mcmc.list(samplers,filter=stage), return_summary = FALSE,print_summary = FALSE,filter=stage,mapped=FALSE)
  n_remove <- round(chain_n(samplers)[,stage][1]/3)
  samplers_short <- try(lapply(samplers,remove_iterations,select=n_remove,filter=stage),silent=TRUE)
  if (is(samplers_short,"try-error")) gd_short <- Inf else
    gd_short <- gd_pmwg(as_mcmc.list(samplers_short,filter=stage), return_summary = FALSE,
                        print_summary = FALSE, filter=stage,mapped=FALSE)
  if (mean(gd_short) < mean(gd)) {
    gd <- gd_short
    samplers <- samplers_short
  }
  if(!is.null(max_gd)){
    ok_max_gd <- ifelse(all(is.finite(gd)), all(gd < max_gd), FALSE)
  } else{
    ok_max_gd <- TRUE
  }
  if(!is.null(mean_gd)){
    ok_mean_gd <- ifelse(all(is.finite(gd)), mean(gd) < mean_gd, FALSE)
  } else{
    ok_mean_gd <- TRUE
  }
  ok_gd <- ok_max_gd & ok_mean_gd

  if(verbose){
    message("Mean mpsrf= ",round(mean(gd),3),", Max alpha mpsrf/psrf = ",round(max(gd),3))
  }
  return(list(gd = gd, gd_done = ok_gd, samplers = samplers))
}


create_eff_proposals <- function(samplers, n_cores){
  samples_merged <- merge_samples(samplers)
  test_samples <- extract_samples(samples_merged, stage = c("adapt", "sample"), max_n_sample = 1000)
  variant_funs <- attr(samplers[[1]], "variant_funs")
  components <- attr(samplers[[1]]$data, "components")
  for(i in 1:length(samplers)){
    iteration = round(test_samples$iteration * i/length(samplers))
    n_pars <- samplers[[1]]$n_pars
    n_subjects <- samplers[[1]]$n_subjects
    eff_mu <- matrix(0, nrow = n_pars, ncol = n_subjects)
    eff_var <- array(0, dim = c(n_pars, n_pars, n_subjects))
    for(comp in unique(components)){
      idx <- comp == components
      conditionals <- parallel::mclapply(X = 1:n_subjects,
                                         FUN = variant_funs$get_conditionals,samples = test_samples,
                                         n_pars = sum(idx), iteration =  iteration, idx = idx,
                                         mc.cores = n_cores)
      conditionals <- array(unlist(conditionals), dim = c(sum(idx), sum(idx) + 1, n_subjects))
      eff_mu[idx,] <- conditionals[,1,]
      eff_var[idx,idx,] <- conditionals[,2:(sum(idx)+1),]
    }

    attr(samplers[[i]], "eff_mu") <- eff_mu
    attr(samplers[[i]], "eff_var") <- eff_var
  }
  return(samplers)
}

sub_blocking <- function(samplers, n_blocks){
  covs <- lapply(samplers, FUN = function(x){return(attr(x, "chains_cov"))})
  out <- array(0, dim = dim(covs[[1]][,,1]))
  for(i in 1:length(covs)){
    cov_tmp <- covs[[1]]
    for(j in 1:dim(cov_tmp)[3]){
      out <- out + cov2cor(cov_tmp[,,j])
    }
  }
  shared_ll_idx <- attr(samplers[[1]]$data, "shared_ll_idx")
  min_comp <- 0
  components <- c()
  for(ll in unique(shared_ll_idx)){
    idx <- ll == shared_ll_idx
    distance <-as.dist(1- abs(out[idx, idx]/out[1,1]))
    clusts <- hclust(distance)
    sub_comps <- min_comp + cutree(clusts, k = n_blocks) # This could go wrong if one group has just one member
    min_comp <- max(sub_comps)
    components <- c(components, sub_comps)
  }
  return(components)
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
  components <- attr(samplers[[1]]$data, "components")
  block_idx <- block_variance_idx(components)
  for(j in 1:n_chains){
    chains_cov <- array(NA_real_, dim = c(n_pars, n_pars, n_subjects))
    for(sub in 1:n_subjects){
      mean_covs <- get_covs(samplers[[j]], samples_idx, sub)
      mean_covs[block_idx] <- 0
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

get_attributes <- function(samplers, attributes = NULL){
  if(is.null(attributes)) {
    return(list(  data_list = attr(samplers,"data_list"),
                  design_list = attr(samplers,"design_list"),
                  model_list = attr(samplers,"model_list")))
  } else{
    attr(samplers,"data_list") <- attributes$data_list
    attr(samplers,"design_list") <- attributes$design_list
    attr(samplers,"model_list") <- attributes$model_list
    return(samplers)
  }
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
  components <- attr(sampler$data, "components")

  if (all(n_unique_sub > min_unique)) {
    if(verbose){
      message("Enough unique values detected: ", min_unique)
      message("Testing proposal distribution creation")
    }
    attempt <- tryCatch({
        for(comp in unique(components)){
          idx <- comp == components
          parallel::mclapply(X = 1:sampler$n_subjects,FUN = variant_funs$get_conditionals,samples = test_samples,
                         n_pars = sum(idx), idx = idx, mc.cores = n_cores_conditional)
        }
    },error=function(e) e, warning=function(w) w)
    if (any(class(attempt) %in% c("warning", "error", "try-error"))) {
      if(verbose){
        message("Can't create efficient distribution yet")
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

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#' Runs burn-in for samplers.
#'
#' Special instance of `run_samplers`, with default arguments specified for completing burn_in.
#' Will run both preburn and burn.
#'
#' @param samplers A list of samplers, could be in any stage, as long as they've been initialized with make_samplers
#' @param max_gd A double. The maximum gelman diagnostic convergence allowed. Will stop if below this number.
#' @param mean_gd A double. The mean gelman diagnostic convergence allowed. Will stop if below this number.
#' @param min_es An integer. The minimal effective size required.
#' @param preburn An integer. The number of iterations run for preburn stage.
#' @param p_accept A double. The target acceptance probability of the MCMC process. This will fine tune the width of the search space. Default = .8
#' @param step_size An integer. After each of these steps, the requirements will be checked if they are met and proposal distributions will be updated. Default = 100.
#' @param verbose Logical. Whether to print emc related messages
#' @param verboseProgress Logical. Whether to print sampling related messages
#' @param fileName A string. If specified will autosave samplers at this location.
#' @param particles An integer. How many particles to use, default is NULL and particle_factor is used. If specified will override particle_factor
#' @param particle_factor An integer. Particle factor multiplied by the square root of the number of sampled parameters will determine the number of particles used.
#' @param cores_per_chain An integer. How many cores to use per chain. Parallelizes across participant calculations.
#' @param cores_for_chains An integer. How many cores to use across chains. Default is the number of chains.
#' @param max_trys An integer. How many times it will try to meet the finish conditions. Default is 50.
#'
#' @return A list of samplers
#' @export
#'
#' @examples
auto_burn <- function(samplers, max_gd = NULL, mean_gd = 1.1, min_es = 0, preburn = 150,
                      p_accept = .8, step_size = 100, verbose = FALSE, verboseProgress = FALSE,
                      fileName = NULL,
                      particles = NULL, particle_factor = 50, cores_per_chain = 1,
                      cores_for_chains = length(samplers), max_trys = 50, useC = FALSE, n_blocks = NULL){
  samplers <- run_samplers(samplers, stage = "preburn", iter = preburn, cores_for_chains = cores_for_chains, p_accept = p_accept,
                           step_size = step_size,  verbose = verbose, verboseProgress = verboseProgress,
                           fileName = fileName,
                           particles = particles, particle_factor =  particle_factor,
                           cores_per_chain = cores_per_chain, max_trys = max_trys, useC = useC, n_blocks = n_blocks)
  samplers <-  run_samplers(samplers, stage = "burn", max_gd = max_gd, mean_gd = mean_gd, min_es = min_es, cores_for_chains = cores_for_chains, p_accept = p_accept,
                            step_size = step_size,  verbose = verbose, verboseProgress = verboseProgress,
                            fileName = fileName,
                            particles = particles, particle_factor =  particle_factor,
                            cores_per_chain = cores_per_chain, max_trys = max_trys, useC = useC, n_blocks = n_blocks)
  return(samplers)
}
#' Runs adapt stage for samplers.
#'
#' Special instance of `run_samplers`, with default arguments specified for completing adaptation.
#'
#' @param samplers A list of samplers, could be in any stage, as long as they've been initialized with make_samplers
#' @param max_gd A double. The maximum gelman diagnostic convergence allowed. Will stop if below this number.
#' @param mean_gd A double. The mean gelman diagnostic convergence allowed. Will stop if below this number.
#' @param min_es An integer. The minimal effective size required.
#' @param min_unique An integer. The minimal number of samples required.
#' @param p_accept A double. The target acceptance probability of the MCMC process. This will fine tune the width of the search space. Default = .8
#' @param step_size An integer. After each of these steps, the requirements will be checked if they are met and proposal distributions will be updated. Default = 100.
#' @param verbose Logical. Whether to print emc related messages
#' @param verboseProgress Logical. Whether to print sampling related messages
#' @param fileName A string. If specified will autosave samplers at this location.
#' @param particles An integer. How many particles to use, default is NULL and particle_factor is used. If specified will override particle_factor
#' @param particle_factor An integer. Particle factor multiplied by the square root of the number of sampled parameters will determine the number of particles used.
#' @param cores_per_chain An integer. How many cores to use per chain. Parallelizes across participant calculations.
#' @param cores_for_chains An integer. How many cores to use across chains. Default is the number of chains.
#' @param max_trys An integer. How many times it will try to meet the finish conditions. Default is 50.
#'
#' @return A list of samplers.
#' @export
#'
#' @examples
run_adapt <- function(samplers, max_gd = NULL, mean_gd = NULL, min_es = 0, min_unique = 600,
                      p_accept = .8, step_size = 100, verbose = FALSE, verboseProgress = FALSE,
                      fileName = NULL,
                      particles = NULL, particle_factor = 50, cores_per_chain = 1,
                      cores_for_chains = length(samplers), max_trys = 50, useC = FALSE, n_blocks = NULL){
  samplers <- run_samplers(samplers, stage = "adapt",  max_gd = max_gd, mean_gd = mean_gd, min_es = min_es, min_unique = min_unique,
                           cores_for_chains = cores_for_chains, p_accept = p_accept,
                           step_size = step_size,  verbose = verbose, verboseProgress = verboseProgress,
                           fileName = fileName,
                           particles = particles, particle_factor =  particle_factor,
                           cores_per_chain = cores_per_chain, max_trys = max_trys, useC = useC, n_blocks = n_blocks)
  return(samplers)
}
#' Runs sample stage for samplers.
#'
#' Special instance of `run_samplers`, with default arguments specified for running sample stage.
#'
#' @param samplers A list of samplers, could be in any stage, as long as they've been initialized with make_samplers
#' @param iter An integer. Indicates how many iterations to run,
#' @param max_gd A double. The maximum gelman diagnostic convergence allowed. Will stop if below this number.
#' @param mean_gd A double. The mean gelman diagnostic convergence allowed. Will stop if below this number.
#' @param min_es An integer. The minimal effective size required.
#' @param p_accept A double. The target acceptance probability of the MCMC process. This will fine tune the width of the search space. Default = .8
#' @param step_size An integer. After each of these steps, the requirements will be checked if they are met and proposal distributions will be updated. Default = 100.
#' @param verbose Logical. Whether to print emc related messages
#' @param verboseProgress Logical. Whether to print sampling related messages
#' @param fileName A string. If specified will autosave samplers at this location.
#' @param particles An integer. How many particles to use, default is NULL and particle_factor is used. If specified will override particle_factor
#' @param particle_factor An integer. Particle factor multiplied by the square root of the number of sampled parameters will determine the number of particles used.
#' @param cores_per_chain An integer. How many cores to use per chain. Parallelizes across participant calculations.
#' @param cores_for_chains An integer. How many cores to use across chains. Default is the number of chains.
#' @param max_trys An integer. How many times it will try to meet the finish conditions. Default is 50.
#' @export
#'
#' @return A list of samplers
run_sample <- function(samplers, iter = 1000, max_gd = 1.1, mean_gd = NULL, min_es = 0,
                       p_accept = .8, step_size = 100, verbose = FALSE, verboseProgress = FALSE,
                       fileName = NULL,
                       particles = NULL, particle_factor = 50, cores_per_chain = 1,
                       cores_for_chains = length(samplers), max_trys = 50, useC = FALSE, n_blocks = NULL){
  samplers <- run_samplers(samplers, stage = "sample", iter = iter, max_gd = max_gd, mean_gd = mean_gd, min_es = min_es, cores_for_chains = cores_for_chains, p_accept = p_accept,
                           step_size = step_size,  verbose = verbose, verboseProgress = verboseProgress,
                           fileName = fileName,
                           particles = particles, particle_factor =  particle_factor,
                           cores_per_chain = cores_per_chain, max_trys = max_trys, useC = useC, n_blocks = n_blocks)
  return(samplers)
}

#' Creates a sampler from data and design combined.
#'
#' This function is used to initialize samplers using the data and the prespecified design.
#'
#' @param data_list A dataframe of data, or a list of a dataframe. Should have the column subjects as identifiers.
#' @param design_list A list with a prespecified design made with make_design.
#' @param model_list A model list, if empty will use the model specified in the design_list.
#' @param type A string indicating whether to run a standard group-level, or blocked, diagonal, factor, or single.
#' @param n_chains An integer. Specifies the amount of mcmc chains to be run. Should be more than 1 to get gelman diagnostics.
#' @param rt_resolution A double. Used for compression, rts will be binned based on this resolution.
#' @param prior_list A list of priors for the group level. Prior distributions should match the type argument.
#' @param par_groups A vector. Only to be specified with type blocked `c(1,1,1,2,2)` means first three parameters first block, last two parameters in the second block
#' @param n_factors An integer. Only to be specified with type factor.
#' @param constraintMat A matrix of rows equal to the number of estimated parameters, and columns equal to the number of factors, only to be specified with type factor.
#' If null will use default settings as specified in Innes et al. 2022
#'
#' @return a list of samplers
#' @export
#'
#' @examples
#'
make_samplers <- function(data_list,design_list,model_list=NULL,
                          type=c("standard","diagonal","blocked","factor","single", "lm", "infnt_factor")[1],
                          n_chains=3,rt_resolution=0.02,
                          prior_list = NULL,
                          par_groups=NULL,
                          n_factors=NULL,constraintMat = NULL, formula = NULL)

{
  if (!(type %in% c("standard","diagonal","blocked","factor","single", "lm", "infnt_factor")))
    stop("type must be one of: standard,diagonal,blocked,factor,infnt_factor", "lm","single")
  if (is(data_list, "data.frame")) data_list <- list(data_list)
  # Sort subject together and add unique trial within subject integer
  # create overarching data list with one list per subject
  if(type == "lm"){
    if(length(data_list) > 1) stop("no joint models for lm yet")
    vars <- c()
    for(form in formula){
      vars <- c(vars, split_form(form)$dep)
    }
    tmp <- data_list[[1]]
    aggr_data <- tmp[cumsum(table(tmp$subjects)),c("subjects", unique(vars))]
    for(i in 1:ncol(aggr_data)){
      if(colnames(aggr_data)[i] != "subjects" & is.factor(aggr_data[,i])){
        aggr_data[,i] <- factor(aggr_data[,i], levels = unique(aggr_data[,i]))
      }
    }
  }
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
    if(is.null(attr(design_list[[i]], "custom_ll"))){
      dadm_list[[i]] <- design_model(data=data_list[[i]],design=design_list[[i]],
                                   model=model_list[[i]],rt_resolution=rt_resolution[i],prior=prior_list[[i]])
    } else{
      dadm_list[[i]] <- design_model_custom_ll(data = data_list[[i]], design = design_list[[i]],
                                               model=model_list[[i]], prior=prior_list[[i]])
    }
  }

  # if(!is.null(subject_covariates)) attr(dadm_list, "subject_covariates") <- subject_covariates
  variant_funs <- get_variant_funs(type = type)
  if (type %in% c("standard", "single", "diagonal", "infnt_factor")) {
    out <- pmwgs(dadm_list, variant_funs)
  } else if (type == "blocked") {
    if (is.null(par_groups)) stop("Must specify par_groups for blocked type")
    out <- pmwgs(dadm_list,par_groups=par_groups, variant_funs)
  } else if (type == "lm") {
    out <- pmwgs(dadm_list,variant_funs, formula = formula, aggr_data = aggr_data)
  }
  # replicate chains
  dadm_lists <- rep(list(out),n_chains)
  # For post predict
  attr(dadm_lists,"data_list") <- data_list
  attr(dadm_lists,"design_list") <- design_list
  attr(dadm_lists,"model_list") <- model_list
  return(dadm_lists)
}

fix_fileName <- function(x){
  ext <- substr(x, nchar(x)-5, nchar(x))
  if(ext != ".RData" & ext != ".Rdata"){
    return(paste0(x, ".RData"))
  } else{
    return(x)
  }
}



extractDadms <- function(dadms, names = 1:length(dadms)){
  N_models <- length(dadms)
  pars <- attr(dadms[[1]], "sampled_p_names")
  prior <- attr(dadms[[1]], "prior")
  ll_func <- attr(dadms[[1]], "model")()$log_likelihood
  subjects <- unique(factor(sapply(dadms, FUN = function(x) levels(x$subjects))))
  dadm_list <- dm_list(dadms[[1]])
  components <- rep(1, length(pars))
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
      components <- c(components, rep(k, length(curr_pars)))
      pars <- c(pars, paste(names[k], curr_pars, sep = "|"))
      prior$theta_mu_mean <- c(prior$theta_mu_mean, attr(dadm, "prior")$theta_mu_mean)
      if(is.matrix(prior$theta_mu_var)){
        prior$theta_mu_var <- adiag(prior$theta_mu_var, attr(dadm, "prior")$theta_mu_var)
      } else{
        prior$theta_mu_var <- c(prior$theta_mu_var, attr(dadm, "prior")$theta_mu_var)
      }
    }
    ll_func <- log_likelihood_joint
    dadm_list <- do.call(mapply, c(list, total_dadm_list, SIMPLIFY = F))
  }
  # subject_covariates_ok <- unlist(lapply(subject_covariates, FUN = function(x) length(x) == length(subjects)))
  # if(!is.null(subject_covariates_ok)) if(any(!subject_covariates_ok)) stop("subject_covariates must be as long as the number of subjects")
  attr(dadm_list, "components") <- components
  attr(dadm_list, "shared_ll_idx") <- components
  return(list(ll_func = ll_func, pars = pars, prior = prior,
              dadm_list = dadm_list, subjects = subjects))
}

#' Runs IS2 from Tran et al. 2021 on a list of samplers
#'
#' Runs IS2 on a list of samplers, only works for types standard, factor and diagonal yet.
#'
#' @param samplers A list of samplers
#' @param filter A string. Indicates which stage to take samples from
#' @param subfilter An integer or vector. If integer specifies how many samples to remove from within that stage. If vector used as index for samples to keep.
#' @param IS_samples An integer. Specifies how many IS2 samples to collect
#' @param max_particles An integer. Specifies the maximum number of particles to collect before stopping one IS iteration.
#' @param stepsize_particles An integer. It will increase particles till optimal variance with this stepsize.
#' @param n_cores An integer. Specifies how many cores to run IS_2 on.
#' @param df An integer. The degrees of freedom used in the t-distribution used as IS distribution for the group-level proposals.
#'
#' @return Samplers, with IS2 attribute
#' @export
#'
#' @examples
run_IS2 <- function(samplers, filter = "sample", subfilter = 0, IS_samples = 1000,
                    stepsize_particles = 500, max_particles = 5000, n_cores = 1, df = 5, useC = TRUE){
  samples_merged <- merge_samples(samplers)
  IS_samples <- IS2(samples_merged, filter, subfilter = subfilter, IS_samples, stepsize_particles, max_particles, n_cores, df, useC)
  attr(samplers, "IS_samples") <- IS_samples
  return(samplers)
}
