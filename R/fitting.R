#### Fitting automation
#' Model estimation in EMC2
#'
#' General purpose function to estimate models specified in EMC2.
#'
#' @details
#'
#' ``stop_criteria`` is either a list of lists with names of the stages,
#' or a single list in which case its assumed to be for the sample `stage` (see examples).
#' The potential stop criteria to be set are:
#'
#' ``selection`` (character vector): For which parameters the ``stop_criteria`` should hold
#'
#' ``mean_gd`` (numeric): The mean Gelman-Rubin diagnostic across all parameters in the selection
#'
#' ``max_gd`` (numeric): The max Gelman-Rubin diagnostic across all parameters in the selection
#'
#' ``min_unique`` (integer): The minimum number of unique samples in the MCMC chains across all parameters in the selection
#'
#' ``min_es`` (integer): The minimum number of effective samples across all parameters in the selection
#'
#' ``omit_mpsrf`` (Boolean): Whether to include the multivariate point-scale reduction factor in the Gelman-Rubin diagnostic. Default is ``FALSE``.
#'
#' ``iter`` (integer): The number of MCMC samples to collect.
#'
#' The estimation is performed using particle-metropolis within-Gibbs sampling.
#' For sampling details see:
#'
#' Gunawan, D., Hawkins, G. E., Tran, M.-N., Kohn, R., & Brown, S. (2020).
#' New estimation approaches for the hierarchical linear ballistic accumulator model.
#' *Journal of Mathematical Psychology* ,96, 102368. doi.org/10.1016/j.jmp.2020.102368
#'
#' Stevenson, N., Donzallaz, M. C., Innes, R. J., Forstmann, B., Matzke, D., & Heathcote, A. (2024).
#' EMC2: An R Package for cognitive models of choice. doi.org/10.31234/osf.io/2e4dq
#'
#'
#' @param samplers An EMC2 samplers object created with ``make_samplers``,
#' or a path to where the samplers object is stored.
#' @param stage A string. Indicates which stage to start the run from, either ``preburn``, ``burn``, ``adapt`` or ``sample``.
#' If unspecified, it will run the subsequent stage (if there is one).
#' @param iter An integer. Indicates how many iterations to run in the sampling stage.
#' @param p_accept A double. The target acceptance probability of the MCMC process.
#' This fine-tunes the width of the search space to obtain the desired acceptance probability. Defaults to .8
#' @param step_size An integer. After each step, the stopping requirements as specified
#' by ``stop_criteria`` are checked and proposal distributions are updated. Defaults to 100.
#' @param verbose Logical. Whether to print messages between each step with the current status regarding the ``stop_criteria``.
#' @param verboseProgress Logical. Whether to print a progress bar within each step or not.
#' Will print one progress bar for each chain and only if ``cores_for_chains = 1``.
#' @param fileName A string. If specified, will auto-save samplers object at this location on every iteration.
#' @param particles An integer. How many particles to use, default is `NULL` and
#' ``particle_factor`` is used instead. If specified, ``particle_factor`` is overwritten.
#' @param particle_factor An integer. ``particle_factor`` multiplied by the square
#' root of the number of sampled parameters determines the number of particles used.
#' @param cores_per_chain An integer. How many cores to use per chain. Parallelizes across
#' participant calculations. Only available on Linux or Mac OS. For Windows, only
#' parallelization across chains (``cores_for_chains``) is available.
#' @param cores_for_chains An integer. How many cores to use across chains.
#' Defaults to the number of chains. The total number of cores used is equal to ``cores_per_chain`` * ``cores_for_chains``.
#' @param max_tries An integer. How many times should it try to meet the finish
#' conditions as specified by ``stop_criteria``? Defaults to 20. ``max_tries`` is
#' ignored if the required number of iterations has not been reached yet.
#' @param n_blocks An integer. Number of blocks. Will block the parameter chains such that they are
#' updated in blocks. This can be helpful in extremely tough models with a large number of parameters.
#' @param stop_criteria A list. Defines the stopping criteria and for which types
#' of parameters these should hold. See the details and examples section.
#'
#' @return A list of samplers
#' @examples \dontrun{
#' # First define a design
#' design_DDMaE <- make_design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#' # Then make the samplers, we've omitted a prior here for brevity so default priors will be used.
#' samplers <- make_samplers(forstmann, design)
#'
#' # With the samplers object we can start sampling by simply calling run_emc
#' samplers <- run_emc(samplers, fileName = "intermediate_save_location.RData")
#'
#' # For particularly hard models it pays off to increase the ``particle_factor``
#' # and, although to a lesser extent, lower ``p_accept``.
#' samplers <- run_emc(samplers, particle_factor = 100, p_accept = .6)
#'
#' # Example of how to use the stop_criteria:
#' samplers <- run_emc(samplers, stop_criteria = list(mean_gd = 1.1, max_gd = 1.5,
#'             selection = c('alpha', 'sigma2'), omit_mpsrf = TRUE, min_es = 1000).
#' # In this case the stop_criteria are set for the sample stage, which will be
#' # run until the mean_gd < 1.1, the max_gd < 1.5 (omitting the multivariate psrf)
#' # and the effective sample size > 1000,
#' # for both the individual-subject parameters ("alpha")
#' # and the group-level variance parameters.
#'
#' # For the unspecified stages in the ``stop_criteria`` the default values
#' # are assumed which are found in Stevenson et al. 2024 <doi.org/10.31234/osf.io/2e4dq>
#'
#' # Alternatively, you can also specify the stop_criteria for specific stages by creating a
#' # list of lists:
#' samplers <- run_emc(samplers, stop_criteria = list("burn" = list(mean_gd = 1.1, max_gd = 1.5,
#'             selection = c('alpha')), "adapt" = list(min_unique = 300)))
#'}
#' @export

run_emc <- function(samplers, stage = NULL, iter = 1000, stop_criteria = NULL,
                    p_accept = .8, step_size = 100, verbose = TRUE, verboseProgress = FALSE, fileName = NULL,
                    particles = NULL, particle_factor=50, cores_per_chain = 1,
                    cores_for_chains = length(samplers), max_tries = 20, n_blocks = 1){
  if(!is.null(stop_criteria) & length(stop_criteria) == 1){
    stop_criteria[["sample"]] <- stop_criteria
  }
  stages_names <- c("preburn", "burn", "adapt", "sample")
  if(is.null(stop_criteria)){
    stop_criteria <- vector("list", length = 4)
    names(stop_criteria) <- stages_names
  }
  for(stage_name in stages_names){
    if(!stage_name %in% names(stop_criteria)) stop_criteria[stage_name] <- list(NULL)
  }

  stop_criteria <- stop_criteria[stages_names]
  stop_criteria <- mapply(get_stop_criteria, stages_names, stop_criteria, MoreArgs = list(type = attr(samplers[[1]], "variant_funs")$type))
  if(is.null(stop_criteria[["sample"]]$iter)) stop_criteria[["sample"]]$iter <- iter
  names(stop_criteria) <- stages_names
  if (is.character(samplers)) {
    samplers <- fix_fileName(samplers)
    if(is.null(fileName)) fileName <- samplers
    samplers <- loadRData(samplers)
  }
  if(is.null(stage)){
    nstage <- colSums(chain_n(samplers))
    if (nstage["preburn"]==0) {
      stage <- "preburn"
    } else{
      if (nstage["sample"]>0){
        stage <- "sample"
      } else{
        stage <- names(nstage)[which(nstage==0)[1] - 1]
      }
    }
  }
  if(stage == "preburn"){
    samplers <- run_samplers(samplers, stage = "preburn", stop_criteria[['preburn']], cores_for_chains = cores_for_chains, p_accept = p_accept,
                             step_size = step_size,  verbose = verbose, verboseProgress = verboseProgress,
                             fileName = fileName,
                             particles = particles, particle_factor =  particle_factor,
                             cores_per_chain = cores_per_chain, max_tries = max_tries, n_blocks = n_blocks)
  }

  if(any(stage %in% c("preburn", "burn"))){
    samplers <-  run_samplers(samplers, stage = "burn", stop_criteria[['burn']], cores_for_chains = cores_for_chains, p_accept = p_accept,
                              step_size = step_size,  verbose = verbose, verboseProgress = verboseProgress,
                              fileName = fileName,
                              particles = particles, particle_factor =  particle_factor,
                              cores_per_chain = cores_per_chain, max_tries = max_tries, n_blocks = n_blocks)
  }
  if(any(stage %in% c("preburn", "burn", "adapt"))){
    samplers <-  run_samplers(samplers, stage = "adapt", stop_criteria[['adapt']],cores_for_chains = cores_for_chains, p_accept = p_accept,
                              step_size = step_size,  verbose = verbose, verboseProgress = verboseProgress,
                              fileName = fileName,
                              particles = particles, particle_factor =  particle_factor,
                              cores_per_chain = cores_per_chain, max_tries = max_tries, n_blocks = n_blocks)
  }
  if(any(stage %in% c("preburn", "burn", "adapt", "sample")) ){
    samplers <-  run_samplers(samplers, stage = "sample",stop_criteria[['sample']],cores_for_chains = cores_for_chains, p_accept = p_accept,
                              step_size = step_size,  verbose = verbose, verboseProgress = verboseProgress,
                              fileName = fileName,
                              particles = particles, particle_factor = particle_factor,
                              cores_per_chain = cores_per_chain, max_tries = max_tries, n_blocks = n_blocks)
  }
  return(samplers)
}

get_stop_criteria <- function(stage, stop_criteria, type){
  if(is.null(stop_criteria)){
    if(stage == "preburn"){
      stop_criteria$iter <- 150
    }
    if(stage == "burn"){
      stop_criteria$mean_gd <- 1.1
      stop_criteria$omit_mpsrf <- TRUE
      if(type != "single"){
        stop_criteria$selection <- c("alpha", "mu")
      } else{
        stop_criteria$selection <- c("alpha")
      }

    }
    if(stage == "adapt"){
      stop_criteria$min_unique <- 600
    }
    if(stage == "sample"){
      stop_criteria$max_gd <- 1.1
      stop_criteria$omit_mpsrf <- TRUE
      if(type != "single"){
        stop_criteria$selection <- c("alpha", "mu")
      } else{
        stop_criteria$selection <- c("alpha")
      }
    }
  }
  if(!is.null(stop_criteria$max_gd) || !is.null(stop_criteria$mean_gd)){
    if(is.null(stop_criteria$selection)) stop_criteria$selection <- c('alpha', 'mu')
  }
  if(stage == "adapt" & is.null(stop_criteria$min_unique)) stop_criteria$min_unique <- 600
  if(stage != "adapt" & !is.null(stop_criteria$min_unique)) stop("min_unique only applicable for adapt stage, try min_es instead.")
  return(stop_criteria)
}

#' Custom function for more controlled model estimation
#'
#' Although typically users will rely on ``run_emc``, this function can be used for more fine-tuned specification of estimation needs.
#' The function will throw an error if a stage is skipped,
#' the stages have to be run in order ("preburn", "burn", "adapt", "sample").
#' More details can be found in the ``run_emc`` help files (``?run_emc``).
#'
#' @param samplers A list of samplers, could be in any stage, as long as they've been initialized with make_samplers
#' @param stage A string. Indicates which stage is to be run, either `preburn`, `burn`, `adapt` or `sample`
#' @param p_accept A double. The target acceptance probability of the MCMC process.
#' This fine-tunes the width of the search space to obtain the desired acceptance probability. Defaults to .8
#' @param step_size An integer. After each step, the stopping requirements as
#' specified by `stop_criteria` are checked and proposal distributions are updated. Defaults to 100.
#' @param verbose Logical. Whether to print messages between each step with the current status regarding the stop_criteria.
#' @param verboseProgress Logical. Whether to print a progress bar within each step or not. Will print one progress bar for each chain and only if cores_for_chains = 1.
#' @param fileName A string. If specified will autosave samplers at this location on every iteration.
#' @param particles An integer. How many particles to use, default is `NULL` and ``particle_factor`` is used instead.
#' If specified will override ``particle_factor``.
#' @param particle_factor An integer. `particle_factor` multiplied by the square root of the number of sampled parameters determines the number of particles used.
#' @param cores_per_chain An integer. How many cores to use per chain.
#' Parallelizes across participant calculations. Only available on Linux or Mac OS.
#' For Windows, only parallelization across chains (``cores_for_chains``) is available.
#' @param cores_for_chains An integer. How many cores to use across chains.
#' Defaults to the number of chains. the total number of cores used is equal to ``cores_per_chain`` * ``cores_for_chains``.
#' @param max_tries An integer. How many times should it try to meet the finish
#' conditions as specified by stop_criteria? Defaults to 20. max_tries is ignored if the required number of iterations has not been reached yet.
#' @param n_blocks An integer. Number of blocks. Will block the parameter chains such that they are updated in blocks. This can be helpful in extremely tough models with a large number of parameters.
#' @param stop_criteria A list. Defines the stopping criteria and for which types of parameters these should hold. See ``?run_emc``.
#' @export
#' @return A list of samplers
#' @examples \dontrun{
#' # First define a design
#' design_DDMaE <- make_design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#' # Then make the samplers, we've omitted a prior here for brevity so default priors will be used.
#' samplers <- make_samplers(forstmann, design)
#'
#' # Now for example we can specify that we only want to run the "preburn" phase
#' # for MCMC 200 iterations
#' samplers <- run_samplers(samplers, stage = "preburn", stop_criteria = list(iter = 200))
#'}

run_samplers <- function(samplers, stage, stop_criteria,
                         p_accept = .8, step_size = 100, verbose = FALSE, verboseProgress = FALSE,
                         fileName = NULL,
                         particles = NULL, particle_factor=50, cores_per_chain = 1,
                         cores_for_chains = length(samplers), max_tries = 20, n_blocks = 1){
  if(Sys.info()[1] == "Windows" & cores_per_chain > 1) stop("only cores_for_chains can be set on Windows")
  if (verbose) message(paste0("Running ", stage, " stage"))
  attributes <- get_attributes(samplers)
  total_iters_stage <- chain_n(samplers)[,stage][1]
  if(stage != "preburn"){
    iter <- stop_criteria[["iter"]] + total_iters_stage
  } else{
    iter <- stop_criteria[["iter"]]
  }
  progress <- check_progress(samplers, stage, iter, stop_criteria, max_tries, step_size, cores_per_chain*cores_for_chains, verbose, n_blocks = n_blocks)
  samplers <- progress$samplers
  progress <- progress[!names(progress) == 'samplers'] # Frees up memory, courtesy of Steven
  particle_factor_in <- particle_factor; p_accept_in <- p_accept
  while(!progress$done){
    if(!is.null(progress$n_blocks)) n_blocks <- progress$n_blocks
    samplers <- add_proposals(samplers, stage, cores_per_chain*cores_for_chains, n_blocks)
    if(!is.null(progress$gds_bad)){
      particle_factor_in <- particle_factor_in + .05 * progress$gds_bad * particle_factor_in
      p_accept_in <- pmax(0.4, p_accept - progress$gds_bad*.3)
      particle_factor_in[!progress$gds_bad] <- particle_factor
    }
    samplers <- auto_mclapply(samplers,run_stages, stage = stage, iter= progress$step_size,
                              verbose=verbose,  verboseProgress = verboseProgress,
                              particles=particles,particle_factor=particle_factor_in,
                              p_accept=p_accept_in, n_cores=cores_per_chain, mc.cores = cores_for_chains)
    for(i in 2:length(samplers)){ # Frees up memory, courtesy of Steven
      samplers[[i]]$data <- samplers[[1]]$data
    }
    progress <- check_progress(samplers, stage, iter, stop_criteria, max_tries, step_size, cores_per_chain*cores_for_chains,verbose, progress,n_blocks)
    samplers <- progress$samplers
    progress <- progress[!names(progress) == 'samplers'] # Frees up memory, courtesy of Steven
    if(!is.null(fileName)){
      fileName <- fix_fileName(fileName)
      attr(samplers,"data_list") <- attributes$data_list
      attr(samplers,"design_list") <- attributes$design_list
      attr(samplers,"model_list") <- attributes$model_list
      class(samplers) <- "emc"
      save(samplers, file = fileName)
    }
  }
  samplers <- get_attributes(samplers, attributes)
  class(samplers) <- "emc"
  return(samplers)
}

run_stages <- function(sampler, stage = "preburn", iter=0, verbose = TRUE, verboseProgress = TRUE,
                       particles=NULL,particle_factor=50, p_accept= NULL, n_cores=1)
{

  if (is.null(particles)){
    # if(!is.null(sampler$g_map_fixed)){
    #   particles <- round(particle_factor*sqrt(length(sampler$pars_random)/sampler$n_subjects))
    # } else{
    #   particles <- round(particle_factor*sqrt(length(sampler$par_names)))
    # }
    max_pars <- max(table(attr(sampler[[1]], "components")))
    particles <- round(particle_factor*sqrt(max_pars))

  }
  if (!sampler$init) {
    # if(!is.null(sampler$g_map_fixed)){
    #   sampler <- init_lm(sampler, n_cores = n_cores)
    # } else{
    #   sampler <- init(sampler, n_cores = n_cores)
    # }
    sampler <- init(sampler, n_cores = n_cores)
  }
  if (iter == 0) return(sampler)
  # if(!is.null(sampler$g_map_fixed)){
  #   sampler <- run_stage_lm(sampler, stage = stage,iter = iter, particles = particles,
  #                           n_cores = n_cores, p_accept = p_accept, verbose = verbose, verboseProgress = verboseProgress)
  # } else{
  #     sampler <- run_stage(sampler, stage = stage,iter = iter, particles = particles,
  #                          n_cores = n_cores, p_accept = p_accept, verbose = verbose, verboseProgress = verboseProgress)
  # }
  sampler <- run_stage(sampler, stage = stage,iter = iter, particles = particles,
                       n_cores = n_cores, p_accept = p_accept, verbose = verbose, verboseProgress = verboseProgress)
  return(sampler)
}

add_proposals <- function(samplers, stage, n_cores, n_blocks){
  if(stage != "preburn"){
    # if(!is.null(samplers[[1]]$g_map_fixed)){
    #   samplers <- create_cov_proposals_lm(samplers)
    # } else{    }
    samplers <- create_cov_proposals(samplers, do_block = stage != "sample")
    if(!is.null(n_blocks)){
      if(n_blocks > 1){
        components <- sub_blocking(samplers, n_blocks)
        for(i in 1:length(samplers)){
          attr(samplers[[i]]$data, "components") <- components
        }
      }
    }
  }
  if(stage == "sample"){
    # if(!is.null(samplers[[1]]$g_map_fixed)){
    #   samplers <- create_eff_proposals_lm(samplers, n_cores)
    # } else{    }
    samplers <- create_eff_proposals(samplers, n_cores)
  }
  return(samplers)
}

check_progress <- function (samplers, stage, iter, stop_criteria,
                            max_tries, step_size, n_cores, verbose, progress = NULL,
                            n_blocks)
{
  min_es <- stop_criteria$min_es
  if(is.null(min_es)) min_es <- 0
  selection <- stop_criteria$selection
  min_unique <- stop_criteria$min_unique
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
  gd <- check_gd(samplers, stage, stop_criteria[["max_gd"]], stop_criteria[["mean_gd"]], trys, verbose,
                 iter = total_iters_stage, selection, omit_mpsrf = stop_criteria[["omit_mpsrf"]],
                 n_blocks)
  iter_done <- ifelse(is.null(iter) || length(iter) == 0, TRUE, total_iters_stage >= iter)
  if (min_es == 0) {
    es_done <- TRUE
  }
  else if (iters_total != 0) {
    curr_min_es <- Inf
    for(select in selection){
      curr_min_es <- min(c(es_summary_new(samplers, selection = select,
                                                filter = stage, stat_only = TRUE), curr_min_es))
    }
    if (verbose)
      message("Smallest effective size = ", round(curr_min_es))
    es_done <- ifelse(!samplers[[1]]$init, FALSE, curr_min_es >
                        min_es)
  }
  else {
    es_done <- FALSE
  }
  trys_done <- ifelse(is.null(max_tries), FALSE, trys >= max_tries)
  if (stage == "adapt") {
    samples_merged <- merge_samples(samplers)
    test_samples <- extract_samples(samples_merged, stage = "adapt",
                                    samples_merged$samples$idx, n_chains = length(samplers))
    # if(!is.null(samplers[[1]]$g_map_fixed)){
    #   adapted <- test_adapted_lm(samplers[[1]], test_samples, min_unique, n_cores, verbose)
    # } else{    }
    adapted <- test_adapted(samplers[[1]], test_samples,
                            min_unique, n_cores, verbose)

  }
  else {
    adapted <- TRUE
  }
  done <- (es_done & iter_done & gd$gd_done & adapted) | (trys_done & iter_done)
  if(es_done & gd$gd_done & adapted & !iter_done){
    step_size <- min(step_size, abs(iter - total_iters_stage))[1]
  }
  if (trys_done & iter_done) {
    if(!(es_done & gd$gd_done & adapted)){
      warning("Max tries reached. If this happens in burn-in while trying to get
            gelman diagnostics small enough, you might have a particularly hard model.
            Make sure your model is well specified. If so, you can run adapt and
            sample, if run for long enough, sample usually converges eventually.")
    }
  }
  return(list(samplers = gd$samplers, done = done, step_size = step_size,
              trys = trys, iters_total = iters_total, n_blocks = gd$n_blocks,
              gds_bad = gd$gds_bad))
}

check_gd <- function(samplers, stage, max_gd, mean_gd, omit_mpsrf, trys, verbose,
                     selection, iter, n_blocks = 1)
{
  get_gds <- function(samplers,omit_mpsrf, selection, stage) {
    gd_out <- c()
    for(select in selection){
      gd <- unlist(gd_summary_new(samplers, selection = select, filter = stage,
                                  omit_mpsrf = omit_mpsrf, stat = NULL))
      gd_out <- c(gd_out, c(gd))
    }
    return(gd_out)
  }

  if(is.null(max_gd) & is.null(mean_gd)) return(list(gd_done = TRUE, samplers = samplers))
  if(!samplers[[1]]$init | !stage %in% samplers[[1]]$samples$stage)
    return(list(gd_done = FALSE, samplers = samplers))
  if(is.null(omit_mpsrf)) omit_mpsrf <- TRUE
  gd <- get_gds(samplers,omit_mpsrf,selection, stage)
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

  ok_gd <- ok_mean_gd & ok_max_gd
  if(!ok_gd) {
    n_remove <- round(chain_n(samplers)[,stage][1]/3)
    samplers_short <- try(lapply(samplers,remove_iterations,subfilter=n_remove,filter=stage),silent=TRUE)
    if (is(samplers_short,"try-error")){
      gd_short <- Inf
    } else{
      gd_short <- get_gds(samplers_short,omit_mpsrf,selection, stage)

    }
    if (is.null(max_gd) & (mean(gd_short) < mean(gd)) | (!is.null(max_gd) & (max(gd_short) < max(gd)))) {
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
    ok_gd <- ok_mean_gd & ok_max_gd
  }


  n_blocks_old <- n_blocks
  # if(iter > 1000 & stage == "sample" & !ok_gd) {
  #   n_blocks <- floor(iter/1000) + 1
  #   n_blocks <- max(n_blocks_old, n_blocks)
  # }
  if(stage == "sample" & !ok_gd & "alpha" %in% selection) {
    gds <- gd_pmwg(samplers, selection = 'alpha', return_summary = F, mapped = F, print_summary = F, filter = stage, omit_mpsrf = omit_mpsrf)
    if(!is.null(mean_gd)){
      gds_bad <- (rowMeans(gds) > mean_gd)
    } else{
      gds_bad <- (apply(gds, 1, max) > max_gd)
    }
  } else{
    gds_bad <- NULL
  }
  if(verbose) {
    if(n_blocks_old != n_blocks) {
      message(paste0("More than ", floor(iter/1000)*1000, " sample iterations past without convergence. Now block updating parameters with ", n_blocks, " blocks."))
    }
  }
  if(verbose) {
    if (omit_mpsrf) type <- "psrf" else type <- "m/psrf"
    if (!is.null(mean_gd)) message("Mean ",type," = ",round(mean(gd),3)) else
      if (!is.null(max_gd)) message("Max ",type," = ",round(max(gd),3))
  }
  return(list(gd = gd, gd_done = ok_gd, samplers = samplers, n_blocks = n_blocks, gds_bad = gds_bad))
}


create_eff_proposals <- function(samplers, n_cores){
  samples_merged <- merge_samples(samplers)
  test_samples <- extract_samples(samples_merged, stage = c("adapt", "sample"), max_n_sample = 750, n_chains = length(samplers))
  variant_funs <- attr(samplers[[1]], "variant_funs")
  components <- attr(samplers[[1]]$data, "components")[!samplers[[1]]$grouped]
  for(i in 1:length(samplers)){
    iteration = round(test_samples$iteration * i/length(samplers))
    n_pars <- sum(!samplers[[1]]$grouped)
    nuisance <- samplers[[1]]$nuisance[!samplers[[1]]$grouped]
    n_subjects <- samplers[[1]]$n_subjects
    eff_mu <- matrix(0, nrow = n_pars, ncol = n_subjects)
    eff_var <- array(0, dim = c(n_pars, n_pars, n_subjects))
    for(comp in unique(components)){
      idx <- comp == components
      nuis_idx <- nuisance[idx]
      if(any(nuis_idx)){
        type <- samples_merged$sampler_nuis$type
        conditionals <- auto_mclapply(X = 1:n_subjects,
                                      FUN = variant_funs$get_conditionals,samples = test_samples,
                                      n_pars = sum(idx[!nuisance]), iteration =  iteration, idx = idx[!nuisance],
                                      mc.cores = n_cores)
        conditionals_nuis <- auto_mclapply(X = 1:n_subjects,
                                           FUN = get_variant_funs(type)$get_conditionals,samples = test_samples$nuisance,
                                           n_pars = sum(idx[nuisance]), iteration =  iteration, idx = idx[nuisance],
                                           mc.cores = n_cores)
        conditionals <- array(unlist(conditionals), dim = c(sum(idx[!nuisance]), sum(idx[!nuisance]) + 1, n_subjects))
        conditionals_nuis <- array(unlist(conditionals_nuis), dim = c(sum(idx[nuisance]), sum(idx[nuisance]) + 1, n_subjects))
        eff_mu[idx & !nuisance,] <- conditionals[,1,]
        eff_var[idx & !nuisance, idx & !nuisance,] <- conditionals[,2:(sum(idx[!nuisance])+1),]
        eff_mu[idx & nuisance,] <- conditionals_nuis[,1,]
        eff_var[idx & nuisance,idx & nuisance,] <- conditionals_nuis[,2:(sum(idx[nuisance])+1),]
      } else{
        conditionals <- auto_mclapply(X = 1:n_subjects,
                                      FUN = variant_funs$get_conditionals,samples = test_samples,
                                      n_pars = sum(idx[!nuisance]), iteration =  iteration, idx = idx[!nuisance],
                                      mc.cores = n_cores)
        conditionals <- array(unlist(conditionals), dim = c(sum(idx[!nuisance]), sum(idx[!nuisance]) + 1, n_subjects))
        eff_mu[idx & !nuisance,] <- conditionals[,1,]
        eff_var[idx & !nuisance,idx & !nuisance,] <- conditionals[,2:(sum(idx[!nuisance])+1),]
      }

    }
    # eff_mu <- lapply(conditionals, FUN = function(x) x$eff_mu)
    # eff_var <- lapply(conditionals, FUN = function(x) x$eff_var)
    # eff_alpha <- lapply(conditionals, FUN = function(x) x$eff_alpha)
    # eff_tau <- lapply(conditionals, FUN = function(x) x$eff_tau)

    eff_mu <- split(eff_mu, col(eff_mu))
    eff_var <- apply(eff_var, 3, identity, simplify = F)
    attr(samplers[[i]], "eff_mu") <- eff_mu
    attr(samplers[[i]], "eff_var") <- eff_var
    # attr(samplers[[i]], "eff_alpha") <- eff_alpha
    # attr(samplers[[i]], "eff_tau") <- eff_tau
  }
  return(samplers)
}


sub_blocking <- function(samplers, n_blocks){
  covs <- lapply(samplers, FUN = function(x){return(attr(x, "chains_cov"))})
  out <- array(0, dim = dim(covs[[1]][[1]]))
  for(i in 1:length(covs)){
    cov_tmp <- covs[[1]]
    for(j in 1:length(cov_tmp)){
      out <- out + cov2cor(cov_tmp[[1]])
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

create_cov_proposals <- function(samplers, samples_idx = NULL, do_block = TRUE){
  get_covs <- function(sampler, samples_idx, sub){
    return(var(t(sampler$samples$alpha[,sub, samples_idx])))
  }
  n_subjects <- samplers[[1]]$n_subjects
  n_chains <- length(samplers)
  grouped <- samplers[[1]]$grouped
  n_pars <- sum(!grouped)

  if(is.null(samples_idx)){
    idx_subtract <- min(250, samplers[[1]]$samples$idx/2)
    samples_idx <- round(samplers[[1]]$samples$idx - idx_subtract):samplers[[1]]$samples$idx
  }
  components <- attr(samplers[[1]]$data, "components")
  block_idx <- block_variance_idx(components[!grouped])
  for(j in 1:n_chains){
    chains_cov <- array(NA_real_, dim = c(n_pars, n_pars, n_subjects))
    for(sub in 1:n_subjects){
      mean_covs <- get_covs(samplers[[j]], samples_idx, sub)
      if(do_block) mean_covs[block_idx] <- 0
      if(is.negative.semi.definite(mean_covs)){
        if(!is.null(attr(samplers[[j]], "chains_cov")[[sub]])){
          chains_cov[,,sub] <- attr(samplers[[j]], "chains_cov")[[sub]]
        } else{
          chains_cov[,,sub] <- diag(nrow(mean_covs))
        }
      } else{
        chains_cov[,,sub] <-  as.matrix(nearPD(mean_covs)$mat)
      }
    }
    if(any(grouped)){
      if(sum(grouped) == 1){
        chains_cov_grouped <- as.matrix(var(samplers[[j]]$samples$grouped_pars[, samples_idx]))
      } else{
        chains_cov_grouped <- var(t(samplers[[j]]$samples$grouped_pars[, samples_idx]))
        if(is.negative.semi.definite(chains_cov_grouped)){
          chains_cov_grouped <- attr(samplers[[j]], "chains_cov_grouped")
        } else{
          chains_cov_grouped <-  as.matrix(nearPD(chains_cov_grouped)$mat)
        }
      }

      attr(samplers[[j]], "chains_cov_grouped") <- chains_cov_grouped
    }
    chains_cov <- apply(chains_cov, 3, identity, simplify = F)
    attr(samplers[[j]], "chains_cov") <- chains_cov
  }
  return(samplers)
}

# create_eff_proposals_lm <- function(samplers, n_cores){
#   samples_merged <- merge_samples(samplers)
#   test_samples <- extract_samples(samples_merged, stage = c("adapt", "sample"), max_n_sample = 1000)
#   for(i in 1:length(samplers)){
#     iteration = round(test_samples$iteration * i/length(samplers))
#     attr(samplers[[i]], "eff_props") <- auto_mclapply(X = samplers[[1]]$subjects,FUN = get_conditionals_lm,samples = test_samples, iteration = iteration,
#                                                       mc.cores = n_cores)
#     attr(samplers[[i]], "eff_props_between") <- get_conditionals_lm_between(samples = test_samples, iteration = iteration)
#   }
#   return(samplers)
# }

# create_cov_proposals_lm <- function(samplers, samples_idx = NULL){
#   get_covs <- function(sampler, samples_idx, par_idx){
#     return(var(t(sampler$samples$random[par_idx, samples_idx])))
#   }
#   n_subjects <- samplers[[1]]$n_subjects
#   n_chains <- length(samplers)
#
#   if(is.null(samples_idx)){
#     idx_subtract <- min(250, samplers[[1]]$samples$idx/2)
#     samples_idx <- round(samplers[[1]]$samples$idx - idx_subtract):samplers[[1]]$samples$idx
#   }
#   pars_per_subject <- lapply(samplers[[1]]$subjects, FUN = function(x, pars) pars[get_sub_idx(x, pars)], samplers[[1]]$pars_random)
#   between_idx <- !grepl("subjects", samplers[[1]]$pars_random)
#
#   chains_cov_sub <- vector("list", length = n_subjects)
#   for(j in 1:n_chains){
#     for(sub in 1:n_subjects){
#       mean_covs <- get_covs(samplers[[j]], samples_idx, par_idx = pars_per_subject[[sub]])
#       if(is.negative.semi.definite(mean_covs)){
#         chains_cov_sub[[sub]] <- attr(samplers[[j]], "chains_cov_sub")[[sub]]
#       } else{
#         chains_cov_sub[[sub]] <-  as.matrix(nearPD(mean_covs)$mat)
#       }
#     }
#     mean_covs_between <- var(t(rbind(samplers[[j]]$samples$fixed[, samples_idx],
#                                      samplers[[j]]$samples$random[between_idx,samples_idx])))
#     if(is.negative.semi.definite(mean_covs_between)){
#       mean_covs_between <- attr(samplers[[j]], "chains_cov_between")
#     } else{
#       mean_covs_between <-  as.matrix(nearPD(mean_covs_between)$mat)
#     }
#     attr(samplers[[j]], "chains_cov_sub") <- chains_cov_sub
#     attr(samplers[[j]], "chains_cov_between") <- mean_covs_between
#
#   }
#   return(samplers)
# }

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

# test_adapted_lm <- function(samples, test_samples, min_unique, n_cores_conditional = 1,
#                             verbose = FALSE){
#   unq_per_sub <- sapply(samples$subjects, FUN = function(x, samples){
#     pars <- samples[get_sub_idx(x, rownames(samples)),]
#     return(length(unique(pars[1,])))}, test_samples$random)
#   if(all(unq_per_sub > min_unique)){
#     if(verbose){
#       message("Enough unique values detected: ", min_unique)
#       message("Testing proposal distribution creation")
#     }
#     attempt <- tryCatch({
#       auto_mclapply(X = samples$subjects,FUN = get_conditionals_lm,samples = test_samples,
#                     mc.cores = n_cores_conditional)
#       get_conditionals_lm_between(samples = test_samples)
#     },error=function(e) e, warning=function(w) w)
#     if (any(class(attempt) %in% c("warning", "error", "try-error"))) {
#       if(verbose){
#         message("Can't create efficient distribution yet")
#         message("Increasing required unique values and continuing adaptation")
#       }
#       return(FALSE)
#     }
#     else {
#       if(verbose) message("Successfully adapted after ", test_samples$iteration, "iterations - stopping adaptation")
#       return(TRUE)
#     }
#
#   } else{
#     return(FALSE)
#   }
# }

test_adapted <- function(sampler, test_samples, min_unique, n_cores_conditional = 1,
                         verbose = FALSE)
{
  # Function used by run_adapt to check whether we can create the conditional.

  # Only need to check uniqueness for one parameter
  first_par <- as.matrix(test_samples$alpha[1, , ])
  if(ncol(first_par) == 1) first_par <- t(first_par)
  # Split the matrix into a list of vectors by subject
  # all subjects is greater than unq_vals
  n_unique_sub <- apply(first_par, 1, FUN = function(x) return(length(unique(x))))
  n_pars <- sampler$n_pars
  variant_funs <- attr(sampler, "variant_funs")
  grouped <- sampler$grouped
  components <- attr(sampler$data, "components")[!grouped]
  nuisance <- sampler$nuisance[!grouped]
  if (length(n_unique_sub) != 0 & all(n_unique_sub > min_unique)) {
    if(verbose){
      message("Enough unique values detected: ", min_unique)
      message("Testing proposal distribution creation")
    }
    attempt <- tryCatch({
      for(comp in unique(components)){
        idx <- comp == components
        nuis_idx <- nuisance[idx]
        if(any(nuis_idx)){
          type <- sampler$sampler_nuis$type
          auto_mclapply(X = 1:sampler$n_subjects,
                        FUN = get_variant_funs(type)$get_conditionals,samples = test_samples$nuisance,
                        n_pars = sum(idx[nuisance]), idx = idx[nuisance],
                        mc.cores = n_cores_conditional)
        }
        auto_mclapply(X = 1:sampler$n_subjects,FUN = variant_funs$get_conditionals,samples = test_samples,
                      n_pars = sum(idx[!nuisance]), idx = idx[!nuisance], mc.cores = n_cores_conditional)
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
      if(verbose) message("Successfully adapted - stopping adaptation")
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
#' @param max_tries An integer. How many times it will try to meet the finish conditions. Default is 50.
#' @param n_blocks An integer. Will block the parameter chains such that they are updated in blocks. This can be helpful in extremely tough models with large number of parameters.
#' @param stop_criteria A list. Defines the stopping criteria and for which types of parameters these should hold. See run_emc
#'
#' @return A list of samplers

auto_burn <- function(samplers, preburn = 150,
                      p_accept = .8, step_size = 100, verbose = FALSE, verboseProgress = FALSE,
                      fileName = NULL, stop_criteria = NULL,
                      particles = NULL, particle_factor=50, cores_per_chain = 1,
                      cores_for_chains = length(samplers), max_tries = 20, n_blocks = 1){
  if(!is.null(stop_criteria) & length(stop_criteria) == 1){
    stop_criteria[["burn"]] <- stop_criteria
  }
  if(is.null(stop_criteria)){
    stop_criteria <- vector("list", length = 2)
    names(stop_criteria) <- c("preburn", "burn")
  }
  for(stage_name in c("preburn", "burn")){
    if(!stage_name %in% names(stop_criteria)) stop_criteria[stage_name] <- list(NULL)
  }
  stop_criteria <- stop_criteria[c("preburn", "burn")]
  stop_criteria <- mapply(get_stop_criteria, c("preburn", "burn"), stop_criteria, MoreArgs = list(type = attr(samplers[[1]], "variant_funs")$type))
  stop_criteria[["preburn"]]$iter <- preburn
  names(stop_criteria) <- c("preburn", "burn")

  samplers <- run_samplers(samplers, stage = "preburn", stop_criteria = stop_criteria[["preburn"]],
                           cores_for_chains = cores_for_chains, p_accept = p_accept,
                           step_size = step_size,  verbose = verbose, verboseProgress = verboseProgress,
                           fileName = fileName,
                           particles = particles, particle_factor =  particle_factor,
                           cores_per_chain = cores_per_chain, max_tries = max_tries, n_blocks = n_blocks)
  samplers <-  run_samplers(samplers, stage = "burn",  stop_criteria = stop_criteria[["preburn"]], cores_for_chains = cores_for_chains, p_accept = p_accept,
                            step_size = step_size,  verbose = verbose, verboseProgress = verboseProgress,
                            fileName = fileName,
                            particles = particles, particle_factor =  particle_factor,
                            cores_per_chain = cores_per_chain, max_tries = max_tries, n_blocks = n_blocks)
  return(samplers)
}
#' Runs adapt stage for samplers.
#'
#' Special instance of `run_samplers`, with default arguments specified for completing adaptation.
#'
#' @param samplers A list of samplers, could be in any stage, as long as they've been initialized with make_samplers
#' @param p_accept A double. The target acceptance probability of the MCMC process. This will fine tune the width of the search space. Default = .8
#' @param step_size An integer. After each of these steps, the requirements will be checked if they are met and proposal distributions will be updated. Default = 100.
#' @param verbose Logical. Whether to print emc related messages
#' @param verboseProgress Logical. Whether to print sampling related messages
#' @param fileName A string. If specified will autosave samplers at this location.
#' @param particles An integer. How many particles to use, default is NULL and particle_factor is used. If specified will override particle_factor
#' @param particle_factor An integer. Particle factor multiplied by the square root of the number of sampled parameters will determine the number of particles used.
#' @param cores_per_chain An integer. How many cores to use per chain. Parallelizes across participant calculations.
#' @param cores_for_chains An integer. How many cores to use across chains. Default is the number of chains.
#' @param max_tries An integer. How many times it will try to meet the finish conditions. Default is 20.
#' @param n_blocks An integer. Will block the parameter chains such that they are updated in blocks. This can be helpful in extremely tough models with large number of parameters.
#' @param stop_criteria A list. Defines the stopping criteria and for which types of parameters these should hold. See run_emc.
#'
#' @return A list of samplers.

run_adapt <- function(samplers, stop_criteria = NULL,
                      p_accept = .8, step_size = 100, verbose = FALSE, verboseProgress = FALSE,
                      fileName = NULL,
                      particles = NULL, particle_factor=50, cores_per_chain = 1,
                      cores_for_chains = length(samplers), max_tries = 20, n_blocks = 1)
{
  if(is.null(stop_criteria)){
    stop_criteria <- list()
    stop_criteria[["adapt"]] <-list(NULL)
  }
  stop_criteria <- get_stop_criteria("adapt", stop_criteria[["adapt"]], type = attr(samplers[[1]], "variant_funs")$type)

  samplers <- run_samplers(samplers, stage = "adapt",  stop_criteria = stop_criteria[[1]],
                           cores_for_chains = cores_for_chains, p_accept = p_accept,
                           step_size = step_size,  verbose = verbose, verboseProgress = verboseProgress,
                           fileName = fileName,
                           particles = particles, particle_factor =  particle_factor,
                           cores_per_chain = cores_per_chain, max_tries = max_tries, n_blocks = n_blocks)
  return(samplers)
}
#' Runs sample stage for samplers.
#'
#' Special instance of `run_samplers`, with default arguments specified for running sample stage.
#'
#' @param samplers A list of samplers, could be in any stage, as long as they've been initialized with make_samplers
#' @param iter An integer. Indicates how many iterations to run,
#' @param p_accept A double. The target acceptance probability of the MCMC process. This will fine tune the width of the search space. Default = .8
#' @param step_size An integer. After each of these steps, the requirements will be checked if they are met and proposal distributions will be updated. Default = 100.
#' @param verbose Logical. Whether to print emc related messages
#' @param verboseProgress Logical. Whether to print sampling related messages
#' @param fileName A string. If specified will autosave samplers at this location.
#' @param particles An integer. How many particles to use, default is NULL and particle_factor is used. If specified will override particle_factor
#' @param particle_factor An integer. Particle factor multiplied by the square root of the number of sampled parameters will determine the number of particles used.
#' @param cores_per_chain An integer. How many cores to use per chain. Parallelizes across participant calculations.
#' @param cores_for_chains An integer. How many cores to use across chains. Default is the number of chains.
#' @param n_blocks An integer. Will block the parameter chains such that they are updated in blocks. This can be helpful in extremely tough models with large number of parameters.
#' @param max_tries An integer. How many times it will try to meet the finish conditions. Default is 20.
#' @param stop_criteria A list. Defines the stopping criteria and for which types of parameters these should hold. See run_emc.
#' @return A list of samplers

run_sample <- function(samplers, iter = 1000, stop_criteria = NULL,
                       p_accept = .8, step_size = 100, verbose = FALSE, verboseProgress = FALSE,
                       fileName = NULL,
                       particles = NULL, particle_factor=50, cores_per_chain = 1,
                       cores_for_chains = length(samplers), max_tries = 20, n_blocks = 1)
{
  if(is.null(stop_criteria)){
    stop_criteria <- list()
  }
  stop_criteria <- get_stop_criteria("sample", stop_criteria, type = attr(samplers[[1]], "variant_funs")$type)
  stop_criteria[["sample"]] <-iter
  samplers <- run_samplers(samplers, stage = "sample", stop_criteria[[1]], cores_for_chains = cores_for_chains, p_accept = p_accept,
                           step_size = step_size,  verbose = verbose, verboseProgress = verboseProgress,
                           fileName = fileName,
                           particles = particles, particle_factor =  particle_factor,
                           cores_per_chain = cores_per_chain, max_tries = max_tries, n_blocks = n_blocks)
  return(samplers)
}

#' `make_samplers()`
#'
#' initializes the sampling process by combining the data, prior,
#' and model specification into a `samplers` object that is needed in `run_emc()`.
#'
#' @param data A data frame, or a list of data frames. Needs to have the variable `subjects` as participant identifier.
#' @param design A list with a pre-specified design, the output of `make_design()`.
#' @param model A model list. If none is supplied, the model specified in `make_design()` is used.
#' @param type A string indicating whether to run a `standard` group-level, `blocked`, `diagonal`, `factor`, or `single` (i.e., non-hierarchical) model.
#' @param n_chains An integer. Specifies the number of mcmc chains to be run (has to be more than 1 to compute `rhat`).
#' @param compress A Boolean, if `TRUE` (i.e., the default), the data is compressed to speed up likelihood calculations.
#' @param rt_resolution A double. Used for compression, response times will be binned based on this resolution.
#' @param par_groups A vector. Only to be specified with type `blocked`, e.g., `c(1,1,1,2,2)` means the covariances
#' of the first three and of the last two parameters are estimated as two separate blocks.
#' @param n_factors An integer. Only to be specified with type `factor`.
#' @param constraintMat A matrix of rows equal to the number of estimated parameters, and columns equal to the number of factors. Only to be specified with type factor.
#' If `NULL`, default settings as specified in Innes et al. 2022 will be used.
#' @param prior_list A named list containing the prior. Default prior created if `NULL`. For the default priors, see `?get_prior_{type}`.
#' @param grouped_pars An integer vector. Parameters on this location of the vector of parameters are treated as constant across subjects
#' @param ... Additional, optional arguments.
#' @return A list
#' @examples dat <- forstmann
#'
#' # function that takes the lR factor (named diff in the following function) and
#' # returns a logical defining the correct response for each stimulus. In this
#' # case the match is simply such that the S factor equals the latent response factor.
#' matchfun <- function(d)d$S==d$lR
#'
#' # design an "average and difference" contrast matrix
#' ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"diff"))
#'
#' # specify design
#' design_LBABE <- make_design(data = dat,model=LBA,matchfun=matchfun,
#' formula=list(v~lM,sv~lM,B~E+lR,A~1,t0~1),
#' contrasts=list(v=list(lM=ADmat)),constants=c(sv=log(1)))
#'
#' # specify priors
#' pmean <- c(v=1,v_lMdiff=1,sv_lMTRUE=log(.5), B=log(.5),B_Eneutral=log(1.5),
#'            B_Eaccuracy=log(2),B_lRright=0, A=log(0.25),t0=log(.2))
#' psd <- c(v=1,v_lMdiff=0.5,sv_lMTRUE=.5,
#'          B=0.3,B_Eneutral=0.3,B_Eaccuracy=0.3,B_lRright=0.3,A=0.4,t0=.5)
#' prior_LBABE <- make_prior(design_LBABE, type = 'standard',pmean=pmean,psd=psd)
#'
#' # create samplers object
#' LBABE <- make_samplers(dat,design_LBABE,type="standard",  prior=prior_LBABE)
#'
#' @export

make_samplers <- function(data,design,model=NULL,
                          type="standard",
                          n_chains=3,compress=TRUE,rt_resolution=0.02,
                          prior_list = NULL,
                          grouped_pars = NULL,
                          par_groups=NULL,
                          n_factors=NULL,constraintMat = NULL, ...){

  # arguments for future compatibility
  formula <- NULL
  Lambda_mat <- NULL
  B_mat <- NULL
  K_mat <- NULL
  G_mat <- NULL
  xy <- NULL
  xeta <- NULL
  nuisance <- NULL
  nuisance_non_hyper <- NULL
  # overwrite those that were supplied
  optionals <- list(...)
  for (name in names(optionals) ) {
    assign(name, optionals[[name]])
  }

  if (!(type %in% c("standard","diagonal","blocked","factor","single", "lm", "infnt_factor", "SEM")))
    stop("type must be one of: standard,diagonal,blocked,factor,infnt_factor, lm, single")

  if(!is.null(nuisance) & !is.null(nuisance_non_hyper)){
    stop("You can only specify nuisance OR nuisance_non_hyper")
  }
  if (is(data, "data.frame")) data <- list(data)
  if(!is.null(prior_list) & !is.null(prior_list$theta_mu_mean)){
    prior_list <- list(prior_list)
  }
  # Sort subject together and add unique trial within subject integer
  # create overarching data list with one list per subject
  if(type == "lm"){
    if(length(data) > 1) stop("no joint models for lm yet")
    vars <- c()
    if(!is.null(formula)){
      for(form in formula){
        vars <- c(vars, split_form(form)$dep)
      }
      tmp <- data[[1]]
      aggr_data <- tmp[cumsum(table(tmp$subjects)),c("subjects", unique(vars))]
      for(i in 1:ncol(aggr_data)){
        if(colnames(aggr_data)[i] != "subjects" & is.factor(aggr_data[,i])){
          aggr_data[,i] <- factor(aggr_data[,i], levels = unique(aggr_data[,i]))
        }
      }
    } else{
      aggr_data <- NULL
    }

  }
  data <- lapply(data,function(d){
    if (!is.factor(d$subjects)) d$subjects <- factor(d$subjects)
    d <- d[order(d$subjects),]
    LC <- attr(d,"LC")
    UC <- attr(d,"UC")
    LT <- attr(d,"LT")
    UT <- attr(d,"UT")
    d <- add_trials(d)
    attr(d,"LC") <- LC
    attr(d,"UC") <- UC
    attr(d,"LT") <- LT
    attr(d,"UT") <- UT
    d
  })
  if (!is.null(names(design)[1]) && names(design)[1]=="Flist")
    design <- list(design)
  if (length(design)!=length(data))
    design <- rep(design,length(data))
  if (is.null(model)) model <- lapply(design,function(x){x$model})
  if (any(unlist(lapply(model,is.null))))
    stop("Must supply model if model is not in all design components")
  if (!is.null(names(model)[1]) && names(model)[1]=="type")
    model <- list(model)
  if (length(model)!=length(data))
    model <- rep(model,length(data))

  dadm_list <- vector(mode="list",length=length(data))
  rt_resolution <- rep(rt_resolution,length.out=length(data))
  for (i in 1:length(dadm_list)) {
    message("Processing data set ",i)
    # if (!is.null(design[[i]]$Ffunctions)) {
    #   pars <- attr(data[[i]],"pars")
    #   data[[i]] <- cbind.data.frame(data[[i]],data.frame(lapply(
    #     design[[i]]$Ffunctions,function(f){f(data[[i]])})))
    #   if (!is.null(pars)) attr(data[[i]],"pars") <- pars
    # }
    if(is.null(attr(design[[i]], "custom_ll"))){
      dadm_list[[i]] <- design_model(data=data[[i]],design=design[[i]],
                                     compress=compress,model=model[[i]],rt_resolution=rt_resolution[i])
      sampled_p_names <- names(attr(design[[i]],"p_vector"))
    } else{
      dadm_list[[i]] <- design_model_custom_ll(data = data[[i]],
                                               design = design[[i]],model=model[[i]])
      sampled_p_names <- attr(design[[i]],"sampled_p_names")
    }

    if(!is.null(prior_list[[i]])){
      prior_list[[i]] <- check_prior(prior_list[[i]], sampled_p_names)
    }
    # create a design model

  }
  prior <- merge_priors(prior_list)

  attr(dadm_list[[1]], "prior") <- prior

  # if(!is.null(subject_covariates)) attr(dadm_list, "subject_covariates") <- subject_covariates
  variant_funs <- get_variant_funs(type = type)
  if (type %in% c("standard", "single", "diagonal", "infnt_factor")) {
    out <- pmwgs(dadm_list, variant_funs, nuisance = nuisance, nuisance_non_hyper =
                   nuisance_non_hyper, grouped_pars = grouped_pars, n_factors = n_factors)
  } else if (type == "blocked") {
    if (is.null(par_groups)) stop("Must specify par_groups for blocked type")
    out <- pmwgs(dadm_list,par_groups=par_groups, variant_funs, nuisance = nuisance,
                 nuisance_non_hyper = nuisance_non_hyper, grouped_pars = grouped_pars)
  } else if (type == "lm") {
    out <- pmwgs(dadm_list,variant_funs, formula = formula, aggr_data = aggr_data,
                 nuisance = nuisance, nuisance_non_hyper = nuisance_non_hyper, grouped_pars = grouped_pars)
  } else if (type == "factor") {
    if (is.null(n_factors)) stop("Must specify n_factors for factor type")
    out <- pmwgs(dadm_list,variant_funs, n_factors = n_factors, nuisance = nuisance,
                 nuisance_non_hyper = nuisance_non_hyper, grouped_pars = grouped_pars,
                 constraintMat = constraintMat)
  } else if (type == "SEM"){
    out <- pmwgs(dadm_list,variant_funs, Lambda_mat = Lambda_mat, B_mat = B_mat, K_mat = K_mat, G_mat = G_mat,
                 xeta = xeta, xy = xy, nuisance = nuisance,
                 nuisance_non_hyper = nuisance_non_hyper, grouped_pars = grouped_pars)
  }
  # replicate chains
  dadm_lists <- rep(list(out),n_chains)
  # For post predict
  attr(dadm_lists,"data_list") <- data
  attr(dadm_lists,"design_list") <- design
  attr(dadm_lists,"model_list") <- model
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

run_IS2 <- function(samplers, filter = "sample", subfilter = 0, IS_samples = 1000,
                    stepsize_particles = 500, max_particles = 5000, n_cores = 1, df = 5){
  samples_merged <- merge_samples(samplers)
  IS_samples <- IS2(samples_merged, filter, subfilter = subfilter, IS_samples, stepsize_particles, max_particles, n_cores, df)
  attr(samplers, "IS_samples") <- IS_samples
  return(samplers)
}


auto_mclapply <- function(X, FUN, mc.cores, ...){
  if(Sys.info()[1] == "Windows"){
    cluster <- parallel::makeCluster(mc.cores)
    list_out <- parallel::parLapply(cl = cluster, X,FUN, ...)
    parallel::stopCluster(cluster)
  } else{
    list_out <- parallel::mclapply(X, FUN, mc.cores = mc.cores, ...)
  }
  return(list_out)
}
