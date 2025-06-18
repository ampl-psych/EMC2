get_stop_criteria <- function(stage, stop_criteria, type){
  if(is.null(stop_criteria)){
    if(stage == "preburn"){
      stop_criteria$iter <- 50
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
      stop_criteria$min_unique <- 150
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
    if(is.null(stop_criteria$omit_mpsrf)) stop_criteria$omit_mpsrf <- TRUE
  }
  if(stage == "adapt" & is.null(stop_criteria$min_unique)) stop_criteria$min_unique <- 600
  if(stage != "adapt" & !is.null(stop_criteria$min_unique)) stop("min_unique only applicable for adapt stage, try min_es instead.")
  return(stop_criteria)
}

#' Fine-Tuned Model Estimation
#'
#' Although typically users will rely on ``fit``, this function can be used for more fine-tuned specification of estimation needs.
#' The function will throw an error if a stage is skipped,
#' the stages have to be run in order ("preburn", "burn", "adapt", "sample").
#' More details can be found in the ``fit`` help files (``?fit``).
#'
#' @param emc An emc object
#' @param stage A string. Indicates which stage is to be run, either `preburn`, `burn`, `adapt` or `sample`
#' @param search_width A double. Tunes target acceptance probability of the MCMC process.
#' This fine-tunes the width of the search space to obtain the desired acceptance probability.
#' 1 is the default width, increases lead to broader search.
#' @param step_size An integer. After each step, the stopping requirements as
#' specified by `stop_criteria` are checked and proposal distributions are updated. Defaults to 100.
#' @param verbose Logical. Whether to print messages between each step with the current status regarding the stop_criteria.
#' @param verboseProgress Logical. Whether to print a progress bar within each step or not. Will print one progress bar for each chain and only if cores_for_chains = 1.
#' @param fileName A string. If specified will autosave emc at this location on every iteration.
#' @param particle_factor An integer. `particle_factor` multiplied by the square root of the number of sampled parameters determines the number of particles used.
#' @param cores_per_chain An integer. How many cores to use per chain.
#' Parallelizes across participant calculations. Only available on Linux or Mac OS.
#' For Windows, only parallelization across chains (``cores_for_chains``) is available.
#' @param cores_for_chains An integer. How many cores to use across chains.
#' Defaults to the number of chains. the total number of cores used is equal to ``cores_per_chain`` * ``cores_for_chains``.
#' @param max_tries An integer. How many times should it try to meet the finish
#' conditions as specified by stop_criteria? Defaults to 20. max_tries is ignored if the required number of iterations has not been reached yet.
#' @param n_blocks An integer. Number of blocks. Will block the parameter chains such that they are updated in blocks. This can be helpful in extremely tough models with a large number of parameters.
#' @param stop_criteria A list. Defines the stopping criteria and for which types of parameters these should hold. See ``?fit``.
#' @param thin A boolean. If `TRUE` will automatically thin the MCMC samples, closely matched to the ESS.
#' Can also be set to a double, in which case 1/thin of the chain will be removed (does not have to be an integer).
#' @param trim A boolean. If `TRUE` will automatically remove redundant samples (i.e. from preburn, burn, adapt).
#' @export
#' @return An emc object
#' @examples \donttest{
#' # First define a design
#' design_in <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1),
#'                            constants=c(s=log(1)))
#' # Then make the emc, we've omitted a prior here for brevity so default priors will be used.
#' emc <- make_emc(forstmann, design_in)
#'
#' # Now for example we can specify that we only want to run the "preburn" phase
#' # for MCMC 10 iterations
#' emc <- run_emc(emc, stage = "preburn", stop_criteria = list(iter = 10), cores_for_chains = 1)
#'}

run_emc <- function(emc, stage, stop_criteria,
                    search_width = 1, step_size = 100, verbose = FALSE, verboseProgress = FALSE,
                    fileName = NULL,particle_factor=50, cores_per_chain = 1,
                    cores_for_chains = length(emc), max_tries = 20, n_blocks = 1,
                    thin = FALSE, trim = TRUE){
  emc <- restore_duplicates(emc)
  if(Sys.info()[1] == "Windows" & cores_per_chain > 1) stop("only cores_for_chains can be set on Windows")
  if (verbose) message(paste0("Running ", stage, " stage"))
  total_iters_stage <- chain_n(emc)[,stage][1]
  if(stage != "preburn"){
    iter <- stop_criteria[["iter"]] + total_iters_stage
  } else{
    iter <- stop_criteria[["iter"]]
  }
  progress <- check_progress(emc, stage, iter, stop_criteria, max_tries, step_size, cores_per_chain*cores_for_chains, verbose, n_blocks = n_blocks)
  emc <- progress$emc
  progress <- progress[!names(progress) == 'emc']
  # We need to multiply step_size by thin to make an accurate guess for good step_size.
  cur_thin <- ifelse(is.numeric(thin), thin, 1)
  while(!progress$done){
    emc <- reset_pm_settings(emc, stage)
    # Remove redundant samples
    if(trim){
      emc <- fit_remove_samples(emc)
    }
    if(!is.null(progress$n_blocks)) n_blocks <- progress$n_blocks
    emc <- add_proposals(emc, stage, cores_per_chain*cores_for_chains, n_blocks)
    last_stage <- get_last_stage(emc)
    if(stage == "preburn"){
      sub_emc <- emc
    } else if(stage != last_stage){
      sub_emc <- subset(emc, filter = chain_n(emc)[1,last_stage]-1, stage = last_stage)
    } else{
      sub_emc <- subset(emc, filter = chain_n(emc)[1,stage] - 1, stage = stage)
    }
    # Actual sampling
    sub_emc <- auto_mclapply(sub_emc,run_stages, stage = stage, iter= progress$step_size*max(1,cur_thin),
                             verbose=verbose,  verboseProgress = verboseProgress,
                             particle_factor=particle_factor,search_width=search_width,
                             n_cores=cores_per_chain, mc.cores = cores_for_chains)

    class(sub_emc) <- "emc"
    if(stage != 'preburn'){
      if(is.numeric(thin)){
        sub_emc <- subset(sub_emc, stage = c("preburn", "burn", "adapt", "sample"), thin = thin)
      } else if(thin){
        sub_emc <- auto_thin(sub_emc, stage = c("preburn", "burn", "adapt", "sample"))
        # Update current rough guess for thinning:
        cur_thin <- progress$step_size/chain_n(sub_emc)[1,stage]
      }
      emc <- concat_emc(emc, sub_emc, progress$step_size, stage)
    } else{
      emc <- sub_emc
    }
    progress <- check_progress(emc, stage, iter, stop_criteria, max_tries, step_size, cores_per_chain*cores_for_chains,
                               verbose, progress,n_blocks)
    emc <- progress$emc
    progress <- progress[!names(progress) == 'emc']
    if(!is.null(fileName)){
      emc <- strip_duplicates(emc)
      fileName <- fix_fileName(fileName)
      class(emc) <- "emc"
      save(emc, file = fileName)
      emc <- restore_duplicates(emc)
    }
  }
  emc <- strip_duplicates(emc)
  class(emc) <- "emc"
  return(emc)
}

run_stages <- function(sampler, stage = "preburn", iter=0, verbose = TRUE, verboseProgress = TRUE,
                       particle_factor=50, search_width= NULL, n_cores=1)
{
  particles <- round(particle_factor*sqrt(sampler$n_pars))
  if (!sampler$init) {
    sampler <- init(sampler, n_cores = n_cores)
  }
  if (iter == 0) return(sampler)
  tune <- list(search_width = search_width)
  sampler <- run_stage(sampler, stage = stage,iter = iter, particles = particles,
                       n_cores = n_cores, tune = tune, verbose = verbose, verboseProgress = verboseProgress)
  return(sampler)
}

add_proposals <- function(emc, stage, n_cores, n_blocks){
  if(stage != "preburn"){
    # if(!is.null(emc[[1]]$g_map_fixed)){
    #   emc <- create_chain_proposals_lm(emc)
    # } else{    }
    emc <- create_chain_proposals(emc, do_block = stage != "sample")
    if(!is.null(n_blocks)){
      if(n_blocks > 1){
        components <- sub_blocking(emc, n_blocks)
        for(i in 1:length(emc)){
          attr(emc[[i]]$data, "components") <- components
        }
      }
    }
  }
  if(stage == "sample"){
    # if(!is.null(emc[[1]]$g_map_fixed)){
    #   emc <- create_eff_proposals_lm(emc, n_cores)
    # } else{    }
    emc <- create_eff_proposals(emc, n_cores)
  }
  return(emc)
}

check_progress <- function (emc, stage, iter, stop_criteria,
                            max_tries, step_size, n_cores, verbose, progress = NULL,
                            n_blocks)
{
  min_es <- stop_criteria$min_es
  if(is.null(min_es)) min_es <- 0
  selection <- stop_criteria$selection
  min_unique <- stop_criteria$min_unique
  total_iters_stage <- chain_n(emc)[, stage][1]
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
  gd <- check_gd(emc, stage, stop_criteria[["max_gd"]], stop_criteria[["mean_gd"]], trys, verbose,
                 iter = total_iters_stage, selection, omit_mpsrf = stop_criteria[["omit_mpsrf"]],
                 n_blocks)
  iter_done <- ifelse(is.null(iter) || length(iter) == 0, TRUE, total_iters_stage >= iter)
  if (min_es == 0) {
    es_done <- TRUE
  } else if (total_iters_stage != 0) {
    class(emc) <- "emc"
    curr_min_es <- Inf
    for(select in selection){
      curr_min_es <- min(c(ess_summary(emc, selection = select,
                                                stage = stage, stat_only = TRUE), curr_min_es))
    }
    if (verbose)
      message("Smallest effective size = ", round(curr_min_es))
    es_done <- ifelse(!emc[[1]]$init, FALSE, curr_min_es >
                        min_es)
  }
  else {
    es_done <- FALSE
  }
  trys_done <- ifelse(is.null(max_tries), FALSE, trys >= max_tries)
  if (stage == "adapt") {
    samples_merged <- merge_chains(emc)
    test_samples <- extract_samples(samples_merged, stage = "adapt",
                                    samples_merged$samples$idx, n_chains = length(emc))
    # if(!is.null(emc[[1]]$g_map_fixed)){
    #   adapted <- test_adapted_lm(emc[[1]], test_samples, min_unique, n_cores, verbose)
    # } else{    }
    adapted <- test_adapted(emc[[1]], test_samples,
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
  return(list(emc = gd$emc, done = done, step_size = step_size,
              trys = trys, n_blocks = gd$n_blocks))
}

check_gd <- function(emc, stage, max_gd, mean_gd, omit_mpsrf, trys, verbose,
                     selection, iter, n_blocks = 1)
{
  get_gds <- function(emc,omit_mpsrf, selection, stage) {
    gd_out <- c()
    alpha <- NULL
    for(select in selection){
      gd <- unlist(gd_summary.emc(emc, selection = select, stage = stage,
                                  omit_mpsrf = omit_mpsrf, stat = NULL))
      gd_out <- c(gd_out, c(gd))
      if(select == "alpha"){
        alpha <- gd
      }
    }
    gd_out[is.na(gd_out)] <- Inf
    return(list(gd = gd_out, alpha = alpha))
  }
  if(is.null(max_gd) & is.null(mean_gd)) return(list(gd_done = TRUE, emc = emc))
  if(!emc[[1]]$init | !stage %in% emc[[1]]$samples$stage)
    return(list(gd_done = FALSE, emc = emc))
  if(is.null(omit_mpsrf)) omit_mpsrf <- TRUE
  gd <- get_gds(emc,omit_mpsrf,selection, stage)
  alpha_gd <- gd$alpha
  gd <- gd$gd
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
    n_remove <- round(chain_n(emc)[,stage][1]/3)
    samplers_short <- subset.emc(emc, filter=n_remove,stage=stage, keep_stages = TRUE)
    if (is(samplers_short,"try-error")){
      gd_short <- Inf
    } else{
      gd_short <- get_gds(samplers_short,omit_mpsrf,selection, stage)
      alpha_gd_short <- gd_short$alpha
      gd_short <- gd_short$gd
    }
    if (is.null(max_gd) & (mean(gd_short) < mean(gd)) | (!is.null(max_gd) & (max(gd_short) < max(gd)))) {
      gd <- gd_short
      alpha_gd <- alpha_gd_short
      emc <- samplers_short
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
  if(verbose) {
    type <- "Rhat"
    if (!is.null(mean_gd)) message("Mean ",type," = ",round(mean(gd),3)) else
      if (!is.null(max_gd)) message("Max ",type," = ",round(max(gd),3))
  }
  emc <- set_tune_ess(emc, alpha_gd, mean_gd, max_gd)
  class(emc) <- "emc"
  return(list(gd = gd, gd_done = ok_gd, emc = emc))
}

set_tune_ess <- function(emc, alpha_gd = NULL, mean_gd = NULL, max_gd = NULL){
  if(is.null(alpha_gd)){
    return(emc)
  }
  alpha_gd <- matrix(alpha_gd, emc[[1]]$n_subjects)
  if(!is.null(mean_gd)){
    mean_alpha_ok <- rowMeans(alpha_gd) < 1+(mean_gd-1)*.5
  } else{
    mean_alpha_ok <- rep(T, nrow(alpha_gd))
  }
  if(!is.null(max_gd)){
    max_alpha_ok <- apply(alpha_gd, 1, max) < 1+(max_gd-1)*.5
  } else{
    max_alpha_ok <- rep(T, nrow(alpha_gd))
  }
  alpha_ok <- max_alpha_ok & mean_alpha_ok
  emc <- lapply(emc, function(x){
    pm_settings <- attr(x$samples, "pm_settings")
    pm_settings <- mapply(pm_settings, alpha_ok, FUN = function(y, z){
      for(i in 1:length(y)){
        y[[i]]$gd_good <- z
      }
      return(list(y))
    })
    attr(x$samples, "pm_settings") <- pm_settings
    return(x)
  })
  return(emc)
}


create_eff_proposals <- function(emc, n_cores){
  samples_merged <- merge_chains(emc)
  test_samples <- extract_samples(samples_merged, stage = c("adapt", "sample"), max_n_sample = 750, n_chains = length(emc))
  type <- emc[[1]]$type
  components <- attr(emc[[1]]$data, "components")
  for(i in 1:length(emc)){
    iteration = round(test_samples$iteration * i/length(emc))
    n_pars <- emc[[1]]$n_pars
    nuisance <- emc[[1]]$nuisance
    n_subjects <- emc[[1]]$n_subjects
    eff_mu <- matrix(0, nrow = n_pars, ncol = n_subjects)
    eff_var <- array(0, dim = c(n_pars, n_pars, n_subjects))
    for(comp in unique(components)){
      idx <- comp == components
      nuis_idx <- nuisance[idx]
      if(any(nuis_idx)){
        type <- samples_merged$sampler_nuis$type
        conditionals <- auto_mclapply(X = 1:n_subjects,
                                      FUN = get_conditionals, samples = test_samples,
                                      n_pars = sum(idx[!nuisance]), iteration =  iteration, idx = idx[!nuisance],
                                      type = type ,
                                      mc.cores = n_cores)
        conditionals_nuis <- auto_mclapply(X = 1:n_subjects,
                                           FUN = get_conditionals, samples = test_samples$nuisance,
                                           n_pars = sum(idx[nuisance]), iteration =  iteration, idx = idx[nuisance],
                                           type = type,
                                           mc.cores = n_cores)
        conditionals <- array(unlist(conditionals), dim = c(sum(idx[!nuisance]), sum(idx[!nuisance]) + 1, n_subjects))
        conditionals_nuis <- array(unlist(conditionals_nuis), dim = c(sum(idx[nuisance]), sum(idx[nuisance]) + 1, n_subjects))
        eff_mu[idx & !nuisance,] <- conditionals[,1,]
        eff_var[idx & !nuisance, idx & !nuisance,] <- conditionals[,2:(sum(idx[!nuisance])+1),]
        eff_mu[idx & nuisance,] <- conditionals_nuis[,1,]
        eff_var[idx & nuisance,idx & nuisance,] <- conditionals_nuis[,2:(sum(idx[nuisance])+1),]
      } else{
        conditionals <- auto_mclapply(X = 1:n_subjects,
                                      FUN = get_conditionals, samples = test_samples,
                                      n_pars = sum(idx[!nuisance]), iteration =  iteration, idx = idx[!nuisance],
                                      type = type,
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
    emc[[i]]$eff_mu <- eff_mu
    emc[[i]]$eff_var <- eff_var
    # attr(emc[[i]], "eff_alpha") <- eff_alpha
    # attr(emc[[i]], "eff_tau") <- eff_tau
  }
  return(emc)
}


sub_blocking <- function(emc, n_blocks){
  covs <- lapply(emc, FUN = function(x){return(x$chains_var)})
  out <- array(0, dim = dim(covs[[1]][[1]]))
  for(i in 1:length(covs)){
    cov_tmp <- covs[[1]]
    for(j in 1:length(cov_tmp)){
      out <- out + cov2cor(cov_tmp[[1]])
    }
  }
  shared_ll_idx <- attr(emc[[1]]$data, "shared_ll_idx")
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

create_chain_proposals <- function(emc, samples_idx = NULL, do_block = TRUE){
  n_subjects <- emc[[1]]$n_subjects
  n_chains <- length(emc)
  n_pars <- emc[[1]]$n_pars
  stage <- emc[[1]]$samples$stage[length(emc[[1]]$samples$stage)]
  if(is.null(samples_idx)){
    idx_subtract <- min(250, emc[[1]]$samples$idx/1.5)
    samples_idx <- round(emc[[1]]$samples$idx - idx_subtract):emc[[1]]$samples$idx
  }
  LL <- get_pars(emc, filter = samples_idx-1, selection = "LL",
                 stage = c('preburn', 'burn', 'adapt', 'sample'),
                 merge_chains = T, return_mcmc = F, remove_constants = F,
                 remove_dup = F)
  alpha <- get_pars(emc, filter = samples_idx-1, selection = "alpha",
                    stage = c('preburn', 'burn', 'adapt', 'sample'),
                    by_subject = T, merge_chains = T, return_mcmc = F,
                    remove_dup = F, remove_constants = F)

  components <- attr(emc[[1]]$data, "components")
  block_idx <- block_variance_idx(components)
  for(j in 1:n_chains){
    # Take only a random half of all the samples for each chain
    rnd_index <- sample(1:ncol(LL), round(ncol(LL)/2))
    chains_var <- vector("list", n_subjects)
    chains_mu <- vector("list", n_subjects)
    for(sub in 1:n_subjects){
      moments <- weighted_moments(alpha[,sub,rnd_index], LL[sub, rnd_index])
      emp_covs <- moments$w_cov
      chains_mu[[sub]] <- moments$w_mu
      if(do_block) emp_covs[block_idx] <- 0
      if(!is.positive.definite(emp_covs)){
        # If not positive definite, do not use it
        next
      } else{
        chains_var[[sub]] <- emp_covs
      }
    }
    # Instead use the mean of all other subjects
    null_idx <- sapply(chains_var, is.null)
    if(all(null_idx)){
      # If all subjects are NULL just come up with something
      mean_chains_var <- diag(n_pars) * .5
    } else{
      mean_chains_var <- Reduce(`+`, chains_var[!null_idx]) / sum(!null_idx)
    }
    if(any(null_idx)){
      for(q in 1:n_subjects){
        if(null_idx[q]) chains_var[[q]] <- mean_chains_var
      }
    }

    new_prop_var <- mean(diag(mean_chains_var))
    if(is.null(attr(emc[[j]], "prop_var"))){
      # This is in preburn case, group-level proposals are wider
      # So we scale the epsilon a bit to account for narrow individual proposals
      prop_var_ratio <- 2
    } else{
      prop_var_ratio <- attr(emc[[j]], "prop_var")/new_prop_var
    }
    if(stage != "sample"){
      emc[[j]] <- update_epsilon_scale(emc[[j]], prop_var_ratio)
    }
    attr(emc[[j]], "prop_var") <- new_prop_var
    emc[[j]]$chains_var <- chains_var
    emc[[j]]$chains_mu <- chains_mu
  }
  return(emc)
}

reset_pm_settings <- function(emc, stage){
  if(stage != get_last_stage(emc) || stage == "burn"){ # In this case we're running a new stage
    for(i in 1:length(emc)){
      pm_settings <- attr(emc[[i]]$samples, "pm_settings")
      attr(emc[[i]]$samples, "pm_settings") <- lapply(pm_settings, function(x){
        for(i in 1:length(x)){
          x[[i]]$proposal_counts <- rep(0, length(x[[i]]$proposal_counts))
          x[[i]]$acc_counts <- rep(0, length(x[[i]]$proposal_counts))
          if(stage != get_last_stage(emc)){
            x[[i]]$iter <- 25
            x[[i]]$mix <- NULL
          }
        }
        return(x)
      })
    }
  }
  return(emc)
}

update_epsilon_scale <- function(pmwgs, prop_var_ratio){
  pm_settings <- attr(pmwgs$samples, "pm_settings")
  pm_settings <- lapply(pm_settings, function(x){
    for(i in 1:length(x)){
      x[[i]]$epsilon <- x[[i]]$epsilon * prop_var_ratio
    }
    return(x)
  })
  attr(pmwgs$samples, "pm_settings") <- pm_settings
  return(pmwgs)
}

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
  components <- attr(sampler$data, "components")
  nuisance <- sampler$nuisance
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
                        FUN = get_conditionals, samples = test_samples$nuisance,
                        n_pars = sum(idx[nuisance]), idx = idx[nuisance],
                        type = type,
                        mc.cores = n_cores_conditional)
        }
        auto_mclapply(X = 1:sampler$n_subjects,FUN = get_conditionals,samples = test_samples,
                      n_pars = sum(idx[!nuisance]), idx = idx[!nuisance], type = sampler$type,
                      mc.cores = n_cores_conditional)
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


#' Make an emc Object
#'
#' Creates an emc object by combining the data, prior,
#' and model specification into a `emc` object that is needed in `fit()`.
#'
#' @param data A data frame, or a list of data frames. Needs to have the variable `subjects` as participant identifier.
#' @param design A list with a pre-specified design, the output of `design()`.
#' @param model A model list. If none is supplied, the model specified in `design()` is used.
#' @param type A string indicating whether to run a `standard` group-level, `blocked`, `diagonal`, `factor`, or `single` (i.e., non-hierarchical) model.
#' @param n_chains An integer. Specifies the number of mcmc chains to be run (has to be more than 1 to compute `rhat`).
#' @param compress A Boolean, if `TRUE` (i.e., the default), the data is compressed to speed up likelihood calculations.
#' @param rt_resolution A double. Used for compression, response times will be binned based on this resolution.
#' @param group_design A design for group-level mappings, made using `group_design()`.
#' @param par_groups A vector. Indicates which parameters are allowed to correlate. Could either be a list of character vectors of covariance blocks. Or
#' a numeric vector, e.g., `c(1,1,1,2,2)` means the covariances
#' of the first three and of the last two parameters are estimated as two separate blocks.
#' @param prior_list A named list containing the prior. Default prior created if `NULL`. For the default priors, see `?get_prior_{type}`.
#' @param ... Additional, optional arguments.
#' @return An uninitialized emc object
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
#' design_LBABE <- design(data = dat,model=LBA,matchfun=matchfun,
#' formula=list(v~lM,sv~lM,B~E+lR,A~1,t0~1),
#' contrasts=list(v=list(lM=ADmat)),constants=c(sv=log(1)))
#'
#' # specify priors
#' pmean <- c(v=1,v_lMdiff=1,sv_lMTRUE=log(.5), B=log(.5),B_Eneutral=log(1.5),
#'            B_Eaccuracy=log(2),B_lRright=0, A=log(0.25),t0=log(.2))
#' psd <- c(v=1,v_lMdiff=0.5,sv_lMTRUE=.5,
#'          B=0.3,B_Eneutral=0.3,B_Eaccuracy=0.3,B_lRright=0.3,A=0.4,t0=.5)
#' prior_LBABE <- prior(design_LBABE, type = 'standard',pmean=pmean,psd=psd)
#'
#' # create emc object
#' LBABE <- make_emc(dat,design_LBABE,type="standard",  prior=prior_LBABE)
#'
#' @export

make_emc <- function(data,design,model=NULL,
                          type="standard",
                          n_chains=3,compress=TRUE,rt_resolution=0.02,
                          prior_list = NULL, group_design = NULL,
                          par_groups=NULL, ...){
  # arguments for future compatibility
  n_factors <- NULL
  formula <- NULL
  Lambda_mat <- NULL
  B_mat <- NULL
  K_mat <- NULL
  G_mat <- NULL
  covariates <- NULL
  nuisance <- NULL
  nuisance_non_hyper <- NULL
  # overwrite those that were supplied
  optionals <- list(...)
  for (name in names(optionals) ) {
    assign(name, optionals[[name]])
  }
  if(!is.null(prior_list) & !is.null(prior_list$theta_mu_mean)){
    prior_list <- list(prior_list)
  }
  if(!is.null(prior_list)){
    type <- attr(prior_list[[1]], "type")
  }
  if(!is.null(group_design) && !type %in% c("standard", "diagonal", "blocked")){
    stop("group_design can only be used with standard, blocked or diagonal type")
  }
  if(type != "single" && length(unique(data$subjects)) == 1){
    stop("can only use type = `single` when there's only one subject in the data")
  }

  if (!(type %in% c("standard","diagonal","blocked","factor","single", "infnt_factor", "SEM", "diagonal-gamma")))
    stop("type must be one of: standard,diagonal,blocked,factor,infnt_factor, single")

  if(!is.null(nuisance) & !is.null(nuisance_non_hyper)){
    stop("You can only specify nuisance OR nuisance_non_hyper")
  }
  if (is(data, "data.frame")) data <- list(data)
  data <- lapply(data,function(d){
    d$subjects <- factor(d$subjects)
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
  if (!is.null(names(design)[1]) && names(design)[1]=="Flist"){
    design <- list(design)
  }
  checks <- sapply(design, function(x) is(x, "emc.design"))
  if (!all(checks)) stop("design must be a list of emc.design objects")
  if (length(design)!=length(data)){
    design <- rep(design,length(data))
  }
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
    if(is.null(attr(design[[i]], "custom_ll"))){
      dadm_list[[i]] <- design_model(data=data[[i]],design=design[[i]],
                                     compress=compress,model=model[[i]],rt_resolution=rt_resolution[i])
      sampled_p_names <- names(attr(design[[i]],"p_vector"))
    } else{
      dadm_list[[i]] <- design_model_custom_ll(data = data[[i]],
                                               design = design[[i]],model=model[[i]])
      sampled_p_names <- attr(design[[i]],"sampled_p_names")
    }
    if(length(prior_list) == length(data)){
      if(!is.null(prior_list[[i]])){
        prior_list[[i]] <- check_prior(prior_list[[i]], sampled_p_names, group_design)
      }
    }
  }
  # Make sure class retains following changes
  class(design) <- "emc.design"
  prior_in <- merge_priors(prior_list)

  prior_in <- prior(design, type, update = prior_in, group_design = group_design, ...)
  attr(dadm_list[[1]], "prior") <- prior_in

  # if(!is.null(subject_covariates)) attr(dadm_list, "subject_covariates") <- subject_covariates
  if (type %in% c("single", "infnt_factor", "diagonal-gamma")) {
    out <- pmwgs(dadm_list, type, nuisance = nuisance,
                 nuisance_non_hyper = nuisance_non_hyper,
                 n_factors = n_factors)
  } else if (type %in% c("standard", "blocked", "diagonal")) {
    if(type == "blocked"){
      if(is.null(par_groups)) stop("par_groups must be specified for blocked models")
    }
    if(type == "diagonal" && is.null(par_groups)){
      par_groups <- 1:length(sampled_pars(design))
    }
    if(type == "standard" && is.null(par_groups)){
      par_groups <- rep(1, length(sampled_pars(design)))
    }
    if(!is.null(par_groups)){
      if(is.character(par_groups)){
        par_groups <- list(par_groups)
      }
      if(is.list(par_groups)){
        par_names <- names(sampled_pars(design))
        new_par_groups <- rep(NA, length(par_names))
        for(i in 1:length(par_groups)){
          if(any(!par_groups[[i]] %in% par_names)) stop("Make sure you specified parameter names in par_groups correctly")
          new_par_groups[par_names %in% par_groups[[i]]] <- i
        }
        new_par_groups[is.na(new_par_groups)] <- (i+1):(i+sum(is.na(new_par_groups)))
        par_groups <- new_par_groups
      }
      if(length(par_groups) != length(sampled_pars(design))){
        stop("par_groups length does not match number of sampled parameters, make sure you specified par_groups correctly")
      }
    }
    if(type %in% c("diagonal", "blocked")) type <- "standard"
    out <- pmwgs(dadm_list, type, par_groups=par_groups,
                 group_design = group_design,
                 nuisance = nuisance,
                 nuisance_non_hyper = nuisance_non_hyper)
  } else if (type == "factor") {
    out <- pmwgs(dadm_list, type, n_factors = n_factors,
                 nuisance = nuisance,
                 nuisance_non_hyper = nuisance_non_hyper,
                 Lambda_mat = Lambda_mat)
  } else if (type == "SEM"){
    out <- pmwgs(dadm_list, type, Lambda_mat = Lambda_mat,
                 B_mat = B_mat, K_mat = K_mat, G_mat = G_mat,
                 covariates = covariates, nuisance = nuisance,
                 nuisance_non_hyper = nuisance_non_hyper)
  }
  out$model <- lapply(design, function(x) x$model)
  # Only for joint models we need to keep a list of functions
  if(length(out$model) == 1) out$model <- out$model[[1]]
  # out <- check_duplicate_designs(out)
  # replicate chains
  dadm_lists <- rep(list(out),n_chains)
  # For post predict
  class(dadm_lists) <- "emc"
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

check_duplicate_designs <- function(out){
  if(is.data.frame(out$data[[1]])) return(out)
  for(i in 1:length(out$data)){ # loop over subjects
    designs <- lapply(out$data[[i]], function(y) attr(y, "designs"))
    duplicacy <- duplicated(designs)
    for(j in 1:length(out$data[[i]])){# Loop over data sets in this sub
      if(duplicacy[j]){
        attr(out$data[[i]][[j]], "designs") <- unq_idx[j]
      }
    }
  }
  return(out)
}

extractDadms <- function(dadms, names = NULL){
  if(is.null(names)) names <- 1:length(dadms)
  N_models <- length(dadms)
  pars <- attr(dadms[[1]], "sampled_p_names")
  prior <- attr(dadms[[1]], "prior")
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
    }
    dadm_list <- do.call(mapply, c(list, total_dadm_list, SIMPLIFY = F))
  }
  # subject_covariates_ok <- unlist(lapply(subject_covariates, FUN = function(x) length(x) == length(subjects)))
  # if(!is.null(subject_covariates_ok)) if(any(!subject_covariates_ok)) stop("subject_covariates must be as long as the number of subjects")
  attr(dadm_list, "components") <- components
  attr(dadm_list, "shared_ll_idx") <- components
  return(list(prior = prior,
              dadm_list = dadm_list, subjects = subjects))
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

#' Strip all entries except samples from EMC list entries
#' @param emc A list of EMC objects
#' @return The same list with everything but samples removed from all but first entry
#' @noRd
strip_duplicates <- function(emc, incl_props = TRUE) {
  # Keep only samples in non-first entries
  for (i in 2:length(emc)) {
    samples <- emc[[i]]$samples
    prop_var <- attr(emc[[i]], "prop_var")
    emc[[i]] <- list(samples = samples)
    attr(emc[[i]], "prop_var") <- prop_var
  }
  if(incl_props){
    # Also remove eff_mu, eff_var, chains_cov
    for (i in 1:length(emc)) {
      emc[[i]]$eff_mu <- NULL
      emc[[i]]$eff_var <- NULL
      emc[[i]]$chains_cov <- NULL
      emc[[i]]$chains_mu <- NULL
    }
  }
  return(emc)
}

get_posterior_weights <- function(ll){
  max_ll <- max(ll)
  weights <- exp(ll - max_ll)
  weights <- weights / sum(weights)
  return(weights)
}

weighted_moments <- function(chain, ll = NULL) {
  # chain: a matrix with each col as a sample, and each row as a parameter
  # weights: an optional vector of weights. If not provided, equal weighting is assumed.
  n <- ncol(chain)
  d <- nrow(chain)

  # Use equal weights if none are provided.
  if (is.null(ll)) {
    weights <- rep(1 / n, n)
  } else {
    weights <- get_posterior_weights(ll)
  }

  # Compute the weighted mean of the chain.
  weighted_mean <- colSums(apply(chain, 1, function(x) x*weights))
  # NIEK THIS CAUSES ERRORS
  # # Compute the weighted covariance matrix.
  # cov_matrix <- matrix(0, nrow = d, ncol = d)
  # for (i in 1:n) {
  #   diff <- chain[,i] - weighted_mean
  #   cov_matrix <- cov_matrix + weights[i] * (diff %*% t(diff))
  # }
  cov_matrix <- cov(t(chain))
  return(list(w_cov = cov_matrix, w_mu = weighted_mean))
}



#' Restore all entries to EMC list entries from first entry
#' @param emc A list of EMC objects with stripped duplicates
#' @return The same list with all fields restored from first entry except samples
#' @noRd
restore_duplicates <- function(emc) {
  # Restore everything except samples from first entry
  for (i in 2:length(emc)) {
    samples <- emc[[i]]$samples
    emc[[i]] <- emc[[1]]
    emc[[i]]$samples <- samples
  }
  return(emc)
}
