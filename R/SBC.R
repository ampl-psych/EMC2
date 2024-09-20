#' Hierarchical simulation-based calibration
#'
#' Runs hierarchical SBC for an EMC2 model and associated design. Note that the
#' normalized rank [0,1] of the group-level mean and the (implied) group-level (co-)variance are returned
#'
#' @param design_in An emc design list. The design of the model to be used in SBC
#' @param prior_in An emc prior list. The prior for the design to be used in SBC
#' @param replicates Integer. The number of samples to draw from the prior
#' @param trials Integer. The number of trials of the simulated data (per subject)
#' @param n_subjects Integer. The number of subjects to be used in data generation of each replicate
#' @param plot_data Boolean. Whether to plot the data simulated (aggregated across subjects)
#' @param verbose Verbose. Whether to print progress related messages
#' @param n_post Integer. The number of posterior samples to be taken and calculated the rank across
#' @param fileName Character. Highly recommended, saves temporary results to the fileName
#' @param ... A list of optional additional arguments that can be passed to `fit` and `make_emc`
#'
#' @return A list with the rank of the group-level mean and (co-)variance, the prior samples and the generated subject-level parameters
#' @export
SBC_hierarchical <- function(design_in, prior_in, replicates = 250, trials = 100, n_subjects = 30,
                             plot_data = FALSE, verbose = TRUE,  n_post = 1000,
                             fileName = NULL, ...){
  dots <- add_defaults(list(...), max_tries = 50, compress = FALSE, rt_resolution = 1e-12)
  dots$verbose <- verbose
  type <- attr(prior_in, "type")
  if(type != "diagonal-gamma") warning("SBC in EMC2 is frequently biased for prior structures with heavy variance mass around 0 \n
                                       please consider using 'type = diagonal-gamma' to test your model. See also BLOGPOST")
  if(is.null(fileName)) warning("Since SBC can take a while it's highly recommended to specify a fileName to save temporary results in case of crashes")
  # Draw prior samples
  prior_mu <- plot_prior(prior_in, design_in, do_plot = F, N = replicates, selection = "mu", return_mcmc = FALSE, map = FALSE)[[1]]
  prior_var <- plot_prior(prior_in, design_in, do_plot = F, N = replicates, selection = "Sigma", return_mcmc = FALSE,
                          remove_constants = FALSE)[[1]]
  rank_mu <- vector("list", length = replicates)
  rank_var <- vector("list", length = replicates)
  attr(rank_mu, "selection") <- "mu"
  attr(rank_var, "selection") <- "Sigma"

  if(!is.null(fileName)) save(prior_mu, prior_var, file = fileName)
  all_rand_effects <- vector("list", length = replicates)
  for(i in 1:replicates){
    if(verbose) print(paste0("Sample ", i, " out of ", replicates))
    rand_effects <- make_random_effects(design_in, prior_mu[,i], n_subj = n_subjects, covariances = prior_var[,,i])
    all_rand_effects[[i]] <- rand_effects
    # if(hyper_only){
    #   rand_effects <- as.data.frame(rand_effects)
    #   rand_effects$subjects <- 1:nrow(rand_effects)
    #   emc <- EMC2:::run_hyper(type, rand_effects, prior = prior_in, iter = 5000)
    #   emc <- list(emc)
    #   attr(emc[[1]], "variant_funs") <- EMC2:::get_variant_funs(type)
    #   class(emc) <- "emc"
    # } else{
    data <- make_data(rand_effects,design_in, trials, model = design_in$model)
    if(plot_data){
      plot_defective_density(data, factors = names(design_in$Ffactors)[names(design_in$Ffactors) != "subjects"])
    }
    emc <-  do.call(make_emc, c(list(data = data, design = design_in, prior_list = prior_in, type = type), fix_dots(dots, make_emc)))
    emc <-  do.call(fit, c(list(emc = emc), fix_dots(dots, fit)))
    n_post_pc <- round((n_post-1)/length(emc))
    mu_rec <- get_pars(emc, selection = "mu", return_mcmc = F, merge_chains = T, length.out = n_post_pc, flatten = T)
    var_rec <- get_pars(emc, selection = "Sigma", return_mcmc = F, merge_chains = T, length.out = n_post_pc, flatten = T)
    prior_var_input <- get_pars(emc, selection = "Sigma", flatten = T, true_pars = prior_var[,,i])[[1]][[1]]
    actual_n_post <- (n_post_pc*length(emc))+1
    rank_mu[[i]] <- (apply(cbind(prior_mu[,i], mu_rec), 1, rank)[1,])/actual_n_post
    rank_var[[i]] <- (apply(cbind(t(prior_var_input), var_rec), 1, rank)[1,])/actual_n_post
    if(!is.null(fileName)) save(rank_mu, rank_var, emc, prior_mu, prior_var, all_rand_effects, file = fileName)
  }
  return(list(rank_mu = rank_mu, rank_var = rank_var,
              prior_mu = prior_mu, prior_var = prior_var,
              rand_effects = rand_effects))
}


#' Simulation-based calibration
#'
#' Runs single subject SBC for an EMC2 model and associated design. Note that the
#' normalzided rank [0,1] of the subject-level parameters is returned.
#' @inheritParams SBC_hierarchical

#' @return A list with the rank of the subject-level parameters and the prior-samples
#' @export
SBC_single <- function(design_in, prior_in, replicates = 250, trials = 100,
                             plot_data = FALSE, verbose = TRUE,  n_post = 1000,
                             fileName = NULL, ...){
  dots <- add_defaults(list(...), max_tries = 50, compress = FALSE, rt_resolution = 1e-12)
  dots$verbose <- verbose
  type <- attr(prior_in, "type")
  if(type != "single") stop("can only use `type = single`")
  if(is.null(fileName)) warning("Since SBC can take a while it's highly recommended to specify a fileName to save temporary results in case of crashes")
  # Draw prior samples
  prior_alpha <- plot_prior(prior_in, design_in, do_plot = F, N = replicates, selection = "mu", return_mcmc = FALSE, map = FALSE)[[1]]
  rank_alpha <- vector("list", length = replicates)
  attr(rank_alpha, "selection") <- "alpha"
  if(!is.null(fileName)) save(prior_alpha, file = fileName)
  for(i in 1:replicates){
    if(verbose) print(paste0("Sample ", i, " out of ", replicates))
    data <- make_data(t(prior_alpha[,,i]),design_in, trials, model = design_in$model)
    if(plot_data){
      plot_defective_density(data, factors = names(design_in$Ffactors)[names(design_in$Ffactors) != "subjects"])
    }
    emc <-  do.call(make_emc, c(list(data = data, design = design_in, prior_list = prior_in, type = type), fix_dots(dots, make_emc)))
    emc <-  do.call(fit, c(list(emc = emc), fix_dots(dots, fit)))
    n_post_pc <- round((n_post-1)/length(emc))
    alpha_rec <- get_pars(emc, selection = "alpha", return_mcmc = F, merge_chains = T, length.out = n_post_pc, flatten = T)[,1,]
    actual_n_post <- (n_post_pc*length(emc))+1
    rank_alpha[[i]] <- (apply(cbind(prior_alpha[,,i], alpha_rec), 1, rank)[1,])/actual_n_post
    if(!is.null(fileName)) save(rank_alpha, emc, prior_alpha, file = fileName)
  }
  return(list(rank_alpha = rank_alpha, prior_alpha = prior_alpha))
}

plot_sbc <- function(ranks, bins, layout = NA, ...){
  selection <- attr(ranks, "selection")
  if(!is.data.frame(ranks) && !is.matrix(ranks)) ranks <- do.call(rbind, ranks)
  oldpar <- par(no.readonly = TRUE) # code line i
  on.exit(par(oldpar)) # code line i + 1
  if(any(is.na(layout))){
    par(mfrow = coda_setmfrow(Nchains = 1, Nparms = ncol(ranks),
                              nplots = 1))
  } else{par(mfrow=layout)}
  n_sample <- nrow(ranks)
  low <- qbinom(0.025, n_sample, 1/bins)
  mid <- qbinom(0.5, n_sample, 1/bins)
  high <- qbinom(0.975, n_sample, 1/bins)
  par_names <- colnames(ranks)
  for(i in 1:ncol(ranks)){
    hist(ranks[,i], main = paste0(selection, " - ", par_names[i]), breaks = bins, ylim = c(0, high + 2))
    abline(h = low, lty = 2)
    abline(h = mid, lty = 2)
    abline(h = high, lty = 2)
  }
}








