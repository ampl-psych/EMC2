#' Simulation-Based Calibration
#'
#' Runs SBC for an EMC2 model and associated design. Returns
#' normalized rank (between 0 and 1) and prior samples. For hierarchical models the group-level mean and
#' the (implied) group-level (co-)variance are returned.
#' For non-hierarchical models only the subject-level parameters rank is returned.
#'
#' @param design_in An emc design list. The design of the model to be used in SBC
#' @param prior_in An emc prior list. The prior for the design to be used in SBC
#' @param replicates Integer. The number of samples to draw from the prior
#' @param trials Integer. The number of trials of the simulated data (per subject)
#' @param n_subjects Integer. Only used for hierarchical models. The number of subjects to be used in data generation of each replicate
#' @param plot_data Boolean. Whether to plot the data simulated (aggregated across subjects)
#' @param verbose Verbose. Whether to print progress related messages
#' @param fileName Character. Highly recommended, saves temporary results to the fileName
#' @param ... A list of optional additional arguments that can be passed to `fit` and `make_emc`
#'
#' @return The ranks and prior samples. For hierarchical models also the prior-generated subject-level parameters.
#' @export
run_sbc <- function(design_in, prior_in, replicates = 250, trials = 100, n_subjects = 30,
                    plot_data = FALSE, verbose = TRUE,
                    fileName = NULL, ...){
  if(is.null(fileName)) message("Since SBC can take a while it's highly recommended to specify a fileName to save temporary results in case of crashes")
  type <- attr(prior_in, "type")
  if(type == "single"){
    out <- SBC_single(design_in, prior_in, replicates, trials,
                      plot_data, verbose, fileName, ...)
  } else{
    out <- SBC_hierarchical(design_in, prior_in, replicates, trials, n_subjects,
                    plot_data, verbose, fileName, ...)
  }
  return(out)
}



SBC_hierarchical <- function(design_in, prior_in, replicates = 250, trials = 100, n_subjects = 30,
                             plot_data = FALSE, verbose = TRUE,
                             fileName = NULL, ...){
  dots <- add_defaults(list(...), max_tries = 50, compress = FALSE, rt_resolution = 1e-12,
                       stop_criteria = list(min_es = 100, max_gd = 1.1,
                                            selection = c("alpha", "mu", "Sigma")))
  dots$verbose <- verbose
  type <- attr(prior_in, "type")
  if(type != "diagonal-gamma") warning("SBC in EMC2 is frequently biased for prior structures with heavy variance mass around 0 \n
                                       please consider using 'type = diagonal-gamma' to test your model. See also vignettes")
  # Draw prior samples
  # Should at a later point go to predict
  prior_mu <- plot(prior_in, design_in, do_plot = F, N = replicates, selection = "mu", return_mcmc = FALSE, map = FALSE)[[1]]
  prior_var <- plot(prior_in, design_in, do_plot = F, N = replicates, selection = "Sigma", return_mcmc = FALSE,
                          remove_constants = FALSE)[[1]]
  rank_mu <- data.frame()
  rank_var <- data.frame()
  par_names <- names(sampled_pars(design_in))
  if(!is.null(fileName)) save(prior_mu, prior_var, file = fileName)
  all_rand_effects <- vector("list", length = replicates)
  for(i in 1:replicates){
    if(verbose) print(paste0("Sample ", i, " out of ", replicates))
    rand_effects <- make_random_effects(design_in, prior_mu[,i], n_subj = n_subjects, covariances = prior_var[,,i])
    all_rand_effects[[i]] <- rand_effects
    data <- make_data(rand_effects,design_in, trials, model = design_in$model)
    if(plot_data){
      plot_density(data, factors = names(design_in$Ffactors)[names(design_in$Ffactors) != "subjects"])
    }
    emc <-  do.call(make_emc, c(list(data = data, design = design_in, prior_list = prior_in, type = type), fix_dots(dots, make_emc)))
    emc <-  do.call(fit, c(list(emc = emc), fix_dots(dots, fit)))

    ESS_mu <- round(pmin(unlist(ess_summary(emc, selection = "mu", stat = NULL)), chain_n(emc)[1,"sample"]))
    mu_rec <- get_pars(emc, selection = "mu", return_mcmc = F, merge_chains = T,flatten = T)
    rank_mu <- rbind(rank_mu, mapply(get_ranks_ESS, split(mu_rec, row(mu_rec)), t(ESS_mu), prior_mu[,i]))

    ESS_var <- round(pmin(unlist(ess_summary(emc, selection = "Sigma", stat = NULL)), chain_n(emc)[1,"sample"]))
    var_rec <- get_pars(emc, selection = "Sigma", return_mcmc = F, merge_chains = T, flatten = T)
    prior_var_input <- get_pars(emc, selection = "Sigma", flatten = T, true_pars = prior_var[,,i])[[1]][[1]]
    rank_var <- rbind(rank_var, mapply(get_ranks_ESS, split(var_rec, row(var_rec)), t(ESS_var), prior_var_input))

    colnames(rank_mu) <- par_names
    colnames(rank_var) <- colnames(prior_var_input)
    if(!is.null(fileName)){
      SBC_temp <- list(rank = list(mu = rank_mu, var =  rank_var),
                     prior = list(mu = prior_mu, var = prior_var),
                     rand_effects = rand_effects, emc = emc)
      save(SBC_temp, file = fileName)
    }
  }
  return(list(rank = list(mu = rank_mu, var =  rank_var),
              prior = list(mu = prior_mu, var = prior_var),
              rand_effects = rand_effects))
}


SBC_single <- function(design_in, prior_in, replicates = 250, trials = 100,
                             plot_data = FALSE, verbose = TRUE,
                             fileName = NULL, ...){
  dots <- add_defaults(list(...), max_tries = 50, compress = FALSE, rt_resolution = 1e-12,
                       stop_criteria = list(min_es = 100, max_gd = 1.1,
                                            selection = c("alpha", "mu", "Sigma")))
  dots$verbose <- verbose
  type <- attr(prior_in, "type")
  if(type != "single") stop("can only use `type = single`")
  # Draw prior samples
  prior_alpha <- parameters(prior_in, N = replicates, selection = "alpha")
  rank_alpha <- data.frame()

  # DEBUG
  #  if(!is.null(fileName)) save(prior_alpha, file = fileName)
  save(prior_alpha,file=paste0("alpha",fileName))

  i <- 1
  if(dots$cores_per_chain > 1 & verbose) print("Since cores_per_chain > 1, estimating multiple data sets simultaneously")
  par_names <- names(sampled_pars(design_in))
  while(i < replicates){
    if(verbose) print(paste0("Sample ", i, " out of ", replicates))
    data <- data.frame()
    start <- i
    for(j in 1:dots$cores_per_chain){
      if(i > replicates) next
      design_in$Ffactors$subjects <- j
      tmp <- make_data(prior_alpha[i,],design_in, trials, model = design_in$model, ...)
      tmp$subjects <- j
      data <- rbind(data, tmp)
      i <- i + 1
    }

    # DEBUG
    save(data,file=paste0("dat",i,fileName))

    if(plot_data){
      plot_density(data, factors = names(design_in$Ffactors)[names(design_in$Ffactors) != "subjects"])
    }
    emc <-  do.call(make_emc, c(list(data = data, design = design_in, prior_list = prior_in, type = type), fix_dots(dots, make_emc)))
    emc <-  do.call(fit, c(list(emc = emc), fix_dots(dots, fit)))
    ESS <- round(pmin(do.call(cbind, ess_summary(emc, selection = "alpha", stat = NULL)), chain_n(emc)[1,"sample"]))
    alpha_rec <- get_pars(emc, selection = "alpha", return_mcmc = F, merge_chains = T, flatten = T)
    for(j in 1:dots$cores_per_chain){
      if(start > replicates) next
      tmp_rec <- alpha_rec[,j,]
      rank_alpha <- rbind(rank_alpha, mapply(get_ranks_ESS, split(tmp_rec, row(tmp_rec)), t(ESS[j,]), prior_alpha[start,]))
      start <- start + 1
    }
    colnames(rank_alpha) <- par_names
    if(!is.null(fileName)){
      SBC_temp <- list(rank = list(alpha = rank_alpha),
                       prior = list(alpha = prior_alpha), emc = emc)
      save(SBC_temp, file = fileName)
    }
  }
  return(list(rank = list(alpha = rank_alpha),
              prior = list(alpha = prior_alpha)))
}

#' Plot the Histogram of the Observed Rank Statistics of SBC
#'
#' Note that this plot is dependent on the number of bins, and a more general
#' visualization is to use `plot_sbc_ecdf`
#'
#' @param ranks A list of named dataframes of the rank statistic
#' @param bins An integer specifying the number of bins to use when plotting the histogram
#' @param layout Optional. A numeric vector specifying the layout using `par(mfrow = layout)`
#'
#' @return No returns
#' @export
plot_sbc_hist <- function(ranks, bins = 10, layout = NA){
  if(!is.null(ranks$rank)) ranks <- ranks$rank
  selects <- names(ranks)
  oldpar <- par(no.readonly = TRUE) # code line i
  on.exit(par(oldpar)) # code line i + 1

  n_sample <- nrow(ranks[[1]])
  low <- qbinom(0.025, n_sample, 1/bins)
  mid <- qbinom(0.5, n_sample, 1/bins)
  high <- qbinom(0.975, n_sample, 1/bins)
  par_names <- colnames(ranks[[1]])
  for(j in 1:length(ranks)){
    if(any(is.na(layout))){
      par(mfrow = coda_setmfrow(Nchains = 1, Nparms = ncol(ranks[[1]]),
                                       nplots = 1))
    } else{par(mfrow=layout)}
    rank <- ranks[[j]]
    for(i in 1:ncol(rank)){
      hist(rank[,i], main = paste0(selects[j], " - ", par_names[i]), breaks = bins, ylim = c(0, high + 2))
      abline(h = low, lty = 2)
      abline(h = mid, lty = 2)
      abline(h = high, lty = 2)
    }
  }
}



get_gamma <- function (N, K, conf_level = 0.95)
{
  p_interior <- function (p_int, x1, x2, z1, z2, gamma, N)
  {
    z_tilde <- (z2 - z1)/(1 - z1)
    N_tilde <- rep(N - x1, each = length(x2))
    p_int <- rep(p_int, each = length(x2))
    x_diff <- outer(x2, x1, "-")
    p_x2_int <- p_int * dbinom(x_diff, N_tilde, z_tilde)
    list(p_int = rowSums(p_x2_int), x1 = x2)
  }
  target <- function(gamma, conf_level, N, K) {
    z <- 1:(K - 1)/K
    z1 <- c(0, z)
    z2 <- c(z, 1)
    x2_lower <- qbinom(gamma/2, N, z2)
    x2_upper <- c(N - rev(x2_lower)[2:K], N)
    x1 <- 0
    p_int <- 1
    for (i in seq_along(z1)) {
      tmp <- p_interior(p_int, x1 = x1, x2 = x2_lower[i]:x2_upper[i],
                        z1 = z1[i], z2 = z2[i], gamma = gamma, N = N)
      x1 <- tmp$x1
      p_int <- tmp$p_int
    }
    abs(conf_level - sum(p_int))
  }
  optimize(target, c(0, 1 - conf_level), conf_level, N = N,
           K = K)$minimum
}

get_lims <- function (N, K, gamma)
{
  lims <- list()
  z <- seq(0, 1, length.out = K)
  lims$lower <- qbinom(gamma/2, N, z)/N - z
  lims$upper <- qbinom(1 - gamma/2, N, z)/N - z
  lims$z <- z
  lims
}

make_smooth <- function(x, y, N = 1000){
  lo <- smooth.spline(x, y, spar=0.5)
  xl <- seq(0, 1, 1/N)
  return(predict(lo,xl)$y)
}

#' Plot the ECDF Difference in SBC Ranks
#'
#' Plots the difference in observed cumulative rank statistics and the
#' expected cumulative distribution of a uniform distribution. The blue
#' shaded areas indicate the 95% credible interval.
#'
#' @param ranks A list of named dataframes of the rank statistic
#' @param layout Optional. A numeric vector specifying the layout using `par(mfrow = layout)`
#'
#' @return No returns
#' @export
plot_sbc_ecdf <- function(ranks, layout = NA){
  if(!is.null(ranks$rank)) ranks <- ranks$rank
  selects <- names(ranks)
  oldpar <- par(no.readonly = TRUE) # code line i
  on.exit(par(oldpar)) # code line i + 1

  K <- N <- nrow(ranks[[1]])
  gamma <- get_gamma(N, K)
  res <- get_lims(N, K, gamma)
  for(j in 1:length(ranks)){
    if(any(is.na(layout))){
      par(mfrow = coda_setmfrow(Nchains = 1, Nparms = ncol(ranks[[1]]),
                                nplots = length(ranks)))
    } else{par(mfrow=layout)}

    rank <- ranks[[j]]
    par_names <- colnames(rank)
    res$x <- apply(rank, 2, function(x) sort(x) - res$z)
    for(i in 1:ncol(rank)){
      plot(res$z, res$x[,i], type = "l", ylim = c(min(res$lower, res$x[,i]) - 0.01, max(res$upper,res$x[,i]) + 0.01), xlim = c(0, 1),
           lwd = 2, ylab = "ECDF Difference", xlab = "Normalized Rank Statistic", main = paste0(selects[j], " - ", par_names[i]))
      polygon(c(res$z, rev(res$z)), c(res$lower, rev(res$upper)), col = adjustcolor("cornflowerblue", 0.2))
    }
  }
}

get_ranks_ESS <- function(posterior, ESS, prior){
  posterior <- posterior[seq(1, length(posterior), length.out = ESS)]
  return(rank(c(prior, posterior))[1]/(ESS+1))
}
