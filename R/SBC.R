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
    data <- make_data(rand_effects,design_in, trials, model = design_in$model, ...)
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


run_SBC_subject <- function(rep, design_in, prior_alpha, trials, prior_in, dots){
  p_vector <- prior_alpha[rep,]
  data <- do.call(make_data, c(list(parameters = p_vector, design = design_in, n_trials = trials), fix_dots(dots, make_data)))
  emc <-  do.call(make_emc, c(list(data = data, design = design_in, prior_list = prior_in, type = "single"), fix_dots(dots, make_emc)))

  # Protection wrapper for the fit call
  fit_result <- tryCatch({
    do.call(fit, c(list(emc = emc), fix_dots(dots, fit)))
  }, warning = function(w) {
    # Save the problematic p_vector
    filename <- paste0("p_vector_rep", rep, ".Rdata")
    save(p_vector, file = filename)
    warning("A warning occurred during fitting for replication ", rep,
            ". The input parameters have been saved as ", filename,
            ". Warning message: ", w$message)
    # Re-throw the warning
    warning(w)
    # Return NULL to indicate failure
    return(NULL)
  }, error = function(e) {
    # Save the problematic p_vector
    filename <- paste0("p_vector_rep", rep, ".Rdata")
    save(p_vector, file = filename)
    warning("An error occurred during fitting for replication ", rep,
            ". The input parameters have been saved as ", filename,
            ". Error message: ", e$message)
    # Return NULL to indicate failure
    return(NULL)
  })

  # Check if fitting failed
  if (is.null(fit_result)) {
    # Return a structure that indicates failure but doesn't break the overall process
    return(list(rank = NULL, med = NULL, bias = NULL, coverage = NULL, failed = TRUE))
  }

  emc <- fit_result

  # Return ESS
  ESS <- ess_summary(emc, stat = NULL)[[1]]
  # Rank
  alpha_rec <- get_pars(emc, selection = "alpha", return_mcmc = F, merge_chains = T, flatten = T)[,1,]
  rank <- mapply(get_ranks_ESS, split(alpha_rec, row(alpha_rec)), t(ESS), p_vector)
  names(rank) <- names(sampled_pars(design_in))
  # For Bias
  CI <- credint(emc)[[1]]
  med <- CI[,2] # And return this for precision
  bias <- med - p_vector
  # For coverage
  coverage <- p_vector > CI[,1] & p_vector < CI[,3]
  return(list(rank = rank, med = med, bias = bias, coverage = coverage))
}

split_list_to_dfs <- function(lst, type = "alpha") {
  comps <- names(lst[[1]])   # get component names from the first element
  out <- lapply(comps, function(nm) {
    res <- list()
    res[[type]] <- do.call(rbind, lapply(lst, `[[`, nm))
    return(res)
  })
  names(out) <- comps
  out
}



SBC_single <- function(
    design_in,
    prior_in,
    replicates = 250,
    trials = 100,
    plot_data = FALSE,
    verbose = TRUE,
    fileName = NULL,
    ...
) {
  if (attr(prior_in, "type") != "single") {
    stop("can only use `type = single`")
  }
  dots <- add_defaults(
    list(...),
    max_tries = 50,
    compress = FALSE, rt_resolution = 1e-12,
    stop_criteria = list(
      min_es = 100, max_gd = 1.1, selection = c("alpha", "mu", "Sigma")
    ),
    cores_per_chain = 1
  )
  dots[["verbose"]] <- verbose
  # Draw prior samples
  prior_alpha <- parameters(prior_in, N = replicates, selection = "alpha")
  rank_alpha <- data.frame()
  if (!is.null(fileName)) {
    save(prior_alpha, file = fileName)
  }
  i <- 1
  if (dots[["cores_per_chain"]] > 1 && verbose) {
    print("Since cores_per_chain > 1, estimating multiple data sets simultaneously")
  }
  par_names <- names(sampled_pars(design_in))
  res <- auto_mclapply(
    X = 1:replicates,
    FUN = run_SBC_subject,
    design_in, prior_alpha, trials, prior_in, dots,
    mc.cores = dots[["cores_per_chain"]]
  )
  SBC <- split_list_to_dfs(res)
  if(!is.null(fileName)) {
    save(SBC, prior_alpha, file = fileName)
  }
  return(SBC)
}

calc_sbc_stats <- function(stats){
  if(is.null(stats$coverage)) return(NULL)
  out <- list()
  out_names <- names(stats[[1]])
  for(i in 1:length(stats[[1]])){
    out[[out_names[i]]] <- list(
      coverage = apply(stats$coverage[[i]], 2, mean),
      # precision = apply(stats$med[[i]], 2, sd),
      bias = apply(stats$bias[[i]], 2, mean)
    )
  }
  return(out)
}


#' Plot the Histogram of the Observed Rank Statistics of SBC
#'
#' Note that this plot is dependent on the number of bins, and a more general
#' visualization is to use `plot_sbc_ecdf`
#'
#' @param ranks A list of named dataframes of the rank statistic
#' @param bins An integer specifying the number of bins to use when plotting the histogram
#' @param layout Optional. A numeric vector specifying the layout using `par(mfrow = layout)`
#' @param add_stats Boolean. Should coverage, bias and precision be included in the figure.
#' @return No returns
#' @export
plot_sbc_hist <- function(ranks, bins = 10, layout = NA, add_stats = TRUE){
  if (!is.null(ranks[["rank"]])) {
    stats <- calc_sbc_stats(ranks)
  }
  if (!is.null(ranks[["rank"]])) {
    ranks <- ranks[["rank"]]
  }
  selects <- names(ranks)
  oldpar <- par(no.readonly = TRUE) # code line i
  on.exit(par(oldpar)) # code line i + 1

  n_sample <- nrow(ranks[[1]])
  low <- qbinom(0.025, n_sample, 1/bins)
  mid <- qbinom(0.5, n_sample, 1/bins)
  high <- qbinom(0.975, n_sample, 1/bins)
  par_names <- colnames(ranks[[1]])
  for (j in seq_along(ranks)) {
    if (any(is.na(layout))) {
      par(mfrow = coda_setmfrow(Nparms = ncol(ranks[[1]])))
    } else {
      par(mfrow = layout)
    }
    rank <- ranks[[j]]
    stat <- stats[[j]]
    for (i in 1:ncol(rank)) {
      hist(
        x = rank[ , i],
        main = paste0(selects[j], " - ", par_names[i]),
        breaks = bins,
        ylim = c(0, high + 2),
        xlab = "rank"
      )
      abline(h = low, lty = 2)
      abline(h = mid, lty = 2)
      abline(h = high, lty = 2)
      if (!is.null(stat) && add_stats) {
        coverage_print <- paste0("coverage : ", round(stat[["coverage"]][i], 2))
        bias_print <- paste0("bias : ", round(stat[["bias"]][i], 3))
        # precision_print <- paste0("precision : ", round(stat[["precision"]][i], 3))
        legend(x = "topleft", legend = coverage_print, bty = "n")
        legend(x = "topright", legend = bias_print, bty = "n")
        # legend(x = "topright", legend = precision_print, bty = "n")
      }
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
#' Plots  F_hat(z) - z  where F_hat is the
#' empirical CDF of the (normalized) ranks and z is the Uniform(0,1) CDF.
#' The shaded band is the simultaneous (1 - gamma) envelope for F_hat(z) - z.
#'
#' @param ranks A list of named dataframes of the rank statistic (raw or normalized)
#' @param layout Optional. A numeric vector specifying the layout using `par(mfrow = layout)`
#' @param add_stats Boolean. Should coverage, bias and precision be included in the figure.
#' @param main Optional. A character specifying plot title.
#' @param K Optional. Effective sample size of the MCMC that produced the ranks.
#' @return No returns
#' @export
plot_sbc_ecdf <- function(ranks, layout = NA, add_stats = TRUE, main = NULL, K = 500) {

  # stats extraction (keep existing behaviour)
  stats <- NULL
  if (!is.null(ranks[["rank"]])) {
    stats <- calc_sbc_stats(ranks)
    ranks <- ranks[["rank"]]
  }

  selects <- names(ranks)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  # N = number of SBC simulations (rows)
  N <- nrow(ranks[[1]])

  # Envelope for F_hat(z) - z
  gamma <- get_gamma(N, K)
  res <- get_lims(N, K, gamma)  # res$z, res$lower, res$upper

  for (j in seq_along(ranks)) {

    if (any(is.na(layout))) {
      par(mfrow = coda_setmfrow(Nparms = ncol(ranks[[j]]), nplots = length(ranks)))
    } else {
      par(mfrow = layout)
    }

    rank <- ranks[[j]]
    stat <- if (!is.null(stats)) stats[[j]] else NULL
    par_names <- colnames(rank)

    for (i in seq_len(ncol(rank))) {

      r <- rank[, i]
      r <- r[is.finite(r)]

      # Convert to normalized ranks in [0,1] if needed
      r_norm <- r
      r_min <- suppressWarnings(min(r_norm, na.rm = TRUE))
      r_max <- suppressWarnings(max(r_norm, na.rm = TRUE))

      if (!(r_max <= 1 && r_min >= 0)) {
        # If looks 1-based integer ranks: 1..K
        if (r_min >= 1 && r_max <= K) {
          r_norm <- (r_norm - 1) / K
        } else if (r_min >= 0 && r_max <= (K - 1)) {
          # If looks 0-based: 0..(K-1)
          r_norm <- r_norm / K
        } else {
          # Fallback: scale by K (keeps it in a reasonable range if ranks are 0..ESS)
          r_norm <- r_norm / K
        }
      }

      # ECDF difference: F_hat(z) - z
      Fn <- ecdf(r_norm)
      z <- res[["z"]]
      y <- Fn(z) - z

      main_i <- if (is.null(main)) paste0(selects[j], " - ", par_names[i]) else main

      ylim_lo <- min(res[["lower"]], y, na.rm = TRUE) - 0.01
      ylim_hi <- max(res[["upper"]], y, na.rm = TRUE) + 0.03

      plot(
        x = z,
        y = y,
        type = "s",
        lwd = 2,
        ylim = c(ylim_lo, ylim_hi),
        xlim = c(0, 1),
        ylab = "ECDF Difference  F(z) - z",
        xlab = "Normalized Rank Statistic (z)",
        main = main_i
      )

      polygon(
        x = c(z, rev(z)),
        y = c(res[["lower"]], rev(res[["upper"]])),
        col = adjustcolor("cornflowerblue", 0.2),
        border = NA
      )

      # redraw ECDF diff on top of ribbon
      lines(z, y, type = "s", lwd = 2)

      abline(h = 0, lty = 2, col = "gray50")

      if (!is.null(stat) && add_stats) {
        coverage_print <- paste0("coverage : ", round(stat[["coverage"]][i], 2))
        bias_print <- paste0("bias : ", round(stat[["bias"]][i], 3))
        legend(x = "topleft", legend = coverage_print, bty = "n")
        legend(x = "topright", legend = bias_print, bty = "n")
      }
    }
  }

  invisible(NULL)
}

get_ranks_ESS <- function(posterior, ESS, prior){
  posterior <- posterior[seq(1, length(posterior), length.out = ESS)]
  return(rank(c(prior, posterior))[1]/(ESS+1))
}
