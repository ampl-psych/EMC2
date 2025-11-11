#' Information Criteria and Marginal Likelihoods
#'
#' Returns the BPIC/DIC or marginal deviance (-2*marginal likelihood) for a list of samples objects.
#'
#' @param sList List of samples objects
#' @param stage A string. Specifies which stage the samples are to be taken from `"preburn"`, `"burn"`, `"adapt"`, or `"sample"`
#' @param filter An integer or vector. If it's an integer, iterations up until the value set by `filter` will be excluded.
#' If a vector is supplied, only the iterations in the vector will be considered.
#' @param use_best_fit Boolean, defaults to `TRUE`, uses the minimal or mean likelihood (whichever is better) in the
#' calculation, otherwise always uses the mean likelihood.
#' @param BayesFactor Boolean, defaults to `TRUE`. Include marginal likelihoods as estimated using WARP-III bridge sampling.
#' Usually takes a minute per model added to calculate
#' @param cores_for_props Integer, how many cores to use for the Bayes factor calculation, here 4 is the default for the 4 different proposal densities to evaluate, only 1, 2 and 4 are sensible.
#' @param cores_per_prop Integer, how many cores to use for the Bayes factor calculation if you have more than 4 cores available. Cores used will be cores_for_props * cores_per_prop. Best to prioritize cores_for_props being 4 or 2
#' @param print_summary Boolean (default `TRUE`), print table of results
#' @param digits Integer, significant digits in printed table for information criteria
#' @param digits_p Integer, significant digits in printed table for model weights
#' @param ... Additional, optional arguments
#'
#' @return Matrix of effective number of parameters, mean deviance, deviance of
#' mean, DIC, BPIC, Marginal Deviance (if `BayesFactor=TRUE`) and associated weights.
#' @examples \donttest{
#' compare(list(samples_LNR), cores_for_props = 1)
#' # Typically we would define a list of two (or more) different models:
#' # # Here the full model is an emc object with the hypothesized effect
#' # # The null model is an emc object without the hypothesized effect
#' # design_full <- design(data = forstmann,model=DDM,
#' #                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#' #                            constants=c(s=log(1)))
#' # # Now without a ~ E
#' # design_null <- design(data = forstmann,model=DDM,
#' #                            formula =list(v~0+S,a~1, t0~1, s~1, Z~1, sv~1, SZ~1),
#' #                            constants=c(s=log(1)))
#' #
#' # full_model <- make_emc(forstmann, design_full)
#' # full_model <- fit(full_model)
#' #
#' # null_model <- make_emc(forstmann, design_null)
#' # null_model <- fit(null_model)
#' # sList <- list(full_model, null_model)
#' # # By default emc uses 4 cores to parallelize marginal likelihood estimation across proposals
#' # # So cores_per_prop = 3 results in 12 cores used.
#' # compare(sList, cores_per_prop = 3)
#' }
#' @export
compare <- function(
    sList,
    stage = "sample",
    filter = NULL,
    use_best_fit = TRUE,
    BayesFactor = TRUE,
    cores_for_props = 4,
    cores_per_prop = 1,
    print_summary = TRUE,
    digits = 0,
    digits_p = 3,
    ...
) {
  if(is(sList, "emc")) {
    sList <- list(sList)
  }
  defaultsf <- 0
  if (is.numeric(filter)) {
    defaultsf <- filter[1]
  }
  sflist <- as.list(
    setNames(
      rep(defaultsf, length(sList)),
      names(sList)
    )
  )
  if (is.list(filter)) {
    for (i in names(filter)) {
      if (i %in% names(sflist)) {
        sflist[[i]] <- filter[[i]]
      }
    }
  }
  dots <- add_defaults(list(...), group_only = FALSE)
  ICs <- setNames(
    vector(mode = "list", length = length(sList)),
    names(sList)
  )
  for (i in seq_along(ICs)) {
    ICs[[i]] <- IC(
      sList[[i]],
      stage = stage,
      filter = sflist[[i]],
      use_best_fit = use_best_fit,
      subject = dots[["subject"]],
      print_summary = FALSE,
      group_only = dots[["group_only"]])
  }
  ICs <- data.frame(do.call(rbind, ICs))
  DICp <- model_weights(ICs[["DIC"]])
  BPICp <- model_weights(ICs[["BPIC"]])
  out <- cbind.data.frame(
    DIC = ICs[["DIC"]],
    wDIC = DICp,
    BPIC = ICs[["BPIC"]],
    wBPIC = BPICp,
    ICs[ , -c(1:2)]
  )
  if (BayesFactor) {
    MLLs <- numeric(length(sList))
    for (i in seq_along(MLLs)) {
      MLLs[i] <- run_bridge_sampling(
        sList[[i]],
        stage = stage,
        filter = sflist[[i]],
        both_splits = FALSE,
        cores_for_props = cores_for_props,
        cores_per_prop = cores_per_prop
      )
    }
    MD <- -2 * MLLs
    modelProbability <- model_weights(MD)
    out <- cbind.data.frame(
      MD = MD,
      wMD = modelProbability,
      out
    )
  }
  if (print_summary) {
    tmp <- out
    digits_p_cols <- c("wDIC", "wBPIC")
    if (BayesFactor) {
      digits_p_cols <- c("wMD", digits_p_cols)
    }
    digits_cols <- setdiff(colnames(tmp), digits_p_cols)
    tmp[ , digits_cols] <- round(tmp[ , digits_cols], digits = digits)
    tmp[ , digits_p_cols] <- round(tmp[ , digits_p_cols], digits = digits_p)
    print(tmp)
  }
  invisible(out)
}

#' @noRd
model_weights <- function(IC) {
  IC_delta <- -(IC - min(IC)) / 2
  log_denom <- log_sum_exp(IC_delta)
  result <- exp(IC_delta - log_denom)
  return(result)
}

# numerically stable version of log(sum(exp(x))), based on matrixStats::logSumExp
#' @noRd
log_sum_exp <- function(x) {
  if (length(x) == 0) return(-Inf)
  if (any(is.na(x))) return(NA_real_)
  if (any(x == Inf)) return(Inf)
  if (all(x == -Inf)) return(-Inf)
  max_x <- max(x)
  # accumulate only the terms smaller than max_x
  nomax_sum <- sum(exp(x[x != max_x] - max_x))
  # add back the "1" for the max element via log1p
  result <- max_x + log1p(nomax_sum)
  return(result)
}



#' Calculate information criteria (DIC, BPIC), effective number of parameters and
#' constituent posterior deviance (D) summaries (meanD = mean of D, Dmean = D
#' for mean of posterior parameters and minD = minimum of D).
#'
#' @param emc emc object or list of these
#' @param stage A string. Specifies which stage you want to plot.
#' @param filter An integer or vector. If it's an integer, iterations up until the value set by `filter` will be excluded.
#' If a vector is supplied, only the iterations in the vector will be considered.
#' @param use_best_fit Boolean, default TRUE use best of minD and Dmean in
#' calculation otherwise always use Dmean
#' @param print_summary Boolean (default TRUE) print table of results
#' @param digits Integer, significant digits in printed table
#' @param subject Integer or string selecting a single subject, default NULL
#' returns sums over all subjects
#' @param group_only Boolean. If `TRUE` will calculate the IC for the group-level only
#'
#' @return Table of DIC, BPIC, EffectiveN, meanD, Dmean, and minD
#' @noRd
IC <- function(
    emc,
    stage = "sample",
    filter = 0,
    use_best_fit = TRUE,
    print_summary = TRUE,
    digits = 0,
    subject = NULL,
    group_only = FALSE
) {
  # Mean log-likelihood for each subject
  ll <- get_pars(
    emc,
    stage = stage,
    filter = filter,
    selection = "LL",
    merge_chains = TRUE
  )
  minDs <- -2 * apply(ll[[1]][[1]], 2, min)
  mean_lls <- apply(ll[[1]][[1]], 2, mean)
  alpha <- get_pars(
    emc,
    selection = "alpha",
    stage = stage,
    filter = filter,
    by_subject = TRUE,
    merge_chains = TRUE
  )
  mean_pars <- lapply(
    alpha,
    function(x) {return(apply(do.call(rbind, x), 2, mean))}
  )
  # log-likelihood for each subject using their mean parameter vector
  data <- emc[[1]][["data"]]
  mean_pars_lls <- setNames(
    numeric(length(mean_pars)),
    names(mean_pars)
  )
  for (sub in names(mean_pars)) {
    mean_pars_lls[sub] <- calc_ll_manager(
      t(mean_pars[[sub]]),
      dadm = data[[sub]],
      emc[[1]][["model"]]
    )
  }
  Dmeans <- -2 * mean_pars_lls

  if (!is.null(subject)) {
    Dmeans <- Dmeans[subject[1]]
    mean_lls <- mean_lls[subject[1]]
    minDs <- minDs[subject[1]]
  } else {
    group_stats <- group_IC(
      emc,
      stage = stage,
      filter = filter,
      type = emc[[1]][["type"]]
    )
    if (group_only) {
      mean_lls <- group_stats[["mean_ll"]]
      minDs <- group_stats[["minD"]]
      Dmeans <- group_stats[["Dmean"]]
    } else {
      mean_lls <- c(mean_lls, group_stats[["mean_ll"]])
      minDs <- c(minDs, group_stats[["minD"]])
      Dmeans <- c(Dmeans, group_stats[["Dmean"]])
    }
  }
  if (use_best_fit) {
    minDs <- pmin(minDs, Dmeans)
  }

  # mean deviance(-2*ll of all data)
  mD <- sum(-2 * mean_lls)
  # Deviance of mean
  Dmean <- sum(Dmeans)
  # mimimum Deviance
  minD <- sum(minDs)

  # Use deviance of mean as best fit or use actual best fit
  if (!use_best_fit) {
    Dm <- Dmean
  } else {
    Dm <- minD
  }

  # effective number of parameters
  pD <- mD - Dm
  # DIC = mean deviance + effective number of parameters
  DIC <- mD + pD
  # BPIC = mean deviance + 2*effective number of parameters
  # Note this is the "easy" BPIC, instead of the complex 2007 one
  BPIC <- mD + 2 * pD

  out <- c(
    DIC = DIC,
    BPIC = BPIC,
    EffectiveN = pD,
    meanD = mD,
    Dmean = Dmean,
    minD = minD
  )
  names(out) <- c("DIC", "BPIC", "EffectiveN", "meanD", "Dmean", "minD")
  if (print_summary) {
    print(round(out, digits))
  }
  invisible(out)
}


#' Bayes Factors
#'
#' returns the Bayes Factor for two models
#'
#' @param MLL1 Numeric. Marginal likelihood of model 1. Obtained with `run_bridge_sampling()`
#' @param MLL2 Numeric. Marginal likelihood of model 2. Obtained with `run_bridge_sampling()`
#'
#' @return The BayesFactor for model 1 over model 2
#' @examples \donttest{
#' # Normally one would compare two different models
#' # Here we use two times the same model:
#' M1 <- M0 <- run_bridge_sampling(samples_LNR, both_splits = FALSE, cores_for_props = 1)
#' get_BayesFactor(M1, M0)
#' }
#' @export
get_BayesFactor <- function(MLL1, MLL2) {
  exp(MLL1 - MLL2)
}


#' Information Criteria For Each Participant
#'
#' Returns the BPIC/DIC based model weights for each participant in a list of samples objects
#'
#' @param sList List of samples objects
#' @param stage A string. Specifies which stage the samples are to be taken from `"preburn"`, `"burn"`, `"adapt"`, or `"sample"`
#' @param filter An integer or vector. If it's an integer, iterations up until the value set by `filter` will be excluded.
#' If a vector is supplied, only the iterations in the vector will be considered.
#' @param use_best_fit Boolean, defaults to `TRUE`, use minimal likelihood or mean likelihood
#' (whichever is better) in the calculation, otherwise always uses the mean likelihood.
#' @param print_summary Boolean (defaults to `TRUE`) print table of results
#' @param digits Integer, significant digits in printed table
#'
#' @return List of matrices for each subject of effective number of parameters,
#' mean deviance, deviance of mean, DIC, BPIC and associated weights.
#' @examples
#' # For a broader illustration see `compare`.
#' # Here we just take two times the same model, but normally one would compare
#' # different models
#' compare_subject(list(m0 = samples_LNR, m1 = samples_LNR))
#' @export
compare_subject <- function(
    sList,
    stage = "sample",
    filter = 0,
    use_best_fit = TRUE,
    print_summary = TRUE,
    digits = 3
) {
  if (is(sList, "emc")) {
    sList <- list(sList)
  }
  subjects <- names(sList[[1]][[1]][["data"]])
  is_single <- sapply(
    sList,
    function(x) {return(x[[1]][["type"]] == "single")}
  )
  if (any(!is_single)) {
    warning("subject-by-subject comparison is best done with models of type `single`")
  }
  out <- setNames(
    vector(mode = "list", length = length(subjects)),
    subjects
  )
  for (i in subjects) {
    out[[i]] <- compare(
      sList,
      subject = i,
      BayesFactor = FALSE,
      stage = stage,
      filter = filter,
      use_best_fit = use_best_fit,
      print_summary = FALSE
    )
  }
  if (print_summary) {
    pDIC <- curate_model_weights(out, "wDIC")
    pBPIC <- curate_model_weights(out, "wBPIC")
    print(
      round(cbind(pDIC, pBPIC), digits = digits)
    )
    mnams <- unlist(
      lapply(
        strsplit(dimnames(pDIC)[[2]], "_"),
        function(x) {return(x[[2]])}
      )
    )
    cat("\nWinners\n")
    print(
      rbind(
        DIC = table(mnams[apply(pDIC, 1, which.max)]),
        BPIC = table(mnams[apply(pBPIC, 1, which.max)])
      )
    )
  }
  invisible(out)
}

#' @noRd
curate_model_weights <- function(out, type = c("wDIC", "wBPIC")) {
  result <- lapply(out, function(x) {return(x[type])})
  result <- lapply(
    result,
    function(x) {
      return(
        setNames(
          data.frame(t(x)),
          paste(type, rownames(x), sep = "_")
        )
      )
    }
  )
  result <- do.call(rbind, result)
  return(result)
}


# #' Calculate a table of model probabilities based for a list of samples objects
# #' based on samples of marginal log-likelihood (MLL) added to these objects by
# #' run_IS2. Probabilities estimated by a bootstrap ath picks a vector of MLLs,
# #' one for each model in the list randomly with replacement nboot times,
# #' calculates model probabilities and averages
# #'
# #'
# #' @param mll List of samples objects with IS_samples attribute added by by run_IS2
# #' @param nboot Integer number of bootstrap samples, the default (1e5) usually
# #' gives stable results at 2 decimal places.
# #' @param print_summary Boolean (default TRUE) print table of results
# #' @param digits Integer, significant digits in printed table
# #'
# #' @return Vector of model probabilities with names from samples list.

# mll is a list of vectors of marginal log-likelihoods for a set of models
# picks a vector of mlls for each model in the list randomly with replacement
# nboot times, calculates model probabilities and averages.
compare_MLL <- function(
    mll,
    nboot = 1e5,
    digits = 2,
    print_summary = TRUE
) {
  mll_samples <- lapply(
    mll,
    function(x) {
      IS_samples <- attr(x, "IS_samples")
      result <- IS_samples[sample(length(IS_samples), nboot, replace = TRUE)]
      return(result)
    }
  )
  mll_samples <- do.call(rbind, mll_samples)
  model_probs <- apply(
    mll_samples, 2, pmp
  )
  mean_model_probs <- apply(
    model_probs, 1, mean
  )
  out <- sort(mean_model_probs, decreasing = TRUE)
  print(round(out, digits = digits))
  invisible(out)
}

# posterior model probability for a vector of marginal log-likelihoods
#' @noRd
pmp <- function(x) {
  log_denom <- log_sum_exp(x)
  result <- exp(x - log_denom)
  return(result)
}


#' Model Averaging
#'
#' Computes model weights and a Bayes factor by comparing two groups of models based on their
#' Information Criterion (IC) values. The function works with either numeric vectors or data
#' frames containing multiple IC measures (e.g., MD, BPIC, DIC).
#'
#' When provided with numeric vectors, it computes the weights for the two groups by first
#' converting the IC values into relative weights and then normalizing them. When provided with
#' a data frame, it assumes that the data frame is the output of a call to `compare`
#' and applies averaging to each IC metric
#'
#' @param IC_for A numeric vector or the output of `compare`
#' @param IC_against A numeric vector or the output of `compare`
#'
#' @return A \code{data.frame} with the following columns:
#'   \describe{
#'     \item{\code{wFor}}{The aggregated weight of the models in favor.}
#'     \item{\code{wAgainst}}{The aggregated weight of the models against.}
#'     \item{\code{Factor}}{The Bayes factor (ratio of \code{wFor} to \code{wAgainst}).}
#'   }
#'   If \code{IC_for} is a data frame, a matrix with rows corresponding to each IC measure is returned.
#'
#' @examples
#' # First set up some example models (normally these would be alternative models)
#' samples_LNR2 <- subset(samples_LNR, length.out = 45)
#' samples_LNR3 <- subset(samples_LNR, length.out = 40)
#' samples_LNR4 <- subset(samples_LNR, length.out = 35)
#'
#' # Run compare on them, BayesFactor = F is set for speed.
#' ICs <- compare(list(S1 = samples_LNR, S2 = samples_LNR2,
#'                     S3 = samples_LNR3, S4 = samples_LNR4), BayesFactor = FALSE)
#'
#' # Model averaging can either be done with a vector of ICs:
#' model_averaging(ICs$BPIC[1:2], ICs$BPIC[2:4])
#'
#' # Or the output of compare:
#' model_averaging(ICs[1:2,], ICs[3:4,])
#'
#' @export
model_averaging <- function(IC_for, IC_against) {
  if(is.null(IC_for)) {
    return(NULL)
  }

  # case: data frame input (from `compare`) - extract vectors, then recursive call
  if(is.data.frame(IC_for)) {
    MD <- model_averaging(IC_for[["MD"]], IC_against[["MD"]])
    BPIC <- model_averaging(IC_for[["BPIC"]], IC_against[["BPIC"]])
    DIC <- model_averaging(IC_for[["DIC"]], IC_against[["DIC"]])
    result <- rbind(MD = MD, BPIC = BPIC, DIC = DIC)
    return(result)
  }

  # case: numeric vector input
  all_IC <- c(IC_for, IC_against)
  logw <- -0.5 * (all_IC - min(all_IC))
  log_denom <- log_sum_exp(logw)
  weights <- exp(logw - log_denom)

  weight_for <- sum(weights[seq_along(IC_for)])
  weight_against <- sum(weights[(length(IC_for) + 1):length(all_IC)])
  bayes_factor <- weight_for / weight_against

  result <- data.frame(
    wFor = weight_for,
    wAgainst = weight_against,
    Factor = bayes_factor
  )
  return(result)
}


# -----------------------------------------------------------------------------
# Summary statistic functions used internally ---------------------------------

get_summary_stat <- function(
    emc,
    selection = "mu",
    fun,
    stat = NULL,
    stat_only = FALSE,
    stat_name = NULL,
    digits = 3,
    ...
) {
  dots <- list(...)
  if (
    is.null(emc[[1]][["n_subjects"]]) || length(dots[["subject"]]) == 1 ||
    emc[[1]][["n_subjects"]] == 1
  ) {
    dots[["by_subject"]] <- TRUE
  }
  MCMC_samples <- do.call(
    get_pars,
    c(list(emc = emc, selection = selection), fix_dots(dots, get_pars))
  )
  out <- vector("list", length = length(MCMC_samples))
  for (i in seq_along(MCMC_samples)) {
    if (length(fun) > 1) {
      outputs <- vector("list", length = length(fun))
      for (j in seq_along(fun)) {
        outputs[[j]] <- do.call(
          fun[[j]],
          c(list(MCMC_samples[[i]]), fix_dots(dots, fun[[j]]))
        )
      }
      out[[i]] <- do.call(cbind, outputs)
      if (!is.null(stat_name)) {
        if (ncol(out[[i]]) != length(stat_name)) {
          stop("`stat_name` must be the same length as function output")
        }
        colnames(out[[i]]) <- stat_name
      }
    } else {
      out[[i]] <- do.call(
        fun,
        c(list(MCMC_samples[[i]]), fix_dots(dots, fun))
      )
    }
  }
  names(out) <- names(MCMC_samples)
  if (length(fun) == 1 & !is.matrix(out[[i]]) & !is.null(stat)) {
    out <- make_nice_summary(
      object = out, stat = stat, stat_only = stat_only, stat_name = stat_name
    )
    out <- round(out, digits)
  } else {
    out <- lapply(out, round, digits)
  }
  return(out)
}


make_nice_summary <- function(
    object,
    stat = "max",
    stat_only = FALSE,
    stat_name = NULL,
    ...
) {
  if (is.null(stat_name)) {
    stat_name <- stat
  }
  row_names <- names(object)
  col_names <- unique(unlist(lapply(object, names)))
  if (all(row_names %in% col_names)) {
    col_names <- row_names
  }
  out_mat <- matrix(NA, nrow = length(row_names), ncol = length(col_names))
  for (i in seq_along(object)) {
    idx <- col_names %in% names(object[[i]])
    out_mat[i, idx] <- object[[i]]
  }
  row_stat <- apply(
    out_mat, 1, get(stat), na.rm = TRUE
  )
  out_mat <- cbind(out_mat, row_stat)

  if (nrow(out_mat) > 1) {
    col_stat <- apply(
      out_mat, 2, get(stat), na.rm = TRUE
    )
    col_stat[length(col_stat)] <- get(stat)(unlist(object))
    out_mat <- rbind(out_mat, c(col_stat))
    rownames(out_mat) <- c(row_names, stat_name)
  } else {
    rownames(out_mat) <- row_names
  }
  colnames(out_mat) <- c(col_names, stat_name)
  if (stat_only) {
    out_mat <- out_mat[nrow(out_mat), ncol(out_mat)]
  }
  return(out_mat)
}


get_posterior_quantiles <- function(x, probs = c(0.025, .5, .975)) {
  result <- summary(x, probs)[["quantiles"]]
  return(result)
}


# -----------------------------------------------------------------------------
# Various distribution functions used internally ------------------------------

std_error_IS2 <- function(IS_samples, n_bootstrap = 5e4) {
  log_marglik_boot <- array(dim = n_bootstrap)
  for (i in 1:n_bootstrap) {
    #resample with replacement from the lw
    log_weight_boot <- sample(IS_samples, length(IS_samples), replace = TRUE)
    log_marglik_boot[i] <- stats::median(log_weight_boot)
  }
  result <- stats::sd(log_marglik_boot)
  return(result)
}

robust_diwish <- function (W, v, S) {
  if (!is.matrix(S)) {
    S <- matrix(S)
  }
  if (!is.matrix(W)) {
    W <- matrix(W)
  }
  p <- nrow(S)
  gammapart <- sum(lgamma((v + 1 - 1:p) / 2))
  ldenom <- gammapart + 0.5 * v * p * log(2) + 0.25 * p * (p - 1) * log(pi)
  cholW <- base::chol(Matrix::nearPD(W)[["mat"]])
  cholS <- base::chol(Matrix::nearPD(S)[["mat"]])
  halflogdetS <- sum(log(diag(cholS)))
  halflogdetW <- sum(log(diag(cholW)))
  invW <- chol2inv(cholW)
  exptrace <- sum(S * invW)
  lnum <- v * halflogdetS - (v + p + 1) * halflogdetW - 0.5 * exptrace
  lpdf <- lnum - ldenom
  out <- exp(lpdf)
  if (!is.finite(out)) {
    return(1e-100)
  }
  return(out)
}

dhalft <- function (x, scale = 25, nu = 1, log = FALSE) {
  x <- as.vector(x)
  scale <- as.vector(scale)
  nu <- as.vector(nu)
  if (any(scale <= 0)) {
    stop("The scale parameter must be positive.")
  }
  NN <- max(length(x), length(scale), length(nu))
  x <- rep(x, len = NN)
  scale <- rep(scale, len = NN)
  nu <- rep(nu, len = NN)
  dens <- log(2) - log(scale) + lgamma((nu + 1)/2) - lgamma(nu/2) -
    0.5 * log(pi * nu) - (nu + 1)/2 * log(1 + (1/nu) * (x/scale) * (x/scale))
  if (!log) {
    return(exp(dens))
  }
  return(dens)
}

rwish <- function(v, S) {
  if (!is.matrix(S)) {
    S <- matrix(S)
  }
  if (nrow(S) != ncol(S)) {
    stop("S not square in rwish().\n")
  }
  if (v < nrow(S)) {
    stop("v is less than the dimension of S in rwish().\n")
  }
  p <- nrow(S)
  CC <- chol(S)
  Z <- matrix(0, p, p)
  diag(Z) <- sqrt(stats::rchisq(p, v:(v - p + 1)))
  if (p > 1) {
    pseq <- 1:(p - 1)
    Z_idx <- rep(p * pseq, pseq) + unlist(lapply(pseq, seq))
    Z[Z_idx] <- stats::rnorm(p * (p - 1)/2)
  }
  result <- crossprod(Z %*% CC)
  return(result)
}

riwish <- function(v, S){
  result <- solve(rwish(v, solve(S)))
  return(result)
}

logdinvGamma <- function(x, shape, rate){
  result <- stats::dgamma(1/x, shape, rate, log = TRUE) - 2 * log(x)
  return(result)
}

condMVN <- function(
    mean, sigma, dependent.ind, given.ind, X.given, check.sigma = TRUE
) {
  if (missing(dependent.ind)) {
    stop("You must specify the indices of dependent random variables in `dependent.ind'")
  }
  if (missing(given.ind) & missing(X.given)) {
    result <- list(
      condMean = mean[dependent.ind],
      condVar = as.matrix(sigma[dependent.ind, dependent.ind])
    )
    return(result)
  }
  if (length(given.ind) == 0) {
    result <- list(
        condMean = mean[dependent.ind],
        condVar = as.matrix(sigma[dependent.ind, dependent.ind])
    )
    return(result)
  }
  if (length(X.given) != length(given.ind)) {
    stop("lengths of `X.given' and `given.ind' must be same")
  }
  if (check.sigma) {
    if (!isSymmetric(sigma)) {
      stop("sigma is not a symmetric matrix")
    }
    eigenvalues <- eigen(sigma, only.values = TRUE)[["values"]]
    if (any(eigenvalues < 1e-08)) {
      sigma <- sigma + abs(diag(stats::rnorm(nrow(sigma), sd = 1e-3)))
    }
  }
  B <- sigma[dependent.ind, dependent.ind]
  C <- sigma[dependent.ind, given.ind, drop = FALSE]
  D <- sigma[given.ind, given.ind]
  CDinv <- C %*% chol2inv(chol(D))
  cMu <- c(mean[dependent.ind] + CDinv %*% (X.given - mean[given.ind]))
  cVar <- B - CDinv %*% t(C)
  result <- list(condMean = cMu, condVar = cVar)
  return(result)
}

