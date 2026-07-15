#' Summarize posterior SSRT moments for stop-signal models
#'
#' Computes posterior summaries of the stop-time distribution mean and standard
#' deviation for EMC2 stop-signal models. For individual-level summaries, SSRT
#' moments are computed from each subject's posterior `alpha` draws. For
#' population-level summaries, implied individuals are sampled from each
#' posterior draw of the hierarchical multivariate normal distribution before
#' transforming to seconds.
#'
#' @param emc A fitted EMC2 object using a stop-signal model.
#' @param level Character. `"population"` summarizes the marginal population
#' distribution implied by the hierarchical mean and covariance draws.
#' `"individual"` summarizes each subject's posterior `alpha` draws.
#' `"condition"` maps stop parameters to design cells before computing SSRT
#' moments. For `type = "single"` fits, the default is resolved to
#' `"individual"` because no population distribution is sampled.
#' @param probs Numeric vector of quantiles to return.
#' @param details Logical. If `FALSE`, return a compact summary with only the
#' SSRT mean and standard deviation rows, requested posterior summary columns,
#' necessary subject/condition identifiers, and `n_population` where relevant.
#' The full detailed summary is retained in the `"ssrt_details"` attribute. If
#' `TRUE`, return the detailed summary directly, including technical metadata
#' and population variance-decomposition rows.
#' @param stage Sampling stage passed to [get_pars()].
#' @param n_population Number of implied population members sampled per
#' posterior draw when `level = "population"`.
#' @param seed Optional random seed used for population sampling.
#' @param data Optional data frame used to define condition cells when
#' `level = "condition"`. Defaults to the fitted data.
#' @param condition_vars Optional character vector of columns defining
#' condition cells. By default, variables appearing in stop-parameter formulas
#' are used.
#' @param condition_level Character. Only used when `level = "condition"`.
#' `"population"` samples implied individuals from each hierarchical
#' population draw and summarizes each condition. `"individual"` maps each
#' observed subject's `alpha` draws and returns subject-by-condition summaries.
#' @param thin,filter,length.out,subject Passed to [get_pars()].
#' @param ... Additional arguments passed to [get_pars()].
#'
#' @return A data frame with posterior summaries. The main rows are:
#' \describe{
#'   \item{`ssrt_mean`}{The mean stop finishing time on the seconds scale.}
#'   \item{`ssrt_sd`}{The total standard deviation of stop finishing times on
#'   the seconds scale. For population summaries this combines the average
#'   within-person stop-time variance and the between-person variance in SSRT
#'   means.}
#'   \item{`between_subject_sd`}{Returned for population summaries when
#'   `details = TRUE`. The
#'   standard deviation of individual SSRT means implied by the hierarchical
#'   population distribution.}
#'   \item{`mean_within_subject_sd`}{Returned for population summaries when
#'   `details = TRUE`.
#'   The average individual stop-time standard deviation.}
#' }
#' Output columns identify the summary target (`level`, `subject`,
#' `condition`, and condition columns when relevant), the stop distribution
#' (`stop_family`), the moment formula used (`method`), the Monte Carlo
#' population size for population summaries (`n_population`), and posterior
#' summary columns. With `details = FALSE`, technical metadata columns are
#' omitted, empty identifier columns are dropped, and the unabridged data frame
#' is stored in `attr(x, "ssrt_details")`. `mean` and `sd` are the posterior mean
#' and posterior standard deviation of the requested SSRT quantity. The remaining columns are
#' the requested quantiles from `probs`, named by their percentages, for example
#' `2.5%`, `50%`, and `97.5%`.
#' @export
ssrt_summary <- function(emc,
                         level = c("population", "individual", "condition"),
                         probs = c(0.025, 0.5, 0.975),
                         details = FALSE,
                         stage = get_last_stage(emc),
                         n_population = 1000,
                         seed = NULL,
                         data = NULL,
                         condition_vars = NULL,
                         condition_level = c("population", "individual"),
                         thin = 1,
                         filter = 0,
                         length.out = NULL,
                         subject = NULL,
                         ...) {
  level_missing <- missing(level)
  condition_level_missing <- missing(condition_level)
  level <- match.arg(level)
  condition_level <- match.arg(condition_level)
  single_level <- ssrt_is_single_level(emc)
  if (single_level && identical(level, "population")) {
    if (level_missing) {
      level <- "individual"
    } else {
      stop(
        "`level = \"population\"` is not defined for `type = \"single\"` fits. ",
        "Use `level = \"individual\"`.",
        call. = FALSE
      )
    }
  }
  if (single_level && identical(level, "condition") &&
      identical(condition_level, "population")) {
    if (condition_level_missing) {
      condition_level <- "individual"
    } else {
      stop(
        "`condition_level = \"population\"` is not defined for ",
        "`type = \"single\"` fits. Use `condition_level = \"individual\"`.",
        call. = FALSE
      )
    }
  }
  if (!is.numeric(probs) || any(probs < 0 | probs > 1)) {
    stop("`probs` must contain probabilities between 0 and 1.", call. = FALSE)
  }
  if (!is.logical(details) || length(details) != 1L || is.na(details)) {
    stop("`details` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  if (!is.numeric(n_population) || length(n_population) != 1L ||
      is.na(n_population) || n_population < 2) {
    stop("`n_population` must be a single number greater than or equal to 2.",
         call. = FALSE)
  }
  n_population <- as.integer(n_population)
  info <- ssrt_model_info(emc)
  has_stop_effect <- ssrt_has_stop_parameter_effects(emc, info)
  if (has_stop_effect && !identical(level, "condition")) {
    stop(
      "`ssrt_summary()` detected condition effects on stop parameters. ",
      "Use `level = \"condition\"` so parameters are mapped to design cells ",
      "before computing SSRT moments.",
      call. = FALSE
    )
  }

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      get(".Random.seed", envir = .GlobalEnv)
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
  }

  if (identical(level, "individual")) {
    out <- ssrt_summary_individual(
      emc = emc, info = info, probs = probs, stage = stage,
      thin = thin, filter = filter, length.out = length.out,
      subject = subject, ...
    )
    return(ssrt_format_output(out, details = details))
  }

  if (identical(level, "condition")) {
    out <- ssrt_summary_condition(
      emc = emc, info = info, probs = probs, stage = stage,
      data = data, condition_vars = condition_vars,
      condition_level = condition_level,
      n_population = n_population, thin = thin, filter = filter, length.out = length.out,
      subject = subject, ...
    )
    return(ssrt_format_output(out, details = details))
  }

  out <- ssrt_summary_population(
    emc = emc, info = info, probs = probs, stage = stage,
    n_population = n_population, thin = thin, filter = filter,
    length.out = length.out, ...
  )
  ssrt_format_output(out, details = details)
}

ssrt_format_output <- function(out, details) {
  if (isTRUE(details)) {
    return(out)
  }

  detailed_out <- out
  out <- out[out$parameter %in% c("ssrt_mean", "ssrt_sd"), , drop = FALSE]
  technical_columns <- c("level", "stop_family", "method")
  out <- out[, setdiff(names(out), technical_columns), drop = FALSE]

  empty_identifier <- vapply(out, function(x) {
    all(is.na(x)) || all(!nzchar(as.character(x)))
  }, logical(1))
  protected_columns <- c("parameter", "mean", "sd", "n_population")
  out <- out[, !empty_identifier | names(out) %in% protected_columns, drop = FALSE]
  attr(out, "ssrt_details") <- detailed_out
  out
}

ssrt_is_single_level <- function(emc) {
  is.list(emc) && length(emc) >= 1L && identical(emc[[1]]$type, "single")
}

ssrt_model_info <- function(emc) {
  design <- get_design(emc)[[1]]
  model <- design$model()
  if (!is.null(model$ss_info) && !is.null(model$ss_info$stop_family)) {
    info <- model$ss_info
  } else if (identical(model$c_name, "SSEXG") ||
             identical(model$c_name, "SSRDEX")) {
    info <- list(
      go_family = if (identical(model$c_name, "SSEXG")) {
        "exgaussian"
      } else {
        "racing_diffusion"
      },
      stop_family = "exgaussian",
      go_parameters = setdiff(
        names(model$p_types),
        c("muS", "sigmaS", "tauS", "exgS_lb", "tf", "gf")
      ),
      stop_parameters = c("muS", "sigmaS", "tauS"),
      stop_mean = "muS + tauS",
      stop_sd = "sqrt(sigmaS^2 + tauS^2)"
    )
  } else {
    stop("`ssrt_summary()` currently supports EMC2 stop-signal models only.",
         call. = FALSE)
  }
  info$exgS_lb <- if ("exgS_lb" %in% names(design$constants)) {
    unname(design$constants[["exgS_lb"]])
  } else if ("exgS_lb" %in% names(model$p_types)) {
    unname(model$p_types[["exgS_lb"]])
  } else {
    0
  }
  info
}

ssrt_required_parameters <- function(stop_family) {
  switch(stop_family,
         exgaussian = c("muS", "sigmaS", "tauS"),
         lognormal = c("meanlogS", "sdlogS"),
         weibull = c("shapeS", "scaleS"),
         stop("Unsupported stop distribution: ", stop_family, call. = FALSE))
}

ssrt_stop_formula_variables <- function(emc, info) {
  design <- get_design(emc)[[1]]
  formulas <- design$Flist
  stop_parameters <- ssrt_required_parameters(info$stop_family)

  unique(unlist(lapply(stop_parameters, function(parameter) {
    form <- formulas[[parameter]]
    if (is.null(form)) return(character(0))
    all.vars(form[[3]])
  }), use.names = FALSE))
}

ssrt_has_stop_parameter_effects <- function(emc, info) {
  design <- get_design(emc)[[1]]
  formulas <- design$Flist
  stop_parameters <- ssrt_required_parameters(info$stop_family)

  has_effect <- vapply(stop_parameters, function(parameter) {
    form <- formulas[[parameter]]
    if (is.null(form)) return(FALSE)
    terms <- stats::terms(form)
    length(attr(terms, "term.labels")) > 0L || attr(terms, "intercept") == 0L
  }, logical(1))

  any(has_effect)
}

ssrt_compute_moments <- function(pars, info) {
  pars <- as.data.frame(pars, check.names = FALSE)
  stop_family <- info$stop_family
  missing <- setdiff(ssrt_required_parameters(stop_family), names(pars))
  if (length(missing)) {
    stop("Missing stop parameter(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  if (identical(stop_family, "exgaussian")) {
    mu <- exp(pars$muS)
    sigma <- exp(pars$sigmaS)
    tau <- exp(pars$tauS)
    lb <- if ("exgS_lb" %in% names(pars)) pars$exgS_lb else info$exgS_lb
    out <- ssrt_exgaussian_lower_truncated_moments(mu, sigma, tau, lb)
    out$method <- "lower_truncated_exgaussian_analytic"
    return(out)
  }

  if (identical(stop_family, "lognormal")) {
    meanlog <- pars$meanlogS
    sdlog <- exp(pars$sdlogS)
    mean <- exp(meanlog + sdlog^2 / 2)
    sd <- sqrt((exp(sdlog^2) - 1) * exp(2 * meanlog + sdlog^2))
    return(data.frame(
      ssrt_mean = mean,
      ssrt_sd = sd,
      method = "lognormal_analytic",
      stringsAsFactors = FALSE
    ))
  }

  shape <- exp(pars$shapeS)
  scale <- exp(pars$scaleS)
  mean <- scale * gamma(1 + 1 / shape)
  second <- scale^2 * gamma(1 + 2 / shape)
  data.frame(
    ssrt_mean = mean,
    ssrt_sd = sqrt(pmax(second - mean^2, 0)),
    method = "weibull_analytic",
    stringsAsFactors = FALSE
  )
}

ssrt_exgaussian_lower_truncated_moments <- function(mu, sigma, tau, lb) {
  n <- max(length(mu), length(sigma), length(tau), length(lb))
  mu <- rep(mu, length.out = n)
  sigma <- rep(sigma, length.out = n)
  tau <- rep(tau, length.out = n)
  lb <- rep(lb, length.out = n)

  mean <- rep(NA_real_, n)
  sd <- rep(NA_real_, n)

  ok <- is.finite(mu) & is.finite(sigma) & is.finite(tau) &
    sigma > 0 & tau > 0 & !is.na(lb)

  untruncated <- ok & is.infinite(lb) & lb < 0
  mean[untruncated] <- mu[untruncated] + tau[untruncated]
  sd[untruncated] <- sqrt(sigma[untruncated]^2 + tau[untruncated]^2)

  finite_lb <- ok & is.finite(lb)
  if (any(finite_lb)) {
    idx <- which(finite_lb)
    b <- (lb[idx] - mu[idx]) / sigma[idx]
    k <- sigma[idx] / tau[idx]

    # Raw tail moments for X = Normal(mu, sigma) + Exp(mean = tau),
    # conditional on X > lb. Split over the latent normal value at lb.
    q <- stats::pnorm(b, lower.tail = FALSE)
    phi <- stats::dnorm(b)
    log_i0 <- -(lb[idx] - mu[idx]) / tau[idx] +
      0.5 * k^2 +
      stats::pnorm(b - k, log.p = TRUE)
    i0 <- ifelse(log_i0 < log(.Machine$double.xmin), 0, exp(pmin(log_i0, 0)))

    survival <- q + i0
    m1_high <- (mu[idx] + tau[idx]) * q + sigma[idx] * phi
    m2_high <- (mu[idx]^2 + 2 * mu[idx] * tau[idx] + 2 * tau[idx]^2) * q +
      2 * sigma[idx] * (mu[idx] + tau[idx]) * phi +
      sigma[idx]^2 * (q + b * phi)
    m1_low <- (lb[idx] + tau[idx]) * i0
    m2_low <- (lb[idx]^2 + 2 * lb[idx] * tau[idx] + 2 * tau[idx]^2) * i0

    valid <- survival > 0 & is.finite(survival)
    idx_valid <- idx[valid]
    raw_mean <- (m1_high[valid] + m1_low[valid]) / survival[valid]
    raw_second <- (m2_high[valid] + m2_low[valid]) / survival[valid]

    mean[idx_valid] <- raw_mean
    sd[idx_valid] <- sqrt(pmax(raw_second - raw_mean^2, 0))
  }

  data.frame(
    ssrt_mean = mean,
    ssrt_sd = sd,
    stringsAsFactors = FALSE
  )
}

ssrt_summarize_draws <- function(draws, probs) {
  rows <- lapply(c("ssrt_mean", "ssrt_sd"), function(parameter) {
    qs <- stats::quantile(draws[[parameter]], probs = probs, na.rm = TRUE)
    data.frame(
      parameter = parameter,
      mean = mean(draws[[parameter]], na.rm = TRUE),
      sd = stats::sd(draws[[parameter]], na.rm = TRUE),
      t(qs),
      row.names = NULL,
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

ssrt_summary_individual <- function(emc, info, probs, stage, thin, filter,
                                    length.out, subject, ...) {
  alpha <- get_pars(
    emc,
    selection = "alpha",
    stage = stage,
    thin = thin,
    filter = filter,
    length.out = length.out,
    map = FALSE,
    by_subject = TRUE,
    return_mcmc = FALSE,
    merge_chains = TRUE,
    remove_constants = FALSE,
    subject = subject,
    ...
  )

  if (length(dim(alpha)) != 3L) {
    stop("Expected `alpha` posterior draws with dimensions parameter x subject x iteration.",
         call. = FALSE)
  }

  subjects <- dimnames(alpha)[[2]]
  out <- lapply(seq_along(subjects), function(i) {
    pars <- t(alpha[, i, , drop = FALSE][, 1, ])
    moments <- ssrt_compute_moments(pars, info)
    summary <- ssrt_summarize_draws(moments, probs)
    cbind(
      level = "individual",
      subject = subjects[[i]],
      stop_family = info$stop_family,
      method = unique(moments$method),
      summary,
      row.names = NULL
    )
  })

  do.call(rbind, out)
}

ssrt_summary_condition <- function(emc, info, probs, stage, data,
                                   condition_vars, condition_level,
                                   n_population, thin, filter, length.out,
                                   subject, ...) {
  design <- get_design(emc)[[1]]
  data <- ssrt_condition_data_frame(emc, data)
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame when `level = \"condition\"`.",
         call. = FALSE)
  }

  if (is.null(condition_vars)) {
    condition_vars <- ssrt_stop_formula_variables(emc, info)
    condition_vars <- setdiff(condition_vars, "subjects")
  }
  stop_formula_vars <- setdiff(ssrt_stop_formula_variables(emc, info), "subjects")
  if (!length(stop_formula_vars)) {
    stop("Stop parameters are not modeled as a function of any condition variable; ",
         "condition-level SSRT summaries are not defined for this model.",
         call. = FALSE)
  }
  condition_vars <- intersect(condition_vars, names(data))
  if (!length(condition_vars)) {
    stop("Could not determine condition variables. Supply `condition_vars`.",
         call. = FALSE)
  }
  unused_condition_vars <- setdiff(condition_vars, stop_formula_vars)
  if (length(unused_condition_vars)) {
    stop("Condition variable(s) not used in stop-parameter formulas: ",
         paste(unused_condition_vars, collapse = ", "), ". ",
         "SSRT only varies by variables included in the stop-parameter formulas.",
         call. = FALSE)
  }

  if (identical(condition_level, "population")) {
    return(ssrt_summary_condition_population_hierarchical(
      emc = emc, design = design, data = data, info = info,
      probs = probs, stage = stage, condition_vars = condition_vars,
      n_population = n_population, thin = thin, filter = filter,
      length.out = length.out, ...
    ))
  }

  alpha <- get_pars(
    emc,
    selection = "alpha",
    stage = stage,
    thin = thin,
    filter = filter,
    length.out = length.out,
    map = FALSE,
    by_subject = TRUE,
    return_mcmc = FALSE,
    merge_chains = TRUE,
    remove_constants = FALSE,
    subject = subject,
    ...
  )
  if (length(dim(alpha)) != 3L) {
    stop("Expected `alpha` posterior draws with dimensions parameter x subject x iteration.",
         call. = FALSE)
  }

  mapped <- par_data_map(
    alpha,
    design = design,
    data = data,
    functions = design$Ffunctions,
    add_recalculated = FALSE,
    group_design = get_group_design(emc)
  )
  mapped_data <- mapped$data
  mapped_pars <- mapped$pars

  missing_vars <- setdiff(condition_vars, names(mapped_data))
  if (length(missing_vars)) {
    stop("Mapped data are missing condition variable(s): ",
         paste(missing_vars, collapse = ", "), call. = FALSE)
  }

  key_vars <- condition_vars
  if ("subjects" %in% names(mapped_data)) {
    key_vars <- c("subjects", key_vars)
  }
  key_vars <- unique(key_vars)
  keep <- !duplicated(mapped_data[, key_vars, drop = FALSE])
  mapped_data <- mapped_data[keep, , drop = FALSE]
  mapped_pars <- mapped_pars[keep, , , drop = FALSE]

  ssrt_summary_condition_by_subject(
    mapped_data = mapped_data, mapped_pars = mapped_pars, info = info,
    probs = probs, condition_vars = condition_vars
  )
}

ssrt_condition_data_frame <- function(emc, data = NULL) {
  if (is.null(data)) {
    data <- get_data(emc)
  }
  if (is.data.frame(data)) {
    return(data)
  }
  if (is.list(data) && length(data) == 1L && is.data.frame(data[[1]])) {
    return(data[[1]])
  }
  data
}

ssrt_condition_label <- function(data, condition_vars) {
  if (length(condition_vars) == 1L) {
    return(as.character(data[[condition_vars]]))
  }
  apply(data[, condition_vars, drop = FALSE], 1, function(x) {
    paste(paste(condition_vars, x, sep = "="), collapse = ", ")
  })
}

ssrt_population_condition_data <- function(data, condition_vars, n_population) {
  if (!"subjects" %in% names(data)) {
    stop("Condition-level population summaries require a `subjects` column.",
         call. = FALSE)
  }

  row_key <- condition_vars
  template <- data[!duplicated(data[, row_key, drop = FALSE]), , drop = FALSE]
  subjects <- sprintf("population_%04d", seq_len(n_population))
  out <- template[rep(seq_len(nrow(template)), times = n_population), ,
                  drop = FALSE]
  out$subjects <- factor(rep(subjects, each = nrow(template)), levels = subjects)
  rownames(out) <- NULL
  out
}

ssrt_population_alpha_draws <- function(emc, stage, n_population, thin, filter,
                                        length.out, ...) {
  mu <- get_pars(
    emc,
    selection = "mu",
    stage = stage,
    thin = thin,
    filter = filter,
    length.out = length.out,
    map = FALSE,
    return_mcmc = FALSE,
    merge_chains = TRUE,
    remove_constants = FALSE,
    ...
  )
  Sigma <- get_pars(
    emc,
    selection = "Sigma",
    stage = stage,
    thin = thin,
    filter = filter,
    length.out = length.out,
    map = FALSE,
    return_mcmc = FALSE,
    merge_chains = TRUE,
    remove_constants = FALSE,
    ...
  )

  if (is.null(rownames(mu))) {
    stop("Population mean draws have no parameter names.", call. = FALSE)
  }
  if (length(dim(Sigma)) != 3L || is.null(dimnames(Sigma)[[1]])) {
    stop("Population covariance draws must be a named parameter x parameter x iteration array.",
         call. = FALSE)
  }

  par_names <- dimnames(Sigma)[[1]]
  n_iter <- min(ncol(mu), dim(Sigma)[3])
  subjects <- sprintf("population_%04d", seq_len(n_population))
  alpha <- array(
    NA_real_,
    dim = c(length(par_names), n_population, n_iter),
    dimnames = list(par_names, subjects, NULL)
  )

  for (i in seq_len(n_iter)) {
    alpha[, , i] <- t(mvtnorm::rmvnorm(
      n_population,
      mean = mu[par_names, i],
      sigma = Sigma[, , i]
    ))
  }

  alpha
}

ssrt_summary_condition_population_hierarchical <- function(emc, design, data,
                                                           info, probs, stage,
                                                           condition_vars,
                                                           n_population,
                                                           thin, filter,
                                                           length.out, ...) {
  if (!is.null(get_group_design(emc))) {
    stop("`ssrt_summary(level = \"condition\", condition_level = \"population\")` ",
         "does not yet support group_design models.",
         call. = FALSE)
  }

  alpha <- ssrt_population_alpha_draws(
    emc = emc,
    stage = stage,
    n_population = n_population,
    thin = thin,
    filter = filter,
    length.out = length.out,
    ...
  )
  population_data <- ssrt_population_condition_data(
    data = data,
    condition_vars = condition_vars,
    n_population = n_population
  )

  mapped <- par_data_map(
    alpha,
    design = design,
    data = population_data,
    functions = design$Ffunctions,
    add_recalculated = FALSE,
    group_design = get_group_design(emc)
  )
  mapped_data <- mapped$data
  mapped_pars <- mapped$pars

  key_vars <- c("subjects", condition_vars)
  keep <- !duplicated(mapped_data[, key_vars, drop = FALSE])
  mapped_data <- mapped_data[keep, , drop = FALSE]
  mapped_pars <- mapped_pars[keep, , , drop = FALSE]

  ssrt_summary_condition_population(
    mapped_data = mapped_data, mapped_pars = mapped_pars, info = info,
    probs = probs, condition_vars = condition_vars,
    n_population = n_population
  )
}

ssrt_mapped_pars_frame <- function(mapped_pars, rows, iter = NULL) {
  par_names <- dimnames(mapped_pars)[[3]]
  if (is.null(iter)) {
    x <- mapped_pars[rows, , , drop = FALSE]
    out <- matrix(x, nrow = dim(mapped_pars)[2], ncol = dim(mapped_pars)[3])
  } else {
    x <- mapped_pars[rows, iter, , drop = FALSE]
    out <- matrix(x, nrow = length(rows), ncol = dim(mapped_pars)[3])
  }
  colnames(out) <- par_names
  as.data.frame(out, check.names = FALSE)
}

ssrt_summary_condition_by_subject <- function(mapped_data, mapped_pars, info,
                                              probs, condition_vars) {
  condition <- ssrt_condition_label(mapped_data, condition_vars)
  subject <- if ("subjects" %in% names(mapped_data)) {
    as.character(mapped_data$subjects)
  } else {
    rep(NA_character_, nrow(mapped_data))
  }

  out <- lapply(seq_len(nrow(mapped_data)), function(i) {
    pars <- ssrt_mapped_pars_frame(mapped_pars, rows = i)
    moments <- ssrt_compute_moments(pars, info)
    summary <- ssrt_summarize_draws(moments, probs)
    cbind(
      level = "condition",
      subject = subject[[i]],
      condition = condition[[i]],
      mapped_data[i, condition_vars, drop = FALSE],
      stop_family = info$stop_family,
      method = unique(moments$method),
      summary,
      row.names = NULL
    )
  })

  do.call(rbind, out)
}

ssrt_summary_condition_population <- function(mapped_data, mapped_pars, info,
                                              probs, condition_vars,
                                              n_population = NA_integer_) {
  condition <- ssrt_condition_label(mapped_data, condition_vars)
  condition_levels <- unique(condition)
  n_iter <- dim(mapped_pars)[2]

  out <- lapply(condition_levels, function(cond) {
    idx <- which(condition == cond)
    draws <- matrix(NA_real_, nrow = n_iter, ncol = 4,
                    dimnames = list(NULL, c("ssrt_mean", "ssrt_sd",
                                            "between_subject_sd",
                                            "mean_within_subject_sd")))
    methods <- character(n_iter)

    for (iter in seq_len(n_iter)) {
      pars <- ssrt_mapped_pars_frame(mapped_pars, rows = idx, iter = iter)
      moments <- ssrt_compute_moments(pars, info)
      within_var <- moments$ssrt_sd^2
      between_var <- if (length(idx) > 1L) stats::var(moments$ssrt_mean) else 0
      draws[iter, "ssrt_mean"] <- mean(moments$ssrt_mean, na.rm = TRUE)
      draws[iter, "ssrt_sd"] <- sqrt(mean(within_var, na.rm = TRUE) + between_var)
      draws[iter, "between_subject_sd"] <- sqrt(between_var)
      draws[iter, "mean_within_subject_sd"] <- mean(moments$ssrt_sd, na.rm = TRUE)
      methods[[iter]] <- moments$method[[1]]
    }

    summary <- ssrt_summarize_draws(as.data.frame(draws), probs)
    extras <- lapply(c("between_subject_sd", "mean_within_subject_sd"), function(parameter) {
      qs <- stats::quantile(draws[, parameter], probs = probs, na.rm = TRUE)
      data.frame(
        parameter = parameter,
        mean = mean(draws[, parameter], na.rm = TRUE),
        sd = stats::sd(draws[, parameter], na.rm = TRUE),
        t(qs),
        row.names = NULL,
        check.names = FALSE
      )
    })
    summary <- rbind(summary, do.call(rbind, extras))

    condition_data <- mapped_data[idx[[1]], condition_vars, drop = FALSE]
    cbind(
      level = "condition",
      subject = NA_character_,
      condition = cond,
      condition_data,
      stop_family = info$stop_family,
      method = unique(methods),
      n_population = n_population,
      summary,
      row.names = NULL
    )
  })

  do.call(rbind, out)
}

ssrt_summary_population <- function(emc, info, probs, stage, n_population,
                                    thin, filter, length.out, ...) {
  if (!is.null(get_group_design(emc))) {
    stop("`ssrt_summary(level = \"population\")` does not yet support ",
         "group_design models. Use `level = \"individual\"` for now.",
         call. = FALSE)
  }

  mu <- get_pars(
    emc,
    selection = "mu",
    stage = stage,
    thin = thin,
    filter = filter,
    length.out = length.out,
    map = FALSE,
    return_mcmc = FALSE,
    merge_chains = TRUE,
    remove_constants = FALSE,
    ...
  )
  Sigma <- get_pars(
    emc,
    selection = "Sigma",
    stage = stage,
    thin = thin,
    filter = filter,
    length.out = length.out,
    map = FALSE,
    return_mcmc = FALSE,
    merge_chains = TRUE,
    remove_constants = FALSE,
    ...
  )

  if (is.null(rownames(mu))) {
    stop("Population mean draws have no parameter names.", call. = FALSE)
  }
  if (length(dim(Sigma)) != 3L || is.null(dimnames(Sigma)[[1]])) {
    stop("Population covariance draws must be a named parameter x parameter x iteration array.",
         call. = FALSE)
  }

  required <- ssrt_required_parameters(info$stop_family)
  missing <- setdiff(required, intersect(rownames(mu), dimnames(Sigma)[[1]]))
  if (length(missing)) {
    stop("Population SSRT summaries require stop parameters to be directly ",
         "sampled population parameters. Missing: ",
         paste(missing, collapse = ", "), call. = FALSE)
  }

  par_names <- dimnames(Sigma)[[1]]
  n_iter <- min(ncol(mu), dim(Sigma)[3])
  draws <- matrix(NA_real_, nrow = n_iter, ncol = 4,
                  dimnames = list(NULL, c("ssrt_mean", "ssrt_sd",
                                          "between_subject_sd",
                                          "mean_within_subject_sd")))
  methods <- character(n_iter)

  for (i in seq_len(n_iter)) {
    pop_pars <- mvtnorm::rmvnorm(
      n_population,
      mean = mu[par_names, i],
      sigma = Sigma[, , i]
    )
    colnames(pop_pars) <- par_names
    moments <- ssrt_compute_moments(pop_pars, info)
    within_var <- moments$ssrt_sd^2
    between_var <- stats::var(moments$ssrt_mean)
    draws[i, "ssrt_mean"] <- mean(moments$ssrt_mean)
    draws[i, "ssrt_sd"] <- sqrt(mean(within_var) + between_var)
    draws[i, "between_subject_sd"] <- sqrt(between_var)
    draws[i, "mean_within_subject_sd"] <- mean(moments$ssrt_sd)
    methods[[i]] <- moments$method[[1]]
  }

  summary <- ssrt_summarize_draws(as.data.frame(draws), probs)
  extras <- lapply(c("between_subject_sd", "mean_within_subject_sd"), function(parameter) {
    qs <- stats::quantile(draws[, parameter], probs = probs, na.rm = TRUE)
    data.frame(
      parameter = parameter,
      mean = mean(draws[, parameter], na.rm = TRUE),
      sd = stats::sd(draws[, parameter], na.rm = TRUE),
      t(qs),
      row.names = NULL,
      check.names = FALSE
    )
  })
  summary <- rbind(summary, do.call(rbind, extras))

  cbind(
    level = "population",
    subject = NA_character_,
    stop_family = info$stop_family,
    method = unique(methods),
    n_population = n_population,
    summary,
    row.names = NULL
  )
}
