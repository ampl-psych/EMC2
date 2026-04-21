.ordered_interval_prob <- function(lower, upper, pars, cdf_fun) {
  cdf_fun(upper, pars) - cdf_fun(lower, pars)
}

.ordered_cut_transform <- function(raw_cut, lR) {
  nr <- length(levels(lR))
  if (length(raw_cut) %% nr != 0) {
    stop("cut vector length must be divisible by the number of response levels")
  }

  cut <- matrix(raw_cut, nrow = nr)
  if (nr == 2L) {
    cut[2, ] <- cut[1, ]
    return(as.vector(cut))
  }

  if (nr > 2L) {
    increments <- apply(exp(cut[2:(nr - 1), , drop = FALSE]), 2, cumsum)
    if (!is.matrix(increments)) {
      increments <- matrix(increments, nrow = 1L)
    }
    cut[2:(nr - 1), ] <- sweep(increments, 2, cut[1, ], "+")
  }
  cut[nr, ] <- cut[nr - 1, ]
  as.vector(cut)
}

.ordered_cut_formula <- function(form, da) {
  lhs <- as.character(stats::terms(form)[[2]])
  rhs_terms <- attr(stats::terms(form), "term.labels")
  has_intercept <- attr(stats::terms(form), "intercept") != 0
  if (length(rhs_terms) == 0 && has_intercept && nlevels(da$lR) > 1L) {
    new_form <- stats::as.formula(paste0(lhs, " ~ 0 + lR"), env = environment(form))
    attr(new_form, "Clist") <- attr(form, "Clist")
    return(new_form)
  }
  form
}

.ordered_cut_data <- function(da) {
  lR_levels <- levels(da$lR)
  if (length(lR_levels) < 2) {
    stop("Ordered models require at least two response levels.")
  }

  cut_levels <- lR_levels[-length(lR_levels)]
  lR_values <- as.character(da$lR)
  lR_values[lR_values == lR_levels[length(lR_levels)]] <- cut_levels[length(cut_levels)]
  da$lR <- factor(lR_values, levels = cut_levels)
  da
}

.ordered_prepare_dm <- function(p_name, form, da) {
  if (!identical(p_name, "cut")) {
    return(list(form = form, da = da))
  }

  da <- .ordered_cut_data(da)
  form <- .ordered_cut_formula(form, da)
  list(form = form, da = da)
}

.ordered_rng <- function(lR, pars, latent_rng,
                         p_types = c("location", "scale", "cut_expanded"), lower = -Inf) {
  if (!all(p_types %in% colnames(pars))) {
    stop("pars must have columns ", paste(p_types, collapse = " "))
  }

  nr <- length(levels(lR))
  n <- nrow(pars) / nr
  first <- seq.int(1L, by = nr, length.out = n)
  cut <- matrix(pars[, "cut_expanded"], nrow = nr)
  cut[nrow(cut), ] <- Inf

  bins <- rbind(
    latent_rng(pars[first, , drop = FALSE]),
    rep(lower, ncol(cut)),
    cut
  )
  bins[nrow(bins), ] <- Inf

  R <- factor(
    apply(bins, 2, function(x) .bincode(x[1], x[-1])),
    levels = seq_len(length(levels(lR))),
    labels = levels(lR)
  )

  data.frame(R = R)
}

.ordered_probit_cdf <- function(x, pars) {
  pnorm(x, mean = pars[, "location"], sd = pars[, "scale"])
}

.ordered_logit_cdf <- function(x, pars) {
  stats::plogis(x, location = pars[, "location"], scale = pars[, "scale"])
}

.ordered_probit_rng <- function(pars) {
  rnorm(nrow(pars), pars[, "location"], pars[, "scale"])
}

.ordered_logit_rng <- function(pars) {
  eps <- .Machine$double.eps
  u <- stats::runif(nrow(pars), min = eps, max = 1 - eps)
  pars[, "location"] + pars[, "scale"] * stats::qlogis(u)
}

.ordered_model <- function(link = c("probit", "logit")) {
  link <- match.arg(link)

  cdf_fun <- switch(link, probit = .ordered_probit_cdf, logit = .ordered_logit_cdf)
  rng_fun <- switch(link, probit = .ordered_probit_rng, logit = .ordered_logit_rng)
  qfun <- switch(link, probit = qnorm, logit = stats::qlogis)

  list(
    type = "ORDERED",
    c_name = switch(link, probit = "ORDERED_PROBIT", logit = "ORDERED_LOGIT"),
    geometry = "ordered",
    link = link,
    p_types = c("location" = 0, "scale" = log(1), "cut" = 0),
    sampled_by_default = "cut",
    transform = list(func = c(location = "identity", scale = "exp", cut = "identity")),
    bound = list(minmax = cbind(location = c(-Inf, Inf), scale = c(0, Inf), cut = c(-Inf, Inf))),
    prepare_dm = .ordered_prepare_dm,
    Ttransform = function(pars, dadm) {
      cut_expanded <- .ordered_cut_transform(pars[, "cut"], dadm$lR)
      pars <- cbind(pars, cut_expanded = cut_expanded)
      pars
    },
    rfun = function(data = NULL, pars) {
      .ordered_rng(data$lR, pars, rng_fun)
    },
    pfun = function(lower, upper, pars) {
      .ordered_interval_prob(lower, upper, pars, cdf_fun)
    },
    qfun = function(p) {
      qfun(p)
    },
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      log_likelihood_ordered(pars = pars, dadm = dadm, model = model, min_ll = min_ll, lower = -Inf)
    }
  )
}

#' Ordered Probit Response Model
#'
#' Ordered responses based on a latent Gaussian evidence variable. Binary
#' response models are the 2-category special case.
#'
#' Model parameters are:
#'    location (unbounded),
#'    scale (log scale), and
#'    cut (first threshold free, remaining thresholds modeled as positive increments).
#'
#' The `cut ~ 1` specification yields the standard flexible `K - 1` threshold
#' parameterization from the ordinal regression literature. The final response
#' category has an implicit upper cut at `Inf`. Internally, `Ttransform` leaves
#' sampled `cut` values unchanged and returns
#' the ordered thresholds in a derived `cut_expanded` column used by the
#' likelihood and random-number generator.
#'
#' @return A model list with all the necessary functions to sample
#' @examples
#' dord <- design(
#'   Rlevels = c("left", "right"),
#'   factors = list(subjects = 1, S = c("left", "right")),
#'   formula = list(location ~ 0 + S, scale ~ 1, cut ~ 1),
#'   matchfun = function(d) d$S == d$lR,
#'   constants = c(scale = log(1)),
#'   model = ordered_probit
#' )
#' @export
ordered_probit <- function() {
  .ordered_model("probit")
}

#' Ordered Logit Response Model
#'
#' Ordered responses based on a latent logistic evidence variable. Binary
#' response models are the 2-category special case.
#'
#' The `scale` parameter is the latent logistic scale, and thresholds follow the
#' same ordered `K - 1` parameterization as `ordered_probit()`, with derived
#' ordered thresholds returned in `cut_expanded`.
#'
#' @return A model list with all the necessary functions to sample
#' @examples
#' dord <- design(
#'   Rlevels = c("left", "right"),
#'   factors = list(subjects = 1, S = c("left", "right")),
#'   formula = list(location ~ 0 + S, scale ~ 1, cut ~ 1),
#'   matchfun = function(d) d$S == d$lR,
#'   constants = c(scale = log(1)),
#'   model = ordered_logit
#' )
#' @export
ordered_logit <- function() {
  .ordered_model("logit")
}
