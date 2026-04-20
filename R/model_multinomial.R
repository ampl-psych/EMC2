.multinomial_softmax_probs <- function(trials, utility) {
  max_u <- ave(utility, trials, FUN = max)
  exp_u <- exp(utility - max_u)
  denom <- ave(exp_u, trials, FUN = sum)
  exp_u / denom
}

.multinomial_probit_trial_probs <- function(utility) {
  n_alt <- length(utility)

  if (n_alt == 1L) return(1)
  if (n_alt == 2L) {
    p_first <- pnorm((utility[1] - utility[2]) / sqrt(2))
    return(c(p_first, 1 - p_first))
  }
  if (n_alt > 21L) {
    stop("multinomial_probit currently supports at most 21 response levels")
  }

  sigma <- matrix(1, nrow = n_alt - 1L, ncol = n_alt - 1L)
  diag(sigma) <- 2
  probs <- numeric(n_alt)
  algorithm <- mvtnorm::Miwa()

  for (i in seq_len(n_alt)) {
    delta <- utility[i] - utility[-i]
    probs[i] <- mvtnorm::pmvnorm(
      lower = -delta,
      upper = rep(Inf, n_alt - 1L),
      sigma = sigma,
      algorithm = algorithm
    )[1]
  }

  probs[!is.finite(probs) | probs < 0] <- 0
  total <- sum(probs)
  if (total <= 0) {
    rep(1 / n_alt, n_alt)
  } else {
    probs / total
  }
}

pMULTINOMIAL_LOGIT <- function(trials, pars) {
  .multinomial_softmax_probs(trials, pars[, "utility"])
}

pMULTINOMIAL_PROBIT <- function(trials, pars) {
  probs <- numeric(nrow(pars))
  trial_groups <- .by_trial_index(trials)

  for (idx in trial_groups) {
    probs[idx] <- .multinomial_probit_trial_probs(pars[idx, "utility"])
  }

  probs
}

rMULTINOMIAL_LOGIT <- function(lR, pars, p_types = "utility") {
  if (!all(p_types %in% colnames(pars))) {
    stop("pars must have columns ", paste(p_types, collapse = " "))
  }

  nr <- length(levels(lR))
  n <- nrow(pars) / nr
  trial_id <- rep(seq_len(n), each = nr)
  p <- pMULTINOMIAL_LOGIT(trial_id, pars)

  p_mat <- matrix(p, nrow = nr, ncol = n)
  u <- runif(n)
  cp_mat <- apply(p_mat, 2, cumsum)
  chosen <- colSums(cp_mat < rep(u, each = nr)) + 1L
  chosen <- pmin(chosen, nr)

  data.frame(R = factor(levels(lR)[chosen], levels = levels(lR)))
}

rMULTINOMIAL_PROBIT <- function(lR, pars, p_types = "utility") {
  if (!all(p_types %in% colnames(pars))) {
    stop("pars must have columns ", paste(p_types, collapse = " "))
  }

  nr <- length(levels(lR))
  n <- nrow(pars) / nr
  utility_mat <- matrix(pars[, "utility"], nrow = nr, ncol = n)
  latent <- utility_mat + matrix(rnorm(length(utility_mat)), nrow = nr)
  chosen <- max.col(t(latent), ties.method = "first")

  data.frame(R = factor(levels(lR)[chosen], levels = levels(lR)))
}

.multinomial_model <- function(link = c("logit", "probit")) {
  link <- match.arg(link)

  pfun <- switch(link, logit = pMULTINOMIAL_LOGIT, probit = pMULTINOMIAL_PROBIT)
  rfun <- switch(link, logit = rMULTINOMIAL_LOGIT, probit = rMULTINOMIAL_PROBIT)

  list(
    type = "MULTINOMIAL",
    c_name = switch(link, logit = "MULTINOMIAL_LOGIT"),
    geometry = "multinomial",
    link = link,
    p_types = c("utility" = 0),
    transform = list(func = c(utility = "identity")),
    bound = list(minmax = cbind(utility = c(-Inf, Inf))),
    Ttransform = function(pars, dadm) {
      pars
    },
    rfun = function(data = NULL, pars) {
      rfun(data$lR, pars)
    },
    pfun = function(trials, pars) {
      pfun(trials, pars)
    },
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      log_likelihood_multinomial(pars = pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}

#' Multinomial Logit Response Model
#'
#' Unordered choice among response alternatives using a softmax over latent
#' utilities.
#'
#' Model parameters are:
#'    utility (unbounded).
#'
#' @return A model list with all the necessary functions to sample
#' @examples
#' dmnl <- design(
#'   Rlevels = c("left", "right"),
#'   factors = list(subjects = 1, S = c("left", "right")),
#'   formula = list(utility ~ 0 + S),
#'   matchfun = function(d) d$S == d$lR,
#'   model = multinomial_logit
#' )
#' @export
multinomial_logit <- function() {
  .multinomial_model("logit")
}

#' Multinomial Probit Response Model
#'
#' Unordered choice among response alternatives using latent Gaussian utility
#' noise. Assumes independent Gaussian shocks across alternatives,
#' which implies correlated pairwise differences in the implied comparison
#' space.
#'
#' Model parameters are:
#'    utility (unbounded).
#'
#' @return A model list with all the necessary functions to sample
#' @examples
#' dmnp <- design(
#'   Rlevels = c("left", "right"),
#'   factors = list(subjects = 1, S = c("left", "right")),
#'   formula = list(utility ~ 0 + S),
#'   matchfun = function(d) d$S == d$lR,
#'   model = multinomial_probit
#' )
#' @export
multinomial_probit <- function() {
  .multinomial_model("probit")
}
