# Correlated multiple-threshold log-normal race (MTLNR; Reynolds, Kvam, Osth &
# Heathcote, 2020): two racing log-normal accumulators whose log finishing times
# are bivariate normal (correlation rho), plus a rating read from the losing
# accumulator's evidence at the decision (rating layer, R/model_rating.R).
#
# With winner w at decision time dt = rt - t0 and loser l, the loser's log
# finishing time given the winner's is normal with
#   m_l|w = m_l + rho * (s_l / s_w) * (log(dt) - m_w),  s_l|w = s_l * sqrt(1 - rho^2).
# The LNR is deterministic, so "loser's evidence below proportion d at dt" is
# "loser finishes after dt / d", and
#   L(R = w, RR = c, rt) = dlnorm(dt; m_w, s_w) *
#       [ S_l|w(dt / d_{K-c+1}) - S_l|w(dt / d_{K-c}) ],   d_0 = 0, d_K = 1.
# Everything is computed on the log scale.

# Log density for vectors over trials: winner (mw, sw), loser (ml, sl), trial
# correlation rho and the loser's rating interval (dlo, dhi) (rating_interval).
mtlnr_logdens <- function(dt, mw, sw, ml, sl, rho, dlo, dhi) {
  out <- rep(-Inf, length(dt))
  ok <- !is.na(dt) & dt > 0 & is.finite(dt)
  if (!any(ok)) return(out)
  ldt <- log(dt[ok])
  sc <- sl[ok] * sqrt((1 - rho[ok]) * (1 + rho[ok]))
  mc <- ml[ok] + rho[ok] * (sl[ok] / sw[ok]) * (ldt - mw[ok])
  z_hi <- (ldt - log(dhi[ok]) - mc) / sc
  z_lo <- (ldt - log(dlo[ok]) - mc) / sc
  lb <- log_pnorm_interval(z_hi, z_lo)
  lb[is.na(lb) | !(sc > 0)] <- -Inf
  out[ok] <- stats::dlnorm(dt[ok], mw[ok], sw[ok], log = TRUE) + lb
  out
}

# Split a two-rows-per-trial parameter matrix (rows ordered by trial, then
# accumulator) into winner and loser rows. winner: logical per row, or an
# integer 1/2 per trial giving the winning accumulator.
mtlnr_rows <- function(pars, winner) {
  n <- nrow(pars) / 2
  first <- seq(1, 2 * n, by = 2)
  if (is.logical(winner)) {
    wi <- matrix(winner, nrow = 2)[1, ]  # TRUE when accumulator 1 wins
    w_idx <- ifelse(wi, first, first + 1)
  } else {
    w_idx <- first + as.integer(winner) - 1L
  }
  l_idx <- ifelse(w_idx == first, first + 1, first)
  list(w = w_idx, l = l_idx, first = first)
}

# Density of the MTLNR. rt, R (factor or integer 1/2) and RR (1 = lowest) are
# per trial; pars has two rows per trial (accumulators in level order) with
# columns m, s, t0, rho and the natural-scale thresholds d1 .. d{K-1}. rho is
# read from the trial's first row, t0 from the winner's, thresholds from the
# loser's.
dMTLNR <- function(rt, R, RR, pars, n_ratings = 3, log = FALSE) {
  K <- rating_check_n(n_ratings)
  ix <- mtlnr_rows(pars, as.integer(R))
  pw <- pars[ix$w, , drop = FALSE]
  pl <- pars[ix$l, , drop = FALSE]
  iv <- rating_interval(pl, RR, K)
  ld <- mtlnr_logdens(rt - pw[, "t0"], pw[, "m"], pw[, "s"], pl[, "m"], pl[, "s"],
                      pars[ix$first, "rho"], iv$lower, iv$upper)
  if (log) ld else exp(ld)
}

# Simulate from the MTLNR: lR the accumulator factor (two rows per trial),
# pars as for dMTLNR. Returns R, rt and RR.
rMTLNR <- function(lR, pars, n_ratings = 3, ok = rep(TRUE, nrow(pars))) {
  K <- rating_check_n(n_ratings)
  if (length(levels(lR)) != 2) stop("MTLNR simulates two accumulators only")
  if (is.null(ok)) ok <- rep(TRUE, nrow(pars))
  n <- nrow(pars) / 2
  first <- seq(1, 2 * n, by = 2)
  p1 <- pars[first, , drop = FALSE]
  p2 <- pars[first + 1, , drop = FALSE]
  rho <- p1[, "rho"]
  z1 <- stats::rnorm(n)
  z2 <- stats::rnorm(n)
  T1 <- exp(p1[, "m"] + p1[, "s"] * z1)
  T2 <- exp(p2[, "m"] + p2[, "s"] * (rho * z1 + sqrt((1 - rho) * (1 + rho)) * z2))
  win1 <- T1 < T2
  dt <- ifelse(win1, T1, T2)
  Tl <- ifelse(win1, T2, T1)
  pl <- p2
  pl[!win1, ] <- p1[!win1, ]
  RR <- rating_from_evidence(dt / Tl, pl, K)
  rt <- ifelse(win1, p1[, "t0"], p2[, "t0"]) + dt
  R <- ifelse(win1, 1L, 2L)
  tok <- matrix(ok, nrow = 2)
  tok <- tok[1, ] & tok[2, ]
  rt[!tok] <- NA
  RR[!tok] <- NA
  R[!tok] <- NA
  data.frame(R = factor(levels(lR)[R], levels = levels(lR)), rt = rt, RR = RR)
}

# Trial log-likelihoods over a (possibly compressed) dadm with two rows per
# trial, expanded to the uncompressed trials and floored at min_ll; rows of
# pars align with rows of dadm. The R reference for the C++ kernel
# (src/model_mtlnr.h, c_name "MTLNR"); used when c_name is removed.
trial_ll_mtlnr <- function(pars, dadm, K, min_ll = log(1e-10)) {
  ix <- mtlnr_rows(pars, dadm$winner)
  pw <- pars[ix$w, , drop = FALSE]
  pl <- pars[ix$l, , drop = FALSE]
  RR <- dadm$RR[ix$w]
  iv <- rating_interval(pl, RR, K)
  ll <- mtlnr_logdens(dadm$rt[ix$w] - pw[, "t0"], pw[, "m"], pw[, "s"], pl[, "m"], pl[, "s"],
                      pars[ix$first, "rho"], iv$lower, iv$upper)
  # a parameter out of bounds on either accumulator row voids the trial (as
  # the C++ likelihood's apply_bounds does)
  ok <- attr(pars, "ok")
  if (!is.null(ok)) ll[!(ok[ix$w] & ok[ix$l])] <- min_ll
  ll[is.na(ll)] <- min_ll
  pmax(min_ll, ll[attr(dadm, "expand")])
}

log_likelihood_mtlnr <- function(pars, dadm, K, min_ll = log(1e-10))
  sum(trial_ll_mtlnr(pars, dadm, K, min_ll))

# Refuses censoring, truncation and missing responses (the likelihood
# is not yet defined for them).
mtlnr_check_data <- function(data, K) {
  if (length(levels(data$R)) != 2) stop("MTLNR needs exactly two response levels")
  if (any(is.na(data$R)) || any(!is.finite(data$rt)))
    stop("MTLNR does not yet support missing responses or censored rts")
  if (("LT" %in% names(data) && any(data$LT > 0, na.rm = TRUE)) ||
      ("UT" %in% names(data) && any(is.finite(data$UT))))
    stop("MTLNR does not yet support truncation")
  rating_check_data(data, K)
}

#' The Correlated Multiple-Threshold Log-Normal Race Model
#'
#' Model file for the correlated multiple-threshold log-normal race (MTLNR) of
#' Reynolds, Kvam, Osth and Heathcote (2020), a model of choice, response time
#' and a rating (e.g. confidence) given after the choice.
#'
#' The data need a factor `R` with two levels, `rt`, and a numeric column `RR`
#' holding the rating, an integer from 1 (lowest) to `n_ratings` (highest).
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Interpretation** |
#' |---------------|---------------|-------------------|-------------|--------------------|
#' | *m*   | -     | \[-Inf, Inf\] | 1        | Mean log finishing time |
#' | *s*   | log   | \[0, Inf\]    | log(1)   | SD of log finishing time |
#' | *t0*  | log   | \[0, Inf\]    | log(0)   | Non-decision time |
#' | *rho* | probit on (-1, 1) | \[-1, 1\] | 0 | Correlation of the log finishing times (trial level, must not depend on lR or lM) |
#' | *c1*  | log   | \[0, Inf\]    | log(0.5) | First rating criterion, -log of the largest threshold |
#' | *ck*  | log   | \[0, Inf\]    | log(0.5) | Increment to the next criterion |
#'
#' The ratings are read from the losing accumulator's evidence at the decision,
#' as a proportion of its threshold, against thresholds
#' `d_{K-k} = exp(-(c1 + ... + ck))`, so rating `c` is given when that proportion
#' lies between `d_{K-c}` and `d_{K-c+1}`. The thresholds belong to the losing
#' accumulator, so the criteria are usually mapped `~ lR`. The natural-scale
#' thresholds `d1 .. d{K-1}` are reported with `add_recalculated = TRUE`.
#'
#' @param n_ratings Number of rating categories K (>= 1). With K = 1 the model
#'   is the binary correlated log-normal race.
#' @return A model list with all the necessary functions for EMC2 to sample.
#' @examples
#' # three ratings is the default; other numbers through a wrapper
#' MTLNR2 <- function() MTLNR(n_ratings = 2)
#' @export
MTLNR <- function(n_ratings = 3) {
  K <- rating_check_n(n_ratings)
  list(
    type = "RACE",
    c_name = "MTLNR",
    n_ratings = K,
    extra_responses = "RR",
    p_types = c("m" = 1, "s" = log(1), "t0" = log(0), "rho" = 0, rating_p_types(K)),
    transform = list(func = c(m = "identity", s = "exp", t0 = "exp", rho = "pnorm",
                              rating_transform(K)),
                     lower = c(rho = -1), upper = c(rho = 1)),
    bound = list(minmax = cbind(m = c(-Inf, Inf), s = c(0, Inf), t0 = c(0.05, Inf),
                                rho = c(-1, 1), rating_bound(K))),
    Ttransform = function(pars, dadm) rating_add_thresholds(pars, K),
    prepare_design = rating_prepare_design(n_choices = 2, trial_level = "rho"),
    check_data = function(data) mtlnr_check_data(data, K),
    rfun = function(data = NULL, pars) rMTLNR(data$lR, pars, K, ok = attr(pars, "ok")),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10))
      log_likelihood_mtlnr(pars, dadm, K, min_ll = min_ll)
  )
}
