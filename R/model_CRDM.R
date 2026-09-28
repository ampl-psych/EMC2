# The conflict racing diffusion model (CRDM): a racing diffusion in which an
# accumulator's mean path carries a conflict pulse (Lueken's eamax). There is
# no closed-form likelihood: the density is the numerical solution of a Volterra
# integral equation (src/crdm.cpp), which makes this model the exact (slow)
# control of the CRDM neural likelihoods and the twin they are registered with
# (register_nn_model(twin = CRDM, hybrid = "amp")).

# pars: natural-scale rows with columns v, amp, tau, s, b and t0.
crdm_dens <- function(rt, pars, dt, what) {
  out <- numeric(length(rt))
  t <- rt - pars[, "t0"]
  use <- !is.na(t) & t > 0
  if (!any(use)) return(out)
  P <- pars[use, c("v", "amp", "tau", "s", "b"), drop = FALSE]
  out[use] <- crdm_dens_rows(t[use], P, dt)[[what]]
  out
}

rCRDM <- function(lR, pars, ok = rep(TRUE, nrow(pars)), t_max = 3, sim_dt = 5e-4) {
  nr <- length(levels(lR)); n <- nrow(pars)
  if (is.null(ok)) ok <- rep(TRUE, n)
  if (any(pars[ok, "A"] != 0)) stop("CRDM: A must be exactly 0 (no start-point variability)")
  trial <- rep(seq_len(n / nr), each = nr)
  ft <- rep(Inf, n)
  pulsed <- ok & pars[, "amp"] > 0
  wd <- which(ok & !pulsed)
  # unpulsed accumulators: exact Wald draws
  if (length(wd)) ft[wd] <- rWald(length(wd), B = pars[wd, "b"] / pars[wd, "s"],
                                 v = pars[wd, "v"] / pars[wd, "s"], A = rep(0, length(wd)))
  if (any(pulsed)) {
    # a pulsed accumulator matters only if it finishes before the fastest
    # unpulsed accumulator of its trial, so it is simulated that far; with no
    # unpulsed accumulator, to t_max
    best <- as.numeric(tapply(ft, trial, min))[trial]
    ix <- which(pulsed)
    hz <- ifelse(is.finite(best[ix]), best[ix], t_max)
    hz <- pmax(sim_dt, ceiling(hz / sim_dt) * sim_dt)
    ft[ix] <- rcrdm_rows(pars[ix, "v"], pars[ix, "amp"], pars[ix, "tau"], pars[ix, "s"],
                         pars[ix, "b"], hz, sim_dt)
  }
  dtm <- matrix(ft, nrow = nr)
  any_ok <- colSums(matrix(ok, nrow = nr)) > 0
  if (any(!is.finite(apply(dtm, 2, min)) & any_ok))
    stop("CRDM: a trial without an unpulsed accumulator had no finisher by t_max = ", t_max,
         "; raise t_max")
  R <- apply(dtm, 2, which.min)
  pick <- cbind(R, seq_len(ncol(dtm)))
  data.frame(R = factor(levels(lR)[R], levels = levels(lR)),
             rt = matrix(pars[, "t0"], nrow = nr)[pick] + dtm[pick])
}

#' The Conflict Racing Diffusion Model
#'
#' Model file for the conflict racing diffusion model (CRDM): a racing
#' diffusion model ([RDM]) in which the mean path of an accumulator carries a
#' transient conflict pulse,
#' \eqn{M(t) = v t + amp\,(e\,t/\tau)\exp(-t/\tau)}{M(t) = v t + amp (e t / tau) exp(-t / tau)},
#' which rises to its peak height `amp` at time `tau` and then decays. With
#' `amp = 0` an accumulator is that of the [RDM].
#'
#' Model files are almost exclusively used in `design()`.
#'
#' @details
#'
#' Default values are used for all parameters that are not explicitly listed in the `formula`
#' argument of `design()`. They can also be accessed with `CRDM()$p_types`.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default**   | **Mapping**          | **Interpretation**                                                |
#' |-----------|-----------|---------------|-----------|------------------|---------------------------------------------------------------|
#' | *v*       | log       | \[0, Inf\]      | log(1)    |                  | Evidence-accumulation rate (drift rate)                        |
#' | *A*       | log       | 0             | log(0)    |                  | Start-point variability; only 0 is admitted                    |
#' | *B*       | log       | \[0, Inf\]      | log(1)    | *b* = *B* + *A*      | Response threshold                                             |
#' | *t0*      | log       | \[0, Inf\]      | log(0)    |                  | Non-decision time                                             |
#' | *s*       | log       | \[0, Inf\]      | log(1)    |                  | Within-trial standard deviation of drift rate                 |
#' | *amp*     | log       | \[0, Inf\]      | log(0)    |                  | Peak height of the pulse; exactly 0 (the default): no pulse    |
#' | *tau*     | log       | \[0, Inf\]      | log(0.1)  |                  | Time of the pulse's peak                                       |
#'
#' **Which accumulator is pulsed** is a property of the trial, as which
#' accumulator matches the stimulus is. Give `design()` a function that marks
#' the pulsed accumulator and let `amp` depend on it, with the unpulsed level
#' fixed at exactly zero:
#' `functions = list(lD = function(d) factor(d$D == d$lR, levels = c(FALSE, TRUE)))`,
#' `formula = list(..., amp ~ 0 + lD)`, `constants = c(amp_lDFALSE = -Inf)`.
#'
#' **Likelihood.** The finishing-time density of a pulsed accumulator has no
#' closed form. It is the solution of a Volterra integral equation of the second
#' kind (the Fortet recursion), solved on the grid `dt`, `2 dt`, ... and
#' interpolated linearly at the decision times `rt - t0`; the CDF is its
#' trapezoidal integral. One solve serves all the trials that share a parameter
#' row (`v`, `amp`, `tau`, `s`, `b`), adjacent or not, and reaches as far as
#' the largest decision time among them. An accumulator with `amp` exactly 0 is
#' the closed-form Wald. The cost of a solve grows with the square of the number
#' of grid points. For one solve to a decision time of 3 s:
#'
#' | **dt** | **time** | **largest error of the density, relative to its peak** |
#' |--------|----------|--------------------------------------------------------|
#' | 0.002  | 5 ms     | 5.9e-3 |
#' | 0.001  | 22 ms    | 2.1e-3 |
#' | 0.0005 | 86 ms    | 5.6e-4 |
#'
#' The likelihood runs in R with the solver in C++ (there is no compiled
#' likelihood pipeline for this model), so fits are slow; the model is meant
#' as the exact control for neural likelihoods of the CRDM (see
#' [register_nn_model()], `hybrid`, and [nn_cell()]), for which
#' `function() CRDM(dt = 0.002)` is a cheaper control.
#'
#' **Simulation.** Unpulsed accumulators are exact Wald draws. A pulsed
#' accumulator is simulated on the grid `sim_dt` (exact mean path, Brownian
#' increments, a Brownian-bridge test for crossings between grid points) only
#' as far as the fastest unpulsed accumulator of its trial, since it matters
#' only if it finishes first; nothing is redrawn. A trial in which every
#' accumulator is pulsed is simulated to `t_max` and must produce a finisher
#' by then (an error otherwise).
#'
#' The model has no start-point variability: `A` admits only 0 (its default).
#'
#' @param dt Step of the grid on which the density is solved (seconds).
#' @param t_max How far (seconds of decision time) a trial whose accumulators
#'   are all pulsed is simulated.
#' @param sim_dt Step of the simulator's grid (seconds).
#'
#' The model, its solver and its simulator follow Malte Lüken's `eamax`
#' (<https://github.com/maltelueken/eamax>).
#'
#' @return A list defining the cognitive model
#' @examples
#' # The accumulator of the response that the distractor D points to is pulsed
#' des <- design(factors = list(subjects = 1, S = c("left", "right"), D = c("left", "right")),
#'               Rlevels = c("left", "right"), matchfun = function(d) d$S == d$lR,
#'               functions = list(lD = function(d) factor(d$D == d$lR, levels = c(FALSE, TRUE))),
#'               formula = list(v ~ lM, B ~ 1, t0 ~ 1, tau ~ 1, amp ~ 0 + lD),
#'               constants = c(amp_lDFALSE = -Inf), model = CRDM)
#' p <- c(v = log(2), v_lMTRUE = .5, B = log(1), t0 = log(.3), tau = log(.1), amp_lDTRUE = log(.4))
#' mapped_pars(des, p)
#' @export
CRDM <- function(dt = 5e-4, t_max = 3, sim_dt = 5e-4) {
  if (!(is.numeric(dt) && length(dt) == 1L && dt > 0)) stop("dt must be a positive number")
  if (!(is.numeric(sim_dt) && length(sim_dt) == 1L && sim_dt > 0)) stop("sim_dt must be a positive number")
  if (!(is.numeric(t_max) && length(t_max) == 1L && t_max >= sim_dt)) stop("t_max must be at least sim_dt")
  list(
    type = "RACE",
    p_types = c("v" = log(1), "B" = log(1), "A" = log(0), "t0" = log(0), "s" = log(1),
                "amp" = log(0), "tau" = log(0.1)),
    transform = list(func = c(v = "exp", B = "exp", A = "exp", t0 = "exp", s = "exp",
                              amp = "exp", tau = "exp")),
    # RDM()'s bounds; A admits only 0; amp = 0 (no pulse) admitted
    bound = list(minmax = cbind(v = c(1e-3, Inf), B = c(0, Inf), A = c(0, 0), t0 = c(0.05, Inf),
                                s = c(0, Inf), amp = c(0, Inf), tau = c(0, Inf)),
                 exception = c(A = 0, v = 0, amp = 0)),
    # Trial dependent parameter transform
    Ttransform = function(pars, dadm) {
      pars <- cbind(pars, b = pars[, "B"] + pars[, "A"])
      pars
    },
    # Random function for racing accumulators
    rfun = function(data = NULL, pars)
      rCRDM(data$lR, pars, ok = attr(pars, "ok"), t_max = t_max, sim_dt = sim_dt),
    # Density function (PDF) for single accumulator
    dfun = function(rt, pars) crdm_dens(rt, pars, dt, "pdf"),
    # Probability function (CDF) for single accumulator
    pfun = function(rt, pars) crdm_dens(rt, pars, dt, "cdf"),
    # Race likelihood combining pfun and dfun
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10))
      log_likelihood_race(pars = pars, dadm = dadm, model = model, min_ll = min_ll)
  )
}
