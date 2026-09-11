#' The Ex-Gaussian Race Model
#'
#' Model file for a race model where each accumulator follows an ex-Gaussian
#' distribution with free mean, Gaussian scale, and exponential tail.
#'
#' @return A model list compatible with `design()`.
#' @export
REXG <- function() {
  dREXG <- function(rt, pars) {
    out <- numeric(length(rt))
    ok <- is.finite(rt) & pars[, "sigma"] > 0 & pars[, "tau"] > 0
    ok[is.na(ok)] <- FALSE

    if (any(ok)) {
      out[ok] <- dexGaussian(
        rt[ok] - pars[ok, "t0"],  # shift by t0
        pars[ok, c("mu", "sigma", "tau"), drop = FALSE]
      )
    }
    out
  }

  pREXG <- function(rt, pars) {
    out <- numeric(length(rt))
    ok <- is.finite(rt) & pars[, "sigma"] > 0 & pars[, "tau"] > 0
    ok[is.na(ok)] <- FALSE
    if (any(ok)) {
      out[ok] <- pexGaussian(
        rt[ok] - pars[ok, "t0"],                          # shift by t0
        pars[ok, c("mu", "sigma", "tau"), drop = FALSE]
      )
    }
    out[is.finite(rt) & rt == Inf] <- 1
    out
  }

  rREXG <- function(lR, pars, p_types = c("mu", "sigma", "tau", "t0"),
                    ok = rep(TRUE, dim(pars)[1])) {
    if (!all(p_types %in% dimnames(pars)[[2]]))
      stop("pars must have columns ", paste(p_types, collapse = " "))

    nr    <- length(levels(lR))
    n_tri <- nrow(pars) / nr
    dt    <- matrix(Inf, nrow = nr, ncol = n_tri)

    pars_ok  <- pars[ok, , drop = FALSE]
    mu_ok    <- pars_ok[, "mu"]
    sigma_ok <- pars_ok[, "sigma"]
    tau_ok   <- pars_ok[, "tau"]

    ok_draw <- is.finite(mu_ok) & sigma_ok > 0 & tau_ok > 0
    if (any(ok_draw)) {
      dt[ok][ok_draw] <-
        stats::rnorm(sum(ok_draw), mean = mu_ok[ok_draw], sd = sigma_ok[ok_draw]) +
        stats::rexp(sum(ok_draw), rate = 1 / tau_ok[ok_draw])
    }

    R    <- apply(dt, 2, which.min)
    pick <- cbind(R, seq_len(n_tri))
    t0   <- matrix(pars[, "t0"], nrow = nr)   # nr x n_tri, same layout as dt
    rt   <- dt[pick] + t0[pick]               # add winner's t0, same as rRDM

    R  <- factor(levels(lR)[R], levels = levels(lR))
    cbind.data.frame(R = R, rt = rt)
  }

  list(
    type = "RACE",
    c_name = "REXG",
    p_types = c(mu = log(.4), sigma = log(.05), tau = log(.1), t0=log(0)),
    p_types_canonical = c("mu", "sigma", "tau","t0"),
    transform = list(func = c(mu = "exp", sigma = "exp", tau = "exp", t0="exp")),
    bound = list(
      minmax = cbind(mu = c(0, Inf), sigma = c(1e-6, Inf), tau = c(1e-6, Inf), t0=c(0, Inf)),
      exception = c(t0 = 0)
    ),
    Ttransform = function(pars, dadm) pars,
    rfun = function(data = NULL, pars) rREXG(data$lR, pars, ok = attr(pars, "ok")),
    dfun = function(rt, pars) dREXG(rt, pars),
    pfun = function(rt, pars) pREXG(rt, pars),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      log_likelihood_race(pars = pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}
