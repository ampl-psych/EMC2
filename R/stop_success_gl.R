# ---------------------------------------------------------------------------
# Pure-R routes for the stop-success integral, with a method switch matching
# the C++ live dispatchers in src/model_SS_EXG.h / src/model_SS_RDEX.h:
#   "integrate" -> adaptive stats::integrate() via my.integrate() (untouched
#                  original route),
#   "gl"        -> fixed Gauss-Legendre with n_nodes,
#   "analytic"/"auto" (EXG only) -> exact n_go == 1 closed form where valid,
#                  GL (with a tight-stop-density node bump) otherwise.
#
# The quadrature routes reuse the exported C++ integrands stopfn_texg() /
# stopfn_rdex(), so the R and C++ routes integrate exactly the same function.
# ---------------------------------------------------------------------------

# Gauss-Legendre nodes/weights on [-1, 1] (Newton "gauleg"; matches gl_quad.h).
# Cached in a package-private environment keyed by node count.
.gl_cache_R <- new.env(parent = emptyenv())
gl_rule_R <- function(n) {
  key <- as.character(n)
  if (!is.null(.gl_cache_R[[key]])) return(.gl_cache_R[[key]])
  x <- numeric(n); w <- numeric(n)
  m <- (n + 1) %/% 2
  eps <- 1e-15
  for (i in seq_len(m)) {
    z <- cos(pi * (i - 0.25) / (n + 0.5))   # initial guess (1-based i)
    repeat {
      p1 <- 1.0; p2 <- 0.0
      for (j in 0:(n - 1)) {
        p3 <- p2; p2 <- p1
        p1 <- ((2 * j + 1) * z * p2 - j * p3) / (j + 1)
      }
      pp <- n * (z * p1 - p2) / (z * z - 1.0)
      z1 <- z; z <- z1 - p1 / pp
      if (abs(z - z1) <= eps) break
    }
    x[i] <- -z;            x[n + 1 - i] <- z
    w[i] <- 2.0 / ((1.0 - z * z) * pp * pp); w[n + 1 - i] <- w[i]
  }
  rule <- list(x = x, w = w)
  .gl_cache_R[[key]] <- rule
  rule
}

# Generic fixed-GL definite integral of a vectorized fn over [lower, upper].
gl_integrate_R <- function(f, lower, upper, n_nodes = 64L, ...) {
  if (!(upper > lower)) return(0)
  gl <- gl_rule_R(n_nodes)
  c1 <- (upper - lower) / 2
  c2 <- (upper + lower) / 2
  x  <- c1 * gl$x + c2
  c1 * sum(gl$w * f(x, ...))
}

# R mirror of C++ gl_auto_nodes() (gl_quad.h): hold a minimum node DENSITY per
# stop-sigma across the window so a sharp (small-tauS, near-Gaussian) stop peak
# is resolved. n_nodes is a floor; quantise up to GL_NODE_STEP, cap at
# GL_MAX_NODES. Keep these constants in sync with src/gl_quad.h.
GL_NODES_PER_SIG <- 6
GL_NODE_STEP     <- 32L
GL_MAX_NODES     <- 256L
gl_auto_nodes_R <- function(n_nodes, lower, upper, sigS) {
  n_nodes <- as.integer(n_nodes)
  if (!(sigS > 0) || !(upper > lower)) return(n_nodes)
  want <- GL_NODES_PER_SIG * (upper - lower) / sigS
  if (!(want > n_nodes)) return(n_nodes)
  n_eff <- as.integer(ceiling(want / GL_NODE_STEP)) * GL_NODE_STEP
  n_eff <- max(n_eff, n_nodes)
  as.integer(min(n_eff, GL_MAX_NODES))
}

#' Stop-success integral (truncated ex-Gaussian race), R route, method-selectable
#'
#' @param mu,sigma,tau,lb Numeric vectors c(stop, go1, go2, ...) of ex-Gaussian
#'   parameters and truncation lower bounds (stop first, then go accumulators).
#' @param SSD Stop-signal delay.
#' @param upper Upper integration limit; Inf -> full integral ("integrate") or
#'   heuristic muS + k_sigma*sigS + k_tau*tauS window ("gl").
#' @param method "integrate" (adaptive; original route, numerically untouched),
#'   "gl" (fixed Gauss-Legendre), "analytic" or "auto" (exact n_go == 1 closed
#'   form where valid, GL fallback otherwise; see stop_success_texg_analytic1_R).
#' @param n_nodes Number of GL nodes when method = "gl" (and the base node
#'   count for the "auto"/"analytic" GL fallback).
#' @param k_sigma,k_tau Finite-window heuristic for the GL routes.
#' @return The integral value (probability of successful stop, pre-trigger-failure).
#' @keywords internal
stop_success_texg_R <- function(mu, sigma, tau, lb, SSD, upper = Inf,
                                method = c("integrate", "gl", "analytic", "auto"),
                                n_nodes = 64L,
                                # defaults mirror SS_WINDOW_K_SIGMA/K_TAU in
                                # src/gl_quad.h — keep in sync
                                k_sigma = 8, k_tau = 16) {
  method <- match.arg(method)
  lower  <- lb[1]
  if (method == "integrate") {
    if (!(upper > lower)) return(0)
    return(my.integrate(f = stopfn_texg, lower = lower, upper = upper,
                        mu = mu, sigma = sigma, tau = tau, lb = lb, SSD = SSD))
  }
  bump <- FALSE
  if (method %in% c("analytic", "auto")) {
    if (length(mu) == 2L && !is.finite(upper)) {
      p <- stop_success_texg_analytic1_R(mu, sigma, tau, lb, SSD)
      if (!is.na(p)) return(p)
    }
    bump <- TRUE   # GL fallback with the tight-stop-density node bump
  }
  # GL needs a finite window; mirror the C++ heuristic on both sides
  if (!is.finite(lower)) lower <- mu[1] - k_sigma * sigma[1] - k_tau * tau[1]
  ub <- if (is.finite(upper)) upper else mu[1] + k_sigma * sigma[1] + k_tau * tau[1]
  if (!(ub > lower)) return(0)
  if (bump) n_nodes <- gl_auto_nodes_R(n_nodes, lower, ub, sigma[1])
  gl_integrate_R(function(x) stopfn_texg(x, mu, sigma, tau, lb, SSD),
                 lower, ub, n_nodes)
}

#' Stop-success integral (RDEX race), R route, method-selectable
#'
#' No analytic form exists for RDEX (Wald survivor arguments are nonlinear in
#' t), so method "auto" is GL with the tight-stop-density node bump.
#' @keywords internal
stop_success_rdex_R <- function(n_acc, mu, sigma, tau, lb, v, B, A, t0, s, SSD,
                                upper = Inf,
                                method = c("integrate", "gl", "auto"),
                                n_nodes = 64L,
                                # defaults mirror SS_WINDOW_K_SIGMA/K_TAU in
                                # src/gl_quad.h — keep in sync
                                k_sigma = 8, k_tau = 16) {
  method <- match.arg(method)
  lower  <- lb[1]
  intfn <- function(x)
    stopfn_rdex(x, n_acc, mu, sigma, tau, lb, v, B, A, t0, s, SSD)
  if (method == "integrate") {
    if (!(upper > lower)) return(0)
    return(my.integrate(f = intfn, lower = lower, upper = upper))
  }
  if (!is.finite(lower)) lower <- mu[1] - k_sigma * sigma[1] - k_tau * tau[1]
  ub <- if (is.finite(upper)) upper else mu[1] + k_sigma * sigma[1] + k_tau * tau[1]
  if (!(ub > lower)) return(0)
  if (method == "auto") n_nodes <- gl_auto_nodes_R(n_nodes, lower, ub, sigma[1])
  gl_integrate_R(intfn, lower, ub, n_nodes)
}

# ---------------------------------------------------------------------------
# Exact n_go == 1 closed form (R mirror of ss_texg_stop_success_analytic1 in
# src/ss_exg_analytic.h; derivation in WorkingTests/stop_success_methods.R (Section 1)).
# Full-line (lbS = lbG = -Inf): 4 pnorm. Truncated: 4 bivariate Phi2 + 2
# boundary terms. Returns NA when a guard trips (caller falls back to GL):
#   * sigS > ANALYTIC_MAX_SIGS_OVER_TAUS * tauS (10; relaxed from 4 once the C++
#     dexg() tail gained the higher-order Mills corrections, and to bound the
#     logC0 cancellation as tauS -> 0 — see src/ss_exg_analytic.h),
#   * kinked domain (lbS + SSD < lbG, or lbS = -Inf with finite lbG),
#   * any term exponent > 600 (exp overflow) or non-finite intermediate.
# ---------------------------------------------------------------------------

# untruncated ex-Gaussian survivor (log-space Phi terms, stable tails)
.sexg_R <- function(t, mu, sigma, tau) {
  pnorm(-(t - mu) / sigma) +
    exp((mu - t) / tau + sigma^2 / (2 * tau^2) +
          pnorm((t - mu) / sigma - sigma / tau, log.p = TRUE))
}

stop_success_texg_analytic1_R <- function(mu, sigma, tau, lb, SSD) {
  muS <- mu[1]; sigS <- sigma[1]; tauS <- tau[1]; lbS <- lb[1]
  muG <- mu[2]; sigG <- sigma[2]; tauG <- tau[2]; lbG <- lb[2]
  ANALYTIC_MAX_SIGS_OVER_TAUS <- 10   # keep in sync with src/ss_exg_analytic.h
  if (sigS > ANALYTIC_MAX_SIGS_OVER_TAUS * tauS) return(NA_real_)
  full_line <- !is.finite(lbS) && !is.finite(lbG)
  if (!full_line) {
    if (!is.finite(lbS)) return(NA_real_)
    if (lbS + SSD < lbG) return(NA_real_)
  }
  m <- muG - SSD
  a0 <- -muS / sigS - sigS / tauS
  b0 <- 1 / sigS
  logC0 <- -log(tauS) + muS / tauS + sigS^2 / (2 * tauS^2)
  total <- 0
  for (A in 0:1) {
    if (A == 0) {           # go survivor's Gaussian branch
      gam <- -1 / tauS; logC <- logC0
      a1 <- m / sigG;   b1 <- -1 / sigG
    } else {                # go survivor's exponential branch
      gam <- -1 / tauS - 1 / tauG
      logC <- logC0 + m / tauG + sigG^2 / (2 * tauG^2)
      a1 <- -m / sigG - sigG / tauG; b1 <- 1 / sigG
    }
    a <- c(a0, a1); b <- c(b0, b1)
    if (!full_line) {       # boundary term at L = lbS
      e <- logC + gam * lbS
      if (!is.finite(e) || e > 600) return(NA_real_)
      total <- total + exp(e) / abs(gam) *
        pnorm(a0 + b0 * lbS) * pnorm(a1 + b1 * lbS)
    }
    for (k in 1:2) {
      xbar <- -a[k] / b[k]
      mustar <- xbar + gam / b[k]^2
      e <- logC + gam * xbar + gam^2 / (2 * b[k]^2)
      if (!is.finite(e) || e > 600) return(NA_real_)
      j <- 3 - k
      beta <- b[j] / abs(b[k])
      alpha <- a[j] + b[j] * mustar
      s <- sqrt(1 + beta^2)
      Q <- if (full_line) pnorm(alpha / s)
           else pbvn_R(-abs(b[k]) * (lbS - mustar), alpha / s, beta / s)
      total <- total + sign(b[k]) * exp(e) * Q / abs(gam)
    }
  }
  if (!is.finite(total) || total <= 0) return(NA_real_)
  p <- total
  if (!full_line) {
    N <- .sexg_R(lbS, muS, sigS, tauS) *
      (if (is.finite(lbG)) .sexg_R(lbG, muG, sigG, tauG) else 1)
    if (!is.finite(N) || N <= 0) return(NA_real_)
    p <- total / N
  }
  if (!is.finite(p) || p <= 0 || p > 1 + 1e-10) return(NA_real_)
  min(p, 1)
}

# ---------------------------------------------------------------------------
# Bivariate normal CDF P(X <= h, Y <= k), correlation rho: scalar R port of
# Alan Genz's BVND (TVPACK), the same algorithm as bvn_cdf() in
# src/ss_exg_analytic.h, so the R and C++ analytic routes agree to ~1e-15.
# Plain R on purpose: no pbivnorm/mvtnorm dependency.
# ---------------------------------------------------------------------------
.bvn_gl_half <- list(
  list(x = c(0.9324695142031522, 0.6612093864662647, 0.2386191860831970),
       w = c(0.1713244923791705, 0.3607615730481384, 0.4679139345726904)),
  list(x = c(0.9815606342467191, 0.9041172563704750, 0.7699026741943050,
             0.5873179542866171, 0.3678314989981802, 0.1252334085114692),
       w = c(0.04717533638651177, 0.1069393259953183, 0.1600783285433464,
             0.2031674267230659, 0.2334925365383547, 0.2491470458134029)),
  list(x = c(0.9931285991850949, 0.9639719272779138, 0.9122344282513259,
             0.8391169718222188, 0.7463319064601508, 0.6360536807265150,
             0.5108670019508271, 0.3737060887154196, 0.2277858511416451,
             0.07652652113349733),
       w = c(0.01761400713915212, 0.04060142980038694, 0.06267204833410906,
             0.08327674157670475, 0.1019301198172404, 0.1181945319615184,
             0.1316886384491766, 0.1420961093183821, 0.1491729864726037,
             0.1527533871307259))
)

pbvn_R <- function(h, k, rho) {
  # Genz BVND computes the upper tail P(X > dh, Y > dk); flip signs.
  dh <- -h; dk <- -k; r <- rho
  if (is.na(dh) || is.na(dk) || is.na(r)) return(NA_real_)
  if (dh == Inf || dk == Inf) return(0)
  if (dh == -Inf) return(if (dk == -Inf) 1 else pnorm(-dk))
  if (dk == -Inf) return(pnorm(-dh))
  gl <- if (abs(r) < 0.3) .bvn_gl_half[[1]]
        else if (abs(r) < 0.75) .bvn_gl_half[[2]]
        else .bvn_gl_half[[3]]
  # rule is symmetric: use +x and -x halves
  x <- c(-gl$x, gl$x); w <- c(gl$w, gl$w)
  hh <- dh; kk <- dk; hk <- hh * kk; bvn <- 0
  if (abs(r) < 0.925) {
    if (abs(r) > 0) {
      hs <- (hh * hh + kk * kk) / 2
      asr <- asin(r)
      sn <- sin(asr * (x + 1) / 2)
      bvn <- sum(w * exp((sn * hk - hs) / (1 - sn^2))) * asr / (4 * pi)
    }
    bvn <- bvn + pnorm(-hh) * pnorm(-kk)
  } else {
    if (r < 0) { kk <- -kk; hk <- -hk }
    if (abs(r) < 1) {
      as_ <- (1 - r) * (1 + r); a <- sqrt(as_)
      bs <- (hh - kk)^2
      cc <- (4 - hk) / 8; dd <- (12 - hk) / 16
      asr <- -(bs / as_ + hk) / 2
      if (asr > -100)
        bvn <- a * exp(asr) *
          (1 - cc * (bs - as_) * (1 - dd * bs / 5) / 3 + cc * dd * as_^2 / 5)
      if (-hk < 100) {
        b <- sqrt(bs)
        bvn <- bvn - exp(-hk / 2) * sqrt(2 * pi) * pnorm(-b / a) * b *
          (1 - cc * bs * (1 - dd * bs / 5) / 3)
      }
      a <- a / 2
      xs <- (a * (x + 1))^2
      rs <- sqrt(1 - xs)
      asr <- -(bs / xs + hk) / 2
      ok <- asr > -100
      if (any(ok))
        bvn <- bvn + sum(a * w[ok] * exp(asr[ok]) *
          (exp(-hk * (1 - rs[ok]) / (2 * (1 + rs[ok]))) / rs[ok] -
             (1 + cc * xs[ok] * (1 + dd * xs[ok]))))
      bvn <- -bvn / (2 * pi)
    }
    if (r > 0) {
      bvn <- bvn + pnorm(-max(hh, kk))
    } else {
      bvn <- -bvn
      if (kk > hh) {
        # Phi(kk) - Phi(hh); use upper tails when both positive for accuracy
        bvn <- bvn + if (hh >= 0) pnorm(-hh) - pnorm(-kk)
                     else pnorm(kk) - pnorm(hh)
      }
    }
  }
  min(1, max(0, bvn))
}

# ---------------------------------------------------------------------------
# Push a model's stop_method / stop_n_nodes into the process-global C++
# stop-success config (emc2_set_stop_method(); see gl_quad.h for the threading
# rationale). Called once per likelihood call by calc_ll_manager() and the
# direct calc_ll users, BEFORE any C++ likelihood evaluation, so forked
# (mclapply) workers inherit it and PSOCK workers set it themselves.
# Models without the fields (e.g. saved designs predating stop_method) get the
# package default "auto" / 64. No-op for non-stop-signal models.
# ---------------------------------------------------------------------------
set_stop_method_from_model <- function(model) {
  c_name <- model$c_name
  if (is.null(c_name) || !c_name %in% c("SSEXG", "SSRDEX"))
    return(invisible(NULL))
  emc2_set_stop_method(model$stop_method %||% "auto",
                       as.integer(model$stop_n_nodes %||% 64L))
  invisible(NULL)
}
