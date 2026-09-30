# Stage M0 (2026-09-30): shared definitions for the MTLBA golden-fixture generators.
# Sourced by make_*.R in this directory. Not package code, not run by testthat.
#
# THE KERNEL (plan section 1). Two LBA accumulators, drifts bivariate normal with correlation
# rho, truncated to BOTH POSITIVE (norm = P(both > 0)), start points independent U(0, A).
# kernel(t) = density that the STOPPING accumulator s reaches threshold bs at decision time t
#             while the OTHER accumulator o has evidence in the band (dl, du).
# dl = 0 means the lowest band ("below du"); otherwise A_o <= dl < du (all thresholds above A).
# Balance of evidence: bs = the choice threshold b. Threshold count: bs = any threshold.

# ---- parameter sets ---------------------------------------------------------------------
# Accumulator order (1, 2). Thresholds: d_k = A + B * prop_k, k = 1..K-1, and b = A + B.
mtlba_sets <- list(
  sym    = list(A = c(.5, .5), B = c(1, 1),    v = c(2.5, 1.5), sv = c(1, 1),    t0 = .2,
                pow = c(1, 1)),
  asym   = list(A = c(.4, .8), B = c(1.2, .7), v = c(3, 1),     sv = c(.7, 1.2), t0 = .15,
                pow = c(1, 1.5)),
  lowv   = list(A = c(1, .3),  B = c(.5, 1.5), v = c(.8, -.3),  sv = c(1, .6),   t0 = .3,
                pow = c(.7, 1)),
  sharp  = list(A = c(.5, .5), B = c(1, 1.2),  v = c(2.5, 2),   sv = c(.25, .4), t0 = .2,
                pow = c(1, 1))
)
mtlba_rho <- c(0, -.7, -.3, .4, .8)
mtlba_dt  <- c(.1, .2, .35, .5, .8, 1.2, 2, 3.5)
mtlba_K   <- 2:4

# 2 x (K+1) matrix of thresholds per accumulator: columns DT0 = 0, DT1 .. DT{K-1}, b
# (the layout of `tmats` in BETE's log_likelihood_mt)
mtlba_thresholds <- function(set, K) {
  out <- t(sapply(1:2, function(a) {
    prop <- if (K > 1) ((1:(K - 1)) / K)^set$pow[a] else numeric(0)
    c(0, set$A[a] + set$B[a] * prop, set$A[a] + set$B[a])
  }))
  colnames(out) <- c("DT0", if (K > 1) paste0("DT", 1:(K - 1)), "b")
  out
}

# ---- piecewise integrate() --------------------------------------------------------------
pw_breaks <- function(lo, hi, centre = NULL, scale = NULL, extra = NULL) {
  b <- c(lo, hi, extra)
  if (!is.null(centre)) b <- c(b, centre + scale * c(-12, -6, -3, -1, 0, 1, 3, 6, 12))
  b <- sort(unique(pmin(pmax(b[is.finite(b)], lo), hi)))
  b
}
pw_int <- function(f, brk, rel.tol = 1e-13, abs.tol = 1e-22, ...) {
  if (length(brk) < 2) return(0)
  s <- 0
  for (i in seq_len(length(brk) - 1)) {
    if (brk[i + 1] - brk[i] <= 0) next
    s <- s + integrate(f, brk[i], brk[i + 1], rel.tol = rel.tol, abs.tol = abs.tol,
                       subdivisions = 2000L, stop.on.error = FALSE, ...)$value
  }
  s
}

# ---- normaliser P(both drifts > 0) ------------------------------------------------------
# one-dimensional: exact product at rho = 0, else integral of a conditional normal tail
ppos_1d <- function(v, sv, rho) {
  if (rho == 0) return(pnorm(v[1] / sv[1]) * pnorm(v[2] / sv[2]))
  sc <- sv[2] * sqrt(1 - rho^2)
  f <- function(x) dnorm(x, v[1], sv[1]) *
    pnorm((v[2] + rho * sv[2] / sv[1] * (x - v[1])) / sc)
  pw_int(f, pw_breaks(0, max(0, v[1]) + 40 * sv[1], v[1], sv[1]))
}
# bivariate normal density written out (no mvtnorm)
dbvn <- function(x, y, v, sv, rho) {
  zx <- (x - v[1]) / sv[1]; zy <- (y - v[2]) / sv[2]
  exp(-(zx^2 - 2 * rho * zx * zy + zy^2) / (2 * (1 - rho^2))) /
    (2 * pi * sv[1] * sv[2] * sqrt(1 - rho^2))
}
# two-dimensional brute force over the positive quadrant
ppos_2d <- function(v, sv, rho) {
  sc <- sv[2] * sqrt(1 - rho^2)
  outer_f <- function(xs) vapply(xs, function(x) {
    mc <- v[2] + rho * sv[2] / sv[1] * (x - v[1])
    pw_int(function(y) dbvn(x, y, v, sv, rho),
           pw_breaks(0, max(0, mc) + 40 * sc, mc, sc), abs.tol = 1e-25)
  }, 0)
  pw_int(outer_f, pw_breaks(0, max(0, v[1]) + 40 * sv[1], v[1], sv[1]))
}

# ---- kernel, one-dimensional form (plan section 1, "Correlated") -----------------------
# Integration variable = the stopping accumulator's drift x = (bs - start)/t, so
#   kernel = 1/(P++ * A_s) * Int_{(bs-A_s)/t}^{bs/t} x * dnorm(x; v_s, sv_s) * band(x) dx
# band(x) = P(dl < start_o + v_o t < du, v_o > 0 | v_s = x), closed form:
#   v_o | x ~ N(mc, sc);  H(c) = E_start[ pnorm(((c - start)/t - mc)/sc) ]
#   band = H(du) - H(dl)            (dl >= A_o: the v_o < 0 part cancels)
#   band = H(du) - pnorm(-mc/sc)    (dl = 0)
# Arguments: s = c(A, v, sv) of the stopping accumulator, o = same for the other one.
Hfun <- function(c, t, A, mc, sc) {
  G <- function(z) z * pnorm(z) + dnorm(z)
  (t * sc / A) * (G((c - t * mc) / (t * sc)) - G((c - A - t * mc) / (t * sc)))
}
kernel_1d <- function(t, As, bs, vs, svs, Ao, dl, du, vo, svo, rho, ppos = NULL) {
  stopifnot(du >= Ao, dl == 0 || dl >= Ao, du > dl, bs >= As, t > 0)
  if (is.null(ppos)) ppos <- ppos_1d(c(vs, vo), c(svs, svo), rho)
  sc <- svo * sqrt(1 - rho^2)
  f <- function(x) {
    mc <- vo + rho * svo / svs * (x - vs)
    band <- Hfun(du, t, Ao, mc, sc) -
      (if (dl == 0) pnorm(-mc / sc) else Hfun(dl, t, Ao, mc, sc))
    x * dnorm(x, vs, svs) * band
  }
  pw_int(f, pw_breaks((bs - As) / t, bs / t, vs, svs)) / (ppos * As)
}

# ---- kernel, brute-force two-dimensional form -------------------------------------------
# No conditional-normal algebra and no closed-form H: integrate the joint drift density
# directly over (x = stopping drift, y = other drift > 0), times the probability (over the
# other accumulator's uniform start point) that its evidence start + y t lies in (dl, du):
#   Ly(y) = length( [dl - y t, du - y t]  intersect  [0, A_o] ) / A_o
# The integrand is non-negative, so there is no cancellation in the tails.
kernel_2d <- function(t, As, bs, vs, svs, Ao, dl, du, vo, svo, rho, ppos = NULL) {
  stopifnot(du >= Ao, dl == 0 || dl >= Ao, du > dl, bs >= As, t > 0)
  v <- c(vs, vo); sv <- c(svs, svo)
  if (is.null(ppos)) ppos <- ppos_2d(v, sv, rho)
  sc <- svo * sqrt(1 - rho^2)
  ylo <- max(0, (dl - Ao) / t); yhi <- du / t
  kinks <- c((dl - Ao) / t, dl / t, (du - Ao) / t, du / t)
  outer_f <- function(xs) vapply(xs, function(x) {
    mc <- vo + rho * svo / svs * (x - vs)
    g <- function(y) dbvn(x, y, v, sv, rho) *
      pmax(pmin(du - y * t, Ao) - pmax(dl - y * t, 0), 0) / Ao
    x * pw_int(g, pw_breaks(ylo, yhi, mc, sc, kinks), abs.tol = 1e-25)
  }, 0)
  pw_int(outer_f, pw_breaks((bs - As) / t, bs / t, vs, svs)) / (ppos * As)
}
