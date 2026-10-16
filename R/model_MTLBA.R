# Multiple-threshold LBA (MTLBA): two LBA accumulators (start point U(0, A),
# drift rate normal with mean v and sd sv), each with K - 1 confidence
# thresholds below its choice threshold b = A + B,
#   T_j = A + B * d_j,  j = 1 .. K - 1,   d_1 < ... < d_{K-1} < 1,
# (d_j from the rating layer's criteria, R/model_rating.R), T_K = b. A trial
# gives a choice R, its rt and a rating RR (1 = lowest).
#
# Drift convention (drifts = "any" | "positive"):
#   "any" (the default; the model of the source paper and its DMC code): the
#     two drifts are bivariate normal conditioned only on NOT BOTH being
#     negative (a trial on which neither accumulator ever finishes), so
#     norm = P(not both < 0) = 1 - Phi(-v1/sv1) Phi(-v2/sv2) (Genz bvnd when
#     correlated). A losing accumulator may have a negative drift and fall
#     below its start point; its band probabilities come from the untruncated
#     normal.
#   "positive": both drifts truncated to be positive, norm = P(both > 0); at
#     rho = 0 each accumulator is EMC2's LBA(). The lowest band then carries
#     the negative-drift mass away (the Phi(v/sv) sentinels below).
# Lowest threshold (lowest_threshold = "above_A" | "free"): with "above_A" all
# d_j are in (0, 1), so A < T_1; with "free" (balance of evidence only) the
# lowest threshold is T_1 = T_2 exp(-c_{K-1}), anywhere in (0, T_2), so d_1
# may be negative (T_1 below A). The band functions take any threshold.
#
# Parameterization (parameterization = "additive" | "proportional"): the above
# is "additive" (b = A + B, A absolute). "proportional" is the source paper's
# nesting: B holds b itself, T_{K-1} = b c1, A = T_{K-1} A (K = 1: A = b A) and
# each lower threshold a proportion of the one above, T_j = A + (T_{j+1} - A)
# c_{K-j}, or T_1 = T_2 c_{K-1} when the lowest threshold is free (c's and A
# proportions in (0, 1)). Ttransform translates such rows exactly into the
# additive model's natural columns (mtlba_proportional; the compiled twin in
# src/model_mtlba.h), so everything below sees only the additive model.
#
# Every response is one KERNEL: the density that the STOPPING accumulator s
# reaches its threshold T_{s,j} at decision time dt = rt - t0 while the OTHER
# accumulator's evidence e_o = start + drift * dt lies in its band
# (T_{o,lo}, T_{o,lo+1}), lo = 0 meaning "below T_{o,1}":
#   kernel = f(dt; A_s, T_{s,j}, v_s, sv_s) * P(band) / norm,
# with f the LBA density of reaching T_{s,j} from an untruncated normal drift
# (only positive drifts get there), and P(band) the joint probability with
# v_o > 0 under "positive" (so that at rho = 0 kernel = f_LBA * P(band | v_o > 0),
# both terms as EMC2's LBA()), the plain probability under "any".
# Two stopping rules share the kernel and differ in the response lookup
# (mtlba_lookup):
#   BE, balance of evidence: the first accumulator to reach b stops and is the
#     choice; the rating is read from the other's band (as MTLNR).
#   TC, threshold count: the race stops at the K-th threshold crossing summed
#     over both accumulators; the choice is the accumulator ahead in crossings
#     (on a tie the one that got there first, i.e. the non-stopping one).
# Don't know (dont_know = TRUE): rating 1 of both choices is one response, whose
# likelihood is the SUM of the two rating-1 kernels.
# t0 is read from the stopping accumulator's row (= the winner's under BE).
#
# The band probability is formed on the log scale. With H(y) = P(e_o <= y) and
# S(y) = 1 - H(y) for untruncated drift (closed forms through
# G(z) = z Phi(z) + phi(z), valid for any y), under "any"
#   band = H(hi) - H(lo) = S(lo) - S(hi),  H(lowest) = 0, S(lowest) = 1;
# under "positive", thresholds above A make the negative-drift mass cancel
# except in the lowest band, where H(0) := Phi(-v/sv), S(0) := Phi(v/sv)
# (plan M0 results), so
#   band * Phi(v/sv) = S(lo) - S(hi) = H(hi) - H(lo),
# (a threshold below A takes the positive-drift forms Hp, Sp of mtlba_logHSp),
# and whichever route has the smaller leading term is used.
#
# Correlated drift rates (rho != 0): the two drifts are bivariate normal with
# correlation rho. Given the stopping accumulator's drift x (its start point
# is then bs - x * dt) the other drift is normal with
# mean v_o + rho * (sv_o / sv_s) * (x - v_s) and sd sv_o * sqrt(1 - rho^2), so
#   kernel = 1 / (A_s * norm) *
#            integral over x in ((bs - A_s) / dt, bs / dt) of
#              x * dnorm(x; v_s, sv_s) * P(band [, v_o > 0] | x) dx,
# a one-dimensional integral of an analytic integrand whose last factor is the
# band function above at the conditional mean and sd, left unnormalised. It is
# done by Gauss-Legendre quadrature on a window clipped to where the integrand
# lives (mtlba_quad_rule, mtlba_logkernel_rho). rho is a trial-level parameter
# (the same on both accumulator rows); trials with rho = 0 take the closed form.
#
# Confidence guessing (gp > 0; see the "confidence guessing" block below): a
# triggering accumulator ends the race at its guess threshold gb = A + B * gB
# with a random rating, so a response is a weighted sum of up to four kernels.
#
# The compiled likelihood (src/model_mtlba.h, c_name "MTLBA") ports this
# algorithm, including the correlated kernel with the default rule
# mtlba_quad_rule(), the correlated stopping terms and guessing; the R functions
# here are its reference, used when c_name is removed.

# log G(z), G(z) = z * pnorm(z) + dnorm(z) = integral of pnorm up to z, accurate
# in the lower tail: G(-x) = dnorm(x) * (1 - x R(x)), R the Mills ratio, with
# 1 - x R(x) = q / (x + q) from Laplace's continued fraction for large x.
mtlba_logG <- function(z) {
  out <- numeric(length(z))
  nz <- is.na(z)
  out[nz] <- NA
  pos <- !nz & z >= 0
  out[pos] <- log(z[pos] * stats::pnorm(z[pos]) + stats::dnorm(z[pos]))
  x <- -z[!nz & !pos]
  if (length(x)) {
    lg <- numeric(length(x))
    mid <- x < 5
    if (any(mid)) {
      xm <- x[mid]
      xR <- xm * exp(stats::pnorm(-xm, log.p = TRUE) - stats::dnorm(xm, log = TRUE))
      lg[mid] <- stats::dnorm(xm, log = TRUE) + log1p(-xR)
    }
    if (any(!mid)) {
      xt <- x[!mid]
      tt <- xt
      # 40 terms: bit-identical to 120 for x in [5, 1e8] (25 already suffice at x = 5)
      for (k in 40:2) tt <- xt + k / tt
      q <- 1 / tt
      lg[!mid] <- stats::dnorm(xt, log = TRUE) + log(q) - log(xt + q)
    }
    out[!nz & !pos] <- lg
  }
  out
}

# log(exp(la) - exp(lb)), -Inf where la <= lb (rounding) or la = -Inf.
mtlba_ldiff <- function(la, lb) {
  out <- rep(-Inf, length(la))
  ok <- !is.na(la) & !is.na(lb) & la > lb
  if (any(ok)) out[ok] <- log_diff_exp(la[ok], lb[ok])
  out[is.na(la) | is.na(lb)] <- NA
  out
}

# log H(y) and log S(y) for evidence e = U(0, A) + v * t, v ~ N(v, sv)
# untruncated, y >= A: H = mean of pnorm over the window [z0, z1],
# z1 = (y - t v) / (t sv), z0 = z1 - al, al = A / (t sv); S likewise with
# pnorm(-u). In general the mean is the difference of G at the two ends over
# al, which loses digits as al -> 0 (error ~ 2e-14 / al in the log). So a
# narrow window, al * max(1, |z1|) < 0.1, takes the mean by 4-node
# Gauss-Legendre instead: positive terms, error < 1e-13 there
# (rating-work/mtlba/stageM5/smallA_probe.R).
mtlba_logHS <- function(y, t, A, v, sv) {
  n <- length(y)
  z1 <- (y - t * v) / (t * sv)
  lH <- lS <- numeric(n)
  small <- A <= 1e-10  # A = 0: the LBA limit (as EMC2's LBA_A_ASYMPTOTIC)
  if (any(small)) {
    lH[small] <- stats::pnorm(z1[small], log.p = TRUE)
    lS[small] <- stats::pnorm(z1[small], lower.tail = FALSE, log.p = TRUE)
  }
  w <- !small
  if (any(w)) {
    al <- A[w] / (t[w] * sv[w])
    za <- z1[w]
    zb <- za - al
    lH[w] <- mtlba_ldiff(mtlba_logG(za), mtlba_logG(zb)) - log(al)
    lS[w] <- mtlba_ldiff(mtlba_logG(-zb), mtlba_logG(-za)) - log(al)
    nar <- which(al * pmax(1, abs(za)) < 0.1)
    if (length(nar)) {
      g <- mtlba_gl(4)
      U <- zb[nar] + outer(al[nar], (g$x + 1) / 2)
      lw <- rep(log(g$w / 2), each = length(nar))
      lmean <- function(L) {
        M <- L[cbind(seq_len(nrow(L)), max.col(L, ties.method = "first"))]
        out <- M + log(rowSums(exp(L - M)))
        out[M == -Inf] <- -Inf
        out
      }
      lH[w][nar] <- lmean(stats::pnorm(U, log.p = TRUE) + lw)
      lS[w][nar] <- lmean(stats::pnorm(U, lower.tail = FALSE, log.p = TRUE) + lw)
    }
  }
  list(H = lH, S = lS)
}

# log P(e <= y, v > 0) and log P(e > y, v > 0) for e = U(0, A) + v t,
# v ~ N(v, sv), for ANY y >= 0, so that a confidence threshold may lie below A
# (lowest_threshold = "free" under "positive" drifts). For y >= A the values
# are H(y) - Phi(-zv) and S(y) (every negative-drift trial ends below
# A <= y). For 0 <= y < A, with zv = v / sv, z1 = (y - t v) / (t sv),
# al = A / (t sv) (every term positive):
#   Hp(y) = [G(z1) - G(-zv)] / al - (y / A) Phi(-zv)
#   Sp(y) = ((A - y) / A) Phi(zv) + [G(zv) - G(-z1)] / al
# (Hp + Sp = Phi(zv), using G(z) - G(-z) = z). Verified against brute-force
# integration in rating-work/mtlba/stageM8/explore/testB.R. Hp is a double
# difference that cancels as y -> 0 (the start points that can end below y
# are s < y, and the drift interval (0, (y - s) / t) shrinks with it), so a
# narrow one, dz = y / (t sv) with dz max(1, |z1|) < 0.1, takes the mean over
# s = y u, u in (0, 1), of the interval probability P(-zv < Z < z1 - dz u) by
# 4-node Gauss-Legendre (positive terms; the integrand is nearly linear there).
mtlba_logHSp <- function(y, t, A, v, sv) {
  hs <- mtlba_logHS(y, t, A, v, sv)
  zv <- v / sv
  lHp <- mtlba_ldiff(hs$H, stats::pnorm(-zv, log.p = TRUE))
  lSp <- hs$S
  below <- !is.na(y) & !is.na(A) & y < A & A > 1e-10
  if (any(below)) {
    w <- which(below)
    yb <- y[w]; tb <- t[w]; Ab <- A[w]; vb <- v[w]; svb <- sv[w]; zvb <- zv[w]
    al <- Ab / (tb * svb)
    z1 <- (yb - tb * vb) / (tb * svb)
    lHp[w] <- mtlba_ldiff(mtlba_ldiff(mtlba_logG(z1), mtlba_logG(-zvb)) - log(al),
                          log(yb / Ab) + stats::pnorm(-zvb, log.p = TRUE))
    lSp[w] <- mtlba_lse(log((Ab - yb) / Ab) + stats::pnorm(zvb, log.p = TRUE),
                        mtlba_ldiff(mtlba_logG(zvb), mtlba_logG(-z1)) - log(al))
    dz <- yb / (tb * svb)
    nar <- which(dz * pmax(1, abs(z1)) < 0.1)
    if (length(nar)) {
      g <- mtlba_gl(4)
      U <- outer(dz[nar], (g$x + 1) / 2)
      # (the upper limit clamped: rounding can put z1 an ulp below -zv when y = 0)
      L <- suppressWarnings(matrix(log_pnorm_interval(rep(-zvb[nar], 4), pmax(z1[nar] - U, -zvb[nar])), length(nar), 4)) +
        rep(log(g$w / 2), each = length(nar))
      L[is.nan(L)] <- -Inf   # a one-ulp interval (log pnorm is not monotone to the ulp): zero
      M <- L[cbind(seq_along(nar), max.col(L, ties.method = "first"))]
      lm <- M + log(rowSums(exp(L - M)))
      lm[M == -Inf] <- -Inf
      lHp[w][nar] <- log(yb[nar] / Ab[nar]) + lm
    }
  }
  list(H = lHp, S = lSp)
}

# log probability of the band (T_lo, T_hi) of the other accumulator at
# decision time t; lowest = lo is the sentinel "no threshold below" (lo_y is
# then ignored).
#   pos = TRUE  (drifts "positive"): log P(T_lo < e_o < T_hi | v_o > 0); the
#     denominator pnorm(v/sv) is floored at 1e-10 as in EMC2's LBA, and
#     normalise = FALSE leaves it out: log P(T_lo < e_o < T_hi, v_o > 0).
#     Thresholds above A use H and S with the Phi(v/sv) sentinels; a threshold
#     below A takes the positive-drift forms Hp / Sp.
#   pos = FALSE (drifts "any"): log P(T_lo < e_o < T_hi) for the untruncated
#     normal drift, H(lowest) = 0, S(lowest) = 1 (normalise is ignored).
mtlba_logband <- function(t, A, v, sv, lo_y, hi_y, lowest, normalise = TRUE, pos = TRUE) {
  hi <- mtlba_logHS(hi_y, t, A, v, sv)
  if (!pos) {
    lo <- list(H = rep(-Inf, length(t)), S = rep(0, length(t)))
    if (any(!lowest)) {
      l2 <- mtlba_logHS(lo_y[!lowest], t[!lowest], A[!lowest], v[!lowest], sv[!lowest])
      lo$H[!lowest] <- l2$H
      lo$S[!lowest] <- l2$S
    }
    useS <- lo$S < hi$H
    useS[is.na(useS)] <- FALSE
    lb <- mtlba_ldiff(hi$H, lo$H)
    if (any(useS)) lb[useS] <- mtlba_ldiff(lo$S[useS], hi$S[useS])
    return(lb)
  }
  zv <- v / sv
  # rows with a threshold below A take the general Hp / Sp route
  below <- (!is.na(hi_y) & hi_y < A) | (!lowest & !is.na(lo_y) & lo_y < A)
  below[is.na(below)] <- FALSE
  if (any(below)) {
    w <- which(below)
    out <- numeric(length(t))
    hib <- mtlba_logHSp(hi_y[w], t[w], A[w], v[w], sv[w])
    lob <- list(H = rep(-Inf, length(w)), S = stats::pnorm(zv[w], log.p = TRUE))
    nl <- !lowest[w]
    if (any(nl)) {
      l2 <- mtlba_logHSp(lo_y[w][nl], t[w][nl], A[w][nl], v[w][nl], sv[w][nl])
      lob$H[nl] <- l2$H; lob$S[nl] <- l2$S
    }
    useS <- lob$S < hib$H; useS[is.na(useS)] <- FALSE
    lb <- mtlba_ldiff(hib$H, lob$H)
    if (any(useS)) lb[useS] <- mtlba_ldiff(lob$S[useS], hib$S[useS])
    if (normalise) lb <- lb - pmax(stats::pnorm(zv[w], log.p = TRUE), log(1e-10))
    out[w] <- lb
    if (all(below)) return(out)
    k <- which(!below)
    out[k] <- mtlba_logband(t[k], A[k], v[k], sv[k], lo_y[k], hi_y[k], lowest[k], normalise)
    return(out)
  }
  lo <- list(H = stats::pnorm(-zv, log.p = TRUE), S = stats::pnorm(zv, log.p = TRUE))
  lo$H <- rep_len(lo$H, length(t)); lo$S <- rep_len(lo$S, length(t))
  if (any(!lowest)) {
    l2 <- mtlba_logHS(lo_y[!lowest], t[!lowest], A[!lowest], v[!lowest], sv[!lowest])
    lo$H[!lowest] <- l2$H
    lo$S[!lowest] <- l2$S
  }
  useS <- lo$S < hi$H
  useS[is.na(useS)] <- FALSE
  lb <- mtlba_ldiff(hi$H, lo$H)
  if (any(useS)) lb[useS] <- mtlba_ldiff(lo$S[useS], hi$S[useS])
  if (!normalise) return(lb)
  lb - pmax(stats::pnorm(zv, log.p = TRUE), log(1e-10))
}

# log LBA density of reaching threshold bs at t from an untruncated normal
# drift (EMC2's dlba formula, all terms taken on the scale of the largest);
# pos = TRUE divides by P(v > 0) (the positive-drift density, as LBA()).
mtlba_logf <- function(t, A, bs, v, sv, pos = TRUE) {
  out <- rep(-Inf, length(t))
  lden <- if (pos) pmax(stats::pnorm(v / sv, log.p = TRUE), log(1e-10)) else rep(0, length(t))
  small <- A <= 1e-10
  if (any(small)) {
    ts <- t[small]
    out[small] <- stats::dnorm(bs[small] / ts, v[small], sv[small], log = TRUE) +
      log(bs[small]) - 2 * log(ts) - lden[small]
  }
  w <- !small
  if (any(w)) {
    cz <- (bs[w] - t[w] * v[w]) / (t[w] * sv[w])
    czm <- cz - A[w] / (t[w] * sv[w])
    lp1 <- stats::dnorm(cz, log = TRUE)
    lp0 <- stats::dnorm(czm, log = TRUE)
    # pnorm(cz) - pnorm(czm); NaN (e.g. sv = 0, out of bounds) stays NaN -> -Inf
    lpi <- rep(NaN, length(cz))
    okz <- !is.na(cz) & !is.na(czm)
    lpi[okz] <- log_pnorm_interval(czm[okz], cz[okz])
    # common scale = the largest of the three terms, so none overflows (the
    # pnorm difference can be ~1 while both dnorms are ~exp(-1000))
    M <- pmax(lp1, lp0, lpi)
    dP <- exp(lpi - M)
    inner <- v[w] * dP + sv[w] * (exp(lp0 - M) - exp(lp1 - M))
    lo <- rep(-Inf, length(cz))
    pos <- !is.na(inner) & inner > 0
    lo[pos] <- log(inner[pos]) + M[pos]
    out[w] <- lo - log(A[w]) - lden[w]
  }
  out
}

# ---- correlated drift rates (rho != 0) ----

# Gauss-Legendre nodes and weights on (-1, 1): Golub-Welsch eigenvalues polished
# by Newton steps on the Legendre recurrence (nodes to ~1e-16). Cached per n.
mtlba_gl_cache <- new.env(parent = emptyenv())
mtlba_gl <- function(n) {
  key <- as.character(n)
  if (!is.null(mtlba_gl_cache[[key]])) return(mtlba_gl_cache[[key]])
  if (n == 1) return(mtlba_gl_cache[[key]] <- list(x = 0, w = 2))
  i <- seq_len(n - 1)
  J <- matrix(0, n, n)
  J[cbind(i, i + 1)] <- J[cbind(i + 1, i)] <- i / sqrt(4 * i^2 - 1)
  x <- sort(eigen(J, symmetric = TRUE, only.values = TRUE)$values)
  legendre <- function(x) {
    p0 <- 1; p1 <- x
    for (k in 2:n) {
      p2 <- ((2 * k - 1) * x * p1 - (k - 1) * p0) / k
      p0 <- p1; p1 <- p2
    }
    list(p = p1, dp = n * (x * p1 - p0) / (x^2 - 1))
  }
  for (it in 1:3) {
    l <- legendre(x)
    x <- x - l$p / l$dp
  }
  x <- (x - rev(x)) / 2  # exactly symmetric
  mtlba_gl_cache[[key]] <- list(x = x, w = 2 / ((1 - x^2) * legendre(x)$dp^2))
}

# Quadrature rule of the correlated kernel:
#   n      Gauss-Legendre nodes per panel;
#   hw     half-width of the integration window in sds of the stopping drift
#          (mtlba_quad_window);
#   panels number of equal panels, or NULL for one panel per pw widths of the
#          integrand's narrowest feature, sv_s * min(1, sqrt(1 - rho^2) / |rho|),
#          at most maxp.
# The default is the fit rule chosen in Stage M4 (rating-work/mtlba/stageM4/
# ACCURACY.md); mtlba_quad_reference() is the high-accuracy rule the tests and
# the accuracy map compare it with.
mtlba_quad_rule <- function(n = 16, hw = 8, panels = NULL, pw = 4, maxp = 64L)
  list(n = n, hw = hw, panels = panels, pw = pw, maxp = maxp)
mtlba_quad_reference <- function() mtlba_quad_rule(n = 16, hw = 12, pw = 1, maxp = 1024L)

# Integration window (a, b) of the stopping accumulator's drift x: its range
# ((bs - A_s) / dt, bs / dt) clipped to where the integrand is not negligible.
# The integrand is the joint (bivariate normal) density of x and the other
# drift y, summed over the y the band allows, y in (yl, yh) with
# yl = max(0, (T_lo - A_o) / dt), yh = T_hi / dt. Its mass in x sits between
# the marginal mean v_s and m = v_s + rho (sv_s / sv_o) (y* - v_o), the mean of
# x given the allowed y nearest to v_o (m = v_s when v_o is allowed; otherwise,
# e.g. v_o < 0 under the both-positive truncation, x is dragged towards m).
# The window keeps x within r of the interval (v_s, m), r^2 = d^2 +
# (hw sv_s)^2, d the distance from that interval to the range (0 when they
# overlap, giving +- hw sds), so the far tail keeps relative accuracy too.
# Under "any" drifts (pos = FALSE) the other drift is not truncated, so yl is
# the band's own lower limit, -Inf for the lowest band.
mtlba_quad_window <- function(dt, As, bs, vs, svs, Ao, vo, svo, lo_y, hi_y, lowest, rho, quad,
                              pos = TRUE) {
  lo_x <- (bs - As) / dt
  hi_x <- bs / dt
  yl <- if (pos) ifelse(lowest, 0, pmax(0, (lo_y - Ao) / dt)) else
    ifelse(lowest, -Inf, (lo_y - Ao) / dt)
  m <- vs + rho * svs / svo * (pmin(pmax(vo, yl), hi_y / dt) - vo)
  c_lo <- pmin(vs, m)
  c_hi <- pmax(vs, m)
  d <- pmax(lo_x - c_hi, c_lo - hi_x, 0)
  r <- sqrt(d^2 + (quad$hw * svs)^2)
  list(a = pmax(lo_x, c_lo - r), b = pmin(hi_x, c_hi + r))
}

# log of the correlated kernel times the normaliser (P(both drifts > 0) under
# "positive", P(not both < 0) under "any"): per trial the stopping accumulator
# (A_s, threshold bs, v_s, sv_s), the other accumulator (A_o, v_o, sv_o) and
# its band (lo_y, hi_y, lowest), decision time dt and rho. The integrand is
# formed on the log scale (band at the conditional mean and sd, unnormalised)
# and summed by log-sum-exp. A_s = 0: no integral, x = bs / dt.
mtlba_logkernel_rho <- function(dt, As, bs, vs, svs, Ao, vo, svo, lo_y, hi_y, lowest, rho,
                                quad = mtlba_quad_rule(), pos = TRUE) {
  n <- length(dt)
  out <- rep(-Inf, n)
  sc <- svo * sqrt(1 - rho^2)
  # log integrand at X (rows = the trials i, columns = nodes)
  lint <- function(i, X) {
    ii <- rep(i, ncol(X))
    x <- as.vector(X)
    mc <- vo[ii] + rho[ii] * svo[ii] / svs[ii] * (x - vs[ii])
    l <- log(x) + stats::dnorm(x, vs[ii], svs[ii], log = TRUE) +
      mtlba_logband(dt[ii], Ao[ii], mc, sc[ii], lo_y[ii], hi_y[ii], lowest[ii], normalise = FALSE,
                    pos = pos)
    l[is.na(l)] <- -Inf
    matrix(l, length(i), ncol(X))
  }
  small <- As <= 1e-10
  if (any(small)) {
    i <- which(small)
    out[i] <- lint(i, matrix(bs[i] / dt[i], ncol = 1))[, 1] - log(dt[i])
  }
  if (any(!small)) {
    win <- mtlba_quad_window(dt, As, bs, vs, svs, Ao, vo, svo, lo_y, hi_y, lowest, rho, quad, pos)
    width <- win$b - win$a
    P <- if (is.null(quad$panels)) {
      sig <- svs * pmin(1, sqrt(1 - rho^2) / abs(rho))
      pmin(quad$maxp, pmax(1, ceiling(width / (quad$pw * sig))))
    } else rep(as.integer(quad$panels), n)
    P[small | is.na(P)] <- 0
    g <- mtlba_gl(quad$n)
    for (p in setdiff(unique(P), 0)) {
      i <- which(P == p)
      h <- width[i] / p
      u <- as.vector(outer((g$x + 1) / 2, 0:(p - 1), "+"))
      L <- lint(i, win$a[i] + outer(h, u))
      L <- L + rep(log(rep(g$w, p) / 2), each = length(i))
      M <- L[cbind(seq_along(i), max.col(L, ties.method = "first"))]
      s <- M + log(rowSums(exp(L - M))) + log(h) - log(As[i])
      s[is.na(s) | M == -Inf] <- -Inf
      out[i] <- s
    }
  }
  out
}

# P(Z1 > h, Z2 > k) for standard bivariate normal (Z1, Z2) with correlation r,
# |r| < 1: Genz (2004, Statistics and Computing 14, 251-260) BVND, the
# double-precision algorithm (Gauss-Legendre in asin(r) for |r| < .925, the
# expansion about |r| = 1 above), absolute error ~1e-15.
mtlba_bvu <- function(h, k, r) {
  gl <- list(
    list(w = c(0.1713244923791705, 0.3607615730481384, 0.4679139345726904),
         x = c(-0.9324695142031522, -0.6612093864662647, -0.2386191860831970)),
    list(w = c(0.04717533638651177, 0.1069393259953183, 0.1600783285433464,
               0.2031674267230659, 0.2334925365383547, 0.2491470458134029),
         x = c(-0.9815606342467191, -0.9041172563704750, -0.7699026741943050,
               -0.5873179542866171, -0.3678314989981802, -0.1252334085114692)),
    list(w = c(0.01761400713915212, 0.04060142980038694, 0.06267204833410906,
               0.08327674157670475, 0.1019301198172404, 0.1181945319615184,
               0.1316886384491766, 0.1420961093183821, 0.1491729864726037,
               0.1527533871307259),
         x = c(-0.9931285991850949, -0.9639719272779138, -0.9122344282513259,
               -0.8391169718222188, -0.7463319064601508, -0.6360536807265150,
               -0.5108670019508271, -0.3737060887154196, -0.2277858511416451,
               -0.07652652113349733)))
  n <- length(h)
  out <- rep(NaN, n)
  ok <- is.finite(h) & is.finite(k) & is.finite(r) & abs(r) < 1
  grp <- ifelse(abs(r) < 0.3, 1L, ifelse(abs(r) < 0.75, 2L, 3L))
  for (gi in 1:3) {
    w <- rep(gl[[gi]]$w, 2)
    x <- c(gl[[gi]]$x, -gl[[gi]]$x)
    i <- which(ok & grp == gi & abs(r) < 0.925)
    if (length(i)) {
      hk <- h[i] * k[i]
      hs <- (h[i]^2 + k[i]^2) / 2
      asr <- asin(r[i])
      sn <- sin(outer(asr, (x + 1) / 2))
      out[i] <- as.vector(exp((sn * hk - hs) / (1 - sn^2)) %*% w) * asr / (4 * pi) +
        stats::pnorm(-h[i]) * stats::pnorm(-k[i])
    }
    i <- which(ok & grp == gi & abs(r) >= 0.925)
    if (length(i)) {
      hh <- h[i]
      neg <- r[i] < 0
      kk <- ifelse(neg, -k[i], k[i])
      hk <- hh * kk
      as <- (1 - r[i]) * (1 + r[i])
      a <- sqrt(as)
      bs <- (hh - kk)^2
      cc <- (4 - hk) / 8
      d <- (12 - hk) / 16
      asr <- -(bs / as + hk) / 2
      bvn <- ifelse(asr > -100,
                    a * exp(asr) * (1 - cc * (bs - as) * (1 - d * bs / 5) / 3 + cc * d * as^2 / 5), 0)
      b <- sqrt(bs)
      bvn <- bvn - ifelse(-hk < 100, exp(-hk / 2) * sqrt(2 * pi) * stats::pnorm(-b / a) * b *
                            (1 - cc * bs * (1 - d * bs / 5) / 3), 0)
      xs <- outer(a / 2, 1 + x)^2
      rs <- sqrt(1 - xs)
      asr <- -(bs / xs + hk) / 2
      term <- exp(asr) * (exp(-hk * (1 - rs) / (2 * (1 + rs))) / rs - (1 + cc * xs * (1 + d * xs)))
      term[asr <= -100] <- 0
      bvn <- -(bvn + (a / 2) * as.vector(term %*% w)) / (2 * pi)
      out[i] <- ifelse(neg, -bvn + ifelse(kk > hh, stats::pnorm(kk) - stats::pnorm(hh), 0),
                       bvn + stats::pnorm(-pmax(hh, kk)))
    }
  }
  out
}

# log P(both drifts > 0)
mtlba_logppos <- function(v1, sv1, v2, sv2, rho)
  log(pmax(mtlba_bvu(-v1 / sv1, -v2 / sv2, rho), 0))

# log P(not both drifts < 0) = log [P(v1 > 0) + P(v1 < 0, v2 > 0)], two
# positive terms (no 1 - P(both < 0) cancellation); the second is the upper
# orthant of (-Z1, Z2), correlation -rho.
mtlba_lognotneg <- function(v1, sv1, v2, sv2, rho)
  log(pmax(stats::pnorm(v1 / sv1) + mtlba_bvu(v1 / sv1, -v2 / sv2, -rho), 0))

# log normaliser of the drift convention, floored at log(1e-10)
mtlba_lognorm <- function(v1, sv1, v2, sv2, rho, pos = TRUE)
  pmax(if (pos) mtlba_logppos(v1, sv1, v2, sv2, rho) else mtlba_lognotneg(v1, sv1, v2, sv2, rho),
       log(1e-10))

# Both-positive correlated drifts: accumulator 1 from its positive-truncated
# marginal normal, accepted with probability proportional to
# P(v2 > 0 | v1) = pnorm(mc / sc), then accumulator 2 from its positive-truncated
# conditional normal.
mtlba_rdrift_rho <- function(m1, s1, m2, s2, rho) {
  n <- length(m1)
  v1 <- v2 <- rep(NA_real_, n)
  sc <- s2 * sqrt(1 - rho^2)
  # log of the largest acceptance probability over v1 > 0 (v1 = 0 when rho < 0)
  lM <- ifelse(rho > 0, 0, stats::pnorm((m2 - rho * s2 / s1 * m1) / sc, log.p = TRUE))
  todo <- which(is.finite(m1 + s1 + m2 + sc + lM))
  for (it in 1:2000) {
    if (!length(todo)) break
    x <- msm::rtnorm(length(todo), m1[todo], s1[todo], lower = 0)
    mc <- m2[todo] + rho[todo] * s2[todo] / s1[todo] * (x - m1[todo])
    acc <- log(stats::runif(length(todo))) + lM[todo] < stats::pnorm(mc / sc[todo], log.p = TRUE)
    if (any(acc)) {
      v1[todo[acc]] <- x[acc]
      v2[todo[acc]] <- msm::rtnorm(sum(acc), mc[acc], sc[todo[acc]], lower = 0)
    }
    todo <- todo[!acc]
  }
  if (length(todo)) warning("MTLBA: ", length(todo), " correlated drift pairs could not be sampled")
  list(v1 = v1, v2 = v2)
}

# Drift pairs under "any": bivariate normal (rho may be 0) by rejection of the
# pairs with both drifts negative.
mtlba_rdrift_any <- function(m1, s1, m2, s2, rho) {
  n <- length(m1)
  v1 <- v2 <- rep(NA_real_, n)
  sc <- s2 * sqrt(1 - rho^2)
  todo <- which(is.finite(m1 + s1 + m2 + sc) & abs(rho) < 1)
  for (it in 1:2000) {
    if (!length(todo)) break
    x <- stats::rnorm(length(todo), m1[todo], s1[todo])
    y <- stats::rnorm(length(todo), m2[todo] + rho[todo] * s2[todo] / s1[todo] * (x - m1[todo]),
                      sc[todo])
    acc <- x > 0 | y > 0
    v1[todo[acc]] <- x[acc]
    v2[todo[acc]] <- y[acc]
    todo <- todo[!acc]
  }
  if (length(todo)) warning("MTLBA: ", length(todo), " drift pairs could not be sampled")
  list(v1 = v1, v2 = v2)
}

# Response lookup: (choice R in 1/2, rating RR, K, rule) -> stopping
# accumulator s (1/2), its threshold j (1 .. K, K = b), and the other
# accumulator's band lo (0 .. K - 1; the band is (T_lo, T_{lo+1})).
# TC: responses in folded order p = 1 .. 2K (R = 1 rating K .. R = 2 rating K)
# are the stopping events ordered by accumulator 1's lead in crossings, ties
# going to the accumulator that did NOT stop: p odd = accumulator 1 stops at
# threshold K - (p - 1) / 2, p even = accumulator 2 stops at threshold p / 2.
mtlba_lookup <- function(R, RR, K, rule = "BE") {
  R <- as.integer(R)
  if (rule == "BE") {
    s <- R
    j <- rep(K, length(R))
  } else {
    p <- ifelse(R == 1L, K - RR + 1, K + RR)
    odd <- p %% 2 == 1
    s <- ifelse(odd, 1L, 2L)
    j <- ifelse(odd, K - (p - 1) / 2, p / 2)
  }
  list(s = as.integer(s), j = as.integer(j), lo = as.integer(K - if (rule == "BE") RR else j))
}

# Absolute thresholds (columns T_0 = 0 placeholder, T_1 .. T_{K-1}, T_K = b)
# from rows with columns A, b and natural-scale d1 .. d{K-1}.
mtlba_thresholds <- function(pars, K) {
  A <- pars[, "A"]
  b <- pars[, "b"]
  out <- matrix(0, nrow(pars), K + 1)
  if (K > 1) out[, 2:K] <- A + (b - A) * pars[, rating_threshold_names(K), drop = FALSE]
  out[, K + 1] <- b
  out
}

# Proportional parameterization -> natural additive columns, row-wise: pars
# holds B = b, A and c1 .. c{K-1} as proportions (after their transforms);
# returns pars with A the natural start-point range, b and the threshold
# proportions d1 .. d{K-1} of B above A (column B, the criteria and everything
# else unchanged). free: the lowest threshold nests on T_2 alone (K > 2; with
# K = 2 the one threshold is T_1 = b c1, above A).
mtlba_proportional <- function(pars, K, free) {
  ok <- attr(pars, "ok")
  b <- pars[, "B"]
  if (K == 1) {
    A <- b * pars[, "A"]
    d <- NULL
  } else {
    P <- pars[, rating_criterion_names(K), drop = FALSE]
    Ttop <- b * P[, 1]
    A <- Ttop * pars[, "A"]
    d <- matrix(0, nrow(pars), K - 1, dimnames = list(NULL, rating_threshold_names(K)))
    d[, K - 1] <- (Ttop - A) / (b - A)
    if (K > 2) for (j in (K - 2):1) {
      p <- P[, K - j]
      d[, j] <- if (j == 1 && free) ((A + (b - A) * d[, 2]) * p - A) / (b - A) else d[, j + 1] * p
    }
  }
  pars[, "A"] <- A
  pars <- cbind(pars[, !(colnames(pars) %in% c("b", colnames(d))), drop = FALSE], b = b, d)
  if (!is.null(ok)) attr(pars, "ok") <- ok
  pars
}

# log kernel per trial for explicit thresholds: accumulator s (1/2) reaches bs
# at dt = rt - t0_s while the other accumulator's evidence lies in (lo_y, hi_y)
# (lowest: no threshold below, lo_y ignored). rt per trial; pars two rows per
# trial (accumulators in level order). An empty band (hi_y <= lo_y) gives -Inf.
# pos: the drift convention ("positive" TRUE, "any" FALSE).
mtlba_logkernel_y <- function(rt, pars, s, bs, lo_y, hi_y, lowest, quad = mtlba_quad_rule(),
                              pos = TRUE) {
  n <- length(rt)
  first <- seq(1, 2 * n, by = 2)
  si <- first + s - 1L
  oi <- first + 2L - s
  dt <- rt - pars[si, "t0"]
  rho <- if ("rho" %in% colnames(pars)) pars[first, "rho"] else rep(0, n)
  out <- rep(-Inf, n)
  # non-finite parameters (e.g. an overflowing proposal) give -Inf, never an error
  fin <- function(i, cols) rowSums(!is.finite(pars[i, cols, drop = FALSE])) == 0
  ok <- !is.na(dt) & is.finite(dt) & dt > 0 & fin(si, c("v", "sv", "A", "b")) &
    fin(oi, c("v", "sv", "A", "b")) & is.finite(bs) & is.finite(hi_y) &
    (lowest | (is.finite(lo_y) & hi_y > lo_y)) & is.finite(rho) & abs(rho) < 1
  ok[is.na(ok)] <- FALSE
  if (!any(ok)) return(out)
  si <- si[ok]; oi <- oi[ok]; dt <- dt[ok]; rho <- rho[ok]
  bs <- bs[ok]; lo_y <- lo_y[ok]; hi_y <- hi_y[ok]; lowest <- lowest[ok]
  res <- numeric(length(dt))
  u <- rho == 0
  if (any(u)) {
    res[u] <- mtlba_logf(dt[u], pars[si[u], "A"], bs[u], pars[si[u], "v"], pars[si[u], "sv"], pos) +
      mtlba_logband(dt[u], pars[oi[u], "A"], pars[oi[u], "v"], pars[oi[u], "sv"],
                    lo_y[u], hi_y[u], lowest[u], pos = pos)
    # "any": the closed form is a product of untruncated terms over P(not both < 0)
    if (!pos) res[u] <- res[u] - mtlba_lognorm(pars[si[u], "v"], pars[si[u], "sv"],
                                               pars[oi[u], "v"], pars[oi[u], "sv"], rho[u], pos)
  }
  if (any(!u)) {
    k <- which(!u)
    res[k] <- mtlba_logkernel_rho(dt[k], pars[si[k], "A"], bs[k], pars[si[k], "v"],
                                  pars[si[k], "sv"], pars[oi[k], "A"], pars[oi[k], "v"],
                                  pars[oi[k], "sv"], lo_y[k], hi_y[k], lowest[k], rho[k], quad,
                                  pos) -
      mtlba_lognorm(pars[si[k], "v"], pars[si[k], "sv"], pars[oi[k], "v"], pars[oi[k], "sv"],
                    rho[k], pos)
  }
  out[ok] <- res
  out
}

# log kernel per trial for given lookup indices (stopping accumulator s, its
# threshold j, the other's band lo); Tm = mtlba_thresholds(pars, K).
mtlba_logkernel <- function(rt, pars, Tm, s, j, lo, quad = mtlba_quad_rule(), pos = TRUE) {
  first <- seq(1, 2 * length(rt), by = 2)
  si <- first + s - 1L
  oi <- first + 2L - s
  mtlba_logkernel_y(rt, pars, s, Tm[cbind(si, j + 1L)], Tm[cbind(oi, lo + 1L)],
                    Tm[cbind(oi, lo + 2L)], lo == 0L, quad, pos)
}

# log(exp(a) + exp(b))
mtlba_lse <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log1p(exp(-abs(a - b)))
  out[m == -Inf] <- -Inf
  out
}

# ---- confidence guessing ----
# Each accumulator has a guess threshold gb = A + B * gB (0 <= gB <= 1, so
# A <= gb <= b). An accumulator that is TRIGGERING on a trial ends the race when
# it reaches gb: it is then the choice and the rating is drawn at random,
# rating k with probability pg_k prod_{i<k} (1 - pg_i) (k < K; the highest
# rating takes the rest), the pg from the choice's row. Otherwise the race ends
# by the evidence rule (BE or TC), read as usual.
#   guess = "accumulator" (the manuscript's model, first-submission supplement
#     section 1.3): the two accumulators trigger independently with
#     probabilities gp_1, gp_2 (each from its own row). With g_a = gp_a, s the
#     stopping accumulator of the evidence response at threshold T_s, o the
#     other one and its band (lo, hi), w the choice of a guessed response:
#       evidence part = cs * [(1 - g_o) k(T_s; lo, hi) + g_o k(T_s; lo, min(hi, gb_o))],
#         cs = 1 - g_s (BE: a triggering accumulator never reaches b), TC:
#         cs = 1 if T_s < gb_s (it stops by count below its guess threshold);
#       guess part = g_w pg(RR) * [(1 - g_o) k(gb_w; 0, U) + g_o k(gb_w; 0, min(U, gb_o))],
#         U = b_o (BE); TC: U = T_{o, K - n_w}, n_w the number of w's thresholds
#         below gb_w (the count must still be below K when w reaches gb_w).
#   guess = "trial" (BETE): the whole trial is a guess with probability gp
#     (trial-level): then both accumulators race to their gb, under either rule;
#       likelihood = (1 - gp) k(T_s; lo, hi) + gp pg(RR) k(gb_w; 0, gb_o).
# k() is the kernel above (correlated when rho != 0, guesses included). Rows
# with gp = 0 on both accumulators take the plain kernel.
#
# gb_constraint (the manuscript's supplement, equations 17 and 18, used there for
# one data set): in the limit where an accumulator always triggers its choice
# threshold b is never reached, so the data fix gb but not b. The constraint
# moves the guess threshold towards b as triggering becomes frequent,
#   gb = b - B (1 - gB) (1 - g^2) = A + B (gB + g^2 (1 - gB)),
# g the accumulator's own gp ("own") or the larger of the two accumulators' gp
# ("max", the supplement's rule); "none": gb = A + B gB.

# Guessing quantities per parameter row, or NULL when no row guesses: g (trigger
# probability), gb (guess threshold), lpg (log rating-guess probabilities,
# columns = ratings 1 .. K). pars: natural-scale gp, gB, pg1 .. pg{K-1}
# (absent: gB = .5, pg1 = 1, further pg = 0) besides A and b.
mtlba_guess_pars <- function(pars, K, guess = "accumulator", gb_constraint = "none") {
  if (!("gp" %in% colnames(pars))) return(NULL)
  g <- pars[, "gp"]
  if (guess == "trial") g <- rep(g[seq(1, length(g), by = 2)], each = 2)
  if (!anyNA(g) && all(g == 0)) return(NULL)
  n <- nrow(pars)
  gB <- if ("gB" %in% colnames(pars)) pars[, "gB"] else rep(0.5, n)
  gB <- mtlba_gB_constrained(gB, g, gb_constraint)
  lpg <- matrix(0, n, K)
  if (K > 1) {
    lrest <- rep(0, n)
    for (k in 1:(K - 1)) {
      nm <- paste0("pg", k)
      pk <- if (nm %in% colnames(pars)) pars[, nm] else rep(as.numeric(k == 1), n)
      lpg[, k] <- lrest + log(pk)
      lrest <- lrest + log1p(-pk)
    }
    lpg[, K] <- lrest
  }
  list(g = g, gb = pars[, "A"] + (pars[, "b"] - pars[, "A"]) * gB, lpg = lpg)
}

# Binomial guess ratings (guess_ratings = "binomial"): a guessed rating is
# 1 + Binomial(K - 1, q), P(k) = choose(K - 1, k - 1) q^(k-1) (1 - q)^(K-k),
# returned as the sequential pg_k = P(k) / P(>= k), k = 1 .. K - 1 (columns
# pg1 .. pg{K-1}), the form mtlba_guess_pars and rMTLBA read. On the log scale
# (the tails P(>= k) as log-sums of their terms) so q near 0 or 1 loses
# nothing; where P(>= k) is 0 the conditional is irrelevant and set to 1. The
# compiled twin is mtlba_binomial_pg in src/model_mtlba.h.
mtlba_binomial_pg <- function(q, K) {
  lq <- log(q); l1q <- log1p(-q)
  lP <- sapply(1:K, function(k) log(choose(K - 1, k - 1)) + (if (k > 1) (k - 1) * lq else 0) +
                 (if (k < K) (K - k) * l1q else 0))
  lP <- matrix(lP, length(q))
  pg <- matrix(1, length(q), K - 1, dimnames = list(NULL, paste0("pg", seq_len(K - 1))))
  tail <- lP[, K]
  for (k in (K - 1):1) {
    tail <- mtlba_lse(lP[, k], tail)
    pk <- exp(lP[, k] - tail)
    pk[which(tail == -Inf)] <- 1
    pg[, k] <- pk
  }
  pg
}

# Guess-threshold proportion under gb_constraint: gB + gc^2 (1 - gB), gc the
# row's own trigger probability g ("own") or the larger of its trial's two
# ("max"); g has two rows per trial.
mtlba_gB_constrained <- function(gB, g, gb_constraint = "none") {
  if (gb_constraint == "none") return(gB)
  gc <- g
  if (gb_constraint == "max") {
    first <- seq(1, length(g), by = 2)
    gc <- rep(pmax(g[first], g[first + 1]), each = 2)
  }
  gB + gc^2 * (1 - gB)
}

# log(w1 exp(l1) + w2 exp(l2)) for weights w1, w2 >= 0
mtlba_lmix <- function(w1, l1, w2, l2) {
  a <- ifelse(w1 > 0, log(w1) + l1, -Inf)
  b <- ifelse(w2 > 0, log(w2) + l2, -Inf)
  a[is.na(a)] <- -Inf
  b[is.na(b)] <- -Inf
  mtlba_lse(a, b)
}

# log density per trial of the response (R, RR), no don't-know merging. G from
# mtlba_guess_pars (NULL: no guessing).
mtlba_logresp <- function(rt, R, RR, pars, Tm, K, rule, quad, G = NULL, guess = "accumulator",
                          pos = TRUE) {
  n <- length(rt)
  lk <- mtlba_lookup(R, RR, K, rule)
  ld <- mtlba_logkernel(rt, pars, Tm, lk$s, lk$j, lk$lo, quad, pos)
  if (is.null(G)) return(ld)
  first <- seq(1, 2 * n, by = 2)
  act <- which(is.na(G$g[first]) | is.na(G$g[first + 1]) | G$g[first] != 0 | G$g[first + 1] != 0)
  if (!length(act)) return(ld)
  # the guessing trials only
  rows <- as.vector(rbind(2 * act - 1, 2 * act))
  P <- pars[rows, , drop = FALSE]
  TT <- Tm[rows, , drop = FALSE]
  g <- G$g[rows]; gb <- G$gb[rows]; lpg <- G$lpg[rows, , drop = FALSE]
  m <- length(act)
  f1 <- seq(1, 2 * m, by = 2)
  rtt <- rt[act]; Rw <- as.integer(R)[act]; RRw <- RR[act]
  s <- lk$s[act]; j <- lk$j[act]; lo <- lk$lo[act]
  si <- f1 + s - 1L; oi <- f1 + 2L - s
  bs <- TT[cbind(si, j + 1L)]
  lo_y <- TT[cbind(oi, lo + 1L)]; hi_y <- TT[cbind(oi, lo + 2L)]; lowest <- lo == 0L
  kern <- function(s, bs, lo_y, hi_y, lowest, w = rep(TRUE, m)) {
    out <- rep(-Inf, m)
    w <- w & !is.na(w)
    if (any(w)) out[w] <- mtlba_logkernel_y(rtt[w], P[rep(w, each = 2), , drop = FALSE], s[w], bs[w],
                                            lo_y[w], hi_y[w], lowest[w], quad, pos)
    out
  }
  full <- ld[act]
  wi <- f1 + Rw - 1L  # the choice's row: the accumulator that guesses
  wo <- f1 + 2L - Rw
  lpgc <- lpg[cbind(wi, RRw)]
  low <- rep(TRUE, m)
  zero <- rep(0, m)
  if (guess == "trial") {
    ev <- mtlba_lmix(1 - g[f1], full, zero, full)
    gk <- kern(Rw, gb[wi], zero, gb[wo], low)
    gu <- mtlba_lmix(g[f1], lpgc + gk, zero, gk)
  } else {
    cs <- 1 - g[si]
    if (rule == "TC") cs <- ifelse(!is.na(bs) & bs < gb[si], 1, cs)
    cap <- pmin(hi_y, gb[oi])
    capped <- !is.na(cap) & cap < hi_y
    tr <- full
    tr[capped] <- kern(s, bs, lo_y, cap, lowest, capped & g[oi] > 0)[capped]
    ev <- mtlba_lmix(cs * (1 - g[oi]), full, cs * g[oi], tr)
    if (rule == "BE") U0 <- P[wo, "b"] else {
      nw <- rowSums(TT[wi, -1, drop = FALSE] < gb[wi])
      U0 <- TT[cbind(wo, K - nw + 1L)]
    }
    U1 <- pmin(U0, gb[wo])
    capped <- !is.na(U1) & U1 < U0
    gk0 <- kern(Rw, gb[wi], zero, U0, low, g[wi] > 0 & (g[wo] < 1 | !capped))
    gk1 <- gk0
    gk1[capped] <- kern(Rw, gb[wi], zero, U1, low, capped & g[wi] > 0 & g[wo] > 0)[capped]
    gu <- lpgc + mtlba_lmix(g[wi] * (1 - g[wo]), gk0, g[wi] * g[wo], gk1)
  }
  out <- mtlba_lse(ev, gu)
  out[is.na(out)] <- -Inf
  ld[act] <- out
  ld
}

# log likelihood per trial: rt, R (1/2), RR per trial, pars two rows per trial.
mtlba_logdens <- function(rt, R, RR, pars, K, rule, dont_know, quad = mtlba_quad_rule(),
                          guess = "accumulator", gb_constraint = "none", pos = TRUE) {
  Tm <- mtlba_thresholds(pars, K)
  G <- mtlba_guess_pars(pars, K, guess, gb_constraint)
  R <- as.integer(R)
  dk <- if (dont_know) !is.na(RR) & RR == 1 else rep(FALSE, length(R))
  # a don't-know trial is the sum of the two rating-1 terms: R = 1 here, R = 2 below
  R[dk] <- 1L
  ld <- mtlba_logresp(rt, R, RR, pars, Tm, K, rule, quad, G, guess, pos)
  if (any(dk)) {
    w <- which(dk)
    idx <- as.vector(rbind(2 * w - 1, 2 * w))
    G2 <- if (is.null(G)) NULL else list(g = G$g[idx], gb = G$gb[idx], lpg = G$lpg[idx, , drop = FALSE])
    ld2 <- mtlba_logresp(rt[w], rep(2L, length(w)), rep(1, length(w)), pars[idx, , drop = FALSE],
                         Tm[idx, , drop = FALSE], K, rule, quad, G2, guess, pos)
    ld[w] <- mtlba_lse(ld[w], ld2)
  }
  ld[is.na(ld)] <- -Inf
  ld
}

# The drift convention as a flag: TRUE for "positive", FALSE for "any".
mtlba_pos <- function(drifts) match.arg(drifts, c("any", "positive")) == "positive"

# Density of the MTLBA. rt, R (factor or integer 1/2) and RR (1 = lowest) per
# trial; pars has two rows per trial (accumulators in level order) with columns
# v, sv, A, b, t0 and natural-scale d1 .. d{K-1} (proportions of B = b - A
# above A; d1 may be negative, a lowest threshold below A), and optionally rho
# (drift-rate correlation, read from the trial's first row; 0 if absent).
# Don't-know trials: RR = 1 (R is ignored). quad: the quadrature rule used
# where rho != 0. Confidence guessing: natural-scale columns gp, gB and
# pg1 .. pg{K-1} (see mtlba_guess_pars; no gp column or gp = 0: none).
# drifts: "any" (conditioned on not both negative) or "positive" (both > 0).
dMTLBA <- function(rt, R, RR, pars, n_ratings = 3, rule = "BE", dont_know = FALSE,
                   log = FALSE, quad = mtlba_quad_rule(), guess = "accumulator",
                   gb_constraint = "none", drifts = "any") {
  K <- rating_check_n(n_ratings)
  rule <- match.arg(rule, c("BE", "TC"))
  guess <- match.arg(guess, c("accumulator", "trial"))
  rating_check_dont_know(K, dont_know)
  gb_constraint <- match.arg(gb_constraint, c("none", "own", "max"))
  pos <- mtlba_pos(drifts)
  ld <- mtlba_logdens(rt, R, RR, pars, K, rule, dont_know, quad, guess, gb_constraint, pos)
  if (log) ld else exp(ld)
}

# Simulate the MTLBA: lR the accumulator factor (two rows per trial), pars as
# for dMTLBA. Returns R, rt and RR; don't-know responses are written as RR = 1
# with R = the first level. Drifts: "positive" draws each accumulator's drift
# from its positive-truncated normal, trials with rho != 0 from the bivariate
# normal truncated to both positive (mtlba_rdrift_rho); "any" draws the pair
# from the (bivariate) normal, rejecting pairs with both drifts negative, so an
# accumulator with a negative drift never crosses anything (crossing time
# Inf) and a loser may end below its start point. Confidence guessing as for
# dMTLBA: triggering accumulators end the race at their guess threshold with a
# random rating.
rMTLBA <- function(lR, pars, n_ratings = 3, rule = "BE", dont_know = FALSE,
                   ok = rep(TRUE, nrow(pars)), guess = "accumulator", gb_constraint = "none",
                   drifts = "any") {
  K <- rating_check_n(n_ratings)
  rule <- match.arg(rule, c("BE", "TC"))
  guess <- match.arg(guess, c("accumulator", "trial"))
  rating_check_dont_know(K, dont_know)
  pos <- mtlba_pos(drifts)
  if (length(levels(lR)) != 2) stop("MTLBA simulates two accumulators only")
  if (is.null(ok)) ok <- rep(TRUE, nrow(pars))
  n <- nrow(pars) / 2
  first <- seq(1, 2 * n, by = 2)
  Tm <- mtlba_thresholds(pars, K)
  rho <- if ("rho" %in% colnames(pars)) pars[first, "rho"] else rep(0, n)
  # a non-positive drift never reaches a threshold ("any")
  tau_of <- function(T, start, v) {
    tau <- (T - start) / v
    tau[!is.na(v) & v <= 0, ] <- Inf
    tau
  }
  # crossing times of thresholds 1 .. K, one n x K matrix per accumulator
  cross <- lapply(0:1, function(a) {
    p <- pars[first + a, , drop = FALSE]
    start <- p[, "A"] * stats::runif(n)
    # "any": every trial's pair is drawn below
    v <- if (pos) msm::rtnorm(n, p[, "v"], p[, "sv"], lower = 0) else rep(NA_real_, n)
    list(tau = tau_of(Tm[first + a, -1, drop = FALSE], start, v), start = start, v = v)
  })
  cr <- which(is.na(rho) | rho != 0 | !pos)
  if (length(cr)) {
    vv <- if (pos) mtlba_rdrift_rho(pars[first[cr], "v"], pars[first[cr], "sv"],
                                    pars[first[cr] + 1, "v"], pars[first[cr] + 1, "sv"], rho[cr]) else
      mtlba_rdrift_any(pars[first[cr], "v"], pars[first[cr], "sv"],
                       pars[first[cr] + 1, "v"], pars[first[cr] + 1, "sv"], rho[cr])
    for (a in 1:2) {
      cross[[a]]$v[cr] <- vv[[a]]
      cross[[a]]$tau[cr, ] <- tau_of(Tm[first[cr] + a - 1, -1, drop = FALSE],
                                     cross[[a]]$start[cr], vv[[a]])
    }
  }
  t1 <- cross[[1]]$tau
  t2 <- cross[[2]]$tau
  # confidence guessing: which accumulators trigger, and when they reach gb
  G <- mtlba_guess_pars(pars, K, guess, match.arg(gb_constraint, c("none", "own", "max")))
  tg <- matrix(Inf, n, 2)
  if (!is.null(G)) {
    trig <- if (guess == "trial") rep(stats::runif(n) < G$g[first], 2) else
      stats::runif(2 * n) < c(G$g[first], G$g[first + 1])
    trig <- matrix(trig & !is.na(trig), n, 2)
    for (a in 1:2) tg[trig[, a], a] <-
        tau_of(cbind(G$gb[first + a - 1]), cross[[a]]$start, cross[[a]]$v)[trig[, a], 1]
    if (guess == "trial") {
      # a guess trial is a plain race to the guess thresholds
      gt <- trig[, 1]
      t1[gt, ] <- Inf
      t2[gt, ] <- Inf
    } else if (rule == "BE") {
      # a triggering accumulator stops at gb, before it can reach b
      t1[trig[, 1], K] <- Inf
      t2[trig[, 2], K] <- Inf
    }
  }
  if (rule == "BE") {
    win1 <- t1[, K] < t2[, K]
    dt <- ifelse(win1, t1[, K], t2[, K])
    # other accumulator's evidence at dt, counted against its thresholds
    e <- ifelse(win1, cross[[2]]$start + cross[[2]]$v * dt, cross[[1]]$start + cross[[1]]$v * dt)
    To <- Tm[ifelse(win1, first + 1, first), , drop = FALSE]
    RR <- if (K > 1) K - rowSums(To[, 2:K, drop = FALSE] < e) else rep(1, n)
    R <- ifelse(win1, 1L, 2L)
    s <- R
  } else {
    # stop at the K-th crossing overall: min over j of max(t1[j], t2[K - j])
    t1z <- cbind(0, t1)
    t2z <- cbind(0, t2)
    stop_t <- rep(Inf, n)
    for (jj in 0:K) stop_t <- pmin(stop_t, pmax(t1z[, jj + 1], t2z[, K - jj + 1]))
    j1 <- rowSums(t1 <= stop_t)
    s1 <- j1 > 0 & t1z[cbind(seq_len(n), j1 + 1)] == stop_t
    s <- ifelse(s1, 1L, 2L)
    j <- ifelse(s1, j1, K - j1)
    p <- ifelse(s1, 2 * (K - j + 1) - 1, 2 * j)
    R <- ifelse(p <= K, 1L, 2L)
    RR <- ifelse(p <= K, K - p + 1, p - K)
    dt <- stop_t
  }
  if (!is.null(G)) {
    # guessed responses: a trigger reached before the evidence rule stops the race
    tgm <- pmin(tg[, 1], tg[, 2])
    gs <- which(!is.na(tgm) & !is.na(dt) & tgm < dt)
    if (length(gs)) {
      a <- ifelse(tg[gs, 1] <= tg[gs, 2], 1L, 2L)
      cp <- t(apply(exp(G$lpg[first[gs] + a - 1, , drop = FALSE]), 1, cumsum))
      if (K == 1) cp <- matrix(1, length(gs), 1)
      RR[gs] <- pmin(K, 1 + rowSums(stats::runif(length(gs)) > cp))
      R[gs] <- a
      s[gs] <- a
      dt[gs] <- tgm[gs]
    }
  }
  if (dont_know) R[RR == 1] <- 1L
  rt <- pars[first + s - 1, "t0"] + dt
  tok <- matrix(ok, nrow = 2)
  tok <- tok[1, ] & tok[2, ]
  rt[!tok] <- NA
  RR[!tok] <- NA
  R[!tok] <- NA
  data.frame(R = factor(levels(lR)[R], levels = levels(lR)), rt = rt, RR = RR)
}

# log P(dt <= y) and log P(dt > y) for the stopping time of the race, y a
# decision-time bound per trial (y <= 0: nothing has stopped; y = Inf: all
# has). With n_a the number of thresholds accumulator a has crossed at y
# (n_a = K: past b), F_a(n) = P(n_a <= n), Fc_a(n) = 1 - F_a(n), p_a(n) =
# P(n_a = n), all given v_a > 0:
#   BE: survivor = F_1(K-1) F_2(K-1), stopped = Fc_1(K-1) + F_1(K-1) Fc_2(K-1);
#   TC: survivor = sum_{n1 < K} p_1(n1) F_2(K-1-n1),
#       stopped  = Fc_1(K-1) + sum_{n1 < K} p_1(n1) Fc_2(K-1-n1);
# every term is positive, so both are formed on the log scale without
# cancellation. F_a is the lowest band up to T_{n+1} (mtlba_logband) and
# Fc_a(n) = S(T_{n+1}) / Phi(v/sv) (thresholds above A: no negative-drift mass
# above them; a lowest threshold below A takes Sp of mtlba_logHSp). The
# response kernels summed over responses are the density of this stopping
# time. i1, i2: rows of accumulators 1 and 2; Tm from mtlba_thresholds.
# Under "any" drifts (pos = FALSE) the same positive-drift terms, unnormalised
# (Hp = P(e < T, v > 0), S = P(e >= T)), are completed by the trials on which
# one drift is negative (N_a = P(v_a < 0); both negative is excluded):
#   stopped  = [stopped_pos + N_1 S_2(T_{2,K}) + N_2 S_1(T_{1,K})] / P(not both < 0)
#   survivor = [survivor_pos + N_1 Hp_2(T_{2,K}) + N_2 Hp_1(T_{1,K})] / P(not both < 0)
# (an accumulator with a negative drift crosses nothing, so with the other
# one positive the race is decided by that one alone), still all positive.
mtlba_logstop <- function(y, pars, Tm, i1, i2, K, rule, pos = TRUE) {
  n <- length(y)
  lP <- rep(-Inf, n)
  lS <- rep(0, n)
  inf <- !is.na(y) & y == Inf
  lP[inf] <- 0
  lS[inf] <- -Inf
  bad <- is.na(y)
  lP[bad] <- lS[bad] <- NaN
  w <- !bad & !inf & y > 0
  if (!any(w)) return(list(P = lP, S = lS))
  acc <- function(i) {
    i <- i[w]
    t <- y[w]
    A <- pars[i, "A"]; v <- pars[i, "v"]; sv <- pars[i, "sv"]
    lPhi <- stats::pnorm(v / sv, log.p = TRUE)
    lden <- if (pos) pmax(lPhi, log(1e-10)) else rep(0, length(t))
    m <- length(t)
    F <- Fc <- band <- matrix(-Inf, m, K)
    for (nn in 0:(K - 1)) {
      hi <- Tm[cbind(i, nn + 2L)]
      F[, nn + 1] <- mtlba_logband(t, A, v, sv, rep(0, m), hi, rep(TRUE, m), normalise = pos)
      Fc[, nn + 1] <- mtlba_logHSp(hi, t, A, v, sv)$S - lden
      band[, nn + 1] <- if (nn == 0) F[, 1] else
        mtlba_logband(t, A, v, sv, Tm[cbind(i, nn + 1L)], hi, rep(FALSE, m), normalise = pos)
    }
    list(F = F, Fc = Fc, band = band, lPhi = lPhi, lN = stats::pnorm(-v / sv, log.p = TRUE))
  }
  a1 <- acc(i1)
  a2 <- acc(i2)
  # "positive": terms normalised by Phi(v/sv) per accumulator; "any": the
  # positive-drift part of P(stopped, not both < 0), unnormalised, with
  # accumulator 2's positivity explicit in the first term
  l2 <- if (pos) 0 else a2$lPhi
  if (rule == "BE") {
    S <- a1$F[, K] + a2$F[, K]
    P <- mtlba_lse(a1$Fc[, K] + l2, a1$F[, K] + a2$Fc[, K])
  } else {
    S <- a1$band[, 1] + a2$F[, K]
    P <- mtlba_lse(a1$Fc[, K] + l2, a1$band[, 1] + a2$Fc[, K])
    if (K > 1) for (n1 in 1:(K - 1)) {
      S <- mtlba_lse(S, a1$band[, n1 + 1] + a2$F[, K - n1])
      P <- mtlba_lse(P, a1$band[, n1 + 1] + a2$Fc[, K - n1])
    }
  }
  if (!pos) {
    lZ <- mtlba_lse(a1$lPhi, a1$lN + a2$lPhi)   # log P(not both < 0)
    S <- mtlba_lse(mtlba_lse(S, a1$lN + a2$F[, K]), a2$lN + a1$F[, K]) - lZ
    P <- mtlba_lse(mtlba_lse(P, a1$lN + a2$Fc[, K]), a2$lN + a1$Fc[, K]) - lZ
  }
  lP[w] <- P
  lS[w] <- S
  list(P = lP, S = lS)
}

# Correlated twin of mtlba_logstop for ONE parameter row and bound y (finite,
# > 0): c(log P(dt <= y), log P(dt > y)). One integral over accumulator 1's
# drift x > 0 (its marginal normal): given x, accumulator 1 has crossed n1
# thresholds with probability w_n1(x) (its uniform start point: piecewise
# linear in x, kinks at (T - A_1) / y and T / y for each threshold T), and
# accumulator 2 has its conditional normal drift, so with low_2(m | x) =
# P(v2 > 0, e_2 < T_{2,m} | x), up_2(m | x) = P(e_2 >= T_{2,m} | x) and
# pos_2(x) = P(v2 > 0 | x) the integrands are
#   BE: survivor w_{<K} low_2(K), stopped w_K pos_2 + w_{<K} up_2(K);
#   TC: survivor sum_{n1 < K} w_n1 low_2(K - n1),
#       stopped  w_K pos_2 + sum_{n1 < K} w_n1 up_2(K - n1),
# all positive terms, divided by P(both > 0). Composite Gauss-Legendre between
# the kinks, panels no wider than quad$pw of the narrowest feature. The range
# is clipped to hw sds around the drift mean, so the terms are accurate
# absolutely (~1e-14), not relatively, in the far tail.
# T1, T2: thresholds 1 .. K (T_K = b) of the two accumulators.
# Under "any" (pos = FALSE) x runs over the whole real line: for x > 0 the
# other drift is unrestricted (low_2 = P(e_2 < T | x), pos_2 = 1), for x < 0
# it must be positive (low_2 as under "positive"; w_K = 0 there), and the
# integrand jumps at x = 0 (a break point); divided by P(not both < 0).
mtlba_stop_rho1 <- function(y, A1, v1, s1, T1, A2, v2, s2, T2, rho, K, rule,
                            quad = mtlba_quad_reference(), pos = TRUE) {
  if (!all(is.finite(c(y, A1, v1, s1, T1, A2, v2, s2, T2, rho))) || abs(rho) >= 1 ||
      s1 <= 0 || s2 <= 0) return(c(NaN, NaN))
  # window as mtlba_quad_window: v2 > 0 drags x towards m when v2's mean is negative
  m <- v1 + rho * s1 / s2 * (max(v2, 0) - v2)
  r <- sqrt((if (pos) max(-max(v1, m), 0) else 0)^2 + (quad$hw * s1)^2)
  a <- if (pos) max(0, min(v1, m) - r) else min(v1, m) - r
  b <- max(v1, m) + r
  kinks <- c((T1 - A1) / y, T1 / y, if (!pos) 0)
  brk <- sort(unique(c(a, b, kinks[kinks > a & kinks < b])))
  sig <- s1 * min(1, sqrt(1 - rho^2) / abs(rho))
  g <- mtlba_gl(quad$n)
  x <- w <- vector("list", length(brk) - 1)
  for (i in seq_len(length(brk) - 1)) {
    width <- brk[i + 1] - brk[i]
    p <- min(quad$maxp, max(1, ceiling(width / (quad$pw * sig))))
    h <- width / p
    x[[i]] <- brk[i] + h * as.vector(outer((g$x + 1) / 2, 0:(p - 1), "+"))
    w[[i]] <- h * rep(g$w, p) / 2
  }
  x <- unlist(x)
  w <- unlist(w) * stats::dnorm(x, v1, s1)
  nx <- length(x)
  mc <- v2 + rho * s2 / s1 * (x - v1)
  sc <- rep(s2 * sqrt(1 - rho^2), nx)
  yy <- rep(y, nx)
  AA <- rep(A2, nx)
  # P(e_1 < T_{1,n} | x), n = 1 .. K
  Fle <- sapply(1:K, function(n)
    if (A1 <= 1e-10) as.numeric(T1[n] > x * y) else pmin(pmax((T1[n] - x * y) / A1, 0), 1))
  # P(e_2 < T_{2,m}, v_2 > 0 | x) (positive-drift form; under "any" for x < 0)
  low2 <- sapply(1:K, function(m)
    exp(mtlba_logband(yy, AA, mc, sc, rep(0, nx), rep(T2[m], nx), rep(TRUE, nx), normalise = FALSE)))
  if (pos) {
    up2 <- sapply(1:K, function(m) exp(mtlba_logHSp(rep(T2[m], nx), yy, AA, mc, sc)$S))
    pos2 <- stats::pnorm(mc / sc)
  } else {
    # x > 0: the other drift is unrestricted
    up2 <- sapply(1:K, function(m) exp(mtlba_logHS(rep(T2[m], nx), yy, AA, mc, sc)$S))
    lowu <- sapply(1:K, function(m) exp(mtlba_logHS(rep(T2[m], nx), yy, AA, mc, sc)$H))
    dim(low2) <- dim(lowu) <- c(nx, K)
    low2[x > 0, ] <- lowu[x > 0, ]
    pos2 <- rep(1, nx)
  }
  dim(Fle) <- dim(low2) <- dim(up2) <- c(nx, K)
  if (rule == "BE") {
    Sint <- Fle[, K] * low2[, K]
    Pint <- (1 - Fle[, K]) * pos2 + Fle[, K] * up2[, K]
  } else {
    wn <- Fle - cbind(0, Fle[, -K, drop = FALSE])  # w_n1, n1 = 0 .. K - 1
    Sint <- rowSums(wn * low2[, K:1, drop = FALSE])
    Pint <- (1 - Fle[, K]) * pos2 + rowSums(wn * up2[, K:1, drop = FALSE])
  }
  log(c(sum(w * Pint), sum(w * Sint))) -
    if (pos) mtlba_logppos(v1, s1, v2, s2, rho) else mtlba_lognotneg(v1, s1, v2, s2, rho)
}

# mtlba_logstop for any rho: rows with rho = 0 take the closed form, the others
# mtlba_stop_rho1, evaluated once per distinct (bound, parameter row).
mtlba_logstop_any <- function(y, pars, Tm, i1, i2, K, rule, pos = TRUE) {
  rho <- if ("rho" %in% colnames(pars)) pars[i1, "rho"] else rep(0, length(i1))
  cr <- is.na(rho) | rho != 0
  if (!any(cr)) return(mtlba_logstop(y, pars, Tm, i1, i2, K, rule, pos))
  n <- length(y)
  lP <- rep(-Inf, n)
  lS <- rep(0, n)
  if (any(!cr)) {
    o <- mtlba_logstop(y[!cr], pars, Tm, i1[!cr], i2[!cr], K, rule, pos)
    lP[!cr] <- o$P
    lS[!cr] <- o$S
  }
  inf <- cr & !is.na(y) & y == Inf
  lP[inf] <- 0
  lS[inf] <- -Inf
  bad <- cr & is.na(y)
  lP[bad] <- lS[bad] <- NaN
  w <- which(cr & !is.na(y) & y != Inf & y > 0)
  if (!length(w)) return(list(P = lP, S = lS))
  acc <- c("A", "v", "sv")
  M <- cbind(y[w], pars[i1[w], acc, drop = FALSE], Tm[i1[w], -1, drop = FALSE],
             pars[i2[w], acc, drop = FALSE], Tm[i2[w], -1, drop = FALSE], rho[w])
  key <- do.call(paste, as.data.frame(M))
  first <- !duplicated(key)
  val <- vapply(which(first), function(r) {
    m <- M[r, ]
    mtlba_stop_rho1(m[1], m[2], m[3], m[4], m[4 + 1:K], m[5 + K], m[6 + K], m[7 + K],
                    m[7 + K + 1:K], m[8 + 2 * K], K, rule, pos = pos)
  }, numeric(2))
  at <- match(key, key[first])
  lP[w] <- val[1, at]
  lS[w] <- val[2, at]
  list(P = lP, S = lS)
}

# mtlba_logstop_any with confidence guessing: the mixture over which
# accumulators trigger. A triggering accumulator has stopped the race once its
# evidence is past gb, so within a state the stopping time is the plain model's
# with that accumulator's thresholds capped at gb (bands above gb are empty);
# guess = "trial": the guess state is the BE race to the two guess thresholds,
# under either rule. Both mixtures are sums of positive terms.
mtlba_logstop_guess <- function(y, pars, Tm, i1, i2, K, rule, G = NULL, guess = "accumulator",
                                pos = TRUE) {
  if (is.null(G)) return(mtlba_logstop_any(y, pars, Tm, i1, i2, K, rule, pos))
  g1 <- G$g[i1]; g2 <- G$g[i2]
  capT <- function(c1, c2) {
    TT <- Tm
    if (c1) TT[i1, -1] <- pmin(Tm[i1, -1, drop = FALSE], G$gb[i1])
    if (c2) TT[i2, -1] <- pmin(Tm[i2, -1, drop = FALSE], G$gb[i2])
    TT
  }
  states <- if (guess == "trial")
    list(list(w = 1 - g1, c = c(FALSE, FALSE), rule = rule), list(w = g1, c = c(TRUE, TRUE), rule = "BE"))
  else list(list(w = (1 - g1) * (1 - g2), c = c(FALSE, FALSE), rule = rule),
            list(w = g1 * (1 - g2), c = c(TRUE, FALSE), rule = rule),
            list(w = (1 - g1) * g2, c = c(FALSE, TRUE), rule = rule),
            list(w = g1 * g2, c = c(TRUE, TRUE), rule = rule))
  lP <- lS <- rep(-Inf, length(y))
  for (st in states) {
    w <- st$w
    use <- is.na(w) | w > 0
    if (!any(use)) next
    o <- mtlba_logstop_any(y[use], pars, capT(st$c[1], st$c[2]), i1[use], i2[use], K, st$rule, pos)
    lP[use] <- mtlba_lse(lP[use], log(w[use]) + o$P)
    lS[use] <- mtlba_lse(lS[use], log(w[use]) + o$S)
  }
  list(P = lP, S = lS)
}

# log P(a < dt <= b) from the stop/survivor pairs at a and b (a <= b): the route
# (survivors or stopped) with the smaller leading term.
mtlba_logwindow <- function(la, lb) {
  useS <- la$S < lb$P
  useS[is.na(useS)] <- FALSE
  out <- mtlba_ldiff(lb$P, la$P)
  if (any(useS)) out[useS] <- mtlba_ldiff(la$S[useS], lb$S[useS])
  out
}

# Censoring and truncation columns of a dadm, one value per trial (first row),
# with EMC2's defaults for absent columns.
mtlba_bounds <- function(dadm, first) {
  col <- function(x, def) if (x %in% names(dadm)) as.numeric(dadm[[x]][first]) else
    rep(def, length(first))
  miss <- if ("missingness" %in% names(dadm)) as.integer(dadm$missingness[first]) else
    rep(NA_integer_, length(first))
  list(miss = miss, LT = col("LT", 0), UT = col("UT", Inf), LC = col("LC", 0), UC = col("UC", Inf))
}

# Refuses what the censoring/truncation terms do not (yet) define: missing
# responses that are not censored (e.g. pContaminant), response-retained
# censoring, and code 4 (silent accumulators).
mtlba_check_missing <- function(R, rt, miss) {
  cens <- !is.na(miss)
  if (any(miss[cens] == 4L))
    stop("MTLBA: missingness code 4 (silent accumulators) is not supported")
  if (any(cens & !is.na(R)))
    stop("MTLBA: censoring with the response retained (LCresponse/UCresponse = TRUE) is not ",
         "supported; censored trials need R = NA (and RR = NA)")
  if (any(!cens & (is.na(R) | !is.finite(rt))))
    stop("MTLBA: missing responses or non-finite rts must be censored (make_missing() with LC/UC); ",
         "contamination (pContaminant) is not supported")
  invisible(TRUE)
}

# log likelihood per trial with censoring and truncation: observed trials the
# kernel density, censored trials (response not retained) the probability of
# stopping in the censored window(s), (1) (LT, LC], (2) (UC, UT], (3) both,
# and every truncated trial divided by P(LT < rt < UT). The window terms need
# one t0 per trial, so t0 must be equal for the two accumulators there.
mtlba_logdens_missing <- function(rt, R, RR, pars, K, rule, dont_know, bnd,
                                  guess = "accumulator", gb_constraint = "none", pos = TRUE) {
  n <- length(rt)
  mtlba_check_missing(R, rt, bnd$miss)
  cens <- !is.na(bnd$miss)
  trunc <- (bnd$LT > 0 & is.finite(bnd$LT)) | is.finite(bnd$UT)
  ld <- rep(-Inf, n)
  if (any(!cens))
    ld[!cens] <- mtlba_logdens(rt[!cens], R[!cens], RR[!cens],
                               pars[rep(2 * which(!cens), each = 2) - c(1, 0), , drop = FALSE],
                               K, rule, dont_know, guess = guess, gb_constraint = gb_constraint,
                               pos = pos)
  need <- cens | trunc
  if (!any(need)) return(ld)
  w <- which(need)
  i1 <- 2 * w - 1
  i2 <- 2 * w
  t0 <- pars[i1, "t0"]
  if (any(t0 != pars[i2, "t0"], na.rm = TRUE))
    stop("MTLBA: censoring and truncation need t0 equal for both accumulators (t0 must not ",
         "depend on lR or lM)")
  Tm <- mtlba_thresholds(pars, K)
  G <- mtlba_guess_pars(pars, K, guess, gb_constraint)
  at <- function(x) mtlba_logstop_guess(x - t0, pars, Tm, i1, i2, K, rule, G, guess, pos)
  sLT <- at(bnd$LT[w])
  sUT <- at(bnd$UT[w])
  cw <- cens[w]
  m <- bnd$miss[w]
  if (any(cw)) {
    lo <- m %in% c(1L, 3L)
    up <- m %in% c(2L, 3L)
    lc <- rep(-Inf, length(w))
    if (any(lo)) {
      sLC <- at(ifelse(lo, bnd$LC[w], 0))
      lc[lo] <- mtlba_logwindow(sLT, sLC)[lo]
    }
    if (any(up)) {
      sUC <- at(ifelse(up, bnd$UC[w], Inf))
      lc[up] <- mtlba_lse(lc[up], mtlba_logwindow(sUC, sUT)[up])
    }
    ld[w[cw]] <- lc[cw]
  }
  tw <- trunc[w]
  if (any(tw)) {
    lZ <- mtlba_logwindow(sLT, sUT)[tw]
    ld[w[tw]] <- ifelse(lZ == -Inf, -Inf, ld[w[tw]] - lZ)
  }
  ld[is.na(ld)] <- -Inf
  ld
}

# Trial log-likelihoods over a (possibly compressed) dadm with two rows per
# trial, expanded to the uncompressed trials and floored at min_ll. The R
# reference for the C++ kernel (src/model_mtlba.h); used when c_name is removed.
trial_ll_mtlba <- function(pars, dadm, K, rule, dont_know, min_ll = log(1e-10),
                           guess = "accumulator", gb_constraint = "none", drifts = "any") {
  n <- nrow(pars) / 2
  first <- seq(1, 2 * n, by = 2)
  bnd <- mtlba_bounds(dadm, first)
  ll <- mtlba_logdens_missing(dadm$rt[first], as.integer(dadm$R[first]), dadm$RR[first], pars,
                              K, rule, dont_know, bnd, guess, gb_constraint, mtlba_pos(drifts))
  # a parameter out of bounds on either accumulator row voids the trial
  ok <- attr(pars, "ok")
  if (!is.null(ok)) ll[!(ok[first] & ok[first + 1])] <- min_ll
  ll[is.na(ll)] <- min_ll
  pmax(min_ll, ll[attr(dadm, "expand")])
}

log_likelihood_mtlba <- function(pars, dadm, K, rule, dont_know, min_ll = log(1e-10),
                                 guess = "accumulator", gb_constraint = "none", drifts = "any")
  sum(trial_ll_mtlba(pars, dadm, K, rule, dont_know, min_ll, guess, gb_constraint, drifts))

# Data checks: two response levels, the rating column, and censoring coded as
# EMC2's make_missing() does with the response not retained.
mtlba_check_data <- function(data, K, dont_know) {
  if (length(levels(data$R)) != 2) stop("MTLBA needs exactly two response levels")
  miss <- if ("missingness" %in% names(data)) as.integer(data$missingness) else
    rep(NA_integer_, nrow(data))
  # RR on a censored trial is ignored (make_data() writes NA there)
  mtlba_check_missing(data$R, data$rt, miss)
  rating_check_data(data, K, dont_know)
}

# prepare_design: the rating checks; rho is a trial-level parameter (it cannot
# depend on lR or lM).
mtlba_prepare_design <- function(trial_level = "rho")
  function(formula, constants, Rlevels = NULL, ...)
    rating_prepare_design(n_choices = 2, trial_level = trial_level)(formula, constants, Rlevels, ...)

#' The Multiple-Threshold Linear Ballistic Accumulator (Rating) Model
#'
#' Model file for the multiple-threshold LBA (MTLBA), a model of choice,
#' response time and a rating (e.g. confidence). Two LBA accumulators race;
#' each has, besides its choice threshold `b = A + B`, `n_ratings - 1`
#' confidence thresholds below `b`. The two drift rates are (bivariate) normal
#' and may be correlated (`rho`). By default the model is the one of its
#' source paper (the DMC model `MTRp3`): the drifts are conditioned only on
#' not both being negative (`drifts = "any"`), and under the balance-of-evidence
#' rule the lowest confidence threshold may lie below `A`
#' (`lowest_threshold = "free"`). The alternative `drifts = "positive"`
#' (both drifts truncated to be positive) with `lowest_threshold = "above_A"`
#' is the model in which, with `rho = 0` and one rating, each accumulator is
#' the standard [LBA()].
#'
#' The data need a factor `R` with two levels, `rt`, and a numeric column `RR`
#' holding the rating, an integer from 1 (lowest) to `n_ratings` (highest).
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Interpretation** |
#' |---------------|---------------|-------------------|-------------|--------------------|
#' | *v*   | -     | \[-Inf, Inf\] | 1        | Mean drift rate |
#' | *sv*  | log   | \[0, Inf\]    | log(1)   | Between-trial SD of drift rate |
#' | *B*   | log   | \[0, Inf\]    | log(1)   | Choice threshold minus the maximum start point |
#' | *A*   | log   | \[0, Inf\]    | log(0)   | Maximum start point |
#' | *t0*  | log   | \[0, Inf\]    | log(0)   | Non-decision time |
#' | *rho* | probit on (-1, 1) | \[-1, 1\] | 0 | Correlation between the two drift rates |
#' | *c1*  | log   | \[0, Inf\]    | log(0.5) | First rating criterion, -log of the largest threshold proportion |
#' | *ck*  | log   | \[0, Inf\]    | log(0.5) | Increment to the next criterion |
#' | *gp*  | probit | \[0, 1\]     | qnorm(0)   | Probability that an accumulator triggers confidence guessing |
#' | *gB*  | probit | \[0, 1\]     | qnorm(0.5) | Guess threshold as a proportion of `B` above `A` |
#' | *pg1* | probit | \[0, 1\]     | qnorm(1)   | Probability that a guess gives rating 1 (the lowest) |
#' | *pgk* | probit | \[0, 1\]     | qnorm(0)   | Probability that a guess gives rating `k`, given it is not a lower rating |
#' | *pq*  | probit | \[0, 1\]     | qnorm(0)   | With `guess_ratings = "binomial"`, in place of the `pg`: the binomial probability `q` of guessed ratings |
#'
#' The guess-rating probabilities (`pg1 ..`, or `pq`) are bounded to `(0, 1)`
#' but also allowed at their defaults (`pg1 = 1`, the others 0), so a
#' sampled value whose probability rounds to the other end (beyond about
#' 8.3 on the probit scale) is out of bounds.
#'
#' Confidence thresholds are `T_j = A + B * d_j` with proportions
#' `d_{K-k} = exp(-(c1 + ... + ck))`, so `A < T_1 < ... < T_{K-1} < b`. With
#' `lowest_threshold = "free"` the lowest threshold is instead
#' `T_1 = T_2 * exp(-c_{K-1})`: the last criterion scales the absolute
#' threshold rather than its proportion above `A` (the two agree when
#' `A = 0`), so `T_1` lies anywhere in `(0, T_2)`, below `A` when the
#' reported proportion `d1` is negative. The natural-scale proportions
#' `d1 .. d{K-1}` are reported with `add_recalculated = TRUE`.
#'
#' Drift convention (`drifts`):
#' * `"any"` (the default, the source paper's model): the drift pair is
#'   bivariate normal conditioned on not both drifts being negative (on such
#'   a trial no accumulator would ever finish), so the likelihood is divided
#'   by `1 - P(both < 0)`. A losing accumulator may have a negative drift and
#'   end below its start point, which is what makes a lowest threshold below
#'   `A` informative: the two options are meant to be used together, and on
#'   the paper's data they fit far better together than either alone.
#' * `"positive"`: both drifts are truncated to be positive, the likelihood
#'   divided by `P(both > 0)`; with `rho = 0` this is the [LBA()] convention
#'   for each accumulator. With `lowest_threshold = "free"` the band below
#'   `A` uses the positive-drift probabilities.
#'
#' Parameterization (`parameterization`): the table above is `"additive"`
#' (the default), with `b = A + B` and `A` absolute. `"proportional"` is the
#' parameterization of the model's source paper and its DMC code, in which
#' the start-point range and every confidence threshold are proportions, so
#' they scale with `b`. The parameter names are the same; their meaning,
#' transforms and defaults change:
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** | **Interpretation** |
#' |---------------|---------------|-------------------|-------------|--------------------|
#' | *B*   | log    | \[0, Inf\] | log(1)     | The choice threshold `b` itself |
#' | *c1*  | probit | \[0, 1\]   | qnorm(0.5) | Highest confidence threshold `T_{K-1}` as a proportion of `b` |
#' | *A*   | probit | \[0, 1\]   | qnorm(0.5) | Maximum start point as a proportion of `T_{K-1}` (of `b` when `n_ratings = 1`) |
#' | *ck*  | probit | \[0, 1\]   | qnorm(0.5) | Threshold `T_{K-k}` as a proportion of the one above it |
#'
#' Below the highest one, each confidence threshold lies between `A` and the
#' threshold above it, `T_{K-k} = A + (T_{K-k+1} - A) * ck`, except a free
#' lowest threshold, which is a proportion of `T_2` alone,
#' `T_1 = T_2 * c_{K-1}`, and so may lie below `A` (the source paper's model;
#' with `n_ratings = 2` the one confidence threshold is the highest, above
#' `A`, and `lowest_threshold` makes no difference). Every other parameter,
#' `gB` included, keeps its meaning. A formula such as `B ~ SA * lR` with
#' `A ~ 1, c1 ~ lR, c2 ~ lR` therefore keeps the start-point range and the
#' confidence thresholds at the same proportions of `b` across the levels of
#' `SA`, and `A ~ 1` gives a different start-point range to each accumulator
#' whenever `b` or `c1` differ between them. The likelihood is that of the
#' additive model at the corresponding values (the translation is exact).
#' Mapped parameters ([mapped_pars()], `map = TRUE` in [get_pars()] and the
#' plots) report `A` as the natural start-point range, over the design cells
#' it varies in (e.g. `A_lRnew` and `A_lRold`), `B` as `b`, and with
#' `add_recalculated = TRUE` also `b` and the proportions `d1 .. d{K-1}` of
#' `b - A` above `A` (as under `"additive"`). The source paper's sampled
#' parameters map one to one: `B = log(b)`, `c1 = zd2`, `c2 = zd1`,
#' `A = zA`, `gB = -zfb`, `gp = zf` and `rho = zcorr` (with equal `sd_v` for
#' the two accumulators).
#'
#' `lowest_threshold`: `"above_A"` keeps every confidence threshold above `A`,
#' `"free"` lets the lowest one lie below `A` as described above, and
#' `"auto"` (the default) resolves to `"free"` under `rule = "BE"` and to
#' `"above_A"` under `rule = "TC"`: under threshold counting a threshold below
#' the start point would be crossed at time zero, so `"free"` is refused
#' there. `"free"` makes no difference with `n_ratings = 1`.
#'
#' Two stopping rules:
#' * `rule = "BE"` (balance of evidence): the first accumulator to reach `b`
#'   gives the choice, and the rating is `c` when the other accumulator's
#'   evidence lies between its thresholds `T_{K-c}` and `T_{K-c+1}` (`T_0 = 0`),
#'   so rating `K` means the other accumulator is below its lowest threshold.
#'   The thresholds read are the losing accumulator's, so the criteria are
#'   usually mapped `~ lR`.
#' * `rule = "TC"` (threshold count): the race stops at the `K`-th threshold
#'   crossing, counted over both accumulators. The choice is the accumulator
#'   that has crossed more thresholds (on a tie, the one that did not make the
#'   last crossing), and the rating increases with its lead.
#'
#' With `dont_know = TRUE` the lowest rating of both choices is one shared
#' "don't know" response (`2 * n_ratings - 1` responses in all). Code such
#' trials `RR = 1` with `R` equal to the first response level; `simulate`
#' and `make_data` write them that way.
#'
#' `t0` is taken from the accumulator that stops the race.
#'
#' `rho` is a constant 0 unless it is given a formula. It is a trial-level
#' parameter: it may depend on design factors but not on `lR` or `lM`. With
#' `rho != 0` the likelihood is a one-dimensional integral over the stopping
#' accumulator's start point, done by Gauss-Legendre quadrature (16 nodes per
#' panel, the number of panels set by the width of the integrand), so a
#' correlated model is slower than an uncorrelated one, by roughly the number
#' of nodes per trial.
#'
#' Confidence guessing (off unless `gp` is given a formula or a non-zero
#' constant): each accumulator has a guess threshold `A + B * gB`, between `A`
#' and its choice threshold. An accumulator that is *triggering* on a trial
#' ends the race when it reaches its guess threshold: it gives the choice,
#' usually faster than an evidence-based response, and the rating is drawn at
#' random, rating 1 with probability `pg1`, rating 2 with probability
#' `(1 - pg1) * pg2`, and so on, the highest rating taking what is left (there
#' are `n_ratings - 1` such parameters, read from the chosen accumulator; with
#' `dont_know` a guessed rating 1 is the don't-know response). If no
#' triggering accumulator reaches its guess threshold first, the race ends by
#' the stopping rule as usual.
#'
#' With `guess_ratings = "binomial"` the guessed rating is instead
#' `1 + Binomial(n_ratings - 1, q)`, one parameter `pq` (`q = pnorm(pq)`)
#' replacing `pg1 .. pg{K-1}`: rating `k` has probability
#' `choose(K - 1, k - 1) * q^(k - 1) * (1 - q)^(K - k)`, so with three ratings
#' `(1 - q)^2`, `2 q (1 - q)` and `q^2`. `q` near 0 describes a participant
#' whose guesses are mostly low ratings, `q` near 1 one whose guesses are
#' mostly high, and `q = .5` guesses centred on the middle; what a single `q`
#' cannot describe is a preference for the middle beyond the binomial's (the
#' middle of three ratings has probability at most .5) or for both extremes
#' over the middle. The likelihood is exactly that of the free (`"free"`, the
#' default) model at the `pg` the binomial implies, `pg_k = P(k) / P(>= k)`,
#' and these are reported with `add_recalculated = TRUE` (and by
#' [mapped_pars()]). The default `pq = qnorm(0)` guesses the lowest rating,
#' as the free model's defaults do.
#' * `guess = "accumulator"` (the default): the two accumulators trigger
#'   independently, each with its own probability `gp`, so `gp ~ lR` (or
#'   `~ lM`) gives different amounts of guessing for the two choices. Under
#'   `rule = "BE"` a triggering accumulator never reaches its choice threshold;
#'   under `rule = "TC"` its crossings of thresholds below the guess threshold
#'   still count.
#' * `guess = "trial"`: the whole trial is a guess with probability `gp`, a
#'   trial-level parameter (not `~ lR` or `~ lM`). Both accumulators then race
#'   to their guess thresholds, under either rule.
#'
#' `gb_constraint` is the option the model's source paper used for one of its
#' data sets: the guess threshold moves towards the choice threshold as
#' triggering becomes frequent, `A + B * (gB + g^2 * (1 - gB))`, with `g` the
#' accumulator's own `gp` (`"own"`) or the larger of the two accumulators'
#' `gp` (`"max"`, the paper's rule). In the limit where an accumulator always
#' triggers its choice threshold is never reached and `B` is then set by the
#' guess threshold. The price is that frequent guessing at a low threshold
#' cannot be represented (with `g = 0.5` the guess threshold is at least a
#' quarter of the way from `A` to `b`, with `g = 0.9` at least 81%). The
#' resulting proportion is reported as `gBc` with `add_recalculated = TRUE`.
#' The default is no constraint. Whatever the setting, ratings are informative
#' about confidence thresholds only on trials that end by the stopping rule: if
#' one accumulator nearly always triggers, the other accumulator's confidence
#' thresholds (read when the first one wins) are barely identified.
#'
#' The three kinds of guessing parameter inform each other: with `gp` near 0
#' the data say nothing about `gB` and the `pg`, and with `gp` near 1 nothing
#' about the confidence thresholds, so they need informative priors. A model
#' with guessing costs up to four times the likelihood evaluations of the same
#' model without it.
#'
#' Censoring and truncation use EMC2's usual coding ([make_missing()], or the
#' `TC` argument of [make_data()]). A censored trial has `R = NA` and
#' `RR = NA` (the response is not retained) and contributes the probability
#' that the race stops in the censored window; truncated trials are divided by
#' the probability of a response time between `LT` and `UT`. Both need `t0`
#' to be the same for the two accumulators. Censoring with the response
#' retained (`LCresponse`/`UCresponse = TRUE`) and missing responses that
#' are not censored (contamination) are not supported.
#'
#' @param n_ratings Number of rating categories K (>= 1; >= 2 with
#'   `dont_know`). With K = 1 the model is the LBA.
#' @param dont_know Logical: rating 1 is one shared "don't know" response.
#' @param rule Stopping rule, `"BE"` (balance of evidence) or `"TC"` (threshold
#'   count).
#' @param guess How confidence guessing is triggered when `gp` is not zero:
#'   `"accumulator"` (independently per accumulator) or `"trial"` (per trial).
#' @param gb_constraint `"none"` (the guess threshold is `A + B * gB`), or
#'   `"own"` / `"max"`: the guess threshold moves towards the choice threshold
#'   as the trigger probability (the accumulator's own, or the larger of the
#'   two) approaches 1. See Details.
#' @param drifts Drift convention, `"any"` (the drift pair is conditioned on
#'   not both drifts being negative; the default) or `"positive"` (both
#'   drifts truncated to be positive). See Details.
#' @param lowest_threshold `"auto"` (the default: `"free"` under `rule = "BE"`,
#'   `"above_A"` under `rule = "TC"`), `"free"` (the lowest confidence
#'   threshold may lie below `A`; balance of evidence only) or `"above_A"`.
#'   See Details.
#' @param parameterization `"additive"` (the default: `b = A + B`, the
#'   confidence thresholds proportions of `B` above `A`) or `"proportional"`
#'   (the source paper's: `B` is `b`, and `A` and the confidence thresholds are
#'   nested proportions of `b`). See Details.
#' @param guess_ratings How a guessed rating is drawn: `"free"` (the default:
#'   `n_ratings - 1` probabilities `pg1 .. pg{K-1}`) or `"binomial"`
#'   (`1 + Binomial(n_ratings - 1, q)`, one parameter `pq`). See Details.
#' @return A model list with all the necessary functions for EMC2 to sample.
#' @examples
#' # three ratings, balance of evidence (the source paper's model)
#' MTLBA3 <- function() MTLBA(n_ratings = 3)
#' # the same with both drifts positive and every threshold above A
#' MTLBA3p <- function() MTLBA(n_ratings = 3, drifts = "positive",
#'                             lowest_threshold = "above_A")
#' # the source paper's proportional parameterization
#' MTLBA3prop <- function() MTLBA(n_ratings = 3, parameterization = "proportional")
#' # confidence guessing with binomial guessed ratings (one parameter pq)
#' MTLBA3b <- function() MTLBA(n_ratings = 3, guess_ratings = "binomial")
#' # two ratings plus don't know, threshold count
#' MTLBA2dk <- function() MTLBA(n_ratings = 2, dont_know = TRUE, rule = "TC")
#' # confidence guessing, more often on one accumulator than on the other
#' design(factors = list(subjects = 1, S = c("left", "right")),
#'        Rlevels = c("left", "right"), matchfun = function(d) d$S == d$lR,
#'        model = MTLBA3,
#'        formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, c1 ~ 1, c2 ~ 1,
#'                       gp ~ lR, gB ~ 1, pg1 ~ 1, pg2 ~ 1))
#' @export
MTLBA <- function(n_ratings = 3, dont_know = FALSE, rule = "BE", guess = "accumulator",
                  gb_constraint = "none", drifts = "any", lowest_threshold = "auto",
                  parameterization = "additive", guess_ratings = "free") {
  K <- rating_check_n(n_ratings)
  rating_check_dont_know(K, dont_know)
  rule <- match.arg(rule, c("BE", "TC"))
  guess <- match.arg(guess, c("accumulator", "trial"))
  gb_constraint <- match.arg(gb_constraint, c("none", "own", "max"))
  drifts <- match.arg(drifts, c("any", "positive"))
  lowest_threshold <- match.arg(lowest_threshold, c("auto", "free", "above_A"))
  if (lowest_threshold == "auto") lowest_threshold <- if (rule == "BE") "free" else "above_A"
  if (lowest_threshold == "free" && rule == "TC")
    stop("MTLBA: lowest_threshold = \"free\" needs rule = \"BE\": under threshold counting a ",
         "threshold below the start point would be crossed at time zero (use \"above_A\" or \"auto\")")
  free <- lowest_threshold == "free" && K > 1
  parameterization <- match.arg(parameterization, c("additive", "proportional"))
  prop <- parameterization == "proportional"
  guess_ratings <- match.arg(guess_ratings, c("free", "binomial"))
  binom <- guess_ratings == "binomial" && K > 1
  crit <- rating_criterion_names(K)
  pg <- if (K > 1) paste0("pg", seq_len(K - 1)) else character(0)
  # binomial: one p_type pq (q = pnorm(pq)) in place of pg1 .. pg{K-1}, which Ttransform adds
  gr <- if (binom) "pq" else pg
  model <- list(
    type = "RACE",
    c_name = "MTLBA",
    # read by the C++ likelihood through attr(dadm, "c_options")
    c_options = c(n_ratings = K, dont_know = as.numeric(dont_know), TC = as.numeric(rule == "TC"),
                  guess_trial = as.numeric(guess == "trial"),
                  gb_constraint = match(gb_constraint, c("none", "own", "max")) - 1,
                  drifts_any = as.numeric(drifts == "any"), free_lowest = as.numeric(free),
                  proportional = as.numeric(prop), binomial_guess = as.numeric(binom)),
    n_ratings = K,
    dont_know = dont_know,
    rule = rule,
    guess = guess,
    gb_constraint = gb_constraint,
    drifts = drifts,
    lowest_threshold = lowest_threshold,
    parameterization = parameterization,
    guess_ratings = guess_ratings,
    extra_responses = "RR",
    # guessing is off by default: gp = 0, with gB = .5 and the lowest rating guessed
    # (binomial: q = 0)
    # proportional: A and the criteria are proportions (default .5), B = log b
    p_types = c("v" = 1, "sv" = log(1), "B" = log(1), "A" = if (prop) stats::qnorm(0.5) else log(0),
                "t0" = log(0), "rho" = 0,
                if (prop) stats::setNames(rep(stats::qnorm(0.5), length(crit)), crit) else rating_p_types(K),
                "gp" = stats::qnorm(0), "gB" = stats::qnorm(0.5),
                stats::setNames(stats::qnorm(as.numeric(seq_along(gr) == 1 & !binom)), gr)),
    transform = list(func = c(v = "identity", sv = "exp", B = "exp", A = if (prop) "pnorm" else "exp",
                              t0 = "exp", rho = "pnorm",
                              if (prop) stats::setNames(rep("pnorm", length(crit)), crit) else rating_transform(K),
                              gp = "pnorm", gB = "pnorm", stats::setNames(rep("pnorm", length(gr)), gr)),
                     lower = c(rho = -1), upper = c(rho = 1)),
    # (bounds are also checked after Ttransform when simulating, where a proportional
    # model's A column is the natural start-point range: hence (0, Inf) there)
    # the guess-rating probabilities (pg, or binomial pq) in (0, 1), each also allowed at its
    # default (pg1 = 1, the others 0: a guess gives the lowest rating)
    bound = list(minmax = cbind(v = c(-Inf, Inf), sv = c(0, Inf), A = c(if (prop) 0 else 1e-4, Inf),
                                B = c(0, Inf), t0 = c(0.05, Inf), rho = c(-1, 1),
                                if (prop) matrix(rep(c(0, 1), length(crit)), 2, dimnames = list(NULL, crit)) else
                                  rating_bound(K),
                                matrix(rep(c(0, 1), length(gr)), 2, dimnames = list(NULL, gr))),
                 exception = c(A = 0, stats::setNames(as.numeric(seq_along(gr) == 1 & !binom), gr))),
    Ttransform = function(pars, dadm) {
      ok <- attr(pars, "ok")
      if (prop) pars <- mtlba_proportional(pars, K, free) else {
        pars <- cbind(pars[, colnames(pars) != "b", drop = FALSE], b = pars[, "B"] + pars[, "A"])
        pars <- rating_add_thresholds(pars, K)
        # a free lowest threshold: T_1 = T_2 exp(-c_{K-1}), as a proportion of B above A
        if (free) {
          A <- pars[, "A"]; b <- pars[, "b"]
          T2 <- if (K > 2) A + (b - A) * pars[, "d2"] else b
          pars[, "d1"] <- (T2 * exp(-pars[, paste0("c", K - 1)]) - A) / (b - A)
        }
      }
      # the guess threshold as a proportion of B above A, after gb_constraint
      # (rows come two per trial; anything else, e.g. a one-row summary, takes its own gp)
      g <- pars[, "gp"]
      paired <- nrow(pars) %% 2 == 0
      if (guess == "trial" && paired) g <- rep(g[seq(1, length(g), by = 2)], each = 2)
      pars <- cbind(pars[, colnames(pars) != "gBc", drop = FALSE],
                    gBc = mtlba_gB_constrained(pars[, "gB"], g,
                                               if (paired) gb_constraint else sub("max", "own", gb_constraint)))
      # binomial guess ratings: the sequential pg1 .. pg{K-1} implied by q
      if (binom) pars <- cbind(pars[, !(colnames(pars) %in% pg), drop = FALSE],
                               mtlba_binomial_pg(pars[, "pq"], K))
      if (!is.null(ok)) attr(pars, "ok") <- ok
      pars
    },
    prepare_design = mtlba_prepare_design(c("rho", if (guess == "trial") "gp")),
    check_data = function(data) mtlba_check_data(data, K, dont_know),
    rfun = function(data = NULL, pars) rMTLBA(data$lR, pars, K, rule, dont_know,
                                              ok = attr(pars, "ok"), guess = guess,
                                              gb_constraint = gb_constraint, drifts = drifts),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10))
      log_likelihood_mtlba(pars, dadm, K, rule, dont_know, min_ll = min_ll, guess = guess,
                           gb_constraint = gb_constraint, drifts = drifts)
  )
  # the natural start-point range replaces the proportion A: map = TRUE reports it
  # over the cells it varies in (it scales with b and c1), not A's formula's
  if (prop) model$recalculated_p_types <- "A"
  model
}
