# Multiple-threshold LBA (MTLBA): two LBA accumulators (start point U(0, A),
# drift rate normal with mean v and sd sv, truncated to positive), each with
# K - 1 confidence thresholds between A and its choice threshold b = A + B,
#   T_j = A + B * d_j,  j = 1 .. K - 1,   0 < d_1 < ... < d_{K-1} < 1,
# (d_j from the rating layer's criteria, R/model_rating.R), T_K = b. A trial
# gives a choice R, its rt and a rating RR (1 = lowest).
#
# Every response is one KERNEL: the density that the STOPPING accumulator s
# reaches its threshold T_{s,j} at decision time dt = rt - t0 while the OTHER
# accumulator's evidence e_o = start + drift * dt lies in its band
# (T_{o,lo}, T_{o,lo+1}), lo = 0 meaning "below T_{o,1}":
#   kernel = f_LBA(dt; A_s, T_{s,j}, v_s, sv_s) * P(band | v_o > 0),
# with f_LBA the positive-drift LBA density (as EMC2's LBA()).
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
# G(z) = z Phi(z) + phi(z)), thresholds above A make the negative-drift mass
# cancel except in the lowest band, where H(0) := Phi(-v/sv), S(0) := Phi(v/sv)
# (plan M0 results), so
#   band * Phi(v/sv) = S(lo) - S(hi) = H(hi) - H(lo),
# and whichever route has the smaller leading term is used.
#
# Correlated drift rates (rho != 0): the two drifts are bivariate normal with
# correlation rho, truncated to both positive. Given the stopping accumulator's
# drift x (its start point is then bs - x * dt) the other drift is normal with
# mean v_o + rho * (sv_o / sv_s) * (x - v_s) and sd sv_o * sqrt(1 - rho^2), so
#   kernel = 1 / (A_s * P(both > 0)) *
#            integral over x in ((bs - A_s) / dt, bs / dt) of
#              x * dnorm(x; v_s, sv_s) * P(band, v_o > 0 | x) dx,
# a one-dimensional integral of an analytic integrand whose last factor is the
# band function above at the conditional mean and sd, left unnormalised. It is
# done by Gauss-Legendre quadrature on a window clipped to where the integrand
# lives (mtlba_quad_rule, mtlba_logkernel_rho). rho is a trial-level parameter
# (the same on both accumulator rows); trials with rho = 0 take the closed form.
#
# The compiled likelihood (src/model_mtlba.h, c_name "MTLBA") ports the
# uncorrelated algorithm and stops on rho != 0; the R functions here are its
# reference and the only correlated likelihood so far, used when c_name is
# removed.

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
# z1 = (y - t v) / (t sv), z0 = z1 - A / (t sv); S likewise with pnorm(-u).
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
  }
  list(H = lH, S = lS)
}

# log P(T_lo < e_o < T_hi | v_o > 0) at decision time t; lowest = lo is the
# sentinel "no threshold below" (lo_y is then ignored). The denominator
# pnorm(v/sv) is floored at 1e-10 as in EMC2's LBA. normalise = FALSE leaves
# the denominator out: log P(T_lo < e_o < T_hi, v_o > 0).
mtlba_logband <- function(t, A, v, sv, lo_y, hi_y, lowest, normalise = TRUE) {
  zv <- v / sv
  hi <- mtlba_logHS(hi_y, t, A, v, sv)
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

# log positive-drift LBA density of reaching threshold bs at t (EMC2's dlba
# formula, all terms taken on the scale of the largest).
mtlba_logf <- function(t, A, bs, v, sv) {
  out <- rep(-Inf, length(t))
  lden <- pmax(stats::pnorm(v / sv, log.p = TRUE), log(1e-10))
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
mtlba_quad_window <- function(dt, As, bs, vs, svs, Ao, vo, svo, lo_y, hi_y, lowest, rho, quad) {
  lo_x <- (bs - As) / dt
  hi_x <- bs / dt
  yl <- ifelse(lowest, 0, pmax(0, (lo_y - Ao) / dt))
  m <- vs + rho * svs / svo * (pmin(pmax(vo, yl), hi_y / dt) - vo)
  c_lo <- pmin(vs, m)
  c_hi <- pmax(vs, m)
  d <- pmax(lo_x - c_hi, c_lo - hi_x, 0)
  r <- sqrt(d^2 + (quad$hw * svs)^2)
  list(a = pmax(lo_x, c_lo - r), b = pmin(hi_x, c_hi + r))
}

# log of the correlated kernel times P(both drifts > 0): per trial the stopping
# accumulator (A_s, threshold bs, v_s, sv_s), the other accumulator (A_o, v_o,
# sv_o) and its band (lo_y, hi_y, lowest), decision time dt and rho. The
# integrand is formed on the log scale (band at the conditional mean and sd,
# unnormalised) and summed by log-sum-exp. A_s = 0: no integral, x = bs / dt.
mtlba_logkernel_rho <- function(dt, As, bs, vs, svs, Ao, vo, svo, lo_y, hi_y, lowest, rho,
                                quad = mtlba_quad_rule()) {
  n <- length(dt)
  out <- rep(-Inf, n)
  sc <- svo * sqrt(1 - rho^2)
  # log integrand at X (rows = the trials i, columns = nodes)
  lint <- function(i, X) {
    ii <- rep(i, ncol(X))
    x <- as.vector(X)
    mc <- vo[ii] + rho[ii] * svo[ii] / svs[ii] * (x - vs[ii])
    l <- log(x) + stats::dnorm(x, vs[ii], svs[ii], log = TRUE) +
      mtlba_logband(dt[ii], Ao[ii], mc, sc[ii], lo_y[ii], hi_y[ii], lowest[ii], normalise = FALSE)
    l[is.na(l)] <- -Inf
    matrix(l, length(i), ncol(X))
  }
  small <- As <= 1e-10
  if (any(small)) {
    i <- which(small)
    out[i] <- lint(i, matrix(bs[i] / dt[i], ncol = 1))[, 1] - log(dt[i])
  }
  if (any(!small)) {
    win <- mtlba_quad_window(dt, As, bs, vs, svs, Ao, vo, svo, lo_y, hi_y, lowest, rho, quad)
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

# log kernel per trial for given lookup indices. rt per trial; pars two rows
# per trial (accumulators in level order); Tm = mtlba_thresholds(pars, K).
mtlba_logkernel <- function(rt, pars, Tm, s, j, lo, quad = mtlba_quad_rule()) {
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
    fin(oi, c("v", "sv", "A", "b")) & is.finite(Tm[cbind(si, j + 1L)]) &
    is.finite(Tm[cbind(oi, lo + 2L)]) & (lo == 0L | is.finite(Tm[cbind(oi, lo + 1L)])) &
    is.finite(rho) & abs(rho) < 1
  ok[is.na(ok)] <- FALSE
  if (!any(ok)) return(out)
  si <- si[ok]; oi <- oi[ok]; dt <- dt[ok]; j <- j[ok]; lo <- lo[ok]; rho <- rho[ok]
  bs <- Tm[cbind(si, j + 1L)]
  lowest <- lo == 0L
  res <- numeric(length(dt))
  u <- rho == 0
  if (any(u))
    res[u] <- mtlba_logf(dt[u], pars[si[u], "A"], bs[u], pars[si[u], "v"], pars[si[u], "sv"]) +
      mtlba_logband(dt[u], pars[oi[u], "A"], pars[oi[u], "v"], pars[oi[u], "sv"],
                    Tm[cbind(oi[u], lo[u] + 1L)], Tm[cbind(oi[u], lo[u] + 2L)], lowest[u])
  if (any(!u)) {
    k <- which(!u)
    res[k] <- mtlba_logkernel_rho(dt[k], pars[si[k], "A"], bs[k], pars[si[k], "v"],
                                  pars[si[k], "sv"], pars[oi[k], "A"], pars[oi[k], "v"],
                                  pars[oi[k], "sv"], Tm[cbind(oi[k], lo[k] + 1L)],
                                  Tm[cbind(oi[k], lo[k] + 2L)], lowest[k], rho[k], quad) -
      pmax(mtlba_logppos(pars[si[k], "v"], pars[si[k], "sv"], pars[oi[k], "v"],
                         pars[oi[k], "sv"], rho[k]), log(1e-10))
  }
  out[ok] <- res
  out
}

# log(exp(a) + exp(b))
mtlba_lse <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log1p(exp(-abs(a - b)))
  out[m == -Inf] <- -Inf
  out
}

# log likelihood per trial: rt, R (1/2), RR per trial, pars two rows per trial.
mtlba_logdens <- function(rt, R, RR, pars, K, rule, dont_know, quad = mtlba_quad_rule()) {
  Tm <- mtlba_thresholds(pars, K)
  R <- as.integer(R)
  dk <- if (dont_know) !is.na(RR) & RR == 1 else rep(FALSE, length(R))
  # a don't-know trial is the sum of the two rating-1 terms: R = 1 here, R = 2 below
  R[dk] <- 1L
  lk <- mtlba_lookup(R, RR, K, rule)
  ld <- mtlba_logkernel(rt, pars, Tm, lk$s, lk$j, lk$lo, quad)
  if (any(dk)) {
    w <- which(dk)
    l2 <- mtlba_lookup(rep(2L, length(w)), rep(1, length(w)), K, rule)
    idx <- as.vector(rbind(2 * w - 1, 2 * w))
    ld2 <- mtlba_logkernel(rt[w], pars[idx, , drop = FALSE], Tm[idx, , drop = FALSE],
                           l2$s, l2$j, l2$lo, quad)
    ld[w] <- mtlba_lse(ld[w], ld2)
  }
  ld[is.na(ld)] <- -Inf
  ld
}

# Density of the MTLBA. rt, R (factor or integer 1/2) and RR (1 = lowest) per
# trial; pars has two rows per trial (accumulators in level order) with columns
# v, sv, A, b, t0 and natural-scale d1 .. d{K-1} (proportions of B = b - A
# above A), and optionally rho (drift-rate correlation, read from the trial's
# first row; 0 if absent). Don't-know trials: RR = 1 (R is ignored). quad: the
# quadrature rule used where rho != 0.
dMTLBA <- function(rt, R, RR, pars, n_ratings = 3, rule = "BE", dont_know = FALSE,
                   log = FALSE, quad = mtlba_quad_rule()) {
  K <- rating_check_n(n_ratings)
  rule <- match.arg(rule, c("BE", "TC"))
  rating_check_dont_know(K, dont_know)
  ld <- mtlba_logdens(rt, R, RR, pars, K, rule, dont_know, quad)
  if (log) ld else exp(ld)
}

# Simulate the MTLBA: lR the accumulator factor (two rows per trial), pars as
# for dMTLBA. Returns R, rt and RR; don't-know responses are written as RR = 1
# with R = the first level. Trials with rho != 0 draw their two drifts from the
# bivariate normal truncated to both positive (mtlba_rdrift_rho).
rMTLBA <- function(lR, pars, n_ratings = 3, rule = "BE", dont_know = FALSE,
                   ok = rep(TRUE, nrow(pars))) {
  K <- rating_check_n(n_ratings)
  rule <- match.arg(rule, c("BE", "TC"))
  rating_check_dont_know(K, dont_know)
  if (length(levels(lR)) != 2) stop("MTLBA simulates two accumulators only")
  if (is.null(ok)) ok <- rep(TRUE, nrow(pars))
  n <- nrow(pars) / 2
  first <- seq(1, 2 * n, by = 2)
  Tm <- mtlba_thresholds(pars, K)
  # crossing times of thresholds 1 .. K, one n x K matrix per accumulator
  cross <- lapply(0:1, function(a) {
    p <- pars[first + a, , drop = FALSE]
    start <- p[, "A"] * stats::runif(n)
    v <- msm::rtnorm(n, p[, "v"], p[, "sv"], lower = 0)
    list(tau = (Tm[first + a, -1, drop = FALSE] - start) / v, start = start, v = v)
  })
  rho <- if ("rho" %in% colnames(pars)) pars[first, "rho"] else rep(0, n)
  cr <- which(is.na(rho) | rho != 0)
  if (length(cr)) {
    vv <- mtlba_rdrift_rho(pars[first[cr], "v"], pars[first[cr], "sv"],
                           pars[first[cr] + 1, "v"], pars[first[cr] + 1, "sv"], rho[cr])
    for (a in 1:2) {
      cross[[a]]$v[cr] <- vv[[a]]
      cross[[a]]$tau[cr, ] <- (Tm[first[cr] + a - 1, -1, drop = FALSE] - cross[[a]]$start[cr]) /
        vv[[a]]
    }
  }
  t1 <- cross[[1]]$tau
  t2 <- cross[[2]]$tau
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
# above them). The response kernels summed over responses are the density of
# this stopping time. i1, i2: rows of accumulators 1 and 2; Tm from
# mtlba_thresholds.
mtlba_logstop <- function(y, pars, Tm, i1, i2, K, rule) {
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
    lden <- pmax(stats::pnorm(v / sv, log.p = TRUE), log(1e-10))
    m <- length(t)
    F <- Fc <- band <- matrix(-Inf, m, K)
    for (nn in 0:(K - 1)) {
      hi <- Tm[cbind(i, nn + 2L)]
      F[, nn + 1] <- mtlba_logband(t, A, v, sv, rep(0, m), hi, rep(TRUE, m))
      Fc[, nn + 1] <- mtlba_logHS(hi, t, A, v, sv)$S - lden
      band[, nn + 1] <- if (nn == 0) F[, 1] else
        mtlba_logband(t, A, v, sv, Tm[cbind(i, nn + 1L)], hi, rep(FALSE, m))
    }
    list(F = F, Fc = Fc, band = band)
  }
  a1 <- acc(i1)
  a2 <- acc(i2)
  if (rule == "BE") {
    S <- a1$F[, K] + a2$F[, K]
    P <- mtlba_lse(a1$Fc[, K], a1$F[, K] + a2$Fc[, K])
  } else {
    S <- a1$band[, 1] + a2$F[, K]
    P <- mtlba_lse(a1$Fc[, K], a1$band[, 1] + a2$Fc[, K])
    if (K > 1) for (n1 in 1:(K - 1)) {
      S <- mtlba_lse(S, a1$band[, n1 + 1] + a2$F[, K - n1])
      P <- mtlba_lse(P, a1$band[, n1 + 1] + a2$Fc[, K - n1])
    }
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
mtlba_stop_rho1 <- function(y, A1, v1, s1, T1, A2, v2, s2, T2, rho, K, rule,
                            quad = mtlba_quad_reference()) {
  if (!all(is.finite(c(y, A1, v1, s1, T1, A2, v2, s2, T2, rho))) || abs(rho) >= 1 ||
      s1 <= 0 || s2 <= 0) return(c(NaN, NaN))
  # window as mtlba_quad_window: v2 > 0 drags x towards m when v2's mean is negative
  m <- v1 + rho * s1 / s2 * (max(v2, 0) - v2)
  r <- sqrt(max(-max(v1, m), 0)^2 + (quad$hw * s1)^2)
  a <- max(0, min(v1, m) - r)
  b <- max(v1, m) + r
  kinks <- c((T1 - A1) / y, T1 / y)
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
  low2 <- sapply(1:K, function(m)
    exp(mtlba_logband(yy, AA, mc, sc, rep(0, nx), rep(T2[m], nx), rep(TRUE, nx), normalise = FALSE)))
  up2 <- sapply(1:K, function(m) exp(mtlba_logHS(rep(T2[m], nx), yy, AA, mc, sc)$S))
  dim(Fle) <- dim(low2) <- dim(up2) <- c(nx, K)
  pos2 <- stats::pnorm(mc / sc)
  if (rule == "BE") {
    Sint <- Fle[, K] * low2[, K]
    Pint <- (1 - Fle[, K]) * pos2 + Fle[, K] * up2[, K]
  } else {
    wn <- Fle - cbind(0, Fle[, -K, drop = FALSE])  # w_n1, n1 = 0 .. K - 1
    Sint <- rowSums(wn * low2[, K:1, drop = FALSE])
    Pint <- (1 - Fle[, K]) * pos2 + rowSums(wn * up2[, K:1, drop = FALSE])
  }
  log(c(sum(w * Pint), sum(w * Sint))) - mtlba_logppos(v1, s1, v2, s2, rho)
}

# mtlba_logstop for any rho: rows with rho = 0 take the closed form, the others
# mtlba_stop_rho1, evaluated once per distinct (bound, parameter row).
mtlba_logstop_any <- function(y, pars, Tm, i1, i2, K, rule) {
  rho <- if ("rho" %in% colnames(pars)) pars[i1, "rho"] else rep(0, length(i1))
  cr <- is.na(rho) | rho != 0
  if (!any(cr)) return(mtlba_logstop(y, pars, Tm, i1, i2, K, rule))
  n <- length(y)
  lP <- rep(-Inf, n)
  lS <- rep(0, n)
  if (any(!cr)) {
    o <- mtlba_logstop(y[!cr], pars, Tm, i1[!cr], i2[!cr], K, rule)
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
                    m[7 + K + 1:K], m[8 + 2 * K], K, rule)
  }, numeric(2))
  at <- match(key, key[first])
  lP[w] <- val[1, at]
  lS[w] <- val[2, at]
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
mtlba_logdens_missing <- function(rt, R, RR, pars, K, rule, dont_know, bnd) {
  n <- length(rt)
  mtlba_check_missing(R, rt, bnd$miss)
  cens <- !is.na(bnd$miss)
  trunc <- (bnd$LT > 0 & is.finite(bnd$LT)) | is.finite(bnd$UT)
  ld <- rep(-Inf, n)
  if (any(!cens))
    ld[!cens] <- mtlba_logdens(rt[!cens], R[!cens], RR[!cens],
                               pars[rep(2 * which(!cens), each = 2) - c(1, 0), , drop = FALSE],
                               K, rule, dont_know)
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
  at <- function(x) mtlba_logstop_any(x - t0, pars, Tm, i1, i2, K, rule)
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
trial_ll_mtlba <- function(pars, dadm, K, rule, dont_know, min_ll = log(1e-10)) {
  n <- nrow(pars) / 2
  first <- seq(1, 2 * n, by = 2)
  bnd <- mtlba_bounds(dadm, first)
  ll <- mtlba_logdens_missing(dadm$rt[first], as.integer(dadm$R[first]), dadm$RR[first], pars,
                              K, rule, dont_know, bnd)
  # a parameter out of bounds on either accumulator row voids the trial
  ok <- attr(pars, "ok")
  if (!is.null(ok)) ll[!(ok[first] & ok[first + 1])] <- min_ll
  ll[is.na(ll)] <- min_ll
  pmax(min_ll, ll[attr(dadm, "expand")])
}

log_likelihood_mtlba <- function(pars, dadm, K, rule, dont_know, min_ll = log(1e-10))
  sum(trial_ll_mtlba(pars, dadm, K, rule, dont_know, min_ll))

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
mtlba_prepare_design <- function(formula, constants, Rlevels = NULL, ...)
  rating_prepare_design(n_choices = 2, trial_level = "rho")(formula, constants, Rlevels, ...)

#' The Multiple-Threshold Linear Ballistic Accumulator (Rating) Model
#'
#' Model file for the multiple-threshold LBA (MTLBA), a model of choice,
#' response time and a rating (e.g. confidence). Two LBA accumulators race;
#' each has, besides its choice threshold `b = A + B`, `n_ratings - 1`
#' confidence thresholds between `A` and `b`. The two drift rates are normal,
#' truncated to be both positive, and may be correlated (`rho`); with
#' `rho = 0` (the default) and one rating the model is the standard [LBA()].
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
#'
#' Confidence thresholds are `T_j = A + B * d_j` with proportions
#' `d_{K-k} = exp(-(c1 + ... + ck))`, so `A < T_1 < ... < T_{K-1} < b`. The
#' natural-scale proportions `d1 .. d{K-1}` are reported with
#' `add_recalculated = TRUE`.
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
#' accumulator's start point, done by Gauss-Legendre quadrature. At present
#' only the R likelihood has it: the compiled likelihood stops when it meets a
#' non-zero `rho`, so a correlated model must be run with `c_name` removed
#' from the model list (slow; for checking, not for large fits).
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
#' @return A model list with all the necessary functions for EMC2 to sample.
#' @examples
#' # three ratings, balance of evidence
#' MTLBA3 <- function() MTLBA(n_ratings = 3)
#' # two ratings plus don't know, threshold count
#' MTLBA2dk <- function() MTLBA(n_ratings = 2, dont_know = TRUE, rule = "TC")
#' @export
MTLBA <- function(n_ratings = 3, dont_know = FALSE, rule = "BE") {
  K <- rating_check_n(n_ratings)
  rating_check_dont_know(K, dont_know)
  rule <- match.arg(rule, c("BE", "TC"))
  list(
    type = "RACE",
    c_name = "MTLBA",
    # read by the C++ likelihood through attr(dadm, "c_options")
    c_options = c(n_ratings = K, dont_know = as.numeric(dont_know), TC = as.numeric(rule == "TC")),
    n_ratings = K,
    dont_know = dont_know,
    rule = rule,
    extra_responses = "RR",
    p_types = c("v" = 1, "sv" = log(1), "B" = log(1), "A" = log(0), "t0" = log(0), "rho" = 0,
                rating_p_types(K)),
    transform = list(func = c(v = "identity", sv = "exp", B = "exp", A = "exp", t0 = "exp",
                              rho = "pnorm", rating_transform(K)),
                     lower = c(rho = -1), upper = c(rho = 1)),
    bound = list(minmax = cbind(v = c(-Inf, Inf), sv = c(0, Inf), A = c(1e-4, Inf),
                                B = c(0, Inf), t0 = c(0.05, Inf), rho = c(-1, 1),
                                rating_bound(K)),
                 exception = c(A = 0)),
    Ttransform = function(pars, dadm) {
      ok <- attr(pars, "ok")
      pars <- cbind(pars[, colnames(pars) != "b", drop = FALSE], b = pars[, "B"] + pars[, "A"])
      pars <- rating_add_thresholds(pars, K)
      if (!is.null(ok)) attr(pars, "ok") <- ok
      pars
    },
    prepare_design = mtlba_prepare_design,
    check_data = function(data) mtlba_check_data(data, K, dont_know),
    rfun = function(data = NULL, pars) rMTLBA(data$lR, pars, K, rule, dont_know,
                                              ok = attr(pars, "ok")),
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10))
      log_likelihood_mtlba(pars, dadm, K, rule, dont_know, min_ll = min_ll)
  )
}
