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
# This file is the R reference (no c_name); rho (correlated drift rates) is a
# p_type fixed at 0 until the correlated likelihood exists.

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
      for (k in 120:2) tt <- xt + k / tt
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
# pnorm(v/sv) is floored at 1e-10 as in EMC2's LBA.
mtlba_logband <- function(t, A, v, sv, lo_y, hi_y, lowest) {
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
  lb - pmax(stats::pnorm(zv, log.p = TRUE), log(1e-10))
}

# log positive-drift LBA density of reaching threshold bs at t (EMC2's dlba
# formula, the two differences taken on a common scale).
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
    M <- pmax(lp1, lp0)
    dP <- exp(log_pnorm_interval(czm, cz) - M)  # pnorm(cz) - pnorm(czm)
    inner <- v[w] * dP + sv[w] * (exp(lp0 - M) - exp(lp1 - M))
    lo <- rep(-Inf, length(cz))
    pos <- !is.na(inner) & inner > 0
    lo[pos] <- log(inner[pos]) + M[pos]
    out[w] <- lo - log(A[w]) - lden[w]
  }
  out
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
mtlba_logkernel <- function(rt, pars, Tm, s, j, lo) {
  n <- length(rt)
  first <- seq(1, 2 * n, by = 2)
  si <- first + s - 1L
  oi <- first + 2L - s
  dt <- rt - pars[si, "t0"]
  out <- rep(-Inf, n)
  # non-finite parameters (e.g. an overflowing proposal) give -Inf, never an error
  fin <- function(i, cols) rowSums(!is.finite(pars[i, cols, drop = FALSE])) == 0
  ok <- !is.na(dt) & is.finite(dt) & dt > 0 & fin(si, c("v", "sv", "A", "b")) &
    fin(oi, c("v", "sv", "A", "b")) & is.finite(Tm[cbind(si, j + 1L)]) &
    is.finite(Tm[cbind(oi, lo + 2L)]) & (lo == 0L | is.finite(Tm[cbind(oi, lo + 1L)]))
  ok[is.na(ok)] <- FALSE
  if (!any(ok)) return(out)
  si <- si[ok]; oi <- oi[ok]; dt <- dt[ok]; j <- j[ok]; lo <- lo[ok]
  bs <- Tm[cbind(si, j + 1L)]
  lowest <- lo == 0L
  out[ok] <- mtlba_logf(dt, pars[si, "A"], bs, pars[si, "v"], pars[si, "sv"]) +
    mtlba_logband(dt, pars[oi, "A"], pars[oi, "v"], pars[oi, "sv"],
                  Tm[cbind(oi, lo + 1L)], Tm[cbind(oi, lo + 2L)], lowest)
  out
}

# log(exp(a) + exp(b))
mtlba_lse <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log1p(exp(-abs(a - b)))
  out[m == -Inf] <- -Inf
  out
}

mtlba_check_rho <- function(pars) {
  if ("rho" %in% colnames(pars) && any(pars[, "rho"] != 0, na.rm = TRUE))
    stop("MTLBA: correlated drift rates (rho != 0) are not implemented yet")
}

# log likelihood per trial: rt, R (1/2), RR per trial, pars two rows per trial.
mtlba_logdens <- function(rt, R, RR, pars, K, rule, dont_know) {
  mtlba_check_rho(pars)
  Tm <- mtlba_thresholds(pars, K)
  R <- as.integer(R)
  dk <- if (dont_know) !is.na(RR) & RR == 1 else rep(FALSE, length(R))
  # a don't-know trial is the sum of the two rating-1 terms: R = 1 here, R = 2 below
  R[dk] <- 1L
  lk <- mtlba_lookup(R, RR, K, rule)
  ld <- mtlba_logkernel(rt, pars, Tm, lk$s, lk$j, lk$lo)
  if (any(dk)) {
    w <- which(dk)
    l2 <- mtlba_lookup(rep(2L, length(w)), rep(1, length(w)), K, rule)
    idx <- as.vector(rbind(2 * w - 1, 2 * w))
    ld2 <- mtlba_logkernel(rt[w], pars[idx, , drop = FALSE], Tm[idx, , drop = FALSE],
                           l2$s, l2$j, l2$lo)
    ld[w] <- mtlba_lse(ld[w], ld2)
  }
  ld[is.na(ld)] <- -Inf
  ld
}

# Density of the uncorrelated MTLBA. rt, R (factor or integer 1/2) and RR (1 =
# lowest) per trial; pars has two rows per trial (accumulators in level order)
# with columns v, sv, A, b, t0 and natural-scale d1 .. d{K-1} (proportions of
# B = b - A above A). Don't-know trials: RR = 1 (R is ignored).
dMTLBA <- function(rt, R, RR, pars, n_ratings = 3, rule = "BE", dont_know = FALSE,
                   log = FALSE) {
  K <- rating_check_n(n_ratings)
  rule <- match.arg(rule, c("BE", "TC"))
  rating_check_dont_know(K, dont_know)
  ld <- mtlba_logdens(rt, R, RR, pars, K, rule, dont_know)
  if (log) ld else exp(ld)
}

# Simulate the uncorrelated MTLBA: lR the accumulator factor (two rows per
# trial), pars as for dMTLBA. Returns R, rt and RR; don't-know responses are
# written as RR = 1 with R = the first level.
rMTLBA <- function(lR, pars, n_ratings = 3, rule = "BE", dont_know = FALSE,
                   ok = rep(TRUE, nrow(pars))) {
  K <- rating_check_n(n_ratings)
  rule <- match.arg(rule, c("BE", "TC"))
  rating_check_dont_know(K, dont_know)
  if (length(levels(lR)) != 2) stop("MTLBA simulates two accumulators only")
  mtlba_check_rho(pars)
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

# Trial log-likelihoods over a (possibly compressed) dadm with two rows per
# trial, expanded to the uncompressed trials and floored at min_ll.
trial_ll_mtlba <- function(pars, dadm, K, rule, dont_know, min_ll = log(1e-10)) {
  n <- nrow(pars) / 2
  first <- seq(1, 2 * n, by = 2)
  ll <- mtlba_logdens(dadm$rt[first], as.integer(dadm$R[first]), dadm$RR[first], pars, K,
                      rule, dont_know)
  # a parameter out of bounds on either accumulator row voids the trial
  ok <- attr(pars, "ok")
  if (!is.null(ok)) ll[!(ok[first] & ok[first + 1])] <- min_ll
  ll[is.na(ll)] <- min_ll
  pmax(min_ll, ll[attr(dadm, "expand")])
}

log_likelihood_mtlba <- function(pars, dadm, K, rule, dont_know, min_ll = log(1e-10))
  sum(trial_ll_mtlba(pars, dadm, K, rule, dont_know, min_ll))

# Refuses censoring, truncation and missing responses (not yet defined).
mtlba_check_data <- function(data, K, dont_know) {
  if (length(levels(data$R)) != 2) stop("MTLBA needs exactly two response levels")
  if (any(is.na(data$R)) || any(!is.finite(data$rt)))
    stop("MTLBA does not yet support missing responses or censored rts")
  if (("LT" %in% names(data) && any(data$LT > 0, na.rm = TRUE)) ||
      ("UT" %in% names(data) && any(is.finite(data$UT))))
    stop("MTLBA does not yet support truncation")
  rating_check_data(data, K, dont_know)
}

# prepare_design: the rating checks, and rho stays at 0 for now (design()
# passes every p_type left out of the formula as `~ 1` with its default as a
# constant, so rho must be a constant equal to 0).
mtlba_prepare_design <- function(formula, constants, Rlevels = NULL, ...) {
  out <- rating_prepare_design(n_choices = 2, trial_level = "rho")(formula, constants,
                                                                   Rlevels, ...)
  if (is.null(constants) || !"rho" %in% names(constants) || constants[["rho"]] != 0)
    stop("MTLBA: correlated drift rates (rho) are not implemented yet; leave rho out of the design")
  out
}

#' The Multiple-Threshold Linear Ballistic Accumulator (Rating) Model
#'
#' Model file for the multiple-threshold LBA (MTLBA), a model of choice,
#' response time and a rating (e.g. confidence). Two LBA accumulators race;
#' each has, besides its choice threshold `b = A + B`, `n_ratings - 1`
#' confidence thresholds between `A` and `b`. This version has independent
#' drift rates (the correlation `rho` is fixed at 0), both truncated to be
#' positive, so with one rating it is the standard [LBA()].
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
#' | *rho* | probit on (-1, 1) | \[-1, 1\] | 0 | Drift-rate correlation; fixed at 0 in this version |
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
#' `t0` is taken from the accumulator that stops the race. Censoring and
#' truncation are not yet supported.
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
