# Statistics read from the stored draws for hardfit() and diagnose(): the stop
# statistics hardfit() reads next to the Rhat rule at each sample check, on the
# draws it keeps (after any first-third discard), and the drift test between
# the first and last thirds of a stage. diagnose() reads the same ones on a
# finished fit. Each stop statistic can be a criterion (see ?hardfit); rhat_b,
# rhat_ident and the drift are reported only.
# Why an ESS floor: on a ridge the group means mix slowly enough that the Rhat
# rule passes with tens of effective draws. Why a per-subject spread: a subject
# in a different mode in some chains moves the summed log-likelihood by too
# little to see. Why two projected Rhats: along the direction the chains
# disagree on most rhat_b reproduces the largest per-parameter Rhat; along the
# best-identified direction rhat_ident says whether a fit stopped at the limit
# is usable off the ridge.

min_stat_draws <- 20L   # the fewest draws a median comparison is read on
scale_floor <- 1e-8     # guards a zero MAD or IQR in a ratio
spread_note <- 5        # subjects whose chains' medians differ by more than this are counted

# The rule's defaults by type. "single" has one subject, so the rule is the
# maximum over its parameters (no quantile) and the threshold is 1.05: in
# single-subject replays 1.1 passed chains that still sat in different modes.
rule_defaults <- function(type, max_gd = NULL, gd_quantile = NULL){
  single <- type == "single"
  list(max_gd = if(!is.null(max_gd)) max_gd else if(single) 1.05 else 1.1,
       gd_quantile = if(!is.null(gd_quantile)) gd_quantile else if(single) 1 else .99,
       selection = if(single) "alpha" else c("alpha", "mu"))
}

# The three criteria on a set of stop statistics: NULL is off; a statistic that
# cannot be computed never passes
stop_rules <- function(st, min_ess, max_ll_gap, max_subj_spread){
  ok <- function(crit, v, below){
    if(is.null(crit)) return(TRUE)
    if(is.null(v) || !is.finite(v)) return(FALSE)
    if(below) v <= crit else v >= crit
  }
  c(ess = ok(min_ess, st$ess_min, FALSE), gap = ok(max_ll_gap, st$ll_gap, TRUE), spread = ok(max_subj_spread, st$subj_spread, TRUE))
}

stat_or_na <- function(x) if(is.null(x) || !is.finite(x)) NA_real_ else x

# The stored iterations of `stage`, per chain
stage_iters <- function(emc, stage) lapply(emc, function(x) which(x$samples$stage[seq_len(x$samples$idx)] %in% stage))

# The iterations of `stage` a chain comparison uses: the last `window` draws of
# every chain, or with window = NULL the last two thirds of the stage once it
# holds at least `min_window` draws (NULL before that). Chains wander over a
# few hundred iterations; a short window mistakes that for a held-apart mode.
stage_window <- function(emc, stage, window = NULL, min_window = 500){
  its <- stage_iters(emc, stage)
  n <- min(lengths(its)); if(n == 0) return(NULL)
  if(is.null(window)){
    if(n < min_window) return(NULL)
    window <- max(min_window, floor(2 * n / 3))
  }
  n <- min(n, window)
  lapply(its, function(it) it[(length(it) - n + 1):length(it)])
}

# The subject parameters at iterations `its` as one (parameter x subject) by
# draws matrix, rows named par:subject
alpha_matrix <- function(samples, its){
  a <- samples$alpha[, , its, drop = FALSE]
  matrix(a, prod(dim(a)[1:2]),
         dimnames = list(as.vector(outer(dimnames(a)[[1]], dimnames(a)[[2]], paste, sep = ":")), NULL))
}

# Per-chain matrices (subjects x draws) over the window `its` (one index vector
# per chain): the subject log-likelihood `ll` (the comparison) and the subject
# log posterior under the chain's own group parameters `lp` (the guard on a
# cull). Only the types whose group level is N(mu, Sigma) over alpha have the
# group term; SEM and the factor types store an all-NA theta_var, so lp = ll.
chain_window_stats <- function(emc, its){
  if(is.null(its)) return(NULL)
  group <- emc[[1]]$type %in% c("standard", "diagonal", "blocked", "diagonal-gamma")
  lapply(seq_along(emc), function(i){
    s <- emc[[i]]$samples; it <- its[[i]]
    ll <- s$subj_ll[, it, drop = FALSE]; lp <- ll
    if(group){
      p <- dim(s$alpha)[1]
      for(k in seq_along(it)){
        dens <- mvtnorm::dmvnorm(t(matrix(s$alpha[, , it[k]], p)), s$theta_mu[seq_len(p), it[k]],
                                 matrix(s$theta_var[, , it[k]], p, p), log = TRUE, checkSymmetry = FALSE)
        dens[!is.finite(dens)] <- NA   # a group covariance draw that is not positive definite
        lp[, k] <- ll[, k] + dens
      }
    }
    list(it = it, ll = ll, lp = lp)
  })
}

# The chain comparison over a window: per chain the median summed subject
# log-likelihood, its gap to the best chain in units and in within-chain spreads
# (z, the pooled MAD), the median own log posterior and whether that is below
# the other chains' median; per subject the range over chains of the median
# log-likelihood with the chains at either end, widest first. NULL below
# min_stat_draws draws. A non-finite log-likelihood anywhere in the window
# makes the comparison NA: no criterion passes on it and nothing is culled.
chain_table <- function(cw, subjects){
  if(is.null(cw) || length(cw[[1]]$it) < min_stat_draws || length(cw) < 2) return(NULL)
  nch <- length(cw)
  if(!all(vapply(cw, function(s) all(is.finite(s$ll)), TRUE)))
    return(list(chains = data.frame(chain = seq_len(nch), median_ll = NA_real_, gap = NA_real_, z = NA_real_, median_lp = NA_real_, lp_lower = NA),
                subjects = data.frame(subject = subjects, spread = NA_real_, best = NA_integer_, worst = NA_integer_)))
  tot <- sapply(cw, function(s) colSums(s$ll))                 # draws x chains
  med <- apply(tot, 2, stats::median); gap <- max(med) - med
  scale <- max(stats::median(apply(tot, 2, stats::mad)), scale_floor)
  lp <- apply(sapply(cw, function(s) colSums(s$lp)), 2, stats::median)
  lp_lower <- vapply(seq_len(nch), function(i) isTRUE(lp[i] < stats::median(lp[-i])), TRUE)
  sm <- sapply(cw, function(s) apply(s$ll, 1, stats::median)); if(is.null(dim(sm))) sm <- matrix(sm, 1)
  subj <- data.frame(subject = subjects, spread = apply(sm, 1, function(v) diff(range(v))),
                     best = apply(sm, 1, which.max), worst = apply(sm, 1, which.min))
  list(chains = data.frame(chain = seq_len(nch), median_ll = med, gap = gap, z = gap / scale, median_lp = lp, lp_lower = lp_lower),
       subjects = subj[order(-subj$spread), , drop = FALSE])
}

# The stop statistics on the last `n` draws of `stage` (n = the shortest chain):
# the smallest effective size over the group means and log group SDs (type
# "single": the subject parameters), each summed over chains; the projected
# Rhats; and the chain comparison from `tab` (by default the last two thirds).
# NULL below min_stat_draws draws or with one chain.
sample_stop_stats <- function(emc, stage = "sample",
                              tab = chain_table(chain_window_stats(emc, stage_window(emc, stage, NULL, min_stat_draws)),
                                                rownames(emc[[1]]$samples$subj_ll))){
  its <- stage_iters(emc, stage)
  n <- min(lengths(its)); nch <- length(emc)
  if(n < min_stat_draws || nch < 2) return(NULL)
  its <- lapply(its, function(it) utils::tail(it, n))
  s1 <- emc[[1]]$samples
  # a parameter whose series is non-finite or constant in any chain (e.g. the group
  # variances of a type without them, such as SEM's all-NA theta_var) has no
  # effective size and is left out of the minimum
  chain_sum <- function(L, p) { e <- vapply(L, function(m) ess_series(m[p, ]), 0); if(anyNA(e)) NA_real_ else sum(e) }
  # ess_summary() has no log group SD selection; this keeps the three series on one path
  if(!is.null(s1$theta_mu) && !is.null(s1$theta_var)){
    MU <- lapply(seq_len(nch), function(i) emc[[i]]$samples$theta_mu[, its[[i]], drop = FALSE])
    pn <- rownames(MU[[1]])
    LSD <- lapply(seq_len(nch), function(i){
      V <- emc[[i]]$samples$theta_var[, , its[[i]], drop = FALSE]
      matrix(apply(V, 3, function(v) 0.5 * log(diag(v))), nrow = length(pn), dimnames = list(pn, NULL))
    })
    ess <- c(stats::setNames(vapply(pn, function(p) chain_sum(MU, p), 0), pn),
             stats::setNames(vapply(pn, function(p) chain_sum(LSD, p), 0), paste0("lsd:", pn)))
  } else {
    MU <- lapply(seq_len(nch), function(i) alpha_matrix(emc[[i]]$samples, its[[i]]))
    pn <- rownames(MU[[1]])
    ess <- stats::setNames(vapply(pn, function(p) chain_sum(MU, p), 0), pn)
  }
  out <- list()
  if(all(is.na(ess))){ out$ess_min <- NA_real_; out$ess_min_par <- NA_character_ }
  else { out$ess_min <- min(ess, na.rm = TRUE); out$ess_min_par <- names(ess)[which.min(ess)] }
  out <- c(out, projected_rhats(MU))
  gap <- if(is.null(tab)) NA_real_ else tab$chains$gap
  out$ll_gap <- max(gap); out$ll_z <- if(is.null(tab)) NA_real_ else max(tab$chains$z)
  out$ll_worst_chain <- if(all(is.na(gap))) NA_integer_ else which.max(gap)
  sp <- if(is.null(tab)) NA_real_ else tab$subjects$spread
  out$subj_spread <- sp[1]; out$subj_spread_n5 <- sum(sp > spread_note)
  out$subj_spread_subject <- if(is.null(tab)) NA_character_ else tab$subjects$subject[1]
  out
}

# Effective size of one chain's series (coda's spectral estimate); NA for a series
# that is too short, non-finite anywhere or constant (nothing to estimate)
ess_series <- function(x){
  if(length(x) < 4 || !all(is.finite(x))) return(NA_real_)
  s <- stats::sd(x)
  if(!is.finite(s) || s == 0) return(NA_real_)
  unname(coda::effectiveSize(x))
}

# Rhat of the group means projected on the direction of largest between-chain
# disagreement and on the pooled within-chain axis of least variance, both in
# units of the pooled within-chain SD. MU: one parameters x draws matrix per chain.
projected_rhats <- function(MU){
  nch <- length(MU); p <- nrow(MU[[1]]); n <- ncol(MU[[1]])
  if(p < 2) return(list(rhat_b = NA_real_, rhat_ident = NA_real_))
  if(!all(vapply(MU, function(m) all(is.finite(m)), TRUE))) return(list(rhat_b = NA_real_, rhat_ident = NA_real_))
  wsd <- sqrt(Reduce(`+`, lapply(MU, function(m) apply(m, 1, stats::var))) / nch)
  if(any(!is.finite(wsd)) || any(wsd <= 0)) return(list(rhat_b = NA_real_, rhat_ident = NA_real_))
  Z <- lapply(MU, function(m) m / wsd)
  Cw <- Reduce(`+`, lapply(Z, function(z) stats::cov(t(z)))) / nch
  ident <- eigen(Cw, symmetric = TRUE)$vectors[, p]
  means <- t(vapply(Z, rowMeans, numeric(p)))
  b <- eigen(stats::cov(means), symmetric = TRUE)$vectors[, 1]
  proj <- function(v) array(vapply(Z, function(z) as.vector(crossprod(v, z)), numeric(n)), c(n, 1, nch))
  list(rhat_b = unname(split_rhat(proj(b))), rhat_ident = unname(split_rhat(proj(ident))))
}

stop_stats_message <- function(st){
  if(is.null(st)) return(NULL)
  paste0(if(!is.null(st$ess_min) && is.finite(st$ess_min)) sprintf(" | min ESS=%.0f (%s)", st$ess_min, st$ess_min_par),
         if(is.finite(st$ll_gap)) sprintf(" | LL gap=%.0f", st$ll_gap) else " | LL gap=NA",
         if(is.finite(st$subj_spread)) sprintf(" | subject spread=%.1f (%s; %d > %s)", st$subj_spread,
                                               st$subj_spread_subject, st$subj_spread_n5, format(spread_note)))
}

# ---- drift test: is the stage still moving? ----

# One series pooled over `nch` chains iteration-major (chain varies fastest).
# drift: |median(first part) - median(last part)| in IQRs. mcse: the Monte Carlo
# error of that difference in the same units, with the effective sizes summed
# over the chains' own sub-series (an interleaved series misleads the ESS
# estimator). ess: the effective size of the whole window, summed over chains;
# below about 20 per chain the window holds too few autocorrelation times to
# tell a slow wave from a trend, and the mcse says so.
flat_series <- function(x, p1, p2, nch = 1){
  n <- length(x); if(n < 10 * nch) return(c(drift = NA, mcse = NA, ess = NA))   # fewer than 10 draws per chain give no IQR and ESS per third
  a <- x[seq_len(max(1, round(n * p1)))]; b <- x[(n - max(1, round(n * p2)) + 1):n]
  iqr <- max(median(c(stats::IQR(a), stats::IQR(b))), scale_floor)
  ess1 <- function(v) tryCatch(as.numeric(coda::effectiveSize(v)), error = function(e) length(v))
  ess <- function(v) max(sum(vapply(seq_len(nch), function(c) ess1(v[seq(c, length(v), by = nch)]), 0)), 1)
  # 1.25 ~ sqrt(pi / 2): the standard error of a median relative to a mean's, for normal data
  mcse <- sqrt((1.25 * sd(a))^2 / ess(a) + (1.25 * sd(b))^2 / ess(b)) / iqr
  c(drift = abs(median(a) - median(b)) / iqr, mcse = mcse, ess = ess(x))
}

# Series pooled across chains, iteration-major (the first third of every
# chain against the last third of every chain)
pooled_matrix <- function(emc, selection, stage){
  its <- stage_iters(emc, stage); if(min(lengths(its)) == 0) return(NULL)
  if(selection == "mu"){
    if(is.null(emc[[1]]$samples$theta_mu)) return(NULL)
    arr <- sapply(seq_along(emc), function(i) emc[[i]]$samples$theta_mu[, its[[i]], drop = FALSE], simplify = "array") # p x n x chains
    m <- matrix(aperm(arr, c(3, 2, 1)), ncol = dim(arr)[1]); colnames(m) <- rownames(emc[[1]]$samples$theta_mu)
  } else if(selection == "alpha"){
    arr <- sapply(seq_along(emc), function(i) alpha_matrix(emc[[i]]$samples, its[[i]]), simplify = "array") # (p*s) x n x chains
    m <- matrix(aperm(arr, c(3, 2, 1)), ncol = dim(arr)[1]); colnames(m) <- rownames(arr)
  } else if(selection == "LL"){
    arr <- sapply(seq_along(emc), function(i) colSums(emc[[i]]$samples$subj_ll[, its[[i]], drop = FALSE]))
    m <- matrix(t(arr), ncol = 1); colnames(m) <- "LL"
  } else if(selection == "subj_ll"){
    arr <- sapply(seq_along(emc), function(i) emc[[i]]$samples$subj_ll[, its[[i]], drop = FALSE], simplify = "array") # subjects x n x chains
    m <- matrix(aperm(arr, c(3, 2, 1)), ncol = dim(arr)[1]); colnames(m) <- rownames(emc[[1]]$samples$subj_ll)
  } else stop("flat selection must be one or more of: mu, alpha, LL, subj_ll")
  m
}

# Drift of the group means, the subject parameters, the summed log-likelihood
# and the per-subject log-likelihoods between the first and last thirds of the
# stage: per selection the maximum (or `quantile`) over its series, its Monte
# Carlo error, the effective size of the worst series and its name.
flat_diagnostics <- function(emc, stage, selection = c("mu", "alpha", "LL", "subj_ll"), p1 = 1/3, p2 = 1/3,
                             quantile = c(mu = 1, alpha = .99, LL = 1, subj_ll = .99)){
  out <- list(); nch <- length(emc)
  for(sel in selection){
    m <- pooled_matrix(emc, sel, stage); if(is.null(m) || nrow(m) < 10) next
    if(sel %in% c("alpha", "subj_ll")){ # medians and IQRs vectorised; the Monte Carlo error from a sample of columns
      n <- nrow(m); a <- m[seq_len(max(1, round(n * p1))), , drop = FALSE]; b <- m[(n - max(1, round(n * p2)) + 1):n, , drop = FALSE]
      qa <- apply(a, 2, stats::quantile, c(.25, .5, .75)); qb <- apply(b, 2, stats::quantile, c(.25, .5, .75))
      iqr <- pmax(apply(rbind(qa[3, ] - qa[1, ], qb[3, ] - qb[1, ]), 2, median), scale_floor)
      drift <- abs(qa[2, ] - qb[2, ]) / iqr
      # 50 columns bound the cost of the Monte Carlo error over p x subjects series
      cols <- unique(c(which.max(drift), round(seq(1, ncol(m), length.out = min(50, ncol(m))))))
      fs <- apply(m[, cols, drop = FALSE], 2, flat_series, p1 = p1, p2 = p2, nch = nch)
      mcse <- median(fs["mcse", ], na.rm = TRUE); ess <- fs["ess", 1]
    } else {
      fs <- apply(m, 2, flat_series, p1 = p1, p2 = p2, nch = nch); drift <- fs["drift", ]
      mcse <- max(fs["mcse", ], na.rm = TRUE); ess <- fs["ess", which.max(drift)]
    }
    q <- if(!is.na(quantile[sel])) quantile[[sel]] else 1
    stat <- if(q >= 1) max(drift, na.rm = TRUE) else unname(stats::quantile(drift, q, na.rm = TRUE))
    out[[sel]] <- list(drift = stat, mcse = mcse, ess = unname(ess), worst = colnames(m)[which.max(drift)], quantile = q)
  }
  out
}

# The drift test is a report: a failure in it must not lose the fit, but must be heard
drift_report <- function(emc, stage){
  tryCatch(flat_diagnostics(emc, stage), error = function(e){
    warning("the drift test failed: ", conditionMessage(e), call. = FALSE); list()
  })
}
