# Chain diagnostics for the stop rules: chain culling (a chain, or one
# subject in one chain, held in a worse mode than the other chains) and the
# flat test (the stage is still drifting). Both read the stored draws only.
#
# Mode comparisons use the log posterior of a chain's state up to the
# hyperpriors: the subject log-likelihoods plus the group density of every
# subject's parameters under that chain's own group mean and covariance.
# Comparing likelihoods alone picks the mode that fits the data best even
# when it does so with an implausible group covariance; on the hard fits
# used to calibrate this the two criteria chose opposite chains.

cull_defaults <- function(cull){
  defaults <- list(persistence = 2L, window = NULL, min_window = 500L, chain_gap = 10, chain_z = 2.5,
                   max_culls = 2L, reburn = 100L)
  if(isTRUE(cull)) return(defaults)
  if(is.null(cull) || isFALSE(cull)) return(NULL)
  if(!is.list(cull)) stop("`cull` must be TRUE, FALSE or a list")
  if(!is.null(cull$window)) defaults$window <- cull$window
  bad <- setdiff(names(cull), names(defaults))
  if(length(bad)) stop("unknown `cull` field(s): ", paste(bad, collapse = ", "),
                       ". Valid fields: ", paste(names(defaults), collapse = ", "))
  modifyList(defaults, cull)
}

# The iterations of `stage` the comparison uses: the last `window` draws of
# every chain, or with window = NULL the last two thirds of the stage once it
# holds at least `min_window` draws (NULL before that). Chains wander over a
# few hundred iterations; a short window mistakes that for a held-apart mode.
stage_window <- function(emc, stage, window = NULL, min_window = 500){
  its <- lapply(emc, function(x) which(x$samples$stage[seq_len(x$samples$idx)] %in% stage))
  n <- min(lengths(its)); if(n == 0) return(NULL)
  if(is.null(window)){
    if(n < min_window) return(NULL)
    window <- max(min_window, floor(2 * n / 3))
  }
  n <- min(n, window)
  lapply(its, function(it) it[(length(it) - n + 1):length(it)])
}

group_logdens <- function(alpha, mu, Sigma){
  # alpha: parameters x subjects; log N(alpha_s | mu, Sigma) per subject
  R <- tryCatch(chol(Sigma), error = function(e) NULL)
  if(is.null(R)) return(rep(NA_real_, ncol(alpha)))
  z <- backsolve(R, alpha - mu, transpose = TRUE)
  -0.5 * colSums(z^2) - sum(log(diag(R))) - 0.5 * nrow(alpha) * log(2 * pi)
}

uses_group_density <- function(emc){
  emc[[1]]$type %in% c("standard", "diagonal", "blocked", "diagonal-gamma")
}

# Per-chain matrices (subjects x draws) over the window: the subject
# log-likelihood `ll`; the subject log posterior under the chain's own group
# parameters `lp` (for the subject-level comparison); and `common`, the subject
# log posterior under one group model shared by all chains (the mean of their
# window means), whose sum compares chains without the wandering of each
# chain's own group covariance entering the comparison.
chain_window_stats <- function(emc, stage = "sample", window = NULL, min_window = 500){
  its <- stage_window(emc, stage, window, min_window)
  if(is.null(its)) return(NULL)
  group <- uses_group_density(emc)
  mu_star <- Sigma_star <- NULL
  if(group){
    mu_star <- rowMeans(sapply(seq_along(emc), function(i) rowMeans(emc[[i]]$samples$theta_mu[, its[[i]], drop = FALSE])))
    Sigma_star <- Reduce(`+`, lapply(seq_along(emc), function(i)
      apply(emc[[i]]$samples$theta_var[, , its[[i]], drop = FALSE], 1:2, mean))) / length(emc)
  }
  lapply(seq_along(emc), function(i){
    s <- emc[[i]]$samples; it <- its[[i]]
    ll <- s$subj_ll[, it, drop = FALSE]
    lp <- common <- ll
    if(group){
      for(k in seq_along(it)){
        a <- s$alpha[, , it[k]]
        lp[, k] <- ll[, k] + group_logdens(a, s$theta_mu[, it[k]], s$theta_var[, , it[k]])
        common[, k] <- ll[, k] + group_logdens(a, mu_star, Sigma_star)
      }
    }
    list(it = it, ll = ll, lp = lp, common = common)
  })
}

# Robust comparison of chain medians: gap to the best chain in units of the
# pooled within-chain spread (median of the chains' MADs)
robust_gaps <- function(med, mad){
  scale <- pmax(apply(mad, 1, median), 1e-8)
  best <- apply(med, 1, max)
  list(gap = best - med, z = (best - med) / scale, best = apply(med, 1, which.max), scale = scale)
}

#' Chain diagnostics: stuck chains and drift
#'
#' Compares the chains of a fitted `emc` object over the last `window`
#' iterations of a stage and reports (i) chains whose log posterior (subject
#' log-likelihoods plus the group density of the subject parameters, under one
#' group mean and covariance shared by all chains so that a chain's own
#' covariance draws do not enter the comparison) is held below the other chains', (ii)
#' the subject spread: per subject, the range over chains of the median
#' log-likelihood (reported only; nothing acts on it), and (iii) the flat test: the drift of the
#' group means, the subject parameters, the summed log-likelihood and the
#' per-subject log-likelihoods between the first and last thirds of the stage,
#' in units of their interquartile range. These are the statistics the
#' `cull` and `max_flat_loc` stop criteria of [fit()] act on.
#'
#' @param emc An emc object
#' @param stage The stage to diagnose (default the last stage run)
#' @param window Number of iterations at the end of the stage the mode
#'   comparison uses; `NULL` (default) uses the last two thirds of the stage
#'   once it holds 500 iterations
#' @param cull A list of thresholds as for the `cull` stop criterion, or `TRUE`
#'   for the defaults
#' @param flat_selection,flat_p1,flat_p2,flat_quantile The flat test's settings,
#'   as in `stop_criteria`
#' @return A list with elements `chains` (one row per chain: median log
#'   posterior, gap to the best chain, robust z, flagged), `subjects` (one row
#'   per subject, widest first: spread, the chains with the highest and lowest
#'   median log-likelihood), `stop` (the sample stop statistics: min ESS, chain
#'   log-likelihood gap, subject spread, projected Rhats; see [fit()]), `flat` (drift
#'   statistics per selection) and `culls` (the record of culls in the object).
#' @export
chain_diagnostics <- function(emc, stage = get_last_stage(emc), window = NULL, cull = TRUE,
                              flat_selection = c("mu", "alpha", "LL", "subj_ll"), flat_p1 = 1/3, flat_p2 = 1/3,
                              flat_quantile = c(mu = 1, alpha = .99, LL = 1, subj_ll = .99)){
  emc <- restore_duplicates(emc)
  cull <- cull_defaults(cull); if(is.null(cull)) cull <- cull_defaults(TRUE)
  if(!is.null(window)) cull$window <- window
  d <- cull_diagnostics(emc, stage, cull)
  f <- flat_diagnostics(emc, stage, flat_selection, flat_p1, flat_p2, flat_quantile)
  out <- list(chains = d$chains, subjects = d$subjects, flat = f, stop = sample_stop_stats(emc, stage), culls = emc[[1]]$culls)
  class(out) <- "emc.diagnostics"
  out
}

#' @export
print.emc.diagnostics <- function(x, ...){
  if(all(is.na(x$chains$gap))) cat("Chains: the stage is too short for the mode comparison (fewer than 500 iterations).\n") else
  cat("Chains (log posterior under the common group model over the window; gap to the best chain, z in within-chain spread):\n")
  print(x$chains, row.names = FALSE, digits = 3)
  if(nrow(x$subjects)) { s <- x$subjects
    cat(sprintf("\nSubject spread (range over chains of the median log-likelihood): %d subject(s) above 5, widest %s at %.1f\n",
                sum(s$spread > 5), s$subject[1], s$spread[1]))
    print(utils::head(s[s$spread > 5, , drop = FALSE], 10), row.names = FALSE, digits = 3) }
  if(length(x$flat)) { cat("\nFlat test (drift between the first and last thirds, in IQRs; mcse in the same units;\n  ess: effective size of the worst series over the window, summed over chains):\n")
    print(do.call(rbind, lapply(names(x$flat), function(n) data.frame(selection = n, drift = x$flat[[n]]$drift, mcse = x$flat[[n]]$mcse, ess = x$flat[[n]]$ess, worst = x$flat[[n]]$worst))), row.names = FALSE, digits = 3) }
  if(!is.null(x$culls)) { cat("\nCulls applied:\n"); print(do.call(rbind, lapply(x$culls, as.data.frame)), row.names = FALSE) }
  invisible(x)
}

# Chain-level statistics over the window, and the subject spread: per subject,
# the range over chains of the median log-likelihood. The spread is reported,
# not acted on: a subject in a different mode in some chains is invisible in
# the summed statistic (tens of units against a sum over hundreds of
# subjects), and moving it between modes needs a proper Metropolis move, not
# a copy of another chain's parameters.
cull_diagnostics <- function(emc, stage, cull){
  st <- chain_window_stats(emc, stage, cull$window, cull$min_window)
  nch <- length(emc)
  empty <- list(chains = data.frame(chain = seq_len(nch), median_lp = NA, gap = NA, z = NA, flagged = FALSE),
                subjects = data.frame(subject = character(0), spread = numeric(0), best = integer(0), worst = integer(0)))
  if(is.null(st) || length(st[[1]]$it) < 20 || nch < 2) return(empty)
  # chains: summed log posterior per draw under the common group model
  tot <- sapply(st, function(s) colSums(s$common))
  med <- matrix(apply(tot, 2, median), 1); mad <- matrix(apply(tot, 2, mad), 1)
  g <- robust_gaps(med, mad)
  chains <- data.frame(chain = seq_len(nch), median_lp = as.vector(med), gap = as.vector(g$gap), z = as.vector(g$z),
                       flagged = as.vector(g$gap > cull$chain_gap & g$z > cull$chain_z))
  # never cull all chains but one on a single check: at least two survive
  if(sum(chains$flagged) > nch - 2){
    keep <- order(chains$gap)[seq_len(2)]
    chains$flagged[keep] <- FALSE
  }
  list(chains = chains, subjects = subject_spread(st, rownames(emc[[1]]$samples$subj_ll)))
}

# Per subject: the range over chains of the median log-likelihood over the
# window, with the chains at either end; widest first
subject_spread <- function(st, subjects){
  sm <- sapply(st, function(s) apply(s$ll, 1, median))
  if(is.null(dim(sm))) sm <- matrix(sm, 1)
  out <- data.frame(subject = subjects, spread = apply(sm, 1, function(v) diff(range(v))),
                    best = apply(sm, 1, which.max), worst = apply(sm, 1, which.min))
  out[order(-out$spread), , drop = FALSE]
}

# ---- flat test ---------------------------------------------------------------

# Drift of one series: |median(first third) - median(last third)| over the
# median of the two thirds' IQRs, and the Monte Carlo error of that
# difference in the same units
# One series pooled over `nch` chains iteration-major (chain varies fastest).
# drift: |median(first part) - median(last part)| in IQRs. mcse: the Monte Carlo
# error of that difference in the same units, with the effective sizes summed
# over the chains' own sub-series (an interleaved series misleads the ESS
# estimator). ess: the effective size of the whole window, summed over chains;
# below about 20 per chain the window holds too few autocorrelation times to
# tell a slow wave from a trend, and the mcse says so.
flat_series <- function(x, p1, p2, nch = 1){
  n <- length(x); if(n < 10 * nch) return(c(drift = NA, mcse = NA, ess = NA))
  a <- x[seq_len(max(1, round(n * p1)))]; b <- x[(n - max(1, round(n * p2)) + 1):n]
  iqr <- max(median(c(stats::IQR(a), stats::IQR(b))), 1e-8)
  ess1 <- function(v) tryCatch(as.numeric(coda::effectiveSize(v)), error = function(e) length(v))
  ess <- function(v) max(sum(vapply(seq_len(nch), function(c) ess1(v[seq(c, length(v), by = nch)]), 0)), 1)
  mcse <- sqrt((1.25 * sd(a))^2 / ess(a) + (1.25 * sd(b))^2 / ess(b)) / iqr
  c(drift = abs(median(a) - median(b)) / iqr, mcse = mcse, ess = ess(x))
}

# Series pooled across chains, iteration-major (the first third of every
# chain against the last third of every chain)
pooled_matrix <- function(emc, selection, stage){
  its <- stage_window(emc, stage, .Machine$integer.max); if(is.null(its)) return(NULL)
  n <- length(its[[1]])
  if(selection == "mu"){
    if(is.null(emc[[1]]$samples$theta_mu)) return(NULL)
    arr <- sapply(seq_along(emc), function(i) emc[[i]]$samples$theta_mu[, its[[i]], drop = FALSE], simplify = "array") # p x n x chains
    m <- matrix(aperm(arr, c(3, 2, 1)), ncol = dim(arr)[1]); colnames(m) <- rownames(emc[[1]]$samples$theta_mu)
  } else if(selection == "alpha"){
    arr <- sapply(seq_along(emc), function(i){ a <- emc[[i]]$samples$alpha[, , its[[i]], drop = FALSE]; matrix(a, prod(dim(a)[1:2])) }, simplify = "array") # (p*s) x n x chains
    m <- matrix(aperm(arr, c(3, 2, 1)), ncol = dim(arr)[1])
    colnames(m) <- as.vector(outer(dimnames(emc[[1]]$samples$alpha)[[1]], dimnames(emc[[1]]$samples$alpha)[[2]], paste, sep = ":"))
  } else if(selection == "LL"){
    arr <- sapply(seq_along(emc), function(i) colSums(emc[[i]]$samples$subj_ll[, its[[i]], drop = FALSE]))
    m <- matrix(t(arr), ncol = 1); colnames(m) <- "LL"
  } else if(selection == "subj_ll"){
    arr <- sapply(seq_along(emc), function(i) emc[[i]]$samples$subj_ll[, its[[i]], drop = FALSE], simplify = "array") # subjects x n x chains
    m <- matrix(aperm(arr, c(3, 2, 1)), ncol = dim(arr)[1]); colnames(m) <- rownames(emc[[1]]$samples$subj_ll)
  } else stop("flat_selection must be one or more of: mu, alpha, LL, subj_ll")
  m
}

flat_diagnostics <- function(emc, stage, selection, p1, p2, quantile){
  out <- list(); nch <- length(emc)
  for(sel in selection){
    m <- pooled_matrix(emc, sel, stage); if(is.null(m) || nrow(m) < 10) next
    if(sel %in% c("alpha", "subj_ll")){ # medians and IQRs vectorised; the Monte Carlo error from a sample of columns
      n <- nrow(m); a <- m[seq_len(max(1, round(n * p1))), , drop = FALSE]; b <- m[(n - max(1, round(n * p2)) + 1):n, , drop = FALSE]
      qa <- apply(a, 2, quantile, c(.25, .5, .75)); qb <- apply(b, 2, quantile, c(.25, .5, .75))
      iqr <- pmax(apply(rbind(qa[3, ] - qa[1, ], qb[3, ] - qb[1, ]), 2, median), 1e-8)
      drift <- abs(qa[2, ] - qb[2, ]) / iqr
      cols <- unique(c(which.max(drift), round(seq(1, ncol(m), length.out = min(50, ncol(m))))))
      fs <- apply(m[, cols, drop = FALSE], 2, flat_series, p1 = p1, p2 = p2, nch = nch)
      mcse <- median(fs["mcse", ], na.rm = TRUE); ess <- fs["ess", 1]
    } else {
      fs <- apply(m, 2, flat_series, p1 = p1, p2 = p2, nch = nch); drift <- fs["drift", ]
      mcse <- max(fs["mcse", ], na.rm = TRUE); ess <- fs["ess", which.max(drift)]
    }
    q <- if(!is.null(quantile[sel])) quantile[[sel]] else 1
    stat <- if(q >= 1) max(drift, na.rm = TRUE) else unname(stats::quantile(drift, q, na.rm = TRUE))
    out[[sel]] <- list(drift = stat, mcse = mcse, ess = unname(ess), worst = colnames(m)[which.max(drift)], quantile = q, all = drift)
  }
  out
}

# TRUE when every selection's drift is within max_flat_loc, or (flat_mcse > 0)
# within flat_mcse Monte Carlo errors of zero
flat_done <- function(flat, max_flat_loc, flat_mcse){
  if(!length(flat)) return(TRUE)
  all(vapply(flat, function(f) !is.finite(f$drift) || f$drift <= max_flat_loc || f$drift <= flat_mcse * f$mcse, logical(1)))
}

# ---- culling: respawn ----------------------------------------------------------

# Overwrite the last stored draw of chain `to` with the state of chain `from`
# at a draw of its window (a random one at least `back` draws from the end
# where the window allows), so the next run starts the chain there. The
# donor's kernel settings come too (the joint move's step and location scales,
# the particle settings): settings adapted where the culled chain was stuck do
# not fit the region it restarts in. Its random number state stays its own.
respawn_chain <- function(emc, to, from, its_from, back = 100){
  pick <- its_from[seq_len(max(1, length(its_from) - back))]
  t <- pick[sample.int(length(pick), 1)]
  s_to <- emc[[to]]$samples; s_from <- emc[[from]]$samples; last <- s_to$idx
  for(key in setdiff(names(s_from), c("stage", "idx", "last_theta_var_inv"))){
    x <- s_from[[key]]; if(!is.array(x) && !is.matrix(x)) next
    d <- length(dim(x)); if(d < 2 || dim(x)[d] < t) next
    y <- s_to[[key]]
    if(d == 2) y[, last] <- x[, t] else if(d == 3) y[, , last] <- x[, , t] else if(d == 4) y[, , , last] <- x[, , , t]
    s_to[[key]] <- y
  }
  if(!is.null(s_from$theta_var)){
    V <- s_from$theta_var[, , t]
    s_to$last_theta_var_inv <- tryCatch(chol2inv(chol(V)), error = function(e) MASS::ginv(V))
  }
  for(a in c("scale_move", "pm_settings")) attr(s_to, a) <- attr(s_from, a)
  emc[[to]]$samples <- s_to
  emc
}

# Apply the whole-chain culls that have persisted. Returns the object and a
# record of what was done.
apply_culls <- function(emc, diag, stage, cull, state){
  its <- stage_window(emc, stage, cull$window, cull$min_window)
  done <- list()
  culled_chains <- diag$chains$chain[diag$chains$flagged & state$chain_count[diag$chains$chain] >= cull$persistence]
  donors <- setdiff(seq_along(emc), diag$chains$chain[diag$chains$flagged])
  best <- donors[which.max(diag$chains$median_lp[donors])]
  for(cc in culled_chains){
    from <- if(length(donors) > 1) donors[sample.int(length(donors), 1)] else best
    emc <- respawn_chain(emc, cc, from, its[[from]], cull$reburn)
    done[[length(done) + 1]] <- list(iteration = emc[[1]]$samples$idx, stage = stage, chain = cc, subject = NA_character_,
                                     donor = from, gap = diag$chains$gap[cc], z = diag$chains$z[cc])
  }
  list(emc = emc, done = done)
}

# Persistence bookkeeping across checks: consecutive flags per chain
update_cull_state <- function(state, diag, nch){
  if(is.null(state)) state <- list(chain_count = integer(nch))
  f <- diag$chains$flagged
  state$chain_count <- ifelse(f, state$chain_count + 1L, 0L)
  state
}

#' Thin a Fitted Object to a Target Effective Sample Size
#'
#' Keeps ``factor`` times ``target`` draws of ``stage`` pooled over the chains,
#' evenly spaced, so the thinning interval is the pooled stage length over
#' ``factor * target``. A parameter whose effective sample size is below
#' ``factor * target`` keeps nearly all of it (its autocorrelation time exceeds the
#' interval), one at ``target`` keeps about nine tenths for the default factor of
#' 2, and better-mixed parameters are capped near ``factor * target``. Nothing is
#' removed when the stage is already that short. The other stages are kept as they
#' are. The interval and the effective sizes of ``selection`` before and after are
#' recorded in ``emc[[1]]$thinning``.
#'
#' @param emc An emc object
#' @param target Numeric. The effective sample size to retain (default 3000)
#' @param stage Character. The stage to thin (default "sample")
#' @param selection Character vector. Parameter types whose effective sizes are reported
#' @param factor Numeric. Draws kept per unit of ``target`` (default 2)
#' @return The thinned emc object
#' @examples
#' thinned <- thin_to_ess(samples_LNR, target = 20)
#' thinned[[1]]$thinning
#' @export
thin_to_ess <- function(emc, target = 3000, stage = "sample", selection = c("mu", "alpha"), factor = 2){
  if(!is.numeric(target) || length(target) != 1 || target <= 0) stop("`target` must be a single positive number")
  if(!is.numeric(factor) || length(factor) != 1 || factor < 1) stop("`factor` must be a single number of at least 1")
  selection <- selection[vapply(selection, function(s) !is.null(emc[[1]]$samples[[switch(s, mu = "theta_mu", Sigma = "theta_var", s)]]), logical(1))]
  ess_of <- function(x) { v <- unlist(lapply(selection, function(s) ess_summary(x, selection = s, stage = stage, stat = NULL)))
    if(!length(v)) return(c(min = NA, median = NA, max = NA)); c(min = min(v), median = stats::median(v), max = max(v)) }
  n <- unname(chain_n(emc)[1, stage]); nch <- length(emc)
  keep <- min(n, ceiling(factor * target / nch))
  before <- ess_of(emc)
  if(keep < n) emc <- subset(emc, stage = stage, length.out = keep, keep_stages = TRUE)
  emc[[1]]$thinning <- list(stage = stage, target = target, factor = factor, n_before = n, n_after = keep,
                            interval = n / keep, ess_before = before, ess_after = if(keep < n) ess_of(emc) else before)
  emc
}

