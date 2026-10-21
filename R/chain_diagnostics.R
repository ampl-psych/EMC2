# Chain diagnostics for hardfit() and diagnose(): chain culling (a chain held
# in a different region from the others) and the drift test (the stage is
# still moving). Both read the stored draws only.
#
# Chains are compared on the summed subject log-likelihood. The log posterior
# under pooled group parameters is not used: its group-density term differs
# between healthy chains by hundreds of units (it follows each chain's group
# covariance draws; replay of 105 fits: log-posterior gap median 617,
# log-likelihood gap median 15), so culling on it removed wandering chains,
# not stuck ones. Below a likelihood gap of lp_guard_below (default 1,000) a
# chain is culled only when it is also lower on its own log posterior (the
# subject log-likelihoods plus the group density under that chain's own group
# parameters) than the other chains: the likelihood alone can favour a state
# that fits the data with an implausible group covariance. Above it the
# likelihood decides: the own-covariance group term differs between healthy
# chains by 200-1,000 units, so at smaller gaps it could veto a genuine cull by
# noise, while a stuck chain thousands of units down is far lower on both.

cull_defaults <- function(cull){
  defaults <- list(persistence = 2L, window = NULL, min_window = 500L, chain_gap = 200, chain_z = 5,
                   lp_guard_below = 1000, max_culls = 1L, reburn = 100L)
  if(isTRUE(cull)) return(defaults)
  if(is.null(cull) || isFALSE(cull)) return(NULL)
  if(!is.list(cull)) stop("`cull` must be TRUE, FALSE or a list")
  if(!is.null(cull$window)) defaults$window <- cull$window
  bad <- setdiff(names(cull), names(defaults))
  if(length(bad)) stop("unknown `cull` field(s): ", paste(bad, collapse = ", "),
                       ". Valid fields: ", paste(names(defaults), collapse = ", "))
  modifyList(defaults, cull)
}

# The settings per stage: burn compares over a shorter window (its draws are
# discarded, so a wrong cull is cheap) and may re-seed twice; sample restarts
# the stage on a cull, so once.
cull_settings <- function(cull, stage){
  given <- if(is.list(cull)) names(cull) else character(0)
  cull <- cull_defaults(cull)
  if(is.null(cull)) return(NULL)
  if(stage == "burn"){
    if(!"min_window" %in% given && is.null(cull$window)) cull$min_window <- 250L
    if(!"max_culls" %in% given) cull$max_culls <- 2L
  }
  cull
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
# log-likelihood `ll` (the comparison) and the subject log posterior under the
# chain's own group parameters `lp` (the guard on a cull).
chain_window_stats <- function(emc, stage = "sample", window = NULL, min_window = 500){
  its <- stage_window(emc, stage, window, min_window)
  if(is.null(its)) return(NULL)
  group <- uses_group_density(emc)
  lapply(seq_along(emc), function(i){
    s <- emc[[i]]$samples; it <- its[[i]]
    ll <- s$subj_ll[, it, drop = FALSE]
    lp <- ll
    if(group){
      p <- dim(s$alpha)[1]; ns <- dim(s$alpha)[2]
      for(k in seq_along(it))
        lp[, k] <- ll[, k] + group_logdens(matrix(s$alpha[, , it[k]], p, ns), s$theta_mu[seq_len(p), it[k]],
                                           matrix(s$theta_var[, , it[k]], p, p))
    }
    list(it = it, ll = ll, lp = lp)
  })
}

# Robust comparison of chain medians: gap to the best chain in units of the
# pooled within-chain spread (median of the chains' MADs)
robust_gaps <- function(med, mad){
  scale <- pmax(apply(mad, 1, median), 1e-8)
  best <- apply(med, 1, max)
  list(gap = best - med, z = (best - med) / scale, best = apply(med, 1, which.max), scale = scale)
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
  empty <- list(chains = data.frame(chain = seq_len(nch), median_ll = NA, gap = NA, z = NA, median_lp = NA, lp_lower = NA, flagged = FALSE),
                subjects = data.frame(subject = character(0), spread = numeric(0), best = integer(0), worst = integer(0)))
  if(is.null(st) || length(st[[1]]$it) < 20 || nch < 2) return(empty)
  # chains: summed subject log-likelihood per draw; a chain is flagged when it is
  # chain_gap units and chain_z within-chain spreads below the best chain (a
  # different region: fine chains differ by tens of units, a stuck one by
  # hundreds or thousands) and, when the gap is below lp_guard_below, also lower
  # on its own log posterior than the median of the other chains
  tot <- sapply(st, function(s) colSums(s$ll))
  med <- matrix(apply(tot, 2, median), 1); mad <- matrix(apply(tot, 2, mad), 1)
  g <- robust_gaps(med, mad)
  lp <- apply(sapply(st, function(s) colSums(s$lp)), 2, median)
  lp_lower <- vapply(seq_len(nch), function(i) isTRUE(lp[i] < stats::median(lp[-i])), TRUE)
  chains <- data.frame(chain = seq_len(nch), median_ll = as.vector(med), gap = as.vector(g$gap), z = as.vector(g$z),
                       median_lp = lp, lp_lower = lp_lower,
                       flagged = as.vector(g$gap > cull$chain_gap & g$z > cull$chain_z) &
                         (lp_lower | as.vector(g$gap) >= cull$lp_guard_below))
  chains$flagged[is.na(chains$flagged)] <- FALSE
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

# ---- drift test ------------------------------------------------------------

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
      iqr <- pmax(apply(rbind(qa[3, ] - qa[1, ], qb[3, ] - qb[1, ]), 2, median), 1e-8)
      drift <- abs(qa[2, ] - qb[2, ]) / iqr
      cols <- unique(c(which.max(drift), round(seq(1, ncol(m), length.out = min(50, ncol(m))))))
      fs <- apply(m[, cols, drop = FALSE], 2, flat_series, p1 = p1, p2 = p2, nch = nch)
      mcse <- median(fs["mcse", ], na.rm = TRUE); ess <- fs["ess", 1]
    } else {
      fs <- apply(m, 2, flat_series, p1 = p1, p2 = p2, nch = nch); drift <- fs["drift", ]
      mcse <- max(fs["mcse", ], na.rm = TRUE); ess <- fs["ess", which.max(drift)]
    }
    q <- if(!is.na(quantile[sel])) quantile[[sel]] else 1
    stat <- if(q >= 1) max(drift, na.rm = TRUE) else unname(stats::quantile(drift, q, na.rm = TRUE))
    out[[sel]] <- list(drift = stat, mcse = mcse, ess = unname(ess), worst = colnames(m)[which.max(drift)], quantile = q, all = drift)
  }
  out
}

flat_message <- function(flat){
  if(!length(flat)) return(NULL)
  paste0(" | drift(se)=", paste(vapply(names(flat), function(n) sprintf("%s %.2f(%.2f)", n, flat[[n]]$drift, flat[[n]]$mcse), ""), collapse = ", "))
}

# ---- culling: respawn ------------------------------------------------------

# Overwrite the last stored draw of chain `to` with the state of chain `from`
# at a draw of its window (a random one at least `back` draws from the end
# where the window allows), so the next run starts the chain there. The
# donor's kernel settings come too (the joint moves' step and location scales,
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
  best <- donors[which.max(diag$chains$median_ll[donors])]
  for(cc in culled_chains){
    from <- if(length(donors) > 1) donors[sample.int(length(donors), 1)] else best
    emc <- respawn_chain(emc, cc, from, its[[from]], cull$reburn)
    done[[length(done) + 1]] <- list(iteration = emc[[1]]$samples$idx, stage = stage, chain = cc,
                                     donor = from, gap = diag$chains$gap[cc], z = diag$chains$z[cc], exit = FALSE)
  }
  list(emc = emc, done = done)
}

# At the end of burn (by its rule or at max_tries): re-seed every chain held in
# a different region, with the cull thresholds but no persistence (burn draws
# are discarded, so this is free), over the last two thirds of burn however
# short the stage. The re-seeded chain starts adapt from a donor's state, with
# its kernel settings. Recorded in emc[[1]]$culls with exit = TRUE.
# Note: only the chain's last draw is replaced. The first adapt step's
# proposals (create_chain_proposals(), merged over chains) and lik_prec are
# built from the last min(250, idx / 1.5) draws of every chain, so they still
# pool the re-seeded chain's stuck burn draws (1/n_chains of that history).
# The window slides past them within adapt's first 250 iterations, which
# adapt_converge's minimum guarantees before the sample kernel is built from
# adapt's tail; the same holds after a cull at a burn check.
burn_exit_reseed <- function(emc, cull, verbose = FALSE){
  cx <- cull; cx$window <- NULL; cx$min_window <- 20L
  d <- cull_diagnostics(emc, "burn", cx)
  bad <- d$chains$chain[d$chains$flagged]
  if(!length(bad)) return(emc)
  its <- stage_window(emc, "burn", NULL, cx$min_window)
  donors <- setdiff(seq_along(emc), bad)
  for(cc in bad){
    from <- donors[sample.int(length(donors), 1)]
    emc <- respawn_chain(emc, cc, from, its[[from]], min(cull$reburn, floor(length(its[[from]]) / 2)))
    rec <- list(iteration = emc[[1]]$samples$idx, stage = "burn", chain = cc,
                donor = from, gap = d$chains$gap[cc], z = d$chains$z[cc], exit = TRUE)
    emc[[1]]$culls <- c(emc[[1]]$culls, list(rec))
    if(verbose) message(sprintf("Burn exit: chain %d re-seeded from chain %d (log-likelihood %.0f below, z = %.1f)",
                                cc, from, rec$gap, rec$z))
  }
  emc
}

# Persistence bookkeeping across checks: consecutive flags per chain
update_cull_state <- function(state, diag, nch){
  if(is.null(state)) state <- list(chain_count = integer(nch))
  f <- diag$chains$flagged
  state$chain_count <- ifelse(f, state$chain_count + 1L, 0L)
  state
}

# After a cull in the sample stage: keep only the last draw, run `reburn`
# iterations under the frozen sample kernel, and keep only the last draw
# again, so the kept sample draws start after the respawned chain has settled.
cull_restart <- function(emc, reburn, verbose, verboseProgress, particle_factor, search_width,
                         cores_per_chain, cores_for_chains, n_blocks, on_singular, context){
  last_only <- function(emc) subset(emc, stage = "sample", filter = chain_n(emc)[1, "sample"] - 1, keep_stages = TRUE)
  emc <- last_only(emc)
  if(reburn > 0){
    emc <- add_proposals(emc, "sample", cores_per_chain * cores_for_chains, n_blocks)
    sub_emc <- run_native_chains(emc, stage = "sample", iter = reburn, verbose = verbose,
                                 verboseProgress = verboseProgress, particle_factor = particle_factor,
                                 search_width = search_width, n_cores = cores_per_chain,
                                 cores_for_chains = cores_for_chains, on_singular = on_singular, context = context)
    class(sub_emc) <- "emc"
    emc <- concat_emc(emc, sub_emc, reburn, "sample")
    emc <- last_only(emc)
  }
  class(emc) <- "emc"
  emc
}
