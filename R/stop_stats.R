# Per-check statistics for the sample stage's stop rule, read on the draws the
# Rhat read (after any first-third discard). Each can be made a criterion
# (opt-in stop_criteria entries; see ?fit) and all are reported at every check:
#
# * min_ess: the smallest effective sample size over the group means and the
#   log group SDs, each summed over chains. Along a ridge the group means mix
#   slowly enough that the Rhat rule can be met with an effective size near 100.
# * ll_gap: the gap between the best and worst chain in the median, over the
#   last two thirds of the draws, of the summed subject log-likelihood. The log
#   posterior is not used here: its group-density term varies between chains by
#   hundreds of units in fits that are fine by every other measure (it follows
#   each chain's group covariance draws), while the summed log-likelihood varies
#   by tens and a chain in a different region is thousands below.
# * subj_spread: per subject, the range over chains of the median
#   log-likelihood over the same window; its maximum, the subject, and the
#   number of subjects above 5 units. A subject in a different mode in some
#   chains moves the summed statistic by too little to see.
# * rhat_b, rhat_ident (reported only): Rhat of the group means projected on the
#   direction of largest between-chain disagreement (the leading eigenvector of
#   the covariance of the chain means, in within-chain SD units) and on the
#   pooled within-chain axis of least variance (the best-identified combination,
#   e.g. the match rate on a v / v_lMd ridge). The first reproduces the largest
#   per-parameter Rhat; the second tells whether a fit stopped at the limit is
#   usable off the ridge.
sample_stop_stats <- function(emc, stage = "sample"){
  its <- lapply(emc, function(x) which(x$samples$stage[seq_len(x$samples$idx)] == stage))
  n <- min(lengths(its)); nch <- length(emc)
  if(n < 20 || nch < 2) return(NULL)
  its <- lapply(its, function(it) utils::tail(it, n))
  out <- list(draws = n)
  s1 <- emc[[1]]$samples
  if(!is.null(s1$theta_mu) && !is.null(s1$theta_var)){
    MU <- lapply(seq_len(nch), function(i) emc[[i]]$samples$theta_mu[, its[[i]], drop = FALSE])
    pn <- rownames(MU[[1]])
    LSD <- lapply(seq_len(nch), function(i){
      V <- emc[[i]]$samples$theta_var[, , its[[i]], drop = FALSE]
      m <- matrix(apply(V, 3, function(v) 0.5 * log(diag(v))), nrow = length(pn), dimnames = list(pn, NULL))
      m
    })
    ess <- c(stats::setNames(vapply(pn, function(p) sum(vapply(MU, function(m) ess_series(m[p, ]), 0)), 0), pn),
             stats::setNames(vapply(pn, function(p) sum(vapply(LSD, function(m) ess_series(m[p, ]), 0)), 0), paste0("lsd:", pn)))
    out$ess_min <- min(ess, na.rm = TRUE); out$ess_min_par <- names(ess)[which.min(ess)]
    out <- c(out, projected_rhats(MU))
  }
  k <- seq(floor(n / 3) + 1, n)
  ll <- lapply(seq_len(nch), function(i) emc[[i]]$samples$subj_ll[, its[[i]][k], drop = FALSE])
  tot <- vapply(ll, colSums, numeric(length(k)))
  if(is.null(dim(tot))) tot <- matrix(tot, nrow = length(k))
  med <- apply(tot, 2, stats::median); md <- apply(tot, 2, stats::mad)
  out$ll_gap <- max(med) - min(med); out$ll_z <- out$ll_gap / max(stats::median(md), 1e-8)
  out$ll_worst_chain <- which.min(med)
  sm <- vapply(ll, function(x) apply(x, 1, stats::median), numeric(nrow(ll[[1]])))
  if(is.null(dim(sm))) sm <- matrix(sm, nrow = 1)
  sp <- apply(sm, 1, function(v) diff(range(v)))
  out$subj_spread <- max(sp); out$subj_spread_n5 <- sum(sp > 5)
  out$subj_spread_subject <- rownames(s1$subj_ll)[which.max(sp)]
  out
}

# Effective size of one chain's series (coda's spectral estimate); 0 for a constant series
ess_series <- function(x){
  if(length(x) < 4 || stats::sd(x) == 0) return(0)
  unname(coda::effectiveSize(x))
}

# Rhat of the group means projected on the direction of largest between-chain
# disagreement and on the pooled within-chain axis of least variance, both in
# units of the pooled within-chain SD. MU: one parameters x draws matrix per chain.
projected_rhats <- function(MU){
  nch <- length(MU); p <- nrow(MU[[1]]); n <- ncol(MU[[1]])
  if(p < 2) return(list(rhat_b = NA_real_, rhat_ident = NA_real_))
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

# The progress line's part for these statistics
stop_stats_message <- function(st){
  if(is.null(st)) return(NULL)
  paste0(if(!is.null(st$ess_min)) sprintf(" | min ESS=%.0f (%s)", st$ess_min, st$ess_min_par),
         sprintf(" | LL gap=%.0f | subject spread=%.1f (%s; %d > 5)", st$ll_gap, st$subj_spread,
                 st$subj_spread_subject, st$subj_spread_n5))
}
