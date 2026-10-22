#' Diagnose a Fit: Chain Comparison, Drift and the Outcome Class
#'
#' Reads the stored draws of a fitted emc object (from [fit()] or [hardfit()])
#' and reports what the chains agree and disagree on, whether the stage is
#' still moving, and which of five outcome classes the fit is in. A
#' ``stop_record`` written by [hardfit()] is read when present but is not needed.
#'
#' The statistics, all over the stage's draws (the last two thirds for the
#' chain comparison, or the last ``window`` draws):
#'
#' * ``rule``: the stop rule's Rhat (the ``gd_quantile`` quantile of the
#'   subject-level Rhats and every group-level Rhat) and the plain maximum.
#' * ``stop``: min ESS and its parameter; ``rhat_b``, the Rhat of the group
#'   means projected on the direction along which the chains differ most, and
#'   ``rhat_ident``, on the direction the posterior pins down most tightly (the
#'   best-identified combination of parameters); the chain log-likelihood gap
#'   and its worst chain; the subject spread and its subject.
#' * ``chains``: per chain, the median summed subject log-likelihood, its gap to
#'   the best chain in units and in within-chain spreads, its own log posterior
#'   and whether that is below the other chains', and whether it would be culled.
#' * ``subjects``: per subject, the range over chains of the median
#'   log-likelihood and the chains at either end, widest first.
#' * ``flat``: the drift of the group means, the subject parameters, the summed
#'   log-likelihood and the subject log-likelihoods between the first and last
#'   thirds of the stage, in IQR units, with its Monte Carlo error.
#'
#' The classes, in the order they are tested:
#'
#' 1. **aborted**: ``stop_record$aborted`` is set. Do not continue; read the
#'    next two cases.
#' 2. **met**: ``stop_record$met``, or every criterion passes on the draws.
#' 3. **ess_held**: the Rhat rule passes but the effective size is below
#'    ``min_ess``. Continue the fit once or twice.
#' 4. **ridge_only**: the rule fails, but ``rhat_ident`` is within
#'    ``rhat_ident_max``, the likelihood gap and subject spread are within
#'    their thresholds and the log-likelihood drift is within ``drift_se``
#'    Monte Carlo errors: the chains agree on everything the data identify and
#'    disagree along a direction the data do not. Report the identified
#'    quantities; reparameterise or tighten the prior along the ridge.
#' 5. **region_mode**: otherwise. A chain is in a worse region
#'    (``ll_worst_chain``) or a subject has two modes (``subj_spread_subject``);
#'    compare seeds.
#'
#' ``falling`` is ``TRUE`` for an unmet fit whose rule Rhat decreased over its
#' last three checks (from ``stop_record$trace``): one continuation is worth it.
#'
#' ``prior``: each group mean's posterior mean (``type = "single"``: each
#' subject parameter's) in prior SDs from its prior mean, over the stage's draws
#' pooled over chains; a parameter more than ``prior_se`` SDs away is flagged.
#' Along a direction the data leave flat the prior is what sets the posterior,
#' and a likelihood leaning the other way makes two shallow plateaus the chains
#' cannot settle between. The flag is reported with the class, never used as a
#' stop criterion.
#'
#' @param emc An emc object
#' @param stage The stage to diagnose (default the last stage run)
#' @param window Integer. Draws at the end of the stage the chain comparison
#'   uses; ``NULL`` (default) uses the last two thirds of the stage
#' @param max_gd,gd_quantile,min_ess,max_ll_gap,max_subj_spread The stop rule's
#'   thresholds, as in [hardfit()] (``max_gd = NULL``: 1.1, or 1.05 for
#'   ``type = "single"``)
#' @param rhat_ident_max,drift_se The ridge-only class's thresholds (defaults
#'   1.1 and 3)
#' @param prior_se Numeric. A parameter whose posterior mean is more than this
#'   many prior SDs from its prior mean is flagged in ``prior`` (default 3)
#' @return A list of class ``emc.diagnose`` with ``class``, ``falling``,
#'   ``next_step``, ``rule``, ``stop``, ``chains``, ``subjects``, ``flat``,
#'   ``prior``, ``record`` and ``culls``.
#' @examples
#' d <- diagnose(samples_LNR)
#' d$class
#' @export
diagnose <- function(emc, stage = get_last_stage(emc), window = NULL, max_gd = NULL, gd_quantile = .99,
                     min_ess = 100, max_ll_gap = 100, max_subj_spread = 15, rhat_ident_max = 1.1, drift_se = 3,
                     prior_se = 3){
  emc <- restore_duplicates(emc)
  single <- emc[[1]]$type == "single"
  if(is.null(max_gd)) max_gd <- if(single) 1.05 else 1.1
  selection <- if(single) "alpha" else c("alpha", "mu")
  n <- unname(chain_n(emc)[1, stage])
  if(n == 0) stop("the object has no draws of stage ", stage)
  cull <- cull_defaults(TRUE)
  if(!is.null(window)) cull$window <- window else cull$min_window <- 20L   # the last two thirds, however short
  d <- cull_diagnostics(emc, stage, cull)
  st <- sample_stop_stats(emc, stage)
  gd <- check_gd(emc, stage, max_gd, NULL, TRUE, selection, gd_quantile)
  rule <- list(rule = gd$rule, max = if(length(gd$gd)) max(gd$gd) else NA_real_, threshold = max_gd, quantile = gd_quantile,
               met = isTRUE(gd$gd_done))
  flat <- tryCatch(flat_diagnostics(emc, stage), error = function(e) list())
  record <- emc[[1]]$stop_record
  val <- function(x) if(is.null(x) || !is.finite(x)) NA_real_ else x
  ess_ok <- is.null(min_ess) || isTRUE(val(st$ess_min) >= min_ess)
  gap_ok <- is.null(max_ll_gap) || isTRUE(val(st$ll_gap) <= max_ll_gap)
  spread_ok <- single || is.null(max_subj_spread) || isTRUE(val(st$subj_spread) <= max_subj_spread)
  ll_drift_ok <- is.null(flat$LL) || !is.finite(flat$LL$drift) || flat$LL$drift <= drift_se * flat$LL$mcse
  class <- if(!is.null(record$aborted)) "aborted"
    else if(isTRUE(record$met) || (rule$met && ess_ok && gap_ok && spread_ok)) "met"
    else if(rule$met && !ess_ok && gap_ok && spread_ok) "ess_held"
    else if(!rule$met && isTRUE(val(st$rhat_ident) <= rhat_ident_max) && gap_ok && spread_ok && ll_drift_ok) "ridge_only"
    else "region_mode"
  falling <- FALSE
  if(class %in% c("ridge_only", "region_mode", "ess_held") && !is.null(record$trace) && nrow(record$trace) >= 3){
    r <- utils::tail(record$trace$rule, 3); falling <- all(is.finite(r)) && all(diff(r) < 0)
  }
  next_step <- switch(class,
    met = "Converged under every criterion: use the fit.",
    aborted = "Stopped early: do not continue. Read rhat_ident, the log-likelihood gap and the subject spread to tell a ridge from a region/mode problem.",
    ess_held = "The chains agree but the effective sample size is short: continue the fit (hardfit(emc)) once or twice.",
    ridge_only = "The chains agree on every identified quantity and disagree along a direction the data do not fix. Report the identified quantities (rhat_ident near 1); reparameterise or tighten the prior on the weakly identified parameter. Continuing helps only while the rule Rhat is still falling.",
    region_mode = "The chains disagree on something the data identify: a chain in a worse region (ll_worst_chain) or a subject with two modes (subj_spread_subject). Fit several seeds and compare their summed log-likelihoods; hold a seed 10 or more units below the best.")
  if(falling) next_step <- paste(next_step, "The rule Rhat fell over the last three checks: one continuation is worth trying.")
  prior <- prior_tension(emc, stage, prior_se)
  if(!is.null(prior) && any(prior$flagged))
    next_step <- paste(next_step, sprintf("The posterior mean of %s is more than %s prior SDs from the prior mean: the prior and the likelihood pull against each other there; check the prior before continuing.",
                                          paste(prior$parameter[prior$flagged], collapse = ", "), format(prior_se)))
  out <- list(class = class, falling = falling, next_step = next_step, rule = rule, stop = st, chains = d$chains,
              subjects = d$subjects, flat = flat, prior = prior, record = record, culls = emc[[1]]$culls, stage = stage, draws = n,
              thresholds = list(min_ess = min_ess, max_ll_gap = max_ll_gap, max_subj_spread = max_subj_spread,
                                rhat_ident_max = rhat_ident_max, drift_se = drift_se, prior_se = prior_se))
  class(out) <- "emc.diagnose"
  out
}

# Posterior mean of each group mean (no group level: each subject parameter,
# named par:subject) in prior SDs from its prior mean, over the stage's draws
# pooled over chains. NULL when the prior has no mean and variance for them
# (the group means of a type whose prior is on something else).
prior_tension <- function(emc, stage, prior_se = 3){
  pr <- emc[[1]]$prior
  if(is.null(pr$theta_mu_mean) || is.null(pr$theta_mu_var)) return(NULL)
  draws <- lapply(emc, function(e){
    its <- which(e$samples$stage == stage)
    if(!is.null(e$samples$theta_mu)) return(e$samples$theta_mu[, its, drop = FALSE])
    a <- e$samples$alpha[, , its, drop = FALSE]
    nm <- if(dim(a)[2] > 1) as.vector(outer(dimnames(a)[[1]], dimnames(a)[[2]], paste, sep = ":")) else dimnames(a)[[1]]
    matrix(a, prod(dim(a)[1:2]), dimnames = list(nm, NULL))
  })
  post <- rowMeans(do.call(cbind, draws))
  pm <- pr$theta_mu_mean; pv <- pr$theta_mu_var
  psd <- sqrt(if(is.matrix(pv)) diag(pv) else pv)
  if(is.null(names(pm))) names(pm) <- names(post)[seq_along(pm)]
  names(psd) <- names(pm)
  nm <- names(post); key <- sub(":.*$", "", nm)                 # par:subject reads the parameter's prior
  ok <- key %in% names(pm)
  if(!any(ok)) return(NULL)
  z <- (post[ok] - pm[key[ok]]) / psd[key[ok]]
  data.frame(parameter = nm[ok], post_mean = unname(post[ok]), prior_mean = unname(pm[key[ok]]),
             prior_sd = unname(psd[key[ok]]), z = unname(z), flagged = is.finite(z) & abs(unname(z)) > prior_se,
             row.names = NULL)
}

#' @rdname diagnose
#' @param x An ``emc.diagnose`` object
#' @param ... Ignored
#' @export
print.emc.diagnose <- function(x, ...){
  f <- function(v, d = 3) if(is.null(v) || is.na(v)) "NA" else formatC(v, digits = d, format = "f")
  cat(sprintf("Class: %s%s (%s stage, %d draws per chain, %d chains)\n", x$class, if(x$falling) ", rule Rhat still falling" else "",
              x$stage, x$draws, nrow(x$chains)))
  cat(sprintf("Rule Rhat %s (threshold %s; plain max %s)", f(x$rule$rule), format(x$rule$threshold), f(x$rule$max)))
  s <- x$stop
  if(!is.null(s)){
    cat(sprintf(" | min ESS %s (%s) | rhat_b %s, rhat_ident %s\nLL gap %s (worst chain %s) | subject spread %s (%s; %s > 5)\n",
                f(s$ess_min, 0), s$ess_min_par, f(s$rhat_b), f(s$rhat_ident), f(s$ll_gap, 0), s$ll_worst_chain,
                f(s$subj_spread, 1), s$subj_spread_subject, s$subj_spread_n5))
  } else cat("\n(stage too short for the stop statistics: fewer than 20 draws)\n")
  if(length(x$flat)){
    cat("Drift (first vs last third, in IQRs; se):", paste(vapply(names(x$flat), function(n)
      sprintf("%s %.2f (%.2f)", n, x$flat[[n]]$drift, x$flat[[n]]$mcse), ""), collapse = ", "), "\n")
  }
  if(!all(is.na(x$chains$gap))){
    cat("\nChains (median summed log-likelihood; gap to the best chain, z in within-chain spreads; own log posterior lower than the others'; would be culled):\n")
    print(x$chains, row.names = FALSE, digits = 4)
  }
  if(nrow(x$subjects) && any(x$subjects$spread > 5)){
    cat(sprintf("\nSubjects with spread > 5 (%d):\n", sum(x$subjects$spread > 5)))
    print(utils::head(x$subjects[x$subjects$spread > 5, , drop = FALSE], 10), row.names = FALSE, digits = 3)
  }
  if(!is.null(x$prior) && any(x$prior$flagged)){
    p <- x$prior[x$prior$flagged, ]
    cat(sprintf("\nPrior tension (posterior mean in prior SDs from the prior mean, beyond %s): %s\n", format(x$thresholds$prior_se),
                paste(sprintf("%s %+.1f", p$parameter, p$z), collapse = ", ")))
  }
  if(!is.null(x$record)){
    r <- x$record
    cat(sprintf("\nRecord: met %s, tries %d, restarts %d, culls %d, discards %d%s\n", r$met, r$tries, r$restarts, r$culls, length(r$discards),
                if(is.null(r$aborted)) "" else sprintf(", aborted at check %d (rule %.3f > %s)", r$aborted$try, r$aborted$rule, format(r$aborted$abort_gd))))
  }
  if(!is.null(x$culls)){ cat("Culls:\n"); print(do.call(rbind, lapply(x$culls, as.data.frame)), row.names = FALSE, digits = 4) }
  cat("\n", x$next_step, "\n", sep = "")
  invisible(x)
}
