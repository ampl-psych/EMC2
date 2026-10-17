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
                   cell_gap = 5, cell_z = 3, max_culls = 2L, reburn = 100L)
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
#' subject-by-chain cells whose parameters another chain would improve on under
#' this chain's group parameters, and (iii) the flat test: the drift of the
#' group means and subject parameters between the first and last thirds of the
#' stage, in units of their interquartile range. These are the statistics the
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
#'   posterior, gap to the best chain, robust z, flagged), `cells` (the flagged
#'   subject-by-chain cells with the donor chain and the gain), `flat` (drift
#'   statistics per selection) and `culls` (the record of culls in the object).
#' @export
chain_diagnostics <- function(emc, stage = get_last_stage(emc), window = NULL, cull = TRUE,
                              flat_selection = c("mu", "alpha"), flat_p1 = 1/3, flat_p2 = 1/3,
                              flat_quantile = c(mu = 1, alpha = .99)){
  emc <- restore_duplicates(emc)
  cull <- cull_defaults(cull); if(is.null(cull)) cull <- cull_defaults(TRUE)
  if(!is.null(window)) cull$window <- window
  d <- cull_diagnostics(emc, stage, cull)
  f <- flat_diagnostics(emc, stage, flat_selection, flat_p1, flat_p2, flat_quantile)
  out <- list(chains = d$chains, cells = d$cells, flat = f, culls = emc[[1]]$culls)
  class(out) <- "emc.diagnostics"
  out
}

#' @export
print.emc.diagnostics <- function(x, ...){
  if(all(is.na(x$chains$gap))) cat("Chains: the stage is too short for the mode comparison (fewer than 500 iterations).\n") else
  cat("Chains (log posterior under the common group model over the window; gap to the best chain, z in within-chain spread):\n")
  print(x$chains, row.names = FALSE, digits = 3)
  if(nrow(x$cells)) { cat("\nSubject-by-chain cells another chain would improve on:\n"); print(x$cells, row.names = FALSE, digits = 3) }
  else cat("\nNo subject-by-chain cell flagged.\n")
  if(length(x$flat)) { cat("\nFlat test (drift between the first and last thirds, in IQRs; mcse in the same units):\n")
    print(do.call(rbind, lapply(names(x$flat), function(n) data.frame(selection = n, drift = x$flat[[n]]$drift, mcse = x$flat[[n]]$mcse, worst = x$flat[[n]]$worst))), row.names = FALSE, digits = 3) }
  if(!is.null(x$culls)) { cat("\nCulls applied:\n"); print(do.call(rbind, lapply(x$culls, as.data.frame)), row.names = FALSE) }
  invisible(x)
}

# Chain- and cell-level statistics over the window. `cells` holds every
# (subject, chain) whose parameters a donor chain improves on by more than
# cell_gap log posterior units and cell_z within-chain spreads, judged under
# the recipient chain's own group parameters.
cull_diagnostics <- function(emc, stage, cull){
  st <- chain_window_stats(emc, stage, cull$window, cull$min_window)
  nch <- length(emc)
  empty <- list(chains = data.frame(chain = seq_len(nch), median_lp = NA, gap = NA, z = NA, flagged = FALSE),
                cells = data.frame(subject = character(0), chain = integer(0), donor = integer(0), gain = numeric(0), z = numeric(0)))
  if(is.null(st) || length(st[[1]]$it) < 20 || nch < 2) return(empty)
  group <- uses_group_density(emc)
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
  # cells: candidate subjects are those whose per-subject log posterior differs between chains
  smed <- sapply(st, function(s) apply(s$lp, 1, median)); smad <- sapply(st, function(s) apply(s$lp, 1, mad))
  if(is.null(dim(smed))) { smed <- matrix(smed, 1); smad <- matrix(smad, 1) }
  cg <- robust_gaps(smed, smad)
  cand <- which(cg$gap > cull$cell_gap & cg$z > cull$cell_z, arr.ind = TRUE)
  cells <- empty$cells
  subjects <- rownames(emc[[1]]$samples$subj_ll)
  if(nrow(cand)){
    rows <- list()
    for(r in seq_len(nrow(cand))){
      s <- cand[r, 1]; cc <- cand[r, 2]; donor <- cg$best[s]
      if(donor == cc) next
      # the donor's median parameter vector for this subject, judged in the recipient's group
      rec <- emc[[cc]]$samples; dn <- st[[donor]]
      a_donor <- apply(emc[[donor]]$samples$alpha[, s, dn$it, drop = FALSE], 1, median)
      ll_donor <- median(dn$ll[s, ])
      if(group){
        gd <- vapply(st[[cc]]$it, function(t) group_logdens(matrix(a_donor), rec$theta_mu[, t], rec$theta_var[, , t]), 0)
        value_in_recipient <- ll_donor + median(gd)
      } else value_in_recipient <- ll_donor
      own <- st[[cc]]$lp[s, ]
      gain <- value_in_recipient - median(own); z <- gain / max(mad(own), 1e-8)
      if(is.finite(gain) && gain > cull$cell_gap && z > cull$cell_z)
        rows[[length(rows) + 1]] <- data.frame(subject = subjects[s], chain = cc, donor = donor, gain = gain, z = z)
    }
    if(length(rows)) cells <- do.call(rbind, rows)
  }
  # A chain whose gap is mostly explained by its flagged cells gets the cell
  # repairs, not a cull of the whole chain; otherwise its cells are not listed.
  for(cc in chains$chain[chains$flagged]){
    explained <- sum(cells$gain[cells$chain == cc])
    if(explained >= .5 * chains$gap[cc]) chains$flagged[cc] <- FALSE
  }
  cells <- cells[!cells$chain %in% chains$chain[chains$flagged], , drop = FALSE]
  list(chains = chains, cells = cells)
}

# ---- flat test ---------------------------------------------------------------

# Drift of one series: |median(first third) - median(last third)| over the
# median of the two thirds' IQRs, and the Monte Carlo error of that
# difference in the same units
flat_series <- function(x, p1, p2){
  n <- length(x); if(n < 10) return(c(drift = NA, mcse = NA))
  a <- x[seq_len(max(1, round(n * p1)))]; b <- x[(n - max(1, round(n * p2)) + 1):n]
  iqr <- max(median(c(stats::IQR(a), stats::IQR(b))), 1e-8)
  ess <- function(v) tryCatch(as.numeric(coda::effectiveSize(v)), error = function(e) length(v))
  mcse <- sqrt((1.25 * sd(a))^2 / max(ess(a), 1) + (1.25 * sd(b))^2 / max(ess(b), 1)) / iqr
  c(drift = abs(median(a) - median(b)) / iqr, mcse = mcse)
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
  } else stop("flat_selection must be one or more of: mu, alpha, LL")
  m
}

flat_diagnostics <- function(emc, stage, selection, p1, p2, quantile){
  out <- list()
  for(sel in selection){
    m <- pooled_matrix(emc, sel, stage); if(is.null(m) || nrow(m) < 10) next
    if(sel == "alpha"){ # medians and IQRs vectorised; the Monte Carlo error from a sample of columns
      n <- nrow(m); a <- m[seq_len(max(1, round(n * p1))), , drop = FALSE]; b <- m[(n - max(1, round(n * p2)) + 1):n, , drop = FALSE]
      qa <- apply(a, 2, quantile, c(.25, .5, .75)); qb <- apply(b, 2, quantile, c(.25, .5, .75))
      iqr <- pmax(apply(rbind(qa[3, ] - qa[1, ], qb[3, ] - qb[1, ]), 2, median), 1e-8)
      drift <- abs(qa[2, ] - qb[2, ]) / iqr
      cols <- unique(round(seq(1, ncol(m), length.out = min(50, ncol(m)))))
      mcse <- median(apply(m[, cols, drop = FALSE], 2, function(x) flat_series(x, p1, p2)[["mcse"]]), na.rm = TRUE)
    } else {
      fs <- apply(m, 2, flat_series, p1 = p1, p2 = p2); drift <- fs["drift", ]; mcse <- max(fs["mcse", ], na.rm = TRUE)
    }
    q <- if(!is.null(quantile[sel])) quantile[[sel]] else 1
    stat <- if(q >= 1) max(drift, na.rm = TRUE) else unname(stats::quantile(drift, q, na.rm = TRUE))
    out[[sel]] <- list(drift = stat, mcse = mcse, worst = colnames(m)[which.max(drift)], quantile = q, all = drift)
  }
  out
}

# TRUE when every selection's drift is within max_flat_loc, or within
# flat_mcse Monte Carlo errors of zero (so the test cannot fail on noise alone)
flat_done <- function(flat, max_flat_loc, flat_mcse){
  if(!length(flat)) return(TRUE)
  all(vapply(flat, function(f) !is.finite(f$drift) || f$drift <= max_flat_loc || f$drift <= flat_mcse * f$mcse, logical(1)))
}

# ---- culling: respawn ----------------------------------------------------------

# Overwrite the last stored draw of chain `to` with the state of chain `from`
# at a draw of its window (a random one at least `back` draws from the end
# where the window allows), so the next run starts the chain there.
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
  emc[[to]]$samples <- s_to
  emc
}

# Replace subject `s` of chain `to` by the donor's parameters at a random draw of the donor's window
respawn_subject <- function(emc, to, from, s, its_from, back = 100){
  pick <- its_from[seq_len(max(1, length(its_from) - back))]
  t <- pick[sample.int(length(pick), 1)]
  last <- emc[[to]]$samples$idx
  emc[[to]]$samples$alpha[, s, last] <- emc[[from]]$samples$alpha[, s, t]
  emc[[to]]$samples$subj_ll[s, last] <- emc[[from]]$samples$subj_ll[s, t]
  emc
}

# Apply the culls that have persisted: whole chains first, then cells of the
# surviving chains. Returns the object and a record of what was done.
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
  cells <- diag$cells
  if(nrow(cells)){
    key <- paste(cells$subject, cells$chain)
    due <- state$cell_count[key] >= cull$persistence; due[is.na(due)] <- FALSE
    for(r in which(due)){
      if(cells$chain[r] %in% culled_chains) next
      s <- match(cells$subject[r], rownames(emc[[1]]$samples$subj_ll))
      emc <- respawn_subject(emc, cells$chain[r], cells$donor[r], s, its[[cells$donor[r]]], cull$reburn)
      done[[length(done) + 1]] <- list(iteration = emc[[1]]$samples$idx, stage = stage, chain = cells$chain[r], subject = cells$subject[r],
                                       donor = cells$donor[r], gap = cells$gain[r], z = cells$z[r])
    }
  }
  list(emc = emc, done = done)
}

# Persistence bookkeeping across checks: consecutive flags per chain and per cell
update_cull_state <- function(state, diag, nch){
  if(is.null(state)) state <- list(chain_count = integer(nch), cell_count = integer(0))
  f <- diag$chains$flagged
  state$chain_count <- ifelse(f, state$chain_count + 1L, 0L)
  key <- if(nrow(diag$cells)) paste(diag$cells$subject, diag$cells$chain) else character(0)
  new <- setNames(integer(length(key)), key)
  for(k in key) new[k] <- (if(!is.na(state$cell_count[k])) state$cell_count[k] else 0L) + 1L
  state$cell_count <- new
  state
}
