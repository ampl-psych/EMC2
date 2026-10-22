# Chain culling for hardfit() and diagnose(): a chain held in a different
# region from the others is respawned from a healthy one. Chains are compared
# on the summed subject log-likelihood (chain_table()), not the log posterior
# under pooled group parameters: its group-density term differs between healthy
# chains by hundreds of units (it follows each chain's group covariance draws),
# so culling on it removed wandering chains, not stuck ones. Below a gap of
# lp_guard_below the chain must also be lower on its own log posterior (the
# likelihood plus the group density under its own group parameters), because
# the likelihood alone can favour a state with an implausible group covariance;
# above it the likelihood decides, since a stuck chain thousands of units down
# is far lower on both and the own-covariance term is noisy at that scale.

cull_defaults <- function(cull){
  defaults <- list(persistence = 2L, window = NULL, min_window = 500L, chain_gap = 200, chain_z = 5,
                   lp_guard_below = 1000, max_culls = 1L, reburn = 100L)
  if(isTRUE(cull)) return(defaults)
  if(is.null(cull) || isFALSE(cull)) return(NULL)
  if(!is.list(cull)) stop("`cull` must be TRUE, FALSE or a list")
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

# The chain table with the cull decision: a chain is flagged when it is
# chain_gap units and chain_z within-chain spreads below the best chain (a
# different region: fine chains differ by tens of units, a stuck one by
# hundreds or thousands) and, below lp_guard_below, also lower on its own log
# posterior. The subject spread is reported, not acted on: moving a subject
# between modes needs a proper Metropolis move, not a copy of another chain's
# parameters.
cull_diagnostics <- function(emc, stage, cull,
                             tab = chain_table(chain_window_stats(emc, stage_window(emc, stage, cull$window, cull$min_window)),
                                               rownames(emc[[1]]$samples$subj_ll))){
  nch <- length(emc)
  if(is.null(tab))
    return(list(chains = data.frame(chain = seq_len(nch), median_ll = NA, gap = NA, z = NA, median_lp = NA, lp_lower = NA, flagged = FALSE),
                subjects = data.frame(subject = character(0), spread = numeric(0), best = integer(0), worst = integer(0))))
  ch <- tab$chains
  ch$flagged <- ch$gap > cull$chain_gap & ch$z > cull$chain_z & (ch$lp_lower | ch$gap >= cull$lp_guard_below)
  ch$flagged[is.na(ch$flagged)] <- FALSE
  # never cull all chains but one on a single check: at least two survive
  if(sum(ch$flagged) > nch - 2) ch$flagged[order(ch$gap)[1:2]] <- FALSE
  list(chains = ch, subjects = tab$subjects)
}

# One cull check at a burn or sample check: the comparison, the persistence
# bookkeeping and any respawn. `state` carries the consecutive flags per chain
# and the culls made in the stage; nothing more is culled once max_culls is
# reached. `tab` is the chain table when the caller has it already.
cull_step <- function(emc, stage, cull, state, tab = NULL){
  if(is.null(state)) state <- list(chain_count = integer(length(emc)), n_culls = 0L)
  if(state$n_culls >= cull$max_culls) return(list(emc = emc, state = state, done = list(), note = NULL))
  its <- stage_window(emc, stage, cull$window, cull$min_window)
  if(is.null(tab)) tab <- chain_table(chain_window_stats(emc, its), rownames(emc[[1]]$samples$subj_ll))
  diag <- cull_diagnostics(emc, stage, cull, tab)
  state$chain_count <- ifelse(diag$chains$flagged, state$chain_count + 1L, 0L)
  culled <- diag$chains$chain[diag$chains$flagged & state$chain_count >= cull$persistence]
  r <- respawn_flagged(emc, culled, diag, stage, its, cull$reburn, exit = FALSE)
  note <- NULL
  if(length(r$done)){
    state$n_culls <- state$n_culls + 1L; state$chain_count[] <- 0L
    note <- paste(vapply(r$done, function(d) sprintf("chain %d respawned from chain %d (log-likelihood %.0f below, z = %.1f)",
                                                     d$chain, d$donor, d$gap, d$z), ""), collapse = "; ")
  }
  list(emc = r$emc, state = state, done = r$done, note = note)
}

# At the end of burn (by its rule or at max_tries), re-seed every chain held
# apart, with no persistence, over the last two thirds of burn however short:
# burn draws are discarded, so a wrong re-seed is free. Only the chain's last
# draw is replaced; adapt's first proposals (create_chain_proposals(), the last
# min(250, idx / 1.5) draws of every chain) still pool 1/n_chains of the stuck
# history, and slide past it within adapt's first 250 iterations, which
# adapt_converge's minimum guarantees before the sample kernel is built.
burn_exit_reseed <- function(emc, cull, verbose = FALSE){
  cx <- cull; cx$window <- NULL; cx$min_window <- min_stat_draws
  d <- cull_diagnostics(emc, "burn", cx)
  bad <- d$chains$chain[d$chains$flagged]
  if(!length(bad)) return(emc)
  its <- stage_window(emc, "burn", NULL, cx$min_window)
  back <- min(cull$reburn, floor(length(its[[1]]) / 2))   # half the burn window at most: burn may be short at exit
  r <- respawn_flagged(emc, bad, d, "burn", its, back, exit = TRUE)
  if(verbose) for(x in r$done)
    message(sprintf("Burn exit: chain %d re-seeded from chain %d (log-likelihood %.0f below, z = %.1f)", x$chain, x$donor, x$gap, x$z))
  r$emc
}

# Respawn each of `chains` from a random chain that is not flagged, and record
# it on emc[[1]]$culls. Returns the object and this call's records.
respawn_flagged <- function(emc, chains, diag, stage, its, back, exit){
  donors <- setdiff(seq_along(emc), diag$chains$chain[diag$chains$flagged])
  done <- list()
  for(cc in chains){
    from <- donors[sample.int(length(donors), 1)]
    emc <- respawn_chain(emc, cc, from, its[[from]], back)
    done[[length(done) + 1]] <- list(iteration = emc[[1]]$samples$idx, stage = stage, chain = cc, donor = from,
                                     gap = diag$chains$gap[cc], z = diag$chains$z[cc], exit = exit)
  }
  emc[[1]]$culls <- c(emc[[1]]$culls, done)
  list(emc = emc, done = done)
}

# Overwrite the last stored draw of chain `to` with the state of chain `from`
# at a random draw of its window at least `back` draws from the end (where the
# window allows), so the next run starts the chain there and, after `back`
# draws, the two chains are as far apart as the donor is from its own past.
# The donor's kernel settings come too (the joint moves' step and location
# scales, the particle settings): settings adapted where the culled chain was
# stuck do not fit the region it restarts in. Its random number state stays.
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

# After a cull in the sample stage: keep only the last draw, run `reburn` more
# under the frozen sample kernel with `run(emc, n)` (hardfit()'s block), and
# keep only the last draw again, so the kept sample draws all start after the
# respawned chain has settled.
cull_restart <- function(emc, reburn, run){
  last_only <- function(emc) subset(emc, stage = "sample", filter = chain_n(emc)[1, "sample"] - 1, keep_stages = TRUE)
  emc <- last_only(emc)
  if(reburn > 0) emc <- last_only(run(emc, reburn))
  emc
}
