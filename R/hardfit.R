#' Fit a Model with the Hard-Mode Convergence Rules
#'
#' Runs the same stages as [fit()] (preburn, burn, adapt, sample) but drives
#' the burn and sample stages one check at a time and applies, between
#' checks, the rules that hard fits need: a stop rule with an effective-sample
#' floor and two likelihood-based vetoes, replacement of chains held in a
#' different region, an early stop for runs that will not converge, and a
#' record of how the fit ended. Everything is computed on the stored draws; the
#' sampler itself is [fit()]'s.
#'
#' At every check of the sample stage (``step_size`` draws per chain) the fit
#' is **met** when, on one set of kept draws, all of these hold:
#'
#' * **Rhat**: the stop rule's Rhat is below ``max_gd``. For hierarchical fits
#'   the rule reads the ``gd_quantile`` quantile (default .99) of the
#'   subject-level Rhats together with every group-level Rhat, so that with
#'   hundreds of subjects one subject's worst parameter does not set the rule;
#'   for ``type = "single"`` it is the maximum over the subject's parameters
#'   and the default threshold is 1.05.
#' * **Effective size**: the smallest effective sample size over the group
#'   means and the log group SDs (``type = "single"``: over the subject
#'   parameters), summed over chains, is at least ``min_ess``.
#' * **No chain in a different region**: the chains' summed subject
#'   log-likelihoods (medians over the last two thirds of the draws) lie within
#'   ``max_ll_gap`` units of each other.
#' * **No subject in a different region**: for every subject, the chains'
#'   median log-likelihoods lie within ``max_subj_spread`` units (not for
#'   ``type = "single"``).
#' * the stage holds at least ``iter`` new draws per chain.
#'
#' If the rule fails on all the stage's draws but passes with the first third
#' dropped, that shorter history is used only if every criterion then passes;
#' the drop is recorded in ``stop_record$discards``. Otherwise every draw is kept.
#'
#' **Culling** (``cull``): at each burn and sample check, a chain whose summed
#' subject log-likelihood is ``chain_gap`` (200) units and ``chain_z`` (5)
#' within-chain spreads below the best chain's for ``persistence`` (2) checks
#' running, and either ``lp_guard_below`` (1,000) units below or also lower on
#' its own log posterior, has its state replaced by a copy of a healthy chain's
#' (kernel settings included). Burn may do this twice at its checks and once
#' more at its end (``exit = TRUE`` in the record, no persistence needed);
#' sample does it once and then restarts its draws after ``reburn`` (100)
#' iterations, with the check count. Replacing a chain's state during warm-up
#' changes a starting point, not the posterior.
#'
#' **Early stop** (``abort_gd``, ``abort_try``): at check ``abort_try`` of the
#' sample stage, a fit whose rule Rhat is still above ``abort_gd`` ends there
#' with a warning and ``stop_record$aborted``. In replays of 46 hierarchical and
#' single-subject fits above 1.5 at check 20, none met the rule later.
#'
#' **The record** ``emc[[1]]$stop_record``: ``criteria`` (the thresholds that
#' were on), ``met``, ``tries``, ``rule`` (the Rhat the rule read last), ``stats``
#' (min ESS and its parameter, ``rhat_b`` and ``rhat_ident`` — Rhat along the
#' direction the chains disagree on most and along the best-identified
#' direction —, the log-likelihood gap and its worst chain, the subject spread
#' and its subject, the drift statistics ``flat``, and ``prior``, each group
#' mean's posterior mean in prior SDs from its prior mean), ``culls`` (count; the list
#' is ``emc[[1]]$culls``), ``discards``, ``aborted``, ``restarts``, ``trace``
#' (one row per check: draws, rule, plain max Rhat, min ESS, gap, spread) and
#' ``previous`` (the record of an earlier ``hardfit()`` call on the same
#' object). [diagnose()] reads it, and the draws, into one of the outcome
#' classes.
#'
#' Chains: the replacement of a bad chain and the detection of a subject with
#' two posterior modes need healthy chains to compare with; with 3 chains a
#' two-mode subject goes undetected in a quarter to a third of runs, with 6 in
#' about one in ten. A note is printed when a hierarchical object has fewer
#' than 6 chains (``make_emc(n_chains = 6)``).
#'
#' @inheritParams fit
#' @param emc An emc object made by [make_emc()], or the name of a file
#'   holding one.
#' @param stage The stage to start at; by default the last stage present in
#'   the object (so a call on a fitted object continues its sample stage).
#' @param iter Integer. The sample stage needs at least this many new draws
#'   per chain (default 1000).
#' @param step_size Integer. Draws per chain between checks (default 100).
#' @param max_tries Integer. Checks per stage before an unmet stage ends
#'   (default 40; ``fit()``'s default is 20).
#' @param max_gd Numeric. The rule's Rhat threshold; ``NULL`` (default) is 1.1,
#'   or 1.05 for ``type = "single"``.
#' @param gd_quantile Numeric in (0, 1]. The quantile of the subject-level
#'   Rhats the rule reads (default .99; 1 is the maximum).
#' @param min_ess,max_ll_gap,max_subj_spread Numeric thresholds of the three
#'   criteria above (defaults 100, 100, 15); ``NULL`` turns one off.
#' @param cull ``TRUE`` (default), ``FALSE``, or a list of settings
#'   (``persistence``, ``window``, ``min_window``, ``chain_gap``, ``chain_z``,
#'   ``lp_guard_below``, ``max_culls``, ``reburn``).
#' @param abort_gd,abort_try The early stop (defaults 1.5 at check 20);
#'   ``abort_gd = NULL`` turns it off.
#' @param burn A list with the burn stage's own ``mean_gd`` (default 1.1, the
#'   mean Rhat rule of [fit()]) and ``max_tries`` (default ``max_tries``).
#' @param history_proposals Logical. As in [fit()]; changing it on a fitted
#'   object re-runs adapt.
#' @param ... Passed to [run_emc()]: ``n_blocks``, ``trim``, ``verboseProgress``.
#' @return The fitted emc object, with ``emc[[1]]$stop_record`` and, when chains
#'   were replaced, ``emc[[1]]$culls``.
#' @examples
#' \dontrun{
#' emc <- make_emc(forstmann, design_DDMaE, n_chains = 6)
#' emc <- hardfit(emc)
#' diagnose(emc)
#' }
#' @export
hardfit <- function(emc, stage = NULL, iter = 1000, step_size = 100, max_tries = 40,
                    max_gd = NULL, gd_quantile = .99, min_ess = 100, max_ll_gap = 100,
                    max_subj_spread = 15, cull = TRUE, abort_gd = 1.5, abort_try = 20,
                    burn = list(), verbose = TRUE, fileName = NULL, particle_factor = 50,
                    search_width = 1, cores_per_chain = 1, cores_for_chains = length(emc),
                    on_singular = NULL, history_proposals = FALSE, ...){
  dots <- add_defaults(list(...), n_blocks = 1, verboseProgress = FALSE, trim = TRUE)
  if(is.character(emc)){
    emc <- fix_fileName(emc)
    if(is.null(fileName)) fileName <- emc
    emc <- loadRData(emc)
  }
  if(length(emc) < 2) stop("hardfit() needs at least 2 chains; 6 are recommended for hierarchical fits")
  single <- emc[[1]]$type == "single"
  if(is.null(max_gd)) max_gd <- if(single) 1.05 else 1.1
  check_number(max_gd, "max_gd", 1); check_number(gd_quantile, "gd_quantile", 0, 1)
  if(gd_quantile <= 0) stop("gd_quantile must be in (0, 1]")
  for(k in c("min_ess", "max_ll_gap", "max_subj_spread")){
    v <- get(k); if(!is.null(v)) check_number(v, k, 0)
  }
  if(!is.null(abort_gd)){
    check_number(abort_gd, "abort_gd", 1)
    if(abort_gd <= 1) stop("abort_gd must be above 1")
    check_number(abort_try, "abort_try", 1)
    if(abort_try != round(abort_try)) stop("abort_try must be a positive integer")
  }
  check_number(iter, "iter", 0); check_number(step_size, "step_size", 1); check_number(max_tries, "max_tries", 1)
  burn <- add_defaults(burn, mean_gd = 1.1, max_tries = max_tries)
  if(!single && length(emc) < 6 && verbose)
    message(sprintf("Note: %d chains. A subject with two posterior modes goes undetected in about a quarter to a third of 3-chain runs and about one in ten 6-chain runs; make_emc(n_chains = 6) is recommended for hard fits.", length(emc)))
  hp <- set_history_proposals(emc, if(missing(history_proposals)) NULL else history_proposals)
  emc <- hp$emc
  if(is.null(stage)) stage <- get_last_stage(emc)
  if(!is.null(hp$stage)) stage <- hp$stage
  stages <- c("preburn", "burn", "adapt", "sample")
  if(!stage %in% stages) stop("stage must be one of ", paste(stages, collapse = ", "))
  to_run <- stages[which(stages == stage):length(stages)]
  selection <- if(single) "alpha" else c("alpha", "mu")
  type <- emc[[1]]$type
  context <- sampler_contexts(restore_duplicates(emc), cores_per_chain)
  # one block of `n` more draws of `stage`; run_emc() computes no Rhat with these criteria. Burn and sample
  # blocks are silent (hardfit() prints its own line per check); preburn and adapt run whole and print fit()'s lines.
  block <- function(emc, stage, n, criteria = list(iter = n), tries = 1L, say_progress = FALSE){
    run_emc(emc, stage, stop_criteria = criteria, search_width = search_width, step_size = n,
            verbose = say_progress && verbose, verboseProgress = dots$verboseProgress, fileName = NULL,
            particle_factor = particle_factor, cores_per_chain = cores_per_chain,
            cores_for_chains = cores_for_chains, max_tries = tries, n_blocks = dots$n_blocks,
            trim = dots$trim, on_singular = on_singular, .context = context)
  }
  save_fit <- function(emc){
    if(is.null(fileName)) return(invisible(NULL))
    emc <- strip_duplicates(emc); class(emc) <- "emc"
    save(emc, file = fix_fileName(fileName))
  }
  say <- function(...) if(verbose) message(sprintf(...))
  eta <- function(t_check, tries_left) format_duration(t_check * max(0, tries_left))
  t_start <- Sys.time()

  if("preburn" %in% to_run){
    emc <- block(emc, "preburn", step_size, get_stop_criteria("preburn", NULL, type), max_tries, say_progress = TRUE)
    save_fit(emc)
  }

  if("burn" %in% to_run){
    say("Running burn stage")
    crit <- get_stop_criteria("burn", list(mean_gd = burn$mean_gd, selection = selection, omit_mpsrf = TRUE), type)
    cull_b <- cull_settings(cull, "burn"); state <- NULL; n_culls <- 0L; met <- FALSE
    for(try in seq_len(burn$max_tries)){
      t0 <- Sys.time()
      emc <- block(emc, "burn", step_size)
      gd <- check_gd(emc, "burn", NULL, crit$mean_gd, TRUE, selection)
      dropped <- 0L
      if(!isTRUE(gd$gd_done) && !is.null(gd$short)){
        dropped <- unname(chain_n(emc)[1, "burn"] - chain_n(gd$short$emc)[1, "burn"])
        emc <- gd$short$emc; gd <- gd$short
      }
      met <- isTRUE(gd$gd_done)
      cull_note <- NULL
      if(!met && !is.null(cull_b) && n_culls < cull_b$max_culls){
        act <- cull_step(emc, "burn", cull_b, state); emc <- act$emc; state <- act$state
        if(length(act$done)){ n_culls <- n_culls + 1L; state <- NULL; cull_note <- act$note }
      }
      save_fit(emc)
      elapsed <- Sys.time() - t0
      say("[burn | try=%d | iters=%d | Mean Rhat=%.3f%s%s] Duration: %s - ETA (to max_tries): %s", try, chain_n(emc)[1, "burn"],
          mean(gd$gd), if(dropped > 0) sprintf(" (first %d draws dropped)", dropped) else "",
          if(is.null(cull_note)) "" else paste0(" | ", cull_note), format_duration(elapsed), eta(elapsed, burn$max_tries - try))
      if(met) break
    }
    if(!met) warning("Burn stage ended at max_tries with its Rhat rule unmet.", call. = FALSE)
    if(!is.null(cull_b)) emc <- burn_exit_reseed(emc, cull_b, verbose)
    save_fit(emc)
  }

  if("adapt" %in% to_run){
    emc <- block(emc, "adapt", step_size, get_stop_criteria("adapt", NULL, type), max_tries, say_progress = TRUE)
    save_fit(emc)
  }

  if("sample" %in% to_run){
    say("Running sample stage")
    previous <- emc[[1]]$stop_record
    cull_s <- cull_settings(cull, "sample"); state <- NULL; n_culls <- 0L
    target <- chain_n(emc)[1, "sample"] + iter
    discards <- list(); aborted <- NULL; restarts <- 0L; trace <- list(); met <- FALSE; st <- NULL; rule <- NA_real_
    try <- 0L
    while(try < max_tries){
      try <- try + 1L; t0 <- Sys.time()
      emc <- block(emc, "sample", step_size)
      r <- sample_check(emc, target, max_gd, gd_quantile, min_ess, max_ll_gap, if(single) NULL else max_subj_spread, selection)
      dropped <- r$dropped
      if(dropped > 0) discards <- c(discards, list(list(try = try, removed = dropped, kept = r$draws)))
      emc <- r$emc; st <- r$st; rule <- r$rule; met <- r$met
      trace[[length(trace) + 1]] <- data.frame(try = try, draws = r$draws, rule = rule, max_rhat = max(r$gd),
        ess_min = if(is.null(st$ess_min)) NA_real_ else st$ess_min, ll_gap = if(is.null(st$ll_gap)) NA_real_ else st$ll_gap,
        subj_spread = if(is.null(st$subj_spread)) NA_real_ else st$subj_spread, culls = n_culls, dropped = dropped)
      will_abort <- !met && !is.null(abort_gd) && try == abort_try && isTRUE(rule > abort_gd)
      shown_try <- try; cull_note <- NULL
      if(!met && !will_abort && !is.null(cull_s) && n_culls < cull_s$max_culls){
        act <- cull_step(emc, "sample", cull_s, state); emc <- act$emc; state <- act$state
        if(length(act$done)){
          n_culls <- n_culls + 1L; state <- NULL
          emc <- cull_restart(emc, cull_s$reburn, FALSE, dots$verboseProgress, particle_factor, search_width,
                              cores_per_chain, cores_for_chains, dots$n_blocks, on_singular, context)
          target <- chain_n(emc)[1, "sample"] + iter; try <- 0L; restarts <- restarts + 1L
          cull_note <- paste0(act$note, "; sample stage restarted")
        }
      }
      save_fit(emc)
      elapsed <- Sys.time() - t0
      say("[sample | try=%d | iters=%d | %s%s%s%s] Duration: %s - ETA (to max_tries): %s", shown_try,
          r$draws, rhat_message(rule, r$gd),
          if(dropped > 0) sprintf(" (first %d draws dropped: the rule passes without them)", dropped) else "",
          if(is.null(stop_stats_message(st))) "" else stop_stats_message(st),
          if(is.null(cull_note)) "" else paste0(" | ", cull_note), format_duration(elapsed), eta(elapsed, max_tries - try))
      if(met) break
      if(will_abort){
        aborted <- list(try = try, rule = rule, abort_gd = abort_gd)
        warning(sprintf("Sample stage stopped early at check %d: the rule's Rhat %.3f is above abort_gd = %s. The fit is recorded as aborted (stop_record$aborted), not met.",
                        try, rule, format(abort_gd)), call. = FALSE)
        break
      }
    }
    if(!met && is.null(aborted)) warning("Max tries reached before all stopping criteria were met (see diagnose(emc)).", call. = FALSE)
    emc <- restore_duplicates(emc)
    flat <- tryCatch(flat_diagnostics(emc, "sample"), error = function(e) list())
    emc[[1]]$stop_record <- list(
      criteria = list(max_gd = max_gd, gd_quantile = gd_quantile, min_ess = min_ess, max_ll_gap = max_ll_gap,
                      max_subj_spread = if(single) NULL else max_subj_spread, iter = iter, abort_gd = abort_gd,
                      abort_try = if(is.null(abort_gd)) NULL else abort_try, cull = cull_s, max_tries = max_tries),
      met = met, tries = if(length(trace)) trace[[length(trace)]]$try else 0L, rule = rule,
      stats = c(st, list(flat = flat, prior = prior_tension(emc, "sample"))), culls = length(emc[[1]]$culls), discards = discards, aborted = aborted,
      restarts = restarts, trace = if(length(trace)) do.call(rbind, trace) else NULL, previous = previous,
      time = as.numeric(difftime(Sys.time(), t_start, units = "mins")))
    emc <- strip_duplicates(emc); class(emc) <- "emc"
    save_fit(emc)
  }
  emc
}

# Record the history_proposals setting on the object (NULL keeps the saved
# one). A change invalidates the sample kernel and the particle settings, and
# the latest draws become burn history so that adapt rebuilds them: the
# returned `stage` then says where to restart. Shared by fit() and hardfit().
set_history_proposals <- function(emc, history_proposals = NULL){
  previous <- isTRUE(emc[[1]]$history_proposals)
  if(is.null(history_proposals)) history_proposals <- previous
  if(!is.logical(history_proposals) || length(history_proposals) != 1 || is.na(history_proposals))
    stop("history_proposals must be TRUE or FALSE")
  stage <- NULL
  if(history_proposals != previous){
    emc[[1]]$sample_kernel <- NULL
    for(i in seq_along(emc)) attr(emc[[i]]$samples, "pm_settings") <- NULL
    if(get_last_stage(emc) %in% c("adapt", "sample")){
      emc <- subset(emc, stage = get_last_stage(emc))
      for(i in seq_along(emc)) emc[[i]]$samples$stage[] <- "burn"
      stage <- "adapt"
    }
  }
  emc[[1]]$history_proposals <- history_proposals
  list(emc = emc, stage = stage)
}

check_number <- function(x, name, min = -Inf, max = Inf){
  if(!is.numeric(x) || length(x) != 1 || is.na(x) || x < min || x > max)
    stop(sprintf("%s must be a single number%s", name,
                 if(is.finite(min) && is.finite(max)) sprintf(" in [%s, %s]", format(min), format(max))
                 else if(is.finite(min)) sprintf(" of at least %s", format(min)) else ""), call. = FALSE)
  invisible(TRUE)
}

# One sample-stage check: the Rhat rule (check_gd()), the stop statistics and
# the draw target, all on one set of draws. If the rule fails on every draw
# but passes with the first third dropped (check_gd()'s `short`), the shorter
# history is taken only if every criterion then passes; otherwise all draws
# are kept. A criterion given as NULL is off; a statistic that cannot be
# computed (NA) never passes.
sample_check <- function(emc, target, max_gd, gd_quantile, min_ess, max_ll_gap, max_subj_spread, selection){
  gd <- check_gd(emc, "sample", max_gd, NULL, TRUE, selection, gd_quantile)
  evaluate <- function(g){
    e <- g$emc; n <- unname(chain_n(e)[1, "sample"]); s <- sample_stop_stats(e, "sample")
    ok <- function(crit, stat, below){
      if(is.null(crit)) return(TRUE)
      v <- s[[stat]]; if(is.null(v) || !is.finite(v)) return(FALSE)
      if(below) v <= crit else v >= crit
    }
    stats_ok <- ok(min_ess, "ess_min", FALSE) && ok(max_ll_gap, "ll_gap", TRUE) && ok(max_subj_spread, "subj_spread", TRUE)
    list(emc = e, st = s, rule = g$rule, gd = g$gd, draws = n, stats_ok = stats_ok,
         met = isTRUE(g$gd_done) && stats_ok && n >= target, dropped = 0L)
  }
  r <- evaluate(gd)
  if(!r$met && !is.null(gd$short)){
    r_short <- evaluate(gd$short)
    if(r_short$met){ r_short$dropped <- r$draws - r_short$draws; r <- r_short }
  }
  class(r$emc) <- "emc"
  r
}

# One cull check: the comparison, the persistence bookkeeping and any respawn
cull_step <- function(emc, stage, cull, state){
  diag <- cull_diagnostics(emc, stage, cull)
  state <- update_cull_state(state, diag, length(emc))
  act <- apply_culls(emc, diag, stage, cull, state)
  note <- NULL
  if(length(act$done)){
    act$emc[[1]]$culls <- c(act$emc[[1]]$culls, act$done)
    note <- paste(vapply(act$done, function(d) sprintf("chain %d respawned from chain %d (log-likelihood %.0f below, z = %.1f)",
                                                       d$chain, d$donor, d$gap, d$z), ""), collapse = "; ")
  }
  list(emc = act$emc, state = state, done = act$done, note = note)
}
