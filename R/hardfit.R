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
#' changes a starting point, not the posterior. A note is printed when a
#' hierarchical object has fewer than 6 chains: the replacement, and the
#' detection of a subject with two modes, need healthy chains to compare with.
#'
#' **Early stop** (``abort_gd``, ``abort_try``): at check ``abort_try`` of the
#' sample stage, a fit whose rule Rhat is still above ``abort_gd`` ends there
#' with a warning and ``stop_record$aborted``. The defaults were calibrated on
#' replayed hard fits; see ``vignette("convergence")``.
#'
#' **The record** ``emc[[1]]$stop_record``: ``criteria`` (the thresholds that
#' were on), ``met``, ``tries``, ``rule`` (the Rhat the rule read last),
#' ``stats`` (the stop statistics of the last check, as described in
#' [diagnose()]), ``flat`` (the drift statistics) and ``prior`` (each group
#' mean's posterior mean in prior SDs from its prior mean), ``culls`` (count
#' for this call; the list is ``emc[[1]]$culls``), ``discards``, ``aborted``,
#' ``restarts``, ``trace`` (one row per check: draws, rule, plain max Rhat,
#' min ESS, gap, spread, culls, dropped), ``previous`` (the record of an
#' earlier ``hardfit()`` call on the same object) and ``time`` (minutes).
#' [diagnose()] reads it, and the draws, into one of the outcome classes.
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
#'   Rhats the rule reads; ``NULL`` (default) is .99, or 1 (the maximum) for
#'   ``type = "single"``.
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
                    max_gd = NULL, gd_quantile = NULL, min_ess = 100, max_ll_gap = 100,
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
  type <- emc[[1]]$type; single <- type == "single"
  rd <- rule_defaults(type, max_gd, gd_quantile)
  max_gd <- rd$max_gd; gd_quantile <- rd$gd_quantile; selection <- rd$selection
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
    message(sprintf("Note: %d chains; a subject with two modes often goes undetected below 6 (make_emc(n_chains = 6)); see vignette(\"convergence\").", length(emc)))
  hp <- set_history_proposals(emc, if(missing(history_proposals)) NULL else history_proposals)
  emc <- hp$emc
  if(is.null(stage)) stage <- get_last_stage(emc)
  if(!is.null(hp$stage)) stage <- hp$stage
  stages <- c("preburn", "burn", "adapt", "sample")
  if(!stage %in% stages) stop("stage must be one of ", paste(stages, collapse = ", "))
  to_run <- stages[which(stages == stage):length(stages)]
  context <- sampler_contexts(restore_duplicates(emc), cores_per_chain)
  # one block of `n` more draws of `stage`, one exact step of the stage's kernel:
  # run_emc() with an iteration-only criterion computes no Rhat. Burn and sample
  # blocks are silent (hardfit() prints its own line per check); preburn and
  # adapt run whole and print fit()'s lines.
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

  # preburn and adapt run whole: adapt's convergence test keeps a history across
  # its checks (adapt_converged()'s stall count) and builds the sample kernel
  # once at its end, which one block at a time would reset
  if("preburn" %in% to_run){
    emc <- block(emc, "preburn", step_size, get_stop_criteria("preburn", NULL, type), max_tries, say_progress = TRUE)
    save_fit(emc)
  }

  if("burn" %in% to_run){
    say("Running burn stage")
    crit <- get_stop_criteria("burn", list(mean_gd = burn$mean_gd, selection = selection, omit_mpsrf = TRUE), type)
    cull_b <- cull_settings(cull, "burn"); state <- NULL; met <- FALSE
    for(try in seq_len(burn$max_tries)){
      t0 <- Sys.time()
      emc <- block(emc, "burn", step_size)
      gd <- check_gd(emc, "burn", NULL, crit$mean_gd, TRUE, selection)
      dropped <- 0L
      # burn has no other criterion, so the shorter history is taken whenever its mean Rhat passes
      if(!isTRUE(gd$gd_done) && !is.null(gd$short)){
        dropped <- unname(chain_n(emc)[1, "burn"] - chain_n(gd$short$emc)[1, "burn"])
        emc <- gd$short$emc; gd <- gd$short
      }
      met <- isTRUE(gd$gd_done)
      cull_note <- NULL
      if(!met && !is.null(cull_b)){
        act <- cull_step(emc, "burn", cull_b, state); emc <- act$emc; state <- act$state; cull_note <- act$note
      }
      save_fit(emc)
      elapsed <- Sys.time() - t0
      say("[burn | try=%d | iters=%d | Mean Rhat=%.3f%s%s] Duration: %s - ETA (to max_tries): %s", try, chain_n(emc)[1, "burn"],
          mean(gd$gd), if(dropped > 0) sprintf(" (first %d draws dropped)", dropped) else "",
          if(is.null(cull_note)) "" else paste0(" | ", cull_note), format_duration(elapsed), eta(elapsed, burn$max_tries - try))
      if(met) break
    }
    if(!met) warning("Burn stage ended at max_tries with its Rhat rule unmet.", call. = FALSE)
    # a chain still apart when burn ends would seed adapt's pooled proposals with its region
    if(!is.null(cull_b)) emc <- burn_exit_reseed(emc, cull_b, verbose)
    save_fit(emc)
  }

  if("adapt" %in% to_run){
    emc <- block(emc, "adapt", step_size, get_stop_criteria("adapt", NULL, type), max_tries, say_progress = TRUE)
    save_fit(emc)
  }

  if("sample" %in% to_run){
    say("Running sample stage")
    previous <- emc[[1]]$stop_record; culls_before <- length(emc[[1]]$culls)
    cull_s <- cull_settings(cull, "sample"); state <- NULL
    criteria <- list(max_gd = max_gd, gd_quantile = gd_quantile, min_ess = min_ess, max_ll_gap = max_ll_gap,
                     max_subj_spread = if(single) NULL else max_subj_spread, iter = iter, abort_gd = abort_gd,
                     abort_try = if(is.null(abort_gd)) NULL else abort_try, cull = cull_s, max_tries = max_tries)
    target <- chain_n(emc)[1, "sample"] + iter   # iter counts new draws: a call on a fitted object continues its sample stage
    discards <- list(); aborted <- NULL; restarts <- 0L; trace <- list(); met <- FALSE; st <- NULL; rule <- NA_real_
    try <- 0L
    while(try < max_tries){
      try <- try + 1L; t0 <- Sys.time()
      emc <- block(emc, "sample", step_size)
      r <- sample_check(emc, target, max_gd, gd_quantile, min_ess, max_ll_gap, criteria$max_subj_spread, selection)
      emc <- r$emc; st <- r$st; rule <- r$rule; met <- r$met
      if(r$dropped > 0) discards <- c(discards, list(list(try = try, removed = r$dropped, kept = r$draws)))
      # a one-time test at check abort_try: a fit under abort_gd there keeps running
      # whatever its Rhat does later (the rule was calibrated at one check)
      will_abort <- !met && !is.null(abort_gd) && try == abort_try && isTRUE(rule > abort_gd)
      shown_try <- try; cull_note <- NULL
      if(!met && !will_abort && !is.null(cull_s)){   # no respawn in a fit that stops at this check
        act <- cull_step(emc, "sample", cull_s, state); emc <- act$emc; state <- act$state; cull_note <- act$note
        if(length(act$done)){
          # the kept sample draws must all come from the chains' final states: the
          # stage restarts after reburn draws, and the check count with it
          # (abort_try and max_tries count from the restart)
          emc <- cull_restart(emc, cull_s$reburn, function(e, n) block(e, "sample", n))
          target <- chain_n(emc)[1, "sample"] + iter; try <- 0L; restarts <- restarts + 1L
          cull_note <- paste0(cull_note, "; sample stage restarted")
        }
      }
      trace[[length(trace) + 1]] <- cbind(try = shown_try, r$trace, culls = if(is.null(state)) 0L else state$n_culls, dropped = r$dropped)
      save_fit(emc)
      elapsed <- Sys.time() - t0
      stats_note <- stop_stats_message(st)
      say("[sample | try=%d | iters=%d | %s%s%s%s] Duration: %s - ETA (to max_tries): %s", shown_try,
          r$draws, rhat_message(rule, r$gd),
          if(r$dropped > 0) sprintf(" (first %d draws dropped: the rule passes without them)", r$dropped) else "",
          if(is.null(stats_note)) "" else stats_note,
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
    emc[[1]]$stop_record <- list(
      criteria = criteria, met = met,
      tries = if(length(trace)) trace[[length(trace)]]$try else 0L,   # checks since the last restart; the trace holds every check
      rule = rule, stats = st, flat = drift_report(emc, "sample"), prior = prior_tension(emc, "sample"),
      culls = length(emc[[1]]$culls) - culls_before, discards = discards, aborted = aborted, restarts = restarts,
      trace = if(length(trace)) do.call(rbind, trace) else NULL, previous = previous,
      time = as.numeric(difftime(Sys.time(), t_start, units = "mins")))
    emc <- strip_duplicates(emc); class(emc) <- "emc"
    save_fit(emc)
  }
  emc
}

check_number <- function(x, name, min = -Inf, max = Inf){
  if(!is.numeric(x) || length(x) != 1 || is.na(x) || x < min || x > max)
    stop(sprintf("%s must be a single number%s", name,
                 if(is.finite(min) && is.finite(max)) sprintf(" in [%s, %s]", format(min), format(max))
                 else if(is.finite(min)) sprintf(" of at least %s", format(min)) else ""), call. = FALSE)
  invisible(TRUE)
}

# One sample check on one set of draws. check_gd()'s shorter history (`short`)
# is taken only when every criterion then passes: draws dropped from the
# object are gone for good. Returns the object, the statistics, the decision
# and the check's trace row.
sample_check <- function(emc, target, max_gd, gd_quantile, min_ess, max_ll_gap, max_subj_spread, selection){
  gd <- check_gd(emc, "sample", max_gd, NULL, TRUE, selection, gd_quantile)
  evaluate <- function(g){
    e <- g$emc; n <- unname(chain_n(e)[1, "sample"]); s <- sample_stop_stats(e, "sample")
    ok <- stop_rules(s, min_ess, max_ll_gap, max_subj_spread)
    list(emc = e, st = s, rule = g$rule, gd = g$gd, draws = n, met = isTRUE(g$gd_done) && all(ok) && n >= target, dropped = 0L,
         trace = data.frame(draws = n, rule = g$rule, max_rhat = if(length(g$gd)) max(g$gd) else NA_real_,
                            ess_min = stat_or_na(s$ess_min), ll_gap = stat_or_na(s$ll_gap), subj_spread = stat_or_na(s$subj_spread)))
  }
  r <- evaluate(gd)
  if(!r$met && !is.null(gd$short)){
    r_short <- evaluate(gd$short)
    if(r_short$met){ r_short$dropped <- r$draws - r_short$draws; r <- r_short }
  }
  r
}
