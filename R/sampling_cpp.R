# The R boundary resolves fixed likelihood inputs once per fit.
native_models <- function(model) {
  if (is.function(model) || !is.null(model$c_name) || !is.null(model$log_likelihood)) model <- list(model)
  lapply(model, function(x) if (is.function(x)) x() else x)
}

# Share the parameter-to-likelihood mapping between blocks and likelihoods.
shared_likelihood_index <- function(pmwgs, models = native_models(pmwgs$model)) {
  shared <- attr(pmwgs$data, "shared_ll_idx")
  if (!is.null(shared)) return(shared)
  par_names <- pmwgs$par_names
  if (length(models) > 1L) match(sub("[|].*", "", par_names),
                               unique(sub("[|].*", "", par_names)))
  else rep(1L, length(par_names))
}

native_likelihood_inputs <- function(pmwgs) {
  models <- native_models(pmwgs$model)
  if (any(vapply(models, function(model) is.null(model$c_name), logical(1))))
    stop("Sampling requires a compiled likelihood")
  par_names <- pmwgs$par_names
  shared <- shared_likelihood_index(pmwgs, models)
  joint <- length(models) > 1L
  lapply(pmwgs$data, function(data) {
    dadms <- if (is.data.frame(data)) list(data) else data
    specs <- vector("list", length(dadms))
    for (i in seq_along(dadms)) {
      dadm <- dadms[[i]]
      if (!is.data.frame(dadm)) next
      model <- models[[i]]
      designs <- attr(dadm, "designs")
      if (is.numeric(designs)) attr(dadm, "designs") <- attr(dadms[[designs]], "designs")
      columns <- which(shared == i)
      constants <- attr(dadm, "constants")
      # neural likelihoods: the registration and the live evaluator, as calc_ll_manager()
      nn <- if (identical(model$c_name, "NN")) nn_native_args(model$nn, constants)
      if (is.null(constants)) constants <- NA_real_
      specs[[i]] <- list(data = dadm, model = model, constants = constants,
        designs = get_designs_expanded(dadm, model), columns = columns,
        component = i, par_names = if (joint) sub(".*[|]", "", par_names[columns]) else par_names[columns],
        nn = nn)
    }
    Filter(Negate(is.null), specs)
  })
}

prepare_native_stage <- function(pmwgs, stage, iter, particles, tune, on_singular) {
  if (pmwgs$type == "SEM") pmwgs <- normalise_SEM_prior(pmwgs)
  tune$components <- attr(pmwgs$data, "components")
  tune$shared_ll_idx <- shared_likelihood_index(pmwgs)
  kernel <- if (is.null(tune$kernel)) stage else tune$kernel
  tune$history_proposals <- isTRUE(pmwgs$history_proposals) && kernel %in% c("adapt", "sample")
  tune <- check_tune_settings(tune, pmwgs$n_pars, kernel, particles)
  tune$kernel <- kernel
  tune$frozen <- stage == "sample"
  tune$exact <- kernel %in% c("adapt", "sample")
  tune$lik_prec <- pmwgs$lik_prec
  tune$min_local <- if (tune$exact && kernel == "sample" && pmwgs$type != "single") .25 else 0
  tune$do_scale <- pmwgs$type != "single" && kernel %in% scale_move_stages()
  pm_settings <- attr(pmwgs$samples, "pm_settings")
  if (is.null(pm_settings)) {
    pm_settings <- replicate(pmwgs$n_subjects,
      vector("list", length(unique(tune$components))), simplify = FALSE)
  }
  pm_settings <- lapply(pm_settings, check_sampling_settings, stage = kernel,
                        n_pars = pmwgs$n_pars, particles = particles,
                        history_proposals = tune$history_proposals)
  scale_settings <- attr(pmwgs$samples, "scale_move")
  if (tune$do_scale) scale_settings <- scale_move_init(scale_settings, pmwgs$par_names)
  list(sampler = extend_sampler(pmwgs, iter, stage), tune = tune,
       pm_settings = pm_settings, scale_settings = scale_settings,
       on_singular = resolve_on_singular(on_singular))
}

sampler_contexts <- function(emc, n_cores) {
  native_sampler_contexts(native_likelihood_inputs(emc[[1]]), length(emc), n_cores)
}

run_native_chains <- function(emc, stage, iter, particle_factor = 50,
                              search_width = 1, n_cores = 1,
                              cores_for_chains = length(emc), on_singular = NULL,
                              kernel = NULL, verbose = FALSE, verboseProgress = FALSE,
                              context = sampler_contexts(emc, n_cores)) {
  emc <- lapply(emc, function(x) if (x$type == "SEM") normalise_SEM_prior(x) else x)
  if (any(!vapply(emc, function(x) x$init, logical(1))))
    emc <- native_init_chains(emc, context, NULL, NULL, 1000L, FALSE)
  if (!iter) return(emc)
  args <- lapply(emc, function(x) prepare_native_stage(x, stage, iter,
    round(particle_factor * sqrt(x$n_pars)),
    list(search_width = search_width, kernel = kernel), on_singular))
  native_run_chains(args, context, stage, iter, cores_for_chains, verboseProgress)
}
