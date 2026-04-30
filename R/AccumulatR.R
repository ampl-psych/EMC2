`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

accumulatr_dist_param_names <- function(dist) {
  dist <- tolower(dist)
  if (dist == "lognormal") return(c("m", "s"))
  if (dist == "gamma") return(c("shape", "rate"))
  if (dist == "exgauss") return(c("mu", "sigma", "tau"))
  if (dist == "lba") return(c("v", "B", "A", "sv"))
  if (dist == "rdm") return(c("v", "B", "A", "s"))
  character(0)
}

accumulatr_param_profile <- function(param_name, runtime_family) {
  if (identical(runtime_family, "q")) {
    return(list(
      latent_default = qnorm(0),
      actual_default = 0,
      func = "pnorm",
      lower = 0,
      upper = 1,
      bound_min = 0,
      bound_max = 1,
      exception = 0
    ))
  }

  if (identical(runtime_family, "t0")) {
    return(list(
      latent_default = log(0),
      actual_default = 0,
      func = "exp",
      lower = 0,
      upper = Inf,
      bound_min = 0.05,
      bound_max = Inf,
      exception = 0
    ))
  }

  if (identical(runtime_family, "w")) {
    return(list(
      latent_default = qnorm(0.5),
      actual_default = 0.5,
      func = "pnorm",
      lower = 0,
      upper = 1,
      bound_min = 0,
      bound_max = 1,
      exception = 0
    ))
  }

  exp_params <- c("s", "sigma", "tau", "shape", "rate", "B", "A", "sv")
  if (param_name %in% exp_params) {
    return(list(
      latent_default = 0,
      actual_default = 1,
      func = "exp",
      lower = 0,
      upper = Inf,
      bound_min = 0,
      bound_max = Inf,
      exception = 0
    ))
  }

  list(
    latent_default = 0,
    actual_default = 0,
    func = "identity",
    lower = -Inf,
    upper = Inf,
    bound_min = -Inf,
    bound_max = Inf,
    exception = NA_real_
  )
}

accumulatr_profile_key <- function(profile) {
  paste(
    profile$func,
    profile$lower,
    profile$upper,
    profile$bound_min,
    profile$bound_max,
    ifelse(is.na(profile$exception), "NA", as.character(profile$exception)),
    sep = "::"
  )
}

accumulatr_component_members <- function(component) {
  members <- component$members %||% character(0)
  as.character(members)
}

accumulatr_component_weight_param <- function(component) {
  weight_param <- component$weight_param %||% component$attrs$weight_param %||% NULL
  if (is.null(weight_param) || !nzchar(weight_param)) {
    return(NULL)
  }
  as.character(weight_param)
}

accumulatr_model_spec <- function(model_spec) {
  model_spec$model_spec %||% model_spec
}

accumulatr_internal_lookup <- function(model_spec) {
  spec <- accumulatr_model_spec(model_spec)
  AccumulatR:::.parameter_name_lookup(spec)
}

accumulatr_trigger_ids <- function(model_spec) {
  spec <- accumulatr_model_spec(model_spec)
  ids <- vapply(spec$triggers %||% list(), function(trig) {
    trig$id %||% NA_character_
  }, character(1))
  ids[!is.na(ids) & nzchar(ids)]
}

accumulatr_shared_trigger_sources <- function(model_spec) {
  spec <- accumulatr_model_spec(model_spec)
  out <- list()
  for (trig in spec$triggers %||% list()) {
    trig_id <- trig$id %||% NA_character_
    if (is.na(trig_id) || !nzchar(trig_id)) next
    if (!identical(trig$draw %||% "shared", "shared")) next
    for (member in as.character(trig$members %||% character(0))) {
      out[[member]] <- trig_id
    }
  }
  out
}

accumulatr_shared_trigger_member_q_names <- function(model_spec) {
  members <- names(accumulatr_shared_trigger_sources(model_spec))
  if (!length(members)) return(character(0))
  paste(members, "q", sep = ".")
}

accumulatr_internal_runtime_family <- function(model_spec, internal_name) {
  spec <- accumulatr_model_spec(model_spec)
  components <- spec$components %||% list()
  for (component in components) {
    weight_param <- accumulatr_component_weight_param(component)
    if (!is.null(weight_param) && identical(internal_name, weight_param)) {
      return("w")
    }
  }

  if (internal_name %in% accumulatr_trigger_ids(spec)) {
    return("q")
  }

  pieces <- strsplit(internal_name, ".", fixed = TRUE)[[1]]
  if (length(pieces) != 2L) {
    stop(
      "Could not determine runtime family for AccumulatR parameter '",
      internal_name,
      "'",
      call. = FALSE
    )
  }
  acc_id <- pieces[[1]]
  param_name <- pieces[[2]]
  if (identical(param_name, "q")) {
    return("q")
  }
  if (identical(param_name, "t0")) {
    return("t0")
  }

  acc_defs <- spec$accumulators %||% list()
  acc_idx <- match(acc_id, vapply(acc_defs, `[[`, character(1), "id"))
  if (is.na(acc_idx)) {
    stop(
      "Unknown AccumulatR accumulator '",
      acc_id,
      "' in parameter '",
      internal_name,
      "'",
      call. = FALSE
    )
  }
  dist <- acc_defs[[acc_idx]]$dist
  base_names <- get("dist_param_names", envir = asNamespace("AccumulatR"))(dist)
  slot_idx <- match(param_name, base_names)
  if (is.na(slot_idx)) {
    stop(
      "Could not determine runtime slot for AccumulatR parameter '",
      internal_name,
      "'",
      call. = FALSE
    )
  }
  paste0("p", slot_idx)
}

accumulatr_internal_profiles <- function(model_spec, internal_names) {
  runtime_family <- stats::setNames(character(length(internal_names)), internal_names)
  profile <- vector("list", length(internal_names))
  names(profile) <- internal_names
  for (internal_name in internal_names) {
    family <- accumulatr_internal_runtime_family(model_spec, internal_name)
    runtime_family[[internal_name]] <- family
    param_name <- sub("^.*\\.", "", internal_name)
    if (identical(family, "q")) {
      param_name <- "q"
    }
    if (identical(family, "w")) {
      param_name <- "w"
    }
    profile[[internal_name]] <- accumulatr_param_profile(param_name, family)
  }
  list(
    internal_to_runtime = runtime_family,
    profiles = profile
  )
}

AccumulatR_bridge_spec <- function(model_spec) {
  internal_to_public <- accumulatr_internal_lookup(model_spec)
  internal_names <- names(internal_to_public)
  target_info <- accumulatr_internal_profiles(model_spec, internal_names)
  hidden_internal <- accumulatr_shared_trigger_member_q_names(model_spec)

  public_order <- AccumulatR::sampled_pars(model_spec)
  if (length(hidden_internal) > 0L) {
    public_order <- public_order[vapply(public_order, function(nm) {
      targets <- names(internal_to_public)[unname(internal_to_public) %in% nm]
      !(length(targets) > 0L && all(targets %in% hidden_internal))
    }, logical(1))]
  }
  public_to_internal <- lapply(public_order, function(nm) {
    names(internal_to_public)[unname(internal_to_public) %in% nm]
  })
  names(public_to_internal) <- public_order

  public_defaults <- stats::setNames(numeric(length(public_order)), public_order)
  public_funcs <- stats::setNames(character(length(public_order)), public_order)
  public_lower <- stats::setNames(numeric(length(public_order)), public_order)
  public_upper <- stats::setNames(numeric(length(public_order)), public_order)
  public_min <- stats::setNames(numeric(length(public_order)), public_order)
  public_max <- stats::setNames(numeric(length(public_order)), public_order)
  public_exc <- stats::setNames(rep(NA_real_, length(public_order)), public_order)
  public_to_runtime <- stats::setNames(character(length(public_order)), public_order)
  public_profiles <- vector("list", length(public_order))
  names(public_profiles) <- public_order

  for (nm in public_order) {
    targets <- public_to_internal[[nm]]
    families <- unique(unname(target_info$internal_to_runtime[targets]))
    if (length(families) != 1L) {
      stop(
        "AccumulatR public parameter '", nm,
        "' must map to exactly one runtime slot family; got ",
        paste(families, collapse = ", "),
        call. = FALSE
      )
    }
    profiles <- target_info$profiles[targets]
    keys <- unique(vapply(profiles, accumulatr_profile_key, character(1)))
    if (length(keys) > 1L) {
      warning(
        "AccumulatR public parameter '", nm,
        "' maps to different parameter types; using the first profile",
        call. = FALSE
      )
    }
    profile <- profiles[[1]]
    public_to_runtime[[nm]] <- families[[1]]
    public_profiles[[nm]] <- profile
    public_defaults[[nm]] <- profile$latent_default
    public_funcs[[nm]] <- profile$func
    public_lower[[nm]] <- profile$lower
    public_upper[[nm]] <- profile$upper
    public_min[[nm]] <- profile$bound_min
    public_max[[nm]] <- profile$bound_max
    public_exc[[nm]] <- profile$exception
  }

  list(
    public_names = public_order,
    p_types = public_defaults,
    transform = list(
      func = public_funcs,
      lower = public_lower,
      upper = public_upper
    ),
    bound = list(
      minmax = rbind(public_min, public_max),
      exception = public_exc[!is.na(public_exc)]
    ),
    public_to_internal = public_to_internal,
    public_to_runtime = public_to_runtime,
    public_profiles = public_profiles,
    internal_to_runtime = target_info$internal_to_runtime,
    internal_to_public = internal_to_public,
    internal_profiles = target_info$profiles
  )
}

accumulatr_trigger_map <- function(model_spec) {
  trig_map <- list()
  triggers <- model_spec$triggers %||% list()
  for (trig in triggers) {
    members <- trig$members %||% character(0)
    if (!length(members)) next
    q_val <- if (!is.null(trig$q)) as.numeric(trig$q) else NA_real_
    draw_mode <- trig$draw %||% "shared"
    for (member in members) {
      trig_map[[member]] <- list(q = q_val, draw = draw_mode)
    }
  }
  trig_map
}

accumulatr_runtime_columns <- function(model_spec) {
  acc_defs <- accumulatr_model_spec(model_spec)$accumulators %||% list()
  max_p <- max(vapply(acc_defs, function(acc) {
    length(get("dist_param_names", envir = asNamespace("AccumulatR"))(acc$dist))
  }, integer(1)), 0L)
  comp_defs <- accumulatr_model_spec(model_spec)$components %||% list()
  weight_params <- unique(vapply(comp_defs, function(component) {
    accumulatr_component_weight_param(component) %||% NA_character_
  }, character(1)))
  weight_params <- weight_params[!is.na(weight_params)]
  c("q", "t0", paste0("p", seq_len(max_p)), weight_params)
}

accumulatr_actual_default <- function(profile) {
  profile$actual_default %||% 0
}

accumulatr_internal_source_name <- function(bridge, internal_name) {
  public_name <- bridge$internal_to_public[[internal_name]]
  if (is.null(public_name) || !length(public_name) || !nzchar(public_name)) {
    return(NA_character_)
  }
  public_name
}

accumulatr_build_runtime_recipe <- function(model_spec, dadm, bridge) {
  acc_defs <- model_spec$prep$accumulators
  acc_ids <- names(acc_defs)
  acc_by_id <- stats::setNames(acc_defs, acc_ids)
  runtime_cols <- accumulatr_runtime_columns(model_spec)
  defaults <- matrix(0, nrow = nrow(dadm), ncol = length(runtime_cols))
  colnames(defaults) <- runtime_cols
  sources <- matrix(NA_character_, nrow = nrow(dadm), ncol = length(runtime_cols))
  colnames(sources) <- runtime_cols
  comp_defs <- accumulatr_model_spec(model_spec)$components %||% list()
  shared_trigger_sources <- accumulatr_shared_trigger_sources(model_spec)

  for (i in seq_len(nrow(dadm))) {
    acc_id <- as.character(dadm$accumulator[[i]])
    acc <- acc_by_id[[acc_id]]
    if (is.null(acc)) {
      stop("Prepared AccumulatR data contains unknown accumulator '", acc_id, "'", call. = FALSE)
    }

    dist_params <- get("dist_param_names", envir = asNamespace("AccumulatR"))(acc$dist)
    acc_params <- acc$params %||% list()
    for (j in seq_along(dist_params)) {
      runtime_col <- paste0("p", j)
      internal_name <- paste(acc_id, dist_params[[j]], sep = ".")
      profile <- bridge$internal_profiles[[internal_name]] %||%
        accumulatr_param_profile(dist_params[[j]], runtime_col)
      defaults[i, runtime_col] <- acc_params[[dist_params[[j]]]] %||%
        accumulatr_actual_default(profile)
      sources[i, runtime_col] <- accumulatr_internal_source_name(bridge, internal_name)
    }

    q_name <- paste(acc_id, "q", sep = ".")
    q_source <- shared_trigger_sources[[acc_id]] %||% q_name
    q_profile <- bridge$internal_profiles[[q_source]] %||%
      bridge$internal_profiles[[q_name]] %||%
      accumulatr_param_profile("q", "q")
    q_default <- acc$q %||% accumulatr_actual_default(q_profile)
    defaults[i, "q"] <- q_default
    sources[i, "q"] <- accumulatr_internal_source_name(bridge, q_source)

    t0_name <- paste(acc_id, "t0", sep = ".")
    t0_profile <- bridge$internal_profiles[[t0_name]] %||% accumulatr_param_profile("t0", "t0")
    defaults[i, "t0"] <- acc_params$t0 %||% accumulatr_actual_default(t0_profile)
    sources[i, "t0"] <- accumulatr_internal_source_name(bridge, t0_name)
  }

  for (component in comp_defs) {
    weight_param <- accumulatr_component_weight_param(component)
    if (is.null(weight_param) || !weight_param %in% colnames(defaults)) next
    defaults[, weight_param] <- component$weight %||% 0
    sources[, weight_param] <- accumulatr_internal_source_name(bridge, weight_param)
  }

  list(
    defaults = defaults,
    source_names = sources
  )
}

accumulatr_apply_runtime_recipe <- function(public_pars, recipe) {
  runtime <- recipe$defaults
  public_pars <- as.matrix(public_pars)
  for (j in seq_len(ncol(runtime))) {
    src <- recipe$source_names[, j]
    idx <- match(src, colnames(public_pars))
    take <- !is.na(idx)
    if (any(take)) {
      runtime[take, j] <- public_pars[cbind(which(take), idx[take])]
    }
  }
  runtime
}

accumulatr_runtime_params <- function(model_spec, trial_df, public_pars, bridge = NULL) {
  bridge <- bridge %||% AccumulatR_bridge_spec(model_spec)
  recipe <- accumulatr_build_runtime_recipe(model_spec, trial_df, bridge)
  accumulatr_apply_runtime_recipe(public_pars, recipe)
}

accumulatr_sim_runtime_params <- function(model_spec, trial_df, public_pars, bridge = NULL) {
  runtime <- accumulatr_runtime_params(model_spec, trial_df, public_pars, bridge)
  if ("w" %in% colnames(runtime)) return(runtime)

  w <- rep(1, nrow(runtime))
  comp_defs <- accumulatr_model_spec(model_spec)$components %||% list()
  if (length(comp_defs) > 0L && "accumulator" %in% names(trial_df)) {
    acc_weight <- list()
    for (component in comp_defs) {
      weight_param <- accumulatr_component_weight_param(component)
      weight <- component$weight %||% 1
      if (!is.null(weight_param) && weight_param %in% colnames(runtime)) {
        weight <- runtime[, weight_param]
      }
      for (member in accumulatr_component_members(component)) {
        acc_weight[[member]] <- weight
      }
    }
    acc <- as.character(trial_df$accumulator)
    for (i in seq_along(acc)) {
      weight <- acc_weight[[acc[[i]]]]
      if (!is.null(weight)) w[[i]] <- weight[[min(i, length(weight))]]
    }
  }

  cbind(
    runtime[, "q", drop = FALSE],
    w = w,
    runtime[, setdiff(colnames(runtime), "q"), drop = FALSE]
  )
}

AccumulatR_model <- function(model_spec){
  bridge <- AccumulatR_bridge_spec(model_spec)
  model_list <- list(
    c_name = "AccumulatR",
    type = "AccumulatR",
    p_types = bridge$p_types,
    transform = bridge$transform,
    Ttransform = function(pars,dadm) pars,
    bound = bridge$bound,
    rfun = function(trial_df, pars) {
      runtime_pars <- accumulatr_sim_runtime_params(model_spec, trial_df, pars, bridge)
      trial_sim <- trial_df
      if ("trial" %in% names(trial_sim) && nrow(trial_sim) > 0L) {
        trial <- as.integer(trial_sim$trial)
        trial_sim$trial <- cumsum(c(TRUE, trial[-1L] != trial[-length(trial)]))
      }
      AccumulatR::simulate(model_spec, runtime_pars, trial_df = trial_sim, keep_component = TRUE, layout = "long")
    },
    spec = model_spec,
    accumulatr_bridge = bridge
  )
  return(function() {return(model_list)})
}

AccumulatR_expand_data <- function(model_spec, data){
  if(is.null(model_spec)) stop("model_specification needs to be made for AccumulatR models")
  if(is.data.frame(model_spec$components) && nrow(model_spec$components) > 1){
    if(is.null(data$component)){
      warning("The AccumulatR model contains multiple components. Specify the components
              column in your data if you do not want to treat them as mixtures
              marginalized in the likelihood")
      cat("expected components :", model_spec$components$component_id)
    }
  }
  if(is.null(data$rt)) data$rt <- rnorm(nrow(data))
  datar <- prepare_data(model_spec, data)
  return(datar)
}

AccumulatR_add_context <- function(dadm){
  model <- attr(dadm, "model")
  if(is.null(model)) return(dadm)
  if(model()$type != "AccumulatR") return(dadm)
  model_list <- model()
  model_spec <- model_list$spec
  expand <- attr(dadm, "expand", exact = TRUE)
  if(is.null(dadm$accumulator) && !is.null(dadm$lR)){
    dadm$accumulator <- as.character(dadm$lR)
  }
  dadm <- AccumulatR::prepare_data(model_spec, dadm)
  if(!is.null(expand)) attr(dadm, "expand") <- expand
  ctx <- make_context(model_spec)
  bridge <- model_list$accumulatr_bridge
  ctx$native_ctx <- ctx$cpp$native
  ctx$bridge <- accumulatr_build_runtime_recipe(model_spec, dadm, bridge)
  attr(dadm, "AccumulatR_context") <- ctx
  return(dadm)
}

AccumulatR_check_context <- function(emc){
  if(emc[[1]]$model()$type != "AccumulatR") return(emc)
  for(i in 1:length(emc)){
    for(j in 1:length(emc[[i]]$data)){
      dat <- emc[[i]]$data[[j]]
      if(is.data.frame(dat)){
        attr(dat, "model") <- emc[[i]]$model
        dat <- AccumulatR_add_context(dat)
        attr(dat, "model") <- NULL
        emc[[i]]$data[[j]] <- dat
      }
    }
  }
  return(emc)
}
