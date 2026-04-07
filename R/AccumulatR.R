accumulatr_empty_transform <- function() {
  list(
    func = stats::setNames(character(0), character(0)),
    lower = stats::setNames(numeric(0), character(0)),
    upper = stats::setNames(numeric(0), character(0))
  )
}

accumulatr_identity_transform <- function(param_names) {
  list(
    func = stats::setNames(rep("identity", length(param_names)), param_names),
    lower = stats::setNames(rep(-Inf, length(param_names)), param_names),
    upper = stats::setNames(rep(Inf, length(param_names)), param_names)
  )
}

accumulatr_open_bound <- function(param_names) {
  minmax <- rbind(
    stats::setNames(rep(-Inf, length(param_names)), param_names),
    stats::setNames(rep(Inf, length(param_names)), param_names)
  )
  list(
    minmax = minmax,
    exception = stats::setNames(numeric(0), character(0))
  )
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
      default = qnorm(0),
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
      default = log(0),
      func = "exp",
      lower = 0,
      upper = Inf,
      bound_min = 0.05,
      bound_max = Inf,
      exception = 0
    ))
  }

  exp_params <- c("s", "sigma", "tau", "shape", "rate", "B", "A", "sv")
  if (param_name %in% exp_params) {
    return(list(
      default = 0,
      func = "exp",
      lower = 0,
      upper = Inf,
      bound_min = 0,
      bound_max = Inf,
      exception = 0
    ))
  }

  list(
    default = 0,
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

accumulatr_internal_targets <- function(model_spec) {
  acc_defs <- model_spec$prep$accumulators
  internal <- character(0)
  family <- character(0)
  profile <- list()
  owner <- list()

  for (acc in acc_defs) {
    acc_id <- acc$id
    base_names <- accumulatr_dist_param_names(acc$dist)
    for (idx in seq_along(base_names)) {
      nm <- paste(acc_id, base_names[[idx]], sep = ".")
      fam <- paste0("p", idx)
      internal <- c(internal, nm)
      family <- c(family, fam)
      profile[[nm]] <- accumulatr_param_profile(base_names[[idx]], fam)
      owner[[nm]] <- acc_id
    }
    q_nm <- paste(acc_id, "q", sep = ".")
    internal <- c(internal, q_nm)
    family <- c(family, "q")
    profile[[q_nm]] <- accumulatr_param_profile("q", "q")
    owner[[q_nm]] <- acc_id

    t0_nm <- paste(acc_id, "t0", sep = ".")
    internal <- c(internal, t0_nm)
    family <- c(family, "t0")
    profile[[t0_nm]] <- accumulatr_param_profile("t0", "t0")
    owner[[t0_nm]] <- acc_id
  }

  comp_defs <- model_spec$components
  if (is.data.frame(comp_defs) && nrow(comp_defs) > 0) {
    for (i in seq_len(nrow(comp_defs))) {
      weight_param <- NULL
      if ("weight_param" %in% names(comp_defs)) {
        weight_param <- comp_defs$weight_param[[i]]
      }
      if ((is.null(weight_param) || !nzchar(weight_param)) && "attrs" %in% names(comp_defs)) {
        attr_list <- comp_defs$attrs[[i]]
        if (is.list(attr_list) && !is.null(attr_list$weight_param)) {
          weight_param <- attr_list$weight_param
        }
      }
      if (is.null(weight_param) || !nzchar(weight_param)) next
      internal <- c(internal, weight_param)
      family <- c(family, weight_param)
      profile[[weight_param]] <- list(
        default = 0,
        func = "identity",
        lower = -Inf,
        upper = Inf,
        bound_min = -Inf,
        bound_max = Inf,
        exception = NA_real_
      )
      owner[[weight_param]] <- NA_character_
    }
  }

  names(family) <- internal
  list(
    internal_names = unique(internal),
    internal_to_runtime = family[unique(internal)],
    profiles = profile[unique(internal)],
    owners = owner[unique(internal)]
  )
}

accumulatr_public_lookup <- function(model_spec, internal_names) {
  lookup <- stats::setNames(internal_names, internal_names)
  mappings <- model_spec$parameters
  if (is.null(mappings) && !is.null(model_spec$model_spec)) {
    mappings <- model_spec$model_spec$parameters
  }
  if (is.null(mappings) || length(mappings) == 0) {
    return(lookup)
  }

  seen_targets <- character(0)
  for (external_name in names(mappings)) {
    targets <- as.character(mappings[[external_name]])
    unknown <- setdiff(targets, internal_names)
    if (length(unknown) > 0) {
      stop(
        "AccumulatR model references unknown parameter(s): ",
        paste(unknown, collapse = ", "),
        call. = FALSE
      )
    }
    duplicated_targets <- intersect(seen_targets, targets)
    if (length(duplicated_targets) > 0) {
      stop(
        "AccumulatR internal parameters cannot be assigned more than once: ",
        paste(duplicated_targets, collapse = ", "),
        call. = FALSE
      )
    }
    lookup[targets] <- external_name
    seen_targets <- c(seen_targets, targets)
  }
  lookup
}

AccumulatR_bridge_spec <- function(model_spec) {
  target_info <- accumulatr_internal_targets(model_spec)
  internal_names <- target_info$internal_names
  lookup <- accumulatr_public_lookup(model_spec, internal_names)

  public_order <- AccumulatR::sampled_pars(model_spec)
  public_to_internal <- split(names(lookup), unname(lookup))
  public_to_internal <- public_to_internal[public_order]

  public_to_runtime <- stats::setNames(character(length(public_order)), public_order)
  public_defaults <- stats::setNames(numeric(length(public_order)), public_order)
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
    keys <- unique(vapply(target_info$profiles[targets], accumulatr_profile_key, character(1)))
    if (length(keys) != 1L) {
      stop(
        "AccumulatR public parameter '", nm,
        "' maps to internal parameters with incompatible transform/bound profiles",
        call. = FALSE
      )
    }
    public_to_runtime[[nm]] <- families[[1]]
    public_defaults[[nm]] <- target_info$profiles[[targets[[1]]]]$default
    public_profiles[[nm]] <- target_info$profiles[[targets[[1]]]]
  }

  runtime_families <- unique(unname(target_info$internal_to_runtime))
  runtime_defaults <- stats::setNames(numeric(length(runtime_families)), runtime_families)
  runtime_funcs <- stats::setNames(character(length(runtime_families)), runtime_families)
  runtime_lower <- stats::setNames(numeric(length(runtime_families)), runtime_families)
  runtime_upper <- stats::setNames(numeric(length(runtime_families)), runtime_families)
  runtime_min <- stats::setNames(numeric(length(runtime_families)), runtime_families)
  runtime_max <- stats::setNames(numeric(length(runtime_families)), runtime_families)
  runtime_exc <- stats::setNames(rep(NA_real_, length(runtime_families)), runtime_families)

  for (fam in runtime_families) {
    targets <- names(target_info$internal_to_runtime)[target_info$internal_to_runtime == fam]
    keys <- unique(vapply(target_info$profiles[targets], accumulatr_profile_key, character(1)))
    if (length(keys) != 1L) {
      publics <- unique(unname(lookup[targets]))
      stop(
        "AccumulatR runtime family '", fam,
        "' has incompatible internal parameter profiles across public parameters: ",
        paste(publics, collapse = ", "),
        call. = FALSE
      )
    }
    profile <- target_info$profiles[[targets[[1]]]]
    runtime_defaults[[fam]] <- profile$default
    runtime_funcs[[fam]] <- profile$func
    runtime_lower[[fam]] <- profile$lower
    runtime_upper[[fam]] <- profile$upper
    runtime_min[[fam]] <- profile$bound_min
    runtime_max[[fam]] <- profile$bound_max
    runtime_exc[[fam]] <- profile$exception
  }

  runtime_bound <- list(
    minmax = rbind(runtime_min, runtime_max),
    exception = runtime_exc[!is.na(runtime_exc)]
  )

  list(
    public_names = public_order,
    public_defaults = public_defaults,
    public_to_internal = public_to_internal,
    public_to_runtime = public_to_runtime,
    public_profiles = public_profiles,
    internal_to_runtime = target_info$internal_to_runtime,
    internal_owner = target_info$owners,
    runtime_defaults = runtime_defaults,
    runtime_transform = list(
      func = runtime_funcs,
      lower = runtime_lower,
      upper = runtime_upper
    ),
    runtime_bound = runtime_bound,
    runtime_to_public = split(names(public_to_runtime), unname(public_to_runtime))
  )
}

AccumulatR_compile_designs <- function(dadm, public_designs, bridge) {
  n_rows <- nrow(dadm)
  row_acc <- if ("accumulator" %in% names(dadm)) as.character(dadm$accumulator) else rep(NA_character_, n_rows)
  public_names <- bridge$public_names

  ownership <- setNames(vector("list", length(public_names)), public_names)
  for (nm in public_names) {
    targets <- bridge$public_to_internal[[nm]]
    owners <- unique(unname(bridge$internal_owner[targets]))
    owners <- owners[!is.na(owners)]
    if (length(owners) == 0L) {
      ownership[[nm]] <- rep(TRUE, n_rows)
    } else {
      ownership[[nm]] <- row_acc %in% owners
    }
  }

  expand_design <- function(dm) {
    exp_idx <- attr(dm, "expand")
    if (is.null(exp_idx)) exp_idx <- seq_len(nrow(dm))
    dm[exp_idx, , drop = FALSE]
  }

  compress_design <- function(dm) {
    keys <- apply(dm, 1, paste, collapse = "\r")
    keep <- !duplicated(keys)
    out <- dm[keep, , drop = FALSE]
    attr(out, "expand") <- match(keys, keys[keep])
    out
  }

  runtime_designs <- list()
  family_order <- unique(unname(bridge$public_to_runtime))
  for (fam in family_order) {
    fam_public <- bridge$runtime_to_public[[fam]]
    parts <- vector("list", length(fam_public))
    for (i in seq_along(fam_public)) {
      public_name <- fam_public[[i]]
      dm <- public_designs[[public_name]]
      if (is.null(dm)) next
      full_dm <- expand_design(dm)
      mask <- as.numeric(ownership[[public_name]])
      parts[[i]] <- full_dm * mask
    }
    parts <- Filter(Negate(is.null), parts)
    if (length(parts) == 0L) next
    runtime_designs[[fam]] <- compress_design(do.call(cbind, parts))
  }

  passthrough <- setdiff(names(public_designs), public_names)
  for (nm in passthrough) {
    runtime_designs[[nm]] <- public_designs[[nm]]
  }

  bridge$row_ownership <- ownership
  bridge$public_designs <- public_designs
  bridge$runtime_design_names <- names(runtime_designs)
  bridge$runtime_family_sharing <- lapply(bridge$runtime_to_public, unique)
  list(designs = runtime_designs, bridge = bridge)
}

AccumulatR_subset_bridge <- function(bridge, row_idx, public_designs = NULL) {
  if (is.null(bridge)) return(NULL)
  bridge$row_ownership <- lapply(bridge$row_ownership, function(x) x[row_idx])
  if (!is.null(public_designs)) {
    bridge$public_designs <- public_designs
  } else {
    bridge$public_designs <- NULL
  }
  bridge
}

AccumulatR_bridge_reporting_designs <- function(dadm) {
  bridge <- attr(dadm, "accumulatr_bridge")
  if (!is.null(bridge) && !is.null(bridge$public_designs)) {
    return(bridge$public_designs)
  }
  attr(dadm, "designs")
}

AccumulatR_project_public_pars <- function(runtime_pars, dadm, model_list) {
  bridge <- attr(dadm, "accumulatr_bridge")
  if (is.null(bridge)) {
    return(runtime_pars)
  }

  public_names <- bridge$public_names
  out <- matrix(NA_real_, nrow = nrow(runtime_pars), ncol = length(public_names))
  colnames(out) <- public_names
  for (nm in public_names) {
    fam <- bridge$public_to_runtime[[nm]]
    if (!fam %in% colnames(runtime_pars)) next
    mask <- bridge$row_ownership[[nm]]
    out[mask, nm] <- runtime_pars[mask, fam]
  }
  ok <- attr(runtime_pars, "ok")
  out <- as.matrix(out)
  attr(out, "ok") <- ok
  out
}

AccumulatR_validate_trend_safety <- function(trend, model) {
  if (is.null(trend) || is.null(model) || !identical(model$type, "AccumulatR")) {
    return(invisible(NULL))
  }
  bridge <- model$accumulatr_bridge
  if (is.null(bridge)) return(invisible(NULL))

  tnames <- names(trend)
  for (i in seq_along(trend)) {
    phase <- trend[[i]]$phase
    if (identical(phase, "premap")) next
    param_name <- tnames[[i]]
    if (!param_name %in% bridge$public_names) next
    fam <- bridge$public_to_runtime[[param_name]]
    owners <- unique(bridge$runtime_to_public[[fam]])
    if (length(owners) > 1L) {
      stop(
        "AccumulatR ", phase, " trend on '", param_name,
        "' is ambiguous because runtime family '", fam,
        "' is shared by public parameters: ",
        paste(owners, collapse = ", "),
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}

AccumulatR_model <- function(model_spec){
  bridge <- AccumulatR_bridge_spec(model_spec)
  public_transform <- accumulatr_identity_transform(bridge$public_names)
  model_list <- list(
    c_name = "AccumulatR",
    type = "AccumulatR",
    p_types = bridge$public_defaults,
    runtime_p_types = bridge$runtime_defaults,
    transform = public_transform,
    runtime_transform = bridge$runtime_transform,
    pre_transform = accumulatr_empty_transform(),
    Ttransform = function(pars,dadm) pars,
    bound = accumulatr_open_bound(bridge$public_names),
    runtime_bound = bridge$runtime_bound,
    rfun = function(model_spec, trial_df, pars) AccumulatR::simulate(model_spec, pars, trial_df = trial_df, keep_component = TRUE),
    spec = model_spec,
    accumulatr_bridge = bridge
  )
  return(function() {return(model_list)})
}

simulate_AccumulatR <- function(model_spec, component, pars){

}


AccumulatR_expand_data <- function(model_spec, data){
  if(is.null(model_spec)) stop("model_specification needs to be made for AccumulatR models")
  if(nrow(model_spec$component) > 1){
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
  ctx <- make_context(model_spec)
  attr(dadm, "AccumulatR_context") <- ctx
  return(dadm)
}

AccumulatR_check_context <- function(emc){
  if(!emc[[1]]$model()$type == "AccumulatR") return(emc)
  model_spec <- emc[[1]]$model()$spec
  for(i in 1:length(emc)){
    for(j in 1:length(emc[[i]]$data)){
      dat <- emc[[i]]$data[[j]] # Looping over chains, and subjects
      if(is.data.frame(dat)){
        ctx <- make_context(model_spec)
        attr(dat, "AccumulatR_context") <- ctx
        emc[[i]]$data[[j]] <- dat
      }
    }
  }
  return(emc)
}
