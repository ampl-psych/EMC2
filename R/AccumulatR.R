.accumulatr_profile <- function(model, internal) {
  prep <- model$prep
  if (internal %in% names(prep$shared_triggers)) {
    return(c(default = -Inf, actual = 0, lower = 0, upper = 1,
             bound_min = 0, bound_max = 1, exception = 0, transform = "pnorm"))
  }

  weight_names <- vapply(
    prep$components$attrs,
    function(x) x$weight_param %||% NA_character_,
    character(1)
  )
  if (internal %in% weight_names) {
    weight <- 1 / length(prep$components$ids)
    return(c(default = qnorm(weight), actual = weight, lower = 0, upper = 1,
             bound_min = 0, bound_max = 1, exception = 0, transform = "pnorm"))
  }

  pieces <- strsplit(internal, ".", fixed = TRUE)[[1]]
  accumulator <- prep$accumulators[[pieces[[1]]]]
  parameter <- pieces[[2]]
  if (identical(parameter, "t0")) {
    return(c(default = -Inf, actual = 0, lower = 0, upper = Inf,
             bound_min = 0.05, bound_max = Inf, exception = 0, transform = "exp"))
  }

  positive <- parameter %in% c("s", "sigma", "tau", "shape", "rate", "B", "A", "sv") ||
    identical(tolower(accumulator$dist), "rdm") && identical(parameter, "v")
  if (positive) {
    return(c(default = 0, actual = 1, lower = 0, upper = Inf,
             bound_min = 0, bound_max = Inf,
             exception = if (parameter %in% c("A", "v")) 0 else NA,
             transform = "exp"))
  }
  c(default = 0, actual = 0, lower = -Inf, upper = Inf,
    bound_min = -Inf, bound_max = Inf, exception = NA, transform = "identity")
}

.accumulatr_bridge <- function(model) {
  lookup <- model$prep$parameter_lookup
  public <- unique(unname(lookup))
  profiles <- lapply(public, function(name) {
    mapped <- lapply(names(lookup)[lookup == name], function(internal) {
      .accumulatr_profile(model, internal)
    })
    keys <- vapply(mapped, paste, collapse = "\r", character(1))
    if (length(unique(keys)) != 1L) {
      stop("AccumulatR parameter '", name, "' combines incompatible parameter types")
    }
    mapped[[1]]
  })
  names(profiles) <- public

  field <- function(name, numeric = TRUE) {
    out <- vapply(profiles, `[[`, character(1), name)
    setNames(if (numeric) as.numeric(out) else out, public)
  }
  p_types <- field("default")
  names(p_types) <- public
  minmax <- rbind(field("bound_min"), field("bound_max"))
  colnames(minmax) <- public
  exception <- field("exception")
  names(exception) <- public
  exception <- exception[!is.na(exception)]

  list(
    p_types = p_types,
    transform = list(
      func = field("transform", FALSE),
      lower = field("lower"),
      upper = field("upper")
    ),
    bound = list(minmax = minmax, exception = exception),
    profiles = profiles
  )
}

.accumulatr_expand_rows <- function(data, model) {
  if (is.null(data$trials)) data$trials <- seq_len(nrow(data))
  accumulators <- names(model$prep$accumulators)
  rows <- rep(seq_len(nrow(data)), each = length(accumulators))
  out <- data[rows, , drop = FALSE]
  out$trials <- rep(data$trials, each = length(accumulators))
  out$racer <- factor(
    rep(accumulators, times = nrow(data)),
    levels = accumulators
  )
  out$lR <- out$racer
  rownames(out) <- NULL
  out
}

.accumulatr_runtime_columns <- function(model) {
  prep <- model$prep
  slots <- vapply(names(prep$accumulators), function(id) {
    internals <- names(prep$parameter_lookup)
    sum(startsWith(internals, paste0(id, "."))) - 1L
  }, integer(1))
  weights <- vapply(
    prep$components$attrs,
    function(x) x$weight_param %||% NA_character_,
    character(1)
  )
  c("q", "t0", paste0("p", seq_len(max(slots))), unname(weights[!is.na(weights)]))
}

.accumulatr_runtime_recipe <- function(model, data, bridge) {
  prep <- model$prep
  lookup <- prep$parameter_lookup
  columns <- .accumulatr_runtime_columns(model)
  defaults <- matrix(0, nrow(data), length(columns), dimnames = list(NULL, columns))
  sources <- matrix(NA_character_, nrow(data), length(columns), dimnames = list(NULL, columns))

  for (id in names(prep$accumulators)) {
    rows <- which(data$racer == id)
    accumulator <- prep$accumulators[[id]]
    internals <- names(lookup)[startsWith(names(lookup), paste0(id, "."))]
    distribution <- setdiff(internals, paste0(id, ".t0"))
    for (slot in seq_along(distribution)) {
      internal <- distribution[[slot]]
      public <- lookup[[internal]]
      defaults[rows, paste0("p", slot)] <- as.numeric(bridge$profiles[[public]][["actual"]])
      sources[rows, paste0("p", slot)] <- public
    }
    t0 <- paste0(id, ".t0")
    if (t0 %in% names(lookup)) {
      public <- lookup[[t0]]
      defaults[rows, "t0"] <- as.numeric(bridge$profiles[[public]][["actual"]])
      sources[rows, "t0"] <- public
    }
    trigger <- accumulator$shared_trigger_id
    if (!is.null(trigger)) {
      public <- lookup[[trigger]]
      defaults[rows, "q"] <- as.numeric(bridge$profiles[[public]][["actual"]])
      sources[rows, "q"] <- public
    }
    if (!is.null(data$component) && length(accumulator$components)) {
      inactive <- !is.na(data$component[rows]) &
        !data$component[rows] %in% accumulator$components
      # q belongs to the shared trigger, even when its representative is inactive.
      sources[rows[inactive], setdiff(columns, "q")] <- NA_character_
    }
  }

  weights <- intersect(columns, names(lookup))
  for (weight in weights) {
    public <- lookup[[weight]]
    defaults[, weight] <- as.numeric(bridge$profiles[[public]][["actual"]])
    sources[, weight] <- public
  }
  list(defaults = defaults, source_names = sources)
}

.accumulatr_runtime_parameters <- function(parameters, recipe) {
  runtime <- recipe$defaults
  for (column in seq_len(ncol(runtime))) {
    source <- recipe$source_names[, column]
    index <- match(source, colnames(parameters))
    rows <- which(!is.na(index))
    if (length(rows)) {
      runtime[cbind(rows, rep(column, length(rows)))] <-
        parameters[cbind(rows, index[rows])]
    }
  }
  runtime
}

.accumulatr_prepare_subject <- function(data, model, bridge) {
  prepared <- AccumulatR::prepare_data(model, data)
  prepared$lR <- factor(prepared$racer, levels = names(model$prep$accumulators))
  attr(prepared, "AccumulatR_bridge") <- list(
    bridge = .accumulatr_runtime_recipe(model, prepared, bridge)
  )
  prepared
}

#' Use an AccumulatR model in EMC2
#'
#' @param model A finalized AccumulatR model.
#' @return An EMC2 model function suitable for `design()`.
#' @export
AccumulatR_model <- function(model) {
  if (!inherits(model, "model_structure")) {
    stop("AccumulatR_model() requires a finalized AccumulatR model")
  }
  bridge <- .accumulatr_bridge(model)
  context <- new.env(parent = emptyenv())
  native_context <- function() {
    if (!accumulatr_context_valid(context$native)) {
      context$native <- AccumulatR::make_context(model)$cpp
    }
    context$native
  }
  simulate <- function(data, parameters) {
    recipe <- .accumulatr_runtime_recipe(model, data, bridge)
    runtime <- .accumulatr_runtime_parameters(parameters, recipe)
    trial_data <- data[, intersect(c("trials", "racer", "component", "onset"), names(data)), drop = FALSE]
    n_accumulators <- length(model$prep$accumulators)
    trial_data$trials <- rep(
      seq_len(nrow(trial_data) / n_accumulators),
      each = n_accumulators
    )
    AccumulatR::simulate(
      model,
      runtime,
      trial_df = trial_data,
      keep_component = TRUE
    )
  }
  ranks <- model$prep$observation$n_outcomes
  ranked_columns <- if (ranks > 1L) {
    unlist(lapply(2:ranks, function(rank) c(paste0("R", rank), paste0("rt", rank))))
  } else {
    character()
  }
  model_list <- list(
    c_name = "AccumulatR",
    type = "AccumulatR",
    p_types = bridge$p_types,
    transform = bridge$transform,
    pre_transform = NULL,
    Ttransform = function(parameters, data) parameters,
    bound = bridge$bound,
    rfun = simulate,
    spec = model,
    accumulatr_bridge = bridge,
    native_context = native_context,
    compression_columns = c(
      "component", "onset", "LT", "UT", "LC", "UC", "missingness",
      ranked_columns
    )
  )
  function() model_list
}
