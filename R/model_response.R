model_type <- function(x) {
  if (is.function(x)) return(model_type(x()))
  if (is.list(x) && !is.null(x$type)) return(x$type)
  as.character(x)
}

.model_list <- function(x) {
  if (is.function(x)) return(x())
  x
}

is_ordered_response_type <- function(x) {
  identical(model_type(x), "ORDERED")
}

is_multinomial_response_type <- function(x) {
  identical(model_type(x), "MULTINOMIAL")
}

is_choice_accumulator_type <- function(x) {
  model_type(x) %in% c("ORDERED", "MULTINOMIAL", "RACE", "MT", "TC")
}

is_choice_only_model_type <- function(x) {
  model_type(x) %in% c("ORDERED", "MULTINOMIAL")
}

.by_trial_index <- function(trials) {
  split(seq_along(trials), trials)
}

model_sampled_by_default <- function(model) {
  model_list <- .model_list(model)
  sampled <- model_list$sampled_by_default
  if (is.null(sampled)) character(0) else sampled
}

model_prepare_design <- function(model, formula, constants, ...) {
  model_list <- .model_list(model)
  if (is.null(model_list$prepare_design)) {
    return(list(formula = formula, constants = constants))
  }
  model_list$prepare_design(formula = formula, constants = constants, ...)
}

model_prepare_dm <- function(model, p_name, form, da) {
  model_list <- .model_list(model)
  if (is.null(model_list$prepare_dm)) {
    return(list(form = form, da = da))
  }
  model_list$prepare_dm(p_name = p_name, form = form, da = da)
}

# Response columns beyond R and rt that a model observes (e.g. a rating RR).
# They are data, not covariates: never imputed, never predictors, always part
# of the compression key, and returned by the model's rfun.
model_extra_responses <- function(model) {
  model_list <- .model_list(model)
  out <- model_list$extra_responses
  if (is.null(out)) character(0) else out
}

# Model-specific checks of real data, run by design_model() before the data are
# augmented (not when data are being simulated).
model_check_data <- function(model, data) {
  model_list <- .model_list(model)
  if (!is.null(model_list$check_data)) model_list$check_data(data)
  invisible(TRUE)
}
