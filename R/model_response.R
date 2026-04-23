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
