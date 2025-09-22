#' Create a stop-signal delay generator
#'
#' @description
#' `make_ssd()` returns a function that can be supplied to the
#' `functions` argument of [make_data()] to generate the `SSD` column for
#' stop-signal simulations. It supports both fixed SSD values and staircase
#' procedures. Staircases can be configured globally or separately for any
#' combination of grouping factors (e.g., subjects, stimulus categories).
#' By default, a separate staircase is maintained for each subject (using the
#' `subjects` column in the data) so that every participant starts at
#' `SSD0`.
#'
#' @param values Numeric vector of fixed SSD values. Supply this when
#'   `staircase = FALSE`. If `NULL` (the default) the generator will create
#'   staircase trials.
#' @param p Numeric vector of probabilities corresponding to `values`. When
#'   its length equals `length(values)` the entries are interpreted as the
#'   proportion of trials allocated to each SSD level (they need not sum to
#'   1, the remainder is treated as go trials). If `length(p) == 1` the value
#'   specifies the overall stop-trial probability and SSD levels are sampled
#'   uniformly. When `NULL`, the default overall stop probability is taken
#'   from `p_stop`.
#' @param staircase Either a logical flag (default `TRUE`) indicating whether
#'   to use a staircase, or a list describing staircase settings. A single
#'   staircase specification is a list with elements `SSD0`, `stairstep`,
#'   `stairmin`, `stairmax`, and `p`. A list of such lists can be supplied to
#'   configure multiple staircases.
#' @param SSD0,stairstep,stairmin,stairmax Default staircase parameters used
#'   when `staircase = TRUE` or when a group-specific specification omits the
#'   corresponding field.
#' @param p_stop Default probability of a stop trial for staircase modes (and
#'   for fixed SSDs when `p` is omitted). Must lie between 0 and 1.
#' @param factors Character vector of column names that define separate
#'   staircases. When `NULL`, all trials share the same staircase.
#' @param formula Optional formula specifying the grouping structure, e.g.
#'   `SSD ~ S | subjects`. The right-hand side (and, if present, the
#'   conditioning part beyond the `|`) are converted into grouping factor
#'   names. The `factors` argument and the variables implied by the formula
#'   are combined.
#'
#' @return A function of a data frame that returns a numeric vector of SSDs.
#'   The function carries class `emc_ssd_function` so that [make_data()] can
#'   attach the generated staircase specifications to the simulated data.
#'
#' @examples
#' # Fixed SSDs sampled on 25% of trials
#' ssd_fixed <- make_ssd(values = c(.26, .35, .46), p = rep(.25 / 3, 3), staircase = FALSE)
#'
#' # Staircase with default parameters operating per subject and stimulus
#' ssd_stair <- make_ssd(factors = c("subjects", "S"))
#'
#' @export
make_ssd <- function(values = NULL,
                     p = NULL,
                     staircase = TRUE,
                     SSD0 = .25,
                     stairstep = .05,
                     stairmin = 0,
                     stairmax = Inf,
                     p_stop = 0.25,
                     factors = NULL,
                     formula = NULL) {

  if (!is.null(values)) {
    if (!is.numeric(values)) {
      stop("`values` must be numeric when supplied.")
    }
    if (isTRUE(staircase)) {
      staircase <- FALSE
    }
  }

  if (isFALSE(staircase) && is.null(values)) {
    stop("`values` must be supplied when `staircase = FALSE`.")
  }

  if (!is.logical(staircase) && !is.list(staircase)) {
    stop("`staircase` must be TRUE/FALSE or a list of staircase specifications.")
  }

  if (!is.numeric(p_stop) || length(p_stop) != 1 || is.na(p_stop) || p_stop < 0 || p_stop > 1) {
    stop("`p_stop` must be a single numeric value between 0 and 1.")
  }

  group_cols <- unique(c(
    factors %||% character(),
    ssd_parse_formula(formula)
  ))
  base_spec <- list(
    SSD0 = SSD0,
    stairstep = stairstep,
    stairmin = stairmin,
    stairmax = stairmax,
    p = p_stop
  )

  if (!isFALSE(staircase) && is.list(staircase) && !length(group_cols) && is.null(values)) {
    extra_entries <- setdiff(names(staircase), names(base_spec))
    if (length(extra_entries) > 1) {
      stop("When supplying multiple staircase specifications you must specify `factors` or `formula` to identify groups.")
    }
  }

  assign_fun <- function(d) {
    if (!is.data.frame(d)) {
      stop("`make_ssd()` generated functions expect a data frame.")
    }

    n_trial <- nrow(d)
    if (!n_trial) {
      return(numeric(0))
    }

    SSD <- rep(Inf, n_trial)
    group_cols_local <- group_cols
    if ("subjects" %in% names(d)) {
      group_cols_local <- c("subjects", setdiff(group_cols_local, "subjects"))
    }

    if (length(group_cols_local)) {
      missing_cols <- setdiff(group_cols_local, names(d))
      if (length(missing_cols)) {
        stop("Grouping variables not found in data: ", paste(missing_cols, collapse = ", "))
      }
      group_data <- d[group_cols_local]
      group_id <- interaction(group_data, drop = TRUE, sep = "::")
    } else {
      group_id <- factor(rep(".all", n_trial), levels = ".all")
    }

    if (isFALSE(staircase)) {
      assign_fixed_ssd(SSD, values, p, p_stop)
    } else {
      specs <- build_staircase_specs(group_id, d, staircase, base_spec, group_cols_local)
      assign_staircase_ssd(SSD, group_id, d, specs)
    }
  }

  class(assign_fun) <- c("emc_ssd_function", class(assign_fun))
  assign_fun
}


assign_fixed_ssd <- function(SSD, values, p, p_stop) {
  n <- length(SSD)
  if (!length(values)) {
    return(SSD)
  }

  prob <- prepare_value_probabilities(values, p, p_stop)
  stop_trials <- stats::rbinom(n, 1, prob$total_prob) == 1
  if (any(stop_trials)) {
    SSD[stop_trials] <- sample(values, size = sum(stop_trials), replace = TRUE, prob = prob$weights)
  }
  SSD
}


assign_staircase_ssd <- function(SSD, group_id, data, specs) {
  n <- length(SSD)
  if (!length(specs)) {
    return(SSD)
  }

  stop_meta <- list(
    specs = specs,
    group_id = NULL,
    data = data[0, , drop = FALSE],
    group_cols = attr(specs, "group_cols")
  )

  is_stop_all <- rep(FALSE, n)
  has_stop <- FALSE
  for (lvl in names(specs)) {
    idx <- which(group_id == lvl)
    if (!length(idx)) {
      next
    }
    spec <- specs[[lvl]]
    p_group <- spec$p
    if (is.null(p_group)) {
      p_group <- attr(specs, "base_spec")$p
    }
    p_group <- max(min(p_group, 1), 0)
    is_stop <- stats::rbinom(length(idx), 1, p_group) == 1
    if (!any(is_stop)) {
      next
    }
    has_stop <- TRUE
    SSD[idx[is_stop]] <- NA_real_
    is_stop_all[idx[is_stop]] <- TRUE
  }

  if (has_stop) {
    stop_meta$group_id <- factor(as.character(group_id[is_stop_all]), levels = names(specs))
    stop_meta$data <- data[is_stop_all, , drop = FALSE]
    stair_fun <- attr(specs, "staircase_function")
    if (!is.null(stair_fun)) {
      attr(stop_meta, "staircase_function") <- stair_fun
    }
    class(stop_meta) <- c("emc_staircase", "list")
    attr(SSD, "emc_ssd") <- list(staircase = stop_meta)
  }

  SSD
}


prepare_value_probabilities <- function(values, p, p_stop) {
  if (!length(values)) {
    stop("`values` must contain at least one SSD level for fixed sampling.")
  }

  if (is.null(p)) {
    total <- p_stop
    weights <- rep(1 / length(values), length(values))
  } else if (length(p) == length(values)) {
    if (any(p < 0)) {
      stop("Probabilities in `p` must be non-negative.")
    }
    total <- sum(p)
    if (total > 1 + sqrt(.Machine$double.eps)) {
      stop("Sum of `p` cannot exceed 1.")
    }
    if (total > 0) {
      weights <- p / total
    } else {
      weights <- rep(1 / length(values), length(values))
    }
  } else if (length(p) == 1) {
    total <- p
    weights <- rep(1 / length(values), length(values))
  } else {
    stop("`p` must have length 1 or length equal to `values`.")
  }

  total <- max(min(total, 1), 0)
  list(total_prob = total, weights = weights)
}


build_staircase_specs <- function(group_id, data, staircase, base_spec, group_cols) {
  levels_id <- levels(group_id)
  specs <- stats::setNames(vector("list", length(levels_id)), levels_id)

  stair_fun <- attr(staircase, "staircase_function")

  staircase_list <- staircase
  if (is.logical(staircase_list)) {
    staircase_list <- list()
  }

  base_overrides <- staircase_list[names(staircase_list) %in% names(base_spec)]
  if (length(base_overrides)) {
    base_spec <- utils::modifyList(base_spec, base_overrides)
    staircase_list <- staircase_list[setdiff(names(staircase_list), names(base_overrides))]
  }

  for (lvl in levels_id) {
    spec <- base_spec
    matched_override <- FALSE
    level_vals <- strsplit(lvl, "::", fixed = TRUE)[[1]]
    if (length(level_vals) < length(group_cols)) {
      level_vals <- c(level_vals, rep("", length(group_cols) - length(level_vals)))
    }

    res <- resolve_staircase_overrides(spec, staircase_list, group_cols, level_vals)
    spec <- res$spec
    matched_override <- matched_override || res$matched

    for (col in group_cols) {
      if (is.list(spec[[col]])) spec[[col]] <- NULL
    }

    specs[[lvl]] <- spec
  }

  attr(specs, "base_spec") <- base_spec
  attr(specs, "group_cols") <- group_cols
  if (!is.null(stair_fun)) {
    attr(specs, "staircase_function") <- stair_fun
  }
  specs
}


ssd_parse_formula <- function(formula) {
  if (is.null(formula)) {
    return(character())
  }
  if (inherits(formula, "list")) {
    formula <- formula[[1]]
  }
  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula or list containing a formula.")
  }
  if (length(formula) != 3L) {
    stop("`formula` must be of the form SSD ~ factors | groups.")
  }

  rhs <- formula[[3]]
  extract_vars <- function(expr) {
    if (is.null(expr)) return(character())
    all.vars(expr)
  }

  if (is.call(rhs) && identical(rhs[[1]], as.name("|"))) {
    vars <- c(extract_vars(rhs[[2]]), extract_vars(rhs[[3]]))
  } else {
    vars <- extract_vars(rhs)
  }
  unique(vars)
}


`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}


resolve_staircase_overrides <- function(spec, overrides, cols, vals) {
  matched <- FALSE
  if (is.null(overrides) || !is.list(overrides)) {
    return(list(spec = spec, matched = matched))
  }

  full_key <- paste(vals, collapse = "::")
  if (nzchar(full_key) && !is.null(overrides[[full_key]]) && is.list(overrides[[full_key]])) {
    spec <- utils::modifyList(spec, overrides[[full_key]])
    matched <- TRUE
  }

  if (!length(cols)) {
    return(list(spec = spec, matched = matched))
  }

  current_lists <- list(overrides)
  for (i in seq_along(cols)) {
    col <- cols[[i]]
    val <- vals[[i]]
    next_lists <- list()
    for (lst in current_lists) {
      if (!is.list(lst)) next
      by_col <- lst[[col]]
      if (is.list(by_col) && !is.null(by_col[[val]])) {
        next_lists <- c(next_lists, list(by_col[[val]]))
      }
      by_val <- lst[[val]]
      if (is.list(by_val)) {
        next_lists <- c(next_lists, list(by_val))
      }
    }
    if (!length(next_lists)) {
      next_lists <- current_lists
    }
    current_lists <- next_lists
  }

  if (length(current_lists)) {
    for (lst in current_lists) {
      if (is.list(lst)) {
        spec <- utils::modifyList(spec, lst)
        matched <- TRUE
      }
    }
  }

  list(spec = spec, matched = matched)
}
