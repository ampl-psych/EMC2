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
#' @param staircase_up Character vector listing the response labels (or `NA`)
#'   that should increase the staircase on the next stop trial. Defaults to
#'   `NA`, matching the classic rule where successful inhibition yields an
#'   `NA` response.
#' @param staircase_down Character vector listing response labels that should
#'   decrease the staircase. If `NULL` (the default), any stop-trial outcome
#'   not matched by `staircase_up` will decrease the staircase.
#'
#' @param UC Optional response deadline (seconds) used by the trial-by-trial
#'   staircase (see Details): a response slower than `UC` is unobserved in the
#'   experiment, so the staircase treats it as a non-response. Ignored when
#'   `NULL`; when the data carry a `UC` column that is used instead.
#'
#' @details
#' The generator can be used in two places:
#'
#' * `make_data(functions = list(SSD = make_ssd(...)))`: SSDs are assigned
#'   once for the whole data set and a staircase is run inside the model's
#'   random function. This is the vectorised path and requires that no
#'   parameter depends on SSD.
#' * `design(functions = list(SSD = make_ssd(...)))`: the generator is a
#'   pre-trial design function. [make_data()] then simulates trial by trial
#'   (the `conditional_on_data = FALSE` path): on each trial the generator
#'   steps the staircase from the previous stop trial's outcome and assigns
#'   the current trial's SSD before the parameters for that trial are
#'   computed. This is required when a trend makes a parameter depend on SSD
#'   (e.g. a `slin_incr` kernel of `SSD` on `muS`). With fitted data the default (conditional)
#'   simulation keeps the observed SSDs, e.g. `predict(emc, conditional_on_data
#'   = TRUE)`; the unconditional simulation (the default of [predict()] when
#'   the design holds the generator) re-runs the staircase.
#'
#' @return A function of a data frame that returns a numeric vector of SSDs.
#'   The function carries class `emc_ssd_function` (and the `pretrial`
#'   attribute) so that [make_data()] can attach the generated staircase
#'   specifications to the simulated data or run the staircase trial by trial.
#'
#' @examples
#' # Fixed SSDs sampled on 25% of trials
#' ssd_fixed <- make_ssd(values = c(.26, .35, .46), p = rep(.25 / 3, 3), staircase = FALSE)
#'
#' # Staircase with default parameters operating per subject and stimulus
#' ssd_stair <- make_ssd(factors = c("subjects", "S"))
#'
#' # Complex paradigms with stop-triggered responses -----------------
#'
#' # By default `staircase_up = NA` and `staircase_down = NULL`, replicating
#' # the classic rule in which an `NA` response indicates successful inhibition
#' # (increase SSD) and any finite response indicates failure (decrease SSD).
#'
#' # You can override these rules to accommodate stop-triggered accumulators.
#' # For example, in a stop-change paradigm, the stop-triggered responses are
#' # recorded as "down" or "up", and the original go responses remain
#' # "left"/"right". The following applies the staircase rules accordingly:
#'
#' ssd_change <- make_ssd(
#'   staircase_up   = c("down", "up"),
#'   staircase_down = c("left", "right")
#' )
#'
#' # A mixture of NA (complete inhibition) and stop-triggered wins can be
#' # handled by listing both the NA marker and the relevant response labels:
#'
#' ssd_sel <- make_ssd(
#'   staircase_up   = c(NA, "ST_left", "ST_right"),
#'   staircase_down = c("left", "right")
#' )
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
                     formula = NULL,
                     staircase_up = NA,
                     staircase_down = NULL,
                     UC = NULL) {

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

  normalise_rule <- function(x) {
    if (is.null(x)) return(NULL)
    if (!length(x)) return(NULL)
    nas <- is.na(x)
    out <- as.character(x)
    out[nas] <- NA_character_
    out <- unique(out)
    if (!length(out)) NULL else out
  }

  up_rule <- normalise_rule(staircase_up)
  down_rule <- normalise_rule(staircase_down)
  if (!is.null(up_rule) && !is.null(down_rule)) {
    overlap <- intersect(up_rule[!is.na(up_rule)], down_rule[!is.na(down_rule)])
    if (length(overlap)) {
      stop("`staircase_up` and `staircase_down` cannot share response labels: ",
           paste(overlap, collapse = ", "))
    }
    if (any(is.na(up_rule)) && any(is.na(down_rule))) {
      stop("`staircase_up` and `staircase_down` cannot both match NA.")
    }
  }
  staircase_rules <- list(up = up_rule, down = down_rule)

  group_cols <- unique(c(
    factors %||% character(),
    ssd_parse_formula(formula)
  ))
  base_spec <- list(
    SSD0 = SSD0,
    stairstep = stairstep,
    stairmin = stairmin,
    stairmax = stairmax,
    p = p_stop,
    rules = staircase_rules
  )

  if (!isFALSE(staircase) && is.list(staircase) && !length(group_cols) && is.null(values)) {
    extra_entries <- setdiff(names(staircase), names(base_spec))
    if (length(extra_entries) > 1) {
      stop("When supplying multiple staircase specifications you must specify `factors` or `formula` to identify groups.")
    }
  }

  # state for the trial-by-trial mode (design function): the ladder per
  # subject::group, reset whenever a subject's first trial is seen
  ladder_env <- new.env(parent = emptyenv())
  ladder_env$ladder <- list()
  ladder_env$specs  <- list()

  group_cols_for <- function(d) {
    gc <- group_cols
    if ("subjects" %in% names(d)) gc <- c("subjects", setdiff(gc, "subjects"))
    missing_cols <- setdiff(gc, names(d))
    if (length(missing_cols)) {
      stop("Grouping variables not found in data: ", paste(missing_cols, collapse = ", "))
    }
    gc
  }

  group_id_for <- function(d, gc) {
    if (length(gc)) interaction(d[gc], drop = TRUE, sep = "::") else
      factor(rep(".all", nrow(d)), levels = ".all")
  }

  # ---- vectorised assignment (one value per trial) ----
  assign_vectorised <- function(d) {
    n_trial <- nrow(d)
    SSD <- rep(Inf, n_trial)
    gc <- group_cols_for(d)
    group_id <- group_id_for(d, gc)
    if (isFALSE(staircase)) {
      assign_fixed_ssd(SSD, values, p, p_stop)
    } else {
      specs <- build_staircase_specs(group_id, d, staircase, base_spec, gc, staircase_rules)
      assign_staircase_ssd(SSD, group_id, d, specs, staircase_rules)
    }
  }

  # ---- trial-by-trial mode: d is the context window handed to a pre-trial
  # design function (previous trial(s) with their simulated R/rt, then the
  # current trial); SSD already exists, NA on staircase stop trials ----
  assign_trialwise <- function(d) {
    SSD <- d$SSD
    if (!anyNA(SSD)) return(SSD)                          # nothing to resolve
    gc <- group_cols_for(d)
    trial <- if ("trials" %in% names(d)) d$trials else seq_len(nrow(d))
    is_cur <- trial == max(trial)
    has_prev <- any(!is_cur)
    if (!has_prev) {                       # a subject's first trial: fresh ladders
      ladder_env$ladder <- list()
    }
    gid_of <- function(rows) {
      if (!length(gc)) return(".all")
      paste(vapply(gc, function(g) as.character(d[[g]][rows][1]), character(1)), collapse = "::")
    }
    spec_of <- function(gid, rows) {
      if (is.null(ladder_env$specs[[gid]])) {
        sp <- build_staircase_specs(factor(gid, levels = gid), d[rows, , drop = FALSE],
                                    staircase, base_spec, gc, staircase_rules)
        spec <- sp[[gid]]
        if (is.null(spec$SSD0)) spec$SSD0 <- attr(sp, "base_spec")$SSD0
        ladder_env$specs[[gid]] <- spec
      }
      ladder_env$specs[[gid]]
    }
    # step the ladder from the immediately preceding trial if it was a stop trial
    if (has_prev && !isFALSE(staircase)) {
      pr <- which(trial == max(trial[!is_cur]))
      if (is.finite(SSD[pr[1]])) {
        gid <- gid_of(pr)
        spec <- spec_of(gid, pr)
        R_i  <- d$R[pr[1]]
        rt_i <- if ("rt" %in% names(d)) d$rt[pr[1]] else NA_real_
        dl <- if ("UC" %in% names(d)) d$UC[pr[1]] else UC
        late <- !is.null(dl) && is.finite(dl) && !is.na(rt_i) && rt_i > dl
        label <- if (is.na(R_i)) NA_character_ else as.character(R_i)
        step_dir <- staircase_step_dir(spec, label, stopped = is.na(R_i), late = late)
        ladder_env$ladder[[gid]] <- staircase_next_ssd(SSD[pr[1]], step_dir, spec)
      }
    }
    # the current trial: decide stop vs go, then the SSD (ladder or fixed set)
    cur <- which(is_cur)
    if (anyNA(SSD[cur])) {
      if (isFALSE(staircase)) {
        pv <- prepare_value_probabilities(values, p, p_stop)
        SSD[cur] <- if (stats::runif(1) < pv$total_prob)
          sample(values, 1, prob = pv$weights) else Inf
      } else {
        gid <- gid_of(cur)
        spec <- spec_of(gid, cur)
        p_group <- spec$p
        if (is.null(p_group)) p_group <- base_spec$p
        if (stats::runif(1) < max(min(p_group, 1), 0)) {
          if (is.null(ladder_env$ladder[[gid]])) ladder_env$ladder[[gid]] <- spec$SSD0
          ssd_now <- staircase_clamp(ladder_env$ladder[[gid]], spec)
          ladder_env$ladder[[gid]] <- ssd_now
          SSD[cur] <- ssd_now
        } else {
          SSD[cur] <- Inf
        }
      }
    }
    SSD
  }

  assign_fun <- function(d) {
    if (!is.data.frame(d)) {
      stop("`make_ssd()` generated functions expect a data frame.")
    }
    if (!nrow(d)) {
      return(numeric(0))
    }
    if ("SSD" %in% names(d)) {
      return(assign_trialwise(d))
    }
    # vectorised: one draw per trial, replicated over accumulator rows if present
    if (all(c("subjects", "trials") %in% names(d)) && "lR" %in% names(d)) {
      key <- paste(d$subjects, d$trials)
      first <- !duplicated(key)
      out <- assign_vectorised(d[first, , drop = FALSE])
      meta <- attr(out, "emc_ssd")
      res <- as.numeric(out)[match(key, key[first])]
      if (!is.null(meta)) attr(res, "emc_ssd") <- meta
      return(res)
    }
    assign_vectorised(d)
  }

  attr(assign_fun, "pretrial") <- TRUE
  attr(assign_fun, "staircase") <- !isFALSE(staircase)
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


assign_staircase_ssd <- function(SSD, group_id, data, specs, rules) {
  n <- length(SSD)
  if (!length(specs)) {
    return(SSD)
  }

  stop_meta <- list(
    specs = specs,
    group_id = NULL,
    data = data[0, , drop = FALSE],
    group_cols = attr(specs, "group_cols"),
    rules = rules
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


build_staircase_specs <- function(group_id, data, staircase, base_spec, group_cols, rules) {
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

    spec$rules <- rules
    specs[[lvl]] <- spec
  }

  attr(specs, "base_spec") <- base_spec
  attr(specs, "group_cols") <- group_cols
  attr(specs, "rules") <- rules
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
