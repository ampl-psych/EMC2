# Neural-likelihood validation kit: can a registered network be used for
# inference? Generalised from the NLE project's evaluation scripts (Wuth,
# 2026: score_bias_draws.R, score_bias_race.R, test_posterior_shift.R, and
# the checks of METHODS.md section 6 and TRAPS.md in its handover package).
#
# Every check works on a *cell* (nn_cell()): one design whose sampled
# parameters have independent normal priors lying inside the network's
# training region to +-k sd, built identically for the network and for an
# analytic control. The control cell differs from the network's in exactly one
# respect, the likelihood: same design, constants, defaults, prior and
# simulator. (A control with different defaults once certified a broken cell
# for two weeks; TRAPS.md #18.) All likelihoods go through EMC2's own pipeline
# (design_model() + calc_ll_manager(), or the model's dfun/pfun exactly as
# log_likelihood_ddm/race combine them), so bounds, transforms and the min_ll
# floor are those a fit uses.

# ---------------------------------------------------------------------------
# Cells
# ---------------------------------------------------------------------------

#' A Validation Cell for a Neural-Network Likelihood
#'
#' Builds the design and prior on which the neural-likelihood checks
#' ([nn_in_box()], [nn_total_mass()], [nn_score_bias()],
#' [nn_posterior_shift()], [nn_sbc_cell()]) run: one design (by default a
#' single cell, every network input sampled with an intercept-only formula),
#' independent normal priors on the sampled scale that lie inside the
#' network's training region to `k` standard deviations, and the same design
#' and prior for an analytic *control* model. Those functions also accept a
#' model directly and then build the cell with the defaults below.
#'
#' @details
#'
#' **Why a control.** A check is only interpretable next to the same check
#' run on the exact likelihood: a score bias or a posterior shift of the
#' analytic model must come out at about zero, and whatever it does show
#' (finite-difference error, Monte Carlo error, the simulator, the design) is
#' not the network's fault. The control cell uses the control model in an
#' otherwise identical design (same formula, constants, defaults and prior),
#' so the two cells differ in the likelihood only. The control must take the
#' same parameters as the network's model, with the same defaults and
#' transforms; by default it is the analytic twin the network was registered
#' with (e.g. [DDM] for [DDMnn]).
#'
#' **Defaults.** Sampled are the network's inputs (and a parameter used by
#' `pre`), except parameters fixed at registration and a parameter named `s`
#' (the within-trial noise, a scaling parameter: it stays at the model
#' default, which identifies the model). For a race network the rate varies
#' over accumulators (`v ~ 0 + lR`); every other formula is intercept only.
#' A prior mean defaults to the model's default where that lies inside the
#' training region and to the centre of the region otherwise; a prior sd
#' defaults to 99% of the largest value that keeps mean +- `k` sd inside the
#' region.
#' Both are on the sampled scale.
#'
#' **Priors inside the box.** Parameter vectors outside the training region
#' are rejected by the network (likelihood floored), which truncates the
#' posterior to the region. For fitting that is a choice; for SBC it breaks
#' calibration, because a replicate whose true parameters lie outside the
#' region cannot be recovered. `nn_cell()` therefore refuses priors that
#' leave the region within `k` sd ([nn_in_box()]); [run_sbc()] makes the same
#' check for neural-likelihood designs.
#'
#' @param model A neural-likelihood model function ([DDMnn], [RDMnn], or one
#'   returned by [register_nn_model()]).
#' @param mean,sd Prior means and sds on the sampled scale, named by sampled
#'   parameter (e.g. `v_lRa`) or by model parameter (e.g. `v`, then used for
#'   every sampled parameter that is the value of `v` in some cell). Never
#'   matched by position. `sd` may be a single unnamed number. Parameters
#'   that are effects (differences between cells) need explicit values.
#' @param constants Constants on the sampled scale, as in [design()].
#' @param formula,factors,Rlevels,matchfun As in [design()]; `factors`
#'   excludes `subjects` (one subject). Defaults: see Details; `Rlevels`
#'   `c("lower", "upper")` for two-response networks, `c("a", "b")` for races.
#' @param control The analytic control: `NULL` (the model's registered twin),
#'   an EMC2 model function taking the same parameters, or `FALSE` (none; the
#'   checks that need a control then refuse).
#' @param k Half-width of the prior, in sds, that must lie inside the
#'   training region.
#' @return An `nn_cell`: a list with the model, the registration, the design
#'   and prior for the network (`design`, `prior`) and for the control
#'   (`control_design`, `control_prior`), the prior `mean` and `sd` per
#'   sampled parameter, `k`, and the box check (`box`).
#' @seealso [nn_in_box()], [nn_score_bias()], [nn_posterior_shift()],
#'   [nn_total_mass()], [nn_sbc_cell()], [register_nn_model()]
#' @examples
#' cell <- nn_cell(RDMnn, mean = c(v = log(1.5), B = log(1), A = log(.3), t0 = log(.2)),
#'                 sd = .1)
#' cell
#' @export
nn_cell <- function(model, mean = NULL, sd = NULL, constants = NULL, formula = NULL,
                    factors = NULL, Rlevels = NULL, matchfun = NULL, control = NULL, k = 4) {
  if (inherits(model, "nn_cell")) return(model)
  control_label <- if (is.null(control)) NULL else if (isFALSE(control)) "none" else
    paste(deparse(substitute(control)), collapse = "")
  mfun <- nn_model_fun(model)
  ml <- mfun()
  reg <- ml$nn
  what <- paste0("Neural likelihood '", reg$label, "': ")
  joint <- ml$type == "DDM"
  if (!is.numeric(k) || length(k) != 1L || !(k > 0)) stop("k must be a positive number")
  if (is.null(Rlevels)) Rlevels <- if (joint) c("lower", "upper") else c("a", "b")
  if (is.null(formula)) {
    samp <- setdiff(c(reg$pars, reg$pre_pars), c(names(reg$fixed), names(constants), "s"))
    formula <- lapply(samp, function(p)
      stats::as.formula(if (!joint && p == "v") "v ~ 0 + lR" else paste(p, "~ 1"), env = baseenv()))
  }
  des_args <- list(factors = c(list(subjects = 1), factors), Rlevels = Rlevels, formula = formula,
                   constants = constants, matchfun = matchfun, report_p_vector = FALSE)
  des <- suppressMessages(do.call(design, c(des_args, list(model = mfun))))
  sp <- names(sampled_pars(des))
  map <- attr(sampled_pars(des, doMap = TRUE), "map")
  owner <- nn_sampled_owner(map, sp)
  tr <- fill_transform(NULL, des$model)
  box_s <- function(p) {
    b <- nn_param_box(p, reg, ml)
    c(nn_to_sampled(b[1], tr$func[[p]], tr$lower[[p]], tr$upper[[p]]),
      nn_to_sampled(b[2], tr$func[[p]], tr$lower[[p]], tr$upper[[p]]))
  }
  def_mean <- vapply(unique(stats::na.omit(owner)), function(p) {
    d <- ml$p_types[[p]]; b <- box_s(p)
    if (is.finite(d) && d > b[1] && d < b[2]) d else if (all(is.finite(b))) sum(b) / 2 else NA_real_
  }, 0)
  mean <- nn_resolve_prior(mean, sp, owner, def_mean, "mean", what)
  def_sd <- vapply(sp, function(s) {
    if (is.na(owner[[s]])) return(NA_real_)
    b <- box_s(owner[[s]]); 0.99 * min(mean[[s]] - b[1], b[2] - mean[[s]]) / k
  }, 0)
  if (!is.null(sd) && is.null(names(sd)) && length(sd) == 1L) sd <- stats::setNames(rep(sd, length(sp)), sp)
  sd <- nn_resolve_prior(sd, sp, owner, def_sd, "sd", what, by_sampled = TRUE)
  if (any(!(sd > 0))) stop(what, "prior sds must be positive (", paste(sp[!(sd > 0)], collapse = ", "), ")")
  box <- nn_box_table(des, reg, mean, sd, k)
  nn_refuse_box(box, reg, k, "the prior")
  pri <- suppressMessages(prior(des, type = "single", pmean = mean, psd = sd))

  ctl <- nn_control_model(reg, control, ml)
  des_ctl <- pri_ctl <- NULL
  if (!is.null(ctl)) {
    des_ctl <- suppressMessages(do.call(design, c(des_args, list(model = ctl))))
    if (!identical(names(sampled_pars(des_ctl)), sp))
      stop(what, "the control's design samples ", paste(names(sampled_pars(des_ctl)), collapse = ", "),
           ", not ", paste(sp, collapse = ", "))
    pri_ctl <- suppressMessages(prior(des_ctl, type = "single", pmean = mean, psd = sd))
    if (is.null(control_label))
      control_label <- if (is.character(reg$twin)) reg$twin else "twin"
  }
  structure(list(model = mfun, reg = reg, design = des, prior = pri, control = ctl,
                 control_design = des_ctl, control_prior = pri_ctl,
                 control_label = if (is.null(ctl)) "none" else control_label,
                 mean = mean, sd = sd, k = k, box = box), class = "nn_cell")
}

#' @export
print.nn_cell <- function(x, ...) {
  cat("Validation cell for neural likelihood '", x$reg$label, "' (", x$reg$kind, "); control: ",
      x$control_label, "\n", sep = "")
  cat("Prior (sampled scale):\n")
  print(data.frame(mean = round(x$mean, 4), sd = round(x$sd, 4)))
  cat(sprintf("Natural scale at +-%g sd, against the training region:\n", x$k))
  print(nn_format_box(x$box), row.names = FALSE)
  invisible(x)
}

# The model function of a neural-likelihood model given as a function, a model
# list or an emc.design.
nn_model_fun <- function(x) {
  if (inherits(x, "emc.design")) x <- x$model
  f <- if (is.function(x)) x else if (is.list(x) && !is.null(x$nn)) { ml <- x; function() ml } else NULL
  if (is.null(f) || !inherits(f()$nn, "emc_nn"))
    stop("Not a neural-likelihood model; register one with register_nn_model() (or use DDMnn, RDMnn)")
  f
}

nn_as_cell <- function(x, ...) {
  if (inherits(x, "nn_cell")) {
    if (length(list(...))) warning("arguments for nn_cell() are ignored when a cell is given")
    return(x)
  }
  nn_cell(x, ...)
}

nn_need_control <- function(cell, what) {
  if (is.null(cell$control))
    stop("Neural likelihood '", cell$reg$label, "': ", what, " needs an analytic control (the same ",
         "check on an exact likelihood must return about zero, else the estimator is at fault); ",
         "build the cell with control = <an EMC2 model function taking the same parameters>")
  invisible(TRUE)
}

# NULL -> the registered twin; FALSE -> none; else an analytic EMC2 model
# function with the network model's parameters, defaults and transforms.
nn_control_model <- function(reg, control, ml) {
  if (isFALSE(control)) return(NULL)
  what <- paste0("Neural likelihood '", reg$label, "': ")
  if (is.null(control)) {
    control <- reg$twin
    if (is.null(control)) return(NULL)
    if (is.character(control)) control <- get(control, envir = asNamespace("EMC2"))
  }
  if (!is.function(control)) stop(what, "control must be an EMC2 model function, NULL (the twin) or FALSE")
  cl <- control()
  if (!is.null(cl$nn)) stop(what, "the control must be an analytic model, not a neural likelihood")
  if (!identical(cl$type, ml$type)) stop(what, "the control is a ", cl$type, " model; the network's is ", ml$type)
  if (!setequal(names(cl$p_types), names(ml$p_types)))
    stop(what, "the control must take the model's parameters (", paste(names(ml$p_types), collapse = ", "),
         "); it takes ", paste(names(cl$p_types), collapse = ", "))
  if (!isTRUE(all.equal(cl$p_types[names(ml$p_types)], ml$p_types)))
    stop(what, "the control's parameter defaults differ from the model's: the cells would differ in ",
         "more than the likelihood (an un-sampled parameter would mean different things)")
  f1 <- fill_transform(NULL, control); f2 <- fill_transform(NULL, function() ml)
  if (!identical(f1$func[names(f2$func)], f2$func))
    stop(what, "the control's parameter transforms differ from the model's")
  if (is.null(cl$c_name) && !is.function(cl$dfun))
    stop(what, "the control has no likelihood (no c_name and no dfun)")
  control
}

# For each sampled parameter, the model parameter it is the value of in some
# design cell (a cell where it is the only non-zero coefficient, equal to 1);
# NA for an effect (a difference between cells).
nn_sampled_owner <- function(map, sp) {
  own <- stats::setNames(rep(NA_character_, length(sp)), sp)
  for (p in names(map)) {
    X <- nn_map_matrix(map[[p]])
    for (s in intersect(colnames(X), sp)) {
      alone <- X[, s] == 1 & rowSums(X != 0) == 1
      if (any(alone) && is.na(own[[s]])) own[[s]] <- p
    }
  }
  own
}

# numeric design-matrix columns of a map entry (drops data columns)
nn_map_matrix <- function(m) {
  m <- as.data.frame(m)
  keep <- vapply(m, function(z) is.numeric(z) && !is.logical(z), TRUE)
  as.matrix(m[, keep, drop = FALSE])
}

nn_resolve_prior <- function(given, sp, owner, default, what_arg, what, by_sampled = FALSE) {
  if (!is.null(given)) {
    if (!is.numeric(given) || is.null(names(given)) || any(!nzchar(names(given))))
      stop(what, "`", what_arg, "` must be a named numeric vector (sampled or model parameter names); ",
           "priors are never matched by position")
    bad <- setdiff(names(given), c(sp, stats::na.omit(owner)))
    if (length(bad)) stop(what, "`", what_arg, "` names ", paste(bad, collapse = ", "),
                          ", which are neither sampled parameters (", paste(sp, collapse = ", "),
                          ") nor model parameters they stand for")
  }
  out <- stats::setNames(rep(NA_real_, length(sp)), sp)
  for (s in sp) {
    o <- owner[[s]]
    out[[s]] <- if (s %in% names(given)) given[[s]] else
      if (!is.na(o) && o %in% names(given)) given[[o]] else
        if (by_sampled) default[[s]] else if (!is.na(o)) default[[o]] else NA_real_
  }
  miss <- sp[!is.finite(out)]
  if (length(miss))
    stop(what, "give the prior ", what_arg, " of ", paste(miss, collapse = ", "), " (an effect, or a parameter ",
         "without a finite training region, has no default)")
  out
}

# Natural-scale region of a model parameter: the training box for a network
# input, else the model's bound.
nn_param_box <- function(p, reg, ml) {
  if (p %in% reg$pars) return(c(reg$lower[[p]], reg$upper[[p]]))
  mm <- ml$bound$minmax
  if (!is.null(mm) && p %in% colnames(mm)) return(unname(mm[, p]))
  c(-Inf, Inf)
}

# EMC2's transforms (map.R): exp -> lower + exp(x), pnorm -> lower + (upper -
# lower) * pnorm(x); both increasing, so a box maps end to end.
nn_to_sampled <- function(x, func, lower, upper)
  switch(func, identity = x, exp = log(pmax(x - lower, 0)),
         pnorm = stats::qnorm(pmin(pmax((x - lower) / (upper - lower), 0), 1)),
         stop("Unsupported transform '", func, "'"))

nn_to_natural <- function(x, func, lower, upper)
  switch(func, identity = x, exp = lower + exp(x), pnorm = lower + (upper - lower) * stats::pnorm(x),
         stop("Unsupported transform '", func, "'"))

# ---------------------------------------------------------------------------
# Priors inside the training box
# ---------------------------------------------------------------------------

#' Check that a Prior Lies Inside a Neural Likelihood's Training Region
#'
#' For every network input and design cell, maps the prior's mean +- `k` sd
#' (on the sampled scale, through the design matrices) to the natural scale
#' and compares it with the network's training region. Parameter vectors
#' outside the region are rejected by the network, so such a prior is
#' silently truncated in a fit and breaks the calibration of an SBC: every
#' replicate whose true parameters leave the region cannot be recovered.
#' [run_sbc()] runs this check (with `k = getOption("emc.nn_box_k", 4)`) for
#' every neural-likelihood design, and [nn_cell()] for its priors.
#'
#' @param x An [nn_cell()], an `emc.design` of a neural-likelihood model, or
#'   an `emc.prior` made for one (its design is then taken from the prior).
#' @param prior An `emc.prior` (made with [prior()]) when `x` is a design.
#'   For a hierarchical prior the group-mean prior is checked; subjects
#'   spread further, so passing is necessary, not sufficient.
#' @param k Half-width in prior sds.
#' @param refuse If `TRUE`, an error lists every parameter and cell whose
#'   range leaves the region; if `FALSE`, the table is returned regardless.
#' @return (Invisibly) a data frame with one row per network input and
#'   design cell: the natural-scale range of the prior at +- `k` sd, the
#'   training region, and whether the former lies inside the latter.
#' @examples
#' des <- design(factors = list(subjects = 1), Rlevels = c("a", "b"), model = RDMnn,
#'               formula = list(v ~ 0 + lR, B ~ 1, t0 ~ 1, A ~ 1), constants = c(s = log(1)))
#' pri <- prior(des, type = "single",
#'              pmean = c(v_lRa = log(2), v_lRb = log(1), B = 0, t0 = log(.2), A = log(.3)),
#'              psd = c(v_lRa = .1, v_lRb = .1, B = .1, t0 = .1, A = .1))
#' nn_in_box(des, pri)
#' @export
nn_in_box <- function(x, prior = NULL, k = 4, refuse = TRUE) {
  if (inherits(x, "nn_cell")) {
    des <- x$design; mu <- x$mean; sdv <- x$sd
    if (missing(k)) k <- x$k
  } else {
    if (inherits(x, "emc.prior")) { prior <- x; x <- attr(prior, "design") }
    des <- if (inherits(x, "emc.design")) x else
      if (is.list(x) && length(x) == 1L && inherits(x[[1]], "emc.design")) x[[1]] else
        stop("x must be an nn_cell, a design or a prior")
    if (is.null(prior)) stop("nn_in_box() needs the prior (an emc.prior made with prior())")
    mu <- prior$theta_mu_mean
    sdv <- sqrt(diag(as.matrix(prior$theta_mu_var)))
    names(sdv) <- names(mu)
  }
  reg <- nn_model_fun(des)()$nn
  box <- nn_box_table(des, reg, mu, sdv, k)
  if (refuse) nn_refuse_box(box, reg, k, "the prior")
  invisible(box)
}

nn_box_table <- function(des, reg, mu, sdv, k) {
  mp <- sampled_pars(des, doMap = TRUE, add_da = TRUE)
  map <- attr(mp, "map")
  tr <- fill_transform(NULL, des$model)
  consts <- des$constants
  if (length(des$Fcovariates))
    warning("covariate effects are not included in the training-region check")
  rows <- list()
  for (p in intersect(reg$pars, names(map))) {
    m <- as.data.frame(map[[p]])
    X <- nn_map_matrix(m)
    labs <- m[, setdiff(names(m), colnames(X)), drop = FALSE]
    cols <- colnames(X)
    mean_c <- ifelse(cols %in% names(mu), mu[cols], consts[cols])
    sd_c <- ifelse(cols %in% names(sdv), sdv[cols], 0)
    if (anyNA(mean_c)) stop("cannot map ", p, ": ", paste(cols[is.na(mean_c)], collapse = ", "),
                            " is neither sampled nor a constant")
    centre <- drop(X %*% mean_c); half <- k * drop(abs(X) %*% sd_c)
    lo <- nn_to_natural(centre - half, tr$func[[p]], tr$lower[[p]], tr$upper[[p]])
    hi <- nn_to_natural(centre + half, tr$func[[p]], tr$lower[[p]], tr$upper[[p]])
    cell <- if (ncol(labs)) apply(labs, 1, function(r) paste(names(labs), r, sep = "=", collapse = ",")) else ""
    rows[[p]] <- unique(data.frame(parameter = p, cell = cell, prior_lower = lo, prior_upper = hi,
                                   box_lower = reg$lower[[p]], box_upper = reg$upper[[p]],
                                   # a zerobox input's region is closed at 0 (exact zeros)
                                   inside = (if (p %in% reg$zerobox) lo >= reg$lower[[p]] else lo > reg$lower[[p]]) &
                                     hi < reg$upper[[p]],
                                   stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

nn_refuse_box <- function(box, reg, k, what) {
  bad <- box[!box$inside, , drop = FALSE]
  if (!nrow(bad)) return(invisible(TRUE))
  lines <- sprintf("  %s%s: [%s, %s] vs training region [%s, %s]", bad$parameter,
                   ifelse(nzchar(bad$cell), paste0(" (", bad$cell, ")"), ""),
                   signif(bad$prior_lower, 4), signif(bad$prior_upper, 4),
                   signif(bad$box_lower, 4), signif(bad$box_upper, 4))
  stop("Neural likelihood '", reg$label, "': ", what, " leaves the network's training region within +-",
       k, " sd (natural scale):\n", paste(lines, collapse = "\n"),
       "\nParameter vectors outside the region are rejected by the network, so a fit truncates the ",
       "prior there and an SBC cannot be calibrated. Narrow or move the prior (see ?nn_cell).",
       call. = FALSE)
}

nn_format_box <- function(box) {
  f <- function(z) formatC(z, digits = 4, format = "g")
  data.frame(parameter = box$parameter, cell = box$cell,
             prior = paste0("[", f(box$prior_lower), ", ", f(box$prior_upper), "]"),
             region = paste0("[", f(box$box_lower), ", ", f(box$box_upper), "]"),
             inside = box$inside)
}

# run_sbc(): refuse a neural-likelihood design whose prior leaves the box.
nn_sbc_box_check <- function(design_in, prior_in) {
  des <- if (inherits(design_in, "emc.design")) design_in else
    if (is.list(design_in) && length(design_in) == 1L && inherits(design_in[[1]], "emc.design"))
      design_in[[1]] else NULL
  if (is.null(des) || !is.function(des$model) || is.null(des$model()$nn)) return(invisible(NULL))
  nn_in_box(des, prior_in, k = getOption("emc.nn_box_k", 4))
}

# ---------------------------------------------------------------------------
# Evaluation helpers
# ---------------------------------------------------------------------------

nn_dadm <- function(data, des)
  suppressMessages(design_model(data, des, compress = FALSE, rt_resolution = 1e-12,
                                verbose = FALSE, rt_check = FALSE))

# One data row per design cell (factor levels crossed, one subject).
nn_cells <- function(des) {
  fac <- des$Ffactors[setdiff(names(des$Ffactors), "subjects")]
  if (!length(fac)) return(data.frame(row.names = 1L))
  expand.grid(lapply(fac, function(l) factor(l, levels = l)), KEEP.OUT.ATTRS = FALSE)
}

# Data with, for each cell c and each response, the times rts[[c]]; `group`
# numbers the (cell, response) blocks in order, `cell` gives each row's cell.
nn_grid_data <- function(des, cells, rts) {
  Rl <- des$Rlevels
  pieces <- grp <- cel <- list()
  g <- 0L
  for (c in seq_len(max(1L, nrow(cells)))) for (r in Rl) {
    g <- g + 1L
    rt <- rts[[c]]
    d <- data.frame(subjects = factor(1), R = factor(rep(r, length(rt)), levels = Rl), rt = rt)
    if (ncol(cells)) d <- cbind(d, cells[rep(c, length(rt)), , drop = FALSE])
    pieces[[g]] <- d; grp[[g]] <- rep(g, length(rt)); cel[[g]] <- rep(c, length(rt))
  }
  data <- do.call(rbind, pieces)
  rownames(data) <- NULL
  list(data = data, group = unlist(grp), cell = unlist(cel))
}

# A dadm with one trial per design cell (for mapping parameters to cells).
nn_cells_dadm <- function(des, cells) {
  nc <- max(1L, nrow(cells))
  d <- data.frame(subjects = factor(rep(1, nc)), R = factor(rep(des$Rlevels[1], nc), levels = des$Rlevels),
                  rt = rep(1, nc))
  if (ncol(cells)) d <- cbind(d, cells)
  nn_dadm(d, des)
}

# Support start (earliest possible time) per design cell at sampled vector p.
nn_cell_support <- function(p, cells_dadm, model, support) {
  pars <- get_pars_matrix_oo(p, cells_dadm, model)
  s <- support(pars)
  as.numeric(tapply(s, cells_dadm$trials, min))
}

nn_support_fun <- function(support, reg) {
  if (!is.null(support)) {
    if (!is.function(support)) stop("support must be a function(pars) returning each row's earliest time")
    return(support)
  }
  shift <- if (is.character(reg$pre)) reg$pre else "t0"
  function(pars) if (shift %in% colnames(pars)) pars[, shift] else rep(0, nrow(pars))
}

# Trial-wise log-likelihoods (rows = data rows, columns = particles) as a fit
# computes them, floored at min_ll. Compiled models through calc_ll_manager;
# R-path models from dfun/pfun exactly as log_likelihood_ddm/race combine them.
nn_trialwise <- function(P, dadm, model, min_ll = log(1e-10)) {
  ml <- model()
  if (!is.null(ml$c_name))
    return(matrix(calc_ll_manager(P, dadm, model, return_trialwise = TRUE), ncol = nrow(P)))
  out <- matrix(NA_real_, length(attr(dadm, "expand")), nrow(P))
  for (i in seq_len(nrow(P))) {
    pars <- get_pars_matrix_oo(P[i, ], dadm, model)
    ok <- attr(pars, "ok"); if (is.null(ok)) ok <- rep(TRUE, nrow(pars))
    if (ml$type == "DDM") {
      like <- numeric(nrow(dadm))
      if (any(ok)) like[ok] <- ml$dfun(dadm$rt[ok], dadm$R[ok], pars[ok, , drop = FALSE])
      like[is.na(like)] <- 0
      out[, i] <- pmax(min_ll, log(like))[attr(dadm, "expand")]
    } else {
      if (!is.null(dadm$RACE)) pars[as.numeric(dadm$lR) > as.numeric(as.character(dadm$RACE)), ] <- NA
      w <- dadm$winner
      lds <- numeric(nrow(dadm))
      lds[w] <- log(ml$dfun(dadm$rt[w], pars[w, , drop = FALSE]))
      n_acc <- length(levels(dadm$R))
      if (n_acc > 1) lds[!w] <- log(1 - ml$pfun(dadm$rt[!w], pars[!w, , drop = FALSE]))
      lds[is.na(lds) | !ok] <- min_ll
      ll <- lds[w]
      if (n_acc > 1) ll <- ll + colSums(matrix(lds[!w], nrow = n_acc - 1))
      ll[is.na(ll)] <- min_ll
      out[, i] <- pmax(min_ll, ll)[attr(dadm, "expand")]
    }
  }
  out
}

# Is every row of the design within the model's bounds at each particle?
nn_particles_ok <- function(P, dadm, model)
  vapply(seq_len(nrow(P)), function(i) {
    pars <- get_pars_matrix_oo(P[i, ], dadm, model)
    ok <- attr(pars, "ok")
    is.null(ok) || all(ok)
  }, TRUE)

# Draws from the cell prior, truncated at +-k sd.
nn_cell_draws <- function(cell, n) {
  sp <- names(cell$mean)
  out <- matrix(NA_real_, n, length(sp), dimnames = list(NULL, sp))
  for (i in seq_len(n)) {
    repeat { z <- stats::rnorm(length(sp)); if (all(abs(z) <= cell$k)) break }
    out[i, ] <- cell$mean + cell$sd * z
  }
  out
}

nn_theta <- function(cell, theta, n, seed) {
  sp <- names(cell$mean)
  if (is.null(theta)) {
    if (!is.null(seed)) set.seed(seed)
    return(nn_cell_draws(cell, n))
  }
  theta <- as.matrix(theta)
  if (is.null(colnames(theta)) || !all(sp %in% colnames(theta)))
    stop("theta must be a matrix with columns named by the sampled parameters (", paste(sp, collapse = ", "), ")")
  theta[, sp, drop = FALSE]
}

nn_lapply <- function(X, FUN, cores) {
  out <- if (cores > 1) auto_mclapply(X, FUN, mc.cores = cores) else lapply(X, FUN)
  err <- vapply(out, inherits, TRUE, "try-error")
  if (any(err)) stop("a worker failed: ", as.character(out[[which(err)[1]]]))
  out
}

nn_mean_se <- function(v) {
  v <- v[is.finite(v)]
  c(mean = if (length(v)) mean(v) else NA_real_,
    se = if (length(v) > 1) stats::sd(v) / sqrt(length(v)) else NA_real_, n = length(v))
}

# ---------------------------------------------------------------------------
# Total mass
# ---------------------------------------------------------------------------

#' Total Probability Mass of a Neural Likelihood
#'
#' Integrates the network's likelihood over rt and sums over responses,
#' `Z(theta) = sum_R int p(rt, R | theta) drt`, for parameter vectors drawn
#' from the cell prior, and reports how much `log Z` varies over them. A
#' normalizing flow integrates to one by construction; a network that
#' regresses the log-likelihood directly does not, which is harmless while
#' `log Z` is constant in theta (it cancels from the posterior) and biases
#' inference where it is not. The same integral of the analytic control, on
#' the same grid, separates the network from the grid (mass beyond `rt_max`,
#' which is large for slow parameter vectors).
#'
#' @details The grid, per design cell, is the union of a log-spaced grid
#'   from 1 ms to `rt_max` and one from 0.1 ms past the support start (by
#'   default `t0`; see `support`) to `rt_max`, each with `n_grid` points;
#'   the integral is a trapezoid sum. Densities are the ones a fit uses,
#'   floored at EMC2's `min_ll` (`1e-10`), so masses below about `1e-9` read
#'   as zero. `before_support` is the mass the network puts before the
#'   support start (where the exact density is zero).
#'
#' @param x An [nn_cell()], or a neural-likelihood model (then `...` go to
#'   [nn_cell()]).
#' @param theta A matrix of sampled-scale parameter vectors (columns named
#'   by the cell's sampled parameters); default `n` draws from the cell
#'   prior.
#' @param n Number of prior draws when `theta` is `NULL`.
#' @param rt_max Upper end of the rt grid (seconds).
#' @param n_grid Points per log-spaced grid.
#' @param support `function(pars)` giving, for each row of a natural-scale
#'   parameter matrix, the earliest possible time; default `t0` (or the
#'   `pre` parameter), else 0.
#' @param seed Seed for the prior draws (`NULL`: leave the RNG alone).
#' @param cores Parallel workers over draws (forked; 1 = serial).
#' @param ... Arguments for [nn_cell()] when `x` is a model.
#' @return An `nn_total_mass` object: `summary` (per cell: mean and sd of
#'   `log Z`, the same for the control and for the difference, and the largest
#'   `before_support`), and `mass` (one row per draw and cell).
#' @examples
#' \donttest{
#' nn_total_mass(RDMnn, n = 5, n_grid = 300)
#' }
#' @export
nn_total_mass <- function(x, theta = NULL, n = 50, rt_max = 20, n_grid = 1000, support = NULL,
                          seed = 1, cores = 1, ...) {
  cell <- nn_as_cell(x, ...)
  theta <- nn_theta(cell, theta, n, seed)
  support <- nn_support_fun(support, cell$reg)
  des <- cell$design
  cells <- nn_cells(des)
  cells_dadm <- nn_cells_dadm(des, cells)
  a_grid <- exp(seq(log(1e-3), log(rt_max), length.out = n_grid))
  trapz <- function(x, y) if (length(x) < 2) 0 else sum((y[-1] + y[-length(y)]) / 2 * diff(x))
  one <- function(i) {
    th <- theta[i, , drop = FALSE]
    sup <- nn_cell_support(th[1, ], cells_dadm, des$model, support)
    rts <- lapply(sup, function(s) sort(unique(c(a_grid, if (is.finite(s) && s < rt_max)
      s + exp(seq(log(1e-4), log(rt_max - s), length.out = n_grid))))))
    g <- nn_grid_data(des, cells, rts)
    dens_n <- exp(nn_trialwise(th, nn_dadm(g$data, des), des$model)[, 1])
    dens_c <- if (!is.null(cell$control))
      exp(nn_trialwise(th, nn_dadm(g$data, cell$control_design), cell$control_design$model)[, 1])
    do.call(rbind, lapply(seq_along(rts), function(c) {
      z <- function(d, upto = Inf) sum(vapply(unique(g$group[g$cell == c]), function(k) {
        r <- g$group == k & g$data$rt <= upto
        trapz(g$data$rt[r], d[r])
      }, 0))
      data.frame(draw = i, cell = c, support = sup[c], log_Z = log(z(dens_n)),
                 log_Z_control = if (is.null(dens_c)) NA_real_ else log(z(dens_c)),
                 before_support = z(dens_n, sup[c]))
    }))
  }
  mass <- do.call(rbind, nn_lapply(seq_len(nrow(theta)), one, cores))
  d <- mass$log_Z - mass$log_Z_control
  summ <- do.call(rbind, lapply(split(seq_len(nrow(mass)), mass$cell), function(r)
    data.frame(cell = mass$cell[r[1]], n = length(r),
               mean_log_Z = mean(mass$log_Z[r]), sd_log_Z = stats::sd(mass$log_Z[r]),
               mean_log_Z_control = mean(mass$log_Z_control[r]),
               sd_log_Z_control = stats::sd(mass$log_Z_control[r]),
               sd_difference = stats::sd(d[r]), max_abs_difference = max(abs(d[r])),
               max_before_support = max(mass$before_support[r]))))
  structure(list(summary = summ, mass = mass, theta = theta, rt_max = rt_max, label = cell$reg$label,
                 control_label = cell$control_label), class = "nn_total_mass")
}

#' @export
print.nn_total_mass <- function(x, digits = 4, ...) {
  cat("Total mass of neural likelihood '", x$label, "' over ", length(unique(x$mass$draw)),
      " prior draws (rt up to ", x$rt_max, " s); control: ", x$control_label, "\n", sep = "")
  print(format(x$summary, digits = digits), row.names = FALSE)
  cat("log Z constant in theta (sd ~0) is harmless; sd_difference separates the network from the grid.\n")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Score bias
# ---------------------------------------------------------------------------

#' Score Bias of a Neural Likelihood
#'
#' The expected score of the network's log-likelihood under the true model,
#' at the true parameters: for each parameter vector theta* drawn from the
#' cell prior and each sampled parameter j,
#' `g_j = d/d theta_j E_{y ~ p(. | theta*)}[log f(y | theta)]` at
#' `theta = theta*`, per trial. An exact likelihood has `g = 0` (the expected
#' score is zero); a network whose error varies with theta has not, and
#' pulls the posterior in the direction of `g`. The same estimator applied to
#' the analytic control must return about zero: it is part of the estimator,
#' not an optional check (a control that is not zero means the estimator,
#' not the network, is off).
#'
#' @details The expectation is a quadrature over rt: per design cell and
#'   response, `n_grid` log-spaced times from `rt_min` to `rt_max` past the
#'   support start, weighted by the control's likelihood at theta* times the
#'   interval widths (without the widths the control returns -0.25, not 0),
#'   and normalised by the total weight. The derivative is a central
#'   difference with step `h` on the sampled scale (halving `h` should
#'   shrink the control by four; at `h = 0.05` the control can be a sizeable
#'   fraction of the signal). All likelihoods are those a fit uses: for a
#'   race, the winner's density times the losers' survivors; floored at
#'   `min_ll`.
#'
#'   The support start (`t0`, or the `pre` parameter) moves the support, where
#'   this estimator does not apply; parameters entering it are excluded by
#'   default. Score bias is local: it ranks networks of one architecture but
#'   does not predict SBC bias across architectures; use
#'   [nn_posterior_shift()] for that.
#'
#' @param x An [nn_cell()] (with a control), or a neural-likelihood model
#'   (then `...` go to [nn_cell()]).
#' @param n_draws Number of prior draws (when `theta` is `NULL`).
#' @param h Finite-difference step on the sampled scale.
#' @param pars Sampled parameters to differentiate; default all except those
#'   entering the support start.
#' @param n_check Number of draws on which the control's score bias is
#'   computed.
#' @param theta A matrix of sampled-scale parameter vectors; default draws
#'   from the cell prior.
#' @param rt_min,rt_max,n_grid The rt grid past the support start.
#' @param support `function(pars)`: the earliest possible time per row of a
#'   natural-scale parameter matrix; default `t0` (or the `pre` parameter),
#'   else 0.
#' @param seed Seed for the prior draws (`NULL`: leave the RNG alone).
#' @param cores Parallel workers over draws (forked; 1 = serial).
#' @param ... Arguments for [nn_cell()] when `x` is a model.
#' @return An `nn_score_bias` object: `summary` (per parameter: mean `g` and
#'   its standard error over draws, the same for the control), the per-draw
#'   values `g` and `g_control`, and `grid_mass`, the probability (per cell) the
#'   rt grid holds at each draw: the estimates are conditional on the grid, so
#'   for slow parameter vectors (much mass beyond `rt_max`) the control is not
#'   zero until `rt_max` is raised.
#' @examples
#' \donttest{
#' nn_score_bias(RDMnn, n_draws = 5, n_check = 2)
#' }
#' @export
nn_score_bias <- function(x, n_draws = 100, h = 0.01, pars = NULL, n_check = 40, theta = NULL,
                          rt_min = 1e-3, rt_max = 8, n_grid = 120, support = NULL, seed = 1,
                          cores = 1, ...) {
  cell <- nn_as_cell(x, ...)
  nn_need_control(cell, "the score bias")
  des <- cell$design; dc <- cell$control_design
  sp <- names(cell$mean)
  theta <- nn_theta(cell, theta, n_draws, seed)
  support <- nn_support_fun(support, cell$reg)
  if (is.null(pars)) {
    map <- attr(sampled_pars(des, doMap = TRUE), "map")
    moving <- intersect(c("t0", if (is.character(cell$reg$pre)) cell$reg$pre), names(map))
    excl <- unique(unlist(lapply(moving, function(p) {
      X <- nn_map_matrix(map[[p]]); colnames(X)[colSums(X != 0) > 0] })))
    pars <- setdiff(sp, excl)
  }
  if (!length(pars) || !all(pars %in% sp))
    stop("pars must be sampled parameters of the cell (", paste(sp, collapse = ", "), ")")
  cells <- nn_cells(des)
  cells_dadm <- nn_cells_dadm(des, cells)
  off <- exp(seq(log(rt_min), log(rt_max), length.out = n_grid))
  d <- diff(off)
  dr <- c(d[1], (utils::head(d, -1) + utils::tail(d, -1)) / 2, utils::tail(d, 1))
  np <- length(pars)
  one <- function(i) {
    th <- theta[i, ]
    na <- stats::setNames(rep(NA_real_, np), pars)
    sup <- nn_cell_support(th, cells_dadm, des$model, support)
    if (any(!is.finite(sup))) return(list(g = na, gc = na, grid_mass = NA_real_))
    g <- nn_grid_data(des, cells, lapply(sup, function(s) s + off))
    w_dr <- rep(dr, length(unique(g$group)))
    dn <- nn_dadm(g$data, des); dcn <- nn_dadm(g$data, dc)
    P <- matrix(th, 2 * np + 1, length(th), byrow = TRUE, dimnames = list(NULL, sp))
    for (j in seq_len(np)) {
      P[1 + j, pars[j]] <- P[1 + j, pars[j]] + h
      P[1 + np + j, pars[j]] <- P[1 + np + j, pars[j]] - h
    }
    w <- exp(nn_trialwise(P[1, , drop = FALSE], dcn, dc$model)[, 1]) * w_dr
    lbar <- function(lp, ok) {
      out <- apply(lp, 2, function(l) { u <- is.finite(l) & is.finite(w) & w > 0; sum(w[u] * l[u]) / sum(w[u]) })
      out[!ok] <- NA
      out
    }
    slope <- function(L) stats::setNames((L[1 + seq_len(np)] - L[1 + np + seq_len(np)]) / (2 * h), pars)
    gn <- slope(lbar(nn_trialwise(P, dn, des$model), nn_particles_ok(P, dn, des$model)))
    gc <- if (i <= n_check) slope(lbar(nn_trialwise(P, dcn, dc$model), nn_particles_ok(P, dcn, dc$model))) else na
    # probability the grid holds (per cell); the estimator conditions on it
    list(g = gn, gc = gc, grid_mass = sum(w) / length(sup))
  }
  res <- nn_lapply(seq_len(nrow(theta)), one, cores)
  G <- do.call(rbind, lapply(res, `[[`, "g"))
  A <- do.call(rbind, lapply(res, `[[`, "gc"))[seq_len(min(n_check, nrow(theta))), , drop = FALSE]
  sg <- apply(G, 2, nn_mean_se); sa <- apply(A, 2, nn_mean_se)
  summ <- data.frame(parameter = pars, g = sg["mean", ], se = sg["se", ], control = sa["mean", ],
                     control_se = sa["se", ], n = sg["n", ], n_control = sa["n", ], row.names = NULL)
  structure(list(summary = summ, g = G, g_control = A, theta = theta, h = h,
                 grid_mass = vapply(res, `[[`, 0, "grid_mass"), rt_max = rt_max, label = cell$reg$label,
                 control_label = cell$control_label), class = "nn_score_bias")
}

#' @export
print.nn_score_bias <- function(x, ...) {
  s <- x$summary
  cat("Score bias of neural likelihood '", x$label, "' (control: ", x$control_label, ")\n", sep = "")
  cat(sprintf("  per trial, sampled scale, h = %g; mean (se) over %d prior draws, control over %d\n",
              x$h, max(s$n), max(s$n_control)))
  for (r in seq_len(nrow(s)))
    cat(sprintf("  %-10s %+9.5f (%.5f)   control %+9.5f (%.5f)\n", s$parameter[r], s$g[r], s$se[r],
                s$control[r], s$control_se[r]))
  ok <- nn_control_ok(s$control, s$control_se)
  cat(if (all(ok)) "  Control ~ 0: yes.\n" else
    sprintf("  Control ~ 0: NO for %s -- the estimator, not the network, is off (step h, grid, weights).\n",
            paste(s$parameter[!ok], collapse = ", ")))
  gm <- min(x$grid_mass, na.rm = TRUE)
  if (gm < .999)
    cat(sprintf(paste0("  The rt grid (to %g s past the support) holds as little as %.4f of the probability:\n",
                       "  the estimates condition on it; raise rt_max if the control is not ~0.\n"), x$rt_max, gm))
  cat("  Score bias is necessary, not sufficient (local); nn_posterior_shift() predicts SBC.\n")
  invisible(x)
}

nn_control_ok <- function(m, se) !is.na(m) & (abs(m) < 1e-3 | abs(m) < 3 * ifelse(is.na(se), 0, se))

# ---------------------------------------------------------------------------
# Posterior shift
# ---------------------------------------------------------------------------

#' Posterior Shift of a Neural Likelihood Against Its Analytic Control
#'
#' What an SBC cell computes, without the sampler: for `n_datasets` data
#' sets simulated from the cell prior (by the control's simulator), the
#' posterior under the network and under the analytic control, each from a
#' Laplace approximation refined by importance sampling from a defensive
#' mixture of both. Reports, per sampled parameter, the standardised bias of
#' each posterior mean (`mean(post mean - truth) / sd(post mean - truth)`
#' over data sets, SBC's statistic with the mean for the median), the shift
#' of the network's posterior mean from the control's in control posterior
#' sds, and the ratio of posterior sds.
#'
#' @details The control column is the check's own control: it should be
#'   within Monte Carlo error of zero (its standard error is about
#'   `1 / sqrt(n)`; the print flags values beyond `3 / sqrt(n)`). Both
#'   columns share the same truths and data, so their Monte Carlo errors are
#'   strongly correlated: the paired `shift` (and the difference between the
#'   columns) measures the network's contribution far more precisely than
#'   either column does. In the NLE project this predicted the standardised bias of
#'   500-replicate SBC cells within about 0.1 (at ~1 minute per model against
#'   hours per cell). Use it to screen, not to certify: importance sampling
#'   from Laplace proposals can miss mass a sampler finds; confirm the
#'   chosen network with an SBC cell and its matched control
#'   ([nn_sbc_cell()]).
#'
#'   Per data set (seed `seed + i`): truth drawn from the cell prior
#'   (truncated at +- k sd); `n_trials` trials by [make_data()] on the
#'   control's design; a
#'   BFGS mode from the truth (plus `multistart` random starts) with the
#'   Hessian's eigenvalues floored at the smallest prior precision; proposal
#'   an equal mixture of the Laplace normals with covariance `inflate^2`
#'   times theirs; `n_draws` draws with self-normalised weights for each
#'   posterior. Likelihoods as in a fit (floored at `min_ll`).
#'
#' @param x An [nn_cell()] (with a control that can simulate), or a
#'   neural-likelihood model (then `...` go to [nn_cell()]).
#' @param n_datasets Number of simulated data sets.
#' @param n_trials Trials per data set (per design cell).
#' @param n_draws Importance-sampling draws per data set.
#' @param seed Data set `i` uses `set.seed(seed + i)`.
#' @param cores Parallel workers over data sets (forked; 1 = serial).
#' @param multistart Extra random BFGS starts per posterior (0: start at the
#'   truth only); with starts, the local network mode is kept as a third
#'   mixture component.
#' @param inflate Scale factor of the proposal components' sds.
#' @param ... Arguments for [nn_cell()] when `x` is a model.
#' @return An `nn_posterior_shift` object: `summary` (per parameter:
#'   `control` and `nn` standardised bias, `shift`, `sd_ratio`), `ess`
#'   (per data set, both posteriors), `results` (per data set: truth,
#'   posterior means and sds, modes) and the settings.
#' @examples
#' \donttest{
#' nn_posterior_shift(RDMnn, n_datasets = 2, n_trials = 100, n_draws = 200)
#' }
#' @export
nn_posterior_shift <- function(x, n_datasets = 48, n_trials = 400, n_draws = 2000, seed = 7000,
                               cores = 1, multistart = 0, inflate = 1.5, ...) {
  cell <- nn_as_cell(x, ...)
  nn_need_control(cell, "the posterior shift")
  sp <- names(cell$mean)
  start <- proc.time()[["elapsed"]]
  one <- function(i) nn_ps_dataset(cell, i, seed, n_trials, n_draws, multistart, inflate)
  res <- nn_lapply(seq_len(n_datasets), one, cores)
  failed <- vapply(res, function(r) !is.null(r$failed), TRUE)
  if (any(failed))
    warning(sum(failed), " of ", n_datasets, " data sets could not be simulated from the control (first: ",
            res[[which(failed)[1]]]$failed, ")", call. = FALSE)
  res <- res[!failed]
  if (length(res) < 2) stop("fewer than two usable data sets")
  get <- function(k, s) do.call(rbind, lapply(res, function(r) r[[k]][[s]]))
  tr <- do.call(rbind, lapply(res, `[[`, "truth"))
  ea <- get("an", "mean") - tr; ef <- get("fl", "mean") - tr
  sb <- function(e) colMeans(e) / apply(e, 2, stats::sd)
  summ <- data.frame(parameter = sp, control = sb(ea), nn = sb(ef),
                     shift = colMeans((get("fl", "mean") - get("an", "mean")) / get("an", "sd")),
                     sd_ratio = colMeans(get("fl", "sd") / get("an", "sd")), row.names = NULL)
  ess <- data.frame(control = vapply(res, function(r) r$an$ess, 0), nn = vapply(res, function(r) r$fl$ess, 0))
  structure(list(summary = summ, ess = ess, results = res, n_datasets = n_datasets, n_used = length(res),
                 n_trials = n_trials, n_draws = n_draws, minutes = (proc.time()[["elapsed"]] - start) / 60,
                 label = cell$reg$label, control_label = cell$control_label), class = "nn_posterior_shift")
}

#' @export
print.nn_posterior_shift <- function(x, ...) {
  cat(sprintf("Posterior shift of neural likelihood '%s' (control: %s)\n", x$label, x$control_label))
  cat(sprintf("  %d/%d data sets x %d trials, %d draws, %.1f min\n", x$n_used, x$n_datasets,
              x$n_trials, x$n_draws, x$minutes))
  cat(sprintf("  %-10s %10s %10s %12s %10s\n", "parameter", "control", "network", "shift (sd)", "sd ratio"))
  s <- x$summary
  for (r in seq_len(nrow(s)))
    cat(sprintf("  %-10s %+10.3f %+10.3f %+12.3f %10.3f\n", s$parameter[r], s$control[r], s$nn[r],
                s$shift[r], s$sd_ratio[r]))
  cat(sprintf("  ESS control median %.0f (min %.0f), network median %.0f (min %.0f)\n",
              stats::median(x$ess$control), min(x$ess$control), stats::median(x$ess$nn), min(x$ess$nn)))
  lim <- 3 / sqrt(x$n_used)
  bad <- s$parameter[abs(s$control) > lim]
  cat(sprintf("  Control within +-%.2f (3/sqrt(n)): %s\n", lim,
              if (length(bad)) paste("NO for", paste(bad, collapse = ", ")) else "yes"))
  cat("  Screen, not certify: confirm with an SBC cell and its matched control (nn_sbc_cell()).\n")
  invisible(x)
}

nn_ps_dataset <- function(cell, i, seed, n_trials, M, k_ms, inflate) {
  set.seed(seed + i)
  truth <- nn_cell_draws(cell, 1)[1, ]
  # the control is the reference model: data come from its simulator
  why <- NULL
  dat <- withCallingHandlers(
    tryCatch(suppressMessages(make_data(truth, cell$control_design, n_trials = n_trials)),
             error = function(e) { why <<- conditionMessage(e); NULL }),
    warning = function(w) { why <<- conditionMessage(w); invokeRestart("muffleWarning") })
  if (!is.data.frame(dat))
    return(list(failed = if (is.null(why)) "make_data() returned no data" else why))
  dn <- nn_dadm(dat, cell$design); dc <- nn_dadm(dat, cell$control_design)
  lp_n <- function(X) nn_lpost(X, dn, cell$design$model, cell$mean, cell$sd)
  lp_c <- function(X) nn_lpost(X, dc, cell$control_design$model, cell$mean, cell$sd)
  la <- nn_laplace(lp_c, truth, cell$sd, k_ms)
  lf <- nn_laplace(lp_n, truth, cell$sd, k_ms)
  comps <- list(la, lf)
  if (k_ms > 0) comps <- c(comps, list(nn_laplace(lp_n, truth, cell$sd, 0)))
  K <- length(comps)
  X <- do.call(rbind, lapply(comps, function(cc) nn_rmvn(M %/% K, cc$m, inflate^2 * cc$S)))
  colnames(X) <- names(truth)
  Q <- vapply(comps, function(cc) nn_dmvn(X, cc$m, inflate^2 * cc$S), numeric(nrow(X)))
  qm <- apply(Q, 1, max)
  lq <- qm + log(rowMeans(exp(Q - qm)))
  list(truth = truth, an = nn_wstats(lp_c(X) - lq, X), fl = nn_wstats(lp_n(X) - lq, X),
       map_a = la$m, map_f = lf$m)
}

# log posterior (log-likelihood as a fit computes it + independent normal
# prior) for the rows of X
nn_lpost <- function(X, dadm, model, mean, sd) {
  X <- matrix(X, ncol = length(mean), dimnames = list(NULL, names(mean)))
  as.numeric(calc_ll_manager(X, dadm, model)) + colSums(stats::dnorm(t(X), mean, sd, log = TRUE))
}

nn_laplace <- function(f, x0, sd, k) {
  nms <- names(x0)
  nlp <- function(x) { v <- -f(matrix(x, 1, dimnames = list(NULL, nms))); if (is.finite(v)) v else 1e10 }
  gr <- function(x) {                        # central differences, one batched likelihood call
    p <- length(x); eps <- 1e-4
    X <- rbind(matrix(x, p, p, byrow = TRUE) + diag(eps, p), matrix(x, p, p, byrow = TRUE) - diag(eps, p))
    v <- f(X)
    g <- -(v[seq_len(p)] - v[p + seq_len(p)]) / (2 * eps)
    g[!is.finite(g)] <- 0
    g
  }
  fit <- function(s) stats::optim(s, nlp, gr, method = "BFGS", control = list(maxit = 500, reltol = 1e-10))
  o <- fit(x0)
  for (j in seq_len(k)) {
    oj <- fit(x0 + 2 * sd * stats::rnorm(length(x0)))
    if (oj$value < o$value) o <- oj
  }
  H <- stats::optimHess(o$par, nlp, gr)
  e <- eigen((H + t(H)) / 2, symmetric = TRUE)
  ev <- pmax(e$values, min(1 / sd^2))
  list(m = stats::setNames(o$par, nms), S = e$vectors %*% diag(1 / ev, length(ev)) %*% t(e$vectors))
}

nn_dmvn <- function(X, m, S) {
  L <- chol(S)
  z <- backsolve(L, t(X) - m, transpose = TRUE)
  -colSums(z^2) / 2 - sum(log(diag(L))) - ncol(X) / 2 * log(2 * pi)
}

nn_rmvn <- function(n, m, S) sweep(matrix(stats::rnorm(n * length(m)), n) %*% chol(S), 2, m, "+")

nn_wstats <- function(lw, X) {
  w <- exp(lw - max(lw))
  w[!is.finite(w)] <- 0
  w <- w / sum(w)
  mu <- colSums(w * X)
  list(mean = mu, sd = sqrt(colSums(w * sweep(X, 2, mu)^2)), ess = 1 / sum(w^2))
}

# ---------------------------------------------------------------------------
# SBC cells
# ---------------------------------------------------------------------------

#' Run (or Prepare) an SBC Cell for a Neural Likelihood, with Its Control
#'
#' Runs [run_sbc()] on a cell's design and prior for the network and, with
#' `run_control = TRUE`, on the same design and prior for the analytic
#' control: a cell with a matched control that differs from it in exactly one
#' respect, the likelihood. Summarises both with [nn_sbc_summary()] and, given
#' `archive_dir`, writes the run in the layout of the EMC2 SBC archive
#' (`README.md`, `bundle/`, `results/`; see the archive's `ARCHIVING.md`).
#'
#' @details With `run = FALSE` nothing is fitted: the bundle (the cell as
#'   `cell.rds`, a `run_cell.R` that re-runs this call from the bundle
#'   directory, and `session.txt` with the EMC2 build, BLAS and artefact
#'   lines) is written to `archive_dir`, ready to be copied to a cluster. A
#'   card registered by path must exist at the same path there; shipped
#'   artefacts are found by name.
#'
#'   The README is a draft: fields the function cannot know (branch and
#'   commit of the EMC2 checkout, why the cell was run, where) come from
#'   `info` or are marked `TODO`, and the verdict it writes is provisional
#'   (gate: |standardised bias| <= 0.5 for every parameter; ECDF envelopes and
#'   KS p-values reported). File it following `ARCHIVING.md` (add the
#'   `INDEX.md` row, commit).
#'
#'   **Cores.** In [run_sbc()], `cores_per_chain` is the number of replicates
#'   fitted at once, and each replicate's [fit()] also runs its chains in
#'   parallel (`cores_for_chains`, by default the number of chains): the
#'   total is their product. For a cell, parallel replicates are the
#'   efficient choice, e.g. `cores_per_chain = 8, cores_for_chains = 1` for 8
#'   cores in all.
#'
#' @param x An [nn_cell()], or a neural-likelihood model (a cell with the
#'   defaults is built).
#' @param trials,replicates Trials per replicate and number of replicates.
#' @param run_control Also run the control cell.
#' @param archive_dir Directory to write the run into (created); `NULL` for
#'   none. Result files go to `results/`, the SBC temporary files next to
#'   them: repeating an interrupted call resumes it (a model whose run
#'   finished is loaded, not repeated).
#' @param run `FALSE` to only write the bundle (needs `archive_dir`).
#' @param info Optional named list for the README: `branch`, `commit`,
#'   `why`, `where`, `topic`.
#' @param ... Passed to [run_sbc()] (and on to [fit()]), e.g.
#'   `cores_per_chain`, `cores_for_chains`, `stop_criteria`.
#' @return (Invisibly) a list: `nn` and `control` (the SBC objects),
#'   `summary` (a data frame, both models), `cell`.
#' @export
nn_sbc_cell <- function(x, trials = 400, replicates = 500, run_control = TRUE, archive_dir = NULL,
                        run = TRUE, info = list(), ...) {
  cell <- nn_as_cell(x)
  nn_in_box(cell)
  if (run_control) nn_need_control(cell, "an SBC cell with run_control = TRUE")
  dots <- list(...)
  if (!run && is.null(archive_dir)) stop("run = FALSE writes the bundle only; give archive_dir")
  f_nn <- f_ctl <- NULL
  if (!is.null(archive_dir)) {
    for (d in c("bundle", "results")) dir.create(file.path(archive_dir, d), recursive = TRUE, showWarnings = FALSE)
    nn_sbc_bundle(archive_dir, cell, trials, replicates, run_control, dots)
    f_nn <- file.path(archive_dir, "results", "sbc_nn.RData")
    f_ctl <- file.path(archive_dir, "results", "sbc_control.RData")
  }
  out <- list(nn = NULL, control = NULL, summary = NULL, cell = cell)
  if (run) {
    out$nn <- nn_sbc_run(cell$design, cell$prior, replicates, trials, f_nn, dots)
    if (run_control)
      out$control <- nn_sbc_run(cell$control_design, cell$control_prior, replicates, trials, f_ctl, dots)
    out$summary <- rbind(cbind(model = "network", nn_sbc_summary(out$nn)),
                         if (run_control) cbind(model = "control", nn_sbc_summary(out$control)))
  }
  if (!is.null(archive_dir)) nn_sbc_write(archive_dir, out, trials, replicates, dots, info)
  invisible(out)
}

# run_sbc() for one model of a cell. A finished run in `f` (its SBC saved and
# its temporary directory gone) is loaded, not repeated, so repeating an
# interrupted nn_sbc_cell() call resumes where it stopped.
nn_sbc_run <- function(des, pri, replicates, trials, f, dots) {
  if (!is.null(f) && file.exists(f) && !dir.exists(paste0(tools::file_path_sans_ext(f), "_temp"))) {
    env <- new.env(parent = emptyenv())
    load(f, envir = env)
    if (!is.null(env$SBC)) {
      message("Loading the finished SBC run in ", f)
      return(env$SBC)
    }
  }
  do.call(run_sbc, c(list(des, pri, replicates = replicates, trials = trials, fileName = f), dots))
}

#' Calibration Summary of an SBC Run
#'
#' Per parameter: the standardised bias `mean(median - truth) / sd(median -
#' truth)` over replicates (the NLE project's SBC gate is |value| <= 0.5), the
#' mean bias, coverage of the 95% credible interval, the Kolmogorov-Smirnov
#' p-value of the normalised ranks against uniformity, and whether the rank
#' ECDF stays inside the simultaneous 95% envelope of [plot_sbc_ecdf()].
#'
#' @param sbc A single-subject SBC object returned by [run_sbc()] (or
#'   [recover_sbc()]).
#' @param K Evaluation points of the ECDF envelope (as [plot_sbc_ecdf()]).
#' @return A data frame, one row per parameter.
#' @export
nn_sbc_summary <- function(sbc, K = 500) {
  if (is.null(sbc$rank$alpha) || is.null(sbc$bias$alpha))
    stop("nn_sbc_summary() needs a single-subject SBC result (run_sbc() with a 'single' prior)")
  rank <- as.matrix(sbc$rank$alpha); bias <- as.matrix(sbc$bias$alpha)
  cover <- as.matrix(sbc$coverage$alpha)
  N <- nrow(rank)
  lims <- get_lims(N, K, get_gamma(N, K))
  env <- apply(rank, 2, function(r) {
    y <- stats::ecdf(r[is.finite(r)])(lims$z) - lims$z
    all(y >= lims$lower - 1e-12 & y <= lims$upper + 1e-12)
  })
  data.frame(parameter = colnames(rank), std_bias = colMeans(bias) / apply(bias, 2, stats::sd),
             bias = colMeans(bias), coverage = colMeans(cover), ks_p = apply(rank, 2, .sbc_ks_p),
             in_envelope = env, n = N, row.names = NULL)
}

nn_sbc_bundle <- function(dir, cell, trials, replicates, run_control, dots) {
  b <- file.path(dir, "bundle")
  saveRDS(cell, file.path(b, "cell.rds"))
  arg <- function(nm, v) paste0(", ", nm, " = ", paste(deparse(v, width.cutoff = 500L), collapse = ""))
  extra <- paste(vapply(names(dots), function(nm) arg(nm, dots[[nm]]), ""), collapse = "")
  writeLines(c(
    "# Generated by EMC2::nn_sbc_cell(). Re-runs this SBC cell; from this bundle directory:",
    "#   Rscript run_cell.R",
    "# with the EMC2 build named in session.txt (R_LIBS=<its library>). An interrupted run",
    "# resumes from the replicates saved in ../results.",
    "library(EMC2)",
    "cell <- readRDS(\"cell.rds\")",
    "EMC2:::nle_artefact_info(cell$model)",
    "cat(\"BLAS:\", sessionInfo()$BLAS, \"\\n\")",
    sprintf("res <- nn_sbc_cell(cell, trials = %d, replicates = %d, run_control = %s, archive_dir = \"..\"%s)",
            trials, replicates, run_control, extra),
    "print(res$summary)"), file.path(b, "run_cell.R"))
  writeLines(c(utils::capture.output(nle_artefact_info(cell$model)),
               paste("BLAS:", utils::sessionInfo()$BLAS), paste("LAPACK:", utils::sessionInfo()$LAPACK),
               "", utils::capture.output(print(utils::sessionInfo()))), file.path(b, "session.txt"))
  invisible(TRUE)
}

nn_sbc_write <- function(dir, out, trials, replicates, dots, info) {
  cell <- out$cell; reg <- cell$reg
  res_dir <- file.path(dir, "results")
  if (!is.null(out$summary)) {
    utils::write.csv(out$summary, file.path(res_dir, "summary.csv"), row.names = FALSE)
    for (m in c("nn", "control")) if (!is.null(out[[m]])) {
      grDevices::pdf(file.path(res_dir, paste0("ecdf_", m, ".pdf")), width = 10, height = 7)
      tryCatch(plot_sbc_ecdf(out[[m]]), finally = grDevices::dev.off())
    }
  }
  `%||%` <- function(a, b) if (is.null(a)) b else a
  s <- out$summary
  verdict <- if (is.null(s)) "TODO (not run yet)" else {
    sn <- s[s$model == "network", ]
    w <- which.max(abs(sn$std_bias))
    txt <- sprintf("worst abs(std bias) %.3f (%s); ECDFs in envelope %d/%d; KS p %s-%s; coverage %.3f-%.3f",
                   abs(sn$std_bias[w]), sn$parameter[w], sum(sn$in_envelope), nrow(sn),
                   signif(min(sn$ks_p), 2), signif(max(sn$ks_p), 2), min(sn$coverage), max(sn$coverage))
    sc <- s[s$model == "control", ]
    if (nrow(sc)) txt <- paste0(txt, sprintf("; control worst abs(std bias) %.3f, envelope %d/%d",
                                             max(abs(sc$std_bias)), sum(sc$in_envelope), nrow(sc)))
    paste0("**", if (all(abs(sn$std_bias) <= 0.5)) "PASS" else "FAIL",
           "** (provisional, gate abs(std bias) <= 0.5) -- ", txt)
  }
  fmt_f <- function(f) paste(vapply(f, function(z) paste(deparse(z), collapse = ""), ""), collapse = ", ")
  des <- cell$design
  tab <- function(d) c(paste0("| ", paste(names(d), collapse = " | "), " |"),
                       paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"),
                       apply(format(d, digits = 3), 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
  readme <- c(
    sprintf("# SBC: %s vs %s (%s)", reg$label, cell$control_label, info$topic %||% "neural-likelihood cell"), "",
    "| | |", "|---|---|",
    sprintf("| **Date** | %s |", Sys.Date()),
    sprintf("| **Branch / commit** | %s (EMC2 %s @ %s) |",
            if (is.null(info$branch) && is.null(info$commit)) "TODO `<branch>` @ `<hash>`" else
              sprintf("`%s` @ `%s`", info$branch %||% "?", info$commit %||% "?"),
            as.character(utils::packageVersion("EMC2")), find.package("EMC2")),
    sprintf("| **Models** | %s (neural likelihood, %s, sha256 %s) and control %s; single subject |",
            reg$label, reg$kind, reg$sha256, cell$control_label),
    sprintf("| **Why** | %s |", info$why %||% "TODO"),
    sprintf("| **Where run** | %s |", info$where %||% paste(Sys.info()[["nodename"]], "(TODO: cores, jobs)")),
    sprintf("| **Verdict** | %s |", verdict),
    "| **Report** | none |", "",
    "## Design",
    sprintf("Factors: %s; Rlevels: %s; formula: %s; constants (sampled scale): %s. %d trials per replicate.",
            paste(names(des$Ffactors), vapply(des$Ffactors, paste, "", collapse = "/"), sep = " = ", collapse = "; "),
            paste(des$Rlevels, collapse = "/"), fmt_f(des$Flist),
            if (length(des$constants)) paste(names(des$constants), signif(des$constants, 4), sep = " = ", collapse = ", ") else "none",
            trials), "",
    "## Priors",
    "Independent normal on the sampled scale (built by `nn_cell()`); the same prior for the control.", "",
    tab(data.frame(parameter = names(cell$mean), mean = cell$mean, sd = cell$sd)), "",
    sprintf("Natural-scale range at +-%g sd vs the training region:", cell$k), "",
    tab(nn_format_box(cell$box)), "",
    "## Conditions",
    sprintf("One cell: %d replicates x %d trials, network%s. Seeds: run_sbc() defaults.", replicates, trials,
            if (!is.null(out$control) || is.null(s)) " and control" else ""), "",
    "## Fitting",
    if (length(dots)) paste0("run_sbc() arguments: ", paste(names(dots), vapply(dots, function(v)
      paste(deparse(v, width.cutoff = 500L), collapse = ""), ""), sep = " = ", collapse = "; "), ".") else
        "run_sbc() defaults.", "",
    "## Result summary",
    if (is.null(s)) "TODO (bundle prepared, not run yet)." else tab(s), "",
    "## Files",
    "`bundle/`: `cell.rds` (design and prior of both cells), `run_cell.R` (re-runs the cell), `session.txt`",
    "(EMC2 build, BLAS, artefact line, sessionInfo). `results/`: `sbc_nn.RData`, `sbc_control.RData`",
    "(run_sbc() output: SBC, prior_alpha), `summary.csv`, `ecdf_nn.pdf`, `ecdf_control.pdf`.", "",
    "## Not covered",
    "Hierarchical models; parameter regions outside this prior (the rest of the training region); other",
    "designs (factors, covariates); censored or truncated data (refused for neural likelihoods);",
    "identifiability (SBC passes by construction for an unidentified parameter: add parameter recovery",
    "and likelihood profiles).")
  writeLines(readme, file.path(dir, "README.md"))
  invisible(TRUE)
}
