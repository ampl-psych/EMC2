# Golden vectors: a fixed set of inputs with the outputs the network must
# return, computed once from the trained weights by the exporter's float64
# port and shipped beside the card (<tag>_golden.csv). Every evaluator of the
# card (compiled EMC2, the R reference port, an ONNX graph) is tested against
# the same file. nn_golden_check() runs EMC2's evaluator on the rows of a
# golden file and reports the largest difference per output.

#' Check a Neural Likelihood Against Its Golden Vectors
#'
#' Evaluates the network of a registered neural likelihood (see
#' [register_nn_model()]) on the parameter vectors and times of its golden
#' file and compares the outputs with the values recorded there. The golden
#' file is written by the exporter from the trained weights (the NLE project's
#' `export_flow*.py`, or `golden.py` from a card), so this is the test that
#' the installed evaluator is the network the card describes: a wrong input
#' order, transform, constant or spline setting fails it by a large margin.
#'
#' @details The golden file is a CSV with one row per (parameter set, time):
#' for race flows the columns `param_id`, the inputs in `context_names` order
#' on the network's scale, `rt`, `log_pdf`, `cdf`, `log_sf`; for joint flows
#' (`flow_joint`, `flow_scaled`) the inputs on the exporter's storage scale
#' (the sampled scale; for a `zerobox` input the log, or the probit for `SZ`,
#' of the natural value, `-Inf` at an exact zero), `R`, `rt`, `log_pdf_joint`,
#' `cdf_defective` and `log_p_R`; for a LAN (`mlp_joint`) the inputs on the
#' natural scale, `rt`, `choice` or `R`, and the log density. For a
#' `flow_scaled` card `rt` is the unit-scale time `u`.
#'
#' Shipped artefacts carry their golden files in
#' `system.file("extdata", "flownn", "golden", package = "EMC2")`; a card
#' `<tag>_weights.json` is checked against `<tag>_golden.csv` beside it.
#'
#' @param x A model function from [register_nn_model()], its model list, a
#'   registration, a card path, or a shipped artefact's name.
#' @param golden Path of the golden CSV. Default: the file beside the card,
#'   or the shipped one.
#' @param tol The largest absolute difference accepted, per output
#'   (`log_pdf` and the CDF). The exporters' float64 port and EMC2 agree to
#'   about `1e-12`; `1e-9` leaves room for a different BLAS. A LAN converted
#'   from a float32 ONNX file needs about `1e-4`.
#' @param flow_only For a joint flow: compare the flow's part only (the joint
#'   log density minus `log_p_R`, and the CDF divided by `P(R)`), for a golden
#'   file written with another classifier than the card's.
#' @param error If `TRUE` (the default), a failed check is an error; otherwise
#'   the result is returned with `pass = FALSE`.
#' @param twin,... Passed to [register_nn_model()] when `x` is a path or a name.
#' @return (Invisibly) an `nn_golden_check` object: `label`, `kind`, `n` rows,
#'   `max_abs` (the largest difference per output), `tol`, `pass`. Its print
#'   method is one line, `PASS` or `FAIL` with the largest difference.
#' @examples
#' nn_golden_check("rdm_small")
#' @export
nn_golden_check <- function(x, golden = NULL, tol = 1e-9, flow_only = FALSE, error = TRUE,
                            twin = NULL, ...) {
  reg <- nn_golden_reg(x, twin, ...)
  if (is.null(golden)) golden <- nn_golden_path(reg)
  if (!file.exists(golden)) stop("No golden file ", golden, " for neural likelihood '", reg$label, "'")
  g <- utils::read.csv(golden, stringsAsFactors = FALSE)
  ctx <- reg$context_names
  miss <- setdiff(ctx, names(g))
  if (length(miss)) stop("Golden file ", golden, " lacks the input column(s) ", paste(miss, collapse = ", "))
  d <- switch(reg$kind,
    flow_race = nn_golden_race(reg, g),
    flow_joint = , flow_scaled = nn_golden_joint(reg, g, flow_only),
    regression_joint = , mlp_joint = nn_golden_mlp(reg, g),
    stop("No golden check for kind ", reg$kind))
  pass <- all(is.finite(d)) && all(d <= tol)
  out <- structure(list(label = reg$label, kind = reg$kind, golden = golden, n = nrow(g), max_abs = d,
                        tol = tol, pass = pass), class = "nn_golden_check")
  if (error && !pass) stop(format(out), call. = FALSE)
  invisible(out)
}

#' @export
format.nn_golden_check <- function(x, ...) {
  worst <- which.max(x$max_abs)
  sprintf("%s golden check of '%s' (%s, %d rows): max |diff| %.2e (%s) against %.0e%s",
          if (x$pass) "PASS" else "FAIL", x$label, x$kind, x$n, x$max_abs[worst], names(x$max_abs)[worst],
          x$tol, if (x$pass) "" else paste0("; per output: ", paste(names(x$max_abs),
                                                                   formatC(x$max_abs, format = "e", digits = 2),
                                                                   sep = " ", collapse = ", ")))
}

#' @export
print.nn_golden_check <- function(x, ...) { cat(format(x), "\n", sep = ""); invisible(x) }

nn_golden_reg <- function(x, twin, ...) {
  if (inherits(x, "emc_nn")) return(x)
  if (is.function(x)) { ml <- x(); if (!is.null(ml$nn)) return(ml$nn) }
  if (is.list(x) && !is.null(x$nn)) return(x$nn)
  if (is.character(x) && length(x) == 1L) return(register_nn_model(x, twin = twin, ...)()$nn)
  stop("x must be a registered neural likelihood (model function, model list or registration), a card path ",
       "or a shipped artefact's name")
}

# <tag>_golden.csv beside a card file (tag: the file name without
# _weights.json / .json / .rds); the shipped golden for a shipped artefact.
nn_golden_path <- function(reg) {
  if (!is.null(reg$artefact))
    return(file.path(nle_dir(), "golden", paste0(reg$artefact, "_golden.csv")))
  tag <- sub("(_weights)?\\.(json|rds)$", "", basename(reg$path), ignore.case = TRUE)
  file.path(dirname(reg$path), paste0(tag, "_golden.csv"))
}

# |a - b| with equal infinities (and both NA) counting as agreement, as the
# NLE port's validate_ddm.R does for the st0 = 0 support rule.
nn_golden_diff <- function(a, b) {
  d <- abs(a - b)
  same <- (is.infinite(a) & is.infinite(b) & sign(a) == sign(b)) | (is.na(a) & is.na(b))
  d[same] <- 0
  d
}

# Golden columns of a joint flow -> the network's scale. The exporter stores
# the EMC2 sampled scale: a zerobox input as log (probit for SZ) of the natural
# value with -Inf at zero; other inputs as the network takes them.
nn_golden_theta_joint <- function(reg, g) {
  ctx <- reg$context_names
  th <- matrix(0, nrow(g), length(ctx))
  for (j in seq_along(ctx)) {
    x <- as.numeric(g[[ctx[j]]])
    if (reg$transforms[j] == "zerobox") {
      nat <- if (ctx[j] == "SZ") stats::pnorm(x) else exp(x)
      x <- nn_tf(nat, "zerobox", reg$zerobox_c[j])
    }
    th[, j] <- x
  }
  th
}

nn_golden_joint <- function(reg, g, flow_only) {
  need <- c("R", "rt", "log_pdf_joint", "cdf_defective")
  miss <- setdiff(need, names(g))
  if (length(miss)) stop("Golden file lacks column(s) ", paste(miss, collapse = ", "))
  th <- nn_golden_theta_joint(reg, g)
  ev <- ddm_ens_eval_trials_cpp(nn_ptr(reg), th, as.numeric(g$rt), as.integer(g$R))
  if (flow_only) {
    if (is.null(g$log_p_R)) stop("flow_only needs the golden column log_p_R")
    c(log_pdf_flow = max(nn_golden_diff(ev$log_pdf - log(ev$p_R), g$log_pdf_joint - g$log_p_R)),
      cdf_flow = max(nn_golden_diff(ev$cdf / ev$p_R, g$cdf_defective / exp(g$log_p_R))))
  } else {
    d <- c(log_pdf = max(nn_golden_diff(ev$log_pdf, g$log_pdf_joint)),
           cdf = max(nn_golden_diff(ev$cdf, g$cdf_defective)))
    if (!is.null(g$log_p_R)) d <- c(d, log_p_R = max(nn_golden_diff(log(ev$p_R), g$log_p_R)))
    d
  }
}

nn_golden_race <- function(reg, g) {
  need <- c("rt", "log_pdf", "cdf", "log_sf")
  miss <- setdiff(need, names(g))
  if (length(miss)) stop("Golden file lacks column(s) ", paste(miss, collapse = ", "))
  th <- as.matrix(g[, reg$context_names, drop = FALSE]); storage.mode(th) <- "double"
  # a decision-time card's golden file holds the time the network sees
  ev <- flow_eval_trials_cpp(nn_ptr(reg), th, as.numeric(g$rt))
  c(log_pdf = max(nn_golden_diff(ev$log_pdf, g$log_pdf)),
    cdf = max(nn_golden_diff(ev$cdf, g$cdf)),
    log_sf = max(nn_golden_diff(ev$log_sf, g$log_sf)))
}

nn_golden_mlp <- function(reg, g) {
  lp_col <- intersect(c("log_pdf_joint", "log_pdf", "logp", "logp_ort"), names(g))
  if (!length(lp_col)) stop("Golden file has no log-density column (log_pdf_joint, log_pdf, logp or logp_ort)")
  if (!is.null(g$R)) R <- as.integer(g$R)
  else if (!is.null(g$choice)) {
    rv <- if (is.null(reg$card$response_values)) c(1, 2) else unlist(reg$card$response_values)
    R <- ifelse(g$choice == rv[2], 2L, ifelse(g$choice == rv[1], 1L, NA_integer_))
    if (anyNA(R)) stop("Golden 'choice' values are not the card's response_values")
  } else stop("Golden file needs a response column (R, or choice in the card's response_values)")
  pars <- as.matrix(g[, reg$context_names, drop = FALSE]); storage.mode(pars) <- "double"
  colnames(pars) <- reg$pars
  lp <- log(nn_eval_joint(reg, as.numeric(g$rt), R, pars, "pdf"))
  c(log_pdf = max(nn_golden_diff(lp, as.numeric(g[[lp_col[1]]]))))
}
