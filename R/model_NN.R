# Neural-likelihood (normalizing-flow) variants of the DDM and the RDM:
# DDMnn and RDMnn.
#
# The likelihoods come from rational-quadratic spline flows trained on
# simulated data (Wuth, 2026, MSc thesis, UvA; NLE project). The trained
# artefacts ship as inst/extdata/flownn/*.rds, are pinned by sha256 in
# inst/extdata/flownn/MANIFEST, and are evaluated by compiled code
# (src/flow_ddm.cpp, src/flow_race.cpp, derived from the NLE handover package,
# which defines the models; the batched evaluator core is src/nle_flow.h).
#
# The flows are only valid inside their training region. That box is enforced
# twice: as the model's `bound` (EMC2 rejects proposals outside it via the
# usual `ok` mechanism) and inside the C++ evaluator (out-of-box rows return
# pdf = 0 and cdf = 0).
#
# Everything about parameter layout is read from the artefact: the theta
# builder follows `context_names` / `context_transforms` and the bounds come
# from `bounds_natural`, so nothing here is hardcoded to one artefact.

# Per-session cache. External pointers do not survive save/load, so the
# evaluator is rebuilt lazily whenever the cached pointer is invalid.
nle_cache <- new.env(parent = emptyenv())

nle_dir <- function() system.file("extdata", "flownn", package = "EMC2", mustWork = TRUE)

# sha256 where available (R >= 4.5), md5 otherwise; MANIFEST carries both.
nle_file_hash <- function(path) {
  if (exists("sha256sum", envir = asNamespace("tools"), inherits = FALSE))
    unname(tools::sha256sum(path))
  else unname(tools::md5sum(path))
}

nle_manifest <- function() {
  lines <- readLines(file.path(nle_dir(), "MANIFEST"))
  lines <- lines[nzchar(lines)]
  parts <- strsplit(lines, "[[:space:]]+")
  stats::setNames(lapply(parts, function(p) p[1]), vapply(parts, `[`, "", 2))
}

# Load (and cache) an artefact's metadata; verifies the pinned hash once.
nle_meta <- function(artefact) {
  key <- paste0("meta:", artefact)
  if (!is.null(nle_cache[[key]])) return(nle_cache[[key]])
  file <- paste0(artefact, ".rds")
  path <- file.path(nle_dir(), file)
  if (!file.exists(path))
    stop("No neural-likelihood artefact '", artefact, "'; shipped: ",
         paste(sub("\\.rds$", "", names(nle_manifest())), collapse = ", "))
  pinned <- nle_manifest()[[file]]
  actual <- nle_file_hash(path)
  if (is.null(pinned) || !identical(pinned, actual))
    stop("Artefact ", file, " does not match the sha256 pinned in inst/extdata/flownn/MANIFEST")
  x <- readRDS(path)
  # DDM bundles are (possibly single-member) ensembles; race bundles are one flow
  members <- if (!is.null(x$members)) x$members else list(x)
  m1 <- members[[1]]
  for (m in members[-1])
    if (!identical(m$context_names, m1$context_names))
      stop("Ensemble members of '", artefact, "' disagree on context_names")
  meta <- list(
    artefact = artefact, sha256 = actual, is_ddm = !is.null(x$members),
    n_members = length(members),
    budget = m1$budget,
    checkpoint_step = if (!is.null(m1$checkpoint_step)) m1$checkpoint_step else
      c(flow = m1$checkpoint_step_flow, classifier = m1$checkpoint_step_classifier),
    context_names = m1$context_names,
    context_transforms = m1$context_transforms,
    bounds_natural = m1$bounds_natural,
    ll_floor_log = m1$ll_floor_log)
  nle_cache[[key]] <- meta
  meta
}

# The evaluator external pointer for an artefact, (re)built on demand.
nle_get <- function(artefact) {
  meta <- nle_meta(artefact)
  key <- paste0("ptr:", artefact)
  ptr <- nle_cache[[key]]
  valid <- !is.null(ptr) && (if (meta$is_ddm) ddm_ptr_valid(ptr) else flow_ptr_valid(ptr))
  if (!valid) {
    x <- readRDS(file.path(nle_dir(), paste0(artefact, ".rds")))
    ptr <- if (meta$is_ddm) ddm_build_ensemble(x$members) else flow_build(x)
    nle_cache[[key]] <- ptr
  }
  ptr
}

# One line naming the artefact actually loaded (call at the top of every script)
nle_artefact_info <- function(artefact) {
  meta <- nle_meta(artefact)
  cat(sprintf("EMC2 %s @ %s\n  artefact %s: budget=%s checkpoint_step=%s sha256=%s\n",
              as.character(utils::packageVersion("EMC2")), find.package("EMC2"),
              artefact, meta$budget,
              paste(names(meta$checkpoint_step), meta$checkpoint_step, sep = "=", collapse = "/"),
              meta$sha256))
  invisible(meta)
}

# theta (n x length(context_names)) on the flow's sampled scale from the
# natural-scale parameter matrix. The order and the transforms are the
# artefact's; a missing column is an error, never a silent reorder.
nle_theta <- function(pars, meta) {
  nm <- meta$context_names
  miss <- setdiff(nm, colnames(pars))
  if (length(miss)) stop("Neural likelihood '", meta$artefact, "' needs parameter(s) ",
                         paste(miss, collapse = ", "), " which the model does not supply")
  theta <- matrix(0, nrow(pars), length(nm))
  for (j in seq_along(nm)) {
    x <- pars[, nm[j]]
    theta[, j] <- switch(meta$context_transforms[[nm[j]]],
                         identity = x, log = log(x), probit = stats::qnorm(x),
                         stop("Unknown context transform '", meta$context_transforms[[nm[j]]], "'"))
  }
  theta
}

# Parameters the artefact cannot represent at their model default (e.g. sv = 0
# in the full-DDM flow): using the default is a silent mis-specification, so
# refuse. `exceptions` are the values the artefact does support (st0 = 0 for
# the st0-free build).
nle_default_check <- function(pars, defaults_nat, lower, upper, exceptions) {
  bad <- names(defaults_nat)[!(defaults_nat >= lower[names(defaults_nat)] &
                                 defaults_nat <= upper[names(defaults_nat)]) &
                               !(names(defaults_nat) %in% names(exceptions))]
  for (p in bad)
    if (p %in% colnames(pars) && all(pars[, p] == defaults_nat[[p]]))
      stop("Parameter '", p, "' is at the model default (", format(defaults_nat[[p]]),
           "), which is outside this neural likelihood's training region [",
           format(lower[[p]]), ", ", format(upper[[p]]), "]. Sample it, or fix it to a ",
           "value inside the region with constants = c(", p, " = <sampled-scale value>).")
  invisible(TRUE)
}

# bound = list(minmax, exception) from the artefact's natural-scale box.
# A context absent from context_names but present in context_transforms was
# fixed at 0 when the artefact was built (st0 in the st0-free DDM).
nle_bound <- function(meta, pars) {
  lo <- stats::setNames(meta$bounds_natural$lower, names(meta$context_transforms))
  hi <- stats::setNames(meta$bounds_natural$upper, names(meta$context_transforms))
  zero <- setdiff(names(meta$context_transforms), meta$context_names)
  exc <- stats::setNames(rep(0, length(zero)), zero)
  minmax <- vapply(pars, function(p) c(lo[[p]], hi[[p]]), numeric(2))
  colnames(minmax) <- pars
  list(lower = lo, upper = hi, exceptions = exc,
       bound = list(minmax = minmax, exception = exc))
}

#' The Diffusion Decision Model — Neural-Network Likelihood
#'
#' Variant of [DDM] whose joint likelihood of response and rt is computed by
#' trained neural networks (neural likelihood estimation) instead of the
#' numerical-integral density: a normalizing flow for the rt distribution
#' conditional on parameters and response, and a classifier for the choice
#' probability. Evaluation takes microseconds per trial.
#'
#' @details
#'
#' Parameters, transforms, defaults and bounds' meaning match [DDM]; the
#' natural-scale ranges below are the flow's **training region** (a bounding
#' hyper-box, read from the artefact): parameter vectors outside it are
#' rejected (likelihood floored at `min_ll`), both through the model's `bound`
#' and inside the compiled evaluator.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default**    | **Interpretation**            |
#' |-----------|-----------|---------------|------------|---------------------------|
#' | *v*       | -         | \[-6, 6\]       | 1          | Mean drift rate |
#' | *a*       | log       | \[0.05, 5\]     | log(1)     | Boundary separation |
#' | *t0*      | log       | \[0.05, 1\]     | log(0)     | Non-decision time |
#' | *s*       | log       | \[0.1, 2\]      | log(1)     | Within-trial drift SD |
#' | *Z*       | probit    | \[0.1, 0.9\]    | qnorm(0.5) | Relative start point |
#' | *SZ*      | probit    | \[0.05, 0.99\]  | qnorm(0)   | Start-point variability |
#' | *sv*      | log       | \[0.01, 4\]     | log(0)     | Between-trial drift SD |
#' | *st0*     | log       | \[0.01, 0.49\]  | log(0)     | Non-decision variability |
#'
#' The defaults are **identical to [DDM]'s**, including `sv = SZ = st0 = 0`. The
#' full-DDM artefact cannot represent those zeros, so an un-sampled `sv`, `SZ`
#' or `st0` is refused with an error instead of silently changing the model:
#' sample the parameter, or fix it inside the box (e.g.
#' `constants = c(sv = log(0.5))`). The `"ddm_st0zero"` artefact was trained
#' without `st0` and supports `st0 = 0` exactly.
#'
#' Unlike [DDM], `SZ` here is the raw proportion the flow was trained on, not
#' remapped by `2 * SZ * min(Z, 1 - Z)`; data generation applies that remapping
#' internally so `rfun` matches [DDM]'s parameter meaning exactly. Data are
#' generated with the exact sampler (`WienR::rWDM`).
#'
#' Loading a neural likelihood is not the same as being able to infer with it:
#' check calibration and identifiability (parameter recovery, likelihood
#' profiles, SBC) before drawing conclusions.
#'
#' Wuth, J. (2026). *Likelihood approximation in evidence accumulation
#' models: A comparison of kernel density and neural likelihood methods*
#' (MSc thesis, University of Amsterdam).
#'
#' @param artefact Name of a shipped artefact in `inst/extdata/flownn`:
#'   `"ddm_cap256w_c4"` (default, full 8-parameter DDM) or `"ddm_st0zero"`
#'   (7 parameters, `st0 = 0`).
#' @return A model list with all the necessary functions for EMC2 to sample
#' @export
DDMnn <- function(artefact = "ddm_cap256w_c4") {
  meta <- nle_meta(artefact)
  if (!meta$is_ddm) stop("'", artefact, "' is not a DDM artefact")
  p_types <- c("v" = 1,"a" = log(1),"sv" = log(0),"t0" = log(0),"st0" = log(0),
               "s" = log(1),"Z" = qnorm(0.5),"SZ" = qnorm(0))
  pn <- names(p_types)
  if (!setequal(names(meta$context_transforms), pn))
    stop("Artefact '", artefact, "' has parameters ", paste(names(meta$context_transforms), collapse = ", "),
         "; DDMnn expects ", paste(pn, collapse = ", "))
  bd <- nle_bound(meta, c("v","a","Z","t0","sv","s","SZ","st0"))
  defaults_nat <- c(v = 1, a = 1, sv = 0, t0 = 0, st0 = 0, s = 1, Z = 0.5, SZ = 0)
  dfun <- function(rt, R, pars)
    ddm_ens_eval_trials_cpp(nle_get(artefact), nle_theta(pars, meta), rt, as.integer(R))$pdf
  pfun <- function(rt, R, pars)
    ddm_ens_eval_trials_cpp(nle_get(artefact), nle_theta(pars, meta), rt, as.integer(R))$cdf
  list(
    type="DDM",
    c_name = NULL, # R-path DDM likelihood; the evaluator is compiled
    p_types = p_types,
    transform=list(func=c(v = "identity",a = "exp",sv = "exp",t0 = "exp",
                          st0 = "exp",s = "exp",Z = "pnorm",SZ = "pnorm")),
    bound = bd$bound,
    # The flow is trained on the raw Z/SZ (cf. DDM, which rescales SZ and adds
    # z/sz here); only refuse defaults the artefact cannot represent.
    Ttransform = function(pars, dadm) {
      nle_default_check(pars, defaults_nat, bd$lower, bd$upper, bd$exceptions)
      pars
    },
    # Exact data; apply the DDM's SZ remapping locally so SZ means the same
    # thing as in the analytic model
    rfun=function(data=NULL, pars) {
      pars[,"SZ"] <- 2*pars[,"SZ"]*pmin(pars[,"Z"], 1 - pars[,"Z"])
      rDDM(data$R, pars, attr(pars, "ok"))
    },
    dfun = dfun,
    pfun = pfun,
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10))
      log_likelihood_ddm(pars=pars, dadm = dadm, model = model, min_ll = min_ll),
    nn = list(artefact = artefact, sha256 = meta$sha256,
              context_names = meta$context_names, ll_floor = meta$ll_floor_log)
  )
}

#' The Racing Diffusion Model — Neural-Network Likelihood
#'
#' Variant of [RDM] whose single-accumulator density and survivor functions
#' are computed by a trained normalizing flow (neural likelihood estimation)
#' instead of the analytic Wald race expressions.
#'
#' @details
#'
#' Parameters, transforms, and defaults match [RDM]; the natural-scale ranges
#' are the flow's **training region**, read from the artefact:
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default**    | **Interpretation**            |
#' |-----------|-----------|---------------|------------|---------------------------|
#' | *v*       | log       | \[0.001, 5\]    | log(1)     | Evidence-accumulation rate |
#' | *B*       | log       | \[0.1, 3\]      | log(1)     | Threshold gap (b = B + A) |
#' | *A*       | log       | \[0.0001, 2\]   | log(0)     | Start-point variability   |
#' | *t0*      | log       | \[0.05, 1\]     | log(0)     | Non-decision time         |
#' | *s*       | log       | \[0.1, 2\]      | log(1)     | Within-trial noise        |
#'
#' Unlike [RDM], `A = 0` and `v = 0` are not available (the flow trained on
#' `A >= 1e-4`, `v >= 1e-3`). As for [DDMnn], an un-sampled `A` at its default
#' `0` is refused with an error; fix it inside the box instead, e.g.
#' `constants = c(A = log(1e-3))`. Out-of-box parameter vectors are rejected
#' (likelihood floored at `min_ll`).
#'
#' Data generation (`rfun`) uses the analytic RDM, so simulation-based
#' checks compare flow-based inference against exact data. Loading a neural
#' likelihood is not the same as being able to infer with it; validate before
#' use.
#'
#' Wuth, J. (2026). *Likelihood approximation in evidence accumulation
#' models: A comparison of kernel density and neural likelihood methods*
#' (MSc thesis, University of Amsterdam).
#'
#' @param artefact Name of a shipped artefact in `inst/extdata/flownn`
#'   (default `"rdm_small"`).
#' @return A model list with all the necessary functions for EMC2 to sample
#' @export
RDMnn <- function(artefact = "rdm_small") {
  meta <- nle_meta(artefact)
  if (meta$is_ddm) stop("'", artefact, "' is not a race artefact")
  pn <- c("v", "B", "A", "t0", "s")
  if (!setequal(names(meta$context_transforms), pn))
    stop("Artefact '", artefact, "' has parameters ", paste(names(meta$context_transforms), collapse = ", "),
         "; RDMnn expects ", paste(pn, collapse = ", "))
  bd <- nle_bound(meta, pn)
  defaults_nat <- c(v = 1, B = 1, A = 0, t0 = 0, s = 1)
  dfun <- function(rt, pars)
    flow_eval_trials_cpp(nle_get(artefact), nle_theta(pars, meta), rt)$pdf
  pfun <- function(rt, pars)
    flow_eval_trials_cpp(nle_get(artefact), nle_theta(pars, meta), rt)$cdf
  list(
    type="RACE",
    c_name = NULL, # R-path race likelihood; the evaluator is compiled
    p_types=c("v" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0),"s" = log(1)),
    transform=list(func=c(v = "exp", B = "exp", A = "exp",t0 = "exp", s = "exp")),
    bound = bd$bound,
    Ttransform = function(pars, dadm) {
      nle_default_check(pars, defaults_nat, bd$lower, bd$upper, bd$exceptions)
      cbind(pars, b = pars[,"B"] + pars[,"A"])
    },
    rfun=function(data=NULL,pars) rRDM(data$lR,pars,ok=attr(pars, "ok")),
    dfun = dfun,
    pfun = pfun,
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10))
      log_likelihood_race(pars=pars, dadm = dadm, model = model, min_ll = min_ll),
    nn = list(artefact = artefact, sha256 = meta$sha256,
              context_names = meta$context_names, ll_floor = meta$ll_floor_log)
  )
}
