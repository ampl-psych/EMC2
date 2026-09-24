# Neural-likelihood models: the registration contract.
#
# register_nn_model() turns a trained network (an .rds bundle or a JSON model
# card) into an ordinary EMC2 model function. Everything about the network's
# interface is read from the card and checked when it is registered:
#   context_names        the network's parameter inputs, in the network's order
#   context_transforms   the scale each input was trained on (identity/log/probit),
#                        keyed by name; a name listed here but absent from
#                        context_names was fixed when the network was trained
#   bounds_natural       training box on the natural scale, in the order of
#                        names(context_transforms)
#   bounds_sampled       the same box on the network's (sampled) scale, in
#                        context_names order; the two must agree, which also
#                        catches a card whose context_names were permuted
# The parameter types, transforms, defaults and simulator are the analytic
# twin's (the card's `model`: DDM, RDM, ...), so an un-sampled parameter means
# exactly what it means in the analytic model (the defaults trap). A network
# with no analytic model gets its defaults from `p_types` (transforms follow
# the card; no simulator) or from a hand-written model function as `twin`.
#
# Registry: cards are read once per session, keyed by normalised path (re-read
# when the file changes); compiled evaluators are rebuilt lazily because
# external pointers do not survive save/load. A model list carries only the
# registration (path or shipped name, sha256, context spec), never the weights,
# and refuses to evaluate if the file no longer has the sha256 it was
# registered with.

# Per-session cache: "card:<path>" -> list(sig, sha256, card);
# "ptr:<kind>:<sha256>" -> evaluator (external pointer, or the member card for
# regression nets).
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

nn_kinds <- c("flow_joint", "flow_race", "regression_joint", "mlp_joint")
nn_mlp_kinds <- c("regression_joint", "mlp_joint")
nn_transform_codes <- c(identity = 0L, log = 1L, probit = 2L)
# card `model` field -> analytic twin (an EMC2 model function)
nn_twin_names <- c(DDM = "DDM", RDM = "RDM", LNR = "LNR", LBA = "LBA")
# network inputs a regression net may take besides its parameter contexts
nn_data_inputs <- c("rt", "log_rt", "R")

nn_tf <- function(x, tf) switch(tf, identity = x, log = log(x), probit = stats::qnorm(x),
                                stop("Unknown context transform '", tf, "'"))
nn_tf_inv <- function(x, tf) switch(tf, identity = x, log = exp(x), probit = stats::pnorm(x),
                                    stop("Unknown context transform '", tf, "'"))

# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

# `path` -> list(path, artefact). A bare name that is not a file but is a
# shipped artefact resolves to inst/extdata/flownn/<name>.rds; shipped files
# are checked against their MANIFEST pin however they are named.
nn_resolve <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path))
    stop("path must be a single file path or the name of a shipped artefact")
  shipped <- sub("\\.rds$", "", names(nle_manifest()))
  # only an existing regular file is a file: a directory that happens to carry
  # an artefact's name (e.g. an output folder "rdm_small") must not shadow it
  if (!file.exists(path) || dir.exists(path)) {
    if (path %in% shipped) return(nn_shipped(path))
    stop("No neural-likelihood artefact or file '", path, "'; shipped artefacts: ",
         paste(shipped, collapse = ", "))
  }
  path <- normalizePath(path, mustWork = TRUE)
  art <- NULL
  if (identical(dirname(path), normalizePath(nle_dir())) &&
      sub("\\.rds$", "", basename(path)) %in% shipped)
    art <- sub("\\.rds$", "", basename(path))
  list(path = path, artefact = art)
}

nn_shipped <- function(name)
  list(path = normalizePath(file.path(nle_dir(), paste0(name, ".rds")), mustWork = TRUE), artefact = name)

# Read and normalise a card, once per session per file version.
nn_card <- function(res) {
  path <- res$path
  fi <- file.info(path)
  if (is.na(fi$size)) stop("Neural-likelihood file ", path, " does not exist (any more)")
  sig <- paste(fi$size, as.numeric(fi$mtime))
  key <- paste0("card:", path)
  hit <- nle_cache[[key]]
  if (!is.null(hit) && identical(hit$sig, sig)) return(hit)
  sha <- nle_file_hash(path)
  if (!is.null(res$artefact)) {
    pinned <- nle_manifest()[[paste0(res$artefact, ".rds")]]
    if (is.null(pinned) || !identical(pinned, sha))
      stop("Artefact ", res$artefact, ".rds does not match the sha256 pinned in inst/extdata/flownn/MANIFEST")
  }
  raw <- if (grepl("\\.rds$", path, ignore.case = TRUE)) readRDS(path) else
    if (grepl("\\.json$", path, ignore.case = TRUE)) nn_read_json(path) else
      stop("A neural-likelihood file must be an .rds bundle or a .json model card: ", path)
  hit <- list(sig = sig, sha256 = sha, card = nn_normalise_card(raw))
  nle_cache[[key]] <- hit
  hit
}

nn_read_json <- function(path) {
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("Reading a .json model card needs the 'jsonlite' package (install it, or use an .rds bundle)")
  jsonlite::fromJSON(path, simplifyDataFrame = FALSE)
}

# card -> list(members, top, kind). Weight layout as the NLE port's JSON
# loaders produce it (W is n_in x n_out), so .rds and .json cards are the same.
nn_normalise_card <- function(raw) {
  fix_mlp <- function(mlp) {
    mlp$layers <- lapply(mlp$layers, function(l) list(W = as.matrix(l$W), b = as.numeric(l$b)))
    if (isTRUE(mlp$use_norm))
      mlp$norms <- lapply(mlp$norms, function(n)
        list(scale = as.numeric(n$scale), bias = as.numeric(n$bias), eps = as.numeric(n$eps)))
    mlp
  }
  members <- if (!is.null(raw$members)) raw$members else list(raw)
  members <- lapply(members, function(m) {
    for (f in c("mlp", "flow_mlp", "classifier_mlp")) if (!is.null(m[[f]])) m[[f]] <- fix_mlp(m[[f]])
    m
  })
  kinds <- vapply(members, nn_infer_kind, "")
  if (length(unique(kinds)) != 1L) stop("Ensemble members hold different kinds of network")
  declared <- if (!is.null(raw$kind)) raw$kind else members[[1]]$kind
  if (!is.null(declared) && !identical(declared, kinds[1]))
    stop("Card declares kind '", declared, "' but holds a ", kinds[1], " network")
  members <- lapply(members, nn_fill_sampled_box)
  m1 <- members[[1]]
  for (m in members[-1])
    for (f in c("context_names", "context_transforms", "bounds_natural", "bounds_sampled"))
      if (!identical(m[[f]], m1[[f]])) stop("Ensemble members disagree on ", f)
  list(members = members, top = if (!is.null(raw$members)) raw[names(raw) != "members"] else NULL,
       kind = kinds[1])
}

# bounds_sampled from bounds_natural when the card has only the latter (the
# compiled evaluators read the sampled box).
nn_fill_sampled_box <- function(m) {
  bn <- m$bounds_natural; tr <- m$context_transforms; ctx <- m$context_names
  if (!is.null(m$bounds_sampled) || is.null(bn) || is.null(tr) || is.null(ctx) ||
      length(bn$lower) != length(tr) || !all(ctx %in% names(tr))) return(m)
  j <- match(ctx, names(tr))
  tfs <- vapply(ctx, function(p) as.character(tr[[p]]), "", USE.NAMES = FALSE)
  m$bounds_sampled <- list(lower = unname(mapply(nn_tf, bn$lower[j], tfs)),
                           upper = unname(mapply(nn_tf, bn$upper[j], tfs)))
  m
}

nn_infer_kind <- function(m) {
  if (!is.null(m$flow_mlp) && !is.null(m$classifier_mlp)) "flow_joint"
  else if (!is.null(m$mlp) && !is.null(m$spline)) "flow_race"
  else if (!is.null(m$mlp)) if (identical(m$kind, "mlp_joint")) "mlp_joint" else "regression_joint"
  else stop("Cannot tell what kind of network this card holds (neither flow_mlp + classifier_mlp nor mlp)")
}

# The evaluator for a registration, (re)built on demand. Shipped artefacts are
# re-resolved by name so a saved design survives a moved package library.
nn_ptr <- function(reg) {
  key <- paste0("ptr:", reg$kind, ":", reg$sha256)
  ptr <- nle_cache[[key]]
  valid <- !is.null(ptr) && switch(reg$kind, flow_joint = ddm_ptr_valid(ptr),
                                   flow_race = flow_ptr_valid(ptr), mlp_lik_valid(ptr))
  if (!valid) {
    res <- if (!is.null(reg$artefact)) nn_shipped(reg$artefact) else list(path = reg$path)
    hit <- nn_card(res)
    if (!identical(hit$sha256, reg$sha256))
      stop("Neural-likelihood file ", res$path, " has changed since the model was registered ",
           "(sha256 ", reg$sha256, ", now ", hit$sha256, "); re-register it and re-make the design")
    mem <- hit$card$members
    ptr <- switch(reg$kind, flow_joint = ddm_build_ensemble(mem), flow_race = flow_build(mem[[1]]),
                  mlp_lik_build(mem[[1]], reg$lower_s, reg$upper_s, reg$oob))
    nle_cache[[key]] <- ptr
  }
  ptr
}

# ---------------------------------------------------------------------------
# Registration
# ---------------------------------------------------------------------------

#' Register a Neural-Network Likelihood as an EMC2 Model
#'
#' Turns a trained neural likelihood (a normalizing flow or a regression net
#' for the log-likelihood) into an ordinary EMC2 model function, to be passed
#' to [design()] like [DDM] or [RDM]. Everything about the network's interface
#' (input order, input scales, training region) is read from its model card
#' and checked when it is registered; nothing is assumed.
#'
#' @details
#'
#' **What comes from where.** The model's parameter types, transforms,
#' defaults and data-generating function are those of the *analytic twin*, the
#' EMC2 model the network approximates (the card's `model` field: `"DDM"`,
#' `"RDM"`, `"LNR"` or `"LBA"`, or `twin`). An un-sampled parameter therefore
#' takes the same default as in the analytic model (for a network without an
#' analytic model, see below); when that default lies
#' outside the network's training region (e.g. `sv = 0` for a DDM flow
#' trained on `sv >= 0.01`) the model refuses it, at [design()] time and again
#' when the likelihood is evaluated, instead of silently changing the model.
#' The likelihood, the bounds (`bound$minmax` = the training box) and the
#' refusals come from the card. Parameter vectors outside the box are rejected
#' (likelihood floored at `min_ll`) through `bound` and inside the evaluator.
#'
#' **Model card.** The fields read are `context_names` (the network's
#' parameter inputs in order), `context_transforms` (per name: `"identity"`,
#' `"log"` or `"probit"`, the scale the network was trained on),
#' `bounds_natural` and/or `bounds_sampled` (the training box, `lower` and
#' `upper`; natural scale in the order of `names(context_transforms)`, sampled
#' scale in `context_names` order; if both are present they must agree), the
#' weights (`flow_mlp` + `classifier_mlp` + `spline` + `scaler` for
#' `"flow_joint"`, `mlp` + `spline` + `scaler` for `"flow_race"`, `mlp` for
#' `"regression_joint"` and `"mlp_joint"`), and optionally `model`, `kind`, `ll_floor_log`,
#' `context_encoding`. A name in `context_transforms` that is not in
#' `context_names` was fixed when the network was trained; its value must be
#' given in `exceptions` unless the card implies it (`context_encoding =
#' "st0zero"` implies `st0 = 0`). `.rds` bundles may hold an ensemble
#' (`members`, `"flow_joint"` only); `.json` cards need the \pkg{jsonlite}
#' package.
#'
#' **Kinds.**
#' * `"flow_joint"`: a spline flow for rt given parameters and response plus a
#'   classifier for the response; a two-response joint model like [DDM].
#'   `dfun`/`pfun` are the joint density and the defective CDF.
#' * `"flow_race"`: a spline flow for one accumulator's finishing time; a race
#'   model like [RDM] (`dfun`/`pfun` per accumulator; the survivor function
#'   of the losing accumulators needs the CDF, so `cdf = FALSE` is refused).
#' * `"regression_joint"`: an MLP whose single output is the joint log
#'   density log p(rt, R | parameters) of a two-response model. Its card must
#'   declare `input_layout`, the network's inputs in order, each a context
#'   name (entered on its sampled scale) or one of `"rt"`, `"log_rt"`, `"R"`;
#'   optionally `scaler` (`mean`, `scale`; one entry per input, standardising
#'   the whole input vector), `response_values` (the input values coding
#'   responses 1 and 2; default `c(1, 2)`) and `output_scaler` (`mean`,
#'   `scale`: output = raw * scale + mean). It has no CDF. Hidden layers are
#'   GELU (tanh approximation) or tanh; rows outside the training box have log
#'   density `-Inf`.
#' * `"mlp_joint"`: the same evaluator for a likelihood approximation network
#'   (LAN, e.g. the published HSSM networks) converted from ONNX. The card
#'   is written by `inst/scripts/onnx_to_card.py` (in the installed package:
#'   `system.file("scripts", "onnx_to_card.py", package = "EMC2")`); EMC2 does
#'   not run ONNX, a LAN being a plain MLP whose `.onnx` file only holds the
#'   weights. As for `"regression_joint"` the card declares `input_layout` and
#'   `response_values` (HSSM: `c(-1, 1)` for responses 1, 2), and in addition
#'   `ll_floor_log`, the log density that rows outside the training box
#'   return (the floor of the labels the network was trained on) instead of
#'   `-Inf`. A LAN is in its own parameterisation, so give `p_types` (the
#'   defaults for its parameters as it takes them) and read the mapping to
#'   the corresponding analytic model from the card; for HSSM's
#'   `ddm_uniform_st` (inputs `v, a, z, t, st`) it is `a_DDM = 2 * a`,
#'   `Z = z`, `t0 = t - st`, `st0 = 2 * st`, `sv = SZ = 0`, `s = 1`.
#'   **Loading a LAN is not the same as being able to infer with it.** Its
#'   calibration is unknown until checked: an approximation error of the
#'   log density that varies over the parameter space biases the posterior
#'   (for `ddm_uniform_st` the NLE project measured a mean score of -0.122 on
#'   log `a`, where a calibrated likelihood gives 0), and the box, the floor and any
#'   rt range the network was trained on are not enforced beyond what the card
#'   states. Run parameter recovery and a simulation-based calibration cell
#'   for the LAN in your own design before using it for inference.
#'
#' **Evaluation.** Networks whose `pre`
#' is `NULL` or a parameter name run inside EMC2's compiled likelihood (the
#' model's `c_name` is `"NN"`): for every particle the parameters are mapped,
#' passed to the network and the log-likelihood summed in C++, so fits need no
#' R code per particle and can use
#' `options(emc.ll_backend = "multithreaded", emc.n_threads = )`. A `pre`
#' function takes the R path (the model's `dfun`/`pfun` for each particle). Both paths give the same log-likelihood
#' (to rounding); both refuse censored or truncated data.
#'
#' **Networks without an analytic model.** Nothing requires an analytic
#' twin. There are two ways to supply the parameterisation yourself:
#' * `p_types`: a named vector of defaults on the sampled scale, one per
#'   network input (plus any parameter used by `pre` or fixed by
#'   `exceptions`). Each input's transform follows the card's scale (`"log"`
#'   -> `"exp"`, `"probit"` -> `"pnorm"`, `"identity"` -> `"identity"`); any
#'   other parameter gets `"identity"` (override with `design(transform = )`).
#'   The model type follows `kind`. There is no simulator, so [make_data()]
#'   refuses; fit observed data, or use the next route.
#' * `twin = function() list(type = , p_types = , transform = , rfun = )`: a
#'   hand-written model function giving the parameterisation and, if you have
#'   one, a simulator (`rfun(data, pars)` as in any EMC2 model; an optional
#'   `Ttransform` is applied before it, and an optional `bound$minmax` bounds
#'   parameters the network does not take as input).
#'
#' Either way the same checks apply as with an analytic twin: a default
#' outside the training region is refused for an un-sampled parameter, and
#' every parameter must be a network input, used by `pre`, or fixed.
#'
#' **Checks that refuse a registration.** A `pars` that would permute the
#' card's inputs; a card whose `bounds_sampled` disagrees with its transformed
#' `bounds_natural` (e.g. permuted `context_names`); `transforms` that
#' disagree with the card; a network input that is not a parameter of the
#' twin; a twin parameter the network ignores (neither an input, nor consumed
#' by `pre`, nor fixed by `exceptions`); a `pre` that is not a shift of rt.
#'
#' Loading a neural likelihood is not the same as being able to infer with it:
#' check calibration and identifiability (parameter recovery, likelihood
#' profiles, SBC) before drawing conclusions. The validation kit runs the
#' checks against the analytic twin: [nn_cell()] builds matched cells,
#' [nn_total_mass()], [nn_score_bias()] and [nn_posterior_shift()] screen a
#' network in minutes, and [nn_sbc_cell()] runs an SBC cell with its matched
#' analytic control (see `vignette("neural-likelihoods", package = "EMC2")`).
#'
#' **Parallel work on macOS.** Accelerate, the BLAS EMC2 uses there,
#' multithreads large matrix products in a way that crashes forked processes
#' (parallel chains, [run_sbc()]); EMC2 therefore sets
#' `VECLIB_MAXIMUM_THREADS=1` when it is loaded, unless it is already set.
#'
#' @param path An `.rds` bundle or `.json` model card, or the name of an
#'   artefact shipped with EMC2 (`"ddm_cap256w_c4"`, `"ddm_st0zero"`,
#'   `"rdm_small"`; these are checked against their pinned sha256).
#' @param pars EMC2 parameter names for the network's inputs, in the
#'   network's order. Default: the card's `context_names`. Use it to rename
#'   inputs; a name the card also uses must sit at the card's position (a
#'   permutation is refused).
#' @param transforms The scale of each input (`"identity"`, `"log"`,
#'   `"probit"`), in the network's order or named. Default: the card's
#'   `context_transforms`; if both are given they must agree.
#' @param kind `"flow_joint"`, `"flow_race"`, `"regression_joint"` or
#'   `"mlp_joint"`. Default:
#'   inferred from the card; if given it must match.
#' @param cdf Whether the network provides a CDF. Default `TRUE` for flows,
#'   `FALSE` for the MLP kinds. With `cdf = FALSE`, `pfun` refuses.
#' @param exceptions Named natural-scale values of twin parameters that the
#'   network does not take as input because they were fixed in training
#'   (e.g. `c(st0 = 0)`). The model admits exactly that value and refuses
#'   any other.
#' @param pre `NULL`, the name of a parameter subtracted from rt before the
#'   network sees it (`"t0"` for a decision-time network), or a
#'   `function(rt, pars)` returning the network's time; it must be a shift
#'   (no Jacobian is applied). Times <= 0 get zero density and CDF.
#' @param twin The analytic twin: an EMC2 model function supplying
#'   `p_types`, `transform` and `rfun` (and `Ttransform`, applied before
#'   `rfun`). Default: from the card's `model` field. For a network with no
#'   analytic model, a hand-written model function (see Details), or use
#'   `p_types` instead.
#' @param p_types For a network with no analytic model: named defaults on the
#'   sampled scale for the model's parameters (see Details). Give either
#'   `twin` or `p_types`; with `p_types` no twin is used, even if the card
#'   names a model.
#'
#' @return A model function (like [DDM]) whose list carries, besides the usual
#'   elements, `nn`: the registration (`kind`, `sha256`, `pars` and
#'   `context_names` in the network's order, `transforms`, the training box,
#'   `fixed`, `pre`, `cdf`, `ll_floor`, the card's metadata and `twin`: the
#'   name of an EMC2 analytic twin, a hand-written twin function, or `NULL`).
#' @seealso [nn_cell()] and the validation kit; [DDMnn()], [RDMnn()].
#' @examples
#' # A shipped artefact; the card names its analytic twin (RDM)
#' m <- register_nn_model("rdm_small")
#' m()$nn$context_names
#'
#' # A network with no analytic model. Here the shipped RDM flow with its
#' # `model` field removed stands in for your own card: give the defaults
#' # (sampled scale); the transforms follow the card's input scales.
#' card <- readRDS(system.file("extdata", "flownn", "rdm_small.rds", package = "EMC2"))
#' card$model <- NULL
#' path <- tempfile(fileext = ".rds")
#' saveRDS(card, path)
#' m2 <- register_nn_model(path, p_types = c(v = log(1), B = log(1), t0 = log(0.3),
#'                                           s = log(1), A = log(0.3)))
#' m2()$transform$func
#'
#' # The same network with a simulator: a hand-written model function as twin.
#' # rfun(data, pars) is your simulator (here RDM's stands in for it), called
#' # on the parameters after Ttransform.
#' my_model <- function() list(
#'   type = "RACE",
#'   p_types = c(v = log(1), B = log(1), t0 = log(0.3), s = log(1), A = log(0.3)),
#'   transform = list(func = c(v = "exp", B = "exp", t0 = "exp", s = "exp", A = "exp")),
#'   Ttransform = function(pars, dadm) cbind(pars, b = pars[, "B"] + pars[, "A"]),
#'   rfun = RDM()$rfun)
#' m3 <- register_nn_model(path, twin = my_model)
#' @export
register_nn_model <- function(path, pars = NULL, transforms = NULL, kind = NULL,
                              cdf = NULL, exceptions = NULL, pre = NULL, twin = NULL,
                              p_types = NULL) {
  res <- nn_resolve(path)
  hit <- nn_card(res)
  card <- hit$card
  m1 <- card$members[[1]]
  label <- if (!is.null(res$artefact)) res$artefact else basename(res$path)
  what <- paste0("Neural likelihood '", label, "': ")
  if (!is.null(kind)) {
    kind <- match.arg(kind, nn_kinds)
    if (kind != card$kind) stop(what, "the card holds a ", card$kind, " network, not ", kind)
  }
  kind <- card$kind
  if (kind != "flow_joint" && length(card$members) > 1L)
    stop(what, "ensembles are only supported for flow_joint networks")

  # --- inputs: order, names, scales -----------------------------------------
  ctx <- m1$context_names
  if (is.null(ctx)) {
    if (is.null(pars)) stop(what, "the card has no context_names; give the input order with pars")
    ctx <- pars
  }
  ctx <- as.character(ctx)
  k <- length(ctx)
  if (is.null(pars)) pars <- ctx
  pars <- as.character(pars)
  if (length(pars) != k) stop(what, "pars has ", length(pars), " names; the network has ", k, " inputs")
  if (anyDuplicated(pars)) stop(what, "pars has duplicated names")
  moved <- which(pars %in% ctx & pars != ctx)
  if (length(moved))
    stop(what, "pars would permute the network's inputs: '", pars[moved[1]], "' is input ",
         match(pars[moved[1]], ctx), " of the card but position ", moved[1], " of pars. ",
         "The network's order is ", paste(ctx, collapse = ", "), ".")
  tr_card <- m1$context_transforms
  if (!is.null(transforms)) {
    transforms <- stats::setNames(as.character(transforms), names(transforms))
    if (!is.null(names(transforms))) {
      key <- if (setequal(names(transforms), pars)) pars else if (setequal(names(transforms), ctx)) ctx else
        stop(what, "names(transforms) must be the network's inputs (", paste(pars, collapse = ", "), ")")
      transforms <- unname(transforms[key])
    }
    if (length(transforms) != k) stop(what, "transforms needs one entry per network input (", k, ")")
    if (!is.null(tr_card)) {
      card_tf <- vapply(ctx, function(p) as.character(tr_card[[p]]), "")
      bad <- which(card_tf != transforms)
      if (length(bad))
        stop(what, "transforms disagrees with the card's context_transforms for '", ctx[bad[1]],
             "' (", transforms[bad[1]], " vs ", card_tf[bad[1]], ")")
    } else tr_card <- as.list(stats::setNames(transforms, ctx))
  }
  if (is.null(tr_card)) stop(what, "the card has no context_transforms; give them with transforms")
  miss <- setdiff(ctx, names(tr_card))
  if (length(miss)) stop(what, "context_transforms has no entry for input(s) ", paste(miss, collapse = ", "))
  tfs <- vapply(ctx, function(p) as.character(tr_card[[p]]), "", USE.NAMES = FALSE)
  if (!all(tfs %in% names(nn_transform_codes)))
    stop(what, "unsupported context transform(s) ", paste(setdiff(tfs, names(nn_transform_codes)), collapse = ", "))

  # --- training box -----------------------------------------------------------
  box <- nn_box(m1, ctx, tr_card, tfs, what)
  if (!(kind %in% nn_mlp_kinds) && is.null(m1$bounds_sampled))   # read by the compiled evaluators
    stop(what, "a flow card must carry bounds_sampled (or bounds_natural with context_transforms)")

  # --- parameterisation: the analytic twin, or the user's own defaults -----------
  want_type <- if (kind == "flow_race") "RACE" else "DDM"
  if (!is.null(p_types)) {
    if (!is.null(twin))
      stop(what, "give either twin (a model supplies the defaults) or p_types (your own defaults), not both")
    tl <- nn_own_parameterisation(p_types, pars, tfs, want_type, what)
  } else {
    if (is.null(twin) && !is.null(m1$model) && m1$model %in% names(nn_twin_names))
      twin <- get(nn_twin_names[[m1$model]], envir = asNamespace("EMC2"))
    if (is.null(twin))
      stop(what, "no analytic twin (card model '", if (is.null(m1$model)) "" else m1$model, "'). ",
           "Give the defaults with p_types = c(<parameter> = <sampled-scale default>, ...) ",
           "(transforms follow the card; no simulator), or pass twin = <a model function> ",
           "with p_types, transform and, to simulate, rfun")
    if (!is.function(twin)) stop(what, "twin must be an EMC2 model function (e.g. DDM)")
    tl <- twin()
  }
  if (!identical(tl$type, want_type))
    stop(what, "a ", kind, " network needs a model of type ", want_type, "; the twin's type is ",
         if (is.null(tl$type)) "missing" else tl$type)
  pt <- tl$p_types
  if (!is.numeric(pt) || is.null(names(pt))) stop(what, "the twin has no named p_types")
  notp <- setdiff(pars, names(pt))
  if (length(notp)) stop(what, "network input(s) ", paste(notp, collapse = ", "),
                         " are not parameters of the model (", paste(names(pt), collapse = ", "), ")")
  tr_filled <- fill_transform(NULL, function() list(p_types = pt, transform = tl$transform))
  defaults <- nn_natural(pt, tr_filled)

  # --- parameters that are not network inputs ------------------------------------
  dropped <- setdiff(names(tr_card), ctx)          # fixed when the network was trained
  if (!is.null(exceptions)) {
    if (is.null(names(exceptions)) || !is.numeric(exceptions)) stop(what, "exceptions must be a named numeric vector")
    bad <- intersect(names(exceptions), pars)
    if (length(bad)) stop(what, "'", bad[1], "' is an input of the network and cannot be fixed by exceptions; ",
                          "fix it inside the training region with constants in design() instead")
    bad <- setdiff(names(exceptions), names(pt))
    if (length(bad)) stop(what, "exceptions names non-parameter(s) ", paste(bad, collapse = ", "))
  }
  fixed <- exceptions
  for (p in setdiff(dropped, names(fixed))) {
    if (identical(m1$context_encoding, "st0zero") && p == "st0") fixed <- c(fixed, st0 = 0)
    else stop(what, "the network was trained without '", p, "' (it is in context_transforms but not ",
              "context_names); declare the value it was fixed at with exceptions = c(", p, " = <value>)")
  }
  if (is.null(fixed)) fixed <- stats::setNames(numeric(0), character(0))
  pre_pars <- nn_check_pre(pre, pars, defaults, box, what)
  ignored <- setdiff(names(pt), c(pars, names(fixed), pre_pars))
  if (length(ignored))
    stop(what, "the model's parameter(s) ", paste(ignored, collapse = ", "), " are not network inputs, ",
         "not used by pre and not fixed by exceptions: the network would ignore them")

  # --- CDF --------------------------------------------------------------------------
  is_mlp <- kind %in% nn_mlp_kinds
  if (is.null(cdf)) cdf <- !is_mlp
  if (!is.logical(cdf) || length(cdf) != 1L || is.na(cdf)) stop(what, "cdf must be TRUE or FALSE")
  if (is_mlp && cdf) stop(what, "a ", kind, " network has no CDF (cdf must be FALSE)")
  if (kind == "flow_race" && !cdf)
    stop(what, "a race needs the survivor function (1 - CDF) of the losing accumulators, ",
         "so a flow_race network without a CDF cannot be assembled into a race")
  if (is_mlp) nn_check_regression(m1, ctx, kind, what)
  oob <- if (kind == "mlp_joint") m1$ll_floor_log else -Inf

  # --- bounds -------------------------------------------------------------------------
  minmax <- vapply(names(pt), function(p) {
    if (p %in% pars) { j <- match(p, pars); c(box$lower[j], box$upper[j]) }
    else if (p %in% names(fixed)) rep(fixed[[p]], 2)
    else if (!is.null(tl$bound$minmax) && p %in% colnames(tl$bound$minmax)) tl$bound$minmax[, p]
    else c(-Inf, Inf)
  }, numeric(2))
  lower <- stats::setNames(box$lower, pars); upper <- stats::setNames(box$upper, pars)
  refuse_default <- pars[!(defaults[pars] > lower & defaults[pars] < upper)]

  card_meta <- m1[setdiff(names(m1), c("mlp", "flow_mlp", "classifier_mlp", "scaler",
                                        "classifier_scaler", "spline"))]
  card_meta <- c(card_meta, card$top[setdiff(names(card$top), names(card_meta))])
  reg <- list(
    kind = kind, label = label, artefact = res$artefact, path = res$path, sha256 = hit$sha256,
    pars = pars, context_names = ctx, transforms = tfs, transform_codes = unname(nn_transform_codes[tfs]),
    lower = lower, upper = upper, lower_s = box$lower_s, upper_s = box$upper_s,
    fixed = fixed, pre = pre, pre_pars = pre_pars, cdf = cdf,
    defaults = defaults, refuse_default = refuse_default, tr_filled = tr_filled,
    ll_floor = m1$ll_floor_log, oob = oob, n_members = length(card$members),
    input_layout = m1$input_layout, card = card_meta,
    # the analytic twin, the default control of the validation kit
    # (R/nn_validate.R): the name of an EMC2 model, a hand-written model
    # function, or NULL (p_types)
    twin = nn_twin_ref(twin))
  reg$native <- nn_native_spec(reg)
  class(reg) <- "emc_nn"
  nn_model_function(nn_model_list(reg, tl, minmax))
}

# An EMC2 analytic model is kept by name (small, and resolved in the installed
# package); a hand-written twin as the function itself.
nn_twin_ref <- function(twin) {
  if (is.null(twin)) return(NULL)
  for (nm in nn_twin_names)
    if (identical(twin, get(nm, envir = asNamespace("EMC2")))) return(nm)
  twin
}

# The parameterisation of a network with no analytic model: the user's
# defaults, transforms from the card's input scales, no simulator.
nn_own_parameterisation <- function(p_types, pars, tfs, type, what) {
  if (!is.numeric(p_types) || is.null(names(p_types)) || any(!nzchar(names(p_types))) ||
      anyDuplicated(names(p_types)) || anyNA(p_types))
    stop(what, "p_types must be a named numeric vector of defaults on the sampled scale")
  miss <- setdiff(pars, names(p_types))
  if (length(miss)) stop(what, "p_types has no default for network input(s) ", paste(miss, collapse = ", "))
  func <- stats::setNames(rep("identity", length(p_types)), names(p_types))
  func[pars] <- unname(c(identity = "identity", log = "exp", probit = "pnorm")[tfs])
  list(type = type, p_types = p_types, transform = list(func = func))
}

# Training box: natural scale per input (network order), sampled scale per
# input. Cross-checks the card's two boxes when both are present.
nn_box <- function(m, ctx, tr_card, tfs, what) {
  bn <- m$bounds_natural; bs <- m$bounds_sampled
  if (is.null(bn) && is.null(bs)) stop(what, "the card has no training box (bounds_natural / bounds_sampled)")
  k <- length(ctx)
  if (!is.null(bn)) {
    if (length(bn$lower) != length(tr_card) || length(bn$upper) != length(tr_card))
      stop(what, "bounds_natural needs one entry per context_transforms entry (", length(tr_card), ")")
    j <- match(ctx, names(tr_card))
    lo <- bn$lower[j]; hi <- bn$upper[j]
    lo_s <- mapply(nn_tf, lo, tfs); hi_s <- mapply(nn_tf, hi, tfs)
  }
  if (!is.null(bs)) {
    if (length(bs$lower) != k || length(bs$upper) != k)
      stop(what, "bounds_sampled needs one entry per network input (", k, ")")
    if (is.null(bn)) {
      lo_s <- bs$lower; hi_s <- bs$upper
      lo <- mapply(nn_tf_inv, lo_s, tfs); hi <- mapply(nn_tf_inv, hi_s, tfs)
    } else {
      agree <- function(a, b) (is.infinite(a) & a == b) | abs(a - b) <= 1e-8 * pmax(1, abs(b))
      bad <- which(!(agree(lo_s, bs$lower) & agree(hi_s, bs$upper)))
      if (length(bad))
        stop(what, "the card's bounds_sampled disagree with its transformed bounds_natural for input ",
             bad[1], " ('", ctx[bad[1]], "'): its context_names may be permuted or mislabelled")
      lo_s <- bs$lower; hi_s <- bs$upper
    }
  }
  if (any(!(lo < hi))) stop(what, "empty training box for ", paste(ctx[!(lo < hi)], collapse = ", "))
  list(lower = unname(lo), upper = unname(hi), lower_s = unname(lo_s), upper_s = unname(hi_s))
}

# Natural-scale values of sampled-scale parameters; `tr` is a transform
# filled (fill_transform) for all of the model's parameters.
nn_natural <- function(x, tr) {
  x <- x[!is.na(x)]
  if (!length(x)) return(x)
  out <- do_transform(matrix(x, 1L, dimnames = list(NULL, names(x))), tr)
  stats::setNames(as.numeric(out[1L, ]), colnames(out))
}

# pre: NULL, a parameter name, or a shift function(rt, pars). Returns the
# parameters it consumes (probed for a function).
nn_check_pre <- function(pre, pars, defaults, box, what) {
  if (is.null(pre)) return(character(0))
  if (is.character(pre)) {
    if (length(pre) != 1L || !(pre %in% names(defaults)))
      stop(what, "pre must name one parameter of the model (e.g. \"t0\") or be a function(rt, pars)")
    return(pre)
  }
  if (!is.function(pre)) stop(what, "pre must be NULL, a parameter name or a function(rt, pars)")
  probe <- defaults
  mid <- ifelse(is.finite(box$lower) & is.finite(box$upper), (box$lower + box$upper) / 2,
                ifelse(is.finite(box$lower), box$lower + 1, ifelse(is.finite(box$upper), box$upper - 1, 0)))
  probe[pars] <- mid
  P <- matrix(probe, 2L, length(probe), byrow = TRUE, dimnames = list(NULL, names(probe)))
  rt <- c(1, 1.5)
  y0 <- pre(rt, P)
  if (!is.numeric(y0) || length(y0) != 2L) stop(what, "pre(rt, pars) must return one number per rt")
  if (any(abs(pre(rt + 0.25, P) - y0 - 0.25) > 1e-9))
    stop(what, "pre must shift rt (e.g. rt - pars[, \"t0\"]); no Jacobian is applied")
  used <- character(0)
  for (p in setdiff(names(probe), pars)) {
    P2 <- P; P2[, p] <- P2[, p] + 0.1
    if (any(abs(pre(rt, P2) - y0) > 1e-12)) used <- c(used, p)
  }
  used
}

nn_check_regression <- function(m, ctx, kind, what) {
  lay <- m$input_layout
  if (is.null(lay)) stop(what, "a ", kind, " card must declare input_layout (the network's inputs in order)")
  if (kind == "mlp_joint" && !(is.numeric(m$ll_floor_log) && length(m$ll_floor_log) == 1L && is.finite(m$ll_floor_log)))
    stop(what, "a mlp_joint card must give ll_floor_log, the finite log density returned outside the training box")
  lay <- as.character(lay)
  bad <- setdiff(lay, c(ctx, nn_data_inputs))
  if (length(bad)) stop(what, "input_layout entries ", paste(bad, collapse = ", "),
                        " are neither context_names nor one of ", paste(nn_data_inputs, collapse = ", "))
  if (anyDuplicated(lay) || !all(ctx %in% lay))
    stop(what, "input_layout must list every context name exactly once")
  n_in <- nrow(m$mlp$layers[[1]]$W)
  if (n_in != length(lay)) stop(what, "input_layout has ", length(lay), " entries; the network has ", n_in, " inputs")
  if (ncol(m$mlp$layers[[length(m$mlp$layers)]]$W) != 1L)
    stop(what, "a ", kind, " network must have a single output (the joint log density)")
  if (!is.null(m$scaler) && (length(m$scaler$mean) != n_in || length(m$scaler$scale) != n_in))
    stop(what, "a ", kind, " scaler needs one mean and scale per input (", n_in, ")")
  if (!is.null(m$response_values) && length(m$response_values) != 2L)
    stop(what, "response_values must give the input values for responses 1 and 2")
  if (!is.null(m$output_scaler) && !(length(m$output_scaler$mean) == 1L && length(m$output_scaler$scale) == 1L))
    stop(what, "output_scaler needs one mean and one scale")
  nle_mlp_forward(m$mlp, matrix(0, 1L, n_in))     # refuses an unsupported activation now
  invisible(TRUE)
}

# The model list, built in its own small environment: closures that capture
# the registration function's frame would drag the weights into every saved
# design.
nn_model_list <- function(reg, tl, minmax) {
  twin_rfun <- tl$rfun; twin_T <- tl$Ttransform
  joint <- reg$kind != "flow_race"
  ml <- list(
    type = if (joint) "DDM" else "RACE",
    # "NN": the compiled likelihood pipeline evaluates the network directly
    # (calc_ll_manager passes nn_native_args()); NULL: the R path below
    c_name = if (!is.null(reg$native)) "NN",
    p_types = tl$p_types,
    transform = tl$transform,
    bound = list(minmax = minmax, exception = reg$fixed),
    # The network is fed natural-scale parameters as it was trained (e.g. the
    # DDM flows take raw SZ, not DDM's 2 * SZ * min(Z, 1 - Z)); only refuse
    # values it cannot represent.
    Ttransform = function(pars, dadm) {
      nn_check_pars(pars, reg, dadm)
      pars
    },
    prepare_design = function(formula, constants, Rlevels = NULL, ...)
      nn_prepare_design(reg, formula, constants, Rlevels, joint),
    # The twin's simulator on the twin's parameterisation (exact data)
    rfun = if (is.null(twin_rfun)) function(data = NULL, pars)
      stop("Neural likelihood '", reg$label, "' has no simulator; to simulate, register it with ",
           "twin = <a model function with an rfun>") else
        function(data = NULL, pars) {
          tp <- if (is.null(twin_T)) pars else twin_T(pars, data)
          attr(tp, "ok") <- attr(pars, "ok")
          twin_rfun(data, tp)
        },
    nn = reg)
  if (joint) {
    ml$dfun <- function(rt, R, pars) nn_eval_joint(reg, rt, R, pars, "pdf")
    ml$pfun <- function(rt, R, pars) nn_eval_joint(reg, rt, R, pars, "cdf")
    ml$log_likelihood <- function(pars, dadm, model, min_ll = log(1e-10)) {
      nn_check_data(dadm, reg)
      log_likelihood_ddm(pars = pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  } else {
    ml$dfun <- function(rt, pars) nn_eval_race(reg, rt, pars, "pdf")
    ml$pfun <- function(rt, pars) nn_eval_race(reg, rt, pars, "cdf")
    ml$log_likelihood <- function(pars, dadm, model, min_ll = log(1e-10)) {
      nn_check_data(dadm, reg)
      log_likelihood_race(pars = pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  }
  ml
}

nn_model_function <- function(ml) {
  force(ml)
  function() ml
}

# ---------------------------------------------------------------------------
# Evaluation (native: the compiled likelihood pipeline)
# ---------------------------------------------------------------------------

# What calc_ll()/calc_ll_multithreaded() (type "NN", src/model_NN.h) need
# besides the evaluator; NULL when the model must take the R path (a `pre`
# function is R code). The refusals mirror nn_check_pars().
nn_native_spec <- function(reg) {
  if (is.function(reg$pre)) return(NULL)
  refuse <- reg$defaults[reg$refuse_default]
  list(kind = reg$kind, label = reg$label, pars = reg$pars,
       transform_codes = reg$transform_codes,
       pre = if (is.null(reg$pre)) "" else reg$pre,
       fixed = reg$fixed,
       fixed_msg = vapply(names(reg$fixed), nn_fixed_message, "", reg = reg, USE.NAMES = FALSE),
       refuse = refuse,
       refuse_msg = vapply(names(refuse), function(p) nn_default_message(p, refuse[[p]], reg), "",
                           USE.NAMES = FALSE))
}

# The native spec with the live evaluator (rebuilt if the session lost it).
# `constants`: the dadm's constants; the default refusal applies only to them
# (see nn_refused_defaults()).
nn_native_args <- function(reg, constants = NULL) {
  a <- c(reg$native, list(ptr = nn_ptr(reg)))
  keep <- names(a$refuse) %in% names(constants)
  a$refuse <- a$refuse[keep]
  a$refuse_msg <- a$refuse_msg[keep]
  a
}

# ---------------------------------------------------------------------------
# Evaluation (R path)
# ---------------------------------------------------------------------------

# theta (n x inputs) on the network's scale from natural-scale parameters,
# in the network's order; a missing column is an error, never a reorder.
nn_context <- function(pars, reg) {
  miss <- setdiff(reg$pars, colnames(pars))
  if (length(miss)) stop("Neural likelihood '", reg$label, "' needs parameter(s) ",
                         paste(miss, collapse = ", "), " which the model does not supply")
  theta <- matrix(0, nrow(pars), length(reg$pars))
  for (j in seq_along(reg$pars)) theta[, j] <- nn_tf(pars[, reg$pars[j]], reg$transforms[j])
  theta
}

# The network's time for each trial, and which trials it can evaluate.
nn_time <- function(rt, pars, reg) {
  if (is.null(reg$pre)) return(rt)
  if (is.character(reg$pre)) rt - pars[, reg$pre] else reg$pre(rt, pars)
}

nn_eval_joint <- function(reg, rt, R, pars, what) {
  if (what == "cdf" && !reg$cdf)
    stop("Neural likelihood '", reg$label, "' was registered without a CDF (cdf = FALSE)")
  tn <- nn_time(rt, pars, reg)
  R <- as.integer(R)
  use <- !is.na(tn) & tn > 0
  out <- numeric(length(rt))
  if (!any(use)) return(out)
  if (!all(use)) { pars <- pars[use, , drop = FALSE]; tn <- tn[use]; R <- R[use] }
  theta <- nn_context(pars, reg)
  out[use] <- if (reg$kind %in% nn_mlp_kinds) exp(mlp_lik_eval_cpp(nn_ptr(reg), theta, tn, R))
  else ddm_ens_eval_trials_cpp(nn_ptr(reg), theta, tn, R)[[what]]
  out
}

nn_eval_race <- function(reg, rt, pars, what) {
  tn <- nn_time(rt, pars, reg)
  use <- !is.na(tn) & tn > 0
  out <- numeric(length(rt))
  if (!any(use)) return(out)
  if (!all(use)) { pars <- pars[use, , drop = FALSE]; tn <- tn[use] }
  out[use] <- flow_eval_trials_cpp(nn_ptr(reg), nn_context(pars, reg), tn)[[what]]
  out
}

# ---------------------------------------------------------------------------
# Refusals
# ---------------------------------------------------------------------------

nn_default_message <- function(p, value, reg)
  paste0("Parameter '", p, "' is at the model default (", format(value), "), which is outside this ",
         "neural likelihood's training region [", format(reg$lower[[p]]), ", ", format(reg$upper[[p]]),
         "]. Sample it, or fix it to a value inside the region with constants = c(", p,
         " = <sampled-scale value>).")

nn_fixed_message <- function(p, reg)
  paste0("Parameter '", p, "' is fixed at ", format(reg$fixed[[p]]), " in neural likelihood '", reg$label,
         "' (the network was trained without it); it cannot be sampled or set to another value.")

# Run time (every likelihood evaluation): a backstop for what design() refuses.
nn_check_pars <- function(pars, reg, dadm = NULL) {
  for (p in names(reg$fixed))
    if (p %in% colnames(pars) && any(pars[, p] != reg$fixed[[p]], na.rm = TRUE))
      stop(nn_fixed_message(p, reg))
  for (p in nn_refused_defaults(reg, dadm))
    if (p %in% colnames(pars) && all(pars[, p] == reg$defaults[[p]]))
      stop(nn_default_message(p, reg$defaults[[p]], reg))
  invisible(TRUE)
}

# The parameters whose default outside the training region is refused: those
# a design leaves at their default (its constants). A sampled parameter that
# reaches the default value exactly (a probit input underflowing to 0 at an
# extreme proposal) is outside the region like any other value, and the
# bounds reject it. Without a dadm (e.g. simulated data) every such parameter
# is checked.
nn_refused_defaults <- function(reg, dadm) {
  if (is.null(dadm) || is.null(attr(dadm, "p_names"))) return(reg$refuse_default)
  intersect(reg$refuse_default, names(attr(dadm, "constants")))
}

# Censoring and truncation are not wired for neural likelihoods; both paths
# refuse such data (src/model_NN.h makes the same checks).
nn_check_data <- function(dadm, reg) {
  what <- paste0("Neural likelihood '", reg$label, "': ")
  miss <- dadm[["missingness"]]
  if (!is.null(miss) && any(!is.na(miss)))
    stop(what, "censored data (a non-missing 'missingness' code) is not supported")
  lt <- dadm[["LT"]]
  if (!is.null(lt) && any(lt > 0 & is.finite(lt), na.rm = TRUE))
    stop(what, "truncated data (LT > 0) is not supported")
  ut <- dadm[["UT"]]
  if (!is.null(ut) && any(is.finite(ut)))
    stop(what, "truncated data (finite UT) is not supported")
  invisible(TRUE)
}

# design() time: intercept-only constants (explicit, or the twin's defaults
# design() filled in for unspecified parameters) must lie inside the training
# region; fixed parameters must not be sampled.
nn_prepare_design <- function(reg, formula, constants, Rlevels, joint) {
  if (joint && !is.null(Rlevels) && length(Rlevels) != 2L)
    stop("Neural likelihood '", reg$label, "' is a two-response model; Rlevels has ", length(Rlevels), " levels")
  lhs <- vapply(formula, function(f) as.character(stats::terms(f)[[2]]), "")
  nat <- nn_natural(constants[intersect(names(constants), names(reg$defaults))], reg$tr_filled)
  for (p in names(reg$fixed)) {
    if (p %in% lhs && !(p %in% names(constants))) stop(nn_fixed_message(p, reg))
    if (p %in% names(nat) && !isTRUE(all.equal(nat[[p]], reg$fixed[[p]]))) stop(nn_fixed_message(p, reg))
  }
  for (p in intersect(names(nat), reg$pars))
    if (!(nat[[p]] > reg$lower[[p]] && nat[[p]] < reg$upper[[p]])) {
      if (isTRUE(all.equal(nat[[p]], reg$defaults[[p]]))) stop(nn_default_message(p, nat[[p]], reg))
      stop("Constant ", p, " = ", format(constants[[p]]), " (natural scale ", format(nat[[p]]),
           ") is outside neural likelihood '", reg$label, "''s training region [",
           format(reg$lower[[p]]), ", ", format(reg$upper[[p]]), "]")
    }
  list(formula = formula, constants = constants)
}

# One line naming the artefact actually loaded (call at the top of every
# script): a shipped name, a path, a model function or list, or a registration.
nle_artefact_info <- function(x) {
  reg <- if (inherits(x, "emc_nn")) x else if (is.function(x)) x()$nn else if (is.list(x)) x$nn else
    register_nn_model(x)()$nn
  cm <- reg$card
  step <- if (!is.null(cm$checkpoint_step)) cm$checkpoint_step else
    c(flow = cm$checkpoint_step_flow, classifier = cm$checkpoint_step_classifier)
  cat(sprintf("EMC2 %s @ %s\n  artefact %s: budget=%s checkpoint_step=%s sha256=%s\n",
              as.character(utils::packageVersion("EMC2")), find.package("EMC2"),
              reg$label, if (is.null(cm$budget)) "?" else cm$budget,
              if (is.null(names(step))) paste(step, collapse = "/") else
                paste(names(step), step, sep = "=", collapse = "/"),
              reg$sha256))
  invisible(reg)
}
