# Neural-likelihood models: the registration contract.
#
# register_nn_model() turns a trained network (an .rds bundle or a JSON model
# card) into an ordinary EMC2 model function. Everything about the network's
# interface is read from the card and checked when it is registered:
#   context_names        the network's parameter inputs, in the network's order
#   context_transforms   the scale each input was trained on (identity/log/probit/
#                        zerobox; zerobox = asinh(x / c) with c from zerobox_c),
#                        keyed by name; a name listed here but absent from
#                        context_names was fixed when the network was trained
#   bounds_natural       training box on the natural scale, in the order of
#                        names(context_transforms)
#   bounds_sampled       the same box on the network's (sampled) scale, in
#                        context_names order; the two must agree, which also
#                        catches a card whose context_names were permuted
#   context_scale        "natural" (with no context_transforms): every input
#                        is on the natural scale (transform identity), and the
#                        box is bounds_natural in context_names order
#   input_scaling, affine, pre
#                        race flows only: see nn_check_flow_race()
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

# The hash that pins a file: sha256 where R has it (tools::sha256sum, R >= 4.5),
# md5 otherwise. MANIFEST carries a line per algorithm and file,
# "<hash> <file> <algorithm>" (a line without the third field is sha256), so
# that the shipped artefacts load on every R this package supports; an R
# without sha256 and a MANIFEST without an md5 line is refused by name.
nle_hash_algo <- function()
  if (exists("sha256sum", envir = asNamespace("tools"), inherits = FALSE)) "sha256" else "md5"

nle_file_hash <- function(path) {
  if (nle_hash_algo() == "sha256") unname(tools::sha256sum(path)) else unname(tools::md5sum(path))
}

# All pins: list(file = c(sha256 = , md5 = )).
nle_manifest_all <- function() {
  lines <- readLines(file.path(nle_dir(), "MANIFEST"))
  lines <- lines[nzchar(lines)]
  parts <- strsplit(trimws(lines), "[[:space:]]+")
  out <- list()
  for (p in parts) {
    algo <- if (length(p) >= 3L) p[3] else "sha256"
    out[[p[2]]][[algo]] <- p[1]
  }
  out
}

# The pin for the hash this R computes, keyed by file (NULL where the MANIFEST
# has no line for that algorithm).
nle_manifest <- function() {
  algo <- nle_hash_algo()
  all <- nle_manifest_all()
  stats::setNames(lapply(all, function(h) h[[algo]]), names(all))
}

nn_kinds <- c("flow_joint", "flow_race", "regression_joint", "mlp_joint", "flow_scaled", "flow_cat")
nn_mlp_kinds <- c("regression_joint", "mlp_joint")
# card versions this EMC2 reads (card_version; absent = 0). The format is
# described in vignette("neural-likelihood-cards").
nn_card_versions <- c(0, 1)
nn_transform_codes <- c(identity = 0L, log = 1L, probit = 2L, zerobox = 3L)
# card `model` field -> analytic twin (an EMC2 model function). DDMsc is the
# scaled DDM card (kind flow_scaled), whose twin is the DDM.
nn_twin_names <- c(DDM = "DDM", RDM = "RDM", LNR = "LNR", LBA = "LBA", CRDM = "CRDM", DDMsc = "DDM")
# the DDM parameters a scaled card (flow_scaled) computes its inputs from, in
# the order the compiled evaluator takes them (src/model_NN.h, sc_pars)
nn_scaled_pars <- c("v", "a", "sv", "t0", "st0", "s", "Z", "SZ")
nn_scaled_inputs_names <- c("kappa", "eta", "Z", "SZ")
# network inputs a regression net may take besides its parameter contexts
nn_data_inputs <- c("rt", "log_rt", "R")

# Natural scale <-> network scale. `c`: the zerobox constant (zerobox_c) of the
# input, NA for the other transforms. A zerobox input is exact at x = 0
# (asinh(0 / c) = 0), so an exact zero needs no special case.
nn_tf <- function(x, tf, c = NA_real_) switch(tf, identity = x, log = log(x), probit = stats::qnorm(x),
                                              zerobox = asinh(x / c),
                                              stop("Unknown context transform '", tf, "'"))
nn_tf_inv <- function(x, tf, c = NA_real_) switch(tf, identity = x, log = exp(x), probit = stats::pnorm(x),
                                                  zerobox = c * sinh(x),
                                                  stop("Unknown context transform '", tf, "'"))

# The zerobox constant of each input in `ctx` (NA where the transform is not
# zerobox). A zerobox input without a constant, or a constant on an input that
# is not zerobox, is refused.
nn_zerobox_c <- function(m, ctx, tfs, what = "") {
  zc <- m[["zerobox_c", exact = TRUE]]
  zc <- if (is.null(zc)) stats::setNames(numeric(0), character(0)) else unlist(zc)
  zb <- ctx[tfs == "zerobox"]
  miss <- setdiff(zb, names(zc))
  if (length(miss))
    stop(what, "input(s) ", paste(miss, collapse = ", "), " have the zerobox transform but the card's ",
         "zerobox_c gives no constant for them")
  extra <- setdiff(names(zc), zb)
  if (length(extra))
    stop(what, "zerobox_c gives a constant for ", paste(extra, collapse = ", "),
         ", which the card does not mark \"zerobox\" in context_transforms")
  if (length(zc) && !(is.numeric(zc) && all(is.finite(zc) & zc > 0)))
    stop(what, "zerobox_c constants must be positive numbers")
  out <- rep(NA_real_, length(ctx))
  out[tfs == "zerobox"] <- zc[zb]
  out
}

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

nn_shipped <- function(name) {
  fname <- if (grepl("\\.(json|rds)$", name, ignore.case = TRUE)) name else paste0(name, ".rds")
  list(path = normalizePath(file.path(nle_dir(), fname), mustWork = TRUE), artefact = name)
}

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
    algo <- nle_hash_algo()
    mkey <- if (grepl("\\.rds$", res$path, ignore.case = TRUE)) paste0(res$artefact, ".rds") else basename(res$path)
    pinned <- nle_manifest()[[mkey]]
    if (is.null(pinned))
      stop("Artefact ", mkey, " cannot be checked on ", R.version.string, ": this R has no ",
           "tools::sha256sum (R >= 4.5), and inst/extdata/flownn/MANIFEST carries no md5 line for it")
    if (!identical(pinned, sha))
      stop("Artefact ", mkey, " does not match the ", algo, " pinned in inst/extdata/flownn/MANIFEST")
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
    mlp$layers <- lapply(mlp[["layers", exact = TRUE]], function(l) list(W = as.matrix(l[["W", exact = TRUE]]), b = as.numeric(l[["b", exact = TRUE]])))
    if (isTRUE(mlp[["use_norm", exact = TRUE]]))
      mlp$norms <- lapply(mlp[["norms", exact = TRUE]], function(n)
        list(scale = as.numeric(n[["scale", exact = TRUE]]), bias = as.numeric(n[["bias", exact = TRUE]]), eps = as.numeric(n[["eps", exact = TRUE]])))
    mlp
  }
  members <- if (!is.null(raw[["members", exact = TRUE]])) raw[["members", exact = TRUE]] else list(raw)
  version <- nn_card_version(raw, members[[1]])
  members <- lapply(members, function(m) {
    for (f in c("mlp", "flow_mlp", "classifier_mlp")) if (!is.null(m[[f]])) m[[f]] <- fix_mlp(m[[f]])
    m$scaling <- nn_scaling_spec(m[["scaling", exact = TRUE]])
    nn_natural_context(m)
  })
  kinds <- vapply(members, nn_infer_kind, "")
  if (length(unique(kinds)) != 1L) stop("Ensemble members hold different kinds of network")
  declared <- if (!is.null(raw[["kind", exact = TRUE]])) raw[["kind", exact = TRUE]] else members[[1]][["kind", exact = TRUE]]
  if (!is.null(declared) && !identical(declared, kinds[1]))
    stop("Card declares kind '", declared, "' but holds a ", kinds[1], " network")
  members <- lapply(members, nn_fill_sampled_box)
  m1 <- members[[1]]
  for (m in members[-1])
    for (f in c("context_names", "context_transforms", "bounds_natural", "bounds_sampled", "zerobox_c", "pre",
                "scaling"))
      if (!identical(m[[f]], m1[[f]])) stop("Ensemble members disagree on ", f)
  list(members = members, top = if (!is.null(raw[["members", exact = TRUE]])) raw[names(raw) != "members"] else NULL,
       kind = kinds[1], version = version)
}

# The card's format version: `card_version` at the top level (or on the first
# member), absent = 0 (every card written before the field existed). A version
# this EMC2 does not read is refused by name rather than read as if it were
# one it knows.
nn_card_version <- function(raw, m1) {
  v <- if (!is.null(raw[["card_version", exact = TRUE]])) raw[["card_version", exact = TRUE]] else m1[["card_version", exact = TRUE]]
  if (is.null(v)) return(0)
  if (!is.numeric(v) || length(v) != 1L || !is.finite(v))
    stop("The card's card_version must be one number (", paste(nn_card_versions, collapse = ", "), ")")
  if (!(v %in% nn_card_versions))
    stop("Card version ", format(v), " is not supported by this EMC2, which reads card versions ",
         paste(nn_card_versions, collapse = " and "), "; update EMC2, or export the card for a version it reads")
  v
}

# The scaling rule of a scaled card (kind flow_scaled): a list with `kind`
# ("ddm_unit": the network models the unit-scale decision time of a DDM with
# a = s = 1 given kappa = v a / s^2, eta = sv a / s^2, Z and SZ) and, for
# ddm_unit, `sz` ("proportion": the network's SZ is EMC2's start-point
# proportion). Cards exported before the field was structured carry a
# description beginning "ddm_unit"; it is read as that kind.
nn_scaling_spec <- function(sc) {
  if (is.null(sc)) return(NULL)
  if (is.character(sc)) {
    if (length(sc) == 1L && startsWith(sc, "ddm_unit")) sc <- list(kind = "ddm_unit")
    else stop("The card's scaling is not understood: '", sc[1], "'")
  }
  if (!is.list(sc) || !is.character(sc[["kind", exact = TRUE]]) || length(sc[["kind", exact = TRUE]]) != 1L)
    stop("The card's scaling must give a kind (e.g. list(kind = \"ddm_unit\"))")
  if (sc[["kind", exact = TRUE]] != "ddm_unit") stop("Unsupported scaling kind '", sc[["kind", exact = TRUE]], "' (this EMC2 knows ddm_unit)")
  if (is.null(sc[["pars", exact = TRUE]])) sc$pars <- nn_scaled_pars
  if (!identical(as.character(sc[["pars", exact = TRUE]]), nn_scaled_pars))
    stop("A ddm_unit scaling takes the DDM parameters ", paste(nn_scaled_pars, collapse = ", "),
         " in that order")
  if (is.null(sc[["sz", exact = TRUE]])) sc$sz <- "proportion"
  if (!identical(sc[["sz", exact = TRUE]], "proportion"))
    stop("A ddm_unit scaling's SZ convention must be \"proportion\" (see vignette(\"neural-likelihood-cards\"))")
  sc
}

# A card whose inputs are on the natural scale (context_scale = "natural", no
# context_transforms): every input has transform "identity". The parameters
# keep the twin's transforms on EMC2's sampled scale; the network is fed the
# natural values.
nn_natural_context <- function(m) {
  if (is.null(m[["context_transforms", exact = TRUE]]) && identical(m[["context_scale", exact = TRUE]], "natural") && !is.null(m[["context_names", exact = TRUE]])) {
    ctx <- as.character(m[["context_names", exact = TRUE]])
    m$context_transforms <- as.list(stats::setNames(rep("identity", length(ctx)), ctx))
  }
  m
}

# bounds_sampled from bounds_natural when the card has only the latter (the
# compiled evaluators read the sampled box).
nn_fill_sampled_box <- function(m) {
  bn <- m[["bounds_natural", exact = TRUE]]; tr <- m[["context_transforms", exact = TRUE]]; ctx <- m[["context_names", exact = TRUE]]
  if (!is.null(m[["bounds_sampled", exact = TRUE]]) || is.null(bn) || is.null(tr) || is.null(ctx) ||
      length(bn[["lower", exact = TRUE]]) != length(tr) || !all(ctx %in% names(tr))) return(m)
  j <- match(ctx, names(tr))
  tfs <- vapply(ctx, function(p) as.character(tr[[p]]), "", USE.NAMES = FALSE)
  zc <- nn_zerobox_c(m, ctx, tfs)
  m$bounds_sampled <- list(lower = unname(mapply(nn_tf, bn[["lower", exact = TRUE]][j], tfs, zc)),
                           upper = unname(mapply(nn_tf, bn[["upper", exact = TRUE]][j], tfs, zc)))
  m
}

nn_infer_kind <- function(m) {
  if (!is.null(m[["flow_mlp", exact = TRUE]]) && !is.null(m[["classifier_mlp", exact = TRUE]])) {
    if (identical(m[["kind", exact = TRUE]], "flow_cat") ||
        (!is.null(m[["n_categories", exact = TRUE]]) && m[["n_categories", exact = TRUE]] > 2)) "flow_cat"
    else if (!is.null(m[["scaling", exact = TRUE]])) "flow_scaled"
    else "flow_joint"
  }
  else if (!is.null(m[["mlp", exact = TRUE]]) && !is.null(m[["spline", exact = TRUE]])) "flow_race"
  else if (!is.null(m[["mlp", exact = TRUE]])) if (identical(m[["kind", exact = TRUE]], "mlp_joint")) "mlp_joint" else "regression_joint"
  else stop("Cannot tell what kind of network this card holds (neither flow_mlp + classifier_mlp nor mlp)")
}

# The evaluator for a registration, (re)built on demand. Shipped artefacts are
# re-resolved by name so a saved design survives a moved package library.
nn_ptr <- function(reg) {
  key <- paste0("ptr:", reg$kind, ":", reg$sha256)
  ptr <- nle_cache[[key]]
  valid <- !is.null(ptr) && switch(reg$kind, flow_joint = , flow_scaled = ddm_ptr_valid(ptr),
                                   flow_race = flow_ptr_valid(ptr), flow_cat = cat_ptr_valid(ptr),
                                   mlp_lik_valid(ptr))
  if (!valid) {
    res <- if (!is.null(reg$artefact)) nn_shipped(reg$artefact) else list(path = reg$path)
    hit <- nn_card(res)
    if (!identical(hit$sha256, reg$sha256))
      stop("Neural-likelihood file ", res$path, " has changed since the model was registered ",
           "(sha256 ", reg$sha256, ", now ", hit$sha256, "); re-register it and re-make the design")
    mem <- hit$card$members
    ptr <- switch(reg$kind, flow_joint = , flow_scaled = ddm_build_ensemble(mem),
                  flow_race = flow_build(mem[[1]]), flow_cat = cat_build(mem[[1]]),
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
#' `"log"`, `"probit"` or `"zerobox"`, the scale the network was trained on),
#' `bounds_natural` and/or `bounds_sampled` (the training box, `lower` and
#' `upper`; natural scale in the order of `names(context_transforms)`, sampled
#' scale in `context_names` order; if both are present they must agree), the
#' weights (`flow_mlp` + `classifier_mlp` + `spline` + `scaler` for
#' `"flow_joint"`, `mlp` + `spline` + `scaler` for `"flow_race"`, `mlp` for
#' `"regression_joint"` and `"mlp_joint"`), and optionally `model`, `kind`, `ll_floor_log`,
#' `context_encoding`, `zerobox_c`, and for race flows `context_scale`,
#' `input_scaling`, `affine`, `affine_scale`, `pre` (below). A name in
#' `context_transforms` that is not in
#' `context_names` was fixed when the network was trained; its value must be
#' given in `exceptions` unless the card records it (field `fixed`: named
#' natural-scale values, e.g. `{"s": 1, "SZ": 0}`) or implies it
#' (`context_encoding = "st0zero"` implies `st0 = 0`). `.rds` bundles may hold an ensemble
#' (`members`, `"flow_joint"` only); `.json` cards need the \pkg{jsonlite}
#' package.
#'
#' **Natural-scale race flows.** A `"flow_race"` card may describe its inputs
#' on the natural scale and carry further stages (the cards converted from
#' Lüken's flows by the NLE project do):
#' * `context_scale = "natural"` and no `context_transforms`: every input has
#'   transform `"identity"` and the training box is `bounds_natural`, in
#'   `context_names` order. On EMC2's sampled scale the parameters keep the
#'   twin's transforms (`log` for the [RDM]'s); the network is fed the natural
#'   values.
#' * `input_scaling` (`eps`, `loc`, `scale`, one per input): the network sees
#'   `(log(x + eps) - loc) / scale` in place of a `scaler`; a card with both is
#'   refused, and with neither the network sees the natural values as they are.
#' * `affine = TRUE` with `affine_scale` (`min_scale`, `offset`): the
#'   conditioner has two more outputs after the `3K + 1` spline outputs (`K`
#'   bins), `loc` and `raw_scale`. With `scale = softplus(raw_scale + offset) +
#'   min_scale` the spline is inverted at `y = (log t - loc) / scale`, `z` being
#'   its value there (identity tails in `y`), and `log pdf = dnorm(z, log) +
#'   log|dz/dy| - log(scale) - log t`, `CDF = pnorm(z)`.
#' * `pre` (e.g. `"t0"`): the flow models the decision time `rt - t0`; it is the
#'   default of the argument `pre`, and an argument that contradicts the card
#'   is refused.
#' * An input whose training region starts at the twin's default, where the
#'   twin admits that value as a bound exception (the pulse height `amp` of
#'   [CRDM]: default and lower edge 0), has a region closed at that edge, so
#'   that the default means what it means in the twin; the other edges are
#'   open.
#' * A threshold input `b` registered as EMC2's `B` (`pars = c("v", "s", "B")`)
#'   is exact only without start-point variability, so `A` must be fixed at 0
#'   (`exceptions = c(A = 0)`).
#'
#' **Hybrid race.** With `hybrid = "amp"` an accumulator whose natural-scale
#' `amp` is exactly 0 (an unpulsed accumulator of the [CRDM]) is scored by the
#' closed-form Wald of the [RDM] with its `v`, `B` and `s` (no start-point
#' variability) at the decision time, and every other row by the network. The
#' network never sees the unpulsed rows, so they carry no approximation error
#' and cost no conditioner pass. The training box on `v`, `s` and `B` applies to
#' the Wald rows as well, so that the same parameter vectors are rejected with
#' and without `hybrid`. This holds on every path (`dfun`, `pfun`, the compiled
#' likelihood). Without `hybrid` every row goes to the network, `amp = 0`
#' included. See [CRDM] for the design that gives exact zeros.
#'
#' **Zero-box inputs.** An input with transform `"zerobox"` enters the
#' network as `asinh(x / c)`, `c` > 0 from the card's `zerobox_c` (a named list,
#' one constant per zero-box input; a zero-box input without one, or a constant
#' on any other input, is refused). The network is trained at, and exact at,
#' `x = 0`: on EMC2's sampled scale the parameter keeps the twin's transform
#' (for [DDM]: `log` for `sv` and `st0`, `probit` for `SZ`; `"exp"` when you
#' give `p_types`), so `-Inf` is an exact zero, and the network value is built
#' from the natural value (`asinh(0 / c) = 0`). The training region is closed at
#' 0 for these inputs; the bounds are the twin's (e.g. [DDM]'s floor of 0.01
#' for `sv` and `SZ`) with exactly 0 admitted as a bound exception, so a
#' default or constant of 0 is accepted and values between 0 and the floor are
#' rejected as the analytic model rejects them. A zero-box DDM card
#' (`context_encoding = "zerobox"`) needs the twin's defaults `sv = SZ = st0 =
#' -Inf` (exactly zero, as in [DDM]). **Support rule:** for such a card, when
#' `st0` is exactly 0 the density and the defective CDF are 0 (log density
#' `-Inf`) for `rt <= t0`, whatever the flow would give; for `st0 > 0` the flow
#' models the smear itself. The rule is applied by the evaluator on every path
#' (one parameter vector, trial-wise parameters, ensembles, the compiled
#' likelihood); `min_ll` flooring applies after it.
#'
#' **Card versions.** A card states its format version in `card_version`
#' (absent: version 0, every card written before the field existed); this
#' EMC2 reads versions 0 and 1 and refuses any other by name. The format,
#' field by field, is `vignette("neural-likelihood-cards", package =
#' "EMC2")`; a version-1 card may also carry `support`, `mass_tolerance`
#' and `beyond_support` (read by [nn_total_mass()]), `control_limits` (race
#' cards: limits of the analytic control, warned about by [nn_cell()]) and
#' `validation` (what the card passed; a version-1 card without one
#' registers with a message). The golden vectors beside a card are checked
#' by [nn_golden_check()].
#'
#' **Kinds.**
#' * `"flow_joint"`: a spline flow for rt given parameters and response plus a
#'   classifier for the response; a two-response joint model like [DDM].
#'   `dfun`/`pfun` are the joint density and the defective CDF.
#' * `"flow_scaled"`: the scaled DDM card (DDMnle): a joint flow over the
#'   unit-scale decision time `u = (rt - t0) (s / a)^2` of a diffusion with
#'   boundary 1 and noise 1, with inputs `kappa = v a / s^2`, `eta = sv a /
#'   s^2`, `Z` and `SZ` (card field `scaling`, kind `ddm_unit`). The
#'   rescaling by `a` and `s`, the shift by `t0` and the uniform `st0` window
#'   (a difference of two of the network's CDFs) are applied exactly, on the R
#'   path and in the compiled likelihood. The model is [DDM]'s in every other
#'   respect (parameters, transforms, bounds, the exact zeros `sv = SZ = st0 =
#'   0`, the simulator); `pars`, `exceptions`, `pre` and `hybrid` do not
#'   apply, and `cdf` is `FALSE` (no CDF of rt is offered). `SZ` is the
#'   start-point proportion throughout (see [nn_sz_convert()] for the
#'   conventions).
#' * `"flow_cat"`: a K-category joint flow (card fields `n_categories`, even,
#'   and `conv`): a K-way classifier plus a spline flow for the decision time
#'   given the category. The category is `k = rating_fold(R, RR)` over the two
#'   response levels and the rating column `RR` (1..K/2, an extra response of
#'   the model: the first level counts its ratings down, the second up). The
#'   non-decision time is applied outside the network by an exact uniform
#'   convolution whose mean and range are the parameters `conv$mean` and
#'   `conv$range` (e.g. `t0` and `st0`, not network inputs). Lower and upper
#'   censored trials (`missingness` codes 1 and 2, from [make_missing()]) are
#'   scored by the network's CDF. A card with a `mix` block integrates a
#'   latent trial strength by Gauss-Hermite quadrature outside the network;
#'   its `context_names` are then the model's natural parameters, from which
#'   the evaluator derives the network's own inputs (`mix$net_names`). See
#'   [RTCONnn()].
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
#' by `pre`, nor fixed by `exceptions`); a `pre` that is not a shift of rt or
#' contradicts the card's; a card with both `scaler` and `input_scaling`; a
#' `hybrid` that is not a network input.
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
#'   `"probit"`, `"zerobox"`), in the network's order or named. Default: the card's
#'   `context_transforms`; if both are given they must agree.
#' @param kind `"flow_joint"`, `"flow_scaled"`, `"flow_cat"`, `"flow_race"`,
#'   `"regression_joint"` or `"mlp_joint"`. Default:
#'   inferred from the card; if given it must match.
#' @param cdf Whether the network provides a CDF. Default `TRUE` for flows,
#'   `FALSE` for the MLP kinds. With `cdf = FALSE`, `pfun` refuses.
#' @param exceptions Named natural-scale values of twin parameters that the
#'   network does not take as input because they were fixed in training
#'   (e.g. `c(st0 = 0)`). The model admits exactly that value and refuses
#'   any other. Values the card records in its `fixed` field are the
#'   default; an exception that contradicts them is refused.
#' @param pre `NULL`, the name of a parameter subtracted from rt before the
#'   network sees it (`"t0"` for a decision-time network), or a
#'   `function(rt, pars)` returning the network's time; it must be a shift
#'   (no Jacobian is applied). Times <= 0 get zero density and CDF. Default:
#'   the card's `pre`, if it has one; an argument that differs from the card's
#'   is refused.
#' @param twin The analytic twin: an EMC2 model function supplying
#'   `p_types`, `transform` and `rfun` (and `Ttransform`, applied before
#'   `rfun`). Default: from the card's `model` field. For a network with no
#'   analytic model, a hand-written model function (see Details), or use
#'   `p_types` instead.
#' @param p_types For a network with no analytic model: named defaults on the
#'   sampled scale for the model's parameters (see Details). Give either
#'   `twin` or `p_types`; with `p_types` no twin is used, even if the card
#'   names a model.
#' @param hybrid For a race flow: the name of a network input (an EMC2
#'   parameter name, e.g. `"amp"`) such that a row where it is exactly 0 is
#'   scored by the closed-form Wald instead of the network (see Details).
#'   Default `NULL`: every row goes to the network.
#'
#' @return A model function (like [DDM]) whose list carries, besides the usual
#'   elements, `nn`: the registration (`kind`, `sha256`, `pars` and
#'   `context_names` in the network's order, `transforms`, the training box,
#'   `fixed`, `pre`, `hybrid`, `cdf`, `ll_floor`, the card's metadata and `twin`: the
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
                              p_types = NULL, hybrid = NULL) {
  res <- nn_resolve(path)
  hit <- nn_card(res)
  card <- hit$card
  m1 <- card[["members", exact = TRUE]][[1]]
  label <- if (!is.null(res$artefact)) res$artefact else nn_label(res$path, m1)
  what <- paste0("Neural likelihood '", label, "': ")
  if (!is.null(kind)) {
    kind <- match.arg(kind, nn_kinds)
    if (kind != card[["kind", exact = TRUE]]) stop(what, "the card holds a ", card[["kind", exact = TRUE]], " network, not ", kind)
  }
  kind <- card[["kind", exact = TRUE]]
  if (kind != "flow_joint" && length(card[["members", exact = TRUE]]) > 1L)
    stop(what, "ensembles are only supported for flow_joint networks")

  # --- inputs: order, names, scales -----------------------------------------
  ctx <- m1[["context_names", exact = TRUE]]
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
  # A scaled card's inputs are computed from the DDM's parameters (scaling,
  # see nn_scaling_spec); they are not parameters and cannot be renamed, fixed
  # or shifted.
  scaled <- kind == "flow_scaled"
  if (scaled) {
    if (!identical(ctx, nn_scaled_inputs_names))
      stop(what, "a scaled DDM card's inputs must be ", paste(nn_scaled_inputs_names, collapse = ", "),
           " in that order; the card has ", paste(ctx, collapse = ", "))
    if (!identical(pars, ctx)) stop(what, "a scaled DDM card's inputs cannot be renamed (pars)")
    if (!is.null(exceptions)) stop(what, "a scaled DDM card takes every DDM parameter; exceptions do not apply")
    if (!is.null(pre)) stop(what, "a scaled DDM card applies t0 itself; pre does not apply")
    if (!is.null(hybrid)) stop(what, "hybrid is for race flows")
    if (!identical(m1[["context_encoding", exact = TRUE]], "zerobox"))
      stop(what, "a scaled DDM card must have context_encoding \"zerobox\" (eta and SZ exact at 0)")
  }
  tr_card <- m1[["context_transforms", exact = TRUE]]
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

  if (kind == "flow_race") nn_check_flow_race(m1, k, what)

  # --- training box -----------------------------------------------------------
  zc <- nn_zerobox_c(m1, ctx, tfs, what)
  box <- nn_box(m1, ctx, tr_card, tfs, zc, what)
  if (!(kind %in% nn_mlp_kinds) && is.null(m1[["bounds_sampled", exact = TRUE]]))   # read by the compiled evaluators
    stop(what, "a flow card must carry bounds_sampled (or bounds_natural with context_transforms)")

  # --- parameterisation: the analytic twin, or the user's own defaults -----------
  want_type <- if (kind == "flow_race") "RACE" else "DDM"
  if (!is.null(p_types)) {
    if (!is.null(twin))
      stop(what, "give either twin (a model supplies the defaults) or p_types (your own defaults), not both")
    tl <- nn_own_parameterisation(p_types, pars, tfs, want_type, what)
  } else {
    if (is.null(twin) && !is.null(m1[["model", exact = TRUE]]) && m1[["model", exact = TRUE]] %in% names(nn_twin_names))
      twin <- get(nn_twin_names[[m1[["model", exact = TRUE]]]], envir = asNamespace("EMC2"))
    if (is.null(twin))
      stop(what, "no analytic twin (card model '", if (is.null(m1[["model", exact = TRUE]])) "" else m1[["model", exact = TRUE]], "'). ",
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
  notp <- setdiff(if (scaled) nn_scaled_pars else pars, names(pt))
  if (length(notp)) stop(what, if (scaled) "the DDM parameter(s) " else "network input(s) ",
                         paste(notp, collapse = ", "),
                         " are not parameters of the model (", paste(names(pt), collapse = ", "), ")")
  tr_filled <- fill_transform(NULL, function() list(p_types = pt, transform = tl$transform))
  defaults <- nn_natural(pt, tr_filled)
  zb <- pars[tfs == "zerobox"]      # closed at 0; see `closed` below for the other closed edges
  # A zero-box DDM network is exact at sv = SZ = st0 = 0, so an omitted
  # parameter must mean exactly zero, as in the analytic DDM() (a positive
  # default would simulate one model and evaluate another). Checked whatever
  # supplies the defaults (DDM, a hand-written twin, or p_types). A scaled
  # card is exact at eta = 0 (sv = 0), SZ = 0 and, by construction, st0 = 0.
  if (identical(m1[["context_encoding", exact = TRUE]], "zerobox") && kind %in% c("flow_joint", "flow_scaled")) {
    zp <- intersect(c("sv", "SZ", "st0"), names(pt))
    if (!all(is.infinite(pt[zp]) & pt[zp] < 0))
      stop(what, "a zero-box DDM network needs the defaults sv = SZ = st0 = -Inf (exactly zero)")
  }

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
  # the card may record the values (field `fixed`, natural scale): they are the
  # default exceptions, and an exception may not contradict them
  card_fixed <- m1[["fixed", exact = TRUE]]
  if (is.null(card_fixed)) card_fixed <- card[["top", exact = TRUE]][["fixed", exact = TRUE]]
  if (!is.null(card_fixed)) {
    if (!is.list(card_fixed) && !is.numeric(card_fixed)) stop(what, "the card's fixed must be named numbers")
    if (!all(vapply(card_fixed, function(x) is.numeric(x) && length(x) == 1L && is.finite(x), TRUE)))
      stop(what, "the card's fixed must be named finite numbers")
    card_fixed <- unlist(card_fixed)
    if (is.null(names(card_fixed)) || any(names(card_fixed) == "")) stop(what, "the card's fixed must be named numbers")
    bad <- setdiff(names(card_fixed), dropped)
    if (length(bad)) stop(what, "the card's fixed names '", bad[1], "', which is not a parameter the network was ",
                          "trained without (in context_transforms but not context_names)")
    for (p in intersect(names(exceptions), names(card_fixed)))
      if (!isTRUE(all.equal(exceptions[[p]], card_fixed[[p]], tolerance = 0)))
        stop(what, "exceptions = c(", p, " = ", format(exceptions[[p]]), ") contradicts the card's fixed (", p,
             " = ", format(card_fixed[[p]]), "): the network was trained at that value only")
    exceptions <- c(exceptions, card_fixed[setdiff(names(card_fixed), names(exceptions))])
  }
  fixed <- exceptions
  for (p in setdiff(dropped, names(fixed))) {
    if (identical(m1[["context_encoding", exact = TRUE]], "st0zero") && p == "st0") fixed <- c(fixed, st0 = 0)
    else stop(what, "the network was trained without '", p, "' (it is in context_transforms but not ",
              "context_names); declare the value it was fixed at with exceptions = c(", p, " = <value>)")
  }
  if (is.null(fixed)) fixed <- stats::setNames(numeric(0), character(0))
  # the network's threshold b is B + A: registered as B it needs A = 0 exactly
  if (any(ctx == "b" & pars == "B") && "A" %in% setdiff(names(pt), names(fixed)[fixed == 0]))
    stop(what, "the network's input 'b' (the threshold, B + A) is registered as 'B', which is exact ",
         "only without start-point variability: fix it with exceptions = c(A = 0)")
  # flow_cat: the convolution's mean and range (card `conv`) are parameters
  # the network does not take; they are consumed like a `pre` shift
  conv <- NULL
  if (kind == "flow_cat") {
    if (!is.null(pre) || !is.null(m1[["pre", exact = TRUE]]))
      stop(what, "a K-category network applies t0 by its uniform convolution; pre is not allowed")
    cv <- m1[["conv", exact = TRUE]]
    if (is.null(cv) || is.null(cv[["mean", exact = TRUE]]) || is.null(cv[["range", exact = TRUE]]))
      stop(what, "a flow_cat card must name the convolution parameters: conv = list(mean = , range = )")
    conv <- c(mean = as.character(cv[["mean", exact = TRUE]]), range = as.character(cv[["range", exact = TRUE]]))
    bad <- setdiff(conv, names(pt))
    if (length(bad)) stop(what, "convolution parameter(s) ", paste(bad, collapse = ", "), " are not parameters of the model")
    if (any(conv %in% pars)) stop(what, "the convolution parameters cannot also be network inputs")
    nk <- m1[["n_categories", exact = TRUE]]
    if (!is.numeric(nk) || length(nk) != 1L || nk < 2) stop(what, "a flow_cat card must give n_categories (>= 2)")
    # the mix block: a latent trial strength integrated by quadrature outside
    # the network; the network's own inputs (mix.net_names) are derived from
    # the registration-facing inputs (context_names), which stay the model's
    # natural parameters per design row
    mx <- m1[["mix", exact = TRUE]]
    if (!is.null(mx)) {
      need <- c(mx[["mean", exact = TRUE]], mx[["sd", exact = TRUE]], unlist(mx[["map_inputs", exact = TRUE]]),
                unlist(mx[["direct", exact = TRUE]]))
      bad <- setdiff(need, ctx)
      if (length(bad)) stop(what, "mix names ", paste(bad, collapse = ", "), " which are not in context_names")
      nodes <- mx[["nodes", exact = TRUE]]
      if (!is.numeric(nodes) || length(nodes) != 1L || nodes < 1 ||
          length(mx[["gh_x", exact = TRUE]]) != nodes || length(mx[["gh_w", exact = TRUE]]) != nodes)
        stop(what, "mix needs nodes, gh_x and gh_w of that length")
      if (length(mx[["net_names", exact = TRUE]]) != length(m1[["scaler", exact = TRUE]][["mean", exact = TRUE]]))
        stop(what, "mix.net_names must have one entry per network input (the scaler's length)")
    }
  }
  if (!is.null(m1[["pre", exact = TRUE]])) {
    card_pre <- m1[["pre", exact = TRUE]]
    if (!is.character(card_pre) || length(card_pre) != 1L)
      stop(what, "the card's pre must name one parameter (e.g. \"t0\")")
    if (is.null(pre)) pre <- card_pre
    else if (!identical(pre, card_pre))
      stop(what, "pre = ", if (is.character(pre)) paste0("\"", pre[1], "\"") else "<function>",
           " contradicts the card, whose network models rt - ", card_pre, " (pre = \"", card_pre, "\")")
  }
  pre_pars <- nn_check_pre(pre, pars, defaults, box, what)
  ignored <- setdiff(names(pt), c(if (scaled) nn_scaled_pars else pars, names(fixed), pre_pars, conv))
  if (length(ignored))
    stop(what, "the model's parameter(s) ", paste(ignored, collapse = ", "), " are not network inputs, ",
         "not used by pre and not fixed by exceptions: the network would ignore them")

  # --- CDF --------------------------------------------------------------------------
  is_mlp <- kind %in% nn_mlp_kinds
  if (is.null(cdf)) cdf <- !is_mlp && !scaled
  if (!is.logical(cdf) || length(cdf) != 1L || is.na(cdf)) stop(what, "cdf must be TRUE or FALSE")
  if (is_mlp && cdf) stop(what, "a ", kind, " network has no CDF (cdf must be FALSE)")
  if (scaled && cdf)
    stop(what, "a scaled DDM card offers no CDF of rt (its CDF is of the unit-scale time; cdf must be FALSE)")
  if (kind == "flow_race" && !cdf)
    stop(what, "a race needs the survivor function (1 - CDF) of the losing accumulators, ",
         "so a flow_race network without a CDF cannot be assembled into a race")
  if (is_mlp) nn_check_regression(m1, ctx, kind, what)
  oob <- if (kind == "mlp_joint") m1[["ll_floor_log", exact = TRUE]] else -Inf

  # --- edges of the training region that are closed; the hybrid race ---------------
  # zerobox inputs at 0; on a natural-scale card, an input whose region starts
  # at the twin's default where the twin admits that value (a bound exception)
  twin_ex <- tl$bound$exception
  at_default <- vapply(seq_along(pars), function(j) {
    p <- pars[j]
    identical(m1[["context_scale", exact = TRUE]], "natural") && tfs[j] == "identity" && p %in% names(twin_ex) &&
      isTRUE(defaults[[p]] == box$lower[j]) && twin_ex[[p]] == box$lower[j]
  }, TRUE)
  closed <- union(zb, pars[at_default])
  wald <- NULL
  if (!is.null(hybrid)) {
    if (!is.character(hybrid) || length(hybrid) != 1L || is.na(hybrid))
      stop(what, "hybrid must name one network input (e.g. \"amp\")")
    if (kind != "flow_race") stop(what, "hybrid is for race flows; this is a ", kind, " network")
    if (!(hybrid %in% pars))
      stop(what, "hybrid = \"", hybrid, "\" is not an input of the network (", paste(pars, collapse = ", "), ")")
    if (!(hybrid %in% closed) || box$lower[match(hybrid, pars)] != 0)
      stop(what, "hybrid = \"", hybrid, "\" needs a training region that starts at 0, the parameter's ",
           "default in the twin (where the accumulator is the twin's Wald)")
    wp <- c(v = "v", b = "B", s = "s")
    if (!all(wp %in% pars) || hybrid %in% wp)
      stop(what, "a hybrid race scores the rows with ", hybrid, " = 0 by the Wald of v, B and s, ",
           "which must be network inputs")
    if ("A" %in% setdiff(names(pt), names(fixed)[fixed == 0]))
      stop(what, "a hybrid race needs exceptions = c(A = 0): its Wald has no start-point variability")
    j <- match(wp, pars)
    wald <- list(pars = unname(wp), lower = box$lower[j], upper = box$upper[j])
  }

  # --- bounds -------------------------------------------------------------------------
  twin_mm <- tl$bound$minmax
  minmax <- vapply(names(pt), function(p) {
    if (p %in% closed) {
      # closed edge (zerobox, ...): the analytic model's floor, with the edge
      # itself admitted through a bound exception (below), capped at the
      # training region's top
      j <- match(p, pars)
      if (!is.null(twin_mm) && p %in% colnames(twin_mm))
        c(max(twin_mm[1, p], box$lower[j]), min(twin_mm[2, p], box$upper[j]))
      else c(box$lower[j], box$upper[j])
    }
    else if (p %in% pars) { j <- match(p, pars); c(box$lower[j], box$upper[j]) }
    else if (p %in% names(fixed)) rep(fixed[[p]], 2)
    else if (!is.null(tl$bound$minmax) && p %in% colnames(tl$bound$minmax)) tl$bound$minmax[, p]
    else c(-Inf, Inf)
  }, numeric(2))
  lower <- stats::setNames(box$lower, pars); upper <- stats::setNames(box$upper, pars)
  # inputs that are model parameters (all of them, except for a scaled card,
  # whose kappa and eta are computed): a default outside the box is refused
  in_pars <- intersect(pars, names(defaults))
  refuse_default <- in_pars[!nn_inside(in_pars, defaults[in_pars], lower, upper, closed)]
  exception <- reg_exception(reg_fixed = fixed, zb = intersect(closed, names(pt)), tl = tl)
  # the convolution is exact outside the network, so the twin's exceptions on
  # its parameters (e.g. st0 = 0) hold as they do for the analytic model
  if (kind == "flow_cat") {
    te <- tl$bound$exception
    keep <- intersect(names(te), setdiff(conv, names(exception)))
    exception[keep] <- te[keep]
  }
  if (scaled) {
    # the DDM's own bounds and exceptions (sv = SZ = st0 = 0): the network's
    # box on kappa, eta, Z and SZ is applied by the evaluator; st0 is free (the
    # uniform window is exact outside the network, by a CDF difference), its
    # exception admitting st0 = 0 as for the analytic DDM
    if (is.null(twin_mm)) stop(what, "the twin of a scaled DDM card must carry bound$minmax")
    minmax <- twin_mm
    exception <- if (is.null(tl$bound$exception)) stats::setNames(numeric(0), character(0)) else tl$bound$exception
  }

  # (a card's embedded golden vectors stay with the card: every saved design
  # carries the registration)
  card_meta <- m1[setdiff(names(m1), c("mlp", "flow_mlp", "classifier_mlp", "scaler",
                                        "classifier_scaler", "spline", "golden"))]
  card_meta <- c(card_meta, card[["top", exact = TRUE]][setdiff(names(card[["top", exact = TRUE]]), names(card_meta))])
  reg <- list(
    kind = kind, label = label, artefact = res$artefact, path = res$path, sha256 = hit$sha256,
    pars = pars, context_names = ctx, transforms = tfs, transform_codes = unname(nn_transform_codes[tfs]),
    zerobox_c = zc, zerobox = zb, closed = closed, exception = exception,
    hybrid = hybrid, wald = wald,
    lower = lower, upper = upper, lower_s = box$lower_s, upper_s = box$upper_s,
    fixed = fixed, pre = pre, pre_pars = pre_pars, cdf = cdf,
    conv = conv, n_categories = if (kind == "flow_cat") as.integer(m1[["n_categories", exact = TRUE]]),
    defaults = defaults, refuse_default = refuse_default, tr_filled = tr_filled,
    ll_floor = m1[["ll_floor_log", exact = TRUE]], oob = oob, n_members = length(card[["members", exact = TRUE]]),
    input_layout = m1[["input_layout", exact = TRUE]], card = card_meta, card_version = card[["version", exact = TRUE]],
    # a scaled card: the rule and the DDM parameters its inputs are computed from
    scaled = scaled, scaling = if (scaled) m1[["scaling", exact = TRUE]], sc_pars = if (scaled) nn_scaled_pars,
    # race cards may state the limits of the analytic control (vignette)
    control_limits = nn_control_limits(m1, kind, what),
    # the analytic twin, the default control of the validation kit
    # (R/nn_validate.R): the name of an EMC2 model, a hand-written model
    # function, or NULL (p_types)
    twin = nn_twin_ref(twin))
  reg$native <- nn_native_spec(reg)
  class(reg) <- "emc_nn"
  if (card[["version", exact = TRUE]] >= 1 && is.null(m1[["validation", exact = TRUE]]))
    message(what, "the card records no validation (score bias, posterior shift, SBC); loading a network ",
            "is not the same as being able to infer with it (see ?nn_cell)")
  nn_model_function(nn_model_list(reg, tl, minmax))
}

# The limits a race card may state for its analytic control (the RDM's Wald):
# `control_limits = list(s = c(lower, upper), vb_over_s2_max = 354)`: the s
# range the control is exact on and the largest v * b / s^2 its CDF (pWald)
# evaluates without overflow. Recorded on the registration and used by the
# validation kit (nn_cell) to warn; a box that reaches past the limit is noted
# at registration, never refused (it is a property of the control).
nn_control_limits <- function(m, kind, what) {
  cl <- m[["control_limits", exact = TRUE]]
  if (is.null(cl)) return(NULL)
  if (kind != "flow_race") stop(what, "control_limits are for race cards")
  if (!is.list(cl)) stop(what, "control_limits must be a list (s = c(lower, upper), vb_over_s2_max = )")
  if (!is.null(cl[["s", exact = TRUE]]) && !(is.numeric(cl[["s", exact = TRUE]]) && length(cl[["s", exact = TRUE]]) == 2L && cl[["s", exact = TRUE]][1] < cl[["s", exact = TRUE]][2]))
    stop(what, "control_limits$s must be c(lower, upper)")
  if (!is.null(cl[["vb_over_s2_max", exact = TRUE]]) && !(is.numeric(cl[["vb_over_s2_max", exact = TRUE]]) && length(cl[["vb_over_s2_max", exact = TRUE]]) == 1L &&
                                       cl[["vb_over_s2_max", exact = TRUE]] > 0))
    stop(what, "control_limits$vb_over_s2_max must be one positive number")
  cl
}

# The largest v * b / s^2 a natural-scale box (named lower/upper over v, B, A, s)
# can reach; b = B + A. NA when a name is missing.
nn_vb_over_s2_max <- function(lower, upper) {
  if (!all(c("v", "B", "s") %in% names(lower))) return(NA_real_)
  b <- upper[["B"]] + if ("A" %in% names(upper)) upper[["A"]] else 0
  upper[["v"]] * b / lower[["s"]]^2
}

# ---------------------------------------------------------------------------
# The SZ conventions (a contract; see vignette("neural-likelihood-cards"))
# ---------------------------------------------------------------------------

#' Convert Between the Start-Point Variability Conventions
#'
#' Three conventions for the DDM's start-point variability are in use, and
#' every interface states which it takes:
#' * `"proportion"`: EMC2's natural-scale `SZ`, a proportion in `[0, 1)` of
#'   the width available at the start point `Z`; the value a model's
#'   `p_types`/`transform` produce from the sampled probit, what the
#'   simulator takes, and the input of the shipped DDM flows and the scaled
#'   DDM card (`scaling$sz = "proportion"`).
#' * `"width"`: the absolute relative width `sw = 2 * SZ * min(Z, 1 - Z)`
#'   that [DDM()]'s `Ttransform` passes to its likelihood (WienR's `sw`).
#'
#' The two are exact inverses of each other at a given `Z`.
#'
#' @param x Values in the `from` convention.
#' @param Z The relative start point, natural scale in (0, 1), recycled.
#' @param from,to The conventions.
#' @return `x` in the `to` convention.
#' @examples
#' nn_sz_convert(0.2, Z = 0.3, from = "proportion", to = "width")   # 2 * 0.2 * 0.3
#' nn_sz_convert(0.12, Z = 0.3, from = "width", to = "proportion")
#' @export
nn_sz_convert <- function(x, Z, from = c("proportion", "width"), to = c("width", "proportion")) {
  from <- match.arg(from); to <- match.arg(to)
  if (from == to) return(x)
  m <- 2 * pmin(Z, 1 - Z)
  if (from == "proportion") x * m else ifelse(x > 0, x / m, 0)
}

# ---------------------------------------------------------------------------
# Evaluation (R path): the scaled DDM card
# ---------------------------------------------------------------------------

# The network's four natural-scale inputs from the DDM's parameters (SZ as the
# proportion, the registered model's convention).
nn_scaled_inputs <- function(pars) {
  s2 <- pars[, "s"]^2; a <- pars[, "a"]
  cbind(kappa = pars[, "v"] * a / s2, eta = pars[, "sv"] * a / s2, Z = pars[, "Z"],
        SZ = ifelse(pars[, "SZ"] > 0, pars[, "SZ"], 0))
}

# The likelihood of rt given the DDM's parameters through a scaled card: the
# R reference of the compiled path (src/model_NN.h, nn_trial_ll_scaled), as
# the NLE project's port/R/ddm_scaled.R ddm_sc_density(), which takes SZ as the
# width (this function takes the proportion). `evaluator(U, u, R)` returns
# list(pdf, cdf) of the unit-scale time u for the natural inputs U; the
# default is the card's network; an analytic oracle serves the tests.
#   f(rt) = c g(u1),                       c = (s/a)^2, u1 = (rt - t0) c   (st0 = 0)
#   f(rt) = [G(u1) - G(max(u0, 0))] / st0, u0 = (rt - t0 - st0) c          (st0 > 0)
#   f(rt) = 0                                                            (rt <= t0)
# A window too narrow to difference (st0 < 1e-7 d) uses the density at its
# midpoint.
nn_eval_scaled <- function(reg, rt, R, pars, what = "pdf", evaluator = NULL) {
  if (what != "pdf")
    stop("Neural likelihood '", reg$label, "' (a scaled DDM card) offers no CDF of rt")
  miss <- setdiff(nn_scaled_pars, colnames(pars))
  if (length(miss)) stop("Neural likelihood '", reg$label, "' needs parameter(s) ",
                         paste(miss, collapse = ", "), " which the model does not supply")
  if (is.null(evaluator)) evaluator <- function(U, u, R) {
    o <- ddm_ens_eval_trials_cpp(nn_ptr(reg), nn_context(U, reg), u, as.integer(R))
    list(pdf = o[["pdf"]], cdf = o[["cdf"]])
  }
  n <- length(rt); out <- numeric(n); R <- as.integer(R)
  U <- nn_scaled_inputs(pars)
  c2 <- (pars[, "s"] / pars[, "a"])^2
  d <- rt - pars[, "t0"]; st0 <- pars[, "st0"]
  live <- !is.na(d) & d > 0
  if (!any(live)) return(out)
  narrow <- live & (st0 <= 0 | st0 < 1e-7 * d)
  wide <- live & !narrow
  if (any(narrow)) { i <- which(narrow)
    um <- (d[i] - pmax(st0[i], 0) / 2) * c2[i]
    out[i] <- c2[i] * evaluator(U[i, , drop = FALSE], um, R[i])$pdf }
  if (any(wide)) { i <- which(wide)
    u1 <- d[i] * c2[i]; u0 <- (d[i] - st0[i]) * c2[i]
    pos <- u0 > 0
    j <- c(seq_along(i), which(pos))
    G <- evaluator(U[i[j], , drop = FALSE], c(u1, u0[pos]), R[i[j]])$cdf
    G0 <- numeric(length(i)); G0[pos] <- G[length(i) + seq_len(sum(pos))]
    out[i] <- pmax(G[seq_along(i)] - G0, 0) / st0[i] }
  out
}

# The label of a card file: its file name, or for a file called card.* (one
# folder per network) the card's `model` field, else the folder's name.
nn_label <- function(path, m) {
  if (!identical(tools::file_path_sans_ext(basename(path)), "card")) return(basename(path))
  if (is.character(m[["model", exact = TRUE]]) && length(m[["model", exact = TRUE]]) == 1L && nzchar(m[["model", exact = TRUE]])) m[["model", exact = TRUE]] else basename(dirname(path))
}

# The optional stages of a race flow's card (src/flow_race.cpp evaluates them).
nn_check_flow_race <- function(m, k, what) {
  if (!is.null(m[["scaler", exact = TRUE]]) && !is.null(m[["input_scaling", exact = TRUE]]))
    stop(what, "the card has both a scaler and input_scaling; a network has one input standardisation")
  if (!is.null(m[["input_scaling", exact = TRUE]])) {
    isc <- m[["input_scaling", exact = TRUE]]
    if (!all(c("eps", "loc", "scale") %in% names(isc)) ||
        any(lengths(isc[c("eps", "loc", "scale")]) != k) || !is.numeric(unlist(isc[c("eps", "loc", "scale")])))
      stop(what, "input_scaling needs eps, loc and scale, one number per network input (", k, ")")
    if (any(unlist(isc[["scale", exact = TRUE]]) == 0)) stop(what, "input_scaling has a zero scale")
  }
  if (!is.null(m[["scaler", exact = TRUE]]) && (length(m[["scaler", exact = TRUE]][["mean", exact = TRUE]]) != k || length(m[["scaler", exact = TRUE]][["scale", exact = TRUE]]) != k))
    stop(what, "the scaler needs one mean and scale per network input (", k, ")")
  if (nrow(m[["mlp", exact = TRUE]][["layers", exact = TRUE]][[1]][["W", exact = TRUE]]) != k)
    stop(what, "the network has ", nrow(m[["mlp", exact = TRUE]][["layers", exact = TRUE]][[1]][["W", exact = TRUE]]), " inputs; the card names ", k)
  affine <- isTRUE(m[["affine", exact = TRUE]])
  if (!is.null(m[["affine", exact = TRUE]]) && !is.logical(m[["affine", exact = TRUE]])) stop(what, "the card's affine must be true or false")
  if (affine && !(is.numeric(m[["affine_scale", exact = TRUE]][["min_scale", exact = TRUE]]) && is.numeric(m[["affine_scale", exact = TRUE]][["offset", exact = TRUE]])))
    stop(what, "an affine card needs affine_scale with min_scale and offset")
  n_out <- ncol(m[["mlp", exact = TRUE]][["layers", exact = TRUE]][[length(m[["mlp", exact = TRUE]][["layers", exact = TRUE]])]][["W", exact = TRUE]])
  want <- 3L * m[["spline", exact = TRUE]][["num_bins", exact = TRUE]] + if (affine) 3L else 1L
  if (n_out != want)
    stop(what, "the conditioner has ", n_out, " outputs; a", if (affine) "n affine" else "",
         " spline of ", m[["spline", exact = TRUE]][["num_bins", exact = TRUE]], " bins needs ", want)
  if (!is.null(m[["spline", exact = TRUE]][["boundary_slopes", exact = TRUE]]) && !identical(m[["spline", exact = TRUE]][["boundary_slopes", exact = TRUE]], "identity"))
    stop(what, "only boundary_slopes = \"identity\" is supported")
  invisible(TRUE)
}

# An EMC2 analytic model is kept by name (small, and resolved in the installed
# package); a hand-written twin as the function itself.
# Inside the training region: open at both ends, except that a closed lower
# end (`zb`: a zerobox input's 0, where the network is exact; an edge at the
# twin's default) is inside.
nn_inside <- function(p, x, lower, upper, zb)
  ifelse(p %in% zb, x >= lower[p], x > lower[p]) & x < upper[p]

# Bound exceptions: fixed parameters admit exactly their value; an input with a
# closed lower edge (`zb`) admits exactly that edge (0) below the analytic
# model's floor (the twin's exception).
reg_exception <- function(reg_fixed, zb, tl) {
  ex <- reg_fixed
  for (p in zb) {
    te <- tl$bound$exception
    ex[[p]] <- if (!is.null(te) && p %in% names(te)) te[[p]] else 0
  }
  if (is.null(ex)) stats::setNames(numeric(0), character(0)) else ex
}

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
  func[pars] <- unname(c(identity = "identity", log = "exp", probit = "pnorm", zerobox = "exp")[tfs])
  list(type = type, p_types = p_types, transform = list(func = func))
}

# Training box: natural scale per input (network order), sampled scale per
# input. Cross-checks the card's two boxes when both are present.
nn_box <- function(m, ctx, tr_card, tfs, zc, what) {
  bn <- m[["bounds_natural", exact = TRUE]]; bs <- m[["bounds_sampled", exact = TRUE]]
  if (is.null(bn) && is.null(bs)) stop(what, "the card has no training box (bounds_natural / bounds_sampled)")
  k <- length(ctx)
  if (!is.null(bn)) {
    if (length(bn[["lower", exact = TRUE]]) != length(tr_card) || length(bn[["upper", exact = TRUE]]) != length(tr_card))
      stop(what, "bounds_natural needs one entry per context_transforms entry (", length(tr_card), ")")
    j <- match(ctx, names(tr_card))
    lo <- bn[["lower", exact = TRUE]][j]; hi <- bn[["upper", exact = TRUE]][j]
    lo_s <- mapply(nn_tf, lo, tfs, zc); hi_s <- mapply(nn_tf, hi, tfs, zc)
  }
  if (!is.null(bs)) {
    if (length(bs[["lower", exact = TRUE]]) != k || length(bs[["upper", exact = TRUE]]) != k)
      stop(what, "bounds_sampled needs one entry per network input (", k, ")")
    if (is.null(bn)) {
      lo_s <- bs[["lower", exact = TRUE]]; hi_s <- bs[["upper", exact = TRUE]]
      lo <- mapply(nn_tf_inv, lo_s, tfs, zc); hi <- mapply(nn_tf_inv, hi_s, tfs, zc)
    } else {
      agree <- function(a, b) (is.infinite(a) & a == b) | abs(a - b) <= 1e-8 * pmax(1, abs(b))
      bad <- which(!(agree(lo_s, bs[["lower", exact = TRUE]]) & agree(hi_s, bs[["upper", exact = TRUE]])))
      if (length(bad))
        stop(what, "the card's bounds_sampled disagree with its transformed bounds_natural for input ",
             bad[1], " ('", ctx[bad[1]], "'): its context_names may be permuted or mislabelled")
      lo_s <- bs[["lower", exact = TRUE]]; hi_s <- bs[["upper", exact = TRUE]]
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
  lay <- m[["input_layout", exact = TRUE]]
  if (is.null(lay)) stop(what, "a ", kind, " card must declare input_layout (the network's inputs in order)")
  if (kind == "mlp_joint" && !(is.numeric(m[["ll_floor_log", exact = TRUE]]) && length(m[["ll_floor_log", exact = TRUE]]) == 1L && is.finite(m[["ll_floor_log", exact = TRUE]])))
    stop(what, "a mlp_joint card must give ll_floor_log, the finite log density returned outside the training box")
  lay <- as.character(lay)
  bad <- setdiff(lay, c(ctx, nn_data_inputs))
  if (length(bad)) stop(what, "input_layout entries ", paste(bad, collapse = ", "),
                        " are neither context_names nor one of ", paste(nn_data_inputs, collapse = ", "))
  if (anyDuplicated(lay) || !all(ctx %in% lay))
    stop(what, "input_layout must list every context name exactly once")
  n_in <- nrow(m[["mlp", exact = TRUE]][["layers", exact = TRUE]][[1]][["W", exact = TRUE]])
  if (n_in != length(lay)) stop(what, "input_layout has ", length(lay), " entries; the network has ", n_in, " inputs")
  if (ncol(m[["mlp", exact = TRUE]][["layers", exact = TRUE]][[length(m[["mlp", exact = TRUE]][["layers", exact = TRUE]])]][["W", exact = TRUE]]) != 1L)
    stop(what, "a ", kind, " network must have a single output (the joint log density)")
  if (!is.null(m[["scaler", exact = TRUE]]) && (length(m[["scaler", exact = TRUE]][["mean", exact = TRUE]]) != n_in || length(m[["scaler", exact = TRUE]][["scale", exact = TRUE]]) != n_in))
    stop(what, "a ", kind, " scaler needs one mean and scale per input (", n_in, ")")
  if (!is.null(m[["response_values", exact = TRUE]]) && length(m[["response_values", exact = TRUE]]) != 2L)
    stop(what, "response_values must give the input values for responses 1 and 2")
  if (!is.null(m[["output_scaler", exact = TRUE]]) && !(length(m[["output_scaler", exact = TRUE]][["mean", exact = TRUE]]) == 1L && length(m[["output_scaler", exact = TRUE]][["scale", exact = TRUE]]) == 1L))
    stop(what, "output_scaler needs one mean and one scale")
  nle_mlp_forward(m[["mlp", exact = TRUE]], matrix(0, 1L, n_in))     # refuses an unsupported activation now
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
    bound = list(minmax = minmax, exception = reg$exception),
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
  if (reg$kind == "flow_cat") {
    # K categories = the two response levels x K / 2 ratings (rating_fold); the
    # rating column RR is an extra response, censoring codes 1 and 2 are scored
    # by the evaluator (nn_loglik_cat), the convolution parameters t0 / st0 come
    # from the model's parameter matrix
    ml$extra_responses <- "RR"
    ml$check_data <- function(data) nn_cat_check_data(data, reg)
    ml$dfun <- function(rt, R, RR, pars) nn_eval_cat(reg, rt, R, RR, pars)$pdf
    ml$pfun <- function(rt, R, RR, pars) nn_eval_cat(reg, rt, R, RR, pars)$cdf
    ml$log_likelihood <- function(pars, dadm, model, min_ll = log(1e-10)) nn_loglik_cat(pars, dadm, reg, min_ll)
  } else if (isTRUE(reg$scaled)) {
    ml$dfun <- function(rt, R, pars) nn_eval_scaled(reg, rt, R, pars, "pdf")
    ml$pfun <- function(rt, R, pars) nn_eval_scaled(reg, rt, R, pars, "cdf")
    ml$log_likelihood <- function(pars, dadm, model, min_ll = log(1e-10)) {
      nn_check_data(dadm, reg)
      log_likelihood_ddm(pars = pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  } else if (joint) {
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
  # the compiled evaluator of a scaled card is the joint flow's, with the
  # `scaled` rule (src/model_NN.h): it computes the inputs from sc_pars
  list(kind = if (isTRUE(reg$scaled)) "flow_joint" else reg$kind, label = reg$label, pars = reg$pars,
       transform_codes = reg$transform_codes,
       zerobox_c = reg$zerobox_c,
       scaled = isTRUE(reg$scaled), sc_pars = if (isTRUE(reg$scaled)) reg$sc_pars else character(0),
       pre = if (is.null(reg$pre)) "" else reg$pre,
       conv_mean = if (is.null(reg$conv)) "" else unname(reg$conv[["mean"]]),
       conv_range = if (is.null(reg$conv)) "" else unname(reg$conv[["range"]]),
       hybrid = if (is.null(reg$hybrid)) "" else reg$hybrid,
       wald_pars = reg$wald$pars, wald_lower = reg$wald$lower, wald_upper = reg$wald$upper,
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
  zc <- if (is.null(reg$zerobox_c)) rep(NA_real_, length(reg$pars)) else reg$zerobox_c
  for (j in seq_along(reg$pars)) theta[, j] <- nn_tf(pars[, reg$pars[j]], reg$transforms[j], zc[j])
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
  if (!is.null(reg$hybrid)) {
    # rows with the hybrid parameter exactly 0: the Wald, inside the box on its
    # parameters (closed, as the evaluator's); they never reach the network
    wd <- use & !is.na(pars[, reg$hybrid]) & pars[, reg$hybrid] == 0
    if (any(wd)) {
      wp <- pars[wd, reg$wald$pars, drop = FALSE]
      inb <- rowSums(t(t(wp) >= reg$wald$lower & t(wp) <= reg$wald$upper)) == 3L
      inb[is.na(inb)] <- FALSE
      val <- numeric(sum(wd))
      if (any(inb)) val[inb] <- nle_wald_cpp(tn[wd][inb], wp[inb, 1L], wp[inb, 2L], wp[inb, 3L])[[what]]
      out[wd] <- val
      use <- use & !wd
      if (!any(use)) return(out)
    }
  }
  out[use] <- flow_eval_trials_cpp(nn_ptr(reg), nn_context(pars[use, , drop = FALSE], reg), tn[use])[[what]]
  out
}

# K-category flow (flow_cat). Category k = rating_fold(R, RR, K / 2): the first
# response level counts its ratings down, the second up (NEW.3 .. OLD.3 for
# K = 6). Rows with R or RR missing are censored rows (scored from LC / UC).
nn_cat_category <- function(R, RR, K) {
  Kr <- K %/% 2
  if (2L * Kr != K) stop("rating_fold needs an even number of categories")
  R <- as.integer(R); RR <- as.integer(round(RR))
  k <- ifelse(R == 1L, Kr - RR + 1L, Kr + RR)
  k[is.na(R) | is.na(RR)] <- NA_integer_
  k
}

# The evaluator's inputs: for a card with a `mix` block the natural-scale
# sources (the evaluator derives the network inputs itself), else the network
# inputs on the network's scale.
nn_cat_theta <- function(pars, reg) {
  if (!is.null(reg$card$mix)) {
    miss <- setdiff(reg$pars, colnames(pars))
    if (length(miss)) stop("Neural likelihood '", reg$label, "' needs parameter(s) ", paste(miss, collapse = ", "), " which the model does not supply")
    return(unname(as.matrix(pars[, reg$pars, drop = FALSE])))
  }
  nn_context(pars, reg)
}

# dfun / pfun of a flow_cat model on observed rows: the convolved joint density
# of (k, rt) and the convolved defective CDF of category k at rt.
nn_eval_cat <- function(reg, rt, R, RR, pars) {
  k <- nn_cat_category(R, RR, reg$n_categories)
  n <- length(rt)
  out <- cat_eval_trials_cpp(nn_ptr(reg), nn_cat_theta(pars, reg), as.numeric(rt),
                             as.numeric(pars[, reg$conv[["mean"]]]), as.numeric(pars[, reg$conv[["range"]]]),
                             k, integer(n), numeric(n), rep(Inf, n))
  list(pdf = exp(out$log_pdf), cdf = out$cdf)
}

# The R path of a flow_cat model (the reference for the compiled branch in
# src/model_NN.h): trial log-likelihoods over the compressed dadm, censored rows
# from LC / UC, expanded and floored as log_likelihood_ddm does.
nn_loglik_cat <- function(pars, dadm, reg, min_ll = log(1e-10)) {
  n <- nrow(dadm)
  miss <- dadm[["missingness"]]; if (is.null(miss)) miss <- rep(NA_integer_, n)
  miss[is.na(miss)] <- 0L
  LC <- dadm[["LC"]]; if (is.null(LC)) LC <- rep(0, n)
  UC <- dadm[["UC"]]; if (is.null(UC)) UC <- rep(Inf, n)
  k <- nn_cat_category(dadm$R, dadm$RR, reg$n_categories)
  k[miss != 0L] <- NA_integer_
  rt <- dadm$rt; rt[is.na(rt)] <- 0
  ok <- attr(pars, "ok"); if (is.null(ok)) ok <- rep(TRUE, n)
  ll <- rep(min_ll, n)
  if (any(ok)) {
    P <- pars[ok, , drop = FALSE]
    out <- cat_eval_trials_cpp(nn_ptr(reg), nn_cat_theta(P, reg), as.numeric(rt[ok]),
                               as.numeric(P[, reg$conv[["mean"]]]), as.numeric(P[, reg$conv[["range"]]]),
                               k[ok], as.integer(miss[ok]), as.numeric(LC[ok]), as.numeric(UC[ok]))
    ll[ok] <- out$log_pdf
  }
  ll[is.na(ll)] <- min_ll
  sum(pmax(min_ll, ll[attr(dadm, "expand")]))
}

# Data checks of a flow_cat model: two response levels, a rating RR in
# 1..K/2 on every observed row, censoring codes 1 / 2 only (as src/model_NN.h),
# no truncation.
nn_cat_check_data <- function(data, reg) {
  what <- paste0("Neural likelihood '", reg$label, "': ")
  if (length(levels(data$R)) != 2L) stop(what, "needs exactly two response levels")
  if (!"RR" %in% names(data)) stop(what, "needs a rating column RR")
  miss <- data[["missingness"]]; if (is.null(miss)) miss <- rep(NA_integer_, nrow(data))
  if (any(!is.na(miss) & !(miss %in% c(1L, 2L))))
    stop(what, "only missingness codes 1 (lower-censored) and 2 (upper-censored) are supported")
  obs <- is.na(miss)
  if (any(is.na(data$R[obs]))) stop(what, "R is NA on a trial that is not censored")
  RR <- data$RR[obs]
  if (any(is.na(RR)) || any(RR != round(RR)) || any(RR < 1 | RR > reg$n_categories %/% 2))
    stop(what, "RR must be an integer in 1..", reg$n_categories %/% 2, " on every observed row")
  lt <- data[["LT"]]; ut <- data[["UT"]]
  if ((!is.null(lt) && any(lt > 0 & is.finite(lt), na.rm = TRUE)) || (!is.null(ut) && any(is.finite(ut))))
    stop(what, "truncation is not supported")
  invisible(TRUE)
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

# Inputs whose training region is closed at its lower edge (registrations made
# before `closed` existed carry only `zerobox`).
nn_closed <- function(reg) if (is.null(reg$closed)) reg$zerobox else reg$closed

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

# make_emc(): a neural likelihood with t0 sampled and the response times
# floored to a resolution coarser than 1 ms. EMC2's default resolution
# (1/60 s; applied by design_model() whether or not the data are compressed)
# biases t0 by about -1.9 posterior sd at N = 100 even for the analytic
# likelihood (the NLE project's measurement); the networks were validated at
# 1e-3 and finer.
nn_rt_resolution_warning <- function(design, rt_resolution) {
  if (!is.function(design$model)) return(invisible(FALSE))
  reg <- tryCatch(design$model()$nn, error = function(e) NULL)
  if (is.null(reg) || is.null(rt_resolution) || !(rt_resolution > 1e-3 + 1e-12)) return(invisible(FALSE))
  lhs <- vapply(design$Flist, function(f) as.character(stats::terms(f)[[2]]), "")
  if (!("t0" %in% lhs)) return(invisible(FALSE))
  warning("Neural likelihood '", reg$label, "': the response times are floored to a resolution of ",
          format(rt_resolution), " s (make_emc's rt_resolution, ", format(1 / rt_resolution, digits = 3),
          " Hz) while t0 is sampled. That floors every rt by up to ", format(1000 * rt_resolution, digits = 3),
          " ms and biases t0 (about -1.9 posterior sd at N = 100 with the default 1/60 s, even for the ",
          "analytic likelihood; compress = FALSE does not change the resolution). Use rt_resolution = 0.001 ",
          "or finer.", call. = FALSE)
  invisible(TRUE)
}

# Censoring and truncation are not wired for neural likelihoods; both paths
# refuse such data (src/model_NN.h makes the same checks). The K-category flow
# scores lower / upper censoring itself (nn_cat_check_data).
nn_check_data <- function(dadm, reg) {
  what <- paste0("Neural likelihood '", reg$label, "': ")
  miss <- dadm[["missingness"]]
  if (!is.null(miss) && any(!is.na(miss)) && !identical(reg$kind, "flow_cat"))
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
    if (!nn_inside(p, nat[[p]], reg$lower, reg$upper, nn_closed(reg))) {
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
  step <- if (!is.null(cm[["checkpoint_step", exact = TRUE]])) cm[["checkpoint_step", exact = TRUE]] else
    c(flow = cm[["checkpoint_step_flow", exact = TRUE]], classifier = cm[["checkpoint_step_classifier", exact = TRUE]])
  cat(sprintf("EMC2 %s @ %s\n  artefact %s: budget=%s checkpoint_step=%s %s=%s\n",
              as.character(utils::packageVersion("EMC2")), find.package("EMC2"),
              reg$label, if (is.null(cm[["budget", exact = TRUE]])) "?" else cm[["budget", exact = TRUE]],
              if (is.null(names(step))) paste(step, collapse = "/") else
                paste(names(step), step, sep = "=", collapse = "/"),
              nle_hash_algo(), reg$sha256))
  invisible(reg)
}
