# The scaled DDM network as an EMC2 model (handover/SCALED_DDM_PLAN.md).
#
# The network models the unit-scale decision time u = (rt - t0) * (s/a)^2 of a diffusion with
# boundary 1 and noise 1, conditioned on four inputs: kappa = v*a/s^2, eta = sv*a/s^2, Z, SZ.
# Everything else is exact:
#   f(rt) = c * g(u; kappa, eta, Z, SZ),  c = (s/a)^2                      (st0 = 0)
#   f(rt) = [G(u1) - G(u0)] / st0,  u1 = (rt - t0) c,  u0 = max(rt - t0 - st0, 0) c   (st0 > 0)
# with g, G the network's defective density and CDF, G(u) = P(R) * Phi(z(u)), which is exactly the
# integral of g because the flow lives on u > 0. The st0 window is EMC2's lower-edge convention
# [t0, t0 + st0] (sbc/st0conv_screen.R's window check).
#
# DDMsc(card) returns a model function: DDM() with its likelihood replaced. Transforms, bounds,
# the zero exceptions for sv, SZ and st0, and the simulator (rfun) are DDM's own, so designs,
# priors and make_data behave exactly as for the analytic DDM. c_name is dropped, so EMC2 takes
# its R likelihood path (log_likelihood_ddm -> model$dfun). The compiled flow evaluator of EMC2
# dev-nle is called unchanged.
#
#   source("port/R/ddm_scaled.R")
#   m <- DDMsc("port/export/ddm_sc_s1_weights.json")
#   design(model = m, ...)
#
# DDMsc(evaluator = ddm_sc_oracle) swaps the network for the analytic unit-scale DDM: the
# test that the scaling and the convolution are exact (tests/test_ddm_scaled.R).

DDM_SC_INPUTS <- c("kappa", "eta", "Z", "SZ")
DDM_SC_C <- c(eta = 0.01, SZ = 0.05)        # asinh constants (simulation/scaled.R, flow data.py)

# Natural parameters (dfun's scale: after DDM's Ttransform, so pars[, "SZ"] is the absolute
# relative width sw = 2 SZ min(Z, 1 - Z)) -> the four unit-scale parameters, natural scale.
ddm_sc_unit <- function(pars) {
  s2 <- pars[, "s"]^2; a <- pars[, "a"]; Z <- pars[, "Z"]
  m <- 2 * pmin(Z, 1 - Z)
  cbind(kappa = pars[, "v"] * a / s2, eta = pars[, "sv"] * a / s2, Z = Z,
        SZ = ifelse(pars[, "SZ"] > 0, pars[, "SZ"] / m, 0))
}

# Unit-scale natural parameters -> the network's input scale (what bounds_sampled is on).
ddm_sc_net <- function(U)
  cbind(U[, "kappa"], asinh(U[, "eta"] / DDM_SC_C[["eta"]]), qnorm(U[, "Z"]),
        asinh(U[, "SZ"] / DDM_SC_C[["SZ"]]))

# An evaluator is function(U, u, R) -> list(pdf, cdf): the unit-scale defective density and CDF
# at u > 0 for response R (1 lower, 2 upper), one row of U per element of u.
ddm_sc_card_evaluator <- function(card) {
  card <- normalizePath(card, mustWork = TRUE)
  env <- new.env(parent = emptyenv())
  ptr <- function() {     # rebuilt after a fork-unsafe reload or in a new process
    if (is.null(env$ptr) || !EMC2:::ddm_ptr_valid(env$ptr)) {
      raw <- jsonlite::fromJSON(card, simplifyDataFrame = FALSE)
      if (!identical(as.character(raw$context_names), DDM_SC_INPUTS))
        stop("DDMsc: ", card, " is not a scaled-DDM card (context_names ",
             paste(raw$context_names, collapse = ", "), ")")
      env$ptr <- EMC2:::ddm_build_ensemble(EMC2:::nn_normalise_card(raw)$members)
    }
    env$ptr
  }
  f <- function(U, u, R) {
    o <- EMC2:::ddm_ens_eval_trials_cpp(ptr(), ddm_sc_net(U), u, as.integer(R))
    list(pdf = o[["pdf"]], cdf = o[["cdf"]])
  }
  attr(f, "card") <- card
  f
}

# The analytic unit-scale DDM in the same interface (a = s = 1, t0 = st0 = 0).
ddm_sc_oracle <- function(U, u, R, precision = 1e-8) {
  P <- cbind(v = U[, "kappa"], a = 1, sv = U[, "eta"], t0 = 0, st0 = 0, s = 1, Z = U[, "Z"],
             SZ = 2 * U[, "SZ"] * pmin(U[, "Z"], 1 - U[, "Z"]))
  Rf <- factor(R, levels = 1:2)
  list(pdf = EMC2:::dDDM(u, Rf, P, precision = precision),
       cdf = EMC2:::pDDM(u, Rf, P, precision = precision))
}

# The likelihood of rt (seconds) given natural parameters, through an evaluator.
ddm_sc_density <- function(rt, R, pars, evaluator) {
  n <- length(rt); out <- numeric(n)
  U <- ddm_sc_unit(pars)
  c2 <- (pars[, "s"] / pars[, "a"])^2
  d <- rt - pars[, "t0"]; st0 <- pars[, "st0"]
  live <- d > 0
  if (!any(live)) return(out)
  # st0 = 0, or a window too narrow to difference: the density (at the window's midpoint)
  narrow <- live & (st0 <= 0 | st0 < 1e-7 * d)
  wide <- live & !narrow
  if (any(narrow)) { i <- which(narrow)
    um <- (d[i] - pmax(st0[i], 0) / 2) * c2[i]
    out[i] <- c2[i] * evaluator(U[i, , drop = FALSE], um, R[i])$pdf }
  if (any(wide)) { i <- which(wide)
    u1 <- d[i] * c2[i]; u0 <- (d[i] - st0[i]) * c2[i]
    pos <- u0 > 0                                    # window entirely after onset: two CDFs
    j <- c(seq_along(i), which(pos))
    G <- evaluator(U[i[j], , drop = FALSE], c(u1, u0[pos]), R[i[j]])$cdf
    G0 <- numeric(length(i)); G0[pos] <- G[length(i) + seq_len(sum(pos))]
    out[i] <- pmax(G[seq_along(i)] - G0, 0) / st0[i] }
  out
}

DDMsc <- function(card = NULL, evaluator = NULL) {
  if (is.null(evaluator)) evaluator <- ddm_sc_card_evaluator(card)
  force(evaluator)
  function() {
    m <- DDM()
    m$c_name <- NULL                 # EMC2's R likelihood path: log_likelihood_ddm -> dfun
    m$dfun <- function(rt, R, pars) ddm_sc_density(rt, as.integer(R), pars, evaluator)
    m$pfun <- NULL                   # no CDF of rt is offered (not needed to fit)
    m$sc_card <- attr(evaluator, "card")
    m
  }
}

# The compiled path (EMC2 branch dev-nle-scaled, src/model_NN.h nn_trial_ll_scaled): a copy of a
# DDM-parameter registration `reg` (e.g. the scaffold card's) whose native spec hands EMC2 the
# scaled card and the eight DDM columns. Use as model list fields c_name = "NN", nn = <this>.
# The R reference above (ddm_sc_density) stays the definition; the two agree to rounding.
ddm_sc_native_reg <- function(reg, card) {
  card <- normalizePath(card, mustWork = TRUE)
  raw <- jsonlite::fromJSON(card, simplifyDataFrame = FALSE)
  if (!identical(as.character(raw$context_names), DDM_SC_INPUTS))
    stop("ddm_sc_native_reg: ", card, " is not a scaled-DDM card")
  tfs <- vapply(DDM_SC_INPUTS, function(p) as.character(raw$context_transforms[[p]]), "")
  zc <- vapply(DDM_SC_INPUTS, function(p) if (is.null(raw$zerobox_c[[p]])) NA_real_ else raw$zerobox_c[[p]], 0)
  stopifnot(identical(unname(tfs), c("identity", "zerobox", "probit", "zerobox")),
            isTRUE(all.equal(unname(zc[c("eta", "SZ")]), unname(DDM_SC_C[c("eta", "SZ")]))))
  r <- reg
  r$path <- card; r$artefact <- NULL; r$kind <- "flow_joint"
  r$sha256 <- EMC2:::nn_card(list(path = card))$sha256
  nat <- reg$native
  r$native <- list(kind = "flow_joint", label = paste0("DDMnle ", basename(card)), pars = DDM_SC_INPUTS,
                   transform_codes = unname(EMC2:::nn_transform_codes[tfs]), zerobox_c = unname(zc),
                   pre = "", hybrid = "", wald_pars = character(0), wald_lower = numeric(0),
                   wald_upper = numeric(0), fixed = nat$fixed, fixed_msg = nat$fixed_msg,
                   refuse = nat$refuse, refuse_msg = nat$refuse_msg,
                   scaled = TRUE, sc_pars = c("v", "a", "sv", "t0", "st0", "s", "Z", "SZ"))
  r
}
