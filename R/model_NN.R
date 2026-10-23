# Neural-likelihood (normalizing-flow) variants of the DDM and the RDM:
# DDMnn and RDMnn, thin wrappers over register_nn_model() (R/nn_register.R);
# and RTCONnn, the RTCON2 rating model as a K-category flow (src/flow_cat.cpp,
# trained by the RTCON2 project).
#
# The likelihoods come from rational-quadratic spline flows trained on
# simulated data (Wuth, 2026, MSc thesis, UvA; NLE project). The trained
# artefacts ship as inst/extdata/flownn/*.json, are pinned by sha256 in
# inst/extdata/flownn/MANIFEST, and are evaluated by compiled code
# (src/flow_ddm.cpp, src/flow_race.cpp, derived from the NLE handover package,
# which defines the models; the batched evaluator core is src/nle_flow.h).

#' The Diffusion Decision Model — Neural-Network Likelihood (Scaled)
#'
#' Variant of [DDM] whose joint likelihood of response and rt is computed by
#' a scaled normalizing flow (neural likelihood estimation, kind
#' `"flow_scaled"`): the network models the unit-scale decision time given
#' kappa = v*a/s^2, eta = sv*a/s^2, Z and SZ, with an exact zero at
#' sv = SZ = 0. Evaluation takes microseconds per trial.
#'
#' @details
#'
#' Parameters, transforms, and defaults match [DDM]; the natural-scale ranges
#' are the flow's **training region** (a bounding hyper-box read from the
#' artefact): parameter vectors outside it are rejected (likelihood floored at
#' `min_ll`).
#'
#' The defaults are **identical to [DDM]'s**, including `sv = SZ = st0 = 0`.
#' The network is exact at sv = 0 and SZ = 0 (zerobox encoding); an un-sampled
#' `sv` or `SZ` at its default 0 is accepted as an exact zero. `st0` is not a
#' network input: the uniform non-decision-time window is applied exactly
#' outside the network, as a difference of the flow's CDF over the window
#' (`st0 = 0` gives the density itself), so `st0` may be sampled or held at
#' any value, 0 included.
#'
#' Loading a neural likelihood is not the same as being able to infer with it:
#' check calibration and identifiability (parameter recovery, likelihood
#' profiles, SBC) before drawing conclusions; [nn_cell()] and the validation
#' kit run these checks against [DDM].
#'
#' `DDMnn(artefact)` is `register_nn_model(artefact, twin = DDM)()`; see
#' [register_nn_model()] for the contract.
#'
#' Wuth, J. (2026). *Likelihood approximation in evidence accumulation
#' models: A comparison of kernel density and neural likelihood methods*
#' (MSc thesis, University of Amsterdam).
#'
#' @param artefact Name of a shipped artefact in `inst/extdata/flownn`
#'   (default `"ddm_sc_s1106_e245_weights.json"`, the DDMnle scaled card).
#' @return A model list with all the necessary functions for EMC2 to sample
#' @export
DDMnn <- function(artefact = "ddm_sc_s1106_e245_weights.json")
  register_nn_model(artefact, twin = DDM)()

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
#' `A >= 1e-4`, `v >= 1e-3`). An un-sampled `A` at its default `0` is refused
#' with an error; fix it inside the box instead, e.g.
#' `constants = c(A = log(1e-3))`. Out-of-box parameter vectors are rejected
#' (likelihood floored at `min_ll`).
#'
#' `RDMnn(artefact)` is `register_nn_model(artefact, kind = "flow_race",
#' twin = RDM)()`; see [register_nn_model()] for the contract.
#'
#' Data generation (`rfun`) uses the analytic RDM, so simulation-based
#' checks compare flow-based inference against exact data. Loading a neural
#' likelihood is not the same as being able to infer with it; validate before
#' use ([nn_cell()] and the validation kit run the checks against [RDM]).
#'
#' Wuth, J. (2026). *Likelihood approximation in evidence accumulation
#' models: A comparison of kernel density and neural likelihood methods*
#' (MSc thesis, University of Amsterdam).
#'
#' @param artefact Name of a shipped artefact in `inst/extdata/flownn`
#'   (default `"rdm_ss_base_s801_e271_weights.json"`).
#' @return A model list with all the necessary functions for EMC2 to sample
#' @export
RDMnn <- function(artefact = "rdm_ss_base_s801_e271_weights.json")
  register_nn_model(artefact, kind = "flow_race", twin = RDM)()

# RTCON2's parameterisation: types, transforms and bounds as the RTCON2
# project's twin (rating-work/rtcon2, R2/emc2/rtcon2_model.R), without its
# simulator (C code outside EMC2). The bounds are the training box of the
# sources (natural scale); t0 and st0 (seconds) are not network inputs (the
# convolution applies them exactly) and keep the box the project simulated.
rtcon2_twin <- function() {
  ctx <- c("mu", "s", "C1", "g2", "g3", "g4", "g5", "a", "sb", paste0("B", 1:6))
  lower <- c(mu = -2.5, s = .2, C1 = -3, g2 = .25, g3 = .25, g4 = .25, g5 = .25, a = .01, sb = 0,
             B1 = .8, B2 = .8, B3 = .8, B4 = .8, B5 = .8, B6 = .8, t0 = .15, st0 = 0)
  upper <- c(mu = 6, s = 3.5, C1 = .5, g2 = 2.5, g3 = 2.5, g4 = 2.5, g5 = 2.5, a = .3, sb = 8,
             B1 = 6, B2 = 6, B3 = 6, B4 = 6, B5 = 6, B6 = 6, t0 = .7, st0 = .3)
  func <- c(mu = "identity", s = "exp", C1 = "identity", g2 = "identity", g3 = "identity", g4 = "identity",
            g5 = "identity", a = "exp", sb = "exp", stats::setNames(rep("exp", 6), paste0("B", 1:6)),
            t0 = "exp", st0 = "exp")
  list(
    type = "DDM",
    p_types = c(mu = 0, s = log(1), C1 = -1, g2 = 0.8, g3 = 0.8, g4 = 0.8, g5 = 0.8,
                a = log(0.04), sb = log(0.8),
                B1 = log(2), B2 = log(2), B3 = log(2), B4 = log(2), B5 = log(2), B6 = log(2),
                t0 = log(0.38), st0 = log(0)),
    transform = list(func = func),
    bound = list(minmax = rbind(lower, upper), exception = c(st0 = 0, sb = 0)))
}

#' RTCON2 — Neural-Network Likelihood
#'
#' The RTCON2 confidence-rating model (Ratcliff and Starns, 2013) for two
#' responses with three ratings each, whose likelihood is a K-category
#' neural likelihood (kind `"flow_cat"`, K = 6): a classifier for the
#' category plus a spline flow for the decision time, trained on the RTCON2
#' simulator (the RTCON2 project, stage R3).
#'
#' @details
#'
#' The data need the response `R` (two levels), `rt` and the rating `RR`
#' (1..3; an extra response of the model). The category is
#' `k = rating_fold(R, RR)`: the first response level counts its ratings down,
#' the second up (`R1 RR3, R1 RR2, R1 RR1, R2 RR1, R2 RR2, R2 RR3` = k 1..6,
#' e.g. NEW.3 .. OLD.3).
#'
#' The parameters (natural scale) are the design cell's strength mean `mu` and
#' SD `s`, the first criterion `C1` and the criterion gaps `g2..g5` (criteria
#' `cumsum(C1, g2..g5)`, with the simulator's minimum gap 0.25), the drift
#' scale `a`, the boundary variability `sb` and the six boundaries `B1..B6`,
#' and the non-decision mean `t0` and range `st0` (seconds). The trial
#' strength `w ~ N(mu, s^2)` is integrated by 64-node Gauss-Hermite
#' quadrature outside the network, and `t0` / `st0` are applied by an exact
#' uniform convolution, so `st0 = 0` (the default) is exact. Lower and upper
#' censored trials ([make_missing()], codes 1 and 2) are scored by the
#' network's CDF.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default** |
#' |-----------|-----------|---------------|------------|
#' | *mu*      | identity  | \[-2.5, 6\]   | 0          |
#' | *s*       | log       | \[0.2, 3.5\]  | log(1)     |
#' | *C1*      | identity  | \[-3, 0.5\]   | -1         |
#' | *g2..g5*  | identity  | \[0.25, 2.5\] | 0.8        |
#' | *a*       | log       | \[0.01, 0.3\] | log(0.04)  |
#' | *sb*      | log       | \[0, 8\]      | log(0.8)   |
#' | *B1..B6*  | log       | \[0.8, 6\]    | log(2)     |
#' | *t0*      | log       | \[0.15, 0.7\] | log(0.38)  |
#' | *st0*     | log       | \[0, 0.3\]    | log(0)     |
#'
#' The ranges of `mu` to `B6` are the network's training region: parameter
#' vectors outside it are rejected (likelihood floored at `min_ll`); `sb = 0`
#' and `st0 = 0` are exact. There is no simulator in EMC2 (the RTCON2
#' simulator is C code in the RTCON2 project); to simulate, register the card
#' with `register_nn_model(artefact, twin = <a model function with an rfun>)`.
#'
#' Loading a neural likelihood is not the same as being able to infer with it:
#' check parameter recovery and calibration before drawing conclusions.
#'
#' Ratcliff, R., & Starns, J. J. (2013). Modeling confidence judgments,
#' response times, and multiple choices in decision making: Recognition
#' memory and motion discrimination. *Psychological Review, 120*(3), 697-719.
#'
#' @param artefact Name of a shipped artefact in `inst/extdata/flownn`
#'   (default `"rtcon2_r3_s2_e006_weights.json"`), or a card path.
#' @return A model list with all the necessary functions for EMC2 to sample
#' @examples
#' RTCONnn()$nn$kind
#' @export
RTCONnn <- function(artefact = "rtcon2_r3_s2_e006_weights.json")
  register_nn_model(artefact, kind = "flow_cat", twin = rtcon2_twin)()
