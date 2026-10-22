# Neural-likelihood (normalizing-flow) variants of the DDM and the RDM:
# DDMnn and RDMnn, thin wrappers over register_nn_model() (R/nn_register.R).
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
#' The network is exact at sv = 0 and SZ = 0 (zerobox encoding); st0 is fixed
#' at 0 (the card was trained without st0 variability). An un-sampled `sv` or
#' `SZ` at its default 0 is accepted as an exact zero; a sampled or non-zero
#' `st0` is refused with an error.
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
