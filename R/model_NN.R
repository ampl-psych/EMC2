# Neural-likelihood (normalizing-flow) variants of the DDM and the RDM:
# DDMnn and RDMnn, thin wrappers over register_nn_model() (R/nn_register.R).
#
# The likelihoods come from rational-quadratic spline flows trained on
# simulated data (Wuth, 2026, MSc thesis, UvA; NLE project). The trained
# artefacts ship as inst/extdata/flownn/*.rds, are pinned by sha256 in
# inst/extdata/flownn/MANIFEST, and are evaluated by compiled code
# (src/flow_ddm.cpp, src/flow_race.cpp, derived from the NLE handover package,
# which defines the models; the batched evaluator core is src/nle_flow.h).

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
#' `constants = c(sv = log(0.5))`). The refusal happens in [design()] (and again
#' at run time). The `"ddm_st0zero"` artefact was trained without `st0`: it
#' admits exactly `st0 = 0` and refuses a sampled or non-zero `st0`.
#'
#' Unlike [DDM], `SZ` here is the raw proportion the flow was trained on, not
#' remapped by `2 * SZ * min(Z, 1 - Z)`; data generation applies that remapping
#' internally so `rfun` matches [DDM]'s parameter meaning exactly. Data are
#' generated with the exact sampler (`WienR::rWDM`).
#'
#' Loading a neural likelihood is not the same as being able to infer with it:
#' check calibration and identifiability (parameter recovery, likelihood
#' profiles, SBC) before drawing conclusions; [nn_cell()] and the validation
#' kit run these checks against [DDM].
#'
#' **Calibration of `ddm_cap256w_c4`, and the t0 shift.** The artefact passes
#' simulation-based calibration from N = 100 to 800 trials on every parameter
#' (standardised bias within the 0.5 gate, coverage nominal), with one
#' systematic feature: its t0 point estimates sit slightly high, by a shift
#' that is constant on the log scale (about 0.002 at N = 800, i.e. 0.2% of
#' t0; rank uniformity fails from N = 200 while interval coverage stays
#' nominal). **No correction is applied by EMC2.** A fixed multiplicative
#' correction of 0.00276 was proposed by the NLE project and reopened on
#' 2026-10-13, when an exact test with saved draws found it over-corrects at
#' N = 800 (mean rank .534, KS p .007); an out-of-sample test of about 0.0019
#' is pending, and the shift is reported here as unresolved. For a DDM with
#' `sv`, `SZ` or `st0` free the scaled DDM card (DDMnle, a `"flow_scaled"`
#' card registered with [register_nn_model()]) replaces this artefact: it
#' needs no t0 correction (t0 inside the envelope in all six variability
#' models at N = 400). The record of the exact test is the NLE project's
#' `sbc/t0corr_exact.R` and the progress log of
#' `handover/ZERO_BOX_RETRAIN_PLAN.md`.
#'
#' `DDMnn(artefact)` is `register_nn_model(artefact, kind = "flow_joint",
#' twin = DDM)()`; see [register_nn_model()] for the contract.
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
DDMnn <- function(artefact = "ddm_cap256w_c4")
  register_nn_model(artefact, kind = "flow_joint", twin = DDM)()

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
#'   (default `"rdm_small"`).
#' @return A model list with all the necessary functions for EMC2 to sample
#' @export
RDMnn <- function(artefact = "rdm_small")
  register_nn_model(artefact, kind = "flow_race", twin = RDM)()
