# Neural-network (normalizing-flow) likelihood variants of the race models:
# LNRnn and RDMnn.
#
# The single-accumulator density/survivor come from a rational-quadratic
# spline flow trained on simulated data (Wuth, 2026, MSc thesis, UvA;
# likelihood-kde-normflow repo). The trained networks ship with the package
# as inst/extdata/flownn/*.rds and are evaluated by compiled code
# (src/flow_race.cpp). One spline inversion yields pdf, cdf, and survivor.
#
# The flows are only valid inside their training region. That box is
# enforced twice: as the model's `bound` (EMC2 rejects proposals outside it,
# via the usual `ok` mechanism) and inside the C++ evaluator (out-of-box
# parameter vectors return pdf = 0 and cdf = 0).

# Per-session cache of compiled flow models. XPtrs do not survive
# save/load, so entries are rebuilt lazily when invalid.
flownn_cache <- new.env(parent = emptyenv())

flownn_get <- function(name) {
  ptr <- flownn_cache[[name]]
  if (is.null(ptr) || !flow_ptr_valid(ptr)) {
    path <- system.file("extdata", "flownn", paste0(name, ".rds"),
                        package = "EMC2", mustWork = TRUE)
    ptr <- flow_build(readRDS(path))
    flownn_cache[[name]] <- ptr
  }
  ptr
}

# DDM-character models use the two-network evaluator (src/flow_ddm.cpp)
ddmnn_get <- function(name) {
  ptr <- flownn_cache[[name]]
  if (is.null(ptr) || !ddm_ptr_valid(ptr)) {
    path <- system.file("extdata", "flownn", paste0(name, ".rds"),
                        package = "EMC2", mustWork = TRUE)
    ptr <- ddm_build(readRDS(path))
    flownn_cache[[name]] <- ptr
  }
  ptr
}

# Flow evaluators. pars arrives on the natural scale (after the model's
# transform); the flows were trained on the sampled scale, so positive
# parameters are logged (and proportions probit-transformed) here. Column
# order must match the flow's training context (recorded in the .rds under
# $context_names).

dLNRnn <- function(rt, pars) {
  theta <- cbind(pars[, "m"], log(pars[, "s"]), log(pars[, "t0"]))
  flow_eval_trials_cpp(flownn_get("lnr_large"), theta, rt)$pdf
}

pLNRnn <- function(rt, pars) {
  theta <- cbind(pars[, "m"], log(pars[, "s"]), log(pars[, "t0"]))
  flow_eval_trials_cpp(flownn_get("lnr_large"), theta, rt)$cdf
}

dRDMnn <- function(rt, pars) {
  theta <- log(cbind(pars[, "v"], pars[, "B"], pars[, "t0"],
                     pars[, "s"], pars[, "A"]))
  flow_eval_trials_cpp(flownn_get("rdm_large"), theta, rt)$pdf
}

# DDM: two networks (spline flow + choice classifier), joint defective
# density over (rt, R). The flow was trained on the RAW (pre-Ttransform)
# parameters, sampled scale: v identity; a, t0, s, sv, st0 log; Z, SZ
# probit. Context order: v, a, t0, s, Z, SZ, sv, st0. R codes the response
# factor level (1 = first level = lower boundary, as in dDDM).

.theta_DDMnn <- function(pars) {
  cbind(pars[, "v"], log(pars[, "a"]), log(pars[, "t0"]), log(pars[, "s"]),
        qnorm(pars[, "Z"]), qnorm(pars[, "SZ"]),
        log(pars[, "sv"]), log(pars[, "st0"]))
}

dDDMnn <- function(rt, R, pars) {
  ddm_eval_trials_cpp(ddmnn_get("ddm_large"), .theta_DDMnn(pars),
                      rt, as.integer(R))$pdf
}

pDDMnn <- function(rt, R, pars) {
  ddm_eval_trials_cpp(ddmnn_get("ddm_large"), .theta_DDMnn(pars),
                      rt, as.integer(R))$cdf
}

pRDMnn <- function(rt, pars) {
  theta <- log(cbind(pars[, "v"], pars[, "B"], pars[, "t0"],
                     pars[, "s"], pars[, "A"]))
  flow_eval_trials_cpp(flownn_get("rdm_large"), theta, rt)$cdf
}


#' The Log-Normal Race Model — Neural-Network Likelihood
#'
#' Variant of [LNR] whose single-accumulator density and survivor functions
#' are computed by a trained normalizing flow (neural likelihood estimation)
#' instead of the analytic log-normal. Intended for validating the
#' flow-likelihood workflow against a model whose analytic likelihood is
#' available for comparison.
#'
#' @details
#'
#' Parameters, transforms, and interpretation match [LNR]:
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default**   | **Interpretation**            |
#'  |-----------|-----------|---------------|-----------|---------------------------|
#'  | *m*       | -         | \[-3, 3\]       | 1         | Scale parameter           |
#'  | *s*       | log       | \[0.1, 2\]      | log(1)    | Shape parameter           |
#'  | *t0*      | log       | \[0.05, 1\]     | log(0.1)  | Non-decision time         |
#'
#' The natural-scale ranges above are the flow's **training region**
#' (a bounding hyper-box): parameter vectors outside it are rejected
#' (likelihood floored at `min_ll`), both through the model's `bound` and
#' inside the compiled evaluator. Note the box is much tighter than the
#' analytic [LNR] bounds — priors and start points must respect it.
#'
#' Data generation (`rfun`) uses the analytic LNR, so simulation-based
#' checks compare flow-based inference against exact data.
#'
#' Wuth, J. (2026). *Likelihood approximation in evidence accumulation
#' models: A comparison of kernel density and neural likelihood methods*
#' (MSc thesis, University of Amsterdam).
#'
#' @return A model list with all the necessary functions for EMC2 to sample
#' @export
LNRnn <- function() {
  list(
    type="RACE",
    c_name = NULL, # R-path race likelihood; flow evaluator is compiled
    p_types=c("m" = 1,"s" = log(1),"t0" = log(0.1)),
    transform=list(func=c(m = "identity",s = "exp", t0 = "exp")),
    # the flow's training region (natural scale)
    bound=list(minmax=cbind(m=c(-3,3), s = c(0.1,2), t0=c(0.05,1))),
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) pars,
    # Random function for racing accumulators: exact LNR data
    rfun=function(data=NULL,pars) rLNR(data$lR, pars, ok = attr(pars, "ok")),
    # Flow density (PDF) for single accumulator
    dfun=function(rt,pars) dLNRnn(rt,pars),
    # Flow probability (CDF) for single accumulator
    pfun=function(rt,pars) pLNRnn(rt,pars),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(pars,dadm, model, min_ll=log(1e-10))
      log_likelihood_race(pars=pars, dadm = dadm, model = model, min_ll = min_ll)
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
#' Parameters, transforms, and interpretation match [RDM]:
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default**    | **Interpretation**            |
#'  |-----------|-----------|---------------|------------|---------------------------|
#'  | *v*       | log       | \[0.001, 5\]    | log(1)     | Evidence-accumulation rate |
#'  | *B*       | log       | \[0.1, 3\]      | log(1)     | Threshold gap (b = B + A) |
#'  | *A*       | log       | \[0.0001, 2\]   | log(0.001) | Start-point variability   |
#'  | *t0*      | log       | \[0.05, 1\]     | log(0.1)   | Non-decision time         |
#'  | *s*       | log       | \[0.1, 2\]      | log(1)     | Within-trial noise        |
#'
#' The natural-scale ranges above are the flow's **training region**
#' (a bounding hyper-box): parameter vectors outside it are rejected
#' (likelihood floored at `min_ll`), both through the model's `bound` and
#' inside the compiled evaluator. Unlike [RDM], `A = 0` and `v = 0` are
#' **not** available (the flow trained on `A >= 1e-4`, `v >= 1e-3`); for a
#' no-start-point-variability model use a small fixed `A` instead, e.g.
#' `constants = c(A = log(1e-3))`.
#'
#' Data generation (`rfun`) uses the analytic RDM, so simulation-based
#' checks compare flow-based inference against exact data.
#'
#' Wuth, J. (2026). *Likelihood approximation in evidence accumulation
#' models: A comparison of kernel density and neural likelihood methods*
#' (MSc thesis, University of Amsterdam).
#'
#' @return A model list with all the necessary functions for EMC2 to sample
#' @export
RDMnn <- function() {
  list(
    type="RACE",
    c_name = NULL, # R-path race likelihood; flow evaluator is compiled
    p_types=c("v" = log(1),"B" = log(1),"A" = log(0.001),"t0" = log(0.1),"s" = log(1)),
    transform=list(func=c(v = "exp", B = "exp", A = "exp",t0 = "exp", s = "exp")),
    # the flow's training region (natural scale); no A = 0 / v = 0 exceptions
    bound=list(minmax=cbind(v=c(1e-3,5), B=c(0.1,3), A=c(1e-4,2),
                            t0=c(0.05,1), s=c(0.1,2))),
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) pars,
    # Random function for racing accumulators: exact RDM data
    rfun=function(data=NULL,pars) rRDM(data$lR,pars,ok=attr(pars, "ok")),
    # Flow density (PDF) for single accumulator
    dfun=function(rt,pars) dRDMnn(rt,pars),
    # Flow probability (CDF) for single accumulator
    pfun=function(rt,pars) pRDMnn(rt,pars),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10))
      log_likelihood_race(pars=pars, dadm = dadm, model = model, min_ll = min_ll)
  )
}


#' The Diffusion Decision Model — Neural-Network Likelihood
#'
#' Variant of [DDM] whose joint likelihood of response and rt is computed by
#' two trained neural networks (neural likelihood estimation) instead of the
#' numerical-integral density: a normalizing flow for the rt distribution
#' conditional on parameters and response, and a classifier for the choice
#' probability. Evaluation takes microseconds per trial — orders of
#' magnitude faster than the numerical DDM density.
#'
#' @details
#'
#' Parameters, transforms, and interpretation match [DDM]:
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default**    | **Interpretation**            |
#' |-----------|-----------|---------------|------------|---------------------------|
#' | *v*       | -         | \[-6, 6\]       | 1          | Mean drift rate |
#' | *a*       | log       | \[0.05, 5\]     | log(1)     | Boundary separation |
#' | *t0*      | log       | \[0.05, 1\]     | log(0.2)   | Non-decision time |
#' | *s*       | log       | \[0.1, 2\]      | log(1)     | Within-trial drift SD |
#' | *Z*       | probit    | \[0.1, 0.9\]    | qnorm(0.5) | Relative start point |
#' | *SZ*      | probit    | \[0.05, 0.99\]  | qnorm(0.1) | Start-point variability |
#' | *sv*      | log       | \[0.01, 4\]     | log(0.5)   | Between-trial drift SD |
#' | *st0*     | log       | \[0.01, 0.49\]  | log(0.05)  | Non-decision variability |
#'
#' The natural-scale ranges above are the flows' **training region**
#' (a bounding hyper-box): parameter vectors outside it are rejected
#' (likelihood floored at `min_ll`), both through the model's `bound` and
#' inside the compiled evaluator. Unlike [DDM], the simplified-model
#' exceptions `sv = 0`, `SZ = 0`, `st0 = 0` are **not** available — this is
#' intrinsically the *full* DDM; for near-zero variabilities fix them to
#' small constants inside the box (e.g., `constants = c(sv = log(0.01))`).
#'
#' Note that unlike [DDM], `SZ` here is the raw proportion (the flow's
#' training parameter), not remapped by `2 * SZ * min(Z, 1 - Z)`; data
#' generation applies that remapping internally so `rfun` matches [DDM]'s
#' parameter meaning exactly.
#'
#' Data generation (`rfun`) uses the exact sampler (`WienR::rWDM`), so
#' simulation-based checks compare flow-based inference against exact data.
#'
#' Wuth, J. (2026). *Likelihood approximation in evidence accumulation
#' models: A comparison of kernel density and neural likelihood methods*
#' (MSc thesis, University of Amsterdam).
#'
#' @return A model list with all the necessary functions for EMC2 to sample
#' @export
DDMnn <- function() {
  list(
    type="DDM",
    c_name = NULL, # R-path DDM likelihood; flow evaluator is compiled
    p_types=c("v" = 1,"a" = log(1),"sv" = log(0.5),"t0" = log(0.2),
              "st0" = log(0.05),"s" = log(1),"Z" = qnorm(0.5),"SZ" = qnorm(0.1)),
    transform=list(func=c(v = "identity",a = "exp",sv = "exp",t0 = "exp",
                          st0 = "exp",s = "exp",Z = "pnorm",SZ = "pnorm")),
    # the flows' training region (natural scale); no sv/SZ/st0 = 0 exceptions
    bound=list(minmax=cbind(v=c(-6,6),a=c(0.05,5),Z=c(.1,.9),t0=c(0.05,1),
                            sv=c(.01,4),s=c(0.1,2),SZ=c(.05,.99),st0=c(0.01,0.49))),
    # flow trained on raw parameters: keep SZ/Z unmodified (cf. DDM, which
    # rescales SZ and adds z/sz here)
    Ttransform = function(pars,dadm) pars,
    # Random function: exact DDM data; apply the DDM's SZ remapping locally
    # so parameters mean the same thing as in the analytic model
    rfun=function(data=NULL,pars) {
      pars[,"SZ"] <- 2*pars[,"SZ"]*pmin(pars[,"Z"], 1 - pars[,"Z"])
      rDDM(data$R,pars, attr(pars, "ok"))
    },
    # Flow joint density (PDF) of (rt, R)
    dfun=function(rt,R,pars) dDDMnn(rt,R,pars),
    # Flow joint (defective) probability (CDF)
    pfun=function(rt,R,pars) pDDMnn(rt,R,pars),
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10)){
      log_likelihood_ddm(pars=pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}
