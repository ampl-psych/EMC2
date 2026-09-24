# Neural-likelihood (flow) integration: golden vectors, trial-wise evaluation,
# agreement with the reference R port, defaults trap, and a smoke fit.
# Golden vectors come from the NLE project (Python float64); tolerances follow
# its validate_ddm.R: C++ vs golden 1e-9, C++ vs R port 1e-10.

golden_dir <- test_path("golden")
nle <- function(f) getFromNamespace(f, "EMC2")

# reference R implementation of the port (source of truth: NLE handover)
port <- new.env()
sys.source(test_path("port", "flow_race.R"), envir = port)
sys.source(test_path("port", "flow_ddm.R"), envir = port)

# the shipped .rds members have the layout the port's JSON loader produces
port_fix <- function(m) {
  fix <- function(mlp) { mlp$layers <- lapply(mlp$layers, function(l)
    list(W = as.matrix(l$W), b = as.numeric(l$b))); mlp }
  if (!is.null(m$mlp)) m$mlp <- fix(m$mlp)
  if (!is.null(m$flow_mlp)) { m$flow_mlp <- fix(m$flow_mlp); m$classifier_mlp <- fix(m$classifier_mlp) }
  m
}
artefact_member <- function(name) {
  x <- readRDS(file.path(nle("nle_dir")(), paste0(name, ".rds")))
  port_fix(if (!is.null(x$members)) x$members[[1]] else x)
}

test_that("artefacts are pinned and report what was loaded", {
  man <- nle("nle_manifest")()
  expect_setequal(names(man), c("ddm_cap256w_c4.rds", "ddm_st0zero.rds", "rdm_small.rds"))
  m <- nle("nle_meta")("ddm_cap256w_c4")
  expect_identical(m$context_names, c("v", "a", "t0", "s", "Z", "SZ", "sv", "st0"))
  expect_equal(as.numeric(m$checkpoint_step[["flow"]]), 72)
  expect_equal(nle("nle_meta")("rdm_small")$checkpoint_step, 195)
  expect_output(nle("nle_artefact_info")("rdm_small"), "sha256=")
  expect_error(nle("nle_meta")("no_such_artefact"), "No neural-likelihood artefact")
})

# --- golden vectors ---------------------------------------------------------
check_ddm_golden <- function(name, tag, joint = TRUE) {
  golden <- read.csv(file.path(golden_dir, paste0(tag, "_golden.csv")))
  meta <- nle("nle_meta")(name)
  ptr <- nle("nle_get")(name)
  ctx <- meta$context_names
  d_lp <- d_cdf <- d_flow <- d_flowcdf <- 0
  for (pid in unique(golden$param_id)) {
    g <- golden[golden$param_id == pid, ]
    Theta <- matrix(as.numeric(g[1L, ctx]), nrow(g), length(ctx), byrow = TRUE)
    ev <- nle("ddm_ens_eval_trials_cpp")(ptr, Theta, g$rt, as.integer(g$R))
    # flow part only: remove each side's own log P(R)
    d_flow <- max(d_flow, abs((ev$log_pdf - log(ev$p_R)) - (g$log_pdf_joint - g$log_p_R)))
    d_flowcdf <- max(d_flowcdf, abs(ev$cdf / ev$p_R - g$cdf_defective / exp(g$log_p_R)))
    d_lp <- max(d_lp, abs(ev$log_pdf - g$log_pdf_joint))
    d_cdf <- max(d_cdf, abs(ev$cdf - g$cdf_defective))
  }
  expect_lt(d_flow, 1e-9)
  expect_lt(d_flowcdf, 1e-9)
  if (joint) { expect_lt(d_lp, 1e-9); expect_lt(d_cdf, 1e-9) }
  invisible(c(flow = d_flow, joint = d_lp))
}

test_that("DDM cap256w_c4: flow matches golden (classifier checked separately)", {
  # The golden's log_p_R comes from classifier step 135, not c4's "tuned-big"
  # classifier, so only the flow columns are compared; the classifier is
  # checked for normalisation below.
  check_ddm_golden("ddm_cap256w_c4", "ddm_cap256w", joint = FALSE)
  meta <- nle("nle_meta")("ddm_cap256w_c4")
  th <- matrix(c(1, log(1), log(.3), 0, 0, qnorm(.3), log(.5), log(.1)), 1)
  p <- vapply(1:2, function(r)
    nle("ddm_ens_eval_trials_cpp")(nle("nle_get")("ddm_cap256w_c4"), th, 1, r)$p_R, 0)
  expect_equal(sum(p), 1, tolerance = 1e-12)
})

test_that("DDM st0zero (z0_s503_e282): joint density matches golden", {
  check_ddm_golden("ddm_st0zero", "ddm_z0_s503_e282")
})

test_that("RDM rdm_small (e195): pdf, cdf and survivor match golden", {
  golden <- read.csv(file.path(golden_dir, "rdm_rdm_small_e195_golden.csv"))
  meta <- nle("nle_meta")("rdm_small")
  ptr <- nle("nle_get")("rdm_small")
  d <- c(lp = 0, cdf = 0, sf = 0)
  for (pid in unique(golden$param_id)) {
    g <- golden[golden$param_id == pid, ]
    Theta <- matrix(as.numeric(g[1L, meta$context_names]), nrow(g), 5, byrow = TRUE)
    ev <- nle("flow_eval_trials_cpp")(ptr, Theta, g$rt)
    d["lp"] <- max(d["lp"], abs(ev$log_pdf - g$log_pdf))
    d["cdf"] <- max(d["cdf"], abs(ev$cdf - g$cdf))
    d["sf"] <- max(d["sf"], abs(ev$log_sf - g$log_sf))
  }
  expect_lt(max(d), 1e-9)
})

# --- trial-wise evaluation --------------------------------------------------
test_that("trial-wise evaluation: caching and out-of-box rows", {
  meta <- nle("nle_meta")("ddm_cap256w_c4")
  ptr <- nle("nle_get")("ddm_cap256w_c4")
  golden <- read.csv(file.path(golden_dir, "ddm_cap256w_golden.csv"))
  ctx <- meta$context_names
  g1 <- golden[golden$param_id == 1, ]; g2 <- golden[golden$param_id == 2, ]
  mid <- (meta$bounds_natural$lower + meta$bounds_natural$upper) / 2
  th_bad <- as.numeric(g1[1L, ctx]); th_bad[1] <- 7          # v above the box
  rows <- function(g) matrix(as.numeric(g[1L, ctx]), nrow(g), length(ctx), byrow = TRUE)
  Theta <- rbind(rows(g1), matrix(th_bad, 3, length(ctx), byrow = TRUE), rows(g2))
  rt <- c(g1$rt, .3, 1, 3, g2$rt); R <- as.integer(c(g1$R, 1, 2, 1, g2$R))
  tv <- nle("ddm_ens_eval_trials_cpp")(ptr, Theta, rt, R)
  n1 <- nrow(g1); bad <- n1 + 1:3
  expect_true(all(tv$log_pdf[bad] == -Inf))
  expect_true(all(tv$cdf[bad] == 0))
  # per-set evaluation is the same as the interleaved, cached one
  one <- nle("ddm_ens_eval_trials_cpp")(ptr, rows(g2), g2$rt, as.integer(g2$R))
  expect_equal(tv$log_pdf[n1 + 3 + seq_len(nrow(g2))], one$log_pdf, tolerance = 1e-12)
})

# --- agreement with the reference R port (protects the evaluator rewrite) ----
test_that("C++ evaluator agrees with the R port over 40 random thetas (1e-10)", {
  set.seed(42)
  check <- function(name, ddm) {
    fl <- artefact_member(name)
    ptr <- nle("nle_get")(name)
    meta <- nle("nle_meta")(name)
    lo <- fl$bounds_sampled$lower; hi <- fl$bounds_sampled$upper
    ctx <- fl$context_names
    if (!ddm) { lo <- lo; hi <- hi }
    if (ddm) {                      # bounds_sampled has one entry per context
      lo <- lo[seq_along(ctx)]; hi <- hi[seq_along(ctx)]
    }
    worst <- 0
    for (i in 1:40) {
      # interior points only (avoid the exact box faces)
      theta <- lo + (hi - lo) * runif(length(lo), .02, .98)
      if (ddm && name == "ddm_cap256w_c4") theta[6] <- theta[6]  # SZ, Z stay in-box by construction
      rt <- exp(runif(25, log(.06), log(3)))
      if (ddm) {
        R <- sample(1:2, 25, TRUE)
        ref <- port$ddm_eval(fl, theta, rt, R)
        cv <- nle("ddm_ens_eval_trials_cpp")(ptr, matrix(theta, 25, length(theta), byrow = TRUE), rt, R)
        worst <- max(worst, abs(cv$log_pdf - ref$log_pdf), abs(cv$cdf - ref$cdf))
      } else {
        ref <- port$flow_eval(fl, theta, rt)
        cv <- nle("flow_eval_trials_cpp")(ptr, matrix(theta, 25, length(theta), byrow = TRUE), rt)
        worst <- max(worst, abs(cv$log_pdf - ref$log_pdf), abs(cv$cdf - ref$cdf))
      }
    }
    worst
  }
  expect_lt(check("ddm_cap256w_c4", TRUE), 1e-10)
  expect_lt(check("ddm_st0zero", TRUE), 1e-10)
  expect_lt(check("rdm_small", FALSE), 1e-10)
})

# --- model definitions ------------------------------------------------------
test_that("DDMnn / RDMnn defaults are identical to DDM / RDM (defaults trap)", {
  expect_identical(DDMnn()$p_types, DDM()$p_types)
  expect_identical(RDMnn()$p_types, RDM()$p_types)
  expect_identical(DDMnn()$transform, DDM()$transform)
  expect_equal(DDMnn("ddm_st0zero")$bound$exception, c(st0 = 0))
  expect_length(DDMnn()$bound$exception, 0)
})

test_that("an un-sampled parameter the artefact cannot represent is refused", {
  pars <- cbind(v = 1, a = 1, sv = 0, t0 = .2, st0 = .1, s = 1, Z = .5, SZ = .3)
  m <- DDMnn()
  expect_error(m$Ttransform(pars, NULL), "Parameter 'sv' is at the model default")
  pars[, "sv"] <- .5
  expect_silent(m$Ttransform(pars, NULL))
  # st0 = 0 is supported by the st0-free artefact, refused by the full one
  pars[, "st0"] <- 0
  expect_error(m$Ttransform(pars, NULL), "'st0'")
  expect_silent(DDMnn("ddm_st0zero")$Ttransform(pars, NULL))
  # a missing context column is an error, never a silent reorder
  expect_error(nle("nle_theta")(pars[, c("v", "a")], nle("nle_meta")("ddm_cap256w_c4")),
               "needs parameter")
})

# --- likelihood vs the analytic twin ---------------------------------------
test_that("RDMnn single-accumulator density and CDF are close to the analytic RDM", {
  rt <- seq(.4, 2.5, length.out = 60)
  P <- cbind(v = 1.5, B = .8, A = .3, t0 = .25, s = 1, b = 1.1)[rep(1, 60), ]
  m <- RDMnn()
  expect_lt(max(abs(m$dfun(rt, P) - RDM()$dfun(rt, P))), 5e-2)
  expect_lt(max(abs(m$pfun(rt, P) - RDM()$pfun(rt, P))), 5e-3)
})

# --- smoke fit --------------------------------------------------------------
test_that("DDMnn runs through make_emc() and a short fit", {
  skip_on_os("windows")
  skip_on_cran()
  set.seed(7)
  des <- design(factors = list(subjects = 1:2, S = c("a", "b")),
                Rlevels = c("a", "b"), model = DDMnn,
                formula = list(v ~ S, a ~ 1, t0 ~ 1),
                constants = c(s = log(1), sv = log(.5), SZ = qnorm(.3),
                              st0 = log(.1), Z = qnorm(.5)),
                report_p_vector = FALSE)
  p <- sampled_pars(des, doMap = FALSE)
  p[] <- c(0.5, -.5, log(1.2), log(.25))[seq_along(p)]
  dat <- make_data(p, des, n_trials = 40)
  expect_true(all(is.finite(dat$rt)))
  emc <- make_emc(dat, des, n_chains = 2, compress = FALSE, rt_resolution = 0.001)
  emc <- fit(emc, cores_for_chains = 1, stop_criteria = list(
    preburn = list(iter = 10), burn = list(mean_gd = 5), adapt = list(min_unique = 5),
    sample = list(iter = 10)), verbose = FALSE, particle_factor = 10, step_size = 10)
  expect_s3_class(emc, "emc")
})
