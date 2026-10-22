# Neural-likelihood (flow) integration: golden vectors, trial-wise evaluation,
# agreement with the reference R port, defaults trap, and a smoke fit.
# Golden vectors come from the NLE project (Python float64); tolerances follow
# its validate_ddm.R: C++ vs golden 1e-9, C++ vs R port 1e-10.

nle <- function(f) getFromNamespace(f, "EMC2")
rdm_name <- "rdm_ss_base_s801_e271_weights.json"
ddm_name <- "ddm_sc_s1106_e245_weights.json"
rtcon_name <- "rtcon2_r3_s2_e006_weights.json"

# reference R implementation of the port (source of truth: NLE handover)
port <- new.env()
sys.source(test_path("port", "flow_race.R"), envir = port)
sys.source(test_path("port", "flow_ddm.R"), envir = port)

test_that("artefacts are pinned and report what was loaded", {
  man <- nle("nle_manifest")()
  expect_setequal(names(man), c(ddm_name, rdm_name, rtcon_name))
  reg_rdm <- register_nn_model(rdm_name)()$nn
  expect_identical(reg_rdm$kind, "flow_race")
  expect_identical(reg_rdm$context_names, c("v", "B", "t0", "s", "A"))
  reg_ddm <- register_nn_model(ddm_name)()$nn
  expect_identical(reg_ddm$kind, "flow_scaled")
  expect_output(nle("nle_artefact_info")(rdm_name), "(sha256|md5)=")
  expect_error(register_nn_model("no_such_artefact"), "No neural-likelihood artefact")
})

# --- K-category flow (flow_cat): RTCON2 ---------------------------------------
test_that("the RTCON2 card registers as a flow_cat model", {
  # no analytic twin in EMC2: a bare registration says so; RTCONnn supplies the parameterisation
  expect_error(register_nn_model(rtcon_name), "no analytic twin")
  reg <- register_nn_model(rtcon_name, twin = nle("rtcon2_twin"))()$nn
  expect_identical(reg$kind, "flow_cat")
  ml <- RTCONnn()
  expect_identical(ml$nn$kind, "flow_cat")
  expect_identical(ml$nn$n_categories, 6L)
  expect_identical(unname(ml$nn$conv), c("t0", "st0"))
  expect_identical(ml$nn$pars, c("mu", "s", "C1", "g2", "g3", "g4", "g5", "a", "sb", paste0("B", 1:6)))
  expect_identical(ml$extra_responses, "RR")
  expect_identical(ml$c_name, "NN")
  expect_equal(ml$bound$exception[c("sb", "st0")], c(sb = 0, st0 = 0))
  expect_null(ml$nn$card$golden)       # the card's embedded golden vectors are not carried by designs
})

test_that("flow_cat: compiled likelihood = R path = multithreaded, with censoring", {
  set.seed(1)
  n <- 200
  dat <- data.frame(subjects = factor(1), S = factor(sample(c("new", "old"), n, TRUE)),
                    R = factor(sample(c("new", "old"), n, TRUE)), RR = sample(1:3, n, TRUE),
                    rt = 0.45 + rexp(n, 2))
  fl <- c(list(mu ~ S), lapply(c("s", "C1", "g2", "g3", "g4", "g5", "a", "sb", paste0("B", 1:6), "t0", "st0"),
                               function(p) stats::as.formula(paste(p, "~ 1"))))
  utils::capture.output(des <- suppressMessages(design(data = dat, model = RTCONnn, formula = fl)))
  pv <- c(mu = 0, mu_Sold = 1.2, s = 0, C1 = -1, g2 = .8, g3 = .8, g4 = .8, g5 = .8, a = log(.04), sb = log(.8),
          B1 = log(2), B2 = log(2), B3 = log(2), B4 = log(2), B5 = log(2), B6 = log(2), t0 = log(.38))
  for (st0 in c(log(.1), -Inf)) for (cens in c(FALSE, TRUE)) {
    d <- dat
    if (cens) { d <- make_missing(dat, LC = 0.5, UC = 1.6); d$RR[!is.na(d$missingness)] <- NA }
    emc <- suppressMessages(make_emc(d, des, type = "single", n_chains = 1, compress = TRUE,
                                     rt_resolution = .001, verbose = FALSE))
    dadm <- emc[[1]]$data[[1]]; model <- emc[[1]]$model
    P <- matrix(c(pv, st0 = st0), 1, dimnames = list(NULL, c(names(pv), "st0")))
    native <- nle("calc_ll_manager")(P, dadm, model)
    mr <- model(); mr$c_name <- NULL
    rpath <- nle("calc_ll_manager")(P, dadm, function() mr)
    o <- options(emc.ll_backend = "multithreaded", emc.n_threads = 2)
    mt <- nle("calc_ll_manager")(P, dadm, model)
    options(o)
    expect_true(is.finite(native))
    expect_lt(abs(native - rpath), 1e-10)
    expect_lt(abs(native - mt), 1e-10)
  }
  # out of the training box (s below 0.2): every trial at the floor
  P[, "s"] <- log(.1)
  expect_equal(as.numeric(nle("calc_ll_manager")(P, dadm, model)), n * log(1e-10))
})

test_that("flow_cat: data checks", {
  dat <- data.frame(subjects = factor(1), R = factor(c("new", "old", "old")), RR = c(1, 2, 3), rt = c(.6, .7, .8))
  utils::capture.output(des <- suppressMessages(design(data = dat, model = RTCONnn, formula = list(mu ~ 1))))
  bad <- dat; bad$RR[2] <- 4
  expect_error(make_emc(bad, des, type = "single", rt_resolution = .001, verbose = FALSE), "RR must be an integer in 1..3")
  bad <- dat; bad$RR <- NULL
  expect_error(make_emc(bad, des, type = "single", rt_resolution = .001, verbose = FALSE), "needs a rating column RR")
  bad <- make_missing(dat, UC = .75); bad$missingness[3] <- 3L
  expect_error(make_emc(bad, des, type = "single", rt_resolution = .001, verbose = FALSE), "only missingness codes 1")
  expect_error(design(data = dat, model = RTCONnn, formula = list(mu ~ RR)), "cannot be used as predictors")
})

# --- batched MLP forward pass (no artefact needed) --------------------------
test_that("batched MLP forward pass matches R matrix algebra, with and without LayerNorm", {
  set.seed(3)
  mk <- function(dims, norm) {
    layers <- lapply(seq_len(length(dims) - 1), function(i)
      list(W = matrix(rnorm(dims[i] * dims[i + 1], sd = 1 / sqrt(dims[i])), dims[i]),
           b = rnorm(dims[i + 1], sd = .1)))
    mlp <- list(activation = "gelu_tanh", layers = layers, use_norm = norm)
    if (norm) mlp$norms <- lapply(dims[2:(length(dims) - 1)], function(d)
      list(scale = runif(d, .5, 1.5), bias = rnorm(d, sd = .1), eps = 1e-5))
    mlp
  }
  for (norm in c(FALSE, TRUE)) {
    mlp <- mk(c(6, 40, 30, 11), norm)
    X <- matrix(rnorm(300 * 6, sd = 3), 300)
    ref <- port$.mlp_fwd(mlp, X, use_norm = norm, norms = mlp$norms)
    expect_lt(max(abs(nle("nle_mlp_forward")(mlp, X) - ref)), 1e-12)
    expect_lt(max(abs(nle("nle_mlp_forward")(mlp, X[1, , drop = FALSE]) - ref[1, ])), 1e-12)
  }
})

# --- RDM: distinct-row batching and per-trial agreement with R port ----------
test_that("RDM: distinct-row batching equals row-by-row evaluation (interleaved, out-of-box)", {
  set.seed(11)
  rdm_path <- system.file("extdata", "flownn", rdm_name, package = "EMC2")
  fl <- port$flow_load(rdm_path)
  ptr <- nle("nn_ptr")(register_nn_model(rdm_name)()$nn)
  lo <- fl$bounds_sampled$lower; hi <- fl$bounds_sampled$upper
  box_rows <- function(n) t(replicate(n, lo + (hi - lo) * runif(length(lo), .02, .98)))
  D <- box_rows(40)
  D[1:3, 1] <- hi[1] + 1   # out of box
  Th <- D[sample(40, 700, TRUE), ]
  rt <- exp(runif(700, log(.06), log(3)))
  tv <- nle("flow_eval_trials_cpp")(ptr, Th, rt)
  rw <- lapply(1:700, function(t) nle("flow_eval_cpp")(ptr, Th[t, ], rt[t]))
  for (o in c("pdf", "cdf", "log_pdf", "log_sf"))
    expect_equal(tv[[o]], vapply(rw, `[[`, 0, o), tolerance = 1e-10, label = o)
  expect_true(all(tv$log_pdf[Th[, 1] > hi[1]] == -Inf))
})

test_that("RDM: per-trial rows agree with the R port (1e-10)", {
  set.seed(5)
  rdm_path <- system.file("extdata", "flownn", rdm_name, package = "EMC2")
  fl <- port$flow_load(rdm_path)
  ptr <- nle("nn_ptr")(register_nn_model(rdm_name)()$nn)
  lo <- fl$bounds_sampled$lower; hi <- fl$bounds_sampled$upper
  n <- 300
  Th <- t(replicate(n, lo + (hi - lo) * runif(length(lo), .02, .98)))
  rt <- exp(runif(n, log(.06), log(3)))
  cv <- nle("flow_eval_trials_cpp")(ptr, Th, rt)
  ref <- lapply(1:n, function(t) port$flow_eval(fl, Th[t, ], rt[t]))
  expect_lt(max(abs(cv$log_pdf - vapply(ref, `[[`, 0, "log_pdf")),
                abs(cv$cdf   - vapply(ref, `[[`, 0, "cdf"))), 1e-10)
})

# --- model definitions -------------------------------------------------------
test_that("DDMnn / RDMnn defaults are identical to DDM / RDM (defaults trap)", {
  expect_identical(DDMnn()$p_types, DDM()$p_types)
  expect_identical(RDMnn()$p_types, RDM()$p_types)
  expect_identical(DDMnn()$transform, DDM()$transform)
  expect_equal(DDMnn()$bound$exception["st0"], c(st0 = 0))
})

test_that("RDMnn rejects t0 at its natural default (0, below floor 0.05); A = 0 is valid via zerobox", {
  m_rdm <- RDMnn()
  pars_bad <- cbind(v = 1, B = 1, A = 0, t0 = 0,  s = 1)   # t0 = 0 below floor
  pars_ok  <- cbind(v = 1, B = 1, A = 0, t0 = .2, s = 1)   # A = 0 valid (zerobox)
  expect_error(m_rdm$Ttransform(pars_bad, NULL), "Parameter 't0' is at the model default")
  expect_silent(m_rdm$Ttransform(pars_ok, NULL))
  # DDM scaled: sv = 0 and SZ = 0 are exact zeros (OK); st0 fixed by card
  pars_ddm <- cbind(v = 1, a = 1, sv = 0, SZ = 0, t0 = .2, s = 1, Z = .5)
  m_ddm <- DDMnn()
  expect_silent(m_ddm$Ttransform(pars_ddm, NULL))
})

# --- likelihood vs the analytic twin ----------------------------------------
test_that("RDMnn single-accumulator density and CDF are close to the analytic RDM", {
  rt <- seq(.4, 2.5, length.out = 60)
  P <- cbind(v = 1.5, B = .8, A = .3, t0 = .25, s = 1, b = 1.1)[rep(1, 60), ]
  m <- RDMnn()
  expect_lt(max(abs(m$dfun(rt, P) - RDM()$dfun(rt, P))), 5e-2)
  expect_lt(max(abs(m$pfun(rt, P) - RDM()$pfun(rt, P))), 5e-3)
})

# --- smoke fit ---------------------------------------------------------------
test_that("DDMnn runs through make_emc() and a short fit", {
  skip_on_os("windows")
  skip_on_cran()
  set.seed(7)
  des <- design(factors = list(subjects = 1:2, S = c("a", "b")),
                Rlevels = c("a", "b"), model = DDMnn,
                formula = list(v ~ S, a ~ 1, t0 ~ 1),
                constants = c(s = log(1), sv = -Inf, SZ = -Inf, Z = qnorm(.5)),
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
