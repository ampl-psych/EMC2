# The native likelihood branch for neural likelihoods (calc_ll type "NN",
# src/model_NN.h) against the R path (dfun/pfun via log_likelihood_ddm/race),
# which stays the reference. Sums agree to 1e-9 relative (the R path sums in
# long double; the race survivor is log(1 - F) computed directly natively).

nle <- function(f) getFromNamespace(f, "EMC2")
ll_manager <- function(...) nle("calc_ll_manager")(...)

# The same model on the R path
r_path <- function(model) {
  ml <- model()
  ml$c_name <- NULL
  function() ml
}

setup_nn <- function(des, p_vector, n_trials, compress = TRUE, seed = 1) {
  set.seed(seed)
  dat <- make_data(p_vector, des, n_trials = n_trials)
  emc <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 1, compress = compress,
                                   verbose = FALSE, rt_resolution = .001))
  list(dadm = emc[[1]]$data[[1]], model = emc[[1]]$model, dat = dat)
}

proposals_around <- function(p_vector, n = 40, sd = .1, n_out = 3) {
  P <- matrix(p_vector, n, length(p_vector), byrow = TRUE, dimnames = list(NULL, names(p_vector)))
  P <- P + matrix(rnorm(length(P), sd = sd), n)
  P[seq_len(n_out), 1] <- P[seq_len(n_out), 1] + 20        # out of bounds
  P
}

rel_diff <- function(a, b) max(abs(a - b) / pmax(1, abs(b)))

# Trial-wise log-likelihoods of one particle on the R path, from the model's
# dfun/pfun as log_likelihood_ddm/race combine them (the R path's
# calc_ll_manager returns sums only)
r_trialwise <- function(p, dadm, model, min_ll = log(1e-10)) {
  m <- model()
  pars <- nle("get_pars_matrix_oo")(p, dadm, model)
  ok <- attr(pars, "ok")
  if (m$type == "DDM") {
    like <- numeric(nrow(dadm))
    like[ok] <- m$dfun(dadm$rt[ok], dadm$R[ok], pars[ok, , drop = FALSE])
    return(pmax(min_ll, log(like))[attr(dadm, "expand")])
  }
  if (!is.null(dadm$RACE)) pars[as.numeric(dadm$lR) > as.numeric(as.character(dadm$RACE)), ] <- NA
  w <- dadm$winner
  lds <- numeric(nrow(dadm))
  lds[w] <- log(m$dfun(dadm$rt[w], pars[w, , drop = FALSE]))
  lds[!w] <- log(1 - m$pfun(dadm$rt[!w], pars[!w, , drop = FALSE]))
  lds[is.na(lds) | !ok] <- min_ll
  n_acc <- length(levels(dadm$R))
  ll <- lds[w] + colSums(matrix(lds[!w], nrow = n_acc - 1))
  pmax(min_ll, ll)[attr(dadm, "expand")]
}

# native vs R path: sums, trial-wise, and the multithreaded backend
expect_native_matches_r <- function(P, dadm, model) {
  expect_identical(model()$c_name, "NN")
  llN <- ll_manager(P, dadm, model)
  llR <- ll_manager(P, dadm, r_path(model))
  expect_lt(rel_diff(llN, llR), 1e-9)
  twN <- ll_manager(P, dadm, model, return_trialwise = TRUE)
  expect_equal(dim(twN), c(length(attr(dadm, "expand")), nrow(P)))
  expect_lt(rel_diff(colSums(twN), as.numeric(llN)), 1e-12)
  for (i in c(1, nrow(P) - 1, nrow(P)))     # an out-of-bounds particle and two others
    expect_lt(max(abs(twN[, i] - r_trialwise(P[i, ], dadm, model))), 1e-9)
  old <- options(emc.ll_backend = "multithreaded", emc.n_threads = 2)
  on.exit(options(old))
  expect_identical(ll_manager(P, dadm, model), llN)
  expect_identical(ll_manager(P, dadm, model, return_trialwise = TRUE), twN)
  invisible(llN)
}

ddm_design <- function(model = DDMnn, formula = list(v ~ S, a ~ 1, t0 ~ 1, Z ~ 1, SZ ~ 1, sv ~ 1, st0 ~ 1),
                       constants = c(s = 0))
  suppressMessages(design(factors = list(subjects = 1, S = c("a", "b")), Rlevels = c("a", "b"),
                          model = model, formula = formula, constants = constants,
                          report_p_vector = FALSE))

ddm_p <- function(des) {
  p <- sampled_pars(des, doMap = FALSE)
  vals <- c(v = 1, v_Sb = -.5, a = log(1.2), t0 = log(.25), Z = 0, SZ = qnorm(.3),
            sv = log(.5), st0 = log(.1))
  p[] <- vals[names(p)]
  p
}

race_design <- function(levels = c("a", "b"), model = RDMnn)
  suppressMessages(design(factors = list(subjects = 1, S = levels), Rlevels = levels,
                          model = model, formula = list(v ~ lR, B ~ 1, t0 ~ 1, A ~ 1),
                          constants = c(s = log(1)), report_p_vector = FALSE))

race_p <- function(des) {
  p <- sampled_pars(des, doMap = FALSE)
  p[] <- log(1)
  p[grepl("^v_lR", names(p))] <- log(1.8)
  p[c("B", "t0", "A")] <- log(c(1, .25, .3))
  p
}

# --- a. agreement with the R path ---------------------------------------------
test_that("DDMnn (cap256w_c4): native equals the R path, all parameters sampled", {
  des <- ddm_design()
  p <- ddm_p(des)
  s <- setup_nn(des, p, 150)
  expect_native_matches_r(proposals_around(p), s$dadm, s$model)
})

test_that("DDMnn (ddm_st0zero, st0 fixed at 0): native equals the R path", {
  des <- ddm_design(model = function() DDMnn("ddm_st0zero"),
                    formula = list(v ~ S, a ~ 1, t0 ~ 1, Z ~ 1, SZ ~ 1, sv ~ 1))
  p <- ddm_p(des)
  s <- setup_nn(des, p, 120)
  expect_native_matches_r(proposals_around(p), s$dadm, s$model)
})

test_that("RDMnn: native equals the R path with two and three accumulators", {
  for (lev in list(c("a", "b"), c("a", "b", "c"))) {
    des <- race_design(lev)
    p <- race_p(des)
    s <- setup_nn(des, p, 120)
    expect_native_matches_r(proposals_around(p), s$dadm, s$model)
  }
})

test_that("pre = 't0': times at or below t0 get zero density / survivor one, as on the R path", {
  m <- register_nn_model("rdm_small", pre = "t0")
  des <- race_design(model = m)
  p <- race_p(des)
  s <- setup_nn(des, p, 120)
  P <- proposals_around(p)
  P[, "t0"] <- log(runif(nrow(P), .2, .45))   # past the fastest rts
  expect_native_matches_r(P, s$dadm, s$model)
  rt0 <- min(s$dadm$rt)
  expect_true(any(exp(P[, "t0"]) > rt0))
})

test_that("absent accumulators (a RACE column) contribute survivor one, as on the R path", {
  des <- race_design(c("a", "b", "c"))
  p <- race_p(des)
  s <- setup_nn(des, p, 150, compress = FALSE)
  dadm <- s$dadm
  n_acc <- 3
  win_lR <- as.integer(dadm$lR)[dadm$winner]
  race <- ifelse(win_lR <= 2 & runif(length(win_lR)) < .5, 2, 3)
  dadm$RACE <- factor(rep(race, each = n_acc), levels = c(2, 3))
  P <- proposals_around(p)
  llN <- expect_native_matches_r(P, dadm, s$model)
  expect_false(isTRUE(all.equal(llN, ll_manager(P, s$dadm, s$model))))
})

# --- b. which models run natively ----------------------------------------------
test_that("a pre function takes the R path; flows otherwise run natively", {
  m_fun <- register_nn_model("rdm_small", pre = function(rt, pars) rt - pars[, "t0"])()
  expect_null(m_fun$c_name)
  expect_null(m_fun$nn$native)
  expect_identical(register_nn_model("rdm_small", pre = "t0")()$c_name, "NN")
  expect_identical(DDMnn()$c_name, "NN")
  expect_identical(RDMnn()$c_name, "NN")

  m <- register_nn_model("rdm_small", pre = function(rt, pars) rt - pars[, "t0"])
  des <- race_design(model = m)
  p <- race_p(des)
  s <- setup_nn(des, p, 60)
  expect_true(is.finite(ll_manager(matrix(p, 1, dimnames = list(NULL, names(p))), s$dadm, s$model)))
})

# --- c. refusals -------------------------------------------------------------------
test_that("censored or truncated data are refused on both paths", {
  des <- ddm_design()
  p <- ddm_p(des)
  s <- setup_nn(des, p, 40, compress = FALSE)
  P <- matrix(p, 1, dimnames = list(NULL, names(p)))
  d1 <- s$dadm
  d1$missingness <- NA_integer_
  expect_true(is.finite(ll_manager(P, d1, s$model)))
  d1$missingness[3] <- 2L
  expect_error(ll_manager(P, d1, s$model), "censored")
  expect_error(ll_manager(P, d1, r_path(s$model)), "censored")
  d2 <- s$dadm
  d2$LT <- 0; d2$UT <- Inf
  expect_true(is.finite(ll_manager(P, d2, s$model)))
  d2$LT <- .1
  expect_error(ll_manager(P, d2, s$model), "truncated")
  expect_error(ll_manager(P, d2, r_path(s$model)), "truncated")
  d2$LT <- 0; d2$UT <- 3
  expect_error(ll_manager(P, d2, s$model), "truncated")
  expect_error(ll_manager(P, d2, r_path(s$model)), "truncated")
})

test_that("the run-time refusals (fixed parameter, default outside the box) match the R path", {
  # a fixed parameter set to another value behind design()'s back
  des <- ddm_design(model = function() DDMnn("ddm_st0zero"), formula = list(v ~ S, a ~ 1, t0 ~ 1),
                    constants = c(s = 0, sv = log(.5), SZ = qnorm(.3), Z = 0))
  p <- ddm_p(des)
  s <- setup_nn(des, p, 30)
  P <- matrix(p, 1, dimnames = list(NULL, names(p)))
  expect_true(is.finite(ll_manager(P, s$dadm, s$model)))
  d <- s$dadm
  attr(d, "constants")[["st0"]] <- log(.1)
  expect_error(ll_manager(P, d, s$model), "fixed at 0")
  expect_error(ll_manager(P, d, r_path(s$model)), "fixed at 0")
  old <- options(emc.ll_backend = "multithreaded", emc.n_threads = 2)
  on.exit(options(old))
  expect_error(ll_manager(P, d, s$model), "fixed at 0")
  options(old)

  # an un-sampled parameter back at the twin's default (sv = 0, outside the box)
  des2 <- ddm_design(formula = list(v ~ S, a ~ 1, t0 ~ 1),
                     constants = c(s = 0, sv = log(.5), SZ = qnorm(.3), Z = 0, st0 = log(.1)))
  s2 <- setup_nn(des2, ddm_p(des2), 30)
  P2 <- matrix(ddm_p(des2), 1, dimnames = list(NULL, names(ddm_p(des2))))
  d2 <- s2$dadm
  attr(d2, "constants")[["sv"]] <- -Inf
  expect_error(ll_manager(P2, d2, s2$model), "'sv' is at the model default")
  expect_error(ll_manager(P2, d2, r_path(s2$model)), "'sv' is at the model default")
})

# --- d. registry: the evaluator is rebuilt for the native path ------------------
test_that("after save/load and a cleared cache the native path rebuilds the evaluator", {
  des <- race_design()
  p <- race_p(des)
  s <- setup_nn(des, p, 40)
  P <- proposals_around(p, n = 5, n_out = 1)
  ll1 <- ll_manager(P, s$dadm, s$model)
  tmp <- tempfile(fileext = ".rds")
  saveRDS(list(dadm = s$dadm, model = s$model), tmp)
  rm(list = ls(nle("nle_cache")), envir = nle("nle_cache"))
  x <- readRDS(tmp)
  expect_identical(ll_manager(P, x$dadm, x$model), ll1)
})

# --- e. a native fit is the R-path fit -------------------------------------------
test_that("a native fit reproduces the R-path fit draw for draw", {
  skip_on_cran()
  des <- ddm_design(formula = list(v ~ S, a ~ 1, t0 ~ 1),
                    constants = c(s = 0, sv = log(.5), SZ = qnorm(.3), Z = 0, st0 = log(.1)))
  p <- ddm_p(des)
  s <- setup_nn(des, p, 60, seed = 11)
  emc <- suppressMessages(make_emc(s$dat, des, type = "single", n_chains = 2,
                                   compress = TRUE, verbose = FALSE, rt_resolution = .001))
  sc <- list(preburn = list(iter = 10), burn = list(mean_gd = 5), adapt = list(min_unique = 5),
             sample = list(iter = 20))
  fit_with <- function(model) {
    e <- emc
    for (k in seq_along(e)) e[[k]]$model <- model
    set.seed(99)
    suppressMessages(fit(e, stop_criteria = sc, cores_for_chains = 1, verbose = FALSE,
                         particle_factor = 10, step_size = 10))
  }
  fN <- fit_with(emc[[1]]$model)
  fR <- fit_with(r_path(emc[[1]]$model))
  aN <- get_pars(fN, stage = "sample", merge_chains = TRUE, return_mcmc = FALSE)
  aR <- get_pars(fR, stage = "sample", merge_chains = TRUE, return_mcmc = FALSE)
  expect_equal(aN, aR, tolerance = 1e-8)
})
