# The neural-likelihood validation kit (R/nn_validate.R): cells, the
# training-region check, total mass, score bias, posterior shift, SBC cells.
# The shipped RDM flow keeps these fast; the controls (analytic RDM) must
# return about zero, and a control that is deliberately not the truth must
# show up with the right sign.

nle <- function(f) getFromNamespace(f, "EMC2")

rdm_mean <- c(v = log(1.5), B = log(1), A = log(.3), t0 = log(.2))
rdm_cell <- function(...) nn_cell(RDMnn, mean = rdm_mean, sd = .1, ...)

# RDM whose rates are 10 % higher than the parameters say: same parameters,
# defaults and transforms (a valid control), but not the network's model
rdm_fast <- function() {
  m <- RDM()
  T0 <- m$Ttransform
  m$c_name <- NULL
  m$Ttransform <- function(pars, dadm) {
    ok <- attr(pars, "ok")
    pars[, "v"] <- 1.1 * pars[, "v"]
    out <- T0(pars, dadm)
    attr(out, "ok") <- ok
    out
  }
  m
}

# --- cells and the training region ---------------------------------------------------
test_that("nn_cell builds matched network and control cells with named priors", {
  cell <- rdm_cell()
  expect_s3_class(cell, "nn_cell")
  expect_identical(names(cell$mean), c("v_lRa", "v_lRb", "B", "t0", "A"))
  expect_identical(unname(cell$mean[c("v_lRa", "v_lRb")]), rep(log(1.5), 2))
  expect_identical(names(sampled_pars(cell$control_design)), names(cell$mean))
  expect_identical(cell$prior$theta_mu_mean, cell$mean)
  expect_identical(cell$control_prior$theta_mu_mean, cell$mean)
  expect_identical(cell$control_label, "RDM")
  expect_true(all(cell$box$inside))
  # s stays at the model default (it sets the scale)
  expect_false("s" %in% names(cell$mean))
  expect_output(print(cell), "training region")
  # defaults: every sampled network input, strictly inside the region
  d <- nn_cell(DDMnn)
  expect_setequal(names(d$mean), c("v", "a", "t0", "Z", "SZ", "sv", "st0"))
  expect_true(all(d$box$inside))
  # the st0 = 0 artefact: st0 is fixed, so not sampled
  expect_false("st0" %in% names(nn_cell(function() DDMnn("ddm_st0zero"))$mean))
})

test_that("nn_cell and nn_in_box refuse priors that leave the training region", {
  expect_error(nn_cell(RDMnn, mean = rdm_mean, sd = .5), "t0: \\[.*training region")
  expect_error(nn_cell(RDMnn, mean = unname(rdm_mean), sd = .1), "never matched by position")
  expect_error(nn_cell(RDMnn, mean = c(rdm_mean, x = 1), sd = .1), "names x")
  # an effect has no default prior
  expect_error(nn_cell(RDMnn, formula = list(v ~ lR, B ~ 1, t0 ~ 1, A ~ 1), mean = rdm_mean, sd = .1),
               "v_lRb")
  cell <- nn_cell(RDMnn, formula = list(v ~ lR, B ~ 1, t0 ~ 1, A ~ 1),
                  mean = c(rdm_mean, v_lRb = 0), sd = c(v = .1, v_lRb = .05, B = .1, t0 = .1, A = .1))
  expect_true(all(cell$box$inside))
  # the effect enters the second accumulator's range: v + v_lRb +- 4 (sd_v + sd_eff)
  rb <- cell$box[cell$box$parameter == "v" & grepl("lR=b", cell$box$cell), ]
  expect_equal(rb$prior_upper, exp(log(1.5) + 4 * .15), tolerance = 1e-12)

  des <- cell$design
  pri <- prior(des, type = "single", pmean = c(v = log(1.5), v_lRb = 0, B = 0, t0 = log(.2), A = log(.3)),
               psd = c(v = .1, v_lRb = .1, B = .1, t0 = .6, A = .1))
  expect_error(nn_in_box(des, pri), "t0")
  tab <- nn_in_box(pri, refuse = FALSE)                 # the design is taken from the prior
  expect_identical(tab$parameter[!tab$inside], "t0")
  expect_true(all(nn_in_box(des, pri, k = 1)$inside))
  # run_sbc refuses before simulating anything
  expect_error(suppressMessages(run_sbc(des, pri, replicates = 2, trials = 10)), "training region")
})

test_that("the control must take the network's parameters, defaults and transforms", {
  bad_default <- function() { m <- RDM(); m$p_types[["A"]] <- log(.1); m }
  expect_error(rdm_cell(control = bad_default), "defaults differ")
  bad_tf <- function() { m <- RDM(); m$transform$func[["B"]] <- "identity"; m }
  expect_error(rdm_cell(control = bad_tf), "transforms differ")
  expect_error(rdm_cell(control = DDM), "RACE model|DDM model")
  cell <- rdm_cell(control = FALSE)
  expect_null(cell$control)
  expect_error(nn_score_bias(cell, n_draws = 2), "needs an analytic control")
  expect_error(nn_posterior_shift(cell, n_datasets = 2), "needs an analytic control")
})

# --- likelihood plumbing ---------------------------------------------------------------
test_that("R-path trial-wise log-likelihoods equal the compiled ones", {
  cell <- rdm_cell()
  set.seed(4)
  dat <- make_data(cell$mean, cell$design, n_trials = 60)
  for (des in list(cell$design, cell$control_design)) {
    dadm <- nle("nn_dadm")(dat, des)
    P <- nle("nn_cell_draws")(cell, 3)
    native <- nle("nn_trialwise")(P, dadm, des$model)
    ml <- des$model(); ml$c_name <- NULL
    rpath <- nle("nn_trialwise")(P, dadm, function() ml)
    expect_equal(dim(native), c(nrow(dat), 3))
    expect_lt(max(abs(native - rpath)), 1e-8)
    expect_equal(colSums(native), as.numeric(nle("calc_ll_manager")(P, dadm, des$model)), tolerance = 1e-10)
  }
})

# --- total mass --------------------------------------------------------------------------
test_that("total mass: the RDM flow integrates to one like its control", {
  tm <- nn_total_mass(rdm_cell(), n = 4, n_grid = 400)
  s <- tm$summary
  expect_lt(abs(s$mean_log_Z), 1e-3)
  expect_lt(abs(s$mean_log_Z_control), 1e-3)
  expect_lt(s$sd_difference, 1e-3)
  # the flow's leading edge is smeared before t0 (its support is rt > 0), a little
  expect_gt(s$max_before_support, 0)
  expect_lt(s$max_before_support, 1e-3)
  expect_equal(nrow(tm$mass), 4)
  expect_output(print(tm), "Total mass")
})

# --- score bias --------------------------------------------------------------------------
test_that("score bias: the analytic control returns ~0, a wrong truth shows with its sign", {
  sb <- nn_score_bias(rdm_cell(), n_draws = 6, n_check = 3)
  s <- sb$summary
  expect_identical(s$parameter, c("v_lRa", "v_lRb", "B", "A"))    # t0 moves the support
  expect_true(all(abs(s$control) < 1e-3))
  expect_true(all(is.finite(s$g)))
  expect_output(print(sb), "Control ~ 0: yes")
  # data 10 % faster than the network's parameters say: the expected score
  # pulls the rates up, strongly; the control is scored against itself (~0)
  sf <- nn_score_bias(rdm_cell(control = rdm_fast), n_draws = 4, n_check = 2)$summary
  expect_true(all(sf$g[1:2] > 10 * abs(s$g[1:2])))
  expect_true(all(sf$g[1:2] > .03))
  expect_true(all(abs(sf$control) < 1e-3))
  # explicit draws and parameters; the control's error is the central
  # difference's O(h^2): halving h divides it by about four (TRAPS.md #2)
  set.seed(9)
  th <- nle("nn_cell_draws")(rdm_cell(), 3)
  c1 <- nn_score_bias(rdm_cell(), theta = th, pars = c("v_lRa", "B"), h = .02)
  c2 <- nn_score_bias(rdm_cell(), theta = th, pars = c("v_lRa", "B"), h = .01)
  expect_equal(dim(c1$g), c(3, 2))
  ratio <- c1$summary$control / c2$summary$control
  expect_true(all(ratio > 3.5 & ratio < 4.5))
})

# --- posterior shift ---------------------------------------------------------------------
test_that("posterior shift: reproducible, fork-safe, and it sees a wrong likelihood", {
  cell <- rdm_cell()
  a <- nn_posterior_shift(cell, n_datasets = 3, n_trials = 100, n_draws = 300)
  expect_identical(a$summary$parameter, names(cell$mean))
  expect_true(all(is.finite(unlist(a$summary[-1]))))
  expect_true(all(abs(a$summary$shift) < .5))
  expect_true(all(a$ess > 20))
  expect_output(print(a), "Posterior shift")
  skip_on_os("windows")
  b <- nn_posterior_shift(cell, n_datasets = 3, n_trials = 100, n_draws = 300, cores = 2)
  expect_identical(b$summary, a$summary)
  # data from 10 % faster rates: the network's rate posteriors sit well above
  # the control's (and B, which trades off with the rates, below); the paired
  # shift is stable even over few data sets (60 give +0.54, +0.54, -0.51)
  f <- nn_posterior_shift(rdm_cell(control = rdm_fast), n_datasets = 3, n_trials = 200, n_draws = 300)
  expect_true(all(f$summary$shift[1:2] > .3))
  expect_lt(f$summary$shift[3], -.2)
})

# --- SBC cells ---------------------------------------------------------------------------
fake_sbc <- function(n = 200, p = c("a", "b"), skew = FALSE) {
  set.seed(2)
  r <- matrix(runif(n * length(p)), n, dimnames = list(NULL, p))
  if (skew) r[, 2] <- r[, 2]^2
  b <- matrix(rnorm(n * length(p), .1, 1), n, dimnames = list(NULL, p))
  list(rank = list(alpha = r), med = list(alpha = b), bias = list(alpha = b),
       coverage = list(alpha = matrix(runif(n * length(p)) < .95, n, dimnames = list(NULL, p))))
}

test_that("nn_sbc_summary: standardised bias, coverage, KS, envelope", {
  s <- nn_sbc_summary(fake_sbc())
  b <- fake_sbc()$bias$alpha
  expect_equal(s$std_bias, unname(colMeans(b) / apply(b, 2, sd)))
  expect_true(all(s$in_envelope))
  s2 <- nn_sbc_summary(fake_sbc(skew = TRUE))
  expect_identical(s2$in_envelope, c(TRUE, FALSE))
  expect_lt(s2$ks_p[2], 1e-6)
  expect_error(nn_sbc_summary(list(rank = list(mu = 1))), "single-subject")
})

test_that("nn_sbc_cell writes an archive bundle without running, and a README after", {
  dir <- file.path(tempdir(), "nn_sbc_bundle")
  unlink(dir, recursive = TRUE)
  cell <- rdm_cell()
  res <- nn_sbc_cell(cell, trials = 50, replicates = 4, archive_dir = dir, run = FALSE,
                     info = list(why = "test"), cores_per_chain = 2)
  expect_null(res$nn)
  expect_true(all(file.exists(file.path(dir, "bundle", c("cell.rds", "run_cell.R", "session.txt")))))
  expect_s3_class(readRDS(file.path(dir, "bundle", "cell.rds")), "nn_cell")
  scr <- readLines(file.path(dir, "bundle", "run_cell.R"))
  expect_silent(parse(text = scr))
  expect_true(any(grepl("trials = 50, replicates = 4, run_control = TRUE, archive_dir = \"..\", cores_per_chain = 2",
                        scr, fixed = TRUE)))
  rd <- readLines(file.path(dir, "README.md"))
  for (h in c("**Date**", "**Branch / commit**", "**Verdict**", "## Design", "## Priors", "## Conditions",
              "## Fitting", "## Result summary", "## Files", "## Not covered"))
    expect_true(any(grepl(h, rd, fixed = TRUE)), info = h)
  expect_true(any(grepl("| **Why** | test |", rd, fixed = TRUE)))
  # after a run (here: fake SBC results), the verdict and the summary table
  out <- list(nn = fake_sbc(p = names(cell$mean)), control = fake_sbc(p = names(cell$mean)), cell = cell)
  out$summary <- rbind(cbind(model = "network", nn_sbc_summary(out$nn)),
                       cbind(model = "control", nn_sbc_summary(out$control)))
  nle("nn_sbc_write")(dir, out, 50, 4, list(), list())
  rd <- readLines(file.path(dir, "README.md"))
  expect_true(any(grepl("**PASS** (provisional", rd, fixed = TRUE)))
  expect_true(all(file.exists(file.path(dir, "results", c("summary.csv", "ecdf_nn.pdf", "ecdf_control.pdf")))))
  expect_error(nn_sbc_cell(cell, run = FALSE), "archive_dir")
  unlink(dir, recursive = TRUE)
})

# --- fork safety (macOS: Accelerate's threaded BLAS crashes in forked children) --------
test_that("per-trial network evaluations work in forked workers", {
  skip_on_os("windows")
  if (identical(Sys.info()[["sysname"]], "Darwin"))
    expect_true(nzchar(Sys.getenv("VECLIB_MAXIMUM_THREADS")))
  m <- register_nn_model(test_path("golden", "lan_ddm_uniform_st.rds"),
                         p_types = c(v = 0, a = 1, z = .5, t = .5, st = .05))()
  n <- 5000
  set.seed(1)
  th <- cbind(v = runif(n, -1, 1), a = 1, z = .5, t = .5, st = .05)
  rt <- runif(n, .6, 2); R <- sample(1:2, n, TRUE)
  here <- sum(log(m$dfun(rt, R, th)))
  there <- parallel::mclapply(1:2, function(i) sum(log(m$dfun(rt, R, th))), mc.cores = 2)
  expect_identical(unlist(there), rep(here, 2))
})
