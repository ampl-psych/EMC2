trapz1d <- function(x, y) {
  sum((y[-1] + y[-length(y)]) * diff(x) / 2)
}

repeat_pars <- function(pars, n) {
  pars[rep(1, n), , drop = FALSE]
}

blended_cdm_fpt <- function(tt, a = 1, sigma = 1) {
  v_scale <- (sigma * sigma) / (a * a)
  s <- tt * v_scale
  s0 <- 0.002
  s1 <- 0.02
  fpt_series <- series_bessel_fpt(tt, a = a, sigma = sigma)
  fpt_small <- v_scale * small_t_fpt(s, 1e-8 / (a * a))
  w <- pmin(pmax((s - s0) / (s1 - s0), 0), 1)
  (1 - w) * fpt_small + w * fpt_series
}

blended_sdm_fpt <- function(tt, a = 1, sigma = 1) {
  v_scale <- (sigma * sigma) / (a * a)
  s <- tt * v_scale
  s0 <- 0.002
  s1 <- 0.02
  fpt_series <- series_sdm_fpt(tt, a = a, sigma = sigma)
  fpt_small <- v_scale * small_t_fpt_sdm(s, 1e-8 / (a * a))
  w <- pmin(pmax((s - s0) / (s1 - s0), 0), 1)
  (1 - w) * fpt_small + w * fpt_series
}

blended_hsdm_fpt <- function(tt, a = 1, sigma = 1) {
  v_scale <- (sigma * sigma) / (a * a)
  s <- tt * v_scale
  s0 <- 0.002
  s1 <- 0.02
  fpt_series <- series_hsdm_fpt(tt, a = a, sigma = sigma)
  fpt_small <- v_scale * small_t_fpt_hsdm(s, 1e-8 / (a * a))
  w <- pmin(pmax((s - s0) / (s1 - s0), 0), 1)
  (1 - w) * fpt_small + w * fpt_series
}

make_circular_design <- function(model, parameters) {
  design(
    model = model,
    factors = list(subjects = 1),
    Rlevels = "resp",
    formula = lapply(parameters, function(p) as.formula(paste0(p, "~ 1"))),
    report_p_vector = FALSE
  )
}

make_snapshot_emc <- function(data, design) {
  make_emc(data, design, type = "single", compress = FALSE, n_chains = 1, verbose = FALSE)
}

cdm_design <- make_circular_design(CDM, c("v", "theta", "a", "t0", "s", "sv"))
sdm_design <- make_circular_design(SDM, c("v", "theta1", "theta2", "a", "t0", "s", "sv"))
hsdm_design <- make_circular_design(HSDM, c("v", "theta1", "theta2", "theta3", "a", "t0", "s", "sv"))
psdm_design <- make_circular_design(PSDM, c("v", "theta1", "a", "t0", "s", "sv"))
phsdm_design <- make_circular_design(PHSDM, c("v", "theta1", "theta2", "a", "t0", "s", "sv"))

cdm_pars <- c(v = log(1.1), theta = qnorm(0.55), a = log(1.0), t0 = log(0.2), s = log(1.0), sv = log(0.1))
sdm_pars <- c(v = log(1.1), theta1 = qnorm(0.35), theta2 = qnorm(0.6), a = log(1.0), t0 = log(0.2), s = log(1.0), sv = log(0.1))
hsdm_pars <- c(v = log(1.0), theta1 = qnorm(0.32), theta2 = qnorm(0.57), theta3 = qnorm(0.66), a = log(1.0), t0 = log(0.2), s = log(1.0), sv = log(0.1))
psdm_pars <- c(v = log(1.0), theta1 = qnorm(0.38), a = log(1.0), t0 = log(0.2), s = log(1.0), sv = log(0.1))
phsdm_pars <- c(v = log(1.0), theta1 = qnorm(0.33), theta2 = qnorm(0.58), a = log(1.0), t0 = log(0.2), s = log(1.0), sv = log(0.1))

RNGkind("L'Ecuyer-CMRG")
set.seed(123)


test_that("CDM simulates and initializes", {
  cdm_data <- make_data(cdm_pars, cdm_design, n_trials = 4)
  cdm_emc <- make_snapshot_emc(cdm_data, cdm_design)
  expect_snapshot(cdm_data)
  expect_snapshot(init_chains(cdm_emc, particles = 10, cores_per_chain = 1)[[1]]$samples)
})

test_that("SDM simulates and initializes", {
  sdm_data <- make_data(sdm_pars, sdm_design, n_trials = 4)
  sdm_emc <- make_snapshot_emc(sdm_data, sdm_design)
  expect_snapshot(sdm_data)
  expect_snapshot(init_chains(sdm_emc, particles = 10, cores_per_chain = 1)[[1]]$samples)
})

test_that("HSDM simulates and initializes", {
  hsdm_data <- make_data(hsdm_pars, hsdm_design, n_trials = 4)
  hsdm_emc <- make_snapshot_emc(hsdm_data, hsdm_design)
  expect_snapshot(hsdm_data)
  expect_snapshot(init_chains(hsdm_emc, particles = 10, cores_per_chain = 1)[[1]]$samples)
})

test_that("PSDM simulates and initializes", {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(123)
  psdm_data <- make_data(psdm_pars, psdm_design, n_trials = 4)
  psdm_emc <- make_snapshot_emc(psdm_data, psdm_design)
  expect_snapshot(psdm_data)
  expect_snapshot(init_chains(psdm_emc, particles = 10, cores_per_chain = 1)[[1]]$samples)
})

test_that("PHSDM simulates and initializes", {
  phsdm_data <- make_data(phsdm_pars, phsdm_design, n_trials = 4)
  phsdm_emc <- make_snapshot_emc(phsdm_data, phsdm_design)
  expect_snapshot(phsdm_data)
  expect_snapshot(init_chains(phsdm_emc, particles = 10, cores_per_chain = 1)[[1]]$samples)
})

test_that("circular family likelihoods are normalized on observed angle coordinates", {
  rt <- 0.9
  tt <- 0.7

  cdm_nat <- matrix(c(0, 0.55, 1.0, 0.2, 1.0, 0), nrow = 1)
  colnames(cdm_nat) <- c("v", "theta", "a", "t0", "s", "sv")
  cdm_grid <- seq(-pi, pi, length.out = 1000)
  cdm_mass <- trapz1d(cdm_grid, dCDM(rep(rt, length(cdm_grid)), cdm_grid, repeat_pars(cdm_nat, length(cdm_grid))))
  expect_equal(cdm_mass, blended_cdm_fpt(tt, a = 1, sigma = 1), tolerance = 1e-4)

  sdm_nat <- matrix(c(0, 0.35, 0.6, 1.0, 0.2, 1.0, 0), nrow = 1)
  colnames(sdm_nat) <- c("v", "theta1", "theta2", "a", "t0", "s", "sv")
  sdm_r <- seq(0, pi, length.out = 150)
  sdm_r2 <- seq(-pi, pi, length.out = 150)
  sdm_mass <- trapz1d(sdm_r, vapply(sdm_r, function(r1) {
    dvals <- dSDM(rep(rt, length(sdm_r2)), rep(r1, length(sdm_r2)), sdm_r2, repeat_pars(sdm_nat, length(sdm_r2)))
    trapz1d(sdm_r2, dvals)
  }, numeric(1)))
  expect_equal(sdm_mass, blended_sdm_fpt(tt, a = 1, sigma = 1), tolerance = 1e-4)

  hsdm_nat <- matrix(c(0, 0.32, 0.57, 0.66, 1.0, 0.2, 1.0, 0), nrow = 1)
  colnames(hsdm_nat) <- c("v", "theta1", "theta2", "theta3", "a", "t0", "s", "sv")
  hsdm_r <- seq(0, pi, length.out = 100)
  hsdm_r2 <- seq(0, pi, length.out = 100)
  hsdm_r3 <- seq(-pi, pi, length.out = 100)
  hsdm_mass <- trapz1d(hsdm_r, vapply(hsdm_r, function(r1) {
    trapz1d(hsdm_r2, vapply(hsdm_r2, function(r2) {
      dvals <- dHSDM(rep(rt, length(hsdm_r3)), rep(r1, length(hsdm_r3)), rep(r2, length(hsdm_r3)), hsdm_r3, repeat_pars(hsdm_nat, length(hsdm_r3)))
      trapz1d(hsdm_r3, dvals)
    }, numeric(1)))
  }, numeric(1)))
  expect_equal(hsdm_mass, blended_hsdm_fpt(tt, a = 1, sigma = 1), tolerance = 1e-4)

  psdm_nat <- matrix(c(0, 0.38, 1.0, 0.2, 1.0, 0), nrow = 1)
  colnames(psdm_nat) <- c("v", "theta1", "a", "t0", "s", "sv")
  psdm_grid <- seq(0, pi, length.out = 1201)
  psdm_mass <- trapz1d(psdm_grid, dPSDM(rep(rt, length(psdm_grid)), psdm_grid, repeat_pars(psdm_nat, length(psdm_grid))))
  expect_equal(psdm_mass, blended_sdm_fpt(tt, a = 1, sigma = 1), tolerance = 1e-4)

  phsdm_nat <- matrix(c(0, 0.33, 0.58, 1.0, 0.2, 1.0, 0), nrow = 1)
  colnames(phsdm_nat) <- c("v", "theta1", "theta2", "a", "t0", "s", "sv")
  phsdm_r <- seq(0, pi, length.out = 100)
  phsdm_r2 <- seq(0, pi, length.out = 100)
  phsdm_mass <- trapz1d(phsdm_r, vapply(phsdm_r, function(r1) {
    dvals <- dPHSDM(rep(rt, length(phsdm_r2)), rep(r1, length(phsdm_r2)), phsdm_r2, repeat_pars(phsdm_nat, length(phsdm_r2)))
    trapz1d(phsdm_r2, dvals)
  }, numeric(1)))
  expect_equal(phsdm_mass, blended_hsdm_fpt(tt, a = 1, sigma = 1), tolerance = 1e-4)
})

test_that("projected models are marginals of their parent models", {
  sdm_nat <- matrix(c(1.3, 0.37, 0.62, 1.1, 0.2, 1.0, 0.15), nrow = 1)
  colnames(sdm_nat) <- c("v", "theta1", "theta2", "a", "t0", "s", "sv")
  rt <- 0.8
  r1 <- 1.1
  hidden_r2 <- seq(-pi, pi, length.out = 1000)
  sdm_marginal <- trapz1d(hidden_r2, dSDM(rep(rt, length(hidden_r2)), rep(r1, length(hidden_r2)), hidden_r2, repeat_pars(sdm_nat, length(hidden_r2))))
  psdm_nat <- sdm_nat[, c("v", "theta1", "a", "t0", "s", "sv"), drop = FALSE]
  expect_equal(dPSDM(rt, r1, psdm_nat), sdm_marginal, tolerance = 1e-4)

  hsdm_nat <- matrix(c(1.2, 0.31, 0.48, 0.66, 1.05, 0.15, 1.0, 0.12), nrow = 1)
  colnames(hsdm_nat) <- c("v", "theta1", "theta2", "theta3", "a", "t0", "s", "sv")
  rt <- 0.75
  r1 <- 1.0
  r2 <- 1.1
  hidden_r3 <- seq(-pi, pi, length.out = 1000)
  hsdm_marginal <- trapz1d(hidden_r3, dHSDM(rep(rt, length(hidden_r3)), rep(r1, length(hidden_r3)), rep(r2, length(hidden_r3)), hidden_r3, repeat_pars(hsdm_nat, length(hidden_r3))))
  phsdm_nat <- hsdm_nat[, c("v", "theta1", "theta2", "a", "t0", "s", "sv"), drop = FALSE]
  expect_equal(dPHSDM(rt, r1, r2, phsdm_nat), hsdm_marginal, tolerance = 1e-4)
})
