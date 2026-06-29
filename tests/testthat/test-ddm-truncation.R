# Tests for DDM truncation / censoring routed through the shared TruncSpec /
# CensorSpec machinery (the same path the race models use). A constant DDM (~1
# formulas, no between-trial variability) is used so every trial shares the same
# parameters and the truncation/censoring mass is analytic via pDDM.

RNGkind("L'Ecuyer-CMRG")

# Build a constant single-subject DDM design + true parameter vector.
ddm_setup <- function() {
  d0 <- forstmann[forstmann$subjects == unique(forstmann$subjects)[1], ]
  d0$subjects <- droplevels(d0$subjects)
  des <- design(data = d0, model = DDM, formula = list(v ~ 1, a ~ 1, t0 ~ 1, Z ~ 1),
                constants = c(s = log(1), sv = log(0), SZ = qnorm(0), st0 = log(0)),
                report_p_vector = FALSE)
  p <- sampled_pars(des)
  p["v"] <- 1; p["a"] <- log(1.2); p["t0"] <- log(0.2); p["Z"] <- qnorm(0.5)
  list(des = des, p = p)
}

# Summed C++ log-likelihood of `dadm` at the single parameter vector `p`.
ddm_calc_ll <- function(emc, p) {
  dadm  <- emc[[1]]$data[[1]]
  model <- emc[[1]]$model()
  p_types <- names(model$p_types)
  designs <- list()
  for (pp in p_types)
    designs[[pp]] <- attr(dadm, "designs")[[pp]][attr(attr(dadm, "designs")[[pp]], "expand"), , drop = FALSE]
  constants <- attr(dadm, "constants"); if (is.null(constants)) constants <- NA
  p_mat <- matrix(p, nrow = 1, dimnames = list(NULL, names(p)))
  EMC2:::calc_ll(p_mat, dadm, constants = constants, designs = designs, type = model$c_name,
                 model$bound, model$transform, model$pre_transform, p_types = p_types,
                 min_ll = log(1e-10), model$trend)
}

# Analytic in-window mass P(LT < RT < UT) for the constant-DDM natural-scale pars.
ddm_window_mass <- function(LT, UT) {
  pars <- matrix(c(1.2, 1, 0.2, 0.5, 0, 0, 0, 1), nrow = 1,
                 dimnames = list(NULL, c("a", "v", "t0", "Z", "SZ", "sv", "st0", "s")))
  Rl <- factor(c("lo", "hi"), levels = c("lo", "hi"))
  Fhi <- if (is.finite(UT)) pDDM(UT, Rl[1], pars) + pDDM(UT, Rl[2], pars) else 1
  Flo <- if (LT > 0)        pDDM(LT, Rl[1], pars) + pDDM(LT, Rl[2], pars) else 0
  Fhi - Flo
}

test_that("DDM truncation divides the likelihood by the in-window mass", {
  s <- ddm_setup()
  set.seed(1)
  dat <- make_data(s$p, s$des, n_trials = 300)
  LT <- 0.25; UT <- 1.2
  dat_tr <- make_missing(dat, LT = LT, UT = UT, rt_resolution = 1e-10)
  emc <- make_emc(dat_tr, s$des, type = "single", compress = FALSE, n_chains = 1,
                  rt_resolution = 1e-10)

  trunc_ll <- ddm_calc_ll(emc, s$p)
  # same data/rows, truncation switched off -> plain summed log density
  emc_un <- emc
  emc_un[[1]]$data[[1]]$LT <- 0
  emc_un[[1]]$data[[1]]$UT <- Inf
  untrunc_ll <- ddm_calc_ll(emc_un, s$p)

  n_kept <- nrow(emc[[1]]$data[[1]])
  Z <- ddm_window_mass(LT, UT)
  expect_true(is.finite(trunc_ll))
  expect_gt(trunc_ll, untrunc_ll)                       # /Z with Z < 1 raises the density
  expect_equal(trunc_ll, untrunc_ll - n_kept * log(Z), tolerance = 1e-6)
})

test_that("DDM upper censoring contributes the tail mass P(RT > UC)", {
  s <- ddm_setup()
  set.seed(2)
  dat <- make_data(s$p, s$des, n_trials = 300)
  UC <- 1.0
  dat_c <- make_missing(dat, UC = UC, rt_resolution = 1e-10)
  emc <- make_emc(dat_c, s$des, type = "single", compress = FALSE, n_chains = 1,
                  rt_resolution = 1e-10)
  cens_ll <- ddm_calc_ll(emc, s$p)

  # reference: observed-trial density (uncensored subset) + n_cens * log S(UC)
  n_cens  <- sum(dat$rt > UC)
  dat_obs <- dat[dat$rt <= UC, ]
  emc_obs <- make_emc(dat_obs, s$des, type = "single", compress = FALSE, n_chains = 1,
                      rt_resolution = 1e-10)
  obs_ll <- ddm_calc_ll(emc_obs, s$p)
  S_UC   <- ddm_window_mass(UC, Inf)       # mass in (UC, Inf) = P(RT > UC)

  expect_true(is.finite(cens_ll) && n_cens > 0)
  expect_equal(cens_ll, obs_ll + n_cens * log(S_UC), tolerance = 1e-6)
})

test_that("R log_likelihood_ddm truncation twin matches the C++ path", {
  s <- ddm_setup()
  set.seed(4)
  dat <- make_data(s$p, s$des, n_trials = 150)
  dat_tr <- make_missing(dat, LT = 0.25, UT = 1.2, rt_resolution = 1e-10)
  emc <- make_emc(dat_tr, s$des, type = "single", compress = FALSE, n_chains = 1,
                  rt_resolution = 1e-10)
  ll_cpp <- ddm_calc_ll(emc, s$p)
  ll_R   <- EMC2:::calc_ll_R(s$p, emc[[1]]$model(), emc[[1]]$data[[1]])
  expect_equal(ll_cpp, ll_R, tolerance = 1e-8)
})

test_that("DDM without truncation/censoring is unchanged by the new path", {
  s <- ddm_setup()
  set.seed(3)
  dat <- make_data(s$p, s$des, n_trials = 100)
  emc <- make_emc(dat, s$des, type = "single", compress = FALSE, n_chains = 1,
                  rt_resolution = 1e-10)
  ll <- ddm_calc_ll(emc, s$p)
  expect_true(is.finite(ll))
  # no LT/UT/missingness columns -> trunc.any()/censor.any() are FALSE -> plain sum
  dadm <- emc[[1]]$data[[1]]
  expect_false(any(c("LT", "UT") %in% colnames(dadm)) &&
                 (any(dadm$LT != 0) || any(is.finite(dadm$UT))))
})
