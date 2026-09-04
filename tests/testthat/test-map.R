RNGkind("L'Ecuyer-CMRG")
set.seed(123)

dat <- forstmann
dat$covariate <- 1:nrow(forstmann)

des <- design(data = dat, formula = list(v ~ covariate*E, B ~ E, t0 ~ S),
              model = LBA, report_p_vector=FALSE)

set.seed(123)
test_that("mapped_pars", {
  expect_snapshot(mapped_pars(des))
  expect_snapshot(mapped_pars(des, p_vector= rnorm(length(sampled_pars(des)))))
  expect_snapshot(mapped_pars(prior(des, mu_mean = c('v_covariate'  = 1))))
  expect_snapshot(mapped_pars(samples_LNR))
  expect_snapshot(mapped_pars(get_prior(samples_LNR)))
  expect_snapshot(mapped_pars(get_design(samples_LNR)))
})

test_that("map", {
  expect_snapshot(credint(samples_LNR, selection = "mu", map = "E"))
  expect_snapshot(credint(samples_LNR, selection = "mu", map = list(~ E*S)))
  expect_snapshot(credint(samples_LNR, selection = "mu", map = TRUE))
})

# Regression tests for the add_acc / accumulator-collapse handling in make_data
# and mapped_pars. These run without any sampling so they stay CRAN-fast.
test_that("make_data / mapped_pars handle models without an lR column", {
  # DDM-type models never get accumulators, so the collapse must not filter on a
  # (non-existent) lR column and wipe every row.
  ddm_des <- design(factors = list(subjects = 1, S = c("left", "right")),
                    Rlevels = c("left", "right"), model = DDM,
                    formula = list(v ~ S, a ~ 1, Z ~ 1, t0 ~ 1),
                    report_p_vector = FALSE)
  p <- sampled_pars(ddm_des, doMap = FALSE)
  p[] <- c(0, 1, log(1), qnorm(.5), log(.3))

  d <- make_data(p, ddm_des, n_trials = 5)
  expect_gt(nrow(d), 0)
  expect_false("lR" %in% names(d))

  # mapped_pars with a p_vector builds a per-row mapping frame; the DDM path used
  # to collapse on a non-existent lR column and return zero rows.
  mp <- mapped_pars(ddm_des, p_vector = p)
  expect_s3_class(mp, "data.frame")
  expect_gt(nrow(mp), 0)
})

test_that("make_data / mapped_pars expand for functions that reference lR", {
  # A design function referencing an accumulator factor (lR) needs the data
  # expanded to accumulators to be evaluated, then collapsed back to one row per
  # trial for simulation / mapping.
  FF <- function(d) factor(ifelse(as.numeric(d$lR) == 1, "a", "b"),
                           levels = c("a", "b"))
  race_des <- design(factors = list(subjects = 1, S = c("left", "right")),
                     Rlevels = c("left", "right"),
                     matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR),
                     functions = list(FF = FF),
                     formula = list(m ~ FF, s ~ 1, t0 ~ 1), model = LNR,
                     report_p_vector = FALSE)
  # The lR-referencing function must appear in the parameterisation
  expect_true("m_FFb" %in% names(sampled_pars(race_des)))

  p <- sampled_pars(race_des, doMap = FALSE)
  p[] <- c(-0.5, 0.2, log(0.3), log(0.2))

  # One simulated row per trial (2 S x 5 trials), no leftover accumulator columns
  d <- make_data(p, race_des, n_trials = 5)
  expect_identical(nrow(d), 10L)
  expect_false(any(c("lR", "lM", "winner") %in% names(d)))

  mp <- mapped_pars(race_des, p_vector = p)
  expect_s3_class(mp, "data.frame")
  expect_gt(nrow(mp), 0)
})

test_that("par_data_map handles models without an lR column (DDM prior plots)", {
  # par_data_map underlies prior plotting. It expands with add_acc=T but never
  # collapses on lR, so DDM-type models (no accumulators) must still map cleanly.
  ddm_des <- design(factors = list(subjects = 1, S = c("left", "right")),
                    Rlevels = c("left", "right"), model = DDM,
                    formula = list(v ~ S, a ~ 1, Z ~ 1, t0 ~ 1),
                    report_p_vector = FALSE)
  p <- sampled_pars(ddm_des, doMap = FALSE)
  par_mcmc <- array(rnorm(length(p) * 1 * 2), dim = c(length(p), 1, 2),
                    dimnames = list(names(p), "1", NULL))

  res <- par_data_map(par_mcmc, ddm_des, n_trials = 5)
  expect_gt(nrow(res$data), 0)
  expect_false("lR" %in% names(res$data))
  # one mapped-parameter slab per mcmc draw
  expect_identical(dim(res$pars)[2], 2L)
})
