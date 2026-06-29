# End-to-end SBC pipeline check for a TRUNCATED DDM: confirms run_sbc generates
# truncated data, fits, and returns well-formed ranks via the new C++ truncated
# likelihood. This is a fast pipeline smoke test (few replicates, short fit) and
# is NOT a statistical calibration assertion -- uniform-rank calibration needs
# ~250 replicates and is run as a separate (NCI) study. The truncation
# *correction itself* is verified exactly in test-ddm-truncation.R.

test_that("run_sbc runs end-to-end on a truncated DDM and returns valid ranks", {
  skip_on_cran()
  set.seed(1)
  d0 <- forstmann[forstmann$subjects == unique(forstmann$subjects)[1], ]
  d0$subjects <- droplevels(d0$subjects)
  des <- design(data = d0, model = DDM, formula = list(v ~ 1, a ~ 1, t0 ~ 1, Z ~ 1),
                constants = c(s = log(1), sv = log(0), SZ = qnorm(0), st0 = log(0)),
                report_p_vector = FALSE)
  pars <- sampled_pars(des)
  pri  <- prior(des, type = "single", mu_mean = pars, mu_sd = rep(0.3, length(pars)))

  out <- run_sbc(des, pri, replicates = 3, trials = 50, LT = 0.2, UT = 1.3,
                 verbose = FALSE, fileName = tempfile(fileext = ".RData"),
                 cores_per_chain = 1, cores_for_chains = 1,
                 stop_criteria = list(iter = 100, selection = "alpha"))

  expect_type(out, "list")
  expect_true(all(c("rank", "med", "bias", "coverage") %in% names(out)))
  ranks <- out$rank[[1]]
  expect_equal(ncol(ranks), length(pars))
  expect_equal(nrow(ranks), 3L)
  expect_true(all(ranks >= 0 & ranks <= 1, na.rm = TRUE))
  expect_false(anyNA(ranks))
})
