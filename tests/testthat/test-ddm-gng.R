# Tests for the C++ DDM go/nogo path (DDMGNG, c_name = "DDMGNG"). A withheld
# trial (rt = NA) contributes 1 - F_go(TIMEOUT); the C++ path must reproduce the
# R log_likelihood_ddmgng exactly (it is the cens-consistent definition).

ddmgng_design <- function() {
  design(Rlevels = c("left", "right"),
         factors = list(subjects = 1, S = c("left", "right")),
         functions = list(
           TIMEOUT = function(d) rep(2.5, nrow(d)),
           Rnogo   = function(d) factor(rep("left",  nrow(d)), levels = c("left", "right")),
           Rgo     = function(d) factor(rep("right", nrow(d)), levels = c("left", "right"))),
         formula = list(v ~ S, a ~ 1, Z ~ 1, t0 ~ 1), model = DDMGNG)
}

test_that("DDMGNG routes to the C++ path and matches R log_likelihood_ddmgng", {
  set.seed(1)
  des <- ddmgng_design()
  p <- sampled_pars(des); p[] <- 0.3; p["a"] <- log(1.2); p["t0"] <- log(0.2); p["Z"] <- qnorm(0.5)
  dat <- make_data(p, des, n_trials = 200)
  expect_gt(sum(is.na(dat$rt)), 0)                      # some withheld (nogo) trials
  emc <- make_emc(dat, des, type = "single", compress = FALSE, n_chains = 1)
  dadm  <- emc[[1]]$data[[1]]
  model <- emc[[1]]$model()
  expect_identical(model$c_name, "DDMGNG")              # C++ dispatch

  p_types <- names(model$p_types)
  designs <- list()
  for (pp in p_types)
    designs[[pp]] <- attr(dadm, "designs")[[pp]][attr(attr(dadm, "designs")[[pp]], "expand"), , drop = FALSE]
  constants <- attr(dadm, "constants"); if (is.null(constants)) constants <- NA
  p_mat <- matrix(p, nrow = 1, dimnames = list(NULL, names(p)))

  ll_cpp <- EMC2:::calc_ll(p_mat, dadm, constants = constants, designs = designs,
                           type = model$c_name, model$bound, model$transform,
                           model$pre_transform, p_types = p_types, min_ll = log(1e-10), model$trend)
  ll_R <- EMC2:::calc_ll_R(p, model, dadm)
  expect_equal(ll_cpp, ll_R, tolerance = 1e-8)
})

test_that("run_sbc runs end-to-end on a DDM go/nogo model", {
  skip_on_cran()
  set.seed(2)
  des <- ddmgng_design()
  pars <- sampled_pars(des)
  pri  <- prior(des, type = "single", mu_mean = pars, mu_sd = rep(0.3, length(pars)))
  out <- run_sbc(des, pri, replicates = 3, trials = 60, verbose = FALSE,
                 fileName = tempfile(fileext = ".RData"),
                 cores_per_chain = 1, cores_for_chains = 1,
                 stop_criteria = list(iter = 100, selection = "alpha"))
  expect_true(all(c("rank", "med", "bias", "coverage") %in% names(out)))
  ranks <- out$rank[[1]]
  expect_equal(nrow(ranks), 3L)
  expect_true(all(ranks >= 0 & ranks <= 1, na.rm = TRUE))
})
