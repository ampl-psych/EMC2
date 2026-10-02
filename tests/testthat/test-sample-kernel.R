# The sample stage runs one fixed kernel, built and tuned at the end of adapt,
# and the hierarchical local proposal follows the current group level through
# the subject's likelihood precision (R/sampling.R, R/fitting.R).

# Conjugate normal model on sufficient statistics: the likelihood precision of
# parameter j is T_j and its linear term T_j * ybar_j, exactly.
conj_ll <- function(pars, dadm, ...) -0.5 * sum(dadm$T * (dadm$ybar - pars)^2)
conj_pars <- c("a", "b", "c", "d")
conj_T <- c(200, 20, 5, 2)
set.seed(11)
conj_dat <- do.call(rbind, lapply(1:4, function(s) data.frame(
  subjects = s, par = conj_pars, T = conj_T, ybar = rnorm(4, 0, .3) + rnorm(4) / sqrt(conj_T))))
conj_dat$subjects <- factor(conj_dat$subjects)
conj_design <- design(model = conj_ll, custom_p_vector = conj_pars, report_p_vector = FALSE)
conj_emc <- make_emc(conj_dat, conj_design, type = "standard", n_chains = 2, compress = FALSE)
conj_stop <- list(preburn = list(iter = 10), burn = list(mean_gd = 2.5), adapt = list(min_unique = 20),
                  sample = list(iter = 20))

test_that("lik_precision recovers the likelihood's quadratic form", {
  dadm <- conj_emc[[1]]$data[[1]]
  centre <- setNames(c(.1, -.2, .3, 0), conj_pars)
  # starting steps far too small and far too large for the likelihood's scale
  for (h in list(rep(1e-4, 4), rep(5, 4))) {
    lik <- EMC2:::lik_precision(centre, h, dadm, conj_emc[[1]]$model)
    expect_equal(unname(lik$prec), diag(conj_T), tolerance = 1e-6)
    expect_equal(unname(lik$lin), conj_T * dadm$ybar, tolerance = 1e-6)
  }
  cond <- EMC2:::conditional_proposal(lik, diag(1 / .09, 4), setNames(rep(0, 4), conj_pars))
  expect_equal(diag(cond$var), 1 / (conj_T + 1 / .09), tolerance = 1e-6)
  expect_equal(unname(cond$mu), conj_T * dadm$ybar / (conj_T + 1 / .09), tolerance = 1e-6)
})

test_that("the sample-stage kernel is fixed", {
  skip_on_os("windows")
  RNGkind("L'Ecuyer-CMRG")
  set.seed(123)
  emc <- fit(conj_emc, cores_for_chains = 1, stop_criteria = conj_stop, verbose = FALSE,
             particle_factor = 20, step_size = 20)
  kernel <- emc[[1]]$sample_kernel
  expect_length(kernel$chains, 2)
  expect_named(kernel$chains[[1]], c("chains_var", "chains_mu", "eff_mu", "eff_var", "prop_var"))
  expect_equal(unname(kernel$lik_prec[[1]]$prec), diag(conj_T), tolerance = 1e-6)
  # the tail of adapt that tunes the kernel is stored as adapt
  expect_gte(chain_n(emc)[1, "adapt"], 100)
  pm <- lapply(emc, function(x) attr(x$samples, "pm_settings"))
  expect_length(pm[[1]][[1]][[1]]$mix, 4)
  # more sample iterations: same proposals, same step size, mixing weights and particles
  emc2 <- fit(emc, cores_for_chains = 1, iter = 40, verbose = FALSE, particle_factor = 20,
              step_size = 20, stop_criteria = list(sample = list(iter = 40)))
  expect_identical(emc2[[1]]$sample_kernel, kernel)
  keep <- c("epsilon", "mix", "n_particles")
  for (ch in 1:2) for (s in 1:4) {
    expect_identical(attr(emc2[[ch]]$samples, "pm_settings")[[s]][[1]][keep], pm[[ch]][[s]][[1]][keep])
  }
})

test_that("the legacy sampler keeps adapting in the sample stage", {
  skip_on_os("windows")
  op <- options(emc.sampler = "legacy"); on.exit(options(op))
  RNGkind("L'Ecuyer-CMRG")
  set.seed(123)
  emc <- fit(conj_emc, cores_for_chains = 1, stop_criteria = conj_stop, verbose = FALSE,
             particle_factor = 20, step_size = 20)
  expect_null(emc[[1]]$sample_kernel)
  expect_null(emc[[1]]$lik_prec)
})
