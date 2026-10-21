design_in <- get_design(samples_LNR)[[1]]
data_in <- get_data(samples_LNR)
joint <- make_emc(list(data_in, data_in), list(a = design_in, b = design_in),
                  type = "single", n_chains = 1)

test_that("joint likelihood sums its components", {
  centre <- c(m = -1, m_lMd = 2, s = log(.5), t0 = log(.2))
  parameters <- likelihood_parameters(setNames(rep(centre, 2),
    c(paste0("a|", names(centre)), paste0("b|", names(centre)))))
  ll <- model_likelihood(joint, parameters)
  a <- make_emc(data_in, design_in, type = "single", n_chains = 1)
  b <- make_emc(data_in, design_in, type = "single", n_chains = 1)
  pa <- parameters[, paste0("a|", names(centre))]
  pb <- parameters[, paste0("b|", names(centre))]
  colnames(pa) <- colnames(pb) <- names(centre)
  expect_equal(ll, model_likelihood(a, pa) + model_likelihood(b, pb))
  expect_snapshot(round(ll, 6))
})
