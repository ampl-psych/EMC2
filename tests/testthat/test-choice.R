RNGkind("L'Ecuyer-CMRG")
set.seed(123)
matchfun <- function(d) d$S == d$lR

d_ord_probit <- design(Rlevels = c("low", "mid", "high"),
                       factors = list(subjects = 1, S = c("low", "mid", "high")),
                       formula = list(location ~ 0 + S, scale ~ 1, cut ~ 1),
                       matchfun = matchfun,
                       constants = c(scale = log(1)),
                       model = ordered_probit)

p_ord_probit <- c(location_Slow = -1, location_Smid = 0, location_Shigh = 1,
                  cut_lRlow = -0.5, cut_lRmid = log(1.5))

d_ord_logit <- design(Rlevels = c("left", "right"),
                      factors = list(subjects = 1, S = c("left", "right")),
                      formula = list(location ~ 0 + S, scale ~ 1, cut ~ 1),
                      matchfun = matchfun,
                      constants = c(scale = log(1)),
                      model = ordered_logit)

p_ord_logit <- c(location_Sleft = -1, location_Sright = 1, cut = 0)

d_mnl <- design(Rlevels = c("left", "right", "up"),
                factors = list(subjects = 1, S = c("left", "right", "up")),
                formula = list(utility ~ 0 + S),
                matchfun = matchfun,
                model = multinomial_logit)

p_mnl <- c(utility_Sleft = 1, utility_Sright = 0, utility_Sup = -1)

d_mnp <- design(Rlevels = c("left", "right", "up"),
                factors = list(subjects = 1, S = c("left", "right", "up")),
                formula = list(utility ~ 0 + S),
                matchfun = matchfun,
                model = multinomial_probit)

p_mnp <- c(utility_Sleft = 1, utility_Sright = 0, utility_Sup = -1)

test_that("ordered_probit", {
  sdat <- make_data(p_ord_probit, d_ord_probit, n_trials = 10)
  ord_probit_s <- make_emc(sdat, d_ord_probit, type = "single", compress = FALSE, n_chains = 1)
  expect_snapshot(init_chains(ord_probit_s, particles = 10, cores_per_chain = 1)[[1]]$samples)
  expect_snapshot(make_data(p_ord_probit, d_ord_probit, n_trials = 10))
})

test_that("ordered_logit", {
  sdat <- make_data(p_ord_logit, d_ord_logit, n_trials = 10)
  ord_logit_s <- make_emc(sdat, d_ord_logit, type = "single", compress = FALSE, n_chains = 1)
  expect_snapshot(init_chains(ord_logit_s, particles = 10, cores_per_chain = 1)[[1]]$samples)
  expect_snapshot(make_data(p_ord_logit, d_ord_logit, n_trials = 10))
})

test_that("multinomial_logit", {
  sdat <- make_data(p_mnl, d_mnl, n_trials = 10)
  mnl_s <- make_emc(sdat, d_mnl, type = "single", compress = FALSE, n_chains = 1)
  expect_snapshot(init_chains(mnl_s, particles = 10, cores_per_chain = 1)[[1]]$samples)
  expect_snapshot(make_data(p_mnl, d_mnl, n_trials = 10))
})

test_that("multinomial_probit", {
  sdat <- make_data(p_mnp, d_mnp, n_trials = 10)
  mnp_s <- make_emc(sdat, d_mnp, type = "single", compress = FALSE, n_chains = 1)
  expect_snapshot(init_chains(mnp_s, particles = 10, cores_per_chain = 1)[[1]]$samples)
  expect_snapshot(make_data(p_mnp, d_mnp, n_trials = 10))
})
