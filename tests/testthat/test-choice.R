matchfun <- function(d) d$S == d$lR

d_ord_probit <- design(Rlevels = c("low", "mid", "high"),
                       factors = list(subjects = 1, S = c("low", "mid", "high")),
                       formula = list(location ~ 0 + S, scale ~ 1, cut ~ 1),
                       matchfun = matchfun,
                       constants = c(scale = log(1)),
                       model = ordered_probit, report_p_vector = FALSE)

p_ord_probit <- c(location_Slow = -1, location_Smid = 0, location_Shigh = 1,
                  cut_lRlow = -0.5, cut_lRmid = log(1.5))

d_ord_logit <- design(Rlevels = c("left", "right"),
                      factors = list(subjects = 1, S = c("left", "right")),
                      formula = list(location ~ 0 + S, scale ~ 1, cut ~ 1),
                      matchfun = matchfun,
                      constants = c(scale = log(1)),
                      model = ordered_logit, report_p_vector = FALSE)

p_ord_logit <- c(location_Sleft = -1, location_Sright = 1, cut = 0)

d_mnl <- design(Rlevels = c("left", "right", "up"),
                factors = list(subjects = 1, S = c("left", "right", "up")),
                formula = list(utility ~ 0 + lR + lM),
                matchfun = matchfun,
                model = multinomial_logit, report_p_vector = FALSE)

p_mnl <- c(utility_lRleft = .5, utility_lRright = 0, utility_lRup = -.5, utility_lMTRUE = 1)

d_mnp <- design(Rlevels = c("left", "right", "up"),
                factors = list(subjects = 1, S = c("left", "right", "up")),
                formula = list(utility ~ 0 + lR + lM),
                matchfun = matchfun,
                model = multinomial_probit, report_p_vector = FALSE)

p_mnp <- c(utility_lRleft = .5, utility_lRright = 0, utility_lRup = -.5, utility_lMTRUE = 1)

# Every stimulus-response combination, independent of the model's simulator.
choice_data <- function(levels) {
  dat <- expand.grid(S = factor(levels, levels = levels),
                     R = factor(levels, levels = levels))
  dat$subjects <- factor(1)
  dat
}

test_that("ordered_probit likelihood", {
  emc <- make_emc(choice_data(c("low", "mid", "high")), d_ord_probit,
                  type = "single", compress = FALSE, n_chains = 1)
  parameters <- likelihood_parameters(p_ord_probit)
  expected <- apply(parameters, 1, function(p) {
    cuts <- c(-Inf, p[["cut_lRlow"]], p[["cut_lRlow"]] + exp(p[["cut_lRmid"]]), Inf)
    sum(vapply(p[1:3], function(location) sum(log(diff(pnorm(cuts - location)))), numeric(1)))
  })
  expect_equal(c(model_likelihood(emc, parameters)), expected, tolerance = 1e-8)
})

test_that("ordered_logit likelihood", {
  emc <- make_emc(choice_data(c("left", "right")), d_ord_logit,
                  type = "single", compress = FALSE, n_chains = 1)
  parameters <- likelihood_parameters(p_ord_logit)
  expected <- apply(parameters, 1, function(p) {
    sum(plogis(p[["cut"]] - p[1:2], log.p = TRUE) +
          plogis(p[["cut"]] - p[1:2], lower.tail = FALSE, log.p = TRUE))
  })
  expect_equal(c(model_likelihood(emc, parameters)), expected, tolerance = 1e-8)
})

test_that("multinomial_logit likelihood", {
  emc <- make_emc(choice_data(c("left", "right", "up")), d_mnl,
                  type = "single", compress = FALSE, n_chains = 1)
  parameters <- likelihood_parameters(p_mnl)
  # A large common offset must cancel before exponentiation.
  shifted <- parameters
  shifted[, 1:3] <- shifted[, 1:3] + 1000
  parameters <- rbind(parameters, shifted)
  expected <- apply(parameters, 1, function(p) {
    sum(vapply(1:3, function(stimulus) {
      utility <- p[1:3] + (1:3 == stimulus) * p[["utility_lMTRUE"]]
      utility <- utility - max(utility)
      sum(utility - log(sum(exp(utility))))
    }, numeric(1)))
  })
  expect_equal(c(model_likelihood(emc, parameters)), expected, tolerance = 1e-8)
})

test_that("multinomial_probit likelihood", {
  emc <- make_emc(choice_data(c("left", "right", "up")), d_mnp,
                  type = "single", compress = FALSE, n_chains = 1)
  parameters <- likelihood_parameters(p_mnp)
  expect_snapshot(round(model_likelihood(emc, parameters), 6))
})
