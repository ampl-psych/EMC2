# The NLE project's direct-regression baseline reg_s402, registered by path (it is NOT
# shipped: a comparison baseline without a normalisation guarantee). Skipped unless the
# converted card and its source regress.json are given:
#   EMC2_NLE_REG_CARD  = card made by inst/scripts/regress_to_card.R
#   EMC2_NLE_REG_JSON  = the NLE project's reg_s402_regress.json
# Checks the contract gaps of Phase 6a: the card registers with the DDM's defaults, the
# density is the NLE project's DDMreg.R evaluator to 1e-10, there is no CDF, the net has
# no hard support edge below t0, and native and R paths of the likelihood agree.

nle <- function(f) getFromNamespace(f, "EMC2")
reg_card <- Sys.getenv("EMC2_NLE_REG_CARD")
reg_json <- Sys.getenv("EMC2_NLE_REG_JSON")
have_reg <- nzchar(reg_card) && nzchar(reg_json) && file.exists(reg_card) && file.exists(reg_json)

# DDMreg.R's evaluator, restated from the regress.json (GELU tanh, raw R after the scaling)
reg_ref_log_pdf <- function(rt, R, pars, net) {
  gelu <- function(x) 0.5 * x * (1 + tanh(sqrt(2 / pi) * (x + 0.044715 * x^3)))
  raw <- cbind(pars[, "v"], log(pars[, "a"]), log(pars[, "t0"]), log(pars[, "s"]), qnorm(pars[, "Z"]),
               qnorm(pars[, "SZ"]), log(pars[, "sv"]), log(pars[, "st0"]), log(rt))
  x <- cbind(sweep(sweep(raw, 2, net$scaler_mean, "-"), 2, net$scaler_scale, "/"), as.numeric(R))
  nl <- length(net$layers$W)
  for (i in seq_len(nl)) {
    x <- sweep(x %*% net$layers$W[[i]], 2, net$layers$b[[i]], "+")
    if (i < nl) x <- gelu(x)
  }
  as.numeric(x)
}

test_that("reg_s402 registers with the DDM's defaults and no CDF", {
  skip_if_not(have_reg, "EMC2_NLE_REG_CARD / EMC2_NLE_REG_JSON not set")
  m <- register_nn_model(reg_card)()
  expect_identical(m$nn$kind, "regression_joint")
  expect_identical(m$nn$pars, c("v", "a", "t0", "s", "Z", "SZ", "sv", "st0"))
  expect_identical(m$p_types, DDM()$p_types)                 # the defaults trap
  expect_false(m$nn$cdf)
  expect_error(m$pfun(.5, 1L, matrix(1, 1, 8, dimnames = list(NULL, m$nn$pars))), "without a CDF")
  expect_identical(m$nn$sha256, unname(tools::sha256sum(reg_card)))
})

test_that("reg_s402 density = the DDMreg.R evaluator at 1e-10, per-trial parameters", {
  skip_if_not(have_reg, "EMC2_NLE_REG_CARD / EMC2_NLE_REG_JSON not set")
  skip_if_not_installed("jsonlite")
  net <- jsonlite::fromJSON(reg_json, simplifyVector = TRUE)
  m <- register_nn_model(reg_card)()
  reg <- m$nn
  set.seed(6002)
  n <- 300
  pars <- t(replicate(n, reg$lower + (reg$upper - reg$lower) * runif(length(reg$lower), .05, .95)))
  colnames(pars) <- reg$pars
  rt <- exp(runif(n, log(.06), log(3))); R <- sample(1:2, n, TRUE)
  ref <- reg_ref_log_pdf(rt, R, pars, net)
  ok <- ref > -600                                           # below that a density underflows
  expect_gt(sum(ok), 50)
  expect_lt(max(abs(log(m$dfun(rt, R, pars))[ok] - ref[ok])), 1e-10)
})

test_that("reg_s402 has no hard support edge at t0: the density decays steeply below it", {
  skip_if_not(have_reg, "EMC2_NLE_REG_CARD / EMC2_NLE_REG_JSON not set")
  m <- register_nn_model(reg_card)()
  p <- matrix(c(1, .8, .3, 1, .5, .2, .1, .05), 1, dimnames = list(NULL, m$nn$pars))
  ld <- function(rt) log(m$dfun(rt, 1L, p))
  expect_equal(ld(.1), -Inf)                                 # underflows well before t0
  expect_lt(ld(.28), ld(.32) - 30)
  expect_gt(ld(.5), -1)
  # unlike a flow (which puts a little mass before t0), a regression net has no floor: the
  # likelihood pipeline's min_ll is what bounds the log density
})

test_that("reg_s402 cell: native equals the R path in the SBC design (factor S)", {
  skip_if_not(have_reg, "EMC2_NLE_REG_CARD / EMC2_NLE_REG_JSON not set")
  m <- register_nn_model(reg_card)
  des <- design(factors = list(subjects = 1, S = c("left", "right")), Rlevels = c("left", "right"),
                formula = list(v ~ 1, a ~ 1, t0 ~ 1, Z ~ 1, sv ~ 1, SZ ~ 1, st0 ~ 1),
                constants = c(s = log(1)), model = m)
  p <- c(v = 1, a = log(.8), t0 = log(.3), Z = 0, sv = log(.1), SZ = qnorm(.2), st0 = log(.05))
  p <- p[names(sampled_pars(des))]
  set.seed(1)
  dat <- make_data(p, des, n_trials = 60)
  emc <- suppressMessages(make_emc(dat, des, type = "single", n_chains = 1, compress = TRUE,
                                   verbose = FALSE, rt_resolution = .001))
  dadm <- emc[[1]]$data[[1]]; model <- emc[[1]]$model
  P <- matrix(p, 20, length(p), byrow = TRUE, dimnames = list(NULL, names(p))) +
    matrix(rnorm(20 * length(p), sd = .05), 20)
  expect_identical(model()$c_name, "NN")
  llN <- nle("calc_ll_manager")(P, dadm, model)
  ml <- model(); ml$c_name <- NULL
  llR <- nle("calc_ll_manager")(P, dadm, function() ml)
  expect_lt(max(abs(llN - llR) / pmax(1, abs(llR))), 1e-9)
})
