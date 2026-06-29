# Tests for the race go/nogo plumbing (G-1): a "nogo" response level switches a
# race model to its GNG variant, and make_data withholds the inhibited trials in
# the cens censoring representation (rt = NA + missingness = 4) for the
# downstream CensorSpec likelihood (G-3). The withheld integrand itself is tested
# in test-rdm-gng.R; here we only check the design/data plumbing.

gng_rdm_design <- function() {
  design(factors = list(subjects = 1, S = c("left", "right")), Rlevels = c("go", "nogo"),
         formula = list(v ~ lR, B ~ lR, A ~ 1, t0 ~ 1), constants = c(s = log(1)),
         model = RDM)
}

test_that("a nogo response level switches a race model to RACEGNG / RDMGNG", {
  des <- gng_rdm_design()
  expect_identical(des$model()$type,   "RACEGNG")
  expect_identical(des$model()$c_name, "RDMGNG")
  # plain race model (no nogo level) is unchanged
  plain <- design(factors = list(subjects = 1, S = c("left", "right")),
                  Rlevels = c("left", "right"),
                  formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1), constants = c(s = log(1)),
                  matchfun = function(d) d$S == d$lR, model = RDM)
  expect_identical(plain$model()$type,   "RACE")
  expect_identical(plain$model()$c_name, "RDM")
})

test_that("make_data withholds nogo trials as rt = NA + missingness = 4", {
  set.seed(1)
  des <- gng_rdm_design()
  p <- sampled_pars(des); p[] <- log(1); p["v_lRnogo"] <- -0.3; p["v"] <- 0.5
  dat <- make_data(p, des, n_trials = 120)

  expect_true("missingness" %in% names(dat))
  expect_true(all(is.na(dat$rt[dat$R == "nogo"])))       # every nogo withheld
  expect_true(all(is.finite(dat$rt[dat$R == "go"])))     # every go observed
  expect_true(all(dat$missingness[dat$R == "nogo"] == 4L))
  expect_true(all(is.na(dat$missingness[dat$R == "go"])))
})

test_that("withheld nogo trials reach the dadm with winner = nogo accumulator", {
  set.seed(1)
  des <- gng_rdm_design()
  p <- sampled_pars(des); p[] <- log(1); p["v_lRnogo"] <- -0.3; p["v"] <- 0.5
  dat <- make_data(p, des, n_trials = 120)
  emc <- make_emc(dat, des, type = "single", compress = FALSE, n_chains = 1)
  dadm <- emc[[1]]$data[[1]]

  expect_true(all(c("missingness", "winner", "lR") %in% names(dadm)))
  expect_gt(sum(dadm$missingness == 4L, na.rm = TRUE), 0)
  # for a withheld trial, exactly the nogo accumulator row is the winner
  wt <- dadm$trials[which(dadm$missingness == 4L)[1]]
  rows <- dadm[dadm$trials == wt, ]
  expect_true(all(is.na(rows$rt)))
  expect_identical(as.character(rows$lR[rows$winner]), "nogo")
})

test_that("plain race models are unaffected (no missingness, no withholding)", {
  set.seed(1)
  des <- design(factors = list(subjects = 1, S = c("left", "right")),
                Rlevels = c("left", "right"),
                formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1), constants = c(s = log(1)),
                matchfun = function(d) d$S == d$lR, model = RDM)
  p <- sampled_pars(des); p[] <- log(1); p["v"] <- 0.8
  dat <- make_data(p, des, n_trials = 50)
  expect_false("missingness" %in% names(dat))
  expect_equal(sum(is.na(dat$rt)), 0)
})
