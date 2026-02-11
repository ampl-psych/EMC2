RNGkind("L'Ecuyer-CMRG")
set.seed(123)

# Set up fake events
events <- data.frame(
  subjects = rep(1, 10),
  run = rep(1, 10),
  onset = seq(0, 90, by = 10),
  condition = rep(c("A", "B"), 5),
  rt = runif(10, 0.5, 1.5),
  accuracy = sample(0:1, 10, replace = TRUE)
)

# Reshape events to long format
reshaped <- reshape_events(events,
                           event_types = c("condition", "accuracy", "rt"),
                           duration = list(condition = 0.5,
                                          accuracy = 0.2,
                                          rt = function(x) x$rt),
                           modulation = list(rt = 1))

# Set up fake timeseries
ts <- data.frame(
  subjects = rep(1, 100),
  run = rep(1, 100),
  time = seq(0, 99),
  ROI1 = rnorm(100)
)

# Convolve long events
design_matrices <-  convolve_design_matrix(
  timeseries = ts,
  events = reshaped,
  covariates = c('accuracy', 'rt'),
  factors = list(cond = c("condition_A", "condition_B")),
  contrasts = list(cond = matrix(c(-1, 1)))
)

# High pass filter timeseries
ts <- high_pass_filter(ts)

# Design for sampling
des <- design_fmri(design_matrices)

# Test whether prior plotting works
test_that("prior_fmri", {
  vdiffr::expect_doppelganger("prior_fmri",plot(prior(des, type = "single"),map=TRUE, N = 1e2))
})

# Test whether sampling works
fmri_emc <- make_emc(ts, des, type = "single")


test_that("joint", {
  expect_snapshot(init_chains(fmri_emc, particles = 10, cores_for_chains = 1)[[1]]$samples, variant = Sys.info()[1])
})



