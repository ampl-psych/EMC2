# When working with lM it is useful to design  an "average and difference"
# contrast matrix, which for binary responses has a simple canonical from:
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
# We also define a match function for lM
matchfun=function(d)d$S==d$lR

# Drop most subjects
dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2],]
dat$subjects <- droplevels(dat$subjects)

design_LNR <- design(data = dat,model=LNR,matchfun=matchfun,
                          formula=list(m~lM,s~1,t0~1),
                          contrasts=list(m=list(lM=ADmat)))

# create emc objects for the different types
LNR_factor <- make_emc(dat, design_LNR, rt_resolution = 0.05, n_chains = 2, type = "factor",
                       n_factors = 2)
LNR_diag <- make_emc(dat, design_LNR, rt_resolution = 0.05, n_chains = 2, type = "diagonal")
LNR_blocked <- make_emc(dat, design_LNR, rt_resolution = 0.05, n_chains = 2, type = "blocked",
                       par_groups = c(1,2,3,3))
LNR_single <- make_emc(dat, design_LNR, rt_resolution = 0.05, n_chains = 2, type = "single")

# Set seed for mclapply
RNGkind("L'Ecuyer-CMRG")
set.seed(123)
# Run the models for the different variants
# Only preburn, so bit risky since adapt could fail
N <- 50

LNR_factor <- init_chains(LNR_factor, cores_for_chains = 1, particles = 10)
LNR_factor <- run_emc(LNR_factor, cores_for_chains = 1, stop_criteria = list(iter = N), stage = "preburn")

LNR_diag <- init_chains(LNR_diag, cores_for_chains = 1, particles = 10)
LNR_diag <- run_emc(LNR_diag, cores_for_chains = 1, stop_criteria = list(iter = N), stage = "preburn")

LNR_blocked <- init_chains(LNR_blocked, cores_for_chains = 1, particles = 10)
LNR_blocked <- run_emc(LNR_blocked, cores_for_chains = 1, stop_criteria = list(iter = N), stage = "preburn")

LNR_single <- init_chains(LNR_single, cores_for_chains = 1, particles = 10)
LNR_single <- run_emc(LNR_single, cores_for_chains = 1, stop_criteria = list(iter = N), stage = "preburn")


idx <- N + 1

test_that("run_factor", {
  expect_snapshot(
    LNR_factor[[1]]$samples$theta_lambda[,,idx], variant = Sys.info()[1]
  )
  expect_snapshot(
    LNR_factor[[1]]$samples$alpha[,,idx], variant = Sys.info()[1]
  )
  expect_snapshot(
    LNR_factor[[1]]$samples$theta_mu[,idx],variant = Sys.info()[1]
  )
  expect_snapshot(
    LNR_factor[[1]]$samples$theta_var[,,idx], variant = Sys.info()[1]
  )
})

test_that("run_diag", {
  expect_snapshot(
    LNR_diag[[1]]$samples$alpha[,,idx], variant = Sys.info()[1]
  )
  expect_snapshot(
    LNR_diag[[1]]$samples$theta_mu[,idx], variant = Sys.info()[1]
  )
  expect_snapshot(
    LNR_diag[[1]]$samples$theta_var[,,idx], variant = Sys.info()[1]
  )
})

test_that("run_blocked", {
  expect_snapshot(
    LNR_blocked[[1]]$samples$alpha[,,idx], variant = Sys.info()[1]
  )
  expect_snapshot(
    LNR_blocked[[1]]$samples$theta_mu[,idx], variant = Sys.info()[1]
  )
  expect_snapshot(
    LNR_blocked[[1]]$samples$theta_var[,,idx], variant = Sys.info()[1]
  )
})

test_that("run_single", {
  expect_snapshot(
    LNR_single[[1]]$samples$alpha[,,idx], variant = Sys.info()[1]
  )
})

test_that("run_bridge", {
  expect_snapshot(
    compare(list(single = LNR_single, diag = LNR_diag, factor = LNR_factor,
                 blocked = LNR_blocked), stage = "preburn", cores_for_props = 1),
    variant = Sys.info()[1]
  )
})





