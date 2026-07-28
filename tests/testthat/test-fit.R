# When working with lM it is useful to design  an "average and difference"
# contrast matrix, which for binary responses has a simple canonical from:
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
# We also define a match function for lM
matchfun=function(d)d$S==d$lR

# Drop most subjects
dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2],]
dat$subjects <- droplevels(dat$subjects)


RNGkind("L'Ecuyer-CMRG")
set.seed(123)


test_that("fit", {
  skip_on_os("windows")
  design_LNR <- design(data = dat,model=LNR,matchfun=matchfun,
                       formula=list(m~lM,s~1,t0~1),
                       contrasts=list(m=list(lM=ADmat)), report_p_vector = FALSE)

  LNR_s <- make_emc(dat, design_LNR, rt_resolution = 0.05, n_chains = 2)


  LNR_s <- fit(LNR_s, cores_for_chains = 1, stop_criteria = list(
    preburn = list(iter = 10), burn = list(mean_gd = 2.5), adapt = list(min_unique = 20),
    sample = list(iter = 25)), verbose = FALSE, particle_factor = 30, step_size = 25)
  idx <- LNR_s[[1]]$samples$idx
  expect_snapshot(
    LNR_s[[1]]$samples$theta_mu[,idx], variant = Sys.info()[1]
  )
  expect_snapshot(
    LNR_s[[1]]$samples$alpha[,,idx], variant = Sys.info()[1]
  )
  expect_snapshot(
    LNR_s[[1]]$samples$theta_var[,,idx], variant = Sys.info()[1]
  )
})

# Regression: on_singular = carry_forward must recover a mid-stage singular
# group covariance without corrupting sampler state. The bug this guards is a
# carried-forward alpha stored as a 3-D slice instead of the 2-D matrix
# gibbs_step returns, which broke downstream with "incorrect number of dimensions".
test_that("on_singular carry_forward recovers a mid-stage singular covariance", {
  skip_on_os("windows")
  design_LNR <- design(data = dat, model = LNR, matchfun = matchfun,
                       formula = list(m ~ lM, s ~ 1, t0 ~ 1),
                       contrasts = list(m = list(lM = ADmat)), report_p_vector = FALSE)
  emc <- make_emc(dat, design_LNR, n_chains = 2, compress = FALSE)

  # Force a singular covariance on calls 4-5 (mid chain-1 preburn, so
  # last_good_pars is already set and the carry_forward branch actually runs).
  orig <- EMC2:::gibbs_step; cnt <- 0
  assignInNamespace("gibbs_step", function(sampler, alpha, type, ridge = FALSE, ...) {
    cnt <<- cnt + 1
    if (cnt %in% c(4, 5)) stop("system is computationally singular")
    orig(sampler, alpha, type, ridge = ridge, ...)
  }, ns = "EMC2")
  on.exit(assignInNamespace("gibbs_step", orig, ns = "EMC2"), add = TRUE)

  res <- run_emc(emc, "preburn", stop_criteria = list(iter = 8),
                 cores_for_chains = 1, cores_per_chain = 1, verbose = FALSE,
                 on_singular = list(max_retries = 0, on_exhausted = "carry_forward",
                                    max_carry_forward = 50))
  # Completed all preburn iterations, with alpha stored as a proper p x n x iters array.
  expect_equal(res[[1]]$samples$idx, 9)
  expect_length(dim(res[[1]]$samples$alpha), 3L)
})
