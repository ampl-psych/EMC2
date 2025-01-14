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

LNR_s <- make_emc(dat, design_LNR, rt_resolution = 0.05, n_chains = 2)

RNGkind("L'Ecuyer-CMRG")
set.seed(123)
LNR_s <- fit(LNR_s, cores_for_chains = 1, stop_criteria = list(
  preburn = list(iter = 50), burn = list(mean_gd = 1.5), adapt = list(min_unique = 50),
  sample = list(iter = 25)), verbose = FALSE)
idx <- LNR_s[[1]]$samples$idx

test_that("fit", {
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
