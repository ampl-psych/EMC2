matchfun <- function(d)d$S==d$lR
# design an "average and difference" contrast matrix
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"diff"))

# specify design
design_LBABE <- design(data = forstmann,model=LBA,
                formula=list(v~1,sv~1,B~E+lR,A~1,t0~1),
                constants=c(sv=log(1)), report_p_vector=FALSE)

pmean <- c(v=1, B=log(.5),B_Eneutral=log(1.5),
           B_Eaccuracy=log(2),B_lRright=0, A=log(0.25),t0=log(.2))
psd <- c(v=1,
         B=0.3,B_Eneutral=0.3,B_Eaccuracy=0.3,B_lRright=0.3,A=0.4,t0=.5)
prior_LBABE <- prior(design_LBABE, type = 'standard',pmean=pmean,psd=psd)
emc <- make_emc(forstmann,design_LBABE,type="standard",  prior=prior_LBABE,
                          compress = FALSE)

test_that("make_emc", {
  expect_snapshot(str(emc, give.attr = F))
})

# Structural identifiability: a parameter mapped to an all-zero design column has
# no effect on the likelihood and must be flagged before a (potentially long) fit.
test_that("make_emc warns about unidentified (all-zero design column) parameters", {
  ADmat <- matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d"))
  # A 2-level factor contrast whose second column is all zero -> a dead parameter
  dat <- droplevels(forstmann[forstmann$E %in% levels(forstmann$E)[1:2], ])
  dat$E <- droplevels(dat$E)
  Z <- matrix(c(1, -1, 0, 0), ncol = 2, dimnames = list(NULL, c("real", "dead")))
  des <- design(data = dat, model = LBA, matchfun = function(d) d$S == d$lR,
                formula = list(v ~ lM, B ~ E, A ~ 1, t0 ~ 1, sv ~ 1),
                contrasts = list(v = list(lM = ADmat), B = list(E = Z)),
                constants = c(sv = log(1)), report_p_vector = FALSE)

  expect_warning(
    make_emc(dat, des, type = "standard", compress = FALSE),
    "B_Edead"
  )
})

test_that("make_emc does not warn for a well-identified design", {
  ADmat <- matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d"))
  des <- design(data = forstmann, model = LBA, matchfun = function(d) d$S == d$lR,
                formula = list(v ~ lM, B ~ E, A ~ 1, t0 ~ 1, sv ~ 1),
                contrasts = list(v = list(lM = ADmat)),
                constants = c(sv = log(1)), report_p_vector = FALSE)
  expect_no_warning(make_emc(forstmann, des, type = "standard", compress = FALSE))
})

# A per-chain failure inside the parallel sampler must produce an informative
# error naming the chain and reason, not a cryptic downstream crash.
test_that("check_chain_failures reports which chain failed and why", {
  good <- list(list(samples = list(idx = 1)), list(samples = list(idx = 1)))
  expect_null(check_chain_failures(good, "burn", NULL))

  te <- structure("Error\n", class = "try-error",
                  condition = simpleError("matrix not positive definite"))
  expect_error(check_chain_failures(list(good[[1]], te), "burn", NULL),
               "chain 2: matrix not positive definite")
  # a worker killed (NULL slot) is also reported
  expect_error(check_chain_failures(list(good[[1]], NULL), "sample", NULL),
               "Sampling failed in 1 of 2 chain")
})
