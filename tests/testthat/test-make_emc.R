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
                          compress = FALSE, n_chains = 3)   # the snapshot is of 3 chains (the default is now 6)

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

# Collinearity: a linear combination of sampled parameters that maps to zero is
# jointly unidentified even when no single column is all zero. This is the
# generalisation of the all-zero-column check and would not be caught by it.
test_that("make_emc warns about collinear (rank-deficient) design blocks", {
  ADmat <- matrix(c(-1/2, 1/2), ncol = 1, dimnames = list(NULL, "d"))
  # 3-level factor with a rank-2, 3-column contrast: column 'ab' = 'a' + 'b'
  Z <- matrix(c(1, 0, -1,  0, 1, -1,  1, 1, -2), nrow = 3, ncol = 3,
              dimnames = list(NULL, c("a", "b", "ab")))
  des <- design(data = forstmann, model = LBA, matchfun = function(d) d$S == d$lR,
                formula = list(v ~ lM, B ~ E, A ~ 1, t0 ~ 1, sv ~ 1),
                contrasts = list(v = list(lM = ADmat), B = list(E = Z)),
                constants = c(sv = log(1)), report_p_vector = FALSE)
  # names the collinear parameters (B_Ea, B_Eb, B_Eab) in the reported dependency
  expect_warning(make_emc(forstmann, des, type = "standard", compress = FALSE),
                 "B_Eab")
})

# on_singular: recovery from a singular group covariance during sampling.
test_that("resolve_on_singular fills defaults and validates", {
  d <- resolve_on_singular(NULL)
  expect_identical(d$max_retries, 3)
  expect_error(resolve_on_singular(list(bad_field = 1)), "unknown")
})
