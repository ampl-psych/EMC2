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

# on_singular: recovery from a singular group covariance during sampling.
test_that("resolve_on_singular fills defaults and validates", {
  d <- resolve_on_singular(NULL)
  expect_identical(d$max_retries, 0)
  expect_identical(d$on_exhausted, "error")
  expect_error(resolve_on_singular(list(bad_field = 1)), "unknown")
  expect_error(resolve_on_singular(list(on_exhausted = "nope")))
})

test_that("is_singular_error recognises the covariance failure", {
  expect_true(is_singular_error(simpleError("system is computationally singular: reciprocal condition number = 1e-16")))
  expect_true(is_singular_error(simpleError("Lapack routine dgesv: system is exactly singular")))
  # chol() non-PD failure, both the old and the R >= 4.x wording
  expect_true(is_singular_error(simpleError("the leading minor of order 38 is not positive")))
  expect_true(is_singular_error(simpleError("the leading minor of order 3 is not positive definite")))
  expect_false(is_singular_error(simpleError("some unrelated error")))
})

test_that("robust_gibbs_step retries then gives up", {
  n <- 0
  orig <- EMC2:::gibbs_step
  assignInNamespace("gibbs_step", function(sampler, alpha, type, ...) {
    n <<- n + 1
    if (n <= 2) stop("system is computationally singular")
    list(ok = TRUE)
  }, ns = "EMC2")
  on.exit(assignInNamespace("gibbs_step", orig, ns = "EMC2"), add = TRUE)

  n <- 0
  res <- robust_gibbs_step(NULL, NULL, "standard", resolve_on_singular(list(max_retries = 5)))
  expect_true(isTRUE(res$ok))
  expect_identical(n, 3)                       # 2 failures + 1 success
  n <- 0
  expect_null(robust_gibbs_step(NULL, NULL, "standard", resolve_on_singular(list(max_retries = 1))))
})
