# Tests for the RDM go/nogo withheld-response integrand (race_gng.h), the G-2
# building block: P(no-go accumulator finishes first in [lower, upper]) =
#   INT f_nogo(t) * prod_go S_go(t) dt, integrated with hcubature.
# rdm_gng_withheld_R(v, B, A, t0, s, nogo, lower, upper) -- nogo is 1-based.

# Reference: the same integrand built from the model's own dRDM/pRDM (= dWald/
# pWald with /s rescaling), integrated with stats::integrate. This pins both the
# parameter scaling (matches the EMC2 RDM likelihood) and the integration.
ref_withheld <- function(v, B, A, t0, s, nogo, lower, upper) {
  mk <- function(k) matrix(c(v[k], B[k], A[k], t0[k], s[k]), nrow = 1,
                           dimnames = list(NULL, c("v", "B", "A", "t0", "s")))
  go <- setdiff(seq_along(v), nogo)
  f <- function(tt) vapply(tt, function(t) {
    dn <- dRDM(t, mk(nogo))
    sg <- 1
    for (j in go) sg <- sg * (1 - pRDM(t, mk(j)))
    dn * sg
  }, numeric(1))
  stats::integrate(f, lower, upper, rel.tol = 1e-9)$value
}

test_that("rdm_gng_withheld matches the model's dRDM/pRDM integrand", {
  grid <- list(
    list(v=c(2,1.5),   B=c(1,1.2),   A=c(0,0),    t0=c(.2,.2),   s=c(1,1)),     # A=0 path
    list(v=c(2,1.5),   B=c(1,1.2),   A=c(.5,.4),  t0=c(.2,.2),   s=c(1,1)),     # A>0 path
    list(v=c(1,2.5),   B=c(.8,1.1),  A=c(.3,0),   t0=c(.15,.25), s=c(1,1)),     # mixed, asymmetric
    list(v=c(2,1.5,1), B=c(1,1.2,.9),A=c(.5,0,.3),t0=c(.2,.15,.25), s=c(1,1,1)) # 3 accumulators
  )
  for (p in grid) for (nogo in seq_along(p$v)) {
    cpp <- rdm_gng_withheld_R(p$v, p$B, p$A, p$t0, p$s, nogo, 0.2, 1.5)
    ref <- ref_withheld(p$v, p$B, p$A, p$t0, p$s, nogo, 0.2, 1.5)
    expect_equal(cpp, ref, tolerance = 1e-5)
  }
})

test_that("winner-finishes-first integral sums to 1 over a proper RDM race", {
  cfgs <- list(
    list(v=c(2,1.5),    B=c(1,1.2),    A=c(0,0),     t0=c(.2,.2),     s=c(1,1)),
    list(v=c(2,1.5),    B=c(1,1.2),    A=c(.5,.4),   t0=c(.2,.2),     s=c(1,1)),
    list(v=c(2,1.5,1),  B=c(1,1.2,.9), A=c(.5,0,.3), t0=c(.2,.15,.25),s=c(1,1,1))
  )
  for (p in cfgs) {
    tot <- sum(vapply(seq_along(p$v),
                      function(k) rdm_gng_withheld_R(p$v,p$B,p$A,p$t0,p$s, k, 0, 60),
                      numeric(1)))
    expect_equal(tot, 1, tolerance = 1e-5)
  }
})

test_that("withheld probability is a valid, monotone interval mass", {
  v<-c(2,1.5); B<-c(1,1.2); A<-c(.4,.3); t0<-c(.2,.2); s<-c(1,1)
  p_small <- rdm_gng_withheld_R(v,B,A,t0,s, 1, 0.3, 0.8)
  p_big   <- rdm_gng_withheld_R(v,B,A,t0,s, 1, 0.3, 2.0)
  expect_true(p_small >= 0 && p_big <= 1)
  expect_gt(p_big, p_small)                      # wider window -> more mass
  expect_equal(rdm_gng_withheld_R(v,B,A,t0,s, 1, 0.5, 0.5), 0, tolerance = 1e-10)  # empty window
})
