withr::local_preserve_seed()

ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
matchfun=function(d)d$S==d$lR
dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2],]
dat$subjects <- droplevels(dat$subjects)

# LNR ---------------------------------------------------------------------
design_LNR <- design(data = dat,model=LNR,matchfun=matchfun,
                     formula=list(m~lM,s~1,t0~1),
                     contrasts=list(m=list(lM=ADmat)), report_p_vector=FALSE)

# RDM ----------------------------------------------------------
design_RDM <- design(data = dat,model=RDM,matchfun=matchfun,
                     formula=list(v~lM,s~1,t0~1,A~1,B~1),
                     contrasts=list(v=list(lM=ADmat)), report_p_vector=FALSE)

# LBA ----------------------------------------------------------
design_LBA <- design(data = dat,model=LBA,matchfun=matchfun,
                     formula=list(v~lM,sv~1,t0~1,A~1,B~1),
                     contrasts=list(v=list(lM=ADmat)), report_p_vector=FALSE)

# WDM ----------------------------------------------------------
design_WDM <- design(data = dat, model=DDM,
                     formula =list(v~1,a~1, t0~1, s~1, Z~1, sv~1, SZ~1),
                     constants=c(s=1, sv=log(0), SZ=qnorm(0)), report_p_vector=FALSE)

# Tests for likelihoods that rely on survivor functions (1-CDFs) ----------
dat_m <- make_missing(dat, LT=0.15, UC=0.75, rt_resolution=0.05) # truncation requires CDFs on each trial, censoring only for NA trials
LNR_s1_m <- make_emc(dat_m, design_LNR, rt_resolution = 0.05, n_chains = 2, compress=TRUE)
RDM_s1_m <- make_emc(dat_m, design_RDM, rt_resolution = 0.05, n_chains = 2, compress=TRUE)
LBA_s1_m <- make_emc(dat_m, design_LBA, rt_resolution = 0.05, n_chains = 2, compress=TRUE)
WDM_s1_m <- make_emc(dat_m, design_WDM, rt_resolution = 0.05, n_chains = 2, compress=TRUE)

# Explicit coefficients on the sampled scale.
p_LNR <- rbind(
  c(m = -1, m_lMd = -.6, s = log(.5), t0 = log(.1)),
  c(m = -.5, m_lMd = -.3, s = log(.8), t0 = log(.05)))
# A = 0 and A > 0 use separate density/CDF implementations.
p_RDM <- rbind(
  c(v = log(2), v_lMd = .6, s = 0, t0 = log(.1), A = -Inf, B = log(.4)),
  c(v = log(2), v_lMd = .6, s = 0, t0 = log(.1), A = log(.2), B = log(.4)))
p_LBA <- rbind(
  c(v = 2, v_lMd = 1, sv = 0, t0 = log(.1), A = log(.2), B = log(.4)),
  c(v = 1.5, v_lMd = .5, sv = log(.7), t0 = log(.05), A = log(.4), B = log(.6)))
p_WDM <- rbind(
  c(v = 1, a = log(2), t0 = log(.1), Z = qnorm(.4)),
  c(v = -1, a = log(1.5), t0 = log(.05), Z = qnorm(.6)))

test_that("LNR with missing", {
  expect_snapshot(round(model_likelihood(LNR_s1_m, p_LNR), 3))
})
test_that("RDM with missing", {
  expect_snapshot(round(model_likelihood(RDM_s1_m, p_RDM), 3))
})
test_that("LBA with missing", {
  expect_snapshot(round(model_likelihood(LBA_s1_m, p_LBA), 3))
})
test_that("WDM with missing", {
  expect_snapshot(round(model_likelihood(WDM_s1_m, p_WDM), 3))
})

# Racing ExGaussian -------------------------------------------------------
design_REXG <- design(data = dat,model=REXG,matchfun=matchfun,
                      formula=list(mu~lM,tau~1,sigma~1,t0~1),
                      contrasts=list(v=list(lM=ADmat)), report_p_vector=FALSE)
REXG_s <- make_emc(dat, design_REXG, rt_resolution = 0.05, n_chains = 2, compress=FALSE)
p_REXG <- rbind(
  c(mu = log(.25), mu_lMTRUE = log(1.3), tau = log(.15), sigma = log(.08), t0 = log(.05)),
  c(mu = log(.3), mu_lMTRUE = log(1.2), tau = log(.25), sigma = log(.12), t0 = log(.1)))
test_that("REXG", {
  expect_snapshot(round(model_likelihood(REXG_s, p_REXG), 3))
})

REXG_s1_m <- make_emc(dat_m, design_REXG, rt_resolution = 0.05, n_chains = 2, compress=TRUE)
test_that("REXG with missing", {
  expect_snapshot(round(model_likelihood(REXG_s1_m, p_REXG), 3))
})
