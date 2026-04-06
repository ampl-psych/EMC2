RNGkind("L'Ecuyer-CMRG")
set.seed(123)

calc_lls <- function(emc, n_particles=1e3) {
  model <- emc[[1]]$model()
  p_types <- names(model$p_types)
  dadm <- emc[[1]]$data[[1]]

  designs <- list()
  for(p in p_types){
    designs[[p]] <- attr(dadm,"designs")[[p]][attr(attr(dadm,"designs")[[p]],"expand"),,drop=FALSE]
  }
  constants <- attr(dadm, "constants")
  if(is.null(constants)) constants <- NA

  # make p_mat

  p_mat <- matrix(rnorm(n_particles*length(p_types)), ncol=length(p_types))
  colnames(p_mat) <- p_types

  lls_new <- EMC2:::calc_ll_oo(p_mat, dadm, constants = constants, designs = designs, type = model$c_name,
                               model$bound, model$transform, model$pre_transform, p_types = p_types,
                               min_ll = log(1e-10), model$trend)
  lls_new
}

# Simplest design, no trend -----------------------------------------------
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
matchfun=function(d)d$S==d$lR
dat <- forstmann[forstmann$subjects %in% unique(forstmann$subjects)[1:2],]
dat$subjects <- droplevels(dat$subjects)

# LNR ---------------------------------------------------------------------
design_LNR <- design(data = dat,model=LNR,matchfun=matchfun,
                     formula=list(m~lM,s~1,t0~1),
                     contrasts=list(m=list(lM=ADmat)))
LNR_s1 <- make_emc(dat, design_LNR, rt_resolution = 0.05, n_chains = 2, compress=FALSE)

test_that("LNR", {
  expect_snapshot(calc_lls(LNR_s1))
})

LNR_s2 <- make_emc(dat, design_LNR, rt_resolution = 0.05, n_chains = 2, compress=TRUE)
test_that("LNR compressed", {
  expect_snapshot(calc_lls(LNR_s2))
})



# RDM ----------------------------------------------------------
design_RDM <- design(data = dat,model=RDM,matchfun=matchfun,
                     formula=list(v~lM,s~1,t0~1,A~1,B~1),
                     contrasts=list(v=list(lM=ADmat)))
RDM_s <- make_emc(dat, design_RDM, rt_resolution = 0.05, n_chains = 2, compress=FALSE)
test_that("RDM", {
  expect_snapshot(calc_lls(RDM_s))
})
RDM_s2 <- make_emc(dat, design_RDM, rt_resolution = 0.05, n_chains = 2, compress=TRUE)
test_that("RDM compressed", {
  expect_snapshot(calc_lls(RDM_s2))
})


# LBA ----------------------------------------------------------
design_LBA <- design(data = dat,model=LBA,matchfun=matchfun,
                     formula=list(v~lM,sv~1,t0~1,A~1,B~1),
                     contrasts=list(v=list(lM=ADmat)))
LBA_s <- make_emc(dat, design_LBA, rt_resolution = 0.05, n_chains = 2, compress=FALSE)
test_that("LBA", {
  expect_snapshot(calc_lls(LBA_s))
})
LBA_s2 <- make_emc(dat, design_LBA, rt_resolution = 0.05, n_chains = 2, compress=TRUE)
test_that("LBA compressed", {
  expect_snapshot(calc_lls(LBA_s2))
})


# WDM ----------------------------------------------------------
design_WDM <- design(data = dat, model=DDM,
                     formula =list(v~1,a~1, t0~1, s~1, Z~1, sv~1, SZ~1),
                     constants=c(s=1, sv=log(0), SZ=qnorm(0)))
WDM_s <- make_emc(dat, design_WDM, rt_resolution = 0.05, n_chains = 2)
test_that("WDM", {
  expect_snapshot(calc_lls(WDM_s))
})
WDM_s2 <- make_emc(dat, design_WDM, rt_resolution = 0.05, n_chains = 2, compress=TRUE)
test_that("WDM compressed", {
  expect_snapshot(calc_lls(WDM_s2))
})


# DDM ---
design_DDM <- design(data = dat,model=DDM,
                     formula =list(v~1,a~1, t0~1, s~1, Z~1, sv~1, SZ~1, st0~1),
                     constants=c(s=1))
DDM_s <- make_emc(dat, design_DDM, rt_resolution = 0.05, n_chains = 2)
test_that("DDM", {
  expect_snapshot(calc_lls(DDM_s))
})
