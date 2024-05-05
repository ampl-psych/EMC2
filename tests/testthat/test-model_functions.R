# When working with lM it is useful to design  an "average and difference"
# contrast matrix, which for binary responses has a simple canonical from:
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
# We also define a match function for lM
matchfun=function(d)d$S==d$lR

# First race models
design_LNR <- make_design(data = forstmann,model=LNR,matchfun=matchfun,
                            formula=list(m~lM + E,s~1,t0~1),
                            contrasts=list(m=list(lM=ADmat)))

design_LBA <- make_design(data = forstmann,model=LBA,matchfun=matchfun,
                            formula=list(v~lM,sv~lM,B~E+lR,A~1,t0~1),
                            contrasts=list(v=list(lM=ADmat)),constants=c(sv=log(1)))

design_RDM <- make_design(data = forstmann,model=RDM,matchfun=matchfun,
                            formula=list(v~lM,s~lM,B~E+lR,A~1,t0~1),
                            contrasts=list(v=list(lM=ADmat)),constants=c(s=log(1)))
# Also DDM
design_DDM <- make_design(data = forstmann,model=DDM,
                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
                            constants=c(s=log(1)))

LNR_s <- make_samplers(forstmann, design_LNR, compress = F, n_chains = 1)
LBA_s <- make_samplers(forstmann, design_LBA, compress = F, n_chains = 1)
RDM_s <- make_samplers(forstmann, design_RDM, compress = F, n_chains = 1)
DDM_s <- make_samplers(forstmann, design_DDM, compress = F, n_chains = 1)


test_that("model estimation", {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(123)
  expect_snapshot(init_chains(LNR_s, particles = 10, cores_for_chains = 1)[[1]]$samples)
  expect_snapshot(init_chains(LBA_s, particles = 10, cores_for_chains = 1)[[1]]$samples)
  expect_snapshot(init_chains(RDM_s, particles = 10, cores_for_chains = 1)[[1]]$samples)
  expect_snapshot(init_chains(DDM_s, particles = 10, cores_for_chains = 1)[[1]]$samples)
})
