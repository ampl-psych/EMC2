rm(list=ls())
library(devtools)
load_all()
load("~/Documents/UVA/2022/EMC_test/PNAS.RData")
dat <- data[,c("s","E","S","R","RT")]
names(dat)[c(1,5)] <- c("subjects","rt")
levels(dat$R) <- levels(dat$S)
# head(dat)


# Average and difference matrix
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
ADmat

Emat <- matrix(c(0,-1,0,0,0,-1),nrow=3)
dimnames(Emat) <- list(NULL,c("a-n","a-s"))
Emat

#### Models -----

# Mu varies by stimulus, emphasis, and latent match
design_mu <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,E=Emat,S=ADmat),
  Flist=list(m~S*E*lM,s~1,t0~1),
  model=lnrMS)


# Test single subject
dat_single <- dat[which(dat$subjects %in% (unique(dat$subjects)[1])),]
dat_single <- droplevels(dat_single)

samplers <- make_samplers(dat_single, design_mu, type = "single", n_chains = 1)
samplers <- run_emc(samplers, verbose = T, stage = "adapt")

microbenchmark::microbenchmark(
  samplers <- run_samplers(samplers, stage = "preburn", iter = 50, verbose = T, cores_per_chain = 1, cores_for_chains = 1),
  times = 15
)
samplers <- run_samplers(samplers, stage = "burn", iter = 250, verbose = T, cores_per_chain = 1, cores_for_chains = 1)

gd_pmwg(samplers, filter = c("preburn", "burn"), selection = "correlation", mapped = F)

plot_chains(samplers, filter = c("preburn", "burn"))

test_imaginary_function()

