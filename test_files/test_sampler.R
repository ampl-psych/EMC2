rm(list = ls())
library(devtools)
library(EMC2)
print(load("test_files/PNAS.RData"))

dat <- data[,c("s","E","S","R","RT")]
names(dat)[c(1,5)] <- c("subjects","rt")
levels(dat$R) <- levels(dat$S)

Emat <- matrix(c(0,-1,0,0,0,-1),nrow=3)
dimnames(Emat) <- list(NULL,c("a-n","a-s"))
Emat

Vmat <- matrix(c(-1,1),ncol=1,dimnames=list(NULL,""))
Vmat

design_at0_full <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),
  Flist=list(v~S,a~E,sv~1, t0~E, st0~1, s~1, Z~1, SZ~1, DP~1),
  constants=c(s=log(1)),
  model=ddmTZD)

dat2 <- dat[which(dat$subjects %in% (unique(dat$subjects)[1:3])),]
dat2 <- droplevels(dat2)

samplers <- make_samplers(dat2, design_at0_full)
samplers_C <- auto_burn(samplers, useC = T, cores_for_chains = 3, cores_per_chain = 3, verbose = T)
samplers <- auto_burn(samplers, useC = T, cores_for_chains = 3, cores_per_chain = 3, verbose = T)


microbenchmark::microbenchmark(
  run_samplers(samplers_C, stage = "burn", iter = 10, useC = T),
  run_samplers(samplers_C, stage = "burn", iter = 10), times = 3
)
