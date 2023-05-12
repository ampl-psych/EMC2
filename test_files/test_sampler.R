rm(list=ls())
devtools::load_all()

print(load("test_files/PNAS.RData"))
dat <- data[,c("s","E","S","R","RT")]
names(dat)[c(1,5)] <- c("subjects","rt")
levels(dat$R) <- levels(dat$S)
head(dat)
# NB: This data has been truncated at 0.25s and 1.5s

# Average rate = intercept, and rate d = difference (match-mismatch) contrast
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
ADmat

Emat <- matrix(c(0,-1,0,0,0,-1),nrow=3)
dimnames(Emat) <- list(NULL,c("a-n","a-s"))
Emat

# Here we fit a series of models
# We'll first build several plausible models and then do model selection
# Only B affected by E
design_B <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM,sv~1,B~E,A~1,t0~1),
  constants=c(sv=log(1)),
  model=lbaB)


samplers <- make_samplers(dat, design_B, nuisance = c(6,7), grouped_pars = 5)

debug(run_stage)
samplers <- run_adapt(samplers, cores_per_chain = 5, cores_for_chains = 1, verbose = T, min_unique = 100)
debug(create_eff_proposals)
samplers <- run_sample(samplers, cores_per_chain = 5, cores_for_chains = 1, verbose = T)


devtools::load_all()
debug(test_adapted)
samplers <- run_adapt(samplers, cores_per_chain = 8, cores_for_chains = 1, verbose = T, min_unique = 40, step_size = 50)
debug(create_eff_proposals)
samplers <- run_sample(samplers, cores_per_chain = 8, cores_for_chains = 1, verbose = T, iter = 200)

debug(merge_samples)
samples <- merge_samples(samplers)

