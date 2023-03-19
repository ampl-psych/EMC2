rm(list = ls())
library(devtools)
# Rcpp only works as fast with install like this unfortunately, if this throws errors, restarting Rstudio worked for me
# devtools::install("~/Documents/UVA/2022/EMC2")
load_all()
# load_all()
load("test_files/PNAS.RData")

dat <- data[,c("s","E","S","R","RT")]
names(dat)[c(1,5)] <- c("subjects","rt")
levels(dat$R) <- levels(dat$S)

# Average rate = intercept, and rate d = difference (match-mismatch) contrast
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
ADmat

Emat <- matrix(c(0,-1,0,0,0,-1),nrow=3)
dimnames(Emat) <- list(NULL,c("a-n","a-s"))
Emat

design_RDM <- make_design(
  Ffactors=list(subjects=levels(dat$subjects),S=levels(dat$S),E=levels(dat$E)),
  Rlevels=levels(dat$R),matchfun=function(d)d$S==d$lR,
  Clist=list(lM=ADmat,lR=ADmat,S=ADmat,E=Emat),
  Flist=list(v~lM,s~1,B~E,A~1,t0~1),
  constants=c(s=log(1)),
  model=rdmB)


# Test single subject
dat_single <- dat[which(dat$subjects %in% (unique(dat$subjects)[1])),]
dat_single <- droplevels(dat_single)

# first speed test:
samplers <- make_samplers(dat_single, design_RDM, type = "single", n_chains = 1)
# first let's run init separately
load_all()
debug(run_stage)
samplers <- run_samplers(samplers, stage = "preburn", iter = 1000, cores_per_chain =1, useC = T)


# Note that the useC argument is a bit deceptive. It uses the C dists in any case, but if TRUE will also use the C log_likelihood_race
microbenchmark::microbenchmark(
  run_samplers(samplers, stage = "preburn", iter = 25),
  run_samplers(samplers, stage = "preburn", iter = 25, useC = T), times = 10
) # Literally twice as fast!

# Then compare
samplers <- make_samplers(dat_single, design_RDM, type = "single")
samplersC <- auto_burn(samplers, verbose = T, min_es = 1000, useC = T)
samplersR <- auto_burn(samplers, verbose = T, min_es = 1000, useC = F)
# samplers <- run_samplers(samplers, stage = "preburn", iter = 50, verbose = T, cores_per_chain = 1, cores_for_chains = 1, useC = T)
# samplers <- run_samplers(samplers, stage = "burn", iter = 50, verbose = T, cores_per_chain = 1, cores_for_chains = 1, useC = T)


samplersC_merg <- merge_samples(samplersC)
samplersR_merg <- merge_samples(samplersR)



library(ggplot2)
library(cowplot)

N <- 1500
plots <- list()
for(i in 1:nrow(samplersR_merg$samples$alpha)){
  R <- samplersR_merg$samples$alpha[i,1,( samplersR_merg$samples$idx - N):  samplersR_merg$samples$idx]
  C <- samplersC_merg$samples$alpha[i,1,( samplersC_merg$samples$idx - N):  samplersC_merg$samples$idx]
  df <- data.frame(par_val = c(R, C), group = rep(c("R", "C"), each = N + 1))
  ggp <- ggplot(df, aes(x = par_val, colour = group, fill = group)) +
    geom_density(alpha = .1) + ggtitle(samplersR_merg$par_names[i]) + theme_bw()
  plots[[i]] <- ggp
}
plot_grid(plots[[1]], plots[[2]])
plot_grid(plots[[3]], plots[[4]])
plot_grid(plots[[5]], plots[[6]])
plot_grid(plots[[7]])



