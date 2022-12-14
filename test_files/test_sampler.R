library(tidyverse)
library(extraDistr)

I <- 25 # nr of subjects
J <- 400 # trials
K <- 2 # nr of accumulators/response options
Y <- array(NA, c(J,K,I)) # two accumulators
mu_psi <- 0.28 # non-decision time/shift
sigma_psi <- 0.02
mu_drift_match <- 4
sigma_drift_match <- 0.35
mu_drift_mismatch <- 1
sigma_drift_mismatch <- 0.35
mu_sigma <- 0.7
sigma_sigma <- 0.15
boundary <- 1
drift_match <- numeric(I)
drift_mismatch <- numeric(I)
sigma <- numeric(I)
psi <- numeric(I)
for(i in 1:I) {
  drift_match[i] <- rnorm(1, mu_drift_match, sigma_drift_match)
  drift_mismatch[i] <- rnorm(1, mu_drift_mismatch, sigma_drift_mismatch)
  sigma[i] <- rnorm(1, mu_sigma, sigma_sigma)
  psi[i] <- rtnorm(1, mean = mu_psi, sd = sigma_psi, a = 0)
}
dat <- data.frame(I = rep(1:I, each = J),
                      m = rep(NA, I*J),
                      t = rep(NA, I*J))
for(i in 1:I){
  for(j in 1:J){
    Y[j,1,i] <- rwald(1, mu = boundary/drift_match[i], lambda = (boundary/sigma[i])^2)
    Y[j,2,i] <- rwald(1, mu = boundary/drift_mismatch[i], lambda = (boundary/sigma[i])^2)
    dat[dat$I==i,"t"][j] <- psi[i] + min(Y[j,,i])
    dat[dat$I==i,"m"][j] <- which.min(Y[j,,i])
  }
}
colnames(dat) <- c("subjects", "R", "rt")

dat$R <- as.factor(dat$R)
dat$subjects <- as.factor(dat$subjects)
library(devtools)
load_all()

#### Models -----

# Mu varies by stimulus, emphasis, and latent match
design_test <- make_design(
  Ffactors=list(subjects=levels(dat$subjects)),
  Rlevels=levels(dat$R),
  Flist=list(v~R,B~1,s~1,t0~1, A ~ 1),
  constants=c(B=log(1), A = log(0)),
  model=rdmB)


samplers <- make_samplers(dat, design_test, type = "standard", n_chains = 3)
samplers <- run_emc(samplers, verbose = T, cores_for_chains = 3, cores_per_chain = 2)
