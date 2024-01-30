rm(list=ls())
library(EMC2)
# Load and format the PNAS data.
load("test_files/PNAS.RData")
dat <- data[,c("s","E","S","R","RT")]
names(dat)[c(1,5)] <- c("subjects","rt")
levels(dat$R) <- levels(dat$S)

plot_defective_density(dat, factors = c("S", "E"), layout = c(2,3))

#### LBA Models ----

# Average rate = intercept, and rate d = difference (match-mismatch) contrast,
# Note d is analogous to r in the DDM, but use of -1/2 as unlike DDM don't want
# half the difference. In this case the intercept is of interest, it is
# and index of decision urgency and stimulus magnitude effects.
ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d")); ADmat
# Emphasis as before, intercept is overall caution like in DDM
Emat <- matrix(c(0,-1,0,0,0,-1),nrow=3)
dimnames(Emat) <- list(NULL,c("a-n","a-s")); Emat

# Standard assumption, only B affected by E.
# Also allow for response bias, and allow sv to vary by correct incorrect. To fix scale sv=1

debug(make_design)
design_BE <- make_design(data = dat, matchfun=function(d)d$S==d$lR,
                         contrasts=list(lM=ADmat, E=Emat,lR=contr.treatment),
                         formula=list(v~lM,sv~1,B~E*lR,A~1,t0~1),
                         constants=c(sv=log(1)),
                         model=LBA)

load("test_prior.RData")
dadm <- dadm[dadm$subjects == "as1t",]
dadm$v <- 1
dadm$v[dadm$lM == TRUE] <- 2

dadm$b <- log(2)
dadm$b_s <- 1

lm(v ~ lM, dadm, contrasts = list(lM = contr.equalprior))
lm(b_s ~ E*lR, dadm, contrasts = list(E = Emat, lR = contr.treatment))

N <- 1e6
a <- rnorm(N, 0, 1)
b <- rnorm(N, .5, sqrt(2))

var(a + b)

attr(dadm, "designs")

contr.bayes <- function(n, contrasts = TRUE) {
  # This makes a balanced contrasts that assigns equal weight to a fixed effect factor
  if (length(n) <= 1L) {
    if (is.numeric(n) && length(n) == 1L && n > 1L)
      TRUE
    else stop("not enough degrees of freedom to define contrasts")
  } else n <- length(n)
  cont <- diag(n)
  if (contrasts) {
    a <- n
    I_a <- diag(a)
    J_a <- matrix(1, nrow = a, ncol = a)
    Sigma_a <- I_a - J_a/a
    cont <- eigen(Sigma_a)$vectors[,seq_len(a-1), drop = FALSE]
  }
  cont
}

stats::model.matrix(v ~ E, dadm, contrasts.arg = list(E = contr.equalprior))

pars <- rnorm(4)

DM_from <- cbind(1, contr.bayes(4))
DM_to <- cbind(1, contr.treatment(4))
reparameterize <- function(pars, DM_from, DM_to){
  new <- DM_from %*% pars
  new_pars <- solve(t(DM_to) %*% DM_to) %*% t(DM_to) %*% c(new) #OLS
}

new_pars <- reparameterize(pars, DM_from, DM_to)
DM_to %*% c(new_pars)
DM_from %*% pars

prior_predict <- function(DM, means = NULL, N = 1e5){
  if(is.null(means)) means <- rep(0, ncol(DM))
  sds <- 3
  out <- matrix(rnorm(N*ncol(DM), means, sds), nrow = N, ncol = ncol(DM), byrow = T)
  for(i in 1:N){
    out[i,] <- DM %*% out[i,]
  }
  return(out)
}

test <- prior_predict(cbind(1, contr.bayes(5)))
for(i in 1:ncol(test)){
  hist(test[,i])
  print(mean(test[,i]))
  print(var(test[,i]))
}

cbind(1,rbind(0,contr.bayes(3))) %*% c(1, 1.5, 1)

contr.bayes(2)





