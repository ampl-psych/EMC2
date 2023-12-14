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

lm(v ~ lM, dadm, contrasts = list(lM = ADmat))
lm(b_s ~ E*lR, dadm, contrasts = list(E = Emat, lR = contr.treatment))

N <- 1e6
a <- rnorm(N, 0, 1)
b <- rnorm(N, .5, 2)

sd(a + b)







