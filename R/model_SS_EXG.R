
#### ExGaussian ----

dexGaussian <- function(rt,pars)
{
  isexp <- pars[,"sigma"] < 1e-4 # shifted exponential
  rt[isexp] <- dexp(rt[isexp]-pars[isexp,"mu"],1/pars[isexp,"tau"])
  isnorm <- !isexp & pars[,"tau"] < 0.05 * pars[,"sigma"] # normal
  rt[isnorm] <- dnorm(rt[isnorm], mean = pars[isnorm,"mu"], sd = pars[isnorm,"sigma"])
  isexg <- !(isexp | isnorm)
  if (any(isexg)) {
    s2 <- pars[isexg,"sigma"]^2
    z <- rt[isexg] - pars[isexg,"mu"] - (s2/pars[isexg,"tau"])
    rt[isexg] <- exp(
      log(pnorm(z/pars[isexg,"sigma"])) -
      log(pars[isexg,"tau"]) -
      (z + (s2/(2 *  pars[isexg,"tau"])))/pars[isexg,"tau"]
    )
  }
  rt
}

pexGaussian <- function(rt,pars)
  # cumulative density for single accumulator
{
  isexp <- pars[,"sigma"] < 1e-4 # shifted exponential
  rt[isexp] <- pexp(rt[isexp]-pars[isexp,"mu"],1/pars[isexp,"tau"])
  isnorm <- !isexp & pars[,"tau"] < 0.05 * pars[,"sigma"] # normal
  rt[isnorm] <- pnorm(rt[isnorm], mean = pars[isnorm,"mu"], sd = pars[isnorm,"sigma"])
  isexg <- !(isexp | isnorm)
  if (any(isexg)) {
    s2 <- pars[isexg,"sigma"]^2
    z <- rt[isexg] - pars[isexg,"mu"] - (s2/pars[isexg,"tau"])
    rt[isexg] <-
      pnorm((rt[isexg] - pars[isexg,"mu"])/pars[isexg,"sigma"]) -
      exp(log(pnorm(z/pars[isexg,"sigma"])) +
        ((pars[isexg,"mu"] + (s2/pars[isexg,"tau"]))^2 - (pars[isexg,"mu"]^2) -
           2 * rt[isexg] * (s2/pars[isexg,"tau"]))/(2 * s2))
  }
  rt
}

rexGaussian <- function(lR,pars,p_types=c("mu","sigma","tau"))
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars is a matrix of corresponding parameter values named as in p_types
  # pars must be sorted so accumulators and parameter for each trial are in
  # contiguous rows.
  #
  # test
  # pars=cbind(mu=c(.5,.6),sigma=c(.1,.1),tau=c(.2,.2)); lR=factor(c(1))
{
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ",paste(p_types,collapse = " "))
  dt <- matrix(
    rnorm(dim(pars)[1],mean=pars[,"mu"],sd=pars[,"sigma"]) +
    rexp(dim(pars)[1],rate=1/pars[,"tau"]),nrow=length(levels(lR)))
  R <- apply(dt,2,which.min)
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  rt <- dt[pick]
  R <- factor(levels(lR)[R],levels=levels(lR))
  cbind.data.frame(R=R,rt=rt)
}

# Go cdf/pdf versions

dexGaussianG <- function(rt,pars)
{
  out <- numeric(length(rt))
  ok <- !is.na(rt)
  out[ok] <- dexGaussian(rt[ok],pars[ok,,drop=FALSE])
  out
}

pexGaussianG <- function(rt,pars)
{
  out <- numeric(length(rt))
  ok <- !is.na(rt)
  out[ok] <- pexGaussian(rt[ok],pars[ok,,drop=FALSE])
  out
}

# Stop cdf/pdf versions

dexGaussianS <- function(rt,pars)
{

  dexGaussS <- function(rt,pars)
  {
    rt <- rt - pars[,"SSD"]
    dimnames(pars)[[2]][dimnames(pars)[[2]]=="muS"] <- "mu"
    dimnames(pars)[[2]][dimnames(pars)[[2]]=="sigmaS"] <- "sigma"
    dimnames(pars)[[2]][dimnames(pars)[[2]]=="tauS"] <- "tau"
    dexGaussian(rt,pars[,c("mu","sigma","tau")])
  }

  out <- numeric(length(rt))
  ok <- !is.na(rt)
  out[ok] <- dexGaussS(rt[ok],pars[ok,,drop=FALSE])
  out
}


pexGaussianS <- function(rt,pars)
{

  pexGaussS <- function(rt,pars)
  {
    rt <- rt - pars[,"SSD"]
    dimnames(pars)[[2]][dimnames(pars)[[2]]=="muS"] <- "mu"
    dimnames(pars)[[2]][dimnames(pars)[[2]]=="sigmaS"] <- "sigma"
    dimnames(pars)[[2]][dimnames(pars)[[2]]=="tauS"] <- "tau"
    pexGaussian(rt,pars[,c("mu","sigma","tau")])
  }

  out <- numeric(length(rt))
  ok <- !is.na(rt)
  out[ok] <- pexGaussS(rt[ok],pars[ok,])
  out
}


# Stop signal random

rSSexGaussian <- function(lR,pars,staircase0,stairstep,stairmin,stairmax,
                          p_types=c("mu","sigma","tau","muS","sigmaS","tauS","tf","gf"))
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars is a matrix of corresponding parameter values named as in p_types
  # pars must be sorted so accumulators and parameter for each trial are in
  # contiguous rows. staircase0 is first SSD in staircase,stairstep is up or
  # down step, stairmin and stairmax are bounds on SSD
{

  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ",paste(p_types,collapse = " "))
  if (!any(levels(lR)=="stop")) stop("lR must have a \"stop\" level")
  # Go failures
  nacc <- length(levels(lR))
  ntrials <- dim(pars)[1]/nacc
  isstop <- lR=="stop"
  isgf <- pars[isstop,"gf"] > runif(ntrials)
  R <- factor(rep("stop",ntrials),levels=levels(lR))
  rt <- rep(NA,ntrials)
  # Not go failures
  nsample <- sum(!isgf)
  dt <- matrix(Inf,nrow=nacc,ncol=nsample)
  isgo <- rep(!isgf,each=nacc) & !isstop
  ngo <- sum(isgo)
  dt[-1,] <- rnorm(ngo,mean=pars[isgo,"mu"],sd=pars[isgo,"sigma"]) +
             rexp(ngo,rate=1/pars[isgo,"tau"])
  SSD <- pars[isstop,"SSD"]
  finiteSSD <- !is.infinite(SSD)
  gostopruns <-  c(finiteSSD & (pars[isstop,"tf"] < runif(ntrials)))[!isgf] # tfs
  if (any(!gostopruns)) { # Run go and trigger failure trials
    r <- apply(dt[,!gostopruns,drop=FALSE],2,which.min)
    R[!isgf][!gostopruns] <- levels(lR)[r]
    pick <- cbind(r,c(1:dim(dt)[2])[!gostopruns]) # Matrix to pick winner
    rt[!isgf][!gostopruns] <- dt[pick]
  }
  if (any(finiteSSD)) { # Run stop trials (and update staircase for gf and tf)
    nstop <- sum(gostopruns)
    dt[1,gostopruns] <- rexp(nstop,rate=1/pars[isstop,"tauS"][!isgf][gostopruns]) +
        rnorm(nstop,mean=pars[isstop,"muS"][!isgf][gostopruns],sd=pars[isstop,"sigmaS"][!isgf][gostopruns])
    isfixed <- gostopruns & !is.na(SSD[!isgf])
    if (any(isfixed)) { # Run fixed trials
      dt[1,isfixed] <- dt[1,isfixed] + SSD[!isgf][isfixed]
      r <- apply(dt[,isfixed,drop=FALSE],2,which.min)
      R[!isgf][isfixed] <- levels(lR)[r]
      pick <- cbind(r,c(1:dim(dt)[2])[isfixed]) # Matrix to pick winner
      rt[!isgf][isfixed] <- dt[pick]
    }
    isstaircase <- is.na(SSD)
    if (any(isstaircase)) {
      sindex <- c(1:ntrials)[isstaircase]
      stairgf <- sindex %in% c(1:ntrials)[isstaircase & isgf]
      dtindex <- c(1:sum(!isgf))[isstaircase[!isgf]]
      SSD[sindex[1]] <- staircase0
      dti <- 1
      for (i in 1:length(sindex)) {
        if (stairgf[i]) SSD <- update_ssd(TRUE,
          sindex[i],sindex[i+1],SSD,stairstep,stairmin,stairmax) else
        {
          dt[1,dtindex[dti]] <- dt[1,dtindex[dti]] + SSD[sindex[i]]
          pick <- which.min(dt[,dtindex[dti]])
          R[sindex[i]] <- levels(lR)[pick]
          rt[sindex[i]] <- dt[pick,dtindex[dti]]
          if ( i != length(sindex) ) {
            SSD <- update_ssd(R[sindex[i]]=="stop",
              sindex[i],sindex[i+1],SSD,stairstep,stairmin,stairmax)
            dti <- dti+1
          }
        }
      }
    }
  }
  rt[R=="stop"] <- NA
  cbind.data.frame(R=R,rt=rt,SSD=SSD)
}




# Following functions moved to C++ model_SS_EXG.cpp

# pEXG <- function (q, mu = 5, sigma = 1, tau = 1, lower_tail = TRUE, log_p = FALSE)
#     # ex-Gaussian cumulative density
#     # Modified from gamlss.dist to make cdf in tau > 0.05 * sigma case robust,
#     # and robust to -Inf and Inf inputs, returns NA for bad sigma or tau and
#     # robust against small sigma cases.
#   {
#     if (sigma <= 0) return(rep(NA,length(q)))
#       if (tau <= 0) return(rep(NA,length(q)))
#
#       # if (sigma < 0.05*tau)
#       if (sigma < 1e-4)
#        return(pexp(q-mu,1/tau,log.p=log_p,lower.tail=lower_tail)) # shfited exponential
#
#     ly <- length(q)
#     sigma <- rep(sigma, length = ly)
#     mu <- rep(mu, length = ly)
#     tau <- rep(tau, length = ly)
#     index <- seq(along = q)
#     z <- q - mu - ((sigma^2)/tau)
#     cdf <- ifelse(is.finite(q),
#       ifelse(tau > 0.05 * sigma,
#          pnorm((q - mu)/sigma) - exp(log(pnorm(z/sigma)) + ((mu + (sigma^2/tau))^2 -
#            (mu^2) -  2 * q * ((sigma^2)/tau))/(2 * sigma^2)),
#          pnorm(q, mean = mu, sd = sigma)),
#         ifelse(q<0,0,1)
#     )
#     if (lower_tail == TRUE)
#       cdf <- cdf
#     else cdf <- 1 - cdf
#     if (log_p == FALSE)
#       cdf <- cdf
#     else cdf <- log(cdf)
#     cdf
#   }
#
# dEXG <- function (x, mu = 5, sigma = 1, tau = 1, log = FALSE)
#   # ex-Gaussian density
#   # gamlss.dist function, but returns NA for bad sigma or tau, and
#   # robust against small sigma cases.
# {
#     if (sigma <= 0) return(rep(NA,length(x)))
#     if (tau <= 0) return(rep(NA,length(x)))
#
#     # if (sigma < 0.05*tau)
#     if (sigma < 1e-4)
#       return(dexp(x-mu,1/tau,log=log)) # shfited exponential
#
#     ly <- length(x)
#     sigma <- rep(sigma, length = ly)
#     mu <- rep(mu, length = ly)
#     tau <- rep(tau, length = ly)
#     z <- x - mu - ((sigma^2)/tau)
#     logfy <- ifelse(tau > 0.05 * sigma,
#       -log(tau) - (z + (sigma^2/(2 *  tau)))/tau + log(pnorm(z/sigma)),
#       dnorm(x, mean = mu, sd = sigma, log = TRUE))
#     if (log == FALSE)
#       fy <- exp(logfy)
#     else fy <- logfy
#     fy
# }
#
# dEXGrace <- function(dt,mu,sigma,tau)
#   # Generates defective PDF for win by first runner, dt (decison time) is
#   # a matrix with length(mu) rows, one row for each runner, and one column
#   # for each decision time for which a defective density value will be
#   # returned.
# {
#   dt[1,] <- dEXG(dt[1,],mu[1],sigma[1],tau[1])
#   if (length(mu)>1) for (i in 2:length(mu))
#     dt[1,] <- dt[1,]*pEXG(dt[i,],mu[i],sigma[i],tau[i],lower_tail=FALSE)
#   dt[1,]
# }
#
#
# stopfn_exg <- function(t,mu,sigma,tau,SSD)
#   # Used by my.integrate, t = vector of times, SSD is a scalar stop-signal delay.
# {
#   dt <- matrix(rep(t+SSD,each=length(mu)),nrow=length(mu))
#   dt[1,] <- dt[1,]-SSD
#   dEXGrace(dt,mu,sigma,tau)
# }

pstopEXG <- function(parstop,n_acc,
  gpars=c("mu","sigma","tau"),spars=c("muS","sigmaS","tauS"))
{

  sindex <- seq(1,dim(parstop)[1],by=n_acc)
  ps <- parstop[sindex,spars]
  SSDs <- parstop[sindex,"SSD"]
  ntrials <- length(SSDs)
  pgo <- array(parstop[-sindex,gpars],dim=c(n_acc-1,ntrials,length(gpars)),
               dimnames=list(NULL,NULL,gpars))
  cells <- character(ntrials)
  for (i in 1:ntrials)
    cells[i] <- paste(SSDs[i],ps[i,],pgo[,i,],collapse="")
  uniq <- !duplicated(cells)
  ups <- numeric(sum(uniq))
  for (i in 1:length(ups)) {
    ups[i] <- my.integrate(f=stopfn_exg,lower=-Inf,SSD=SSDs[uniq][i],
      mu=c(ps[uniq,"muS"][i],pgo[,uniq,"mu",drop=FALSE][,i,1]),
      sigma=c(ps[uniq,"sigmaS"][i],pgo[,uniq,"sigma",drop=FALSE][,i,1]),
      tau=c(ps[uniq,"tauS"][i],pgo[,uniq,"tau",drop=FALSE][,i,1]))
  }
  ups[as.numeric(factor(cells,levels=cells[uniq]))]
}


#### Model list ----
#' Title
#'
#' @return A model list with all the necessary functions to sample
#' @export

SSexGaussian <- function() {
  list(
    type="RACE",
    p_types=c("mu","sigma","tau","muS","sigmaS","tauS","tf","gf"),
    Ntransform=function(x) {
      # transform parameters back to real line
      isprobit <- dimnames(x)[[2]] %in% c("tf","gf")
      x[,!isprobit] <- exp(x[,!isprobit])
      x[,isprobit] <- pnorm(x[,isprobit])
      x
    },
    # p_vector transform
    transform = function(x) x,
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      if (any(names(dadm)=="SSD")) pars <- cbind(pars,SSD=dadm$SSD) else
                                   pars <- cbind(pars,SSD=rep(Inf,dim(pars)[1]))
      attr(pars,"ok") <- (pars[,"tauS"] < 1) & (pars[,"sigmaS"] < 1) &
        (pars[,"tauS"] > 1e-3) & (pars[,"sigmaS"] > 1e-3) & (pars[,"muS"] > 1e-3) &
        ((pars[,"tf"] > 1e-6) | pars[,"tf"] == 0) & ((pars[,"gf"] > 1e-6) | pars[,"gf"] == 0)
      pars
    },
    # Density function (PDF) for single go racer
    dfunG=function(rt,pars) dexGaussianG(rt,pars),
    # Probability function (CDF) for single go racer
    pfunG=function(rt,pars) pexGaussianG(rt,pars),
    # Density function (PDF) for single stop racer
    dfunS=function(rt,pars) dexGaussianS(rt,pars),
    # Probability function (CDF) for single stop racer
    pfunS=function(rt,pars) pexGaussianS(rt,pars),
      # Stop probability integral
    sfun=function(pars,n_acc) pstopEXG(pars,n_acc),
    # Random function for SS race
    rfun=function(lR=NULL,pars) {
      ok <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        (pars[,"tauS"] > 1e-3) & (pars[,"sigmaS"] > 1e-3) & (pars[,"muS"] > 1e-3) &
        (pars[,"tauS"] < 1) & (pars[,"sigmaS"] < 1) &
        ((pars[,"tf"] > 1e-6) | pars[,"tf"] == 0) & ((pars[,"gf"] > 1e-6) | pars[,"gf"] == 0)
      if (is.null(lR)) ok else rSSexGaussian(lR,pars,ok=ok,
      staircase0=.2,stairstep=.05,stairmin=0,stairmax=Inf)
    },
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
      log_likelihood_race_ss(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
  )
}

