#### Staircase ----
#' Assign Stop-Signal Delays (SSDs) to trials
#'
#' @description
#' This function assigns stop-signal delays (SSDs) to trials based on specified probabilities.
#'
#' @param d A data frame containing trial data with a response factor column 'lR'
#' @param SSD A vector of stop-signal delays to be assigned to trials
#' @param pSSD A vector of probabilities for each SSD value. If length is one less than SSD,
#'             the remaining probability is calculated automatically. Default is 0.25.
#'
#' @return A vector of SSDs with the same length as the number of rows in 'd'.
#'         Trials without a stop signal are assigned Inf.
SSD_function <- function(d,SSD=NA,pSSD=.25) {
  if (sum(pSSD)>1) stop("pSSD sum cannot exceed 1.")
  if (length(pSSD)==length(SSD)-1) pSSD <- c(pSSD,1-sum(pSSD))
  if (length(pSSD)!=length(SSD))
    stop("pSSD must be the same length or 1 less than the length of SSD")
  n_acc <- length(levels(d$lR))
  n_trial <- nrow(d)
  out <- rep(Inf,n_trial)
  trials <- c(1:n_trial)
  for (i in 1:length(pSSD)) {
    pick <- sample(trials,floor(pSSD[i]*n_trial))
    out[pick] <- SSD[i]
    trials <- trials[!(trials %in% pick)]
  }
  return(out)
}

staircase_function <- function(dts,staircase) {
  ns <- ncol(dts)
  SSD <- sR <- srt <- numeric()
  SSD[1] <- staircase$SSD0
  for (i in 1:ns) {
    if (SSD[i]<staircase$stairmin) SSD[i] <- staircase$stairmin
    if (SSD[i]>staircase$stairmax) SSD[i] <- staircase$stairmax
    dts[1,i] <- dts[1,i] + SSD[i]
    if (all(is.infinite(dts[-1,i]))) Ri <- 1 else # GF or GF & TF
      Ri <- which.min(dts[,i])
    if (Ri==1) {
      sR[i] <- srt[i] <- NA
      if (i<ns) SSD[i+1] <- round(SSD[i] + staircase$stairstep,3)
    } else {
      sR[i] <- Ri-1
      srt[i] <- min(dts[-1,i])
    if (i<ns) SSD[i+1] <- round(SSD[i] - staircase$stairstep,3)
    }
  }
  list(sR=sR,srt=srt,SSD=SSD)
}


check_staircase <- function(staircase){
  if (!is.list(staircase)){
    staircase <- list(SSD0=.25,stairstep=.05,stairmin=0,stairmax=Inf)
  }
  return(staircase)
}




ST_staircase_function <- function(dts,staircase) {
  ns <- ncol(dts)
  SSD <- sR <- srt <- numeric()
  SSD[1] <- staircase$SSD0
  if (is.null(staircase$accST))  # Indices for accumulator with SSD
    iSSD <- 1 else               # Only stop accumulator
    iSSD <- c(1,staircase$accST) # NB: accST refers to rows in dts
  for (i in 1:ns) {
    if (SSD[i]<staircase$stairmin) SSD[i] <- staircase$stairmin
    if (SSD[i]>staircase$stairmax) SSD[i] <- staircase$stairmax
    dts[iSSD,i] <- dts[iSSD,i] + SSD[i]
    if (all(is.infinite(dts[-1,i]))) { # GF (no ST) or GF & TF
      sR[i] <- srt[i] <- NA
      SSD[i+1] <- round(SSD[i] + staircase$stairstep,3)
    } else {
      Ri <- which.min(dts[,i])
      if (Ri==1) {
        if (!is.null(staircase$accST)) { # ST present
          sR[i] <- staircase$accST[which.min(dts[staircase$accST,i])]-1
          srt[i] <- dts[sR[i]+1,i]
        } else sR[i] <- srt[i] <- NA
        if (i<ns) SSD[i+1] <- round(SSD[i] + staircase$stairstep,3)
      } else {
        sR[i] <- Ri-1
        srt[i] <- dts[Ri,i]
        if (i<ns) {
          if (Ri %in% iSSD) # Stop-triggered response
            SSD[i+1] <- round(SSD[i] + staircase$stairstep,3) else
            SSD[i+1] <- round(SSD[i] - staircase$stairstep,3)
        }
      }
    }
  }
  list(sR=sR,srt=srt,SSD=SSD)
}

#### Single exGaussian functions ----

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

dexGaussian <- function(rt,pars)
  # exGaussian pdf (returns normal or exponential for small tau/sigma)
{
  isexp <- pars[,"sigma"] < 1e-4 # shifted exponential
  rt[isexp] <- dexp(rt[isexp]-pars[isexp,"mu"],1/pars[isexp,"tau"])
  isnorm <- !isexp & pars[,"tau"] < 0.05 * pars[,"sigma"] # normal
  rt[isnorm] <- dnorm(rt[isnorm], mean = pars[isnorm,"mu"],
                                  sd = pars[isnorm,"sigma"])
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
  # exGaussian cdf (returns normal or exponential for small tau/sigma)
{
  isexp <- pars[,"sigma"] < 1e-4 # shifted exponential
  rt[isexp] <- pexp(rt[isexp]-pars[isexp,"mu"],1/pars[isexp,"tau"])
  isnorm <- !isexp & pars[,"tau"] < 0.05 * pars[,"sigma"] # normal
  rt[isnorm] <- pnorm(rt[isnorm], mean = pars[isnorm,"mu"],
                                  sd = pars[isnorm,"sigma"])
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

# following adapted from code released with Tanis et al. (2024)
# https://osf.io/u3k5f/
# DynamicCognitivePsychometrics/dmc/models/WALD-SSEXG/dists.R: functions dEXG and pEXG
# for simplicity assuming no upper truncation

dtexGaussian <- function(rt,pars) {
  out <- Fa <- numeric(length(rt))
  a <- pars[,"exg_lb"]
  ok <- (rt>a) & (rt<Inf) & (a<Inf)
  a_real <- (a>-Inf)
  Fa[ok & a_real] <- pexGaussian(a[ok & a_real], pars[ok & a_real,,drop=FALSE])
  normaliser <- 1 - Fa
  out[ok] <- dexGaussian(rt[ok], pars[ok,,drop=FALSE]) / normaliser[ok]
  return(out)
}

ptexGaussian <- function(rt,pars) {
  out <- Fa <- numeric(length(rt))
  a <- pars[,"exg_lb"]
  ok <- (rt>a) & (rt<Inf) & (a<Inf)
  out[rt==Inf] <- 1
  a_real <- (a>-Inf)
  Fa[ok & a_real] <- pexGaussian(a[ok & a_real], pars[ok & a_real,,drop=FALSE])
  normaliser <- 1 - Fa
  out[ok] <- (pexGaussian(rt[ok], pars[ok,,drop=FALSE]) - Fa[ok]) / normaliser[ok]
  return(out)
}


#### Go Single ExGaussian ----

# Go cdf/pdf (strips out NAs and calls d/pexGaussian)
# Is stripping out rt NAs really necessary?

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

dtexGaussianG <- function(rt,pars)
{
  out <- numeric(length(rt))
  ok <- !is.na(rt)
  out[ok] <- dtexGaussian(rt[ok],pars[ok,,drop=FALSE])
  out
}

ptexGaussianG <- function(rt,pars)
{
  out <- numeric(length(rt))
  ok <- !is.na(rt)
  out[ok] <- ptexGaussian(rt[ok],pars[ok,,drop=FALSE])
  out
}



#### Stop Single ExGaussian ----

# Stop cdf/pdf, subtracts SSD and uses muS/sigmaS/tauS by renaming

dexGaussianS <- function(rt,pars)
{
  rt <- rt - pars[,"SSD"]
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="muS"] <- "mu"
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="sigmaS"] <- "sigma"
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="tauS"] <- "tau"
  dexGaussian(rt,pars)
}


pexGaussianS <- function(rt,pars)
{
  rt <- rt - pars[,"SSD"]
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="muS"] <- "mu"
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="sigmaS"] <- "sigma"
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="tauS"] <- "tau"
  pexGaussian(rt,pars)
}


dtexGaussianS <- function(rt,pars)
{
  rt <- rt - pars[,"SSD"]
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="muS"] <- "mu"
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="sigmaS"] <- "sigma"
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="tauS"] <- "tau"
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="exgS_lb"] <- "exg_lb"
  dtexGaussian(rt,pars)
}


ptexGaussianS <- function(rt,pars)
{
  rt <- rt - pars[,"SSD"]
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="muS"] <- "mu"
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="sigmaS"] <- "sigma"
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="tauS"] <- "tau"
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="exgS_lb"] <- "exg_lb"
  ptexGaussian(rt,pars)
}


#### ExG Race function ----

# Following functions moved to C++ model_SS_EXG.cpp

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

#### ExGaussian random ----

rexG <- function(n,mu,sigma,tau) {
  rnorm(n,mean=mu,sd=sigma) + rexp(n,rate=1/tau)
}

# following copied directly from Tanis et al. 2024 (https://osf.io/u3k5f/)
rTEXG <- function(n,mu=.5,sigma=.05,tau=.1,a=-Inf,b=Inf)
{
  out <- rnorm(n=n, mean=mu,sd=sigma) + rexp(n=n,rate=1/tau)
  ok <- (out > a) & (out < b)
  if (any(!ok)) {
    nneed <- sum(!ok)
    nget <- ceiling(1.1*nneed/mean(ok))
    rtok <- numeric(nneed)
    full <- logical(nneed)
    repeat {
      tmp <- rnorm(n=nget, mean=mu,sd=sigma) + rexp(n=nget,rate=1/tau)
      okrt <- (tmp > a) & (tmp < b)
      if ( sum(okrt)>=nneed ) {
        rtok[!full] <- tmp[okrt][1:nneed]
        break
      }
      ngot <- sum(okrt)
      if (ngot>0) {
        rtok[!full][1:ngot] <- tmp[okrt]
        full[!full][1:ngot] <- TRUE
      }
      nneed <- nneed - ngot
      nget <- ceiling(nneed*2/mean(ok))
    }
    out[!ok] <- rtok
  }
  out
}



rexGaussian <- function(lR,pars,p_types=c("mu","sigma","tau"),
                        ok=rep(TRUE,dim(pars)[1]))
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
  dt <- matrix(rexG(dim(pars)[1],pars[,"mu"],pars[,"sigma"],pars[,"tau"]),
               nrow=length(levels(lR)))
  R <- apply(dt,2,which.min)
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  rt <- dt[pick]
  R <- factor(levels(lR)[R],levels=levels(lR))
  cbind.data.frame(R=R,rt=rt)
}


#### EXG Stop signal random -----

# FIX ME ok NOT IMPLEMENTED

rSSexGaussian <- function(data,pars,ok=rep(TRUE,dim(pars)[1]))
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars must contain an SSD column and an lI column indicating if an
  # accumulator is triggered by the stop signal (ST = stop-triggered).

  # NB1: Go failures will only apply to accumulators where lI = TRUE
  #      and can still have a stop-triggered response on a go-failure trial.
{
  lR <- data$lR
  nacc <- length(levels(lR))   # Does not include stop runner
  ntrials <- dim(pars)[1]/nacc # Number of trials to simulate
  is1 <- lR==levels(lR)[1]     # First go accumulator
  acc <- 1:nacc                # choice (go and ST) accumulator index
  allSSD <- data$SSD[is1]
  # stop-triggered racers
  isST <- pars[,"lI"]==1              # Boolean for all pars
  accST <- acc[pars[1:nacc,"lI"]==1]  # Index of ST accumulator, from 1st trial

  # Setup race among nacc+1 (includes stop accumulator)
  # Default Inf finishing time so if not changed always looses
  dt <- matrix(Inf,nrow=nacc+1,ncol=ntrials)

  # Go failure trials (rep for each accumulator)
  isgf <- rep(pars[is1,"gf"] > runif(ntrials),each=nacc)
  # Go accumulators that don't fail (removing ST)
  isGO <- !isgf & !isST
  ngo <- sum(isGO)

  # Fill in go accumulators
  if (any(isGO)) dt[-1,][isGO] <- rexG(ngo,pars[isGO,"mu"],pars[isGO,"sigma"],
                                       pars[isGO,"tau"])

  # pick out stop trials and races with SSD that is not Inf (i.e., finite or
  # NA, the latter so a staircase can be filled in)
  isStrial <- !is.infinite(pars[is1,"SSD"])
  isSrace <- rep(isStrial,each=nacc)

  # pick out stop trials that are triggered
  isT <- pars[is1,"tf"][isStrial] < runif(sum(isStrial))

  # Logical to pick stop-triggered accumulators that are triggered
  isSTT <- logical(ntrials*nacc) # Default false

  # Pick out stop-triggered accumulators that are triggered
  # NB: can occur on gf trials
  isSTT[isSrace][rep(isT,each=nacc) & isST[isSrace]] <- TRUE
  nst <- sum(isSTT) # Number of triggered ST accumulators

  # Fill in stop-triggered accumulators
  if (any(isSTT)) dt[-1,][isSTT] <- rexG(nst,pars[isSTT,"mu"],
                                      pars[isSTT,"sigma"],pars[isSTT,"tau"])

  # pick out triggered stop racers
  isTS <- logical(ntrials)
  isTS[isStrial][isT] <- TRUE
  ns <- sum(isTS)

  # Fill in stop accumulators
  if (any(isTS)) dt[1,isTS] <- rexG(ns,pars[is1,"muS"][isTS],
                                 pars[is1,"sigmaS"][isTS],pars[is1,"tauS"][isTS])

  # staircase algorithm
  pstair <- is.na(pars[,"SSD"])
  stair <- pstair[is1]
  if (any(stair)) {
    if (is.null(attr(data,"staircase")))
      stop("When SSD has NAs a staircase list must be supplied!")
    staircase <- attr(data,"staircase")
    if (length(accST)>0)
      staircase$accST <- 1+accST
    if (is.null(attr(staircase,"staircase_function"))) {
      if (length(accST)>0)
        stop("Do not use default starircase function with stop-triggered accumulators")
      attr(staircase,"staircase_function") <- staircase_function
    }

    allR <- allrt <- numeric(ncol(dt))  # to store unified results
    dts <- dt[,stair,drop=F]

    # Non-staircase trials
    dt <- dt[,!stair,drop=F]
    spars <- pars[pstair,,drop=F]
    pars <- pars[!pstair,,drop=F]
    isTS <- isTS[!stair]
    isSTT <- isSTT[!pstair]
    isST <- isST[!pstair]
    ntrials <- sum(!stair)
    is1 <- is1[!pstair]
    stair_res <- attr(staircase,"staircase_function")(dts,staircase)
    allR[stair] <- stair_res$sR
    allrt[stair] <- stair_res$srt
  }


  # All SSD already filled in so return R and rt

  # Add SSD to triggered stop accumulators
  if (any(isTS)) dt[1,isTS] <- dt[1,isTS] + pars[is1,"SSD"][isTS]

  # Add SSD to stop-triggered accumulators that are triggered
  if (any(isSTT)) dt[-1,][isSTT] <- dt[-1,][isSTT] + pars[isSTT,"SSD"]

  # R <- factor(rep(NA,ntrials),levels=levels(lR))
  R <- rt <- rep(NA,ntrials)

  # All accumulators Inf (usually when both gf and tf)
  allinf <- apply(dt,2,function(x)all(is.infinite(x)))

  # get winner of stop and go where there is a race
  r <- c(0, acc)[apply(dt[,!allinf,drop=FALSE],2,which.min)]

  # stop wins
  stopwins <- r==0

  # First fill in cases where stop looses, can include wins by ST
  if (any(!stopwins)) {
    rgo <- r[!stopwins]
    R[!allinf][!stopwins] <- rgo
    pick <- cbind(rgo,c(1:sum(!stopwins))) # Matrix to pick winner
    rt[!allinf][!stopwins] <-
      dt[-1,!allinf,drop=FALSE][,!stopwins,drop=FALSE][pick]
  }

  # then if stop wins and kills go, find ST winner
  if (any(isST) & any(stopwins)) {
    # stop triggered accumulators that are racing
    rst <- dt[-1,!allinf,drop=FALSE][accST,stopwins,drop=FALSE]
    # stop-triggered winners
    rtw <- apply(rst,2,which.min)
    # index for stop-triggered
    R[!allinf][stopwins] <- accST[rtw]
    pick <- cbind(rtw,1:ncol(rst))
    rt[!allinf][stopwins] <- rst[pick]
  }

  rt[is.na(R)] <- NA

  if (any(stair)) {
    allrt[!stair] <- rt
    allR[!stair] <- R
    allSSD[stair] <- stair_res$SSD
    out <- cbind.data.frame(R=factor(allR,levels=1:nacc,labels=levels(lR)),
                            rt=allrt, SSD = allSSD)
    return(out)
  }
  cbind.data.frame(R=factor(R,levels=1:nacc,labels=levels(lR)),rt=rt)
}


#### ExG stop probability (no stop triggered) ----

my.integrate_old <- function(...,upper=Inf,big=10)
  # Avoids bug in integrate upper=Inf that uses only 1  subdivision
  # Use of  big=10 is arbitrary ...
{
  out <- try(integrate(...,upper=upper),silent=TRUE)
  if (inherits(out, "try-error")) 0 else
  {
    if (upper==Inf & out$subdivisions==1)
    {
      out <- try(integrate(...,upper=big),silent=TRUE)
      if (inherits(out, "try-error")) 0 else
      {
        if (out$subdivisions==1) 0 else out$value
      }
    } else out$value
  }
}


pstopEXG_old <- function(parstop,n_acc,upper=Inf,
                     gpars=c("mu","sigma","tau"),spars=c("muS","sigmaS","tauS"))
{
  sindex <- seq(1,nrow(parstop),by=n_acc)  # Stop accumulator index
  ps <- parstop[sindex,spars,drop=FALSE]   # Stop accumulator parameters
  SSDs <- parstop[sindex,"SSD",drop=FALSE] # SSDs
  ntrials <- length(SSDs)
  if (length(upper)==1) upper <- rep(upper,length.out=ntrials)
  pgo <- array(parstop[,gpars],dim=c(n_acc,ntrials,length(gpars)),
               dimnames=list(NULL,NULL,gpars))
  cells <- apply(
    cbind(SSDs,ps,upper,matrix(as.vector(aperm(pgo,c(2,1,3))),nrow=ntrials))
  ,1,paste,collapse="")
  # cells <- character(ntrials)
  # for (i in 1:ntrials)
  #   cells[i] <- paste(SSDs[i],ps[i,],pgo[,i,],upper[i],collapse="")
  uniq <- !duplicated(cells)
  ups <- sapply(1:sum(uniq),function(i){
    my.integrate_old(f=stopfn_exg_old,lower=-Inf,SSD=SSDs[i],upper=upper[i],
                           mu=c(ps[i,"muS"],pgo[,i,"mu"]),
                           sigma=c(ps[i,"sigmaS"],pgo[,i,"sigma"]),
                           tau=c(ps[i,"tauS"],pgo[,i,"tau"]))
  })
  ups[as.numeric(factor(cells,levels=cells[uniq]))]
}

pstopEXG <- function(parstop,n_acc,upper=Inf,
                     gpars=c("mu","sigma","tau"),spars=c("muS","sigmaS","tauS"))
{
  sindex <- seq(1,nrow(parstop),by=n_acc)  # Stop accumulator index
  ps <- parstop[sindex,spars,drop=FALSE]   # Stop accumulator parameters
  SSDs <- parstop[sindex,"SSD",drop=FALSE] # SSDs
  ntrials <- length(SSDs)
  if (length(upper)==1) upper <- rep(upper,length.out=ntrials)
  pgo <- array(parstop[,gpars],dim=c(n_acc,ntrials,length(gpars)),
               dimnames=list(NULL,NULL,gpars))
  cells <- apply(
    cbind(SSDs,ps,upper,matrix(as.vector(aperm(pgo,c(2,1,3))),nrow=ntrials))
  ,1,paste,collapse="")
  # cells <- character(ntrials)
  # for (i in 1:ntrials)
  #   cells[i] <- paste(SSDs[i],ps[i,],pgo[,i,],upper[i],collapse="")
  uniq <- !duplicated(cells)
  ups <- sapply(1:sum(uniq),function(i){
    my.integrate(f=stopfn_exg,lower=-Inf,SSD=SSDs[i],upper=upper[i],
                           mu=c(ps[i,"muS"],pgo[,i,"mu"]),
                           sigma=c(ps[i,"sigmaS"],pgo[,i,"sigma"]),
                           tau=c(ps[i,"tauS"],pgo[,i,"tau"]))
  })
  ups[as.numeric(factor(cells,levels=cells[uniq]))]
}


pstopTEXG <- function(
    parstop, n_acc, upper=Inf,
    gpars=c("mu","sigma","tau","exg_lb"), spars=c("muS","sigmaS","tauS","exgS_lb")
) {
  sindex <- seq(1,nrow(parstop),by=n_acc)  # Stop accumulator index
  ps <- parstop[sindex,spars,drop=FALSE]   # Stop accumulator parameters
  SSDs <- parstop[sindex,"SSD",drop=FALSE] # SSDs
  ntrials <- length(SSDs)
  if (length(upper)==1) upper <- rep(upper,length.out=ntrials)
  pgo <- array(parstop[,gpars],dim=c(n_acc,ntrials,length(gpars)),
               dimnames=list(NULL,NULL,gpars))
  cells <- apply(
    cbind(SSDs,ps,upper,matrix(as.vector(aperm(pgo,c(2,1,3))),nrow=ntrials))
    ,1,paste,collapse="")
  # cells <- character(ntrials)
  # for (i in 1:ntrials)
  #   cells[i] <- paste(SSDs[i],ps[i,],pgo[,i,],upper[i],collapse="")
  uniq <- !duplicated(cells)
  ups <- sapply(1:sum(uniq),function(i){
    my.integrate_old(f=stopfn_texg,lower=ps[i,"exgS_lb"],SSD=SSDs[i],upper=upper[i],
                     mu=c(ps[i,"muS"],pgo[,i,"mu"]),
                     sigma=c(ps[i,"sigmaS"],pgo[,i,"sigma"]),
                     tau=c(ps[i,"tauS"],pgo[,i,"tau"]),
                     lb=c(ps[i,"exgS_lb"],pgo[,i,"exg_lb"]))
  })
  ups[as.numeric(factor(cells,levels=cells[uniq]))]
}

# #### Stop probability stop triggered ----
#
#
# stopfn_exgST <- function(t,mu,sigma,tau,SSD,st=1)
#   # Used by my.integrate, t = vector of times, SSD is a scalar stop-signal delay.
#   # st is a vector of indices for the stop and stop-triggered accumulators
# {
#   dt <- matrix(rep(t+SSD,each=length(mu)),nrow=length(mu))
#   dt[st,] <- dt[st,]-SSD
#   dEXGrace(dt,mu,sigma,tau)
# }
#
#
# pstopEXGST <- function(parstop,n_acc,upper=Inf,st=1,
#                      gpars=c("mu","sigma","tau"),spars=c("muS","sigmaS","tauS"))
# {
#   sindex <- seq(1,nrow(parstop),by=n_acc)
#   ps <- parstop[sindex,spars,drop=FALSE]
#   SSDs <- parstop[sindex,"SSD",drop=FALSE]
#   ntrials <- length(SSDs)
#   if (length(upper)==1) upper <- rep(upper,length.out=ntrials)
#   pgo <- array(parstop[,gpars],dim=c(n_acc,ntrials,length(gpars)),
#                dimnames=list(NULL,NULL,gpars))
#   cells <- apply(cbind(SSDs,ps,upper,
#     matrix(as.vector(aperm(pgo,c(2,1,3))),nrow=ntrials)),1,paste,collapse="")
#   # cells <- character(ntrials)
#   # for (i in 1:ntrials)
#   #   cells[i] <- paste(SSDs[i],ps[i,],pgo[,i,],upper[i],collapse="")
#   uniq <- !duplicated(cells)
#   ups <- sapply(1:sum(uniq),function(i){
#     my.integrate(f=stopfn_exgST,lower=-Inf,SSD=SSDs[i],upper=upper[i],
#                            mu=c(ps[i,"muS"],pgo[,i,"mu"]),
#                            sigma=c(ps[i,"sigmaS"],pgo[,i,"sigma"]),
#                            tau=c(ps[i,"tauS"],pgo[,i,"tau"]),st=st)
#   })
#   ups[as.numeric(factor(cells,levels=cells[uniq]))]
# }



#### SSexG Model list ----

#' Stop-signal ex-Gaussian race
#'
#' Model file to estimate the ex-Gaussian race model for Stop-Signal data
#'
#' Model files are almost exclusively used in `design()`.
#'
#' @details
#'
#' Default values are used for all parameters that are not explicitly listed in the `formula`
#' argument of `design()`.They can also be accessed with `SSexG()$p_types`.
#'
#' @details
#' | **Parameter** | **Transform** | **Natural scale** | **Default**   | **Mapping**                    | **Interpretation**                                            |
#' |-----------|-----------|---------------|-----------|----------------------------|-----------------------------------------------------------|
#' | *mu*       | log         | \[0, Inf\]     | log(.4)         |                            | mu parameter of ex-Gaussian go finishing time distribution              |
#' | *sigma*       | log       | \[0, Inf\]        | log(.05)    |                            | sigma parameter of ex-Gaussian go finishing time distribution                                        |
#' | *tau*      | log       | \[0, Inf\]        | log(.1)    |                            | tau parameter of ex-Gaussian go finishing time distribution                                          |
#' | *muS*       | log       | \[0, Inf\]        | log(.3)    |                            | mu parameter of ex-Gaussian stopping finishing time distribution           |
#' | *sigmaS*       | log    | \[0, Inf\]        | log(.025)|                   | sigma parameter of ex-Gaussian stopping finishing time distribution                              |
#' | *tauS*      | log    | \[0, Inf\]        | log(.05)  |  | tau parameter of ex-Gaussian stopping finishing time distribution       |
#' | *tf*      | probit       | \[0, 1\]        | qnorm(0)    |                            | Trigger failure probability           |
#' | *gf*     | probit       | \[0, 1\]        | qnorm(0)    |                            | Go failure probability    |
#'
#'
#'
#' @return A model list with all the necessary functions to sample
#' @export
SSexG_old <- function() {
  list(
    type="RACE",
    p_types=c(mu=log(.4),sigma=log(.05),tau=log(.1),
              muS=log(.3),sigmaS=log(.025),tauS=log(.05),tf=qnorm(0),gf=qnorm(0)),
    transform=list(func=c( mu = "exp",  sigma = "exp",  tau = "exp",
                          muS = "exp", sigmaS = "exp", tauS = "exp",
                          tf="pnorm",gf="pnorm")),
    bound=list(minmax=cbind( mu=c(0,Inf),  sigma=c(1e-4,Inf),  tau=c(1e-4,Inf),
                            muS=c(0,Inf), sigmaS=c(1e-4,Inf), tauS=c(1e-4,Inf),
                            tf=c(.001,.999),gf=c(.001,.999)),
                            exception=c(tf=0,gf=0)),
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      # if (any(names(dadm)=="SSD")) pars <- cbind(pars,SSD=dadm$SSD) else
      #   pars <- cbind(pars,SSD=rep(Inf,dim(pars)[1]))
      pars <- cbind(pars,SSD=dadm$SSD)
      pars <- cbind(pars,lI=as.numeric(dadm$lI))  # Only necessary for data generation.
      pars
    },
    # Density function (PDF) for single go racer
    dfunG=function(rt,pars) dexGaussianG(rt,pars),
    # Probability function (CDF) for single go racer
    pfunG=function(rt,pars) pexGaussianG(rt,pars),
    # Density function (PDF) for single stop racer
    dfunS=function(rt,pars)
      dexGaussianS(rt,pars[,c("muS","sigmaS","tauS","SSD"),drop=FALSE]),
    # Probability function (CDF) for single stop racer
    pfunS=function(rt,pars)
      pexGaussianS(rt,pars[,c("muS","sigmaS","tauS","SSD"),drop=FALSE]),
    # Stop probability integral
    sfun=function(pars,n_acc,upper=Inf) pstopEXG_old(pars,n_acc,upper=upper),
    # Random function for SS race
    rfun=function(data=NULL,pars) {
      rSSexGaussian(data,pars,ok=attr(pars, "ok"))
    },
    # Race likelihood combining pfun and dfun
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10))
      log_likelihood_race_ss_old(pars, dadm, model, min_ll = min_ll)
  )
}


#' Stop-signal ex-Gaussian race
#'
#' Model file to estimate the ex-Gaussian race model for Stop-Signal data
#'
#' Model files are almost exclusively used in `design()`.
#'
#' @details
#'
#' Default values are used for all parameters that are not explicitly listed in the `formula`
#' argument of `design()`.They can also be accessed with `SSexG()$p_types`.
#'
#' @details
#' | **Parameter** | **Transform** | **Natural scale** | **Default**   | **Mapping**                    | **Interpretation**                                            |
#' |-----------|-----------|---------------|-----------|----------------------------|-----------------------------------------------------------|
#' | *mu*       | log         | \[0, Inf\]     | log(.4)         |                            | mu parameter of ex-Gaussian go finishing time distribution              |
#' | *sigma*       | log       | \[0, Inf\]        | log(.05)    |                            | sigma parameter of ex-Gaussian go finishing time distribution                                        |
#' | *tau*      | log       | \[0, Inf\]        | log(.1)    |                            | tau parameter of ex-Gaussian go finishing time distribution                                          |
#' | *muS*       | log       | \[0, Inf\]        | log(.3)    |                            | mu parameter of ex-Gaussian stopping finishing time distribution           |
#' | *sigmaS*       | log    | \[0, Inf\]        | log(.025)|                   | sigma parameter of ex-Gaussian stopping finishing time distribution                              |
#' | *tauS*      | log    | \[0, Inf\]        | log(.05)  |  | tau parameter of ex-Gaussian stopping finishing time distribution       |
#' | *tf*      | probit       | \[0, 1\]        | qnorm(0)    |                            | Trigger failure probability           |
#' | *gf*     | probit       | \[0, 1\]        | qnorm(0)    |                            | Go failure probability    |
#'
#'
#' @return A model list with all the necessary functions to sample
#' @export
SSexG <- function() {
  list(
    type = "RACE",
    c_name = "SS_EXG",
    p_types = c(
      mu = log(.4), sigma = log(.05), tau = log(.1),
      muS = log(.3), sigmaS = log(.025), tauS = log(.05),
      tf = qnorm(0), gf = qnorm(0)
    ),
    transform = list(
      func = c(
        mu = "exp", sigma = "exp", tau = "exp",
        muS = "exp", sigmaS = "exp", tauS = "exp",
        tf = "pnorm", gf = "pnorm"
      )
    ),
    bound = list(
      minmax = cbind(
        mu = c(0, Inf), sigma = c(1e-4, Inf), tau = c(1e-4, Inf),
        muS = c(0, Inf), sigmaS = c(1e-4, Inf), tauS = c(1e-4, Inf),
        tf = c(.001, .999), gf = c(.001, .999)
      ),
      exception = c(
        tf = 0, gf = 0
      )
    ),
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      # if (any(names(dadm)=="SSD")) pars <- cbind(pars,SSD=dadm$SSD) else
      #   pars <- cbind(pars,SSD=rep(Inf,dim(pars)[1]))
      pars <- cbind(pars, SSD = dadm$SSD)
      pars <- cbind(pars, lI = as.numeric(dadm$lI))  # Only necessary for data generation.
      return(pars)
    },
    # Density function (PDF) for single go racer
    dfunG = function(rt, pars) return(dexGaussianG(rt, pars)),
    # Probability function (CDF) for single go racer
    pfunG = function(rt, pars) return(pexGaussianG(rt, pars)),
    # Density function (PDF) for single stop racer
    dfunS = function(rt, pars) {
      parsS <- pars[ , c("muS", "sigmaS", "tauS", "SSD"), drop = FALSE]
      return(dexGaussianS(rt, parsS))
    },
    # Probability function (CDF) for single stop racer
    pfunS = function(rt, pars) {
      parsS <- pars[ , c("muS", "sigmaS", "tauS", "SSD"), drop = FALSE]
      return(pexGaussianS(rt, parsS))
    },
    # Stop probability integral
    sfun = function(pars, n_acc, upper = Inf) {
      return(pstopEXG(pars, n_acc, upper = upper))
    },
    # Random function for SS race
    rfun = function(data = NULL, pars) {
      return(rSSexGaussian(data, pars, ok = attr(pars, "ok")))
    },
    # Race likelihood combining pfun and dfun
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      return(log_likelihood_race_ss(pars, dadm, model, min_ll = min_ll))
    }
  )
}


#' The ex-Gaussian race model of the stop signal task
#'
#' Model file to estimate the ex-Gaussian race model of stop signal task data in EMC2.
#'
#' Model files are almost exclusively used in `design()`.
#'
#' @details
#'
#' Default values are used for all parameters that are not explicitly listed in the `formula`
#' argument of `design()`.They can also be accessed with `SSexG()$p_types`.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default**   | **Mapping**                    | **Interpretation**                                            |
#' |-----------|-----------|---------------|-----------|----------------------------|-----------------------------------------------------------|
#' | *mu*       | log         | \[0, Inf\]     | log(.4)         |                            | Mean of Gaussian component of ex-Gaussian go finish time distribution              |
#' | *sigma*       | log       | \[0, Inf\]        | log(.05)    |                            | Standard deviation of Gaussian component of ex-Gaussian go finish time distribution                                        |
#' | *tau*      | log       | \[0, Inf\]        | log(.1)    |                            | Mean (inverse rate) of exponential component of ex-Gaussian go finish time distribution                                          |
#' | *muS*       | log       | \[0, Inf\]        | log(.3)    |                            | Mean of Gaussian component of ex-Gaussian stop finish time distribution           |
#' | *sigmaS*       | log    | \[0, Inf\]        | log(.025)|                   | Standard deviation of Gaussian component of ex-Gaussian stop finish time distribution                              |
#' | *tauS*      | log    | \[0, Inf\]        | log(.05)  |  | Mean (inverse rate) of exponential component of ex-Gaussian stop finish time distribution       |
#' | *tf*      | probit       | \[0, 1\]        | qnorm(0)    |                            | Attentional lapse rate for stop process ("trigger failure")           |
#' | *gf*     | probit       | \[0, 1\]        | qnorm(0)    |                            | Attentional lapse rate for go process ("go failure")    |
#' | *exg_lb*      | -       | \[-Inf, Inf\]        | .05    |                            | Lower bound of ex-Gaussian go finish time distribution           |
#' | *exgS_lb*     | -       | \[-Inf, Inf\]        | .05    |                            | Lower bound of ex-Gaussian stop finish time distribution    |
#'
#' The independent race model of the stop signal task was introduced in Logan & Cowan (1984).
#' Several review papers provide further background on the task paradigm and modelling framework (e.g., Matzke et al., 2018; Verbruggen et al., 2019; Colonius et al., 2023).
#'
#' The ex-Gaussian model was introduced as a simple race (i.e., unitary go process racing against stop process) in Matzke et al. (2013).
#' This was extended in Matzke et al. (2019) to allow for an arbitrary number of go processes.
#'
#' Attentional lapse rate parameters for the stop process and go process, termed "trigger failure" and "go failure" respectively, were introduced in Matzke et al. (2017) and Matzke et al. (2019).
#'
#' Lower truncation of the ex-Gaussian finish time distributions was introduced in Tanis et al. (2024).
#'
#' Note that the ex-Gaussian parameters do not have unique psychological interpretations (Matzke & Wagenmakers, 2009). Parameter estimates are typically summarised in terms of the mean of the ex-Gaussian distribution, which is given by `mu + tau`.
#'
#' @references
#' Colonius, H., & Diederich, A. (2023). Modeling response inhibition in the stop signal task. *F. G. Ashby, H. Colonius, & E. N. Dzhafarov (Eds.), New Handbook of Mathematical Psychology*, *3*, 311-356. \doi{10.1017/9781108902724.008}
#'
#' Logan, G. D., & Cowan, W. B. (1984). On the ability to inhibit thought and action: A theory of an act of control. *Psychological Review*, *91*(3), 295. \doi{10.1037/0033-295X.91.3.295}
#'
#' Matzke, D., & Wagenmakers, E. J. (2009). Psychological interpretation of the ex-Gaussian and shifted Wald parameters: A diffusion model analysis. *Psychonomic Bulletin & Review*, *16*, 798-817. \doi{10.3758/PBR.16.5.798}
#'
#' Matzke, D., Dolan, C. V., Logan, G. D., Brown, S. D., & Wagenmakers, E. J. (2013). Bayesian parametric estimation of stop-signal reaction time distributions. *Journal of Experimental Psychology: General*, *142*(4), 1047. \doi{10.1037/a0030543}
#'
#' Matzke, D., Love, J., & Heathcote, A. (2017). A Bayesian approach for estimating the probability of trigger failures in the stop-signal paradigm. *Behavior Research Methods*, *49*, 267-281. \doi{10.3758/s13428-015-0695-8}
#'
#' Matzke, D., Verbruggen, F., & Logan, G. (2018). The stop-signal paradigm. *Stevens’ handbook of experimental psychology and cognitive neuroscience*, *5*, 383-427. \doi{10.1002/9781119170174.epcn510}
#'
#' Matzke, D., Curley, S., Gong, C. Q., & Heathcote, A. (2019). Inhibiting responses to difficult choices. *Journal of Experimental Psychology: General*, *148*(1), 124. \doi{10.1037/xge0000525}
#'
#' Tanis, C. C., Heathcote, A., Zrubka, M., & Matzke, D. (2024). A hybrid approach to dynamic cognitive psychometrics: Dynamic cognitive psychometrics. *Behavior Research Methods*, *56*(6), 5647-5666. \doi{10.3758/s13428-023-02295-y}
#'
#' Verbruggen, F., Aron, A. R., Band, G. P., Beste, C., Bissett, P. G., Brockett, A. T., ... & Boehler, C. N. (2019). A consensus guide to capturing the ability to inhibit actions and impulsive behaviors in the stop-signal task. *elife*, *8*, e46323. \doi{10.7554/eLife.46323}
#'
#' @return A model list with all the necessary functions to sample
#' @export
SStexG <- function() {
  list(
    type = "RACE",
    c_name = "SS_TEXG",
    p_types = c(
      mu = log(.4), sigma = log(.05), tau = log(.1),
      muS = log(.3), sigmaS = log(.025), tauS = log(.05),
      tf = qnorm(0), gf = qnorm(0),
      exg_lb = .05, exgS_lb = .05
    ),
    transform = list(
      func = c(
        mu = "exp", sigma = "exp", tau = "exp",
        muS = "exp", sigmaS = "exp", tauS = "exp",
        tf = "pnorm", gf = "pnorm",
        exg_lb = "identity", exgS_lb = "identity"
      )
    ),
    bound = list(
      minmax = cbind(
        mu = c(0, Inf), sigma = c(1e-4, Inf), tau = c(1e-4, Inf),
        muS = c(0, Inf), sigmaS = c(1e-4, Inf), tauS = c(1e-4, Inf),
        tf = c(.001, .999), gf = c(.001, .999),
        exg_lb = c(-Inf, Inf), exgS_lb = c(-Inf, Inf)
      ),
      exception = c(
        tf = 0, gf = 0,
        exg_lb = -Inf, exgS_lb = -Inf
      )
    ),
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      # if (any(names(dadm)=="SSD")) pars <- cbind(pars,SSD=dadm$SSD) else
      #   pars <- cbind(pars,SSD=rep(Inf,dim(pars)[1]))
      pars <- cbind(pars, SSD = dadm$SSD)
      pars <- cbind(pars, lI = as.numeric(dadm$lI))  # Only necessary for data generation.
      return(pars)
    },
    # Density function (PDF) for single go racer
    dfunG = function(rt, pars) return(dtexGaussianG(rt, pars)),
    # Probability function (CDF) for single go racer
    pfunG = function(rt, pars) return(ptexGaussianG(rt, pars)),
    # Density function (PDF) for single stop racer
    dfunS = function(rt, pars) {
      parsS <- pars[ , c("muS", "sigmaS", "tauS", "SSD", "exgS_lb"), drop=FALSE]
      return(dtexGaussianS(rt, parsS))
    },
    # Probability function (CDF) for single stop racer
    pfunS = function(rt, pars) {
      parsS <- pars[ , c("muS", "sigmaS", "tauS", "SSD", "exgS_lb"), drop=FALSE]
      return(ptexGaussianS(rt, parsS))
    },
    # Stop probability integral
    sfun = function(pars, n_acc, upper = Inf) {
      return(pstopTEXG(pars, n_acc, upper = upper))
    },
    # Random function for SS race
    # TODO
    rfun = function(data = NULL, pars) {
      return(rSSexGaussian(data, pars, ok = attr(pars, "ok")))
    },
    # Race likelihood combining pfun and dfun
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      return(log_likelihood_race_ss_old(pars, dadm, model, min_ll = min_ll))
    }
  )
}


#####################  RDEX ----

#### RDEX SS random ----

rSShybrid <- function(data,pars,ok=rep(TRUE,dim(pars)[1]))
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars must contain an SSD column and an lI column indicating if an
  # accumulator is triggered by the stop signal (ST = stop-triggered).

  # NB1: Go failures will only apply to accumulators where lI = TRUE
  #      and can still have a stop-triggered response on a go-failure trial.
{
  lR <- data$lR # For Michelle as an example
  pars[,c("A","B","v")] <- pars[,c("A","B","v")]/pars[ok,"s"]

  nacc <- length(levels(lR))   # Does not include stop runner
  ntrials <- dim(pars)[1]/nacc # Number of trials to simulate
  is1 <- lR==levels(lR)[1]     # First go accumulator
  acc <- 1:nacc                # choice (go and ST) accumulator index

  # stop-triggered racers
  isST <- pars[,"lI"]==1              # Boolean for all pars
  accST <- acc[pars[1:nacc,"lI"]==1]  # Index of ST accumulator, from 1st trial

  # Setup race among nacc+1 (includes stop accumulator)
  # Default Inf finishing time so if not changed always looses
  dt <- matrix(Inf,nrow=nacc+1,ncol=ntrials)

  # Go failure trials (rep for each accumulator)
  isgf <- rep(pars[is1,"gf"] > runif(ntrials),each=nacc)
  # Go accumulators that don't fail (removing ST)
  isGO <- !isgf & !isST
  ngo <- sum(isGO)

  # Fill in go accumulators
  if (any(isGO)) dt[-1,][isGO] <-
    pars[isGO,"t0"] + rWald(ngo,pars[isGO,"B"],pars[isGO,"v"],pars[isGO,"A"])

  # pick out stop trials and races with SSD that is not Inf (i.e., finite or
  # NA, the latter so a staircase can be filled in)
  isStrial <- !is.infinite(pars[is1,"SSD"])
  isSrace <- rep(isStrial,each=nacc)

  # pick out stop trials that are triggered
  isT <- pars[is1,"tf"][isStrial] < runif(sum(isStrial))

  # Logical to pick stop-triggered accumulators that are triggered
  isSTT <- logical(ntrials*nacc) # Default false

  # Pick out stop-triggered accumulators that are triggered
  # NB: can occur on gf trials
  isSTT[isSrace][rep(isT,each=nacc) & isST[isSrace]] <- TRUE
  nst <- sum(isSTT)

  # Fill in stop-triggered accumulators
  if (any(isSTT)) dt[-1,][isSTT] <-
    pars[isSTT,"t0"] + rWald(ngo,pars[isSTT,"B"],pars[isSTT,"v"],pars[isSTT,"A"])

  # pick out triggered stop racers
  isTS <- logical(ntrials)
  isTS[isStrial][isT] <- TRUE
  ns <- sum(isTS)

  # Fill in stop accumulators
  if (any(isTS)) dt[1,isTS] <- rexG(ns,pars[is1,"muS"][isTS],
    pars[is1,"sigmaS"][isTS],pars[is1,"tauS"][isTS])

  # staircase algorithm
  pstair <- is.na(pars[,"SSD"])
  stair <- pstair[is1]
  if (any(stair)) {
    if (is.null(attr(pars,"staircase")))
      stop("When SSD has NAs a staircase list must be supplied!")

    staircase <- attr(pars,"staircase")
    data <- data[is1,]
    if (length(accST)>0)
      staircase$accST <- 1+accST
    if (is.null(attr(staircase,"staircase_function")))
      attr(staircase,"staircase_function") <- staircase_function

    allR <- allrt <- numeric(ncol(dt))  # to store unified results
    dts <- dt[,stair,drop=F]

    # Non-staircase trials
    dt <- dt[,!stair,drop=F]
    spars <- pars[pstair,,drop=F]
    pars <- pars[!pstair,,drop=F]
    isTS <- isTS[!stair]
    isSTT <- isSTT[!pstair]
    isST <- isST[!pstair]
    ntrials <- sum(!stair)
    is1 <- is1[!pstair]
    data <- data[stair,]

    stair_res <- attr(staircase,"staircase_function")(dts,staircase)
    allR[stair] <- stair_res$sR
    allrt[stair] <- stair_res$srt
  }


  # All SSD already filled in so return R and rt

  # Add SSD to triggered stop accumulators
  if (any(isTS)) dt[1,isTS] <- dt[1,isTS] + pars[is1,"SSD"][isTS]

  # Add SSD to stop-triggered accumulators that are triggered
  if (any(isSTT)) dt[-1,][isSTT] <- dt[-1,][isSTT] + pars[isSTT,"SSD"]

  # R <- factor(rep(NA,ntrials),levels=levels(lR))
  R <- rt <- rep(NA,ntrials)

  # All accumulators Inf (usually when both gf and tf)
  allinf <- apply(dt,2,function(x)all(is.infinite(x)))

  # get winner of stop and go where there is a race
  r <- c(0, acc)[apply(dt[,!allinf,drop=FALSE],2,which.min)]

  # stop wins
  stopwins <- r==0

  # First fill in cases where stop looses, can include wins by ST
  if (any(!stopwins)) {
    rgo <- r[!stopwins]
    R[!allinf][!stopwins] <- rgo
    pick <- cbind(rgo,c(1:sum(!stopwins))) # Matrix to pick winner
    rt[!allinf][!stopwins] <-
      dt[-1,!allinf,drop=FALSE][,!stopwins,drop=FALSE][pick]
  }

  # then if stop wins and kills go, find ST winner
  if (any(isST) & any(stopwins)) {
    # stop triggered accumulators that are racing
    rst <- dt[-1,!allinf,drop=FALSE][accST,stopwins,drop=FALSE]
    # stop-triggered winners
    rtw <- apply(rst,2,which.min)
    # index for stop-triggered
    R[!allinf][stopwins] <- accST[rtw]
    pick <- cbind(rtw,1:ncol(rst))
    rt[!allinf][stopwins] <- rst[pick]
  }

  rt[is.na(R)] <- NA

  if (any(stair)) {
    allrt[!stair] <- rt
    allR[!stair] <- R
    allSSD <- NA
    allSSD[stair] <- stair_res$SSD
    out <- cbind.data.frame(R=factor(allR,levels=1:nacc,labels=levels(lR)),
                            rt=allrt, SSD = allSSD)
    return(out)
  }
  cbind.data.frame(R=factor(R,levels=1:nacc,labels=levels(lR)),rt=rt)
}

#### RDEX stop probability ----

# # NB: these functions are in Rcpp
# dWald_RDEX
# pWald_RDEX
# EMC2:::stopfn_rdex


pstopHybrid <- function(parstop,n_acc,upper=Inf,
  gpars=c("v","B","t0","A"),spars=c("muS","sigmaS","tauS"))
{
  sindex <- seq(1,nrow(parstop),by=n_acc)  # Stop accumulator index
  ps <- parstop[sindex,spars,drop=FALSE]   # Stop accumulator parameters
  SSDs <- parstop[sindex,"SSD",drop=FALSE] # SSDs
  ntrials <- length(SSDs)
  if (length(upper)==1) upper <- rep(upper,length.out=ntrials)
  pgo <- array(parstop[,gpars],dim=c(n_acc,ntrials,length(gpars)),
               dimnames=list(NULL,NULL,gpars))
  cells <- apply(
    cbind(SSDs,ps,upper,matrix(as.vector(aperm(pgo,c(2,1,3))),nrow=ntrials))
  ,1,paste,collapse="")
  uniq <- !duplicated(cells)
  ups <- sapply(1:sum(uniq),function(i){
    my.integrate_old(f=stopfn_rdex,lower=0,upper=upper[i],
      mu=ps[i,"muS"],sigma=ps[i,"sigmaS"],tau=ps[i,"tauS"],
      v=pgo[,i,"v"],B=pgo[,i,"B"],A=pgo[,i,"A"],t0=pgo[,i,"t0"],
      SSD=SSDs[i],n_acc=n_acc)
  })
  ups[as.numeric(factor(cells,levels=cells[uniq]))]
}


#### RDEX model list ----

#' Stop-signal Hybrid (RDM go, ExGaussian stop) race
#'
#' Model file to estimate the Hybrid race model for Stop-Signal data
#'
#' Model files are almost exclusively used in `design()`.
#'
#' @details
#'
#' Default values are used for all parameters that are not explicitly listed in the `formula`
#' argument of `design()`.They can also be accessed with `SShybrid()$p_types`.
#'
#' @return A model list with all the necessary functions to sample
#' @export
SShybrid <- function() {
  list(
    type="RACE",
    p_types=c("v" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0),"s" = log(1),
              muS=log(.3),sigmaS=log(.025),tauS=log(.05),tf=qnorm(0),gf=qnorm(0)),
    transform=list(func=c(v = "exp", B = "exp", A = "exp",t0 = "exp", s = "exp",
                          muS = "exp", sigmaS = "exp", tauS = "exp",
                          tf="pnorm",gf="pnorm")),
    bound=list(minmax=cbind(v=c(1e-3,Inf), B=c(0,Inf), A=c(1e-4,Inf),t0=c(0.05,Inf),
                            s=c(0,Inf),muS=c(0,Inf), sigmaS=c(1e-4,Inf), tauS=c(1e-4,Inf),
                            tf=c(.001,.999),gf=c(.001,.999)),
               exception=c(A=0, v=0,tf=0,gf=0)),
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars <- cbind(pars,b=pars[,"B"] + pars[,"A"])
      pars <- cbind(pars,SSD=dadm$SSD)
      pars <- cbind(pars,lI=as.numeric(dadm$lI))  # Only necessary for data generation.
      pars
    },
    # Density function (PDF) for single go racer
    dfunG=function(rt,pars) dRDM(rt,pars),
    # Probability function (CDF) for single go racer
    pfunG=function(rt,pars) pRDM(rt,pars),
    # Density function (PDF) for single stop racer
    dfunS=function(rt,pars) dexGaussianS(rt,
      pars[,c("muS","sigmaS","tauS","SSD"),drop=FALSE]),
    # Probability function (CDF) for single stop racer
    pfunS=function(rt,pars) pexGaussianS(rt,
      pars[,c("muS","sigmaS","tauS","SSD"),drop=FALSE]),
    # Stop probability integral
    sfun=function(pars,n_acc,upper=Inf) pstopHybrid(pars,n_acc,upper=upper),
    # Random function for SS race
    rfun=function(data=NULL,pars) {
      rSShybrid(data,pars,ok=attr(pars, "ok"))
    },
    # Race likelihood combining pfun and dfun
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10))
      log_likelihood_race_ss(pars, dadm, model, min_ll = min_ll)
  )
}






