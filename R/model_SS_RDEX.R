

dRDEXG <- function(rt,pars)
{
  if (any(dimnames(pars)[[2]]=="s")) # rescale
     pars[,c("A","B","v")] <- pars[,c("A","B","v")]/pars[,"s"]
  out <- numeric(length(rt))
  ok <- !is.na(rt)
  out[ok] <- dRDM(rt[ok],pars[ok,,drop=FALSE])
  out
}


pRDEXG <- function(rt,pars)
{
  if (any(dimnames(pars)[[2]]=="s")) # rescale
     pars[,c("A","B","v")] <- pars[,c("A","B","v")]/pars[,"s"]
  out <- numeric(length(rt))
  ok <- !is.na(rt)
  out[ok] <- pRDM(rt[ok],pars[ok,,drop=FALSE])
  out
}

# NB Uses update_ssd declared in model_EXGSS

rRDEX <- function(lR,pars,staircase0,stairstep,stairmin,stairmax,
                  p_types=c("v","B","A","t0","muS","sigmaS","tauS","tf","gf"))
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
  if (any(dimnames(pars)[[2]]=="s")) # rescale
     pars[isgo,c("A","B","v")] <- pars[isgo,c("A","B","v")]/pars[isgo,"s"]
  ngo <- sum(isgo)
  dt[-1,] <- pars[isgo,"t0"] + rWald(ngo,
      B=pars[isgo,"B"],v=pars[isgo,"v"],A=pars[isgo,"A"])
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
      r <- apply(dt[,isfixed],2,which.min)
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

# NB: Uses pEXG, dEXG, dEXGrace and pstopEXG and stopfn_exg from model_EXGSS

# NB: Uses pWald, pigt and pigt0 from model_RDM and dEXG from model_EXGSS

# Moved to C++ in model_SS_RDEX

# dRDEXrace <- function(dt,mu,sigma,tau,v,B,A,t0)
#     # Generates defective PDF for win by first runner, dt (decison time) is
#     # a matrix with length(mu) rows, one row for each runner, and one column
#     # for each decision time for which a defective density value will be
#     # returned.
# {
#     dt[1,] <- dEXG(dt[1,],mu,sigma,tau)
#     for (i in 1:length(v))
#       dt[1,] <- dt[1,]*(1-pWald(dt[i+1,],v[i],B[i],A[i],t0[i]))
#     dt[1,]
# }
#
# stopfn_rdex <- function(t,n_acc,mu,sigma,tau,v,B,A,t0,SSD)
#     # Used by my.integrate, t = vector of times, SSD is a scalar stop-signal delay.
# {
#     dt <- matrix(rep(t+SSD,each=n_acc),nrow=n_acc)
#     dt[1,] <- dt[1,]-SSD
#     dRDEXrace(dt,mu,sigma,tau,v,B,A,t0)
# }


pstopRDEX <- function(parstop,n_acc,
  gpars=c("v","B","A","t0"),spars=c("muS","sigmaS","tauS"))
{
  if (any(dimnames(parstop)[[2]]=="s")) # rescale
     parstop[,c("A","B","v")] <- parstop[,c("A","B","v")]/parstop[,"s"]
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
  ups <- sapply(1:sum(uniq),function(i){
  	my.integrate(f=stopfn_rdex,lower=0,
      n_acc=n_acc,SSD=SSDs[uniq][i],
      mu=ps[uniq,"muS"][i],
      sigma=ps[uniq,"sigmaS"][i],
      tau=ps[uniq,"tauS"][i],
      v=pgo[,uniq,"v",drop=FALSE][,i,1],
      B=pgo[,uniq,"B",drop=FALSE][,i,1],
      A=pgo[,uniq,"A",drop=FALSE][,i,1],
      t0=pgo[,uniq,"t0",drop=FALSE][,i,1])
  })
  ups[as.numeric(factor(cells,levels=cells[uniq]))]
}


#### Model list ----
SSrdexB <- function() {
    list(
    type="RACE",
    p_types=c("v","B","A","t0","s","muS","sigmaS","tauS","tf","gf"),
    # Transform to natural scale
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
      attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        (pars[,"tauS"] < 1) & (pars[,"sigmaS"] < 1) &
        ((pars[,"tf"] > 1e-6) | pars[,"tf"] == 0) & ((pars[,"gf"] > 1e-6) | pars[,"gf"] == 0)
      pars
    },
    # Density function (PDF) for single go racer
    dfunG=function(rt,pars) dRDEXG(rt,pars),
    # Probability function (CDF) for single go racer
    pfunG=function(rt,pars) pRDEXG(rt,pars),
    # Density function (PDF) for single stop racer
    dfunS=function(rt,pars) dexGaussianS(rt,pars),
    # Probability function (CDF) for single stop racer
    pfunS=function(rt,pars) pexGaussianS(rt,pars),
    # Stop probability integral
    sfun=function(pars,n_acc) pstopRDEX(pars,n_acc),
    # Random function for SS race
    rfun=function(lR,pars) rRDEX(lR,pars,
      staircase0=.2,stairstep=.05,stairmin=0,stairmax=Inf),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
      log_likelihood_race_ss(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
  )
}
