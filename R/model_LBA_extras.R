#### Multiple threshold models ----


rLBA_BE <- function(pars,ok=rep(TRUE,dim(pars)[1]),return_activation = FALSE)
  # Balance of evidence
  # Rating from 1=high confidence left to length(d1)+length(d2)+2 for
  # high confidence right, assumes posdrift=TRUE.
  # pars is a matrix with columns A, b, t0, v, sv and r sorted so parameters for
  # each trial are in contiguous rows. Returns nrow(pars)/2 trials
  # d1 and d2 are matrices of lower thresholds for first and second accumulator
  # with one row per trial.
{

  getq <- function(pmat)
    # sample from mvtnorm
  {
    sv2 <- pmat[,"sv"]^2
    sigma <- matrix(c(sv2[1], rep(prod(pmat[,"sv"]) * pmat[1,"r"], 2), sv2[2]),2, 2, byrow = TRUE)
    tmvtnorm::rtmvnorm(n = 1, mean = pmat[,"v"], sigma = sigma, lower = rep(0, 2))
  }

  get_rate <- function(x){.bincode(x[1],x[-1])}

  bad <- rep(NA, dim(pars)[1]/2)
  out <- data.frame(R = bad, rt = bad)
  pars <- pars[ok,]

  pmats <- array(pars,dim=c(2,nrow(pars)/2,ncol(pars)),
                dimnames=list(NULL,NULL,colnames(pars)))
  n <- dim(pmats)[2] # number of trials
  q <- apply(pmats,2,getq)
  A <- matrix(runif(length(pmats[,,"A"]),min=0,max=pmats[,,"A"]),nrow=2)
  rts <- (pmats[,,"b"]-A) / q
  rts <- rts + pmats[,,"t0"]
  response <- apply(rts, 2, which.min)
  rt <- rts[cbind(response,seq_len(n))]
  # rating
  looser2 <- response==1
  lmat <- cbind(1 + as.numeric(looser2),seq_len(n))
  ylooser <- A[lmat] + q[lmat]*(rt-pmats[1,,"t0"])
  rating <- numeric(length(response))
  # 2=nacc x n x #DT
  d <- pmats[,,substr(dimnames(pmats)[[3]],1,2)=="DT",drop=FALSE]
  nr <- dim(d)[3]+1 # number of ratings = DT1, DT2 ..., b
  if (any(looser2))
    rating[looser2] <- apply(cbind(ylooser[looser2],rep(0,sum(looser2)),d[2,looser2,],pmats[2,looser2,"b"]),1,get_rate)
  if (any(!looser2))
    rating[!looser2] <- (2*nr+1)-apply(cbind(ylooser[!looser2],rep(0,sum(!looser2)),d[1,!looser2,],pmats[2,!looser2,"b"]),1,get_rate)
  ok <- matrix(ok,nrow=2)[1,]
  out[ok,] <- cbind(response = rating,rt = rt)

  if (return_activation) {
    tmp <- cbind(out[ok,],y = t(A) + t(q) * (rt - pmats[1,,"t0"]))
    out <- cbind(out,y1=bad,y2=bad)
    out[ok,] <- tmp
  }
  out$R <- factor(out$R,levels=1:(2*nr))
  out
}

rLBA_TC <- function(pars,ok=rep(TRUE,dim(pars)[1]),return_activation = FALSE)
  # Threshold count
  # Rating from 1=high confidence left to length(d1)+length(d2)+2 for
  # high confidence right, assumes posdrift=TRUE.
  # pars is a matrix with columns A, b, t0, v, sv and r sorted so parameters for
  # each trial are in contiguous rows. Returns nrow(pars)/2 trials
  # d1 and d2 are matrices of lower thresholds for first and second accumulator
  # with one row per trial.
{

  getq <- function(pmat)
    # sample from mvtnorm
  {
    sv2 <- pmat[,"sv"]^2
    sigma <- matrix(c(sv2[1], rep(prod(pmat[,"sv"]) * pmat[1,"r"], 2), sv2[2]),2, 2, byrow = TRUE)
    tmvtnorm::rtmvnorm(n = 1, mean = pmat[,"v"], sigma = sigma, lower = rep(0, 2))
  }

  bad <- rep(NA, dim(pars)[1]/2)
  out <- data.frame(R = bad, rt = bad)
  pars <- pars[ok,]

  Dnams <- dimnames(pars)[[2]][substr(dimnames(pars)[[2]],1,2)=="DT"]
  tmats <- array(cbind(pars[,c(Dnams,"b")]),dim=c(2,nrow(pars)/2,1+length(Dnams)),
                dimnames=list(NULL,NULL,c(Dnams,"b")))
  pnams <- dimnames(pars)[[2]][!(dimnames(pars)[[2]] %in% c(Dnams,"b"))]
  pmats <- array(pars[,pnams],dim=c(2,nrow(pars)/2,length(pnams)),
                dimnames=list(NULL,NULL,pnams))

  nr <- dim(tmats)[3] # number of ratings = triggering count
  nt=nr*2

  # rating lookup
  dL <- matrix(nrow=nt,ncol=2)
  even <- array(c(1:nt),dim=c(2,nr))
  odd <- as.vector(even[1,])
  even <- as.vector(even[2,])
  dL[odd,1] <- 1; dL[even,1] <- 2
  dL[odd,2] <- nr:1; dL[even,2] <- 1:nr
  rate <- dL[,1] + (dL[,2]-1)*2

  n <- dim(pmats)[2]        # number of trials
  q <- apply(pmats,2,getq)  # rates
  A <- matrix(runif(length(pmats[,,"A"]),min=0,max=pmats[,,"A"]),nrow=2)
  rts=array(apply(tmats,3,function(x){(x-A)/q}),dim=dim(tmats))

  rts[rts<0] <- 0

  ok <- matrix(ok,nrow=2)[1,]
  out[ok,] <- t(apply(rts,2,function(x){
    rt <- sort(x)[nr]
    c(rate[x==rt],rt)
  }))
  out[ok,"rt"] <- out[ok,"rt"] + pmats[1,,"t0"]
  if (return_activation) {
    tmp <- cbind(out[ok,],y = t(A) + t(q) * (out[ok,"rt"] - pmats[1,,"t0"]))
    out <- cbind(out,y1=bad,y2=bad)
    out[ok,] <- tmp
  }
  out$R <- factor(out$R,levels=1:(2*nr))
  out
}


n1PDF_MTR_1 <- function(rt, pars,dl,du,b)
  # Works for Balance of Evidence and Threshold count.
  # Assumes single rt value, posdrift = TRUE
  # pars is a 2 row matrix of parameters with columns A, t0, v, sv, r
  # with first row corresponding to the accumulator that stops the race
  # b is stopping threshold, dl and du are the upper and lower thresholds for
  # the accumulator that does not stop the race.
  # r and t0 same for both (second ignored)
{

  if (rt <= pars[1,"t0"]) {
    return(0)
  }
  n <- 2
  t <- rt - pars[1,"t0"]
  sv2 <- pars[,"sv"]^2
  sigma <- matrix(c(sv2[1], rep(prod(pars[,"sv"]) * pars[1,"r"], 2), sv2[2]),2, 2, byrow = TRUE)

  Pqpos <- MomTrunc:::recintab0(kappa = rep(0, n),
                                  mu = pars[,"v"],
                                  S = sigma,
                                  a = rep(0, n),
                                  b = rep(Inf, n))

  # Integral bounds
  lb <- ub <- numeric(4)
  lb[1] <- 0
  lb[2] <- ifelse((du - pars[2,"A"]) / t <= 0, - du / pars[2,"A"], -1)
  lb[3] <- 0
  lb[4] <- ifelse((dl - pars[2,"A"]) / t <= 0, - dl / pars[2,"A"], -1)

  ub[1] <- ifelse((du - pars[2,"A"]) / t <= 0, 0, (du - pars[2,"A"]) / t)
  ub[2] <- ifelse(du / t          <= 0, - du / pars[2,"A"], 0)
  ub[3] <- ifelse((dl - pars[2,"A"]) / t <= 0, 0, (dl - pars[2,"A"]) / t)
  ub[4] <- ifelse(dl / t       <= 0, - dl / pars[2,"A"], 0)

  # signs
  signs <- c(1, -1, -1, 1)

  # Cs
  cs <- cbind(NA, c(1, t / pars[2,"A"]), NA, c(1, t / pars[2,"A"]))

  # ms
  ms <- cbind(NA, c(0, - du / pars[2,"A"]), NA, c(0, - dl / pars[2,"A"]))

  # index for which integrals relevant
  index <- which(lb != ub)

  ints <- sapply(index, function(x) {

    if (x == 1 || x == 3) {

      mom <- MomTrunc:::recintab0(kappa = c(1, 0),
                                  mu = pars[,"v"],
                                  S = sigma,
                                  a = c((b - pars[1,"A"]) / t, lb[x]),
                                  b = c( b / t,         ub[x]))
      return(signs[x] * mom[length(mom)])

    } else if (x == 2 || x == 4) {

      C_j <- diag(cs[,x])
      m_j <- ms[,x]

      m_x <- as.vector(m_j + C_j %*% pars[,"v"])
      sigma_x <- C_j %*% sigma %*% t(C_j)

      mom <- MomTrunc:::recintab0(kappa = rep(1, 2),
                                  mu = m_x,
                                  S = sigma_x,
                                  a = c((b - pars[1,"A"]) / t, lb[x]),
                                  b = c( b / t,         ub[x]))
      return(signs[x] * mom[length(mom)])

    }

  })

  out <- sum(ints) / (Pqpos * pars[1,"A"])

  return(pmax(out,0))

}


#' Mlba_B
#'
#' LBA model accommodating missing values (truncation and censoring) and
#' assuming positive rates (i.e., no intrinsic omissions)
#'
#' @return A model list with all the necessary functions to sample
#' @export
MlbaB <- function(){
  list(
    type="RACE",
    p_types=c("v","sv","B","A","t0","pContaminant"),
    Ntransform=function(x) {
      # Transform to natural scale
      doexp <- !(dimnames(x)[[2]] %in% c("v","pContaminant"))
      x[,doexp] <- exp(x[,doexp])
      doprobit <- dimnames(x)[[2]] == "pContaminant"
      x[,doprobit] <- pnorm(x[,doprobit])
      x
    },
    # p_vector transform
    transform = function(p) p,
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars <- cbind(pars,b=pars[,"B"] + pars[,"A"])
      attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | (pars[,"A"] == 0)) &
        (abs(pars[,"v"])<100) & ((pars[,"sv"] > 1e-3) | pars[,"sv"] == 0)
      pars
    },
    # # Random function for racing accumulator
    # Random function for racing accumulator
    rfun=function(lR=NULL,pars) {
      ok <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        (abs(pars[,"v"])<100) & ((pars[,"sv"] > 1e-3) | pars[,"sv"] == 0)
      if (is.null(lR)) ok else rLBA(lR,pars,posdrift=TRUE)
    },
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dLBA(rt,pars,posdrift = TRUE, robust = FALSE),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pLBA(rt,pars,posdrift = TRUE, robust = FALSE),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm){
      log_likelihood_race_missing(p_vector=p_vector, dadm = dadm, min_ll=log(1e-10))
    }
  )
}

#' MUlbaB
#'
#' LBA model accommodating missing values (truncation and censoring) and
#' assuming unbounded rates (i.e., allows intrinsic omissions)
#'
#' @return A model list with all the necessary functions to sample
#' @export
MUlbaB <- function(){
  list(
    type="RACE",
    p_types=c("v","sv","B","A","t0","pContaminant"),
    Ntransform=function(x) {
      # Transform to natural scale
      doexp <- !(dimnames(x)[[2]] %in% c("v","pContaminant"))
      x[,doexp] <- exp(x[,doexp])
      doprobit <- dimnames(x)[[2]] == "pContaminant"
      x[,doprobit] <- pnorm(x[,doprobit])
      x
    },
    # p_vector transform
    transform = function(p) p,
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars <- cbind(pars,b=pars[,"B"] + pars[,"A"])
      attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | (pars[,"A"] == 0)) &
        (abs(pars[,"v"])<100) & ((pars[,"sv"] > 1e-3) | pars[,"sv"] == 0)
      pars
    },
    # # Random function for racing accumulator
    # Random function for racing accumulator
    rfun=function(lR=NULL,pars) {
      ok <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        (abs(pars[,"v"])<100) & ((pars[,"sv"] > 1e-3) | pars[,"sv"] == 0)
      if (is.null(lR)) ok else rLBA(lR,pars,posdrift=FALSE)
    },
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dLBA(rt,pars,posdrift = FALSE, robust = FALSE),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pLBA(rt,pars,posdrift = FALSE, robust = FALSE),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm){
      log_likelihood_race_missing_LBAU(p_vector=p_vector, dadm = dadm, min_ll=log(1e-10))
    }
  )
}


#' Balance of Evidence 2 Threshold LBA model
#'
#' All thresholds greater than A, DT1 <- A+DT1, b <- DT1 + B
#'
#' @return A model list with all the necessary functions to sample
#' @export
BE2lbaB <- function(){
  list(
    type="MT",
    # p_vector transform, sets sv as a scaling parameter
    p_types=c("v","sv","r","A","DT1","B","t0"),
    transform = function(p) p,
    # Transform to natural scale
    Ntransform=function(x) {
      x[,dimnames(x)[[2]]=="r"] <- 2*pnorm(x[,dimnames(x)[[2]]=="r"])-1
      x[,!(dimnames(x)[[2]] %in% c("v","r"))] <- exp(x[,!(dimnames(x)[[2]] %in% c("v","r"))])
      x
    },
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pnams <- c("A",dimnames(pars)[[2]][substr(dimnames(pars)[[2]],1,2)=="DT"],"B")
      pars[,pnams] <- t(apply(pars[,pnams],1,cumsum))
      dimnames(pars)[[2]][dimnames(pars)[[2]]=="B"] <- "b"
      attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        ( abs(pars[,"v"])< 10) & ((pars[,"sv"] > 1e-3) | pars[,"sv"] == 0) & (abs(pars[,"b"])<10)
      pars
    },
    # Random function for racing accumulator
    rfun=function(lR=NULL,pars) {
      ok <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        ( abs(pars[,"v"])< 10) & ((pars[,"sv"] > 1e-3) | pars[,"sv"] == 0) & (abs(pars[,"b"])<10)
      if (is.null(lR)) ok else rLBA_BE(pars,ok=ok)
    },
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dLBA(rt,pars,posdrift = TRUE, robust = FALSE),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pLBA(rt,pars,posdrift = TRUE, robust = FALSE),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)){
      log_likelihood_mt(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
    }
  )
}


#' Balance of Evidence 2 Threshold LBA model
#'
#' Intermediate thresholds proportional to b, can be less than A
#'
#' @return A model list with all the necessary functions to sample
#' @export
BE2lbaBP <- function(){
  list(
    type="MT",
    # p_vector transform, sets sv as a scaling parameter
    p_types=c("v","sv","r","A","DT1","B","t0"),
    transform = function(p) p,
    # Transform to natural scale
    Ntransform=function(x) {
      x[,dimnames(x)[[2]] %in% c("A","B","t0","sv")] <-
        exp(x[,dimnames(x)[[2]] %in% c("A","B","t0","sv")])
      isDT <- substr(dimnames(x)[[2]],1,2)=="DT"
      x[,isDT] <- pnorm(x[,isDT])
      x[,dimnames(x)[[2]]=="r"] <- 2*pnorm(x[,dimnames(x)[[2]]=="r"])-1
      x
    },
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars <- cbind(pars,b=pars[,"B"] + pars[,"A"])
      isDT <- substr(dimnames(pars)[[2]],1,2)=="DT"
      pars[,isDT] <- pars[,isDT]*pars[,"b"]
      attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        ( abs(pars[,"v"])< 10) & ((pars[,"sv"] > 1e-3) | pars[,"sv"] == 0) & (abs(pars[,"b"])<10)
      pars
    },
    # Random function for racing accumulator
    rfun=function(lR=NULL,pars) {
      ok <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        ( abs(pars[,"v"])< 10) & ((pars[,"sv"] > 1e-3) | pars[,"sv"] == 0) & (abs(pars[,"b"])<10)
      if (is.null(lR)) ok else rLBA_BE(pars,ok=ok)
    },
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dLBA(rt,pars,posdrift = TRUE, robust = FALSE),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pLBA(rt,pars,posdrift = TRUE, robust = FALSE),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)){
      log_likelihood_mt(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
    }
  )
}


#' Threshold Count 2 LBA model
#'
#' All thresholds greater than A, DT1 <- A+DT1, b <- DT1 + B
#'
#' @return A model list with all the necessary functions to sample
#' @export
TC2lbaB <- function(){
  list(
    type="TC",
    # p_vector transform, sets sv as a scaling parameter
    p_types=c("v","sv","r","A","DT1","B","t0"),
    transform = function(p) p,
    # Transform to natural scale
    Ntransform=function(x) {
      x[,dimnames(x)[[2]]=="r"] <- 2*pnorm(x[,dimnames(x)[[2]]=="r"])-1
      x[,!(dimnames(x)[[2]] %in% c("v","r"))] <- exp(x[,!(dimnames(x)[[2]] %in% c("v","r"))])
      x
    },
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pnams <- c("A",dimnames(pars)[[2]][substr(dimnames(pars)[[2]],1,2)=="DT"],"B")
      pars[,pnams] <- t(apply(pars[,pnams],1,cumsum))
      dimnames(pars)[[2]][dimnames(pars)[[2]]=="B"] <- "b"
      attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        ( abs(pars[,"v"])< 10) & ((pars[,"sv"] > 1e-3) | pars[,"sv"] == 0) & (abs(pars[,"b"])<10)
      pars
    },
    # Random function for racing accumulator
    rfun=function(lR=NULL,pars) {
      ok <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        ( abs(pars[,"v"])< 10) & ((pars[,"sv"] > 1e-3) | pars[,"sv"] == 0) & (abs(pars[,"b"])<10)
      if (is.null(lR)) ok else rLBA_TC(pars,ok=ok)
    },
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dLBA(rt,pars,posdrift = TRUE, robust = FALSE),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pLBA(rt,pars,posdrift = TRUE, robust = FALSE),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)){
      log_likelihood_mt(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
    }
  )
}


#' Balance of Evidence 3 Threshold LBA model
#'
#' All thresholds greater than A, DT1 <- A+DT1, DT2 <- DT1 + DT2, b <- DT2 + B
#'
#' @return A model list with all the necessary functions to sample
#' @export
BE3lbaB <- function(){
  list(
    type="MT",
    # p_vector transform, sets sv as a scaling parameter
    p_types=c("v","sv","r","A","DT1","DT2","B","t0"),
    transform = function(p) p,
    # Transform to natural scale
    Ntransform=function(x) {
      x[,dimnames(x)[[2]]=="r"] <- 2*pnorm(x[,dimnames(x)[[2]]=="r"])-1
      x[,!(dimnames(x)[[2]] %in% c("v","r"))] <- exp(x[,!(dimnames(x)[[2]] %in% c("v","r"))])
      x
    },
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pnams <- c("A",dimnames(pars)[[2]][substr(dimnames(pars)[[2]],1,2)=="DT"],"B")
      pars[,pnams] <- t(apply(pars[,pnams],1,cumsum))
      dimnames(pars)[[2]][dimnames(pars)[[2]]=="B"] <- "b"
      attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        ( abs(pars[,"v"])< 10) & ((pars[,"sv"] > 1e-3) | pars[,"sv"] == 0) & (abs(pars[,"b"])<10)
      pars
    },
    # Random function for racing accumulator
    rfun=function(lR=NULL,pars) {
      ok <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        ( abs(pars[,"v"])< 10) & ((pars[,"sv"] > 1e-3) | pars[,"sv"] == 0) & (abs(pars[,"b"])<10)
      if (is.null(lR)) ok else rLBA_BE(pars,ok=ok)
    },
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dLBA(rt,pars,posdrift = TRUE, robust = FALSE),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pLBA(rt,pars,posdrift = TRUE, robust = FALSE),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)){
      log_likelihood_mt(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
    }
  )
}

#' Balance of Evidence 3 Threshold LBA model
#'
#' Intermediate thresholds proportional to b, can be less than A
#'
#' @return A model list with all the necessary functions to sample
#' @export
BE3lbaBP <- function(){
  list(
    type="MT",
    # p_vector transform, sets sv as a scaling parameter
    p_types=c("v","sv","r","A","DT1","DT2","B","t0"),
    transform = function(p) p,
    # Transform to natural scale
    Ntransform=function(x) {
      x[,dimnames(x)[[2]] %in% c("A","B","t0","sv")] <-
        exp(x[,dimnames(x)[[2]] %in% c("A","B","t0","sv")])
      isDT <- substr(dimnames(x)[[2]],1,2)=="DT"
      x[,isDT] <- pnorm(x[,isDT])
      x[,dimnames(x)[[2]]=="r"] <- 2*pnorm(x[,dimnames(x)[[2]]=="r"])-1
      x
    },
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars <- cbind(pars,b=pars[,"B"] + pars[,"A"])
      isDT <- substr(dimnames(pars)[[2]],1,2)=="DT"
      pars[,isDT] <- pars[,isDT]*pars[,"b"]
      attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        ( abs(pars[,"v"])< 10) & ((pars[,"sv"] > 1e-3) | pars[,"sv"] == 0) & (abs(pars[,"b"])<10)
      pars
    },
    # Random function for racing accumulator
    rfun=function(lR=NULL,pars) {
      ok <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        ( abs(pars[,"v"])< 10) & ((pars[,"sv"] > 1e-3) | pars[,"sv"] == 0) & (abs(pars[,"b"])<10)
      if (is.null(lR)) ok else rLBA_BE(pars,ok=ok)
    },
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dLBA(rt,pars,posdrift = TRUE, robust = FALSE),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pLBA(rt,pars,posdrift = TRUE, robust = FALSE),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)){
      log_likelihood_mt(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
    }
  )
}

#' Threshold Count 3 LBA model
#'
#' All thresholds greater than A, DT1 <- A+DT1, DT2 <- DT1 + DT2, b <- DT2 + B
#'
#' @return A model list with all the necessary functions to sample
#' @export
TC3lbaB <- function(){
  list(
    type="TC",
    # p_vector transform, sets sv as a scaling parameter
    p_types=c("v","sv","r","A","DT1","DT2","B","t0"),
    transform = function(p) p,
    # Transform to natural scale
    Ntransform=function(x) {
      x[,dimnames(x)[[2]]=="r"] <- 2*pnorm(x[,dimnames(x)[[2]]=="r"])-1
      x[,!(dimnames(x)[[2]] %in% c("v","r"))] <- exp(x[,!(dimnames(x)[[2]] %in% c("v","r"))])
      x
    },
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pnams <- c("A",dimnames(pars)[[2]][substr(dimnames(pars)[[2]],1,2)=="DT"],"B")
      pars[,pnams] <- t(apply(pars[,pnams],1,cumsum))
      dimnames(pars)[[2]][dimnames(pars)[[2]]=="B"] <- "b"
      attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        ( abs(pars[,"v"])< 10) & ((pars[,"sv"] > 1e-3) | pars[,"sv"] == 0) & (abs(pars[,"b"])<10)
      pars
    },
    # Random function for racing accumulator
    rfun=function(lR=NULL,pars) {
      ok <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        ( abs(pars[,"v"])< 10) & ((pars[,"sv"] > 1e-3) | pars[,"sv"] == 0) & (abs(pars[,"b"])<10)
      if (is.null(lR)) ok else rLBA_TC(pars,ok=ok)
    },
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dLBA(rt,pars,posdrift = TRUE, robust = FALSE),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pLBA(rt,pars,posdrift = TRUE, robust = FALSE),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)){
      log_likelihood_mt(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
    }
  )
}


#' Balance of Evidence 4 Threshold LBA model
#'
#' All thresholds greater than A, DT1 <- A+DT1, DT2 <- DT1 + DT2, DT3 <- DT2 + DT3, b <- DT3 + B
#'
#' @return A model list with all the necessary functions to sample
#' @export
BE4lbaB <- function(){
  list(
    type="MT",
    # p_vector transform, sets sv as a scaling parameter
    p_types=c("v","sv","r","A","DT1","DT2","DT3","B","t0"),
    transform = function(p) p,
    # Transform to natural scale
    Ntransform=function(x) {
      x[,dimnames(x)[[2]]=="r"] <- 2*pnorm(x[,dimnames(x)[[2]]=="r"])-1
      x[,!(dimnames(x)[[2]] %in% c("v","r"))] <- exp(x[,!(dimnames(x)[[2]] %in% c("v","r"))])
      x
    },
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pnams <- c("A",dimnames(pars)[[2]][substr(dimnames(pars)[[2]],1,2)=="DT"],"B")
      pars[,pnams] <- t(apply(pars[,pnams],1,cumsum))
      dimnames(pars)[[2]][dimnames(pars)[[2]]=="B"] <- "b"
      attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        ( abs(pars[,"v"])< 10) & ((pars[,"sv"] > 1e-3) | pars[,"sv"] == 0) & (abs(pars[,"b"])<10)
      pars
    },
    # Random function for racing accumulator
    rfun=function(lR=NULL,pars) {
      ok <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        ( abs(pars[,"v"])< 10) & ((pars[,"sv"] > 1e-3) | pars[,"sv"] == 0) & (abs(pars[,"b"])<10)
      if (is.null(lR)) ok else rLBA_BE(pars,ok=ok)
    },
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dLBA(rt,pars,posdrift = TRUE, robust = FALSE),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pLBA(rt,pars,posdrift = TRUE, robust = FALSE),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)){
      log_likelihood_mt(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
    }
  )
}

#' Balance of Evidence 4 Threshold LBA model
#'
#' Intermediate thresholds proportional to b, can be less than A
#'
#' @return A model list with all the necessary functions to sample
#' @export
BE4lbaBP <- function(){
  list(
    type="MT",
    # p_vector transform, sets sv as a scaling parameter
    p_types=c("v","sv","r","A","DT1","DT2","DT3","B","t0"),
    transform = function(p) p,
    # Transform to natural scale
    Ntransform=function(x) {
      x[,dimnames(x)[[2]] %in% c("A","B","t0","sv")] <-
        exp(x[,dimnames(x)[[2]] %in% c("A","B","t0","sv")])
      isDT <- substr(dimnames(x)[[2]],1,2)=="DT"
      x[,isDT] <- pnorm(x[,isDT])
      x[,dimnames(x)[[2]]=="r"] <- 2*pnorm(x[,dimnames(x)[[2]]=="r"])-1
      x
    },
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars <- cbind(pars,b=pars[,"B"] + pars[,"A"])
      isDT <- substr(dimnames(pars)[[2]],1,2)=="DT"
      pars[,isDT] <- pars[,isDT]*pars[,"b"]
      attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        ( abs(pars[,"v"])< 10) & ((pars[,"sv"] > 1e-3) | pars[,"sv"] == 0) & (abs(pars[,"b"])<10)
      pars
    },
    # Random function for racing accumulator
    rfun=function(lR=NULL,pars) {
      ok <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        ( abs(pars[,"v"])< 10) & ((pars[,"sv"] > 1e-3) | pars[,"sv"] == 0) & (abs(pars[,"b"])<10)
      if (is.null(lR)) ok else rLBA_BE(pars,ok=ok)
    },
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dLBA(rt,pars,posdrift = TRUE, robust = FALSE),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pLBA(rt,pars,posdrift = TRUE, robust = FALSE),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)){
      log_likelihood_mt(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
    }
  )
}


#' Threshold Count 4 LBA model
#'
#' All thresholds greater than A, DT1 <- A+DT1, DT2 <- DT1 + DT2, DT3 <- DT2 + DT3, b <- DT3 + B
#'
#' @return A model list with all the necessary functions to sample
#' @export
TC4lbaB <- function(){
  list(
    type="TC",
    # p_vector transform, sets sv as a scaling parameter
    p_types=c("v","sv","r","A","DT1","DT2","DT3","B","t0"),
    transform = function(p) p,
    # Transform to natural scale
    Ntransform=function(x) {
      # x[,dimnames(x)[[2]] %in% c("A","B","t0","sv")] <-
      #   exp(x[,dimnames(x)[[2]] %in% c("A","B","t0","sv")])
      # isDT <- substr(dimnames(x)[[2]],1,2)=="DT"
      # x[,isDT] <- pnorm(x[,isDT])
      x[,dimnames(x)[[2]]=="r"] <- 2*pnorm(x[,dimnames(x)[[2]]=="r"])-1
      x[,!(dimnames(x)[[2]] %in% c("v","r"))] <- exp(x[,!(dimnames(x)[[2]] %in% c("v","r"))])
      x
    },
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      # pars <- cbind(pars,b=pars[,"B"] + pars[,"A"])
      # isDT <- substr(dimnames(pars)[[2]],1,2)=="DT"
      # pars[,isDT] <- pars[,isDT]*pars[,"b"]
      pnams <- c("A",dimnames(pars)[[2]][substr(dimnames(pars)[[2]],1,2)=="DT"],"B")
      pars[,pnams] <- t(apply(pars[,pnams],1,cumsum))
      dimnames(pars)[[2]][dimnames(pars)[[2]]=="B"] <- "b"
      attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        ( abs(pars[,"v"])< 10) & ((pars[,"sv"] > 1e-3) | pars[,"sv"] == 0) & (abs(pars[,"b"])<10)
      pars
    },
    # Random function for racing accumulator
    rfun=function(lR=NULL,pars) {
      ok <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0) &
        ( abs(pars[,"v"])< 10) & ((pars[,"sv"] > 1e-3) | pars[,"sv"] == 0) & (abs(pars[,"b"])<10)
      if (is.null(lR)) ok else rLBA_TC(pars,ok=ok)
    },
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dLBA(rt,pars,posdrift = TRUE, robust = FALSE),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pLBA(rt,pars,posdrift = TRUE, robust = FALSE),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)){
      log_likelihood_mt(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
    }
  )
}


albaB <- function(){
  list(
    type="RACE",
    p_types=c("v_0","v_S","v_D","sv","B","A","t0"),
    # Transform to natural scale
    Ntransform=function(x) {
      x[,dimnames(x)[[2]] != "v"] <- exp(x[,dimnames(x)[[2]] != "v"])
      x
    },
    # p_vector transform, sets sv as a scaling parameter
    transform = function(p) {p},
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars <- cbind(pars,v = pars[,"v_0"] + pars[,"v_D"]*dadm$SD + pars[,"v_S"]*dadm$SS)
      pars <- cbind(pars,b=pars[,"B"] + pars[,"A"])
      attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0)
      pars

    },
    # Random function for racing accumulator
    rfun=function(lR=NULL,pars) {
      ok <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0)
      if (is.null(lR)) ok else rLBA(lR,pars,posdrift=TRUE,ok=ok)
    },
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dLBA(rt,pars,posdrift = TRUE, robust = FALSE),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pLBA(rt,pars,posdrift = TRUE, robust = FALSE),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)){
      log_likelihood_race(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
    }
  )}
