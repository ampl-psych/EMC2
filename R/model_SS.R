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

# vectorised ex-Gaussian (log) PDF
# written to be equivalent to version in src/exgaussian_functions.h
dexGaussian <- function(rt, pars, log.d = FALSE) {
  log_out <- vector("numeric", length(rt))
  mu <- pars[, "mu"]
  sigma <- pmax(pars[, "sigma"], 1e-12)
  tau <- pmax(pars[, "tau"], 1e-12)
  z <- (rt - mu) / sigma - sigma / tau
  log_phi <- pnorm(z, log.p = TRUE)
  isnan <- is.nan(log_phi)
  if (any(isnan)) {
    log_out[isnan] <- NA_real_
  }
  isinf <- is.infinite(log_phi)
  if (any(isinf)) {
    log_out[isinf] <- log_phi[isinf]
  }
  ok <- !isnan & !isinf
  if (any(ok)) {
    log_exp <- (mu[ok] - rt[ok]) / tau[ok] + sigma[ok]^2 / (2 * tau[ok]^2)
    log_out[ok] <- -log(tau[ok]) + log_exp + log_phi[ok]
  }
  if (!log.d) {
    return(exp(log_out))
  } else {
    return(log_out)
  }
}

# vectorised ex-Gaussian (log) CDF
# written to be equivalent to version in src/exgaussian_functions.h
pexGaussian <- function(rt, pars, lower.tail = TRUE, log.p = FALSE) {
  log_out_lower <- vector("numeric", length(rt))
  mu <- pars[, "mu"]
  sigma <- pmax(pars[, "sigma"], 1e-12)
  tau <- pmax(pars[, "tau"], 1e-12)
  isinf <- is.infinite(rt)
  if (any(isinf)) {
    log_out_lower[isinf] <- ifelse(rt[isinf] < 0, -Inf, 0)
  }
  ok <- !isinf
  z <-  (rt[ok] - mu[ok]) / sigma[ok]
  log_phi_1 <- pnorm(z, log.p = TRUE)
  log_phi_2 <- pnorm(z - sigma[ok] / tau[ok], log.p = TRUE)
  log_exp_term <- (mu[ok] - rt[ok]) / tau[ok] + sigma[ok]^2 / (2 * tau[ok]^2)
  log_out_lower[ok] <- log_diff_exp(log_phi_1, (log_exp_term + log_phi_2))
  if (lower.tail) {
    if (log.p) {
      return(log_out_lower)
    } else {
      return(exp(log_out_lower))
    }
  } else {
    if (log.p) {
      return(log1m_exp(log_out_lower))
    } else {
      return(-expm1(log_out_lower))
    }
  }
}

# wrapper around exG cdf/pdf for go trials (strips out NAs and calls d/pexGaussian)
# Is stripping out rt NAs really necessary?
dexGaussianG <- function(rt, pars, log.d = FALSE) {
  out <- numeric(length(rt))
  ok <- !is.na(rt)
  out[ok] <- dexGaussian(rt[ok], pars[ok, , drop = FALSE], log.d)
  return(out)
}

pexGaussianG <- function(rt, pars, lower.tail = TRUE, log.p = FALSE) {
  out <- numeric(length(rt))
  ok <- !is.na(rt)
  out[ok] <- pexGaussian(rt[ok], pars[ok, , drop = FALSE], lower.tail, log.p)
  return(out)
}


#### Stop Single ExGaussian ----

# wrapper around exG cdf/pdf for stop process
# subtracts SSD and uses muS/sigmaS/tauS by renaming
dexGaussianS <- function(rt, pars, log.d = FALSE) {
  rt <- rt - pars[,"SSD"]
  dimnames(pars)[[2]][dimnames(pars)[[2]] == "muS"] <- "mu"
  dimnames(pars)[[2]][dimnames(pars)[[2]] == "sigmaS"] <- "sigma"
  dimnames(pars)[[2]][dimnames(pars)[[2]] == "tauS"] <- "tau"
  return(dexGaussian(rt, pars, log.d))
}

pexGaussianS <- function(rt, pars, lower.tail = TRUE, log.p = FALSE) {
  rt <- rt - pars[,"SSD"]
  dimnames(pars)[[2]][dimnames(pars)[[2]] == "muS"] <- "mu"
  dimnames(pars)[[2]][dimnames(pars)[[2]] == "sigmaS"] <- "sigma"
  dimnames(pars)[[2]][dimnames(pars)[[2]] == "tauS"] <- "tau"
  return(pexGaussian(rt, pars, lower.tail, log.p))
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

rexG <- function(n,mu,sigma,tau)
  rnorm(n,mean=mu,sd=sigma) + rexp(n,rate=1/tau)

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

pstopEXG <- function(
    parstop, n_acc,
    upper = Inf, gpars = c("mu","sigma","tau"), spars = c("muS","sigmaS","tauS")
) {
  sindex <- seq(1, nrow(parstop), by = n_acc)   # Stop accumulator index
  ps <- parstop[sindex, spars, drop = FALSE]    # Stop accumulator parameters
  SSDs <- parstop[sindex,"SSD", drop = FALSE]   # SSDs
  ntrials <- length(SSDs)
  if (length(upper) == 1) {
    upper <- rep(upper, ntrials)
  }
  # 3D array of go accumulator parameter values:
  # slices = parameters (mu/sigma/tau); rows = accumulators, columns = trials
  pgo <- array(
    data = parstop[ , gpars],
    dim = c(n_acc, ntrials, length(gpars)),
    dimnames = list(NULL, NULL, gpars)
  )
  # transpose so that rows = trials and columns = accumulators, then collapse into vector:
  # [param1_acc1_trial1, param1_acc1_trail2, ..., param1_acc1_trialN, param1_acc2_trial1, ..., param1_accN_trialN, param2_acc1_trial1, ..., paramN_accN_trialN]
  # finally puts that in a matrix where rows = trials and columns = [param_1_acc1, param1_acc2, ..., param1_accN, param2_acc1, ..., paramN_accN]
  pgo_mat <- matrix(
    as.vector(
      aperm(pgo, c(2, 1, 3))),
    nrow = ntrials
  )
  # bring together all relevant variables: SSD, stop parameters, upper limit of
  # integration, and go parameters
  cells <- apply(
    X = cbind(SSDs, ps, upper, pgo_mat),
    MARGIN = 1,
    FUN = paste,
    collapse = ""
  )
  # get indices of unique parameter combinations across trials
  unique_cells <- !duplicated(cells)
  # calculate stop success probability by integration, only for unique
  # parameter combinations
  ss_integral <- sapply(
    X = 1:sum(unique_cells),
    FUN = function(i) {
      # wrapper function for numerical integration, defined in likelihood.R
      my.integrate(
        f = stopfn_exg, # integrand, defined in src/model_SS_EXG.h
        lower = -Inf, # lower limit of integration
        SSD = SSDs[i],
        upper = upper[i],
        mu = c(ps[i, "muS"], pgo[ , i , "mu"]),
        sigma = c(ps[i, "sigmaS"], pgo[ , i, "sigma"]),
        tau = c(ps[i, "tauS"], pgo[ , i, "tau"])
      )
    }
  )
  # expand into trial-wise likelihood by getting indices of unique cells
  out <- ss_integral[as.numeric(factor(cells, levels = cells[unique_cells]))]
  return(out)
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
#' @return A model list with all the necessary functions to sample
#' @export
SSexG <- function() {
  list(
    type = "RACE",
    c_name = "SS_EXG",
    p_types = c(
      mu = log(.4), sigma = log(.05), tau = log(.1),
      muS = log(.3), sigmaS = log(.025), tauS = log(.05),
      tf = qnorm(0), gf = qnorm(0)#,
      # exg_lb = log(.05), exgS_lb = log(.05)
    ),
    transform = list(
      func = c(
        mu = "exp", sigma = "exp", tau = "exp",
        muS = "exp", sigmaS = "exp", tauS = "exp",
        tf = "pnorm", gf = "pnorm"#,
        # exg_lb = "exp", exgS_lb = "exp"
      )
    ),
    bound = list(
      minmax = cbind(
        mu = c(0, Inf), sigma = c(1e-4, Inf), tau = c(1e-4, Inf),
        muS = c(0, Inf), sigmaS = c(1e-4,Inf), tauS = c(1e-4,Inf),
        tf = c(.001, .999), gf = c(.001, .999)#,
        # exg_lb = c(0, Inf), exgS_lb = c(0, Inf)
      ),
      exception = c(tf = 0, gf = 0)
    ),
    # Trial dependent parameter transform
    Ttransform = function(pars, dadm) {
      pars <- cbind(pars, SSD = dadm$SSD)
      pars <- cbind(pars, lI = as.numeric(dadm$lI)) # only necessary for data generation
      return(pars)
    },
    # Density function (PDF) for single go racer
    dfunG = function(rt, pars, log.d = FALSE) {
      return(dexGaussianG(rt, pars, log.d))
    },
    # Probability function (CDF) for single go racer
    pfunG = function(rt, pars, lower.tail = TRUE, log.p = FALSE) {
      return(pexGaussianG(rt, pars, lower.tail, log.p))
    },
    # Density function (PDF) for single stop racer
    dfunS = function(rt, pars, log.d = FALSE) {
      parsS <- pars[ , c("muS", "sigmaS", "tauS", "SSD"), drop = FALSE]
      return(dexGaussianS(rt, parsS, log.d))
    },
    # Probability function (CDF) for single stop racer
    pfunS = function(rt, pars, lower.tail = TRUE, log.p = FALSE) {
      parsS <- pars[ , c("muS", "sigmaS", "tauS", "SSD"), drop = FALSE]
      return(pexGaussianS(rt, parsS, lower.tail, log.p))
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


#####################  RDEX ----

#### RDEX SS random ----

rSShybrid <- function(data,pars,ok=rep(TRUE,dim(pars)[1]))
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars must contain an SSD column and an lI column indicating if an
  # accumulator is triggered by the stop signal (ST = stop-triggered).

  # NB1: Go failures will only apply to accumulators where lI = TRUE
  #      and can still have a stop-triggered response on a go-failure trial.
{
  lR <- data$lR
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
    my.integrate(f=stopfn_rdex,lower=0,upper=upper[i],
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

