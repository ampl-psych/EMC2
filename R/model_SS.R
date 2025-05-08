<<<<<<< HEAD
log_likelihood_race_ss <- function(pars,dadm,model,min_ll=log(1e-10)){

  # pars is a named parameter matrix. nrow = number of trials times the number of go-accumulators (here 2: left and right), ncol = number of parameters plus SSD (stop-signal delay)
  # the columns are: mu, sigma, tau --> ex-Gaussian parameters of the go racer,
  # muS, sigmaS, tauS --> ex-Gaussian parameters of the stop racer
  # tf and gf --> probability that the stop racer is not triggered (tf), probability that the go racer is not triggered (gf)
  # SSD --> stop-signal delay in seconds. is finite on stop-trials and infinite (Inf) on go-trials

  # dadm is a data augmented model. it is a data frame


  # Set up indices:
  # "is" = logical, "isp" pars/dadm index, "t" trials index, "ist" logical on trial index
  # "n_" number of integer

  if (is.null(attr(pars,"ok"))){ # ok is a vector of length
    ok <- !logical(dim(pars)[1]) # creates a vector with TRUE of trial number length
  } else {
    ok <- attr(pars,"ok") # gets the bounds check from the parameter object
  }

    # # spurious go winners on no-response trials
    # dadm$winner[is.na(dadm$R)] <- FALSE
    if (any(is.infinite(dadm$rt))){
      stop("RT contains infinite values")
    }

    # Counts
    # create lI factor that allows for later model flexibility: addint stop-triggered accumulatory ("st")
    # in the standard model, lI is always "go", or TRUE. lI may lateron be user-supplied as a function
    n_accG <- length(levels(dadm$lR))                   # total number of go accumulators

    # check if stop-triggered accumulators have been specified
   if(!is.null(dadm$lI)){
     n_accST <- sum(dadm[1:n_acc,"lI"]==levels(dadm$lI)[1])  # stop-triggered accumulators
     STR <- levels(dadm$lR)[dadm[1:n_acc,"lI"]==levels(dadm$lI)[1]]
   } else {
     n_accST <- 0
     # adds latent inhibitor to dadm, is only "go" in simple stop-signal paradigm
     dadm$lI <- factor(dadm$lR %in% levels(dadm$lR),levels=c(FALSE,TRUE),labels=c("st","go"))
   }

    GOR <- levels(dadm$lR)[dadm[1:n_acc,"lI"]==levels(dadm$lI)[2]]

    # Likelihood for all trials and for ok trials
    allLL <- numeric(nrow(dadm)/n_acc) # placeholder vector of the length of number of trials for participant
    allok <- ok[c(dadm$lR==levels(dadm$lR)[1])] # used to put back into allLL --> ok of length of nr of trials for participant

    # remove bad trials, from here onwards only use the rows where parameters are within bounds
    pars <- pars[ok,,drop=FALSE]
    dadm <- dadm[ok,,drop=FALSE]

    # isp1 indexes the pars/dadm
    isp1 <- dadm$lR==levels(dadm$lR)[1]      # 1st accumulator rows
    ispGOacc <- dadm$lI==levels(dadm$lI)[2] # picks out the go accumulator rows (in simple stop signal, there are only go accumulator rows)
    ispStop <- is.finite(dadm$SSD) # picks out the stop-trial rows
    gf <- pars[isp1,"gf"]
    tf <- pars[isp1,"tf"]

    # first deal with non-responses --------------------------------------------

    # If no response on any of the stop or go trials
    ispNR <- is.na(dadm$R) # no response
    if (any(ispNR)) { # If there are missing responses in the data (to be expected, also in the simple SST)
      ispgoNR <- ispNR & !ispStop # no response on go trial
      tgoNR <- c(1:sum(isp1))[ispgoNR[isp1]] # go trials with no response
      if (any(ispgoNR)) { # any non-responses on go-trials (go failures)
        allLL[allok][tgoNR] <- log(gf[tgoNR]) # put log of go failure probability in there
      }
      ispstopNR <- ispNR & ispStop # successful stop
      tstopNR <- rep(1:(nrow(dadm)/n_acc),each=n_acc)[ispstopNR & isp1] # stop trials
      if (any(ispstopNR)) { # Stop trials
        if (n_accST>0) {
          pStop <- 0 # you don't stop if there is a stop-triggered accumulator, you respond in those cases
          # ...
        } else {
          pStop <- pmax(0,model$sfun(pars[ispStop & ispGOacc & ispNR,,drop=FALSE],n_acc=n_accG)) # sfun computes the probability of stopping
        }
        allLL[allok][tstopNR] <- log(gf[tstopNR] + (1-gf[tstopNR])*(1-tf[tstopNR])*pStop) # something is not right here, gives back -Inf and NAs
      }
    }

    # then deal with responses -------------------------------------------------

    # first remove no response trials
    allr <- !ispNR[isp1] # used to put back into allLL
    pars <- pars[!ispNR,,drop=FALSE]
    dadm <- dadm[!ispNR,,drop=FALSE]
    isp1 <- dadm$lR==levels(dadm$lR)[1]      # 1st accumulator rows
    ispGOacc <- dadm$lI==levels(dadm$lI)[2] # Go accumulator rows
    ispStop <- is.finite(dadm$SSD) # stop-trial rows
    gf <- pars[isp1,"gf"]
    tf <- pars[isp1,"tf"]

    n_trials <- nrow(dadm)/n_acc # number of trial
    trials <- 1:n_trials         # trial number
    ptrials <- rep(trials,each=n_acc) # trials replicated
    like <- numeric(n_trials)
    lds <- numeric(nrow(dadm)) # log density and survivor, used for both go and stop trials; length of dadm, for the winner the density, for the loser the survivor

    # deal with go trials with response
    if (any(!ispStop)) {
      ispGOwin <-  !ispStop & dadm$winner # pick out go winners
      tGO <- trials[!ispStop[isp1]] # pick out go winner trials
      # Winner density, just deal with the standard race
      lds[ispGOwin] <- log(model$dfunG(
        # go win RTs and go win parameters
        rt=dadm[ispGOwin,"rt"],pars=pars[ispGOwin,,drop=FALSE]))
      like[tGO] <- lds[ispGOwin] # density of the winner
      if (n_accG >1) {  # Looser survivor go accumulator(s)
        ispGOloss <- !ispStop & !dadm$winner & ispGOacc # Looser go accumulator rows
        lds[ispGOloss] <- log(1-model$pfunG(
          rt=dadm$rt[ispGOloss],pars=pars[ispGOloss,,drop=FALSE]))
        # add the log of all of the survivors
        like[tGO] <- like[tGO] + apply(matrix(lds[ispGOloss],nrow=n_accG-1),2,sum)
      }
      like[tGO] <- (1-gf[tGO])*exp(like[tGO]) # make a likelihood on the likelihood scale because gf is a probability
    }

    # Deal with Stop trials with a response
    if (any(ispStop)) {
      # Stop looses
      ispSwin <-       ispStop & dadm$winner              # Winner go of stop accumulator rows
      ispSlossGOacc <- ispStop & !dadm$winner & ispGOacc  # Loosing go accumulator rows
      ispSlossSTacc <- ispStop & !dadm$winner & !ispGOacc # Loosing ST accumulator rows

      # pStop at observed rt if ST present (calculate before correcting rt)
      if (n_accST>0) { # probably not working yet
        tST <- trials[ispStop[isp1] & as.numeric(dadm$lI)[dadm$winner] == 1]
        ispST <- ptrials %in% tST
        if (any(ispST)) {
          pStop <- numeric(n_trials)
          upper <- dadm$rt[dadm$winner][tST]
          pStop[tST] <- pmax(0,model$sfun(pars[ispST,,drop=FALSE],
                                                         upper=upper,n_acc=n_acc,st=c(1,1+accST)))
        }
      }

      # For following race calculations correct rt with SSD for ST accumulators
      if (any(ispSlossSTacc)) { # subtract SSD for the stop loosers
        # dadm will not be used later on so can be "corrected" here
        dadm[ispSlossSTacc,"rt"] <- dadm[ispSlossSTacc,"rt"]-dadm[ispSlossSTacc,"SSD"]
      }

      # Fill in lds and sums over survivors for race
      lds[ispSwin] <- log(model$dfunG(rt=dadm[ispSwin,"rt"],pars=pars[ispSwin,,drop=FALSE])) # signal-respond RTs winning go
      if (n_acc >1) {  # Survivor for looser go and/or ST accumulator(s) / signal-respond RTs losing go
        lds[ispSlossGOacc | ispSlossSTacc] <- log(1-model$pfunG(
          rt=dadm[ispSlossGOacc | ispSlossSTacc,"rt"],
          pars=pars[ispSlossGOacc | ispSlossSTacc,,drop=FALSE]))
        # Sum survivor over loosing ST and GO accumulators
        SSTGO <- tapply(lds[ispSlossGOacc | ispSlossSTacc],
                        cbind.data.frame(trials=ptrials[ispSlossGOacc | ispSlossSTacc],
                                         lI=dadm[ispSlossGOacc | ispSlossSTacc,"lI"]),sum)
        SSTGO[is.na(SSTGO)] <- 0 # cases where no ST or GO survivor
      } else {
        SSTGO <- matrix(0,ncol=2)
      }

      # Stop accumulator survivor (no integration)
      sStop <- log(1-model$pfunS(rt=dadm[ispSwin,"rt"],pars=pars[ispSwin,,drop=FALSE]))

      # Put all of it together, get likelihood
      tS <- trials[ispStop[isp1]]
      # Sub-select from tS
      istSgo <- dadm$R[isp1][tS] %in% GOR
      # Stop looses (is the survivor) or not present (trigger failure), can produce ST but only if wins in the absolute race
      like[tS] <- # no failures, works for all responses
        (1-gf[tS])*(1-tf[tS])*exp(lds[ispSwin]+SSTGO[,1]+SSTGO[,2]+sStop)
      # if the trials in which
      if (any(istSgo)){ # stop trial with a go response
        like[tS][istSgo] <- like[tS][istSgo] + # tf (no stop runner), no gf, only produces GO responses
          (1-gf[tS][istSgo])*(tf[tS][istSgo])*exp(lds[ispSwin][istSgo]+SSTGO[,2][istSgo])
        # If both tf and gf then no response, handled previously
      }

      # if there are stop-triggered accumulators --------------------------
      if (n_accST>0) {
        ### Stop wins at some time before rt, and so must be ST response, add to ST race winners
        istSst <- dadm$R[isp1][tS] %in% STR
        like[tS][istSst] <- like[tS][istSst] + # gf (no go runner) no tf, only produces ST responses
          (gf[tS][istSst])*(1-tf[tS][istSst])*exp(lds[ispSwin][istSst]+SSTGO[,1][istSst]+sStop[istSst])
        SST <- numeric(n_trials)
        # Winner and looser ST density already computed.
        SST[tST] <- lds[dadm$winner][tST]
        if (n_accST>1) SST[tST] <- SST[tST] + apply(matrix(
          lds[(ptrials %in% tST) & !dadm$winner & !ispGOacc],nrow=n_accST-1),2,sum)
        like[tS][istSst] <- like[tS][istSst] + pStop[tS][istSst]*exp(SST[tST])
      }
    }

    # put everything together
    allLL[allok][allr] <- log(like)
    return(sum(pmax(min_ll,allLL)))
}



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


# Stop signal random

rSSexGaussian <- function(lR,pars)
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars must contain an SSD column and an lI column. If SSD contains any
  # NAs then return the dt matrix (for use in staircase creation), else return
  # the usual Rrt data frame.
  # NB1: Go failures will only apply to accumulators where lI = TRUE
  #      and can still have a stop-triggered response on a go-failure trial.
{

  nacc <- length(levels(lR)) # Does not include stop runner
  ntrials <- dim(pars)[1]/nacc
  is1 <- lR==levels(lR)[1]
  acc <- 1:nacc

  # stop-triggered racers
  isST <- pars[,"lI"]==1
  accST <- acc[pars[1:nacc,"lI"]==1]

  # Default Inf finishing time so if not changed always looses
  dt <- matrix(Inf,nrow=nacc+1,ncol=ntrials)

  # Go failures
  isgf <- rep(pars[is1,"gf"] > runif(ntrials),each=nacc)
  # Expand to match go accumulators that don't fail
  isGO <- !isgf & !isST
  ngo <- sum(isGO)

  # Fill in go accumulators
  if (any(isGO)) dt[-1,][isGO] <-
    rnorm(ngo,mean=pars[isGO,"mu"],sd=pars[isGO,"sigma"]) +
    rexp(ngo,rate=1/pars[isGO,"tau"])

  # pick out stop trials and races with SSD that is not Inf (i.e., finite or
  # NA, the latter so a staircase can be filled in)
  isStrial <- !is.infinite(pars[is1,"SSD"])
  isSrace <- rep(isStrial,each=nacc)

  # pick out stop trials that are triggered
  isT <- pars[is1,"tf"][isStrial] < runif(sum(isStrial))

  # Logical to pick stop-triggered accumulators that are triggered
  isSTT <- logical(ntrials*nacc)

  # Pick out stop-triggered accumulators that are triggered
  isSTT[isSrace][rep(isT,each=nacc) & isST[isSrace]] <- TRUE
  nst <- sum(isSTT)

  # Fill in stop-triggered accumulators
  if (any(isSTT)) dt[-1,][isSTT] <-
    rnorm(nst,mean=pars[isSTT,"mu"],sd=pars[isSTT,"sigma"]) +
    rexp(nst,rate=1/pars[isSTT,"tau"])

  # pick out triggered stop racers
  isTS <- logical(ntrials)
  isTS[isStrial][isT] <- TRUE
  ns <- sum(isTS)

  # Fill in stop accumulators
  if (any(isTS)) dt[1,isTS] <-
    rnorm(ns,mean=pars[is1,"muS"][isTS],sd=pars[is1,"sigmaS"][isTS]) +
    rexp(ns,rate=1/pars[is1,"tauS"][isTS])

  # return dt to be used by a staircase algorithm
  if (any(is.na(pars[,"SSD"]))) return(dt)

  if (any(isTS)) dt[1,isTS] <- dt[1,isTS] + pars[is1,"SSD"][isTS]
  if (any(isSTT)) dt[-1,][isSTT] <- dt[-1,][isSTT] + pars[isSTT,"SSD"]

  # All SSD already filled in so return R and rt

  # R <- factor(rep(NA,ntrials),levels=levels(lR))
  R <- rt <- rep(NA,ntrials)

  # All accumulators Inf (usually when both go and tf)
  allinf <- apply(dt,2,\(x)all(is.infinite(x)))

  # get winner of stop and go where there is a race
  r <- c(1, 1 + acc)[apply(dt[,!allinf,drop=FALSE],2,which.min)]

  # stop wins
  stopwins <- r==1

  # First fill in cases where stop looses
  if (any(!stopwins)) {
    rgo <- r[!stopwins]-1
    R[!allinf][!stopwins] <- rgo
    pick <- cbind(rgo,c(1:sum(!stopwins))) # Matrix to pick winner
    rt[!allinf][!stopwins] <- dt[-1,!allinf,drop=FALSE][,!stopwins,drop=FALSE][pick]
  }

  # then if stop triggers extra accumulators find their winner
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
  cbind.data.frame(R=factor(R,levels=1:nacc,labels=levels(lR)),rt=rt) #,SSD=SSD)
}

#' Staircase function
#'
#' Staircase function for stop-signal models. Can be used to generate
#' stop-signal data.
#'
#' @param dadm A data-augmented design
#' @param p Proportion of stop-trials
#' @param pars Parameter matrix, needs to have the same number of rows as `dadm`
#' @param SSD0 Minimum SSD in seconds
#' @param stairstep Step size in seconds
#' @param stairmin Minimum SSD in seconds
#' @param stairmax Maximum SSD in seconds
#'
#' @return A data frame with variables `R`, `rt`, and `SSD`
#' @export
#'
staircase_fun <- function(dadm,p=.25,pars=NULL,
                          SSD0=.25,stairstep=.05,stairmin=0,stairmax=Inf)
  # random p of trials get NA, ready to be filled in by a staircase
{
  # if pars not supplied return an SSD column indicating all go trials
  if (is.null(pars)) return(rep(Inf,nrow(dadm)))

  # levels(dadm$lR) <- levels(dadm$lR)

  # if pars supplied potentially run staircase for each participant
  for(i in 1:length(unique(dadm$subjects))){
    dadm_i <- dadm[dadm$subjects==unique(dadm$subjects)[i],]
    pars_i <- pars[dadm$subjects==unique(dadm$subjects)[i],]
    nacc <- length(levels(dadm_i$lR))
    is1 <- dadm_i$lR==levels(dadm_i$lR)[1] # first accumulator
    ntrials <- sum(is1)
    if (!any(colnames(dadm_i)=="SSD")) { # pick staircase trials
      dadm_i$SSD <- rep(Inf, nrow(dadm_i))
      tmp <- matrix(dadm_i$SSD,nrow=nacc)
      tmp[,sample(1:ntrials,round(ntrials*p))] <- NA
      dadm_i$SSD <- as.vector(tmp)
    }
    out_i <- setNames(data.frame(matrix(NA,nrow=ntrials,ncol=3)),c("R","rt","SSD"))
    pick <- is.infinite(dadm_i$SSD)
    if (any(pick)) { # fill in no-stop trials
      out_i$SSD[pick[is1]] <- Inf
      pars_i[pick,"SSD"] <- Inf
      out_i[pick[is1],c("R","rt")] <- attributes(dadm_i)$model()$rfun(dadm_i$lR[pick],pars_i[pick,,drop=FALSE])
    }
    pick <- is.finite(dadm_i$SSD)
    if (any(pick)) { # fill in fixed SSD trials
      out_i$SSD[pick[is1]] <- dadm_i$SSD[pick][is1[pick]]
      pars_i[pick,"SSD"] <- dadm_i$SSD[pick]
      out_i[pick[is1],c("R","rt")] <- attributes(dadm_i)$model()$rfun(dadm_i$lR[pick],pars_i[pick,,drop=FALSE])
    }
    isna <- is.na(dadm_i$SSD)
    if (any(isna)) { # run staircase if any NAs to fill in
      nstair <- sum(isna)/nacc # number of staircase trials
      trials <- rep(0,ntrials*nacc) # used to pick out_i each staircase trial in pars_i
      trials[is.na(dadm_i$SSD)] <- rep(1:nstair,each=nacc)
      for (s in 1:nstair) { # run staircase
        current <-  trials == s # current staircase trial
        if (s==1)  dadm_i$SSD[current] <- out_i$SSD[current[is1]] <- SSD0 # initialize
        p_stairs <- pars_i[current,,drop=FALSE] # parameters for current staircase trial
        # simulate 1 trial, because is.na(SSD) in pars_i rfun returns dt, an nacc x 1 matrix
        dt <- attributes(dadm_i)$model()$rfun(dadm_i$lR[current],p_stairs)
        inhibit <- p_stairs[,"lI"]==1 # inhibition triggered
        # add SSD to stop and inhibition triggered
        dt[c(TRUE,inhibit),] <- dt[c(TRUE,inhibit),] + dadm_i$SSD[current][1]
        winner <- which.min(dt)
        iinhibit <- c(1,1+c(1:nacc)[inhibit]) # stop or inhibition triggered index
        if (s != nstair) { # set SSD for next staircase trial
          nexts <- trials == (s+1)
          if (any(iinhibit==winner))   # round as otherwise get spurious tiny differences
            dadm_i$SSD[nexts] <- round(dadm_i$SSD[current] + stairstep,3) else  # successful stop
              dadm_i$SSD[nexts] <- round(dadm_i$SSD[current] - stairstep,3)       # failed stop
            if ((dadm_i$SSD[nexts][1]<stairmin) | (dadm_i$SSD[nexts][1]>stairmax))
              dadm_i$SSD[nexts] <- dadm_i$SSD[current] # dont step
            out_i$SSD[c(trials == s+1)[is1]] <- dadm_i$SSD[nexts][1]
        }
        if (winner==1) { # stop wins
          if (any(inhibit)) { # inhibition triggered response
            if (all(is.infinite(dt[c(FALSE,inhibit),])))  # tf
              out_i[current[is1],c("R","rt")] <- c(NA,NA) else {
                pick <- which.min(dt[c(FALSE,inhibit),])
                out_i[current[is1],c("R","rt")] <-
                  c(c(1:nacc)[inhibit][pick],dt[c(FALSE,inhibit),][pick])
              }
          } # otherwise no response
        } else { # pick from all except stop
          out_i[current[is1],c("R","rt")] <- c(winner-1,dt[winner,])
        }
      }
    }
    out_i$R <- factor(out_i$R,levels=1:nacc,labels=levels(dadm_i$lR))
    if(i == 1){
      out <- out_i
    } else {
      out <- rbind(out, out_i)
    }
    if(length(unique(dadm$subjects)) == 1){
      return(out)
    }
=======
#### Staircase ----
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
>>>>>>> dev
  }
  return(out)
}


<<<<<<< HEAD

=======
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
>>>>>>> dev

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
<<<<<<< HEAD
=======

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


#### ExG Race function ----

# Following functions moved to C++ model_SS_EXG.cpp

>>>>>>> dev
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
<<<<<<< HEAD
=======

>>>>>>> dev
#
#
# stopfn_exg <- function(t,mu,sigma,tau,SSD)
#   # Used by my.integrate, t = vector of times, SSD is a scalar stop-signal delay.
# {
#   dt <- matrix(rep(t+SSD,each=length(mu)),nrow=length(mu))
#   dt[1,] <- dt[1,]-SSD
#   dEXGrace(dt,mu,sigma,tau)
# }

<<<<<<< HEAD
pstopEXG <- function(parstop,n_acc,upper=Inf,
                     gpars=c("mu","sigma","tau"),spars=c("muS","sigmaS","tauS"))
{
  sindex <- seq(1,nrow(parstop),by=n_acc)
  ps <- parstop[sindex,spars]
  SSDs <- parstop[sindex,"SSD"]
=======
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

rSSexGaussian <- function(lR,pars,ok=rep(TRUE,dim(pars)[1]))
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars must contain an SSD column and an lI column indicating if an
  # accumulator is triggered by the stop signal (ST = stop-triggered).

  # NB1: Go failures will only apply to accumulators where lI = TRUE
  #      and can still have a stop-triggered response on a go-failure trial.
{

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
    if (is.null(attr(pars,"staircase")))
      stop("When SSD has NAs a staircase list must be supplied!")
    staircase <- attr(pars,"staircase")
    staircase$data <- staircase$data[is1,]
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
    if (!is.null(staircase$data))
      staircase$data <- staircase$data[stair,]

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
    out <- cbind.data.frame(R=factor(allR,levels=1:nacc,labels=levels(lR)),
                            rt=allrt)
    attr(out,"SSD") <- stair_res$SSD
    return(out)
  }
  cbind.data.frame(R=factor(R,levels=1:nacc,labels=levels(lR)),rt=rt)
}


#### ExG stop probability (no stop triggered) ----

pstopEXG <- function(parstop,n_acc,upper=Inf,
                     gpars=c("mu","sigma","tau"),spars=c("muS","sigmaS","tauS"))
{
  sindex <- seq(1,nrow(parstop),by=n_acc)  # Stop accumulator index
  ps <- parstop[sindex,spars,drop=FALSE]   # Stop accumulator parameters
  SSDs <- parstop[sindex,"SSD",drop=FALSE] # SSDs
>>>>>>> dev
  ntrials <- length(SSDs)
  if (length(upper)==1) upper <- rep(upper,length.out=ntrials)
  pgo <- array(parstop[,gpars],dim=c(n_acc,ntrials,length(gpars)),
               dimnames=list(NULL,NULL,gpars))
<<<<<<< HEAD
  cells <- apply(cbind(SSDs,ps,upper,
                       matrix(as.vector(aperm(pgo,c(2,1,3))),nrow=ntrials)),1,paste,collapse="")
=======
  cells <- apply(
    cbind(SSDs,ps,upper,matrix(as.vector(aperm(pgo,c(2,1,3))),nrow=ntrials))
  ,1,paste,collapse="")
>>>>>>> dev
  # cells <- character(ntrials)
  # for (i in 1:ntrials)
  #   cells[i] <- paste(SSDs[i],ps[i,],pgo[,i,],upper[i],collapse="")
  uniq <- !duplicated(cells)
  ups <- sapply(1:sum(uniq),function(i){
    my.integrate(f=stopfn_exg,lower=-Inf,SSD=SSDs[i],upper=upper[i],
<<<<<<< HEAD
                 mu=c(ps[i,"muS"],pgo[,i,"mu"]),
                 sigma=c(ps[i,"sigmaS"],pgo[,i,"sigma"]),
                 tau=c(ps[i,"tauS"],pgo[,i,"tau"]))
=======
                           mu=c(ps[i,"muS"],pgo[,i,"mu"]),
                           sigma=c(ps[i,"sigmaS"],pgo[,i,"sigma"]),
                           tau=c(ps[i,"tauS"],pgo[,i,"tau"]))
>>>>>>> dev
  })
  ups[as.numeric(factor(cells,levels=cells[uniq]))]
}

<<<<<<< HEAD
stopfn_exgST <- function(t,mu,sigma,tau,SSD,st=1)
  # Used by my.integrate, t = vector of times, SSD is a scalar stop-signal delay.
  # st is a vector of indices for the stop and stop-triggered accumulators
{
  dt <- matrix(rep(t+SSD,each=length(mu)),nrow=length(mu))
  dt[st,] <- dt[st,]-SSD
  dEXGrace(dt,mu,sigma,tau)
}

pstopEXGST <- function(parstop,n_acc,upper=Inf,st=1,
                       gpars=c("mu","sigma","tau"),spars=c("muS","sigmaS","tauS"))
{
  sindex <- seq(1,nrow(parstop),by=n_acc)
  ps <- parstop[sindex,spars]
  SSDs <- parstop[sindex,"SSD"]
  ntrials <- length(SSDs)
  if (length(upper)==1) upper <- rep(upper,length.out=ntrials)
  pgo <- array(parstop[,gpars],dim=c(n_acc,ntrials,length(gpars)),
               dimnames=list(NULL,NULL,gpars))
  cells <- apply(cbind(SSDs,ps,upper,
                       matrix(as.vector(aperm(pgo,c(2,1,3))),nrow=ntrials)),1,paste,collapse="")
  # cells <- character(ntrials)
  # for (i in 1:ntrials)
  #   cells[i] <- paste(SSDs[i],ps[i,],pgo[,i,],upper[i],collapse="")
  uniq <- !duplicated(cells)
  ups <- sapply(1:sum(uniq),function(i){
    my.integrate(f=stopfn_exgST,lower=-Inf,SSD=SSDs[i],upper=upper[i],
                 mu=c(ps[i,"muS"],pgo[,i,"mu"]),
                 sigma=c(ps[i,"sigmaS"],pgo[,i,"sigma"]),
                 tau=c(ps[i,"tauS"],pgo[,i,"tau"]),st=st)
  })
  ups[as.numeric(factor(cells,levels=cells[uniq]))]
}


#### Model list ----
#' Stop-signal exGaussian race
=======
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
>>>>>>> dev
#'
#' @return A model list with all the necessary functions to sample
#' @export
SSexG <- function() {
  list(
<<<<<<< HEAD
    c_name = NULL,
    type="RACE",
    p_types=c(mu=log(.4),sigma=log(.05),tau=log(.1),
              muS=log(.3),sigmaS=log(.025),tauS=log(.05),tf=qnorm(0),gf=qnorm(0)),
    transform =list(func=c(mu = "exp",sigma = "exp",tau = "exp",
                           muS = "exp",sigmaS = "exp",tauS = "exp",tf = "pnorm", gf = "pnorm")),
    bound=list(minmax=cbind(mu = c(0.001,Inf), sigma = c(0.001,Inf), tau = c(0.001,Inf),
                            muS = c(0.001, Inf), sigmaS = c(0.001, Inf), tauS = c(0.001, Inf),
                            tf = c(-Inf, 1), gf = c(-Inf, 1))),
    Ttransform = function(pars,dadm) {
      if (any(names(dadm)=="SSD")){
        pars <- cbind(pars,SSD=dadm$SSD)
      } else {
        pars <- cbind(pars,SSD=rep(NA,dim(pars)[1]))
      }
      pars <- cbind(pars,lI=as.numeric(dadm$lI))  # Only necessary for data generation.
      return(pars)
    },
    # Random function for SS race
    rfun=function(lR=NULL,pars) {
      rSSexGaussian(lR,pars)
=======
    type="RACE",
    p_types=c(mu=log(.4),sigma=log(.05),tau=log(.1),
              muS=log(.3),sigmaS=log(.025),tauS=log(.05),tf=qnorm(0),gf=qnorm(0)),
    transform=list(func=c( mu = "exp",  sigma = "exp",  tau = "exp",
                          muS = "exp", sigmaS = "exp", tauS = "exp",
                          tf="pnorm",gf="pnorm")),
    bound=list(minmax=cbind( mu=c(0,Inf),  sigma=c(0,Inf),  tau=c(1e-4,Inf),
                            muS=c(0,Inf), sigmaS=c(0,Inf), tauS=c(1e-4,Inf),
                            tf=c(.001,.999),gf=c(.001,.999)),
                            exception=c(tf=0,gf=0)),
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      # if (any(names(dadm)=="SSD")) pars <- cbind(pars,SSD=dadm$SSD) else
      #   pars <- cbind(pars,SSD=rep(Inf,dim(pars)[1]))
      pars <- cbind(pars,SSD=dadm$SSD)
      pars <- cbind(pars,lI=as.numeric(dadm$lI))  # Only necessary for data generation.
      pars
>>>>>>> dev
    },
    # Density function (PDF) for single go racer
    dfunG=function(rt,pars) dexGaussianG(rt,pars),
    # Probability function (CDF) for single go racer
    pfunG=function(rt,pars) pexGaussianG(rt,pars),
    # Density function (PDF) for single stop racer
<<<<<<< HEAD
    dfunS=function(rt,pars) dexGaussianS(rt,pars[,c("muS","sigmaS","tauS","SSD"),drop=FALSE]),
    # Probability function (CDF) for single stop racer
    pfunS=function(rt,pars) pexGaussianS(rt,pars[,c("muS","sigmaS","tauS","SSD"),drop=FALSE]),
    # Stop probability integral
    sfun=function(pars,n_acc,st=1,upper=Inf) pstopEXGST(pars,n_acc,upper=upper,st=st),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10)){
      log_likelihood_race_ss(pars=pars, dadm = dadm, model= model, min_ll = min_ll)
    }
  )
}

my.integrate <- function(...,upper=Inf,big=10)
  # Avoids bug in integrate upper=Inf that uses only 1  subdivision
  # Use of  big=10 is arbitrary ...
{
  out <- try(integrate(...,upper=upper),silent=TRUE)
  if (is(out,"try-error")) 0 else
  {
    if (upper==Inf & out$subdivisions==1)
    {
      out <- try(integrate(...,upper=big),silent=TRUE)
      if (is(out,"try-error")) 0 else
      {
        if (out$subdivisions==1) 0 else out$value
      }
    } else out$value
  }
=======
    dfunS=function(rt,pars)
      dexGaussianS(rt,pars[,c("muS","sigmaS","tauS","SSD"),drop=FALSE]),
    # Probability function (CDF) for single stop racer
    pfunS=function(rt,pars)
      pexGaussianS(rt,pars[,c("muS","sigmaS","tauS","SSD"),drop=FALSE]),
    # Stop probability integral
    sfun=function(pars,n_acc,upper=Inf) pstopEXG(pars,n_acc,upper=upper),
    # Random function for SS race
    rfun=function(lR=NULL,pars) {
      rSSexGaussian(lR,pars,ok=attr(pars, "ok"))
    },
    # Race likelihood combining pfun and dfun
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10))
      log_likelihood_race_ss(pars, dadm, model, min_ll = min_ll)
  )
}


#####################  RDEX ----

#### RDEX SS random ----

rSShybrid <- function(lR,pars,ok=rep(TRUE,dim(pars)[1]))
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars must contain an SSD column and an lI column indicating if an
  # accumulator is triggered by the stop signal (ST = stop-triggered).

  # NB1: Go failures will only apply to accumulators where lI = TRUE
  #      and can still have a stop-triggered response on a go-failure trial.
{

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
    staircase$data <- staircase$data[is1,]
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
    if (!is.null(staircase$data))
      staircase$data <- staircase$data[stair,]

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
    out <- cbind.data.frame(R=factor(allR,levels=1:nacc,labels=levels(lR)),
                            rt=allrt)
    attr(out,"SSD") <- stair_res$SSD
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
                            s=c(0,Inf),muS=c(0,Inf), sigmaS=c(0,Inf), tauS=c(1e-4,Inf),
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
    rfun=function(lR=NULL,pars) {
      rSShybrid(lR,pars,ok=attr(pars, "ok"))
    },
    # Race likelihood combining pfun and dfun
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10))
      log_likelihood_race_ss(pars, dadm, model, min_ll = min_ll)
  )
>>>>>>> dev
}

