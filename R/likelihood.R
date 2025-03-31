calc_ll_R <- function(p_vector, model, dadm){
  if(!is.null(model$transform)){
    pars <- get_pars_matrix(p_vector, dadm, model)
  } else{
    pars <- p_vector
  }
  ll <- model$log_likelihood(pars, dadm, model)
  return(ll)
}


log_likelihood_race <- function(pars,dadm,model,min_ll=log(1e-10))
  # Race model summed log likelihood
{
  if (any(names(dadm)=="RACE")) # Some accumulators not present
    pars[as.numeric(dadm$lR)>as.numeric(as.character(dadm$RACE)),] <- NA

  if (is.null(attr(pars,"ok")))
    ok <- !logical(dim(pars)[1]) else ok <- attr(pars,"ok")

    lds <- numeric(dim(dadm)[1]) # log pdf (winner) or survivor (losers)
    lds[dadm$winner] <- log(model$dfun(rt=dadm$rt[dadm$winner],
                                                      pars=pars[dadm$winner,]))
    n_acc <- length(levels(dadm$R))
    if (n_acc>1) lds[!dadm$winner] <- log(1-model$pfun(rt=dadm$rt[!dadm$winner],pars=pars[!dadm$winner,]))
    lds[is.na(lds) | !ok] <- min_ll
    lds <- lds[attr(dadm,"expand")] # decompress
    if (n_acc>1) {
      winner <- dadm$winner[attr(dadm,"expand")]
      ll <- lds[winner]
      if (n_acc==2) {
        ll <- ll + lds[!winner]
      } else {
        ll <- ll + apply(matrix(lds[!winner],nrow=n_acc-1),2,sum)
      }
      ll[is.na(ll)] <- min_ll
      return(sum(pmax(min_ll,ll)))
    } else return(sum(pmax(min_ll,lds)))
}



log_likelihood_ddm <- function(pars,dadm,model,min_ll=log(1e-10))
  # DDM summed log likelihood, with protection against numerical issues
{
  like <- numeric(dim(dadm)[1])
  if (any(attr(pars,"ok")))
    like[attr(pars,"ok")] <- model$dfun(dadm$rt[attr(pars,"ok")],dadm$R[attr(pars,"ok")],
                                                       pars[attr(pars,"ok"),,drop=FALSE])
  like[attr(pars,"ok")][is.na(like[attr(pars,"ok")])] <- 0
  sum(pmax(min_ll,log(like[attr(dadm,"expand")])))
}


log_likelihood_ddmgng <- function(pars,dadm,model,min_ll=log(1e-10))
  # DDM summed log likelihood for go/nogo model
{
  like <- numeric(dim(dadm)[1])
  if (any(attr(pars,"ok"))) {
    isna <- is.na(dadm$rt)
    ok <- attr(pars,"ok") & !isna
    like[ok] <- model$dfun(dadm$rt[ok],dadm$R[ok],pars[ok,,drop=FALSE])
    ok <- attr(pars,"ok") & isna
    like[ok] <- # dont terminate on go boundary before timeout
      pmax(0,pmin(1,(1-model$pfun(dadm$TIMEOUT[ok],dadm$Rgo[ok],pars[ok,,drop=FALSE]))))

  }
  like[attr(pars,"ok")][is.na(like[attr(pars,"ok")])] <- 0
  sum(pmax(min_ll,log(like[attr(dadm,"expand")])))
}



#### sdt choice likelihoods ----

log_likelihood_sdt <- function(pars,dadm,lb=-Inf, model, min_ll=log(1e-10))
  # probability of ordered discrete choices based on integrals of a continuous
  # distribution between thresholds, with fixed lower bound for first response
  # lb. Upper bound for last response is a fixed value in threshold vector
{
  first <- dadm$lR==levels(dadm$lR)[1]
  last <- dadm$lR==levels(dadm$lR)[length(levels(dadm$lR))]
  pars[last,"threshold"] <- Inf
  # upper threshold
  ut <- pars[dadm$winner,"threshold"]
  # lower threshold fixed at lb for first response
  pars[first &  dadm$winner,"threshold"] <- lb
  # otherwise threshold of response before one made
  notfirst <- !first &  dadm$winner
  pars[notfirst,"threshold"] <- pars[which(notfirst)-1,"threshold"]
  lt <- pars[dadm$winner,"threshold"]
  # fix race expand to suit SDT
  nr <- length(levels(dadm$lR))      # number of responses
  ne <- length(attr(dadm,"expand"))  # length of expand
  # Shorten expand to only one per lR set
  expand <- (attr(dadm,"expand")[(c(1:ne) %% nr)==0] + 1 ) %/% nr
  # log probability
  ll <- numeric(sum(dadm$winner))
  if (!is.null(attr(pars,"ok"))) { # Bad parameter region
    ok <- attr(pars,"ok")
    okw <- ok[dadm$winner]
    ll[ok] <- log(model$pfun(lt=lt[okw],ut=ut[okw],pars=pars[dadm$winner & ok,,drop=FALSE]))
  } else ll <- log(model$pfun(lt=lt,ut=ut,pars=pars[dadm$winner,,drop=FALSE]))
  ll <- ll[expand]
  ll[is.na(ll)] <- 0
  sum(pmax(min_ll,ll))
}

log_likelihood_joint <- function(proposals, dadms, model_list, component = NULL){
  parPreFixs <- unique(gsub("[|].*", "", colnames(proposals)))
  i <- 0
  total_ll <- 0
  if(!is.null(component)) dadms <- dadms[component]
  for(dadm in dadms){
    if(is.data.frame(dadm)){
      i <- i + 1
      parPrefix <- parPreFixs[i]
      columns_to_use <- sapply(strsplit(colnames(proposals), "|", fixed = TRUE), function(x) x == parPrefix)[1,]
      currentPars <- proposals[,columns_to_use, drop = F]
      colnames(currentPars) <- gsub(".*[|]", "", colnames(currentPars))
      total_ll <- total_ll +  calc_ll_manager(currentPars, dadm, model_list[[i]])
    }
  }
  return(total_ll)
}

####  Stop-signal ----

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
}

log_likelihood_race_ss <- function(pars,dadm,min_ll=log(1e-10))
{

  # Set up indices:
  # "is" = logical, "isp" pars/dadm index, "t" trials index, "ist" logical on trial index
  # "n_" number of integer
  n_acc <- length(levels(dadm$lR))                   # total number of accumulators

  if (is.null(attr(pars,"ok")))
    ok <- !logical(dim(pars)[1]) else ok <- attr(pars,"ok")
  if (!any(ok)) return(min_ll*length(attr(dadm, "expand_winner")))

  # # spurious go winners on no-response trials
  # dadm$winner[is.na(dadm$R)] <- FALSE
  # if (any(is.infinite(dadm$rt))) stop("BUGGER!")

  # Counts
  n_accG <- sum(as.numeric(dadm[1:n_acc,"lI"])==2)   # go accumulators
  n_accST <- sum(as.numeric(dadm[1:n_acc,"lI"])==1)  # stop-triggered accumulators
  GOR <- levels(dadm$lR)[as.numeric(dadm[1:n_acc,"lI"])==2] # Go responses
  if (n_accST>0)
    STR <- levels(dadm$lR)[as.numeric(dadm[1:n_acc,"lI"])==1]

  # Likelihood for all trials and for ok trials
  allLL <- rep(-Inf,nrow(dadm)/n_acc)
  allok <- ok[c(dadm$lR==levels(dadm$lR)[1])] # used to put back into allLL

  # remove bad trials
  pars <- pars[ok,,drop=FALSE]
  dadm <- dadm[ok,,drop=FALSE]
  isp1 <- dadm$lR==levels(dadm$lR)[1]      # 1st accumulator rows
  ispGOacc <- dadm$lI==levels(dadm$lI)[2]  # Go accumulator rows
  ispStop <- is.finite(dadm$SSD) # stop-trial rows
  gf <- pars[isp1,"gf"]
  tf <- pars[isp1,"tf"]

  # No response
  ispNR <- is.na(dadm$R)
  if (any(ispNR)) { # Test as some models always respond
    ispgoNR <- ispNR & !ispStop
    tgoNR <- c(1:sum(isp1))[ispgoNR[isp1]]
    if (any(ispgoNR)) allLL[allok][tgoNR] <- log(gf[tgoNR]) # Go trials
    ispstopNR <- ispNR & ispStop
    tstopNR <- (rep(1:(nrow(dadm)/n_acc),each=n_acc))[ispstopNR & isp1]
    if (any(ispstopNR)) { # Stop trials
      if (n_accST>0) pStop <- 0 else
        pStop <- pmax(0,attr(dadm,"model")()$sfun(pars[ispStop & ispGOacc & ispNR,,drop=FALSE],n_acc=n_accG))
      allLL[allok][tstopNR] <- log(gf[tstopNR] + (1-gf[tstopNR])*(1-tf[tstopNR])*pStop)
    }
   }

  # remove no response trials
  allr <- !ispNR[isp1] # used to put back into allLL
  pars <- pars[!ispNR,,drop=FALSE]
  dadm <- dadm[!ispNR,,drop=FALSE]
  isp1 <- dadm$lR==levels(dadm$lR)[1]      # 1st accumulator rows
  ispGOacc <- dadm$lI==levels(dadm$lI)[2] # Go accumulator rows
  ispStop <- is.finite(dadm$SSD) # stop-trial rows with a response
  gf <- pars[isp1,"gf"]
  tf <- pars[isp1,"tf"]

  n_trials <- nrow(dadm)/n_acc # number of trial
  trials <- 1:n_trials         # trial number
  ptrials <- rep(trials,each=n_acc)
  accST <- c(1:n_acc)[pars[1:n_acc,"lI"]==1]

  like <- numeric(n_trials)
  lds <- numeric(nrow(dadm)) # log density and survivor, used for both go and stop trials

  # Go trials with response
  if (any(!ispStop)) {
    ispGOwin <-  !ispStop & dadm$winner # Winner go accumulator rows
    tGO <- trials[!ispStop[isp1]]
    # Winner density
    lds[ispGOwin] <- log(attr(dadm,"model")()$dfunG(
      rt=dadm[ispGOwin,"rt"],pars=pars[ispGOwin,,drop=FALSE]))
    like[tGO] <- lds[ispGOwin]
    if (n_accG >1) {  # Looser survivor go accumulator(s)
      ispGOloss <- !ispStop & !dadm$winner & ispGOacc # Looser go accumulator rows
      lds[ispGOloss] <- log(1-attr(dadm,"model")()$pfunG(
         rt=dadm$rt[ispGOloss],pars=pars[ispGOloss,,drop=FALSE]))
      like[tGO] <- like[tGO] + apply(matrix(lds[ispGOloss],nrow=n_accG-1),2,sum)
    }
    like[tGO] <- (1-gf[tGO])*exp(like[tGO])
  }

  # Stop trials with a response
  if (any(ispStop)) {
    # Stop looses
    ispSwin <-       ispStop & dadm$winner              # Winner go of ST accumulator rows
    ispSlossGOacc <- ispStop & !dadm$winner & ispGOacc  # Loosing go accumulator rows
    ispSlossSTacc <- ispStop & !dadm$winner & !ispGOacc # Loosing ST accumulator rows

    # pStop at observed rt if ST present (calculate before correcting rt)
    if (n_accST>0) {
      tST <- trials[ispStop[isp1] & as.numeric(dadm$lI)[dadm$winner] == 1]
      ispST <- ptrials %in% tST
      if (any(ispST)) {
        pStop <- numeric(n_trials)
        upper <- dadm$rt[dadm$winner][tST]
        pStop[tST] <- pmax(0,attr(dadm,"model")()$sfun(pars[ispST,,drop=FALSE],
          upper=upper,n_acc=n_acc,st=c(1,1+accST)))
      }
    }

    # For following race calculations correct rt with SSD for ST accumulators
    if (any(ispSlossSTacc)) dadm[ispSlossSTacc,"rt"] <-
        dadm[ispSlossSTacc,"rt"]-dadm[ispSlossSTacc,"SSD"]

    # Fill in lds and sums over survivors for race
    lds[ispSwin] <- log(attr(dadm,"model")()$dfunG(
      rt=dadm[ispSwin,"rt"],pars=pars[ispSwin,,drop=FALSE]))
    if (n_acc >1) {  # Survivor for looser go and/or ST accumulator(s)
      lds[ispSlossGOacc | ispSlossSTacc] <- log(1-attr(dadm,"model")()$pfunG(
          rt=dadm[ispSlossGOacc | ispSlossSTacc,"rt"],
          pars=pars[ispSlossGOacc | ispSlossSTacc,,drop=FALSE]))
      # Sum survivor over loosing ST and GO accumulators
      SSTGO <- tapply(lds[ispSlossGOacc | ispSlossSTacc],
        cbind.data.frame(trials=ptrials[ispSlossGOacc | ispSlossSTacc],
                         lI=dadm[ispSlossGOacc | ispSlossSTacc,"lI"]),sum)
      SSTGO[is.na(SSTGO)] <- 0 # cases where no ST or GO survivor
    } else SSTGO <- matrix(0,ncol=2)

    # Stop accumulator survivor
    sStop <- log(1-attr(dadm,"model")()$pfunS(
      rt=dadm[ispSwin,"rt"],pars=pars[ispSwin,,drop=FALSE]))

    # Get like
    tS <- trials[ispStop[isp1]]
    # Sub-select from tS
    istSgo <- dadm$R[isp1][tS] %in% GOR
    # Stop looses or not present, can produce ST but only if wins absolute race
    like[tS] <- # no failures, works for all responses
      (1-gf[tS])*(1-tf[tS])*exp(lds[ispSwin]+SSTGO[,1]+SSTGO[,2]+sStop)
    if (any(istSgo)) like[tS][istSgo] <- like[tS][istSgo] + # tf (no stop runner), no gf, only produces GO responses
          (1-gf[tS][istSgo])*(tf[tS][istSgo])*exp(lds[ispSwin][istSgo]+SSTGO[,2][istSgo])
    # If both tf and gf then no response, handled previously

    ### Stop wins at some time before rt, and so must be ST response, add to ST race winners
    if (n_accST>0) {
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
  allLL[allok][allr] <- log(like)
  sum(pmax(min_ll,allLL[attr(dadm,"expand_winner")]))
}

