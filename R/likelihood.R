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

log_likelihood_race_ss <- function(pars,dadm,model,min_ll=log(1e-10))
{

  # # DELETE me when sure no longer happening
  # # spurious go winners on no-response trials
  # dadm$winner[is.na(dadm$R)] <- FALSE
  # if (any(is.infinite(dadm$rt))) stop("BUGGER!")

  # All bad?
  if (is.null(attr(pars,"ok")))
    ok <- !logical(dim(pars)[1]) else ok <- attr(pars,"ok")
  if (!any(ok)) return(min_ll*length(attr(dadm, "expand_winner")))

  # Nomenclature for indices:
  # "is" = logical, "isp" pars/dadm index,
  # "t" trials index, "ist" logical on trial index
  # "n_" number of integer

  # Counts
  n_acc <- length(levels(dadm$lR))                   # total number of accumulators
  n_trials <- nrow(dadm)/n_acc                       # number of trials
  n_accG <- sum(as.numeric(dadm[1:n_acc,"lI"])==2)   # go accumulators
  n_accST <- sum(as.numeric(dadm[1:n_acc,"lI"])==1)

  # Likelihood for all trials and for ok trials
  allLL <- rep(-Inf,n_trials)
  # used to put results into allLL[allok]
  allok <- ok[dadm$lR==levels(dadm$lR)[1]]

  # Remove trials not ok
  pars <- pars[ok,,drop=FALSE]
  dadm <- dadm[ok,,drop=FALSE]

  # Only keep OK trials for stop trial computations
  n_trials <- nrow(dadm)/n_acc # number of trials
  trials <- 1:n_trials         # trial number
  # ptrials <- rep(trials,each=n_acc) # trial number for pars/dadm

  # Booleans for ok pars/dadm
  isp1 <- dadm$lR==levels(dadm$lR)[1]      # 1st accumulator rows
  ispGOacc <- dadm$lI==levels(dadm$lI)[2]  # Go accumulator rows
  ispStop <- is.finite(dadm$SSD)           # Stop-trial rows

  # Failure parameters for each trial
  gf <- pars[isp1,"gf"]
  tf <- pars[isp1,"tf"]

  # Go response names
  GoR <- as.character(dadm[1:n_acc,"lR"][dadm[1:n_acc,"lI"]==2])

  # No response
  ispNR <- is.na(dadm$R)
  if ( any(ispNR) ) {  # Note by definition no ST present

    # Go trial with no response
    ispgoNR <- ispNR & !ispStop            # boolean
    tgoNR <- c(1:sum(isp1))[ispgoNR[isp1]] # trial number

    if (any(ispgoNR)) # Fill in go failure
      allLL[allok][tgoNR] <- log(gf[tgoNR])

    # Stop trial with no response
    ispstopNR <- ispNR & ispStop        # boolean
    # ptrials[ispstopNR & isp1]
    if ( any(ispstopNR) ) { # Stop trial probability
        # Non-response and stop trial & go/stop accumulator parameters
        pStop <- pmin(1,pmax(0,  # protection to keep in 0-1
          attr(dadm,"model")()$sfun(
            pars[ispNR & ispStop & ispGOacc,,drop=FALSE],n_acc=n_accG)
        ))
      # Fill in stop-trial non-response probabilities, either 1) go failure
      # 2) not go failure and not trigger failure and stop wins
      tstopNR <- trials[ispstopNR[isp1]]  # trial number
      allLL[allok][tstopNR] <- log(gf[tstopNR] + (1-gf[tstopNR])*(1-tf[tstopNR])*pStop)
    }
  }

  # Response made
  if (any(!ispNR)) {
    # Only keep response trials for further computation
    allr <- !ispNR[isp1] # used to put back into allLL[allok]
    pars <- pars[!ispNR,,drop=FALSE]
    dadm <- dadm[!ispNR,,drop=FALSE]
    n_trials <- nrow(dadm)/n_acc # number of trials
    trials <- 1:n_trials
    ptrials <- rep(trials,each=n_acc) # trial number for pars/dadm

    isp1 <- dadm$lR==levels(dadm$lR)[1]      # 1st accumulator rows
    ispGOacc <- dadm$lI==levels(dadm$lI)[2]  # Go accumulator rows
    ispStop <- is.finite(dadm$SSD)           # stop-trial rows with a response
    gf <- pars[isp1,"gf"]
    tf <- pars[isp1,"tf"]

    # likelihoods for trials with a response (eventually, intermediately probabilities)
    like <- numeric(n_trials)
    lds <- numeric(nrow(dadm)) # log density and survivor, used for both go and stop trials

    # Go trials with response
    if (any(!ispStop)) {
      ispGOwin <-  !ispStop & dadm$winner # Winner go accumulator rows
      # tGO <- ptrials[ispGOwin]
      tGO <- trials[!ispStop[isp1]]       # Go trials
      # Winner density
      # lds[ispGOwin] <- log(attr(dadm,"model")()$dfunG(
      #   rt=dadm[ispGOwin,"rt"],pars=pars[ispGOwin,,drop=FALSE]))
      # like[tGO] <- lds[ispGOwin]
      like[tGO] <- log(attr(dadm,"model")()$dfunG(
        rt=dadm[ispGOwin,"rt"],pars=pars[ispGOwin,,drop=FALSE]))

      if (n_accG >1) {  # Looser survivor go accumulator(s)
        ispGOloss <- !ispStop & !dadm$winner & ispGOacc # Looser go accumulator rows
        # lds[ispGOloss] <- log(1-attr(dadm,"model")()$pfunG(
        #    rt=dadm$rt[ispGOloss],pars=pars[ispGOloss,,drop=FALSE]))
        # like[tGO] <- like[tGO] + apply(matrix(lds[ispGOloss],nrow=n_accG-1),2,sum)
        like[tGO] <- like[tGO] + apply(matrix(log(1-attr(dadm,"model")()$pfunG(
          rt=dadm$rt[ispGOloss],pars=pars[ispGOloss,,drop=FALSE])),nrow=n_accG-1),2,sum)
      }
      # Transform back to densities to include go failure
      like[tGO] <- (1-gf[tGO])*exp(like[tGO])
    }

    # Stop trials with a response, occurs if
    # 1) All triggered
    #   a) Stop does not beat go before rt, produces go and ST
    #   b) Stop beats go before rt, produces ST only
    # 2) go failure, stop triggered, produces ST only
    # 3) go triggered, stop failure, produces go only
    # NB: cant have both go and stop failure as then no response

    if (any(ispStop)) {                       # Stop trials
      ispSGO <- ispStop & (dadm$R %in% GoR)   # Go responses
      ispSST <- ispStop & !(dadm$R %in% GoR)

      # Go beats stop and ST (if any)
      if (any(ispSGO)) {

        ispGOwin <-  ispSGO & dadm$winner # Winner go accumulator rows
        tGO <- ptrials[ispGOwin]          # Go trials
        like[tGO] <- log(attr(dadm,"model")()$dfunG(
          rt=dadm[ispGOwin,"rt"],pars=pars[ispGOwin,,drop=FALSE])) +
            log(1-attr(dadm,"model")()$pfunS(
                rt=dadm[ispGOwin,"rt"],pars=pars[ispGOwin,,drop=FALSE]))
        if (n_acc > 1) {  # Looser survivor accumulator(s)
          ispGOloss <- ispSGO & !dadm$winner
          like[tGO] <- like[tGO] + apply(matrix(log(1-attr(dadm,"model")()$pfunG(
            rt=dadm$rt[ispGOloss],pars=pars[ispGOloss,,drop=FALSE])),nrow=n_acc-1),2,sum)
        }
        # Transform back to densities to include go failure
        like[tGO] <- (1-gf[tGO])*exp(like[tGO])
      }

      # ST WINS (never tf)
      if (any(ispSST)) { # Go triggered 1a) with (1-pStop), all race,
                         #              1b) with pStop, only ST races
                         # Go failure    2) Only ST races
        # Stop probability on ST win trials
        pStop <- attr(dadm,"model")()$sfun(
            pars[ispSST,,drop=FALSE],n_acc=n_accG)
        # ST win ll
        ispSSTwin <-  dadm$winner &  ispSST
        tST <- ptrials[ispSSTwin]
        like[tST] <- log(attr(dadm,"model")()$dfunG(
          rt=dadm[ispSSTwin,"rt"]-dadm[ispSSTwin,"SSD"], # correct ST racers for SSD delay
          pars=pars[ispSSTwin,,drop=FALSE]))
        # ST survivor
        if (n_accST > 1) {  # Survivor for looser for ST accumulator(s)
          ispSSTloss <-  !dadm$winner &  ispSST
          llST <-  log(1-attr(dadm,"model")()$pfunG(
                rt=dadm[ispSSTloss,"rt"],pars=pars[ispSSTloss,,drop=FALSE]))
            if (n_accST == 2) like[tST] <- like[tST] + llST else
              like[tST] <- like[tST] + apply(matrix(llST,nrow=n_accS-1),2,sum)
        }
        # Go looser survivor
        ispSGloss <- ispSST & !dadm$winner # Looser go accumulator rows
        llG <- apply(matrix(log(1-attr(dadm,"model")()$pfunG(
            rt=dadm[ispSGloss,"rt"],pars=pars[ispSGloss,,drop=FALSE])),
          nrow=n_accG),2,sum)
        like[tST] <- (1-tf[tST])*(                 # Never trigger failure
              gf[tST]*exp(like[tST]) +             # Case 2
          (1-gf[tST])*(pStop*exp(like[tST]) +      # Case 1b
                    (1-pStop)*exp(like[tST]+llG))) # Case 1a
      }
    }
    allLL[allok][allr] <- log(like)
  }
  sum(pmax(min_ll,allLL[attr(dadm,"expand_winner")]))
}

