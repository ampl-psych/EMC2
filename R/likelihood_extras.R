log_likelihood_race_missing <- function(p_vector,dadm,min_ll=log(1e-10))
  # Race model summed log likelihood
{

  f <- function(t,p)
    # Called by integrate to get race density for vector of times t given
    # matrix of parameters where first row is the winner.
  {
    if (!is.matrix(p)) p <- matrix(p,nrow=1,dimnames=list(NULL,names(p)))
    out <- attr(dadm,"model")()$dfun(t,
                                   matrix(rep(p[1,],each=length(t)),nrow=length(t),
                                          dimnames=list(NULL,dimnames(p)[[2]])))
    if (dim(p)[1]>1) for (i in 2:dim(p)[1])
      out <- out*(1-attr(dadm,"model")()$pfun(t,
        matrix(rep(p[i,],each=length(t)),nrow=length(t),dimnames=list(NULL,dimnames(p)[[2]]))))
    out
  }

  logTrunc_race <- function(LT,UT,dadm,pars)
    # log truncation probability for each rt
  {
    if ( (LT > 0) | (is.finite(UT)) )
      # Get truncation adjustment. same length as winner
    {
      upars <- pars[attr(dadm,"unique_nort"),]
      looser <- matrix(!dadm[attr(dadm,"unique_nort"),"winner"],nrow=n_acc)
      np <- dim(upars)[1]/n_acc
      upars <- array(upars,dim=c(n_acc,np,dim(upars)[2]),
                     dimnames=list(NULL,NULL,dimnames(upars)[[2]]))
      lps <- numeric(np)
      for (i in 1:np) {
        tmp <- try(integrate(f,lower=LT,upper=UT,p=upars[,i,][order(looser[,i]),]),silent=TRUE)
        if (inherits(tmp, "try-error") || tmp$value == 0 ) lps[i] <- Inf else
          lps[i] <- log(tmp$value)
      }
      lps[is.na(lps) | is.nan(lps)] <- Inf
      return(lps[attr(dadm,"expand_nort")])
    } else return(0)
  }

  # Censoring
  logCens_race <- function(LT,UT,dadm,pars,n_trials,n_acc,upper=TRUE)
    # b = censor bound, bT = truncation bound
  {
    if (upper)
      b <- attr(dadm,"UC")[1] else
        b <- attr(dadm,"LC")[1]
      out <- numeric(n_trials)
      if (is.null(b)) return(out)
      if (upper) ok <- dadm$rt==Inf else ok <- dadm$rt==-Inf
      if (!any(ok)) return(out)
      d <- dadm[ok,]
      p <- pars[ok,]
      d[,"winner"][is.na(d[,"R"])] <- NA
      looser <- matrix(!d[,"winner"],nrow=n_acc)
      np <- dim(p)[1]/n_acc
      p <- array(p,dim=c(n_acc,np,dim(p)[2]),
                 dimnames=list(NULL,NULL,dimnames(p)[[2]]))
      if (upper) {
        lo <- b
        up <- UT
      } else {
        lo <- LT
        up <- b
      }
      lps <- numeric(np)
      for (i in 1:np) {
        if ( all(is.na(looser[,i])) ) {
          l <- !logical(length(looser[,i]))
          num <- den <- vector(mode="list",length=np)
          for (j in 1:length(looser[,i])) {
            temp <- l; temp[j] <- FALSE
            num[[j]] <- try(integrate(f,lower=lo,upper=up,p=p[,i,][order(temp),]),silent=TRUE)
            if (inherits(num[[j]], "try-error")) {num <- NA; break}
            den[[j]] <- try(integrate(f,lower=LT,upper=UT,p=p[,i,][order(temp),]),silent=TRUE)
            if (inherits(den[[j]],"try-error") || den[[j]]$value==0) {den <- NA; break}
          }
          if (any(is.na(num)) || any(is.na(den))) lps[i] <- -Inf else {
            num <- unlist(lapply(num,function(x)x$value))
            den <- unlist(lapply(den,function(x)x$value))
            nd <- num/den
            w <- num/sum(num)
            lps[i] <- log(sum(w*nd))
          }
        } else {
          num <- try(integrate(f,lower=lo,upper=up,p=p[,i,][order(looser[,i]),]),silent=TRUE)
          if (inherits(num, "try-error")) lps[i] <- -Inf else {
            den <- try(integrate(f,lower=LT,upper=UT,p=p[,i,][order(looser[,i]),]),silent=TRUE)
            if (inherits(den, "try-error") || den$value == 0) lps[i] <- -Inf else
              lps[i] <- log(num$value/den$value)
          }
        }
      }
      lps[is.na(lps) | is.nan(lps)] <- -Inf
      if (upper) return(c(0,lps)[attr(dadm,"expand_uc")]) else
        return(c(0,lps)[attr(dadm,"expand_lc")])
  }

  pars <- get_pars(p_vector,dadm)

  # Density for winner rows
  lds <- numeric(dim(dadm)[1]) # log pdf (winner) or survivor (losers)
  lds[attr(dadm,"ok_dadm_winner")] <- log(attr(dadm,"model")()$dfun(
    rt=dadm$rt[attr(dadm,"ok_dadm_winner")],pars=pars[attr(dadm,"ok_dadm_winner"),]))

  n_acc <- length(levels(dadm$R))                # number of accumulators
  n_trials <- length(attr(dadm,"expand_winner")) # number of trials

  # Truncation
  LT <- attr(dadm,"LT")[1]
  UT <- attr(dadm,"UT")[1]
  if (is.null(LT)) LT <- 0
  if (is.null(UT)) UT <- Inf
  # log truncation
  lt <- logTrunc_race(LT,UT,dadm,pars)
  if (is.character(lt)) return(min_ll*length(attr(dadm,"expand")))
  if (length(lt)>1) lt[!attr(dadm,"ok_trial")] <- 0
  # log upper censor
  luc <- logCens_race(LT,UT,dadm,pars,n_trials,n_acc)
  if (is.character(luc)) return(min_ll*length(attr(dadm,"expand")))
  # log lower censor
  llc <- logCens_race(LT,UT,dadm,pars,n_trials,n_acc,upper=FALSE)
  if (is.character(llc)) return(min_ll*length(attr(dadm,"expand")))

  # Survivors and truncation
  if (n_acc>1) {
    lds[attr(dadm,"ok_dadm_looser")] <- log(1-attr(dadm,"model")()$pfun(
      rt=dadm$rt[attr(dadm,"ok_dadm_looser")],pars=pars[attr(dadm,"ok_dadm_looser"),]))
  }
  lds <- lds[attr(dadm,"expand")]  # decompress
  lds[is.na(lds)] <- 0
  if (n_acc>1) {
    ll <- numeric(n_trials)
    ll[attr(dadm,"ok_trials")] <- lds[attr(dadm,"ok_da_winner")]
    if (n_acc==2)
      ll[attr(dadm,"ok_trials")] <- ll[attr(dadm,"ok_trials")] + lds[attr(dadm,"ok_da_looser")] else
        ll[attr(dadm,"ok_trials")] <- ll[attr(dadm,"ok_trials")] +
      apply(matrix(lds[attr(dadm,"ok_da_looser")],nrow=n_acc-1),2,sum)
    ll <- ll - lt + luc + llc
    ll[is.na(ll)] <- 0
    sum(pmax(min_ll,ll))
  } else sum(pmax(min_ll,lds - lt + luc + llc))
}



log_likelihood_race_ss <- function(p_vector,dadm,min_ll=log(1e-10))
  # Race model summed log likelihood, one accumulator (lR=="stop") is ExGaussian
  # others are Wald.
{

  pars <- get_pars(p_vector,dadm)

  if (is.null(attr(pars,"ok")))
    ok <- !logical(dim(pars)[1]) else ok <- attr(pars,"ok")

  # standard race, ignores R=NA
  isstop <- dadm[,"lR"]=="stop"
  stopwinner <- dadm$winner & isstop
  gowinner <- dadm$winner & !isstop
  stoplooser <- !dadm$winner & isstop
  golooser <- !dadm$winner & !isstop
  lds <- numeric(dim(dadm)[1]) # log pdf (winner) or survivor (losers)
  lds[gowinner] <- log(attr(dadm,"model")()$dfunG(rt=dadm$rt[gowinner],pars=pars[gowinner,,drop=FALSE]))
  lds[stopwinner] <- log(attr(dadm,"model")()$dfunS(rt=dadm$rt[stopwinner],pars=pars[stopwinner,,drop=FALSE]))
  n_acc <- length(levels(dadm$R))
  if (n_acc>1) {
    lds[golooser] <- log(1-attr(dadm,"model")()$pfunG(rt=dadm$rt[golooser],pars=pars[golooser,,drop=FALSE]))
    lds[stoplooser] <- log(1-attr(dadm,"model")()$pfunS(rt=dadm$rt[stoplooser],pars=pars[stoplooser,,drop=FALSE]))
  }
  lds[is.na(lds) | !ok] <- 0
  lds <- lds[attr(dadm,"expand")] # decompress
  winner <- dadm$winner[attr(dadm,"expand")]
  if (n_acc>1) {
    ll <- lds[winner]
    if (n_acc==2)
      ll <- ll + lds[!winner] else
        ll <- ll + apply(matrix(lds[!winner],nrow=n_acc-1),2,sum)
    ll[is.na(ll)] <- 0
  } else ll <- lds
  like <- exp(ll)

  # stop success
  stopsucess_accs <- is.finite(dadm$SSD[attr(dadm,"expand")]) & (dadm$R[attr(dadm,"expand")] == "stop")
  stoptrial <- is.finite(dadm$SSD[attr(dadm,"expand")][winner])
  goresp <- dadm$R[attr(dadm,"expand")][winner] != "stop"
  if (any(stopsucess_accs)) {
    parstop <- pars[attr(dadm,"expand"),][stopsucess_accs,]
    iss <- stoptrial & !goresp
    like[iss] <- attr(dadm,"model")()$sfun(parstop,n_acc)
    like[iss][is.na(like[iss])] <- 0
  }

  # trigger failures
  tf <- pars[attr(dadm,"expand"),"tf"][winner]
  dotf <- goresp & stoptrial & (tf > 0)
  if (any(dotf)) {
    tflike <- exp(apply(matrix(lds,nrow=n_acc)[-1,dotf,drop=FALSE],2,sum))
    like[dotf] <- tf[dotf]*tflike + (1-tf[dotf])*like[dotf]
    stopsucess <- !dotf & stoptrial
    like[stopsucess] <- like[stopsucess]*(1-tf[stopsucess])
  }

  # go failures
  gf <- pars[,"gf"][attr(dadm,"expand")][winner]
  if (any(gf>0)) {
    like[!goresp] <- gf[!goresp] + (1-gf[!goresp])*like[!goresp]
    like[goresp] <- like[goresp]*(1-gf[goresp])
  }
  like[like<0 | is.na(like)] <- 0
  sum(pmax(min_ll,log(like)))
}


