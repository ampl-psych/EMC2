my_integrate <- function(...,upper=Inf,big=10)
  # Avoids bug in integrate upper=Inf that uses only 1  subdivision
  # Use of  big=10 is arbitrary ...
{
  out <- try(integrate(...,upper=upper),silent=TRUE)
  if (!is(out,"try-error") && upper==Inf && out$subdivisions==1)
      out <- try(integrate(...,upper=big),silent=TRUE)
  out
}


log_likelihood_race_missing <- function(p_vector,dadm,min_ll=log(1e-10))
  # Race model summed log likelihood for models allowing missing values
{

  f <- function(t,p,dfun,pfun) {
    # Called by integrate to get race density for vector of times t given
    # matrix of parameters where first row is the winner.

    out <- dfun(t,
      matrix(rep(p[1,],each=length(t)),nrow=length(t),dimnames=list(NULL,dimnames(p)[[2]])))
    if (dim(p)[1]>1) for (i in 2:dim(p)[1])
      out <- out*(1-pfun(t,
        matrix(rep(p[i,],each=length(t)),nrow=length(t),dimnames=list(NULL,dimnames(p)[[2]]))))
    out
  }

  pr_pt <- function(LT,UT,ps,dadm)
    # p(untrucated response)/p(truncated response), > 1, multiplicative truncation correction
  {
    pr <- my_integrate(f,lower=0,upper=Inf,p=ps,
              dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
    if (inherits(pr, "try-error") || suppressWarnings(is.nan(pr$value))) return(NA)
    if (pr$value==0) return(0)
    pt <- my_integrate(f,lower=LT,upper=UT,p=ps,
              dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
    if (inherits(pt, "try-error") || suppressWarnings(is.nan(pt$value)) || pt$value==0) return(NA)
    out <- pmax(0,pmin(pr$value,1))/pmax(0,pmin(pt$value,1))
    if (is.infinite(out)) return(NA)
    out
  }

  pLU <- function(LT,LC,UC,UT,ps,dadm)
    # Probability from LT-LC + UC-UT
  {
    pL <- my_integrate(f,lower=LT,upper=LC,p=ps,
              dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
   if (inherits(pL,"try-error") || suppressWarnings(is.nan(pL$value))) return(NA)
   pU <- my_integrate(f,lower=UC,upper=UT,p=ps,
            dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
    if (inherits(pU,"try-error") || suppressWarnings(is.nan(pU$value))) return(NA)
   pmax(0,pmin(pL$value,1))+pmax(0,pmin(pU$value,1))
  }


  pars <- get_pars(p_vector,dadm)

  if (any(names(dadm)=="NACC")) # Some accumulators not present
    pars[as.numeric(dadm$lR)>dadm$NACC,] <- NA

  if (is.null(attr(pars,"ok")))
    ok <- !logical(dim(pars)[1]) else ok <- attr(pars,"ok")

  lds <- numeric(dim(dadm)[1]) # log pdf (winner) or survivor (losers)
  lds[dadm$winner] <- log(attr(dadm,"model")()$dfun(rt=dadm$rt[dadm$winner],
                               pars=pars[dadm$winner,]))
  n_acc <- length(levels(dadm$R))
  if (n_acc>1) lds[!dadm$winner] <-
    log(1-attr(dadm,"model")()$pfun(rt=dadm$rt[!dadm$winner],pars=pars[!dadm$winner,]))
  lds[is.na(lds) | !ok] <- -Inf

  # Calculate truncation?
  LT <- attr(dadm,"LT")
  UT <- attr(dadm,"UT")
  dotrunc <- (!is.null(LT) | !is.null(UT))
  if (is.null(LT)) LT <- 0
  if (is.null(UT)) UT <- Inf

  # Calculate censoring
  LC <- attr(dadm,"LC")
  UC <- attr(dadm,"UC")

  # Response known
  # Fast
  nort <- dadm$rt==-Inf; nort[is.na(nort)] <- FALSE; nort <- nort & !is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
          dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofix <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
                              pi <- mpars[,i,][order(!winner[,i]),]
      tmp <- my_integrate(f,lower=LT,upper=LC,p=pi,
        dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
      if ( !inherits(tmp, "try-error") && suppressWarnings(!is.nan(log(tmp$value))) )
        lds[tofix][i] <- log(pmax(0,pmin(tmp$value,1)))
    }
  }
  # Slow
  nort <- dadm$rt==Inf; nort[is.na(nort)] <- FALSE; nort <- nort & !is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
          dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofix <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
                              pi <- mpars[,i,][order(!winner[,i]),]
      tmp <- my_integrate(f,lower=UC,upper=UT,p=pi,
        dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
      if (!inherits(tmp, "try-error") && suppressWarnings(!is.nan(log(tmp$value))))
        lds[tofix][i] <- log(pmax(0,pmin(tmp$value,1)))
    }
  }
  # No direction
  nort <- is.na(dadm$rt) & !is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
          dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofix <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
                              pi <- mpars[,i,][order(!winner[,i]),]
      tmp <- my_integrate(f,lower=LT,upper=LC,p=pi,
        dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
      if ( !inherits(tmp,"try-error") && suppressWarnings(!is.nan(tmp$value)) ) {
        p <- tmp$value
        if (dim(mpars)[[1]]==1) pi <- mpars[,i,,drop=FALSE] else
                                pi <- mpars[,i,][order(!winner[,i]),]
        tmp <- my_integrate(f,lower=UC,upper=UT,p=pi,
          dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
        if ( !inherits(tmp,"try-error") && suppressWarnings(!is.nan(tmp$value)) )
          p <- p + tmp$value else p <- 0
      } else p <- 0
      lds[tofix][i] <- log(pmax(0,pmin(p,1)))
    }
  }

  # Response unknown.
  # Fast
  nort <- dadm$rt==-Inf; nort[is.na(nort)] <- FALSE; nort <- nort & is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
          dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofixfast <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
                              pi <- mpars[,i,]
      pc <- my_integrate(f,lower=LT,upper=LC,p=pi,
        dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
      if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value)))
        p <- NA else p <- pmax(0,pmin(pc$value,1))
      if (!is.na(p)) {
        if (p != 0 && !(LT==0 & UT==Inf))  cf <- pr_pt(LT,UT,mpars[,i,],dadm) else cf <- 1
        if (!is.na(cf)) p <- p*cf
      }
      if (!is.na(p) & n_acc>1) for (j in 2:n_acc) {
        pc <- my_integrate(f,lower=LT,upper=LC,p=mpars[,i,][c(j,c(1:n_acc)[-j]),],
          dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
       if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value))) {
          p <- NA; break
        }
        if (pc$value != 0 & !(LT==0 & UT==Inf))
          cf <- pr_pt(LT,UT,mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm) else cf <- 1
        if (!is.na(cf)) p <- p + pc$value*cf
      }
      lp <- log(p)
      if (!is.nan(lp) & !is.na(lp)) lds[tofixfast][i] <- lp else lds[tofixfast][i] <- -Inf
    }
  } else tofixfast <- NA
  # Slow
  nort <- dadm$rt==Inf; nort[is.na(nort)] <- FALSE; nort <- nort & is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
          dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofixslow <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
                              pi <- mpars[,i,]
      pc <- my_integrate(f,lower=UC,upper=UT,p=pi,
        dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
      if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value)))
        p <- NA else p <- pmax(0,pmin(pc$value,1))
      if (!is.na(p)) {
        if (p != 0 && !(LT==0 & UT==Inf))  cf <- pr_pt(LT,UT,mpars[,i,],dadm) else cf <- 1
        if (!is.na(cf)) p <- p*cf
      }
      if (!is.na(p) & n_acc>1) for (j in 2:n_acc) {
        pc <- my_integrate(f,lower=UC,upper=UT,p=mpars[,i,][c(j,c(1:n_acc)[-j]),],
          dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
       if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value))) {
          p <- NA; break
        }
        if (pc$value != 0 & !(LT==0 & UT==Inf))
          cf <- pr_pt(LT,UT,mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm) else cf <- 1
        if (!is.na(cf)) p <- p + pc$value*cf
      }
      lp <- log(p)
      if (!is.nan(lp) & !is.na(lp)) lds[tofixslow][i] <- lp else lds[tofixslow][i] <- -Inf
    }
  } else tofixslow <- NA
  # no direction
  nort <- is.na(dadm$rt) & is.na(dadm$R)
  nort <- nort & (pars[,"pContaminant"] == 0) # Otherwise not identifiable
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
          dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofix <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
                              pi <- mpars[,i,]
      pc <- pLU(LT,LC,UC,UT,pi,dadm)
      if (is.na(pc)) p <- NA else {
        if (pc!=0 & !(LT==0 & UT==Inf)) cf <- pr_pt(LT,UT,pi,dadm) else cf <- 1
        if (!is.na(cf)) p <- pc*cf else p <- NA
        if (!is.na(p) & n_acc>1) for (j in 2:n_acc) {
          pc <- pLU(LT,LC,UC,UT,mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm)
          if (is.na(pc)) {
            p <- NA; break
          }
          if (pc!=0 & !(LT==0 & UT==Inf))
            cf <- pr_pt(LT,UT,mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm) else cf <- 1
          if (is.na(cf)) {
            p <- NA; break
          } else p <- p +  pc*cf
        }
      }
      lp <- log(pmax(0,pmin(p,1)))
      if (!is.nan(lp) & !is.na(lp)) lds[tofix][i] <- lp
    }
  }


  # Truncation where not censored or censored and response known
  ok <- is.finite(lds[attr(dadm,"unique_nort") & dadm$winner])
  alreadyfixed <- is.na(dadm$R[attr(dadm,"unique_nort") & dadm$winner])
  ok <- ok & !alreadyfixed
  if ( dotrunc & any(ok) ) {
    tpars <- pars[attr(dadm,"unique_nort"),,drop=FALSE]
    tpars <- array(tpars[,,drop=FALSE],dim=c(n_acc,nrow(tpars)/n_acc,ncol(tpars)),
          dimnames = list(NULL,NULL,colnames(tpars)))[,,,drop=FALSE]
    winner <- matrix(dadm$winner[attr(dadm,"unique_nort")],nrow=n_acc)[,,drop=FALSE]
    cf <- rep(NA,length(ok))
    for (i in 1:length(ok)) if (ok[i]) {
      if (dim(tpars)[[1]]==1) pi <- t(as.matrix(tpars[,i,])) else
                              pi <- tpars[,i,][order(!winner[,i]),]
      cf[i] <- pr_pt(LT,UT,pi,dadm)
    }
    cf <- rep(log(cf),each=n_acc)[attr(dadm,"expand_nort")]
    fix <- dadm$winner & !is.na(cf) & !is.nan(cf) & is.finite(cf)
    if (any(fix)) lds[fix] <- lds[fix] + cf[fix]
    badfix <- dadm$winner & (is.na(cf) | is.nan(cf) | is.infinite(cf))
    if (!all(is.na(tofixfast))) badfix <- badfix & !tofixfast
    if (!all(is.na(tofixslow))) badfix <- badfix & !tofixslow
    if (any(badfix)) lds[badfix] <- -Inf
  }

  # Non-process (contaminant) miss.
  ispContaminant <- pars[,"pContaminant"]>0
  if ( any(ispContaminant) ) {
    p <- exp(lds[dadm$winner])
    pc <- pars[dadm$winner,"pContaminant"]
    isMiss <- is.na(dadm$R[dadm$winner])
    p[isMiss] <- pc[isMiss] + (1-pc[isMiss])*p[isMiss]
    p[!isMiss] <- (1-pc[!isMiss])*p[!isMiss]
    lds[dadm$winner] <- log(p)
  }

  lds <- lds[attr(dadm,"expand")] # decompress
  if (n_acc>1) {
    winner <- dadm$winner[attr(dadm,"expand")]
    ll <- lds[winner]
    if (n_acc==2) {
      ll <- ll + lds[!winner]
    } else {
      ll <- ll + apply(matrix(lds[!winner],nrow=n_acc-1),2,sum)
    }
  } else ll <- lds
  ll[is.na(ll) | is.nan(ll)] <- -Inf

  return(sum(pmax(min_ll,ll)))
}


log_likelihood_race_missing_LBAU <- function(p_vector,dadm,min_ll=log(1e-10))
  # Race model summed log likelihood for an LBA allowing negative rates and
  # missing values
{

  f <- function(t,p,dfun,pfun) {
    # Called by integrate to get race density for vector of times t given
    # matrix of parameters where first row is the winner.

    out <- dfun(t,
      matrix(rep(p[1,],each=length(t)),nrow=length(t),dimnames=list(NULL,dimnames(p)[[2]])))
    if (dim(p)[1]>1) for (i in 2:dim(p)[1])
      out <- out*(1-pfun(t,
        matrix(rep(p[i,],each=length(t)),nrow=length(t),dimnames=list(NULL,dimnames(p)[[2]]))))
    out
  }

  pr_pt <- function(LT,UT,ps,dadm)
    # p(untrucated response)/p(truncated response), > 1, multiplicative truncation correction
  {
    pr <- my_integrate(f,lower=0,upper=Inf,p=ps,
              dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
    if (inherits(pr, "try-error") || suppressWarnings(is.nan(pr$value))) return(NA)
    if (pr$value==0) return(0)
    pt <- my_integrate(f,lower=LT,upper=UT,p=ps,
              dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
    if (inherits(pt, "try-error") || suppressWarnings(is.nan(pt$value)) || pt$value==0) return(NA)
    out <- pmax(0,pmin(pr$value,1))/pmax(0,pmin(pt$value,1))
    if (is.infinite(out)) return(NA)
    out
  }

  pLU <- function(LT,LC,UC,UT,ps,dadm)
    # Probability from LT-LC + UC-UT
  {
    pL <- my_integrate(f,lower=LT,upper=LC,p=ps,
              dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
   if (inherits(pL,"try-error") || suppressWarnings(is.nan(pL$value))) return(NA)
   pU <- my_integrate(f,lower=UC,upper=UT,p=ps,
            dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
    if (inherits(pU,"try-error") || suppressWarnings(is.nan(pU$value))) return(NA)
   pmax(0,pmin(pL$value,1))+pmax(0,pmin(pU$value,1))
  }


  pars <- get_pars(p_vector,dadm)

  if (any(names(dadm)=="NACC")) # Some accumulators not present
    pars[as.numeric(dadm$lR)>dadm$NACC,] <- NA

  if (is.null(attr(pars,"ok")))
    ok <- !logical(dim(pars)[1]) else ok <- attr(pars,"ok")

  lds <- numeric(dim(dadm)[1]) # log pdf (winner) or survivor (losers)
  lds[dadm$winner] <- log(attr(dadm,"model")()$dfun(rt=dadm$rt[dadm$winner],
                               pars=pars[dadm$winner,]))
  n_acc <- length(levels(dadm$R))
  if (n_acc>1) lds[!dadm$winner] <-
    log(1-attr(dadm,"model")()$pfun(rt=dadm$rt[!dadm$winner],pars=pars[!dadm$winner,]))
  lds[is.na(lds) | !ok] <- -Inf

  # Calculate truncation?
  LT <- attr(dadm,"LT")
  UT <- attr(dadm,"UT")
  dotrunc <- (!is.null(LT) | !is.null(UT))
  if (is.null(LT)) LT <- 0
  if (is.null(UT)) UT <- Inf

  # Calculate censoring
  LC <- attr(dadm,"LC")
  UC <- attr(dadm,"UC")

  # Response known
  # Fast
  nort <- dadm$rt==-Inf; nort[is.na(nort)] <- FALSE; nort <- nort & !is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
          dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofix <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
                              pi <- mpars[,i,][order(!winner[,i]),]
      tmp <- my_integrate(f,lower=LT,upper=LC,p=pi,
        dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
      if ( !inherits(tmp, "try-error") && suppressWarnings(!is.nan(log(tmp$value))) )
        lds[tofix][i] <- log(pmax(0,pmin(tmp$value,1)))
    }
  }
  # Slow
  nort <- dadm$rt==Inf; nort[is.na(nort)] <- FALSE; nort <- nort & !is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
          dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofix <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
                              pi <- mpars[,i,][order(!winner[,i]),]
      tmp <- my_integrate(f,lower=UC,upper=UT,p=pi,
        dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
      # if (!inherits(tmp, "try-error") && suppressWarnings(!is.nan(log(tmp$value))))
      #   lds[tofix][i] <- log(pmax(0,pmin(tmp$value,1)))
      if (!inherits(tmp, "try-error") && suppressWarnings(!is.nan(log(tmp$value))))
        lds[tofix][i] <- pmax(0,pmin(tmp$value,1))
      lds[tofix][i] <- log(lds[tofix][i] + (1-lds[tofix][i])*prod(pnorm(0,mpars[,i,"v"],mpars[,i,"sv"])))
    }
  }
  # No direction
  nort <- is.na(dadm$rt) & !is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
          dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofix <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
                              pi <- mpars[,i,][order(!winner[,i]),]
      tmp <- my_integrate(f,lower=LT,upper=LC,p=pi,
        dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
      if ( !inherits(tmp,"try-error") && suppressWarnings(!is.nan(tmp$value)) ) {
        p <- tmp$value
        if (dim(mpars)[[1]]==1) pi <- mpars[,i,,drop=FALSE] else
                                pi <- mpars[,i,][order(!winner[,i]),]
        tmp <- my_integrate(f,lower=UC,upper=UT,p=pi,
          dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
        if ( !inherits(tmp,"try-error") && suppressWarnings(!is.nan(tmp$value)) )
          p <- p + tmp$value else p <- 0
      } else p <- 0
      lds[tofix][i] <- log(pmax(0,pmin(p,1)))
    }
  }

  # Response unknown.
  # Fast
  nort <- dadm$rt==-Inf; nort[is.na(nort)] <- FALSE; nort <- nort & is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
          dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofixfast <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
                              pi <- mpars[,i,]
      pc <- my_integrate(f,lower=LT,upper=LC,p=pi,
        dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
      if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value)))
        p <- NA else p <- pmax(0,pmin(pc$value,1))
      if (!is.na(p)) {
        if (p != 0 && !(LT==0 & UT==Inf))  cf <- pr_pt(LT,UT,mpars[,i,],dadm) else cf <- 1
        if (!is.na(cf)) p <- p*cf
      }
      if (!is.na(p) & n_acc>1) for (j in 2:n_acc) {
        pc <- my_integrate(f,lower=LT,upper=LC,p=mpars[,i,][c(j,c(1:n_acc)[-j]),],
          dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
       if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value))) {
          p <- NA; break
        }
        if (pc$value != 0 & !(LT==0 & UT==Inf))
          cf <- pr_pt(LT,UT,mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm) else cf <- 1
        if (!is.na(cf)) p <- p + pc$value*cf
      }
      lp <- log(p)
      if (!is.nan(lp) & !is.na(lp)) lds[tofixfast][i] <- lp else lds[tofixfast][i] <- -Inf
    }
  } else tofixfast <- NA
  # Slow
  nort <- dadm$rt==Inf; nort[is.na(nort)] <- FALSE; nort <- nort & is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
          dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofixslow <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
                              pi <- mpars[,i,]
      pc <- my_integrate(f,lower=UC,upper=UT,p=pi,
        dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
      if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value)))
        p <- NA else p <- pmax(0,pmin(pc$value,1))
      if (!is.na(p)) {
        if (p != 0 && !(LT==0 & UT==Inf))  cf <- pr_pt(LT,UT,mpars[,i,],dadm) else cf <- 1
        if (!is.na(cf)) p <- p*cf
      }
      if (!is.na(p) & n_acc>1) for (j in 2:n_acc) {
        pc <- my_integrate(f,lower=UC,upper=UT,p=mpars[,i,][c(j,c(1:n_acc)[-j]),],
          dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun)
       if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value))) {
          p <- NA; break
        }
        if (pc$value != 0 & !(LT==0 & UT==Inf))
          cf <- pr_pt(LT,UT,mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm) else cf <- 1
        if (!is.na(cf)) p <- p + pc$value*cf
      }
      lp <- log(p + (1-p)*prod(pnorm(0,mpars[,i,"v"],mpars[,i,"sv"])))
      if (!is.nan(lp) & !is.na(lp)) lds[tofixslow][i] <- lp else lds[tofixslow][i] <- -Inf
    }
  } else tofixslow <- NA
  # no direction
  nort <- is.na(dadm$rt) & is.na(dadm$R)
  nort <- nort & (pars[,"pContaminant"] == 0) # Otherwise not identifiable
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
          dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofix <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      if (dim(mpars)[[1]]==1) pi <- t(as.matrix(mpars[,i,])) else
                              pi <- mpars[,i,]
      pc <- pLU(LT,LC,UC,UT,pi,dadm)
      if (is.na(pc)) p <- NA else {
        if (pc!=0 & !(LT==0 & UT==Inf)) cf <- pr_pt(LT,UT,pi,dadm) else cf <- 1
        if (!is.na(cf)) p <- pc*cf else p <- NA
        if (!is.na(p) & n_acc>1) for (j in 2:n_acc) {
          pc <- pLU(LT,LC,UC,UT,mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm)
          if (is.na(pc)) {
            p <- NA; break
          }
          if (pc!=0 & !(LT==0 & UT==Inf))
            cf <- pr_pt(LT,UT,mpars[,i,][c(j,c(1:n_acc)[-j]),],dadm) else cf <- 1
          if (is.na(cf)) {
            p <- NA; break
          } else p <- p +  pc*cf
        }
      }
      lp <- log(pmax(0,pmin(p,1)))
      if (!is.nan(lp) & !is.na(lp)) lds[tofix][i] <- lp
    }
  }


  # Truncation where not censored or censored and response known
  ok <- is.finite(lds[attr(dadm,"unique_nort") & dadm$winner])
  alreadyfixed <- is.na(dadm$R[attr(dadm,"unique_nort") & dadm$winner])
  ok <- ok & !alreadyfixed
  if ( dotrunc & any(ok) ) {
    tpars <- pars[attr(dadm,"unique_nort"),,drop=FALSE]
    tpars <- array(tpars[,,drop=FALSE],dim=c(n_acc,nrow(tpars)/n_acc,ncol(tpars)),
          dimnames = list(NULL,NULL,colnames(tpars)))[,,,drop=FALSE]
    winner <- matrix(dadm$winner[attr(dadm,"unique_nort")],nrow=n_acc)[,,drop=FALSE]
    cf <- rep(NA,length(ok))
    for (i in 1:length(ok)) if (ok[i]) {
      if (dim(tpars)[[1]]==1) pi <- t(as.matrix(tpars[,i,])) else
                              pi <- tpars[,i,][order(!winner[,i]),]
      cf[i] <- pr_pt(LT,UT,pi,dadm)
    }
    cf <- rep(log(cf),each=n_acc)[attr(dadm,"expand_nort")]
    fix <- dadm$winner & !is.na(cf) & !is.nan(cf) & is.finite(cf)
    if (any(fix)) lds[fix] <- lds[fix] + cf[fix]
    badfix <- dadm$winner & (is.na(cf) | is.nan(cf) | is.infinite(cf))
    if (!all(is.na(tofixfast))) badfix <- badfix & !tofixfast
    if (!all(is.na(tofixslow))) badfix <- badfix & !tofixslow
    if (any(badfix)) lds[badfix] <- -Inf
  }

  # Non-process (contaminant) miss.
  ispContaminant <- pars[,"pContaminant"]>0
  if ( any(ispContaminant) ) {
    p <- exp(lds[dadm$winner])
    pc <- pars[dadm$winner,"pContaminant"]
    isMiss <- is.na(dadm$R[dadm$winner])
    p[isMiss] <- pc[isMiss] + (1-pc[isMiss])*p[isMiss]
    p[!isMiss] <- (1-pc[!isMiss])*p[!isMiss]
    lds[dadm$winner] <- log(p)
  }

  lds <- lds[attr(dadm,"expand")] # decompress
  if (n_acc>1) {
    winner <- dadm$winner[attr(dadm,"expand")]
    ll <- lds[winner]
    if (n_acc==2) {
      ll <- ll + lds[!winner]
    } else {
      ll <- ll + apply(matrix(lds[!winner],nrow=n_acc-1),2,sum)
    }
  } else ll <- lds
  ll[is.na(ll) | is.nan(ll)] <- -Inf

  return(sum(pmax(min_ll,ll)))
}

log_likelihood_mt <- function(p_vector,dadm,min_ll=log(1e-10))
  # Multiple threshold (BE or TC) summed log likelihood
  # attr(dadm,"dL")
{

  pars <- get_pars(p_vector,dadm)

  Dnams <- dimnames(pars)[[2]][substr(dimnames(pars)[[2]],1,2)=="DT"]
  tmats <- array(cbind(rep(0,dim(pars)[1]),pars[,c(Dnams,"b")]),dim=c(2,nrow(pars)/2,2+length(Dnams)),
                dimnames=list(NULL,NULL,c("DT0",Dnams,"b")))
  pnams <- dimnames(pars)[[2]][!(dimnames(pars)[[2]] %in% c(Dnams,"b"))]
  pmats <- array(pars[,pnams],dim=c(2,nrow(pars)/2,length(pnams)),
                dimnames=list(NULL,NULL,pnams))
  nt <- dim(pmats)[2]
  like <- numeric(nt)
  for (i in 1:nt) {
    i2 <- i*2
    # Get look up table for current rating (dadm$R[i2])
    pick <- attr(dadm,"dL")[as.numeric(dadm$R[i2]),]
    tmp <- try(n1PDF_MTR_1(rt=dadm$rt[i*2], pars = pmats[c(pick[1],pick[3]),i,],
      dl = tmats[pick[3],i,pick[4]], du = tmats[pick[5],i,pick[6]], b = tmats[pick[1],i,pick[2]]),silent=TRUE)
    if (inherits(tmp,"try-error")) like[i] <- 0 else like[i] <- tmp
  }
  like <- log(like)
  like[is.na(like) | is.nan(like) | like == Inf] <- 0
  return(sum(pmax(min_ll,like[attr(dadm,"expand_winner")])))
}


log_likelihood_mt <- function(p_vector,dadm,min_ll=log(1e-10),n_cores=10)
  # Multiple threshold (BE or TC) summed log likelihood
  # attr(dadm,"dL")
{

  mt <- function(i,dadm,tmats) {
    i2 <- i*2
    # Get look up table for current rating (dadm$R[i2])
    pick <- attr(dadm,"dL")[as.numeric(dadm$R[i2]),]
    tmp <- try(n1PDF_MTR_1(rt=dadm$rt[i*2], pars = pmats[c(pick[1],pick[3]),i,],
      dl = tmats[pick[3],i,pick[4]], du = tmats[pick[5],i,pick[6]],
      b = tmats[pick[1],i,pick[2]]),silent=TRUE)
    if (inherits(tmp,"try-error")) 0 else tmp
  }

  pars <- get_pars(p_vector,dadm)

  Dnams <- dimnames(pars)[[2]][substr(dimnames(pars)[[2]],1,2)=="DT"]
  tmats <- array(cbind(rep(0,dim(pars)[1]),pars[,c(Dnams,"b")]),dim=c(2,nrow(pars)/2,2+length(Dnams)),
                dimnames=list(NULL,NULL,c("DT0",Dnams,"b")))
  pnams <- dimnames(pars)[[2]][!(dimnames(pars)[[2]] %in% c(Dnams,"b"))]
  pmats <- array(pars[,pnams],dim=c(2,nrow(pars)/2,length(pnams)),
                dimnames=list(NULL,NULL,pnams))
  nt <- dim(pmats)[2]
  like <- log(unlist(mclapply(1:nt,mt,dadm=dadm,tmats=tmats,mc.cores=n_cores)))
  like[is.na(like) | is.nan(like) | like == Inf] <- 0
  return(sum(pmax(min_ll,like[attr(dadm,"expand_winner")])))
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
  lds[is.na(lds) | !ok] <- min_ll
  lds <- lds[attr(dadm,"expand")] # decompress
  winner <- dadm$winner[attr(dadm,"expand")]
  if (n_acc>1) {
    ll <- lds[winner]
    if (n_acc==2)
      ll <- ll + lds[!winner] else
        ll <- ll + apply(matrix(lds[!winner],nrow=n_acc-1),2,sum)
    ll[is.na(ll)] <- min_ll
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


