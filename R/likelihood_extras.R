# log_likelihood_race_missing <- function(p_vector,dadm,min_ll=log(1e-10))
#   # Race model summed log likelihood
# {
#
#   f <- function(t,p)
#     # Called by integrate to get race density for vector of times t given
#     # matrix of parameters where first row is the winner.
#   {
#     if (!is.matrix(p)) p <- matrix(p,nrow=1,dimnames=list(NULL,names(p)))
#     out <- attr(dadm,"model")()$dfun(t,
#                                    matrix(rep(p[1,],each=length(t)),nrow=length(t),
#                                           dimnames=list(NULL,dimnames(p)[[2]])))
#     if (dim(p)[1]>1) for (i in 2:dim(p)[1])
#       out <- out*(1-attr(dadm,"model")()$pfun(t,
#         matrix(rep(p[i,],each=length(t)),nrow=length(t),dimnames=list(NULL,dimnames(p)[[2]]))))
#     out
#   }
#
#   logTrunc_race <- function(LT,UT,dadm,pars)
#     # log truncation probability for each rt
#   {
#     if ( (LT > 0) | (is.finite(UT)) )
#       # Get truncation adjustment. same length as winner
#     {
#       upars <- pars[attr(dadm,"unique_nort"),]
#       looser <- matrix(!dadm[attr(dadm,"unique_nort"),"winner"],nrow=n_acc)
#       np <- dim(upars)[1]/n_acc
#       upars <- array(upars,dim=c(n_acc,np,dim(upars)[2]),
#                      dimnames=list(NULL,NULL,dimnames(upars)[[2]]))
#       lps <- numeric(np)
#       for (i in 1:np) {
#         tmp <- try(integrate(f,lower=LT,upper=UT,p=upars[,i,][order(looser[,i]),]),silent=TRUE)
#         if (inherits(tmp, "try-error") || tmp$value == 0 ) lps[i] <- Inf else
#           lps[i] <- log(tmp$value)
#       }
#       lps[is.na(lps) | is.nan(lps)] <- Inf
#       return(lps[attr(dadm,"expand_nort")])
#     } else return(0)
#   }
#
#   # Censoring
#   logCens_race <- function(LT,UT,dadm,pars,n_trials,n_acc,upper=TRUE)
#     # b = censor bound, bT = truncation bound
#   {
#     if (upper)
#       b <- attr(dadm,"UC")[1] else
#         b <- attr(dadm,"LC")[1]
#       out <- numeric(n_trials)
#       if (is.null(b)) return(out)
#       if (upper) ok <- dadm$rt==Inf else ok <- dadm$rt==-Inf
#       if (!any(ok)) return(out)
#       d <- dadm[ok,]
#       p <- pars[ok,]
#       d[,"winner"][is.na(d[,"R"])] <- NA
#       looser <- matrix(!d[,"winner"],nrow=n_acc)
#       np <- dim(p)[1]/n_acc
#       p <- array(p,dim=c(n_acc,np,dim(p)[2]),
#                  dimnames=list(NULL,NULL,dimnames(p)[[2]]))
#       if (upper) {
#         lo <- b
#         up <- UT
#       } else {
#         lo <- LT
#         up <- b
#       }
#       lps <- numeric(np)
#       for (i in 1:np) {
#         if ( all(is.na(looser[,i])) ) {
#           l <- !logical(length(looser[,i]))
#           num <- den <- vector(mode="list",length=np)
#           for (j in 1:length(looser[,i])) {
#             temp <- l; temp[j] <- FALSE
#             num[[j]] <- try(integrate(f,lower=lo,upper=up,p=p[,i,][order(temp),]),silent=TRUE)
#             if (inherits(num[[j]], "try-error")) {num <- NA; break}
#             den[[j]] <- try(integrate(f,lower=LT,upper=UT,p=p[,i,][order(temp),]),silent=TRUE)
#             if (inherits(den[[j]],"try-error") || den[[j]]$value==0) {den <- NA; break}
#           }
#           if (any(is.na(num)) || any(is.na(den))) lps[i] <- -Inf else {
#             num <- unlist(lapply(num,function(x)x$value))
#             den <- unlist(lapply(den,function(x)x$value))
#             nd <- num/den
#             w <- num/sum(num)
#             lps[i] <- log(sum(w*nd))
#           }
#         } else {
#           num <- try(integrate(f,lower=lo,upper=up,p=p[,i,][order(looser[,i]),]),silent=TRUE)
#           if (inherits(num, "try-error")) lps[i] <- -Inf else {
#             den <- try(integrate(f,lower=LT,upper=UT,p=p[,i,][order(looser[,i]),]),silent=TRUE)
#             if (inherits(den, "try-error") || den$value == 0) lps[i] <- -Inf else
#               lps[i] <- log(num$value/den$value)
#           }
#         }
#       }
#       lps[is.na(lps) | is.nan(lps)] <- -Inf
#       if (upper) return(c(0,lps)[attr(dadm,"expand_uc")]) else
#         return(c(0,lps)[attr(dadm,"expand_lc")])
#   }
#
#   pars <- get_pars(p_vector,dadm)
#
#   # Density for winner rows
#   lds <- numeric(dim(dadm)[1]) # log pdf (winner) or survivor (losers)
#   lds[attr(dadm,"ok_dadm_winner")] <- log(attr(dadm,"model")()$dfun(
#     rt=dadm$rt[attr(dadm,"ok_dadm_winner")],pars=pars[attr(dadm,"ok_dadm_winner"),]))
#
#   n_acc <- length(levels(dadm$R))                # number of accumulators
#   n_trials <- length(attr(dadm,"expand_winner")) # number of trials
#
#   # Truncation
#   LT <- attr(dadm,"LT")[1]
#   UT <- attr(dadm,"UT")[1]
#   if (is.null(LT)) LT <- 0
#   if (is.null(UT)) UT <- Inf
#   # log truncation
#   lt <- logTrunc_race(LT,UT,dadm,pars)
#   if (is.character(lt)) return(min_ll*length(attr(dadm,"expand")))
#   if (length(lt)>1) lt[!attr(dadm,"ok_trial")] <- 0
#   # log upper censor
#   luc <- logCens_race(LT,UT,dadm,pars,n_trials,n_acc)
#   if (is.character(luc)) return(min_ll*length(attr(dadm,"expand")))
#   # log lower censor
#   llc <- logCens_race(LT,UT,dadm,pars,n_trials,n_acc,upper=FALSE)
#   if (is.character(llc)) return(min_ll*length(attr(dadm,"expand")))
#
#   # Survivors and truncation
#   if (n_acc>1) {
#     lds[attr(dadm,"ok_dadm_looser")] <- log(1-attr(dadm,"model")()$pfun(
#       rt=dadm$rt[attr(dadm,"ok_dadm_looser")],pars=pars[attr(dadm,"ok_dadm_looser"),]))
#   }
#   lds <- lds[attr(dadm,"expand")]  # decompress
#   lds[is.na(lds)] <- 0
#   if (n_acc>1) {
#     ll <- numeric(n_trials)
#     ll[attr(dadm,"ok_trials")] <- lds[attr(dadm,"ok_da_winner")]
#     if (n_acc==2)
#       ll[attr(dadm,"ok_trials")] <- ll[attr(dadm,"ok_trials")] + lds[attr(dadm,"ok_da_looser")] else
#         ll[attr(dadm,"ok_trials")] <- ll[attr(dadm,"ok_trials")] +
#       apply(matrix(lds[attr(dadm,"ok_da_looser")],nrow=n_acc-1),2,sum)
#     ll <- ll - lt + luc + llc
#     ll[is.na(ll)] <- 0
#     sum(pmax(min_ll,ll))
#   } else sum(pmax(min_ll,lds - lt + luc + llc))
# }


#' Title
#'
#' @param p_vector
#' @param dadm
#' @param min_ll
#'
#' @return
#' @export
#'
#' @examples
log_likelihood_race_missing <- function(p_vector,dadm,min_ll=log(1e-10))
  # Race model summed log likelihood
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
    pr <- try(integrate(f,lower=0,upper=Inf,p=ps,
              dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun),silent=TRUE)
    if (inherits(pr, "try-error") || suppressWarnings(is.nan(pr$value))) return(NA)
    if (pr$value==0) return(0)
    pt <- try(integrate(f,lower=LT,upper=UT,p=ps,
              dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun),silent=TRUE)
    if (inherits(pt, "try-error") || suppressWarnings(is.nan(pt$value)) || pt$value==0) return(NA)
    out <- pmax(0,pmin(pr$value,1))/pmax(0,pmin(pt$value,1))
    if (is.infinite(out)) return(NA)
    out
  }

  pLU <- function(LT,LC,UC,UT,ps,dadm)
    # Probability from LT-LC + UC-UT
  {
    pL <- try(integrate(f,lower=LT,upper=LC,p=ps,
              dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun),silent=TRUE)
   if (inherits(pL,"try-error") || suppressWarnings(is.nan(pL$value))) return(NA)
   pU <- try(integrate(f,lower=UC,upper=UT,p=ps,
            dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun),silent=TRUE)
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
      tmp <- try(integrate(f,lower=LT,upper=LC,p=mpars[,i,][order(!winner[,i]),],
        dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun),silent=TRUE)
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
      tmp <- try(integrate(f,lower=UC,upper=UT,p=mpars[,i,][order(!winner[,i]),],
        dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun),silent=TRUE)
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
      tmp <- try(integrate(f,lower=LT,upper=LC,p=mpars[,i,][order(!winner[,i]),],
        dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun),silent=TRUE)
      if ( !inherits(tmp,"try-error") && suppressWarnings(!is.nan(tmp$value)) ) {
        p <- tmp$value
        tmp <- try(integrate(f,lower=UC,upper=UT,p=mpars[,i,][order(!winner[,i]),],
          dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun),silent=TRUE)
        if ( !inherits(tmp,"try-error") && suppressWarnings(!is.nan(tmp$value)) )
          p <- p + tmp$value else p <- 0
      } else p <- 0
      lds[tofix][i] <- log(pmax(0,pmin(p,1)))
    }
  }

  # Response unknown
  # Fast
  nort <- dadm$rt==-Inf; nort[is.na(nort)] <- FALSE; nort <- nort & is.na(dadm$R)
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
          dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofixfast <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      pc <- try(integrate(f,lower=LT,upper=LC,p=mpars[,i,],
        dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun),silent=TRUE)
      if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value)))
        p <- NA else p <- pmax(0,pmin(pc$value,1))
      if (!is.na(p)) {
        if (p != 0 && !(LT==0 & UT==Inf))  cf <- pr_pt(LT,UT,mpars[,i,],dadm) else cf <- 1
        if (!is.na(cf)) p <- p*cf
      }
      if (!is.na(p) & n_acc>1) for (j in 2:n_acc) {
        pc <- try(integrate(f,lower=LT,upper=LC,p=mpars[,i,][c(j,c(1:n_acc)[-j]),],
          dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun),silent=TRUE)
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
      pc <- try(integrate(f,lower=UC,upper=UT,p=mpars[,i,],
        dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun),silent=TRUE)
      if (inherits(pc, "try-error") || suppressWarnings(is.nan(pc$value)))
        p <- NA else p <- pmax(0,pmin(pc$value,1))
      if (!is.na(p)) {
        if (p != 0 && !(LT==0 & UT==Inf))  cf <- pr_pt(LT,UT,mpars[,i,],dadm) else cf <- 1
        if (!is.na(cf)) p <- p*cf
      }
      if (!is.na(p) & n_acc>1) for (j in 2:n_acc) {
        pc <- try(integrate(f,lower=UC,upper=UT,p=mpars[,i,][c(j,c(1:n_acc)[-j]),],
          dfun=attr(dadm,"model")()$dfun,pfun=attr(dadm,"model")()$pfun),silent=TRUE)
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
  if ( any(nort) ) {
    mpars <- array(pars[nort,,drop=FALSE],dim=c(n_acc,sum(nort)/n_acc,ncol(pars)),
          dimnames = list(NULL,NULL,colnames(pars)))
    winner <- matrix(dadm$winner[nort],nrow=n_acc)
    tofix <- dadm$winner & nort
    for (i in 1:dim(mpars)[2]) {
      pc <- pLU(LT,LC,UC,UT,mpars[,i,],dadm)
      if (is.na(pc)) p <- NA else {
        if (pc!=0 & !(LT==0 & UT==Inf)) cf <- pr_pt(LT,UT,mpars[,i,],dadm) else cf <- 1
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
    tpars <- pars[attr(dadm,"unique_nort"),]
    tpars <- array(tpars[,,drop=FALSE],dim=c(n_acc,nrow(tpars)/n_acc,ncol(tpars)),
          dimnames = list(NULL,NULL,colnames(tpars)))[,,,drop=FALSE]
    winner <- matrix(dadm$winner[attr(dadm,"unique_nort")],nrow=n_acc)[,,drop=FALSE]
    cf <- rep(NA,length(ok))
    for (i in 1:length(ok)) if (ok[i]) {
      cf[i] <- pr_pt(LT,UT,tpars[,i,][order(!winner[,i]),],dadm)
    }
    cf <- rep(log(cf),each=n_acc)[attr(dadm,"expand_nort")]
    fix <- dadm$winner & !is.na(cf) & !is.nan(cf) & is.finite(cf)
    if (any(fix)) lds[fix] <- lds[fix] + cf[fix]
    badfix <- dadm$winner & (is.na(cf) | is.nan(cf) | is.infinite(cf))
    if (!all(is.na(tofixfast))) badfix <- badfix & !tofixfast
    if (!all(is.na(tofixslow))) badfix <- badfix & !tofixslow
    if (any(badfix)) lds[badfix] <- -Inf
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
  }
  ll[is.na(ll) | is.nan(ll)] <- -Inf

  return(sum(pmax(min_ll,ll)))
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


