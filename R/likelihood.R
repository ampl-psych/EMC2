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
  if (any(names(dadm)=="RACE")){# Some accumulators not present
    pars[as.numeric(dadm$lR)>as.numeric(as.character(dadm$RACE)),] <- NA
  }

  if (is.null(attr(pars,"ok"))){
    ok <- !logical(dim(pars)[1])
  } else ok <- attr(pars,"ok")
  lds <- numeric(dim(dadm)[1]) # log pdf (winner) or survivor (losers)
  lds[dadm$winner] <- log(model$dfun(rt=dadm$rt[dadm$winner],
                                                    pars=pars[dadm$winner,]))
  n_acc <- length(levels(dadm$R))
  if (n_acc>1) lds[!dadm$winner] <- log(1-model$pfun(rt=dadm$rt[!dadm$winner],pars=pars[!dadm$winner,]))
  lds[is.na(lds) | !ok] <- min_ll
  if (n_acc>1) {
    ll <- lds[dadm$winner]
    if (n_acc==2) {
      ll <- ll + lds[!dadm$winner]
    } else {
      ll <- ll + apply(matrix(lds[!dadm$winner],nrow=n_acc-1),2,sum)
    }
    ll[is.na(ll)] <- min_ll
    return(sum(pmax(min_ll,ll[attr(dadm,"expand")])))
  } else return(sum(pmax(min_ll,lds[attr(dadm,"expand")])))
}


log_likelihood_race_missing <- function(pars,dadm,model,min_ll=log(1e-10))
  # Race model summed log likelihood for models allowing missing values
{

  if (any(names(dadm)=="RACE")) # Some accumulators not present
    pars[as.numeric(dadm$lR)>as.numeric(as.character(dadm$RACE)),] <- NA

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

  if (n_acc>1) {
    ll <- lds[dadm$winner]
    if (n_acc==2) {
      ll <- ll + lds[!dadm$winner]
    } else {
      ll <- ll + apply(matrix(lds[!dadm$winner],nrow=n_acc-1),2,sum)
    }
  } else ll <- lds
  ll[is.na(ll) | is.nan(ll)] <- -Inf

  return(sum(pmax(min_ll,ll[attr(dadm,"expand_winner")])))
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

log_likelihood_sdt <- function(pars,dadm, model,lb=-Inf, min_ll=log(1e-10))
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
  # log probability
  ll <- numeric(sum(dadm$winner))
  if (!is.null(attr(pars,"ok"))) { # Bad parameter region
    ok <- attr(pars,"ok")
    okw <- ok[dadm$winner]
    ll[ok] <- log(model$pfun(lt=lt[okw],ut=ut[okw],pars=pars[dadm$winner & ok,,drop=FALSE]))
  } else ll <- log(model$pfun(lt=lt,ut=ut,pars=pars[dadm$winner,,drop=FALSE]))
  ll <- ll[attr(dadm,"expand")]
  ll[is.na(ll)] <- 0
  sum(pmax(min_ll,ll))
}

# Two options:
# 1) component = NULL in which case we do all likelihoods in one block
# 2) component = integer, in which case we are blocking the ll and only want that one
log_likelihood_joint <- function(proposals, dadms, model_list, component = NULL){
  parPreFixs <- unique(gsub("[|].*", "", colnames(proposals)))
  i <- 0
  k <- 0
  total_ll <- 0
  for(dadm in dadms){
    i <- i + 1
    if(is.null(component) || component == i){
      if(is.data.frame(dadm)){
        k <- k + 1
        # Sometimes designs are the same across models (typically in fMRI)
        # Instead of storing multiple, we just store a pointer to the first original design
        if(is.numeric(attr(dadm, "designs"))){
          ref_idx <- attr(dadm, "designs")
          attr(dadm, "designs") <- attr(dadms[[ref_idx]], "designs")
        }
        parPrefix <- parPreFixs[k]
        # Unfortunately indexing can't get quicker than this as far as I can tell.
        columns_to_use <- sapply(strsplit(colnames(proposals), "|", fixed = TRUE), function(x) x == parPrefix)[1,]
        currentPars <- proposals[,columns_to_use, drop = F]
        colnames(currentPars) <- gsub(".*[|]", "", colnames(currentPars))
        total_ll <- total_ll +  calc_ll_manager(currentPars, dadm, model_list[[i]])
      }
    }
  }
  return(total_ll)
}

####  Stop-signal ----

# log(1 - x)
log1m <- function(x) {
  ifelse(is.na(x), NA_real_,
         ifelse(x >= 1, -Inf,
                ifelse(x == -Inf, 0,
                       log1p(-x))))
}

# log(1 - exp(x))
log1m_exp <- function(x) {
  ifelse(is.na(x), NA_real_,
         ifelse(x > 0, NA_real_,  # 1 - exp(x) < 0 for x > 0
                ifelse(x == 0, -Inf,  # log(1 - 1) = log(0) = -Inf
                       ifelse(x == -Inf, 0,  # log(1 - 0) = log(1) = 0
                              ifelse(x > -log(2),  # x > -log(2)
                                     log(-expm1(x)),
                                     log1m(exp(x)))))))
}

# log(exp(a) + exp(b))
log_sum_exp <- function(a, b = NULL) {
  if (is.null(b)) {
    # Vector version: log(sum(exp(x)))
    x <- a
    if (length(x) == 0) return(-Inf)
    if (length(x) == 1) return(x[1])

    max_val <- max(x, na.rm = TRUE)
    if (is.infinite(max_val) && max_val < 0) return(-Inf)

    # Subtract max for numerical stability
    exp_terms <- exp(x - max_val)
    exp_terms[is.infinite(x) & x < 0] <- 0  # Handle -Inf cases

    return(max_val + log(sum(exp_terms, na.rm = TRUE)))
  } else {
    # Two-argument version: log(exp(a) + exp(b))
    ifelse(a == -Inf, b,
           ifelse(b == -Inf, a,
                  ifelse(a > b,
                         a + log1p(exp(b - a)),
                         b + log1p(exp(a - b)))))
  }
}

# log(exp(a) - exp(b))
log_diff_exp <- function(a, b) {
  # Handle edge cases
  result <- ifelse(is.na(a) | is.na(b), NA_real_,
                   ifelse(is.infinite(a) & a > 0, NA_real_,  # +Inf case
                          ifelse(a < b, NA_real_,  # log of negative
                                 ifelse(a == b, -Inf,  # log(0)
                                        ifelse(b == -Inf, a, NA_real_)))))  # Will be overwritten

  # For valid cases where a > b and b != -Inf
  valid <- !is.na(result) & a > b & is.finite(b)

  if (any(valid)) {
    diff <- b[valid] - a[valid]

    # For numerical stability
    stable_case <- diff > -log(2)  # exp(b)/exp(a) > 0.5

    result[valid] <- ifelse(stable_case,
                            a[valid] + log1m(exp(diff)),
                            a[valid] + log(-expm1(diff)))
  }

  # Handle b == -Inf case
  result[b == -Inf & !is.na(a)] <- a[b == -Inf & !is.na(a)]

  return(result)
}

# log mixture density: log(theta * exp(lambda1) + (1-theta) * exp(lambda2))
log_mix <- function(theta, lambda1, lambda2) {
  # Input validation and edge cases
  result <- ifelse(is.na(theta) | is.na(lambda1) | is.na(lambda2), NA_real_,
                   ifelse(theta < 0 | theta > 1, NA_real_,
                          ifelse(theta == 0, lambda2,
                                 ifelse(theta == 1, lambda1, NA_real_))))  # Will be overwritten

  # Handle cases where both lambdas are -Inf
  both_neginf <- lambda1 == -Inf & lambda2 == -Inf
  result[both_neginf & !is.na(result)] <- -Inf

  # Handle cases where one lambda is -Inf
  lambda1_neginf <- lambda1 == -Inf & lambda2 != -Inf & !both_neginf
  lambda2_neginf <- lambda2 == -Inf & lambda1 != -Inf & !both_neginf

  result[lambda1_neginf & !is.na(result)] <- log1m(theta[lambda1_neginf]) + lambda2[lambda1_neginf]
  result[lambda2_neginf & !is.na(result)] <- log(theta[lambda2_neginf]) + lambda1[lambda2_neginf]

  # General case: neither lambda is -Inf and theta is in (0,1)
  general_case <- !is.na(result) & theta > 0 & theta < 1 &
    is.finite(lambda1) & is.finite(lambda2)

  if (any(general_case)) {
    result[general_case] <- log_sum_exp(
      log(theta[general_case]) + lambda1[general_case],
      log1m(theta[general_case]) + lambda2[general_case]
    )
  }

  return(result)
}

# wrapper function for numerical integration of race likelihood function
# retries with arbitrary large but finite upper limit, if initially failed with Inf upper
my.integrate <- function(..., upper = Inf, big = 10) {
  out <- try(
    integrate(..., upper = upper, rel.tol = 1e-6, abs.tol = 1e-8),
    silent = TRUE
  )
  if (inherits(out, "try-error")) {
    return(0)
  }
  if (upper == Inf && out$subdivisions == 1) {
    out <- try(
      integrate(..., upper = big, rel.tol = 1e-6, abs.tol = 1e-8),
      silent = TRUE
    )
    if (inherits(out, "try-error") || out$subdivisions == 1) {
      return(0)
    }
  }
  return(out$value)
}


log_likelihood_race_ss <- function(pars, dadm, model, min_ll = log(1e-10)) {
  # All bad?
  if (is.null(attr(pars, "ok"))) {
    ok <- !logical(dim(pars)[1])
  } else {
    ok <- attr(pars, "ok")
  }
  if (!any(ok)) {
    return(min_ll * length(attr(dadm, "expand")))
  }

  # Nomenclature for indices:
  # "is" = logical, "isp" pars/dadm index,
  # "t" trials index, "ist" logical on trial index
  # "n_" number of integer

  # Counts
  n_acc <- length(levels(dadm$lR))                    # total number of accumulators
  n_trials <- nrow(dadm) / n_acc                      # number of trials
  n_accG <- sum(as.numeric(dadm[1:n_acc, "lI"]) == 2) # go accumulators
  n_accST <- sum(as.numeric(dadm[1:n_acc, "lI"]) == 1)

  # Likelihood for all trials and for ok trials
  allLL <- rep(min_ll, n_trials)
  # used to put results into allLL[allok]
  allok <- ok[dadm$lR == levels(dadm$lR)[1]]

  # Remove trials not ok
  pars <- pars[ok, , drop = FALSE]
  dadm <- dadm[ok, , drop = FALSE]

  # Only keep OK trials for stop trial computations
  n_trials <- nrow(dadm) / n_acc    # number of trials
  trials <- 1:n_trials              # trial number

  # Booleans for ok pars/dadm
  isp1 <- dadm$lR == levels(dadm$lR)[1]     # 1st accumulator rows
  ispGOacc <- dadm$lI == levels(dadm$lI)[2] # Go accumulator rows
  ispStop <- is.finite(dadm$SSD)            # Stop-trial rows

  # Failure parameters for each trial
  gf <- pars[isp1, "gf"]
  tf <- pars[isp1, "tf"]

  # Go response names
  GoR <- as.character(dadm[1:n_acc, "lR"][dadm[1:n_acc, "lI"] == 2])

  # No response
  ispNR <- is.na(dadm$R)
  if (any(ispNR)) {  # Note by definition no ST present

    # Go failures
    ispgoNR <- ispNR & !ispStop            # No response and no stop signal
    tgoNR <- c(1:sum(isp1))[ispgoNR[isp1]] # trial number
    if (any(ispgoNR)) {
      # go trial with no response: log likelihood is simply log go failure prob
      allLL[allok][tgoNR] <- log(gf[tgoNR])
    }

    # Stop trial with no response and no ST accumulator
    ispstopNR <- ispNR & ispStop & (n_accST == 0)
    if (any(ispstopNR)) {
      # Stop trial probability
      # Non-response and stop trial & go/stop accumulator parameters
      # protection to keep in 0-1
      pStop <- pmin(
        1,
        pmax(
          0,
          model$sfun(
            pars = pars[ispNR & ispStop & ispGOacc, , drop = FALSE],
            n_acc = n_accG
            # ,min_ll = min_ll
          )
        )
      )
      # Fill in stop-trial non-response probabilities, either 1) go failure
      # 2) not go failure and not trigger failure and stop wins
      tstopNR <- trials[ispstopNR[isp1]]  # trial number
      # likelihood = gf + [(1-gf) x (1-tf) x pStop]
      stop_success_lprob <- log1m(gf[tstopNR]) + log1m(tf[tstopNR]) + log(pStop)
      allLL[allok][tstopNR] <- log_sum_exp(log(gf[tstopNR]), stop_success_lprob)
    }
  } # end of bracket for trials with no response

  # Response made
  if (any(!ispNR)) {
    # Only keep response trials for further computation
    allr <- !ispNR[isp1] # used to put back into allLL[allok]
    pars <- pars[!ispNR, , drop = FALSE]
    dadm <- dadm[!ispNR, , drop = FALSE]
    n_trials <- nrow(dadm) / n_acc # number of trials
    trials <- 1:n_trials
    ptrials <- rep(trials, each = n_acc) # trial number for pars/dadm

    isp1 <- dadm$lR == levels(dadm$lR)[1]      # 1st accumulator rows
    ispGOacc <- dadm$lI == levels(dadm$lI)[2]  # Go accumulator rows
    ispStop <- is.finite(dadm$SSD)             # stop-trial rows with a response
    gf <- pars[isp1, "gf"]
    tf <- pars[isp1, "tf"]

    # likelihoods for trials with a response (eventually, intermediately probabilities)
    loglike <- numeric(n_trials)
    lds <- numeric(nrow(dadm)) # log density and survivor, used for both go and stop trials

    # Go trials with response
    if (any(!ispStop)) {
      ispGOwin <- !ispStop & dadm$winner # Winner go accumulator rows
      tGO <- ptrials[ispGOwin]  # Go trials
      # Winner density
      loglike[tGO] <- log(model$dfunG(
        rt = dadm$rt[ispGOwin],
        pars = pars[ispGOwin, , drop = FALSE]
        #,log.d = TRUE,min_ll = min_ll
      ))
      if (n_accG > 1) {
        # Looser survivor go accumulator(s)
        ispGOloss <- !ispStop & !dadm$winner & ispGOacc # Looser go accumulator rows
        # add (summed) log survivor probability to winner log density
        loglike[tGO] <- loglike[tGO] + apply(
          X = matrix(
            data = log(1-model$pfunG(
              rt = dadm$rt[ispGOloss],
              pars = pars[ispGOloss, , drop = FALSE]
              # lower.tail = FALSE,log.p = TRUE,min_ll = min_ll
            )),
            nrow = n_accG - 1
          ),
          MARGIN = 2,
          FUN = sum
        )
      }
      # likelihood = (1-gf) x go_race_prob
      loglike[tGO] <- log1m(gf[tGO]) + loglike[tGO]
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
      ispSST <- ispStop & !(dadm$R %in% GoR)  # ST responses
      # Go beats stop and ST (if any)
      if (any(ispSGO)) {
        ispGOwin <- ispSGO & dadm$winner # Winner go accumulator rows
        tGO <- ptrials[ispGOwin]          # Go trials
        loglike[tGO] <- log(model$dfunG(
          rt = dadm$rt[ispGOwin],
          pars = pars[ispGOwin, , drop = FALSE]
          #,log.d = TRUE,min_ll = min_ll
        ))
        if (n_accG > 1) {
          # Loser survivor gp accumulators
          ispGOloss <- ispSGO & !dadm$winner & ispGOacc
          # add (summed) log survivor probability to winner log density
          loglike[tGO] <- loglike[tGO] + apply(
            X = matrix(
              data = log(1-model$pfunG(
                rt = dadm$rt[ispGOloss],
                pars = pars[ispGOloss, , drop = FALSE]
                #,lower.tail = FALSE,log.p = TRUE,min_ll = min_ll
              )),
              nrow = n_accG - 1
            ),
            MARGIN = 2,
            FUN = sum
          )
        }
        # trigger stop, add in stop survivor
        ts <- loglike[tGO] + log(1-model$pfunS(
          rt = dadm$rt[ispGOwin],
          pars = pars[ispGOwin, , drop = FALSE]
          #,lower.tail = FALSE,log.p = TRUE,min_ll = min_ll
        ))
        # ST loosers
        if (n_accST == 0) {
          stl <- 0
        } else {
          ispSTloss <- ispSGO & !ispGOacc
          stl <- apply(
            X = matrix(
              data = log(1-model$pfunG(
                rt = dadm$rt[ispSTloss] - pars[ispSTloss, "SSD"], # correct for SSD
                pars=pars[ispSTloss, , drop = FALSE]
                # lower.tail = FALSE,log.p = TRUE,min_ll = min_ll
              )),
              nrow = n_accST
            ),
            MARGIN = 2,
            FUN = sum
          )
        }
        # likelihood = (1-gf) x [tf x go_prob + (1-tf) x stop_fail_prob]
        loglike[tGO] <- log1m(gf[tGO]) + log_mix(tf[tGO], loglike[tGO], ts + stl)
      }

      # ST WINS (never tf)
      if (any(ispSST)) { # Go triggered 1a) with (1-pStop), all race,
                         #              1b) with pStop, only ST races
                         # Go failure    2) Only ST races
        # ST winner rows
        ispSSTwin <- dadm$winner & ispSST
        # Stop probability on ST win trials, only needs go accumulators
        pStop <- model$sfun(
          pars = pars[ispSST & ispGOacc, , drop = FALSE],
          n_acc = n_accG,
          upper = dadm$rt[ispSSTwin]
          # ,min_ll = min_ll
        ) # pStop before observed RT
        # ST win ll
        tST <- ptrials[ispSSTwin]
        loglike[tST] <- log(model$dfunG(
          rt = dadm$rt[ispSSTwin] - dadm$SSD[ispSSTwin], # correct ST racers for SSD delay
          pars = pars[ispSSTwin, , drop = FALSE]
          #,log.d = TRUE,min_ll = min_ll
        ))
        # ST looser survivors
        if (n_accST > 1) {
          # Survivor for loser for ST accumulator(s)
          ispSSTloss <- !dadm$winner & ispSST & !ispGOacc
          llST <- log(1-model$pfunG(
            rt = dadm$rt[ispSSTloss] - dadm$SSD[ispSSTloss],
            pars = pars[ispSSTloss, , drop = FALSE]
            #,lower.tail = FALSE,log.p = TRUE,min_ll = min_ll
          ))
          if (n_accST == 2) {
            # Could remove branch, maybe faster as no matrix sum?
            loglike[tST] <- loglike[tST] + llST
          } else {
            loglike[tST] <- loglike[tST] + apply(
              X = matrix(data = llST, nrow = n_accST - 1),
              MARGIN = 2,
              FUN = sum
            )
          }
        }
        # Go loser survivor
        ispSGloss <- ispSST & ispGOacc
        llG <- apply(
          X = matrix(
            data = log(1-model$pfunG(
              rt = dadm$rt[ispSGloss],
              pars = pars[ispSGloss,,drop=FALSE]
              #,lower.tail = FALSE,log.p = TRUE,min_ll = min_ll
            )),
            nrow = n_accG
          ),
          MARGIN = 2,
          FUN = 2
        )

        # like[tST] <- (1 - tf[tST]) * (            # Never trigger failure
        #   gf[tST] * exp(like[tST]) +              # Case 2, gf only ST race
        #   (1 - gf[tST]) * (                       # Case 1, no gf...
        #     pStop * exp(like[tST]) +              # Case 1b, no gf, stop beats go, ST race
        #     (1 - pStop) * exp(like[tST] + llG)    # Case 1a, no gf, all race (no stop win)
        #   )
        # )
        loglike[tST] <- log1m(tf[tST]) + log_mix(
          gf[tST],
          loglike[tST],
          log_mix(pStop, loglike[tST], loglike[tST] + llG)
        )
      }
    }

    # end of bracket for trials with a response; assign log likelihoods to
    # overall log likelihood vector
    allLL[allok][allr] <- loglike
  }

  # protect against non-finity / tiny log likelihoods
  allLL[is.na(allLL) | is.nan(allLL)] <- min_ll
  allLL <- pmax(min_ll, allLL)
  # de-compress
  allLL <- allLL[attr(dadm, "expand")]

  llR <<- allLL

  # return summed log likelihood
  return(sum(allLL))

}


log_likelihood_race_ss_old <- function(pars,dadm,model,min_ll=log(1e-10))
{
  # All bad?
  if (is.null(attr(pars,"ok")))
    ok <- !logical(dim(pars)[1]) else ok <- attr(pars,"ok")
    if (!any(ok)) return(min_ll*length(attr(dadm, "expand")))

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
    allLL <- rep(min_ll,n_trials)
    # used to put results into allLL[allok]
    allok <- ok[dadm$lR==levels(dadm$lR)[1]]

    # Remove trials not ok
    pars <- pars[ok,,drop=FALSE]
    dadm <- dadm[ok,,drop=FALSE]

    # Only keep OK trials for stop trial computations
    n_trials <- nrow(dadm)/n_acc # number of trials
    trials <- 1:n_trials         # trial number

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

      # Go failures
      ispgoNR <- ispNR & !ispStop            # No response and no stop signal
      tgoNR <- c(1:sum(isp1))[ispgoNR[isp1]] # trial number
      if (any(ispgoNR))
        allLL[allok][tgoNR] <- log(gf[tgoNR])

      # Stop trial with no response and no ST accumulator
      ispstopNR <- ispNR & ispStop & (n_accST == 0)
      if ( any(ispstopNR) ) { # Stop trial probability
        # Non-response and stop trial & go/stop accumulator parameters
        pStop <- pmin(1,pmax(0,  # protection to keep in 0-1
                             model$sfun(
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
        tGO <- ptrials[ispGOwin]  # Go trials
        # Winner density
        like[tGO] <- log(model$dfunG(
          rt=dadm$rt[ispGOwin],pars=pars[ispGOwin,,drop=FALSE]))
        if (n_accG >1) {  # Looser survivor go accumulator(s)
          ispGOloss <- !ispStop & !dadm$winner & ispGOacc # Looser go accumulator rows
          like[tGO] <- like[tGO] + apply(matrix(log(1-model$pfunG(
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
        ispSST <- ispStop & !(dadm$R %in% GoR)  # ST responses
        # Go beats stop and ST (if any)
        if (any(ispSGO)) {
          ispGOwin <-  ispSGO & dadm$winner # Winner go accumulator rows
          tGO <- ptrials[ispGOwin]          # Go trials
          like[tGO] <- log(model$dfunG(
            rt=dadm$rt[ispGOwin],pars=pars[ispGOwin,,drop=FALSE]))
          if (n_accG > 1) {  # Looser survivor gp accumulators
            ispGOloss <- ispSGO & !dadm$winner & ispGOacc
            like[tGO] <- like[tGO] + apply(matrix(log(1-model$pfunG(
              rt=dadm$rt[ispGOloss],pars=pars[ispGOloss,,drop=FALSE])),
              nrow=n_accG-1),2,sum)
          }
          # trigger stop, add in stop survivor
          ts <- like[tGO] + log(1-model$pfunS(
            rt=dadm$rt[ispGOwin],pars=pars[ispGOwin,,drop=FALSE]))
          # ST loosers
          if (n_accST == 0) stl <- 0 else {
            ispSTloss <- ispSGO & !ispGOacc
            stl <- apply(matrix(log(1-model$pfunG(
              rt=dadm$rt[ispSTloss]-pars[ispSTloss,"SSD"], # correct for SSD
              pars=pars[ispSTloss,,drop=FALSE])),
              nrow=n_accST),2,sum)
          }

          # Transform back to densities to include failures
          like[tGO] <- (1-gf[tGO])*(tf[tGO]*exp(like[tGO]) +
                                      (1-tf[tGO])*exp(ts+stl))
        }

        # ST WINS (never tf)
        if (any(ispSST)) { # Go triggered 1a) with (1-pStop), all race,
          #              1b) with pStop, only ST races
          # Go failure    2) Only ST races
          # ST winner rows
          ispSSTwin <-  dadm$winner &  ispSST
          # Stop probability on ST win trials, only needs go accumulators
          pStop <- model$sfun(pars[ispSST & ispGOacc,,drop=FALSE],n_acc=n_accG,
                              upper=dadm$rt[ispSSTwin]) # pStop before observed RT
          # ST win ll
          tST <- ptrials[ispSSTwin]
          like[tST] <- log(model$dfunG(
            rt=dadm$rt[ispSSTwin]-dadm$SSD[ispSSTwin], # correct ST racers for SSD delay
            pars=pars[ispSSTwin,,drop=FALSE]))
          # ST looser survivors
          if (n_accST > 1) {  # Survivor for looser for ST accumulator(s)
            ispSSTloss <-  !dadm$winner &  ispSST & !ispGOacc
            llST <-  log(1-model$pfunG(
              rt=dadm$rt[ispSSTloss]-dadm$SSD[ispSSTloss],
              pars=pars[ispSSTloss,,drop=FALSE]))
            if (n_accST == 2) # Could remove branch, maybe faster as no matrix sum?
              like[tST] <- like[tST] + llST else
                like[tST] <- like[tST] + apply(matrix(llST,nrow=n_accST-1),2,sum)
          }
          # Go looser survivor
          ispSGloss <- ispSST & ispGOacc
          llG <- apply(matrix(log(1-model$pfunG(
            rt=dadm$rt[ispSGloss],pars=pars[ispSGloss,,drop=FALSE])),
            nrow=n_accG),2,sum)
          like[tST] <- (1-tf[tST])*(                 # Never trigger failure
            gf[tST]*exp(like[tST]) +             # Case 2, gf only ST race
              (1-gf[tST])*(pStop*exp(like[tST]) +      # Case 1b, no gf, stop beats go, ST race
                             (1-pStop)*exp(like[tST]+llG))) # Case 1a, no  gf, all race (no stop win)
        }
      }
      allLL[allok][allr] <- log(like)
    }

    allLL[is.na(allLL)|is.nan(allLL)] <- min_ll
    allLL <- pmax(min_ll,allLL)

    llR_old <<- allLL

    sum(allLL[attr(dadm,"expand")])
}

