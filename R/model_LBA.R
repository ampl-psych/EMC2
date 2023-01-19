dLBA <- function (rt, pars, posdrift = TRUE, robust = FALSE)
  # posdrift = truncated positive normal rates
  # robust slower, deals with extreme rate values
{

  dlba_norm <- function (dt,A,b,v,sv,posdrift=TRUE,robust=FALSE)
    # like dlba_norm_core but t0 dealt with outside (removed from dt)
  {

    pnormP <- function (x, mean = 0, sd = 1, lower.tail = TRUE)
      ifelse(abs(x) < 7, pnorm(x, mean = mean, sd = sd, lower.tail = lower.tail),
             ifelse(x < 0, 0, 1))

    dnormP <- function (x, mean = 0, sd = 1)
      ifelse(abs(x) < 7, dnorm(x, mean = mean, sd = sd), 0)

    if (robust) {
      pnorm1 <- pnormP
      dnorm1 <- dnormP
    } else {
      pnorm1 <- pnorm
      dnorm1 <- dnorm
    }

    if (posdrift)
      denom <- pmax(pnorm1(v/sv), 1e-10) else
        denom <- rep(1, length(t))

    A_small <- A < 1e-10
    if (any(A_small)) {
      out <- numeric(length(dt))
      out[A_small] <- pmax(0, ((b[A_small]/dt[A_small]^2) *
                                 dnorm1(b[A_small]/dt[A_small],v[A_small], sd = sv[A_small]))/denom[A_small])
      zs <- dt[!A_small] * sv[!A_small]
      zu <- dt[!A_small] * v[!A_small]
      chiminuszu <- b[!A_small] - zu
      chizu <- chiminuszu/zs
      chizumax <- (chiminuszu - A[!A_small])/zs
      out[!A_small] <- pmax(0, (v[!A_small] * (pnorm1(chizu) -
                                                 pnorm1(chizumax)) + sv[!A_small] * (dnorm1(chizumax) -
                                                                                       dnorm1(chizu)))/(A[!A_small] * denom[!A_small]))
      return(out)
    } else {
      zs <- dt * sv
      zu <- dt * v
      chiminuszu <- b - zu
      chizu <- chiminuszu/zs
      chizumax <- (chiminuszu - A)/zs
      return(pmax(0, (v * (pnorm1(chizu) - pnorm1(chizumax)) +
                        sv * (dnorm1(chizumax) - dnorm1(chizu)))/(A * denom)))
    }
  }


  dt <- rt - pars[,"t0"]
  tpos <- (dt>0) & (pars[,"b"] >= pars[,"A"])
  out <- numeric(length(dt))
  out[tpos] <- dlba_norm(dt = dt[tpos], A = pars[tpos,"A"], b = pars[tpos,"b"],
                         v = pars[tpos,"v"], sv = pars[tpos,"sv"],
                         posdrift = posdrift, robust = robust)
  out
}


pLBA <- function (rt, pars, posdrift = TRUE, robust = FALSE)
  # posdrift = truncated positive normal rates
  # robust slower, deals with extreme rate values
{

  plba_norm <- function (dt,A,b,v,sv,posdrift=TRUE,robust=FALSE)
    # like plba_norm_core but t0 dealt with outside (removed from dt)
  {

    pnormP <- function (x, mean = 0, sd = 1, lower.tail = TRUE)
      ifelse(abs(x) < 7, pnorm(x, mean = mean, sd = sd, lower.tail=lower.tail),
             ifelse(x < 0, 0, 1))

    dnormP <- function (x, mean = 0, sd = 1)
      ifelse(abs(x) < 7, dnorm(x, mean = mean, sd = sd), 0)

    if (robust) {
      pnorm1 <- pnormP
      dnorm1 <- dnormP
    } else {
      pnorm1 <- pnorm
      dnorm1 <- dnorm
    }
    if (posdrift)
      denom <- pmax(pnorm1(v/sv), 1e-10) else
        denom <- 1
    A_small <- A < 1e-10
    if (any(A_small)) {
      out <- numeric(length(dt))
      out[A_small] <- pmin(1, pmax(0, (pnorm1(b[A_small]/dt[A_small],
                                              mean = v[A_small], sd = sv[A_small],
                                              lower.tail = FALSE))/denom[A_small]))
      zs <- dt[!A_small] * sv[!A_small]
      zu <- dt[!A_small] * v[!A_small]
      chiminuszu <- b[!A_small] - zu
      xx <- chiminuszu - A[!A_small]
      chizu <- chiminuszu/zs
      chizumax <- xx/zs
      tmp1 <- zs * (dnorm1(chizumax) - dnorm1(chizu))
      tmp2 <- xx * pnorm1(chizumax) - chiminuszu * pnorm1(chizu)
      out[!A_small] <- pmin(pmax(0, (1 + (tmp1 + tmp2)/A[!A_small])/denom[!A_small]),1)
      return(out)
    } else {
      zs <- dt * sv
      zu <- dt * v
      chiminuszu <- b - zu
      xx <- chiminuszu - A
      chizu <- chiminuszu/zs
      chizumax <- xx/zs
      tmp1 <- zs * (dnorm1(chizumax) - dnorm1(chizu))
      tmp2 <- xx * pnorm1(chizumax) - chiminuszu * pnorm1(chizu)
      return(pmin(pmax(0, (1 + (tmp1 + tmp2)/A)/denom), 1))
    }
  }


  dt <- rt - pars[,"t0"]
  tpos <- (dt>0) & (pars[,"b"] >= pars[,"A"])
  out <- numeric(length(dt))
  out[tpos] <- plba_norm(dt = dt[tpos], A = pars[tpos,"A"], b = pars[tpos,"b"],
                         v = pars[tpos,"v"], sv = pars[tpos,"sv"],
                         posdrift = posdrift, robust = robust)
  out
}

#### random

rLBA <- function(lR,pars,p_types=c("v","sv","b","A","t0"),posdrift = TRUE)
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars is a matrix of corresponding parameter values named as in p_types
  # pars must be sorted so accumulators and parameter for each trial are in
  # contiguous rows.
{
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ",paste(p_types,collapse = " "))
  dt <- matrix((pars[,"b"]-pars[,"A"]*runif(dim(pars)[1]))/
                 msm::rtnorm(dim(pars)[1],pars[,"v"],pars[,"sv"],ifelse(posdrift,0,-Inf)),
               nrow=length(levels(lR)))
  R <- apply(dt,2,which.min)
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  # Any t0 difference with lR due to response production time (no effect on race)
  rt <- matrix(pars[,"t0"],nrow=length(levels(lR)))[pick] + dt[pick]
  R <- factor(levels(lR)[R],levels=levels(lR))
  bad <- !is.finite(rt)
  R[bad] <- NA
  rt[bad] <- NA
  cbind.data.frame(R=R,rt=rt)
}

# lba_B parameterization
lbaB <- function(){
  list(
    type="RACE",
    # p_vector transform, sets sv as a scaling parameter
    p_types=c("v","sv","B","A","t0"),
    transform = function(p) p,
    # Transform to natural scale
    Ntransform=function(x) {
      x[,dimnames(x)[[2]] != "v"] <- exp(x[,dimnames(x)[[2]] != "v"])
      x
    },
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars <- cbind(pars,b=pars[,"B"] + pars[,"A"])
      attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0)
      pars
    },
    # Random function for racing accumulator
    rfun=function(lR,pars) rLBA(lR,pars,posdrift=TRUE),
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dLBA(rt,pars,posdrift = TRUE, robust = FALSE),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pLBA(rt,pars,posdrift = TRUE, robust = FALSE),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)){
      log_likelihood_race(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
    }
  )
}

