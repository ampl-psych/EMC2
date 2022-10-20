rDDM <- function(lR,pars,precision=3,ok=NULL)
  # lR is an empty latent response factor lR with one level for each boundary
  # pars is a matrix of parameter values named as in p_types
  # lower is mapped to first level of lR and upper to second
  # test
  # pars=cbind.data.frame(a=c(1,1),v=c(-1,1),t0=c(.2,.2),z=c(.5,.5),d=c(0,0),
  #                       sz=c(0,0),sv=c(0,0),st0=c(0,0),s=c(1,1))
  # lR <- factor(c("left","right"))

{
  if (is.null(ok)) ok <- rep(TRUE,length(lR))
  bad <- rep(NA,length(lR))
  out <- data.frame(response=bad,rt=bad)
  out[ok,2:1] <- rtdists::rdiffusion(length(lR[ok]), a = pars[ok,"a"], v = pars[ok,"v"], t0 = pars[ok,"t0"],
                            z = pars[ok,"z"], d = pars[ok,"d"], sz = pars[ok,"sz"], sv = pars[ok,"sv"],
                            st0 = pars[ok,"st0"], s = pars[ok,"s"])
  # cbind.data.frame(R=factor(out[,"response"],levels=c("lower","upper"),labels=levels(lR)),rt=out[,"rt"])
  cbind.data.frame(R=factor(out[,"response"],levels=1:2,labels=levels(lR)),rt=out[,"rt"])
}

dDDM <- function(rt,R,pars,precision=3)
  # DDM density for response factor R with rt
  # lower is mapped to first level of R and upper to second
  # test
  # pars=cbind.data.frame(a=c(1,1),v=c(-1,1),t0=c(.2,.2),z=c(.5,.5),d=c(0,0),
  #                       sz=c(0,0),sv=c(0,0),st0=c(0,0),s=c(1,1))
  # R <- factor(c("left","right")); rt=c(1,1)
{
  levels(R) <- c("lower","upper")
  rtdists::ddiffusion(rt,response=as.character(R),
             a = pars[,"a"], v = pars[,"v"], t0 = pars[,"t0"], z = pars[,"z"],d = pars[,"d"],
             sz = pars[,"sz"], sv = pars[,"sv"],st0 = pars[,"st0"], s = pars[,"s"])
}


pDDM <- function(rt,R,pars,precision=3)
  # DDM cdf for response factor R with rt
  # lower is mapped to first level of R and upper to second
{
  levels(R) <- c("lower","upper")
  rtdists::pdiffusion(rt,response=as.character(R),
             a = pars[,"a"], v = pars[,"v"], t0 = pars[,"t0"], z = pars[,"z"],d = pars[,"d"],
             sz = pars[,"sz"], sv = pars[,"sv"],st0 = pars[,"st0"], s = pars[,"s"])
}

ddmTZD <- list(
  type="DDM",
  p_types=c("v","a","sv","t0","st0","s","Z","SZ","DP"),
  # The "TZD" parameterization defined relative to the "rtdists" package is:
  # natural scale
  #   v = rtdists rate v (positive favors upper)
  # log scale
  #   t0 > 0: lower bound of non-decision time
  #   st0 > 0: rtdists width of non-decision time distribution
  #   a > 0: rtdists upper threshold, a
  #   sv > 0: rtdists v standard deviation sv
  #   s > 0: rtdists moment-to-moment standard deviation, s
  # probit scale
  #   0 < Z < 1: rtdists start point z = Z*a
  #   0 < SZ < 1: rtdists start-point variability, sz = 2*SZ*min(c(a*Z,a*(1-Z))
  #   0 < DP < 1: rtdists d = t0(upper)-t0(lower) = (2*DP-1)*t0  #
  #
  # Transform to natural scale
  Ntransform=function(x) {
    islog <- dimnames(x)[[2]] %in% c("a","sv","t0","st0","s")
    isprobit <- dimnames(x)[[2]] %in% c("Z","SZ","DP")
    x[,islog] <- exp(x[,islog])
    x[,isprobit] <- pnorm(x[,isprobit])
    x
  },
  # Trial dependent parameter transform
  Ttransform = function(pars,dadm) {
    pars <- cbind(pars,z=pars[,"a"]*pars[,"Z"],
                  sz = 2*pars[,"SZ"]*pars[,"a"]*apply(cbind(pars[,"Z"],1-pars[,"Z"]),1,min))
    pars <- cbind(pars, d = pars[,"t0"]*(2*pars[,"DP"]-1))
    attr(pars,"ok") <-
      !( abs(pars[,"v"])> 20 | pars[,"a"]> 10 | pars[,"sv"]> 10 | pars[,"SZ"]> .999 | pars[,"st0"]>.2)
    if (pars[1,"sv"] !=0) attr(pars,"ok") <- attr(pars,"ok") & pars[,"sv"] > .001
    if (pars[1,"SZ"] !=0) attr(pars,"ok") <- attr(pars,"ok") & pars[,"SZ"] > .001
    pars
  },
  # p_vector transform, sets s as a scaling parameter
  transform = function(p) p,
  # Random function
  rfun=function(lR,pars) {
    ok <- !( abs(pars[,"v"])> 20 | pars[,"a"]> 10 | pars[,"sv"]> 10 | pars[,"SZ"]> .999 | pars[,"st0"]>.2)
    if (pars[1,"sv"] !=0) attr(pars,"ok") <- attr(pars,"ok") & pars[,"sv"] > .001
    if (pars[1,"SZ"] !=0) attr(pars,"ok") <- attr(pars,"ok") & pars[,"SZ"] > .001
    rDDM(lR,pars,precision=3,ok)
  },
  # Density function (PDF)
  dfun=function(rt,R,pars) dDDM(rt,R,pars,precision=3),
  # Probability function (CDF)
  pfun=function(rt,R,pars) pDDM(rt,R,pars,precision=3),
  log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)){
    log_likelihood_ddm(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
  }
)


ddmTZDt0natural <- list(
  type="DDM",
  p_types=c("v","a","sv","t0","st0","s","Z","SZ","DP"),
  # The "TZD" parameterization defined relative to the "rtdists" package is:
  # natural scale
  #   v = rtdists rate v (positive favors upper)
  # log scale
  #   t0 > 0: lower bound of non-decision time
  #   st0 > 0: rtdists width of non-decision time distribution
  #   a > 0: rtdists upper threshold, a
  #   sv > 0: rtdists v standard deviation sv
  #   s > 0: rtdists moment-to-moment standard deviation, s
  # probit scale
  #   0 < Z < 1: rtdists start point z = Z*a
  #   0 < SZ < 1: rtdists start-point variability, sz = 2*SZ*min(c(a*Z,a*(1-Z))
  #   0 < DP < 1: rtdists d = t0(upper)-t0(lower) = (2*DP-1)*t0  #
  #
  Ntransform=function(x) {
    islog <- dimnames(x)[[2]] %in% c("a","sv","st0","s")
    isprobit <- dimnames(x)[[2]] %in% c("Z","SZ","DP")
    x[,islog] <- exp(x[,islog])
    x[,isprobit] <- pnorm(x[,isprobit])
    x
  },
  # p_vector transform, sets s as a scaling parameter
  transform = function(p) p,
  # Trial dependent parameter transform
  # Trial dependent parameter transform
  Ttransform = function(pars,dadm) {
    pars <- cbind(pars,z=pars[,"a"]*pars[,"Z"],
                  sz = 2*pars[,"SZ"]*pars[,"a"]*apply(cbind(pars[,"Z"],1-pars[,"Z"]),1,min))
    pars <- cbind(pars, d = pars[,"t0"]*(2*pars[,"DP"]-1))
    attr(pars,"ok") <-
      !( abs(pars[,"v"])> 20 | pars[,"a"]> 10 | pars[,"sv"]> 10 | pars[,"SZ"]> .999 | pars[,"st0"]>.2)
    if (pars[1,"sv"] !=0) attr(pars,"ok") <- attr(pars,"ok") & pars[,"sv"] > .001
    if (pars[1,"SZ"] !=0) attr(pars,"ok") <- attr(pars,"ok") & pars[,"SZ"] > .001
    pars
  },
  # Random function
  rfun=function(lR,pars) {
    ok <- !( abs(pars[,"v"])> 20 | pars[,"a"]> 10 | pars[,"sv"]> 10 | pars[,"SZ"]> .999 | pars[,"st0"]>.2)
    if (pars[1,"sv"] !=0) attr(pars,"ok") <- attr(pars,"ok") & pars[,"sv"] > .001
    if (pars[1,"SZ"] !=0) attr(pars,"ok") <- attr(pars,"ok") & pars[,"SZ"] > .001
    rDDM(lR,pars,precision=3,ok)
  },
  # Density function (PDF)
  dfun=function(rt,R,pars) dDDM(rt,R,pars,precision=3),
  # Probability function (CDF)
  pfun=function(rt,R,pars) pDDM(rt,R,pars,precision=3),
  log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)){
    log_likelihood_ddm(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
  }
)
