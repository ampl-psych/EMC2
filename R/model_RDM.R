

# n-choice uniformly varying start point (0-A) Wald
#    race, with t0, v, A, b (boundary) parameterixaiton

# pigt, digt, rwaldt Copyright (C) 2013  Trisha Van Zandt distributed with:
# Logan, Van Zandt, Verbruggen, and Wagenmakers (2014).  On the ability to
# inhibit thought and action: General and special theories of an act of control.
# Psychological Review. Comments and changes added by Andrew Heathcote. Trish's
# code is for k = threshold, a = half width of uniform threshold variability,
# l = rate of accumulation. Note that Wald mean = k/l and shape = k^2.
#
# I use a different parameterization in terms of v=l (rate),
# uniform start point variability from 0-A (A>=0), threshold b (>0) and hence
# B=b-A (>=0) as a threshold gap. Hence k = b-A/2 = B + A/2 and a=A
#
# Added ability to set s = diffusive standard deviation, by scaling A, B and v
# i.e., given s, same result for A/s, B/s, v/s with s=1

#### distribution functions

# # Moved to C++ model_RDM.cpp
#
# digt0 <- function(t,k=1,l=1) {
#       # pdf of inverse gaussian at t with no k variability
#       # much faster than statmod's dinvgauss funciton
#
#       lambda <- k^2
#       l0 <- l==0
#       e <- numeric(length(t))
#       if ( any(!l0) ) {
#         mu <- k[!l0]/l[!l0]
#         e[!l0] <- -(lambda[!l0]/(2*t[!l0])) * (t[!l0]^2/mu^2 - 2*t[!l0]/mu  + 1)
#       }
#       if ( any(l0) )  e[l0] <- -.5*lambda[l0]/t[l0]
#       x <- exp(e + .5*log(lambda) - .5*log(2*t^3*pi))
#       x[t<=0] <- 0
#       x
# }
#
#
# digt <- function(t,k=1,l=1,a=.1,tiny=1e-10) {
#     # pdf of inverse gaussian at t with k +/- a/2 uniform variability
#     # returns digt0 if a<1e-10
#
#
#     options(warn=-1)
#     if(length(k)!=length(t)) k <- rep(k,length.out=length(t))
#     if(length(l)!=length(t)) l <- rep(l,length.out=length(t))
#     if(length(a)!=length(t)) a <- rep(a,length.out=length(t))
#
#     tpos <- t<=0
#
#     atiny <- a<=tiny & !tpos
#     a[atiny] <- 0
#
#     ltiny <- (l<=tiny) & !atiny & !tpos
#     notltiny <- (l>tiny) & !atiny & !tpos
#     l[l<=tiny] <- 0
#
#     x <- numeric(length(t))
#
#     # No threshold variability
#     if ( any(atiny) )
#       x[atiny] <- digt0(t=t[atiny],k=k[atiny],l=l[atiny])
#
#     # Threshold variability
#     if ( any(!atiny) ) {
#
#       if ( any(notltiny) ) { # rate non-zero
#
#         sqr.t <- sqrt(t[notltiny])
#
#         term.1a <- -(a[notltiny]-k[notltiny]+t[notltiny]*l[notltiny])^2/(2*t[notltiny])
#         term.1b <- -(a[notltiny]+k[notltiny]-t[notltiny]*l[notltiny])^2/(2*t[notltiny])
#         term.1 <- (exp(term.1a) - exp(term.1b))/sqrt(2*pi*t[notltiny])
#
#         term.2a <- log(.5)+log(l[notltiny])
#         term.2b <- 2*pnorm((-k[notltiny]+a[notltiny])/sqr.t+sqr.t*l[notltiny])-1
#         term.2c <- 2*pnorm((k[notltiny]+a[notltiny])/sqr.t-sqr.t*l[notltiny])-1
#         term.2d <- term.2b+term.2c
#         term.2 <- exp(term.2a)*term.2d
#
#         term.3 <- term.1+term.2
#         term.4 <- log(term.3)-log(2)-log(a[notltiny])
#         x[notltiny] <- exp(term.4)
#       }
#
#       if ( any(ltiny) ) {  # rate zero
#         log.t <- log(t[ltiny])
#         term.1 <- -.5*(log(2)+log(pi)+log.t)
#         term.2 <- (k[ltiny]-a[ltiny])^2/(2*t[ltiny])
#         term.3 <- (k[ltiny]+a[ltiny])^2/(2*t[ltiny])
#         term.4 <- (exp(-term.2)-exp(-term.3))
#         term.5 <- term.1+log(term.4) - log(2) - log(a[ltiny])
#         x[ltiny] <- exp(term.5)
#       }
#
#     }
#
#     x[x<0 | is.nan(x) ] <- 0
#     x
# }
#
#
# pigt0 <- function(t,k=1,l=1) {
#       # cdf of inverse gaussian at t with no k variability
#       # much faster than statmod's pinvgauss funciton
#
#       mu <- k/l
#       lambda <- k^2
#
#       e <- exp(log(2*lambda) - log(mu))
#       add <- sqrt(lambda/t) * (1 + t/mu)
#       sub <- sqrt(lambda/t) * (1 - t/mu)
#
#       p.1 <- 1 - pnorm(add)
#       p.2 <- 1 - pnorm(sub)
#       x <- exp(e + log(p.1)) + p.2
#
#       x[t<0] <- 0
#       x
# }
#
#
# pigt <- function(t,k=1,l=1,a=.1,tiny=1e-10) {
#     # cdf of inverse gaussian at t with k +/- a/2 uniform variability
#     # returns pigt0 if a<=0
#
#     options(warn=-1)
#     if(length(k)!=length(t)) k <- rep(k,length.out=length(t))
#     if(length(l)!=length(t)) l <- rep(l,length.out=length(t))
#     if(length(a)!=length(t)) a <- rep(a,length.out=length(t))
#
#     tpos <- t<=0
#
#     atiny <- a<=tiny & !tpos
#     a[atiny] <- 0
#
#     ltiny <- (l<=tiny) & !atiny & !tpos
#     notltiny <- (l>tiny) & !atiny & !tpos
#     l[l<=tiny] <- 0
#
#     x <- numeric(length(t))
#
#     # No threshold variability
#     if ( any(atiny) )
#       x[atiny] <- pigt0(t[atiny],k[atiny],l[atiny])
#
#     # Threshold variability
#     if ( any(!atiny) ) {
#
#       if ( any(notltiny) ) { # rate non-zero
#
#         log.t <- log(t[notltiny])
#         sqr.t <- sqrt(t[notltiny])
#
#         term.1a <- .5*log.t-.5*log(2*pi)
#         term.1b <- exp(-((k[notltiny]-a[notltiny]-t[notltiny]*l[notltiny])^2/t[notltiny])/2)
#         term.1c <- exp(-((k[notltiny]+a[notltiny]-t[notltiny]*l[notltiny])^2/t[notltiny])/2)
#         term.1 <- exp(term.1a)*(term.1b-term.1c)
#
#         term.2a <- exp(2*l[notltiny]*(k[notltiny]-a[notltiny]) +
#                          log(pnorm(-(k[notltiny]-a[notltiny]+t[notltiny]*l[notltiny])/sqr.t)))
#         term.2b <- exp(2*l[notltiny]*(k[notltiny]+a[notltiny]) +
#                          log(pnorm(-(k[notltiny]+a[notltiny]+t[notltiny]*l[notltiny])/sqr.t)))
#         term.2 <- a[notltiny] + (term.2b-term.2a)/(2*l[notltiny])
#
#         term.4a <- 2*pnorm((k[notltiny]+a[notltiny])/sqr.t-sqr.t*l[notltiny])-1
#         term.4b <- 2*pnorm((k[notltiny]-a[notltiny])/sqr.t-sqr.t*l[notltiny])-1
#         term.4c <- .5*(t[notltiny]*l[notltiny] - a[notltiny] - k[notltiny] + .5/l[notltiny])
#         term.4d <- .5*(k[notltiny] - a[notltiny] - t[notltiny]*l[notltiny] - .5/l[notltiny])
#         term.4 <- term.4c*term.4a + term.4d*term.4b
#
#         x[notltiny] <- (term.4 + term.2 + term.1)/(2*a[notltiny])
#       }
#
#       if ( any(ltiny) ) {  # rate zero
#         sqr.t <- sqrt(t[ltiny])
#         log.t <- log(t[ltiny])
#         term.5a <- 2*pnorm((k[ltiny]+a[ltiny])/sqr.t)-1
#         term.5b <- 2*pnorm(-(k[ltiny]-a[ltiny])/sqr.t)-1
#         term.5 <- (-(k[ltiny]+a[ltiny])*term.5a - (k[ltiny]-a[ltiny])*term.5b)/(2*a[ltiny])
#
#         term.6a <- -.5*(k[ltiny]+a[ltiny])^2/t[ltiny] - .5*log(2) -.5*log(pi) + .5*log.t - log(a[ltiny])
#         term.6b <- -.5*(k[ltiny]-a[ltiny])^2/t[ltiny] - .5*log(2) -.5*log(pi) + .5*log.t - log(a[ltiny])
#         term.6 <- 1 + exp(term.6b) - exp(term.6a)
#
#         x[ltiny] <- term.5 + term.6
#       }
#
#     }
#
#     x[x<0 | is.nan(x) ] <- 0
#     x
# }
#
#

dRDM <- function(rt,pars)
  # density for single accumulator
{
  out <- numeric(length(rt))
  ok <- rt > pars[,"t0"] & !pars[,"v"] < 0  # code handles rate zero case
  ok[is.na(ok)] <- FALSE
  if (any(dimnames(pars)[[2]]=="s")) # rescale
    pars[ok,c("A","B","v")] <- pars[ok,c("A","B","v")]/pars[ok,"s"]
  out[ok] <- dWald(rt[ok],v=pars[ok,"v"],B=pars[ok,"B"],A=pars[ok,"A"],t0=pars[ok,"t0"])
  out
}


pRDM <- function(rt,pars)
  # cumulative density for single accumulator
{
  out <- numeric(length(rt))
  ok <- rt > pars[,"t0"] & !pars[,"v"] < 0  # code handles rate zero case
  ok[is.na(ok)] <- FALSE
  if (any(dimnames(pars)[[2]]=="s")) # rescale
    pars[ok,c("A","B","v")] <- pars[ok,c("A","B","v")]/pars[ok,"s"]
  out[ok] <- pWald(rt[ok],v=pars[ok,"v"],B=pars[ok,"B"],A=pars[ok,"A"],t0=pars[ok,"t0"])
  out
}

#### random

rWald <- function(n,B,v,A)
  # random function for single accumulator
{

  rwaldt <- function(n,k,l,tiny=1e-6) {
    # random sample of n from a Wald (or Inverse Gaussian)
    # k = criterion, l = rate, assumes sigma=1 Browninan motion
    # about same speed as statmod rinvgauss

    rlevy <- function(n=1, m=0, c=1) {
      if (any(c<0)) stop("c must be positive")
      c/qnorm(1-runif(n)/2)^2+m
    }

    flag <- l>tiny
    x <- rep(NA,times=n)

    x[!flag] <- rlevy(sum(!flag),0,k[!flag]^2)
    mu <- k/l
    lambda <- k^2

    y <- rnorm(sum(flag))^2
    mu.0 <- mu[flag]
    lambda.0 <- lambda[flag]

    x.0 <- mu.0 + mu.0^2*y/(2*lambda.0) -
      sqrt(4*mu.0*lambda.0*y + mu.0^2*y^2)*mu.0/(2*lambda.0)

    z <- runif(length(x.0))
    test <- mu.0/(mu.0+x.0)
    x.0[z>test] <- mu.0[z>test]^2/x.0[z>test]
    x[flag] <- x.0
    x[x<0] <- max(x)
    x
  }

  # Kluge to return Inf for negative rates
  out <- numeric(n)
  ok <- !v<0
  nok <- sum(ok)
  bs <- B[ok]+runif(nok,0,A[ok])
  out[ok] <- rwaldt(nok,k=bs,l=v[ok])
  out[!ok] <- Inf
  out
}

rRDM <- function(lR,pars,p_types=c("v","B","A","t0"),ok=rep(TRUE,dim(pars)[1]))
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars is a matrix of corresponding parameter values named as in p_types
  # pars must be sorted so accumulators and parameter for each trial are in
  # contiguous rows. "s" parameter will be used but can be ommitted
  #
  # test
  # pars=cbind(B=c(1,2),v=c(1,1),A=c(0,0),t0=c(.2,.2)); lR=factor(c(1,2))
{
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ",paste(p_types,collapse = " "))
  if (any(dimnames(pars)[[2]]=="s")) # rescale
    pars[,c("A","B","v")] <- pars[,c("A","B","v")]/pars[,"s"]
  pars[,"B"][pars[,"B"]<0] <- 0 # Protection for negatives
  pars[,"A"][pars[,"A"]<0] <- 0
  bad <- rep(NA, length(lR)/length(levels(lR)))
  out <- data.frame(R = bad, rt = bad)
  nr <- length(levels(lR))
  dt <- matrix(Inf,nrow=nr,ncol=nrow(pars)/nr)
  t0 <- pars[,"t0"]
  pars <- pars[ok,]
  dt[ok] <- rWald(sum(ok),B=pars[,"B"],v=pars[,"v"],A=pars[,"A"])
  R <- apply(dt,2,which.min)
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  # Any t0 difference with lR due to response production time (no effect on race)
  rt <- matrix(t0,nrow=nr)[pick] + dt[pick]
  out$R <- levels(lR)[R]
  out$R <- factor(out$R,levels=levels(lR))
  out$rt <- rt
  out
}

#' The Racing Diffusion Model
#'
#' Model file to estimate the Racing Diffusion Model (RDM), also known as the Racing Wald Model.
#'
#' Model files are almost exclusively used in `design()`.
#'
#' @details
#'
#' Default values are used for all parameters that are not explicitly listed in the `formula`
#' argument of `design()`.They can also be accessed with `RDM()$p_types`.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default**   | **Mapping**          | **Interpretation**                                                |
#' |-----------|-----------|---------------|-----------|------------------|---------------------------------------------------------------|
#' | *v*       | log       | \[0, Inf\]      | log(1)    |                  | Evidence-accumulation rate (drift rate)                        |
#' | *A*       | log       | \[0, Inf\]      | log(0)    |                  | Between-trial variation (range) in start point                 |
#' | *B*       | log       | \[0, Inf\]      | log(1)    | *b* = *B* + *A*      | Distance from *A* to *b* (response threshold)                  |
#' | *t0*      | log       | \[0, Inf\]      | log(0)    |                  | Non-decision time                                             |
#' | *sv*      | log       | \[0, Inf\]      | log(1)    |                  | Within-trial standard deviation of drift rate                 |
#'
#'
#' All parameters are estimated on the log scale.
#'
#' The parameterization *b* = *B* + *A* ensures that the response threshold is
#' always higher than the between trial variation in start point.
#'
#' Conventionally, `s` is fixed to 1 to satisfy scaling constraints.
#'
#' Because the RDM is a race model, it has one accumulator per response option.
#' EMC2 automatically constructs a factor representing the accumulators `lR` (i.e., the
#' latent response) with level names taken from the `R` column in the data.
#'
#' The `lR` factor is mainly used to allow for response bias, analogous to *Z* in the
#' DDM. For example, in the RDM, response thresholds are determined by the *B*
#' parameters, so `B~lR` allows for different thresholds for the accumulator
#' corresponding to "left" and "right" stimuli, for example, (e.g., a bias to respond left occurs
#' if the left threshold is less than the right threshold).
#'
#' For race models in general, the argument `matchfun` can be provided in `design()`.
#' One needs to supply a function that takes the `lR` factor (defined in the augmented data (d)
#' in the following function) and returns a logical defining the correct
#' response. In the example below, this is simply whether the `S` factor equals the
#' latent response factor: `matchfun=function(d)d$S==d$lR`. Using `matchfun` a latent match factor (`lM`) with
#' levels `FALSE` (i.e., the stimulus does not match the accumulator) and `TRUE`
#' (i.e., the stimulus does match the accumulator). This is added internally
#' and can also be used in model formula, typically for parameters related to
#' the rate of accumulation.
#'
#' Tillman, G., Van Zandt, T., & Logan, G. D. (2020). Sequential sampling models
#' without random between-trial variability: The racing diffusion model of speeded
#' decision making. *Psychonomic Bulletin & Review, 27*(5), 911-936.
#' https://doi.org/10.3758/s13423-020-01719-6
#'
#' @return A list defining the cognitive model
#' @examples

#' # When working with lM it is useful to design  an "average and difference"
#' # contrast matrix, which for binary responses has a simple canonical from:

#' ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
#' # We also define a match function for lM
#' matchfun=function(d)d$S==d$lR
#' # We now construct our design, with v ~ lM and the contrast for lM the ADmat.
#' design_RDMBE <- design(data = forstmann,model=RDM,matchfun=matchfun,
#'                        formula=list(v~lM,s~lM,B~E+lR,A~1,t0~1),
#'                        contrasts=list(v=list(lM=ADmat)),constants=c(s=log(1)))
#' # For all parameters that are not defined in the formula, default values are assumed
#' # (see Table above).
#' @export

RDM <- function(){
  list(
    type="RACE",
    c_name = "RDM",
    p_types=c("v" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0),"s" = log(1)),
    transform=list(func=c(v = "exp", B = "exp", A = "exp",t0 = "exp", s = "exp")
                   , lower=c(t0=.05)),
    bound=list(minmax=cbind(v=c(1e-3,Inf), B=c(0,Inf), A=c(1e-4,Inf),t0=c(0,Inf), s=c(0,Inf)),
               exception=c(A=0, v=0)),
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars <- cbind(pars,b=pars[,"B"] + pars[,"A"])
      pars
    },
    # Random function for racing accumulators
    rfun=function(lR=NULL,pars)  rRDM(lR,pars,ok=attr(pars, "ok")),
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dRDM(rt,pars),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pRDM(rt,pars),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
      log_likelihood_race(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
  )
}
