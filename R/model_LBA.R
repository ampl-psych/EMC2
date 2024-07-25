#### Standard LBA ----
# # Moved to C+_ in model_LBA.cpp
#
# pnormP <- function (x, mean = 0, sd = 1, lower.tail = TRUE)
#       ifelse(abs(x) < 7, pnorm(x, mean = mean, sd = sd, lower.tail = lower.tail),
#              ifelse(x < 0, 0, 1))
#
# dnormP <- function (x, mean = 0, sd = 1)
#       ifelse(abs(x) < 7, dnorm(x, mean = mean, sd = sd), 0)
#
#
# dlba_norm <- function (dt,A,b,v,sv,posdrift=TRUE,robust=FALSE)
#     # like dlba_norm_core but t0 dealt with outside (removed from dt)
# {
#
#
#     if (robust) {
#       pnorm1 <- pnormP
#       dnorm1 <- dnormP
#     } else {
#       pnorm1 <- pnorm
#       dnorm1 <- dnorm
#     }
#
#     if (posdrift)
#       denom <- pmax(pnorm1(v/sv), 1e-10) else
#         denom <- rep(1, length(t))
#
#     A_small <- A < 1e-10
#     if (any(A_small)) {
#       out <- numeric(length(dt))
#       out[A_small] <- pmax(0, ((b[A_small]/dt[A_small]^2) *
#                                  dnorm1(b[A_small]/dt[A_small],v[A_small], sd = sv[A_small]))/denom[A_small])
#       zs <- dt[!A_small] * sv[!A_small]
#       zu <- dt[!A_small] * v[!A_small]
#       chiminuszu <- b[!A_small] - zu
#       chizu <- chiminuszu/zs
#       chizumax <- (chiminuszu - A[!A_small])/zs
#       out[!A_small] <- pmax(0, (v[!A_small] * (pnorm1(chizu) -
#                                                  pnorm1(chizumax)) + sv[!A_small] * (dnorm1(chizumax) -
#                                                                                        dnorm1(chizu)))/(A[!A_small] * denom[!A_small]))
#       return(out)
#     } else {
#       zs <- dt * sv
#       zu <- dt * v
#       chiminuszu <- b - zu
#       chizu <- chiminuszu/zs
#       chizumax <- (chiminuszu - A)/zs
#       return(pmax(0, (v * (pnorm1(chizu) - pnorm1(chizumax)) +
#                         sv * (dnorm1(chizumax) - dnorm1(chizu)))/(A * denom)))
#     }
# }
#
#
# plba_norm <- function (dt,A,b,v,sv,posdrift=TRUE,robust=FALSE)
#     # like plba_norm_core but t0 dealt with outside (removed from dt)
# {
#
#     if (robust) {
#       pnorm1 <- pnormP
#       dnorm1 <- dnormP
#     } else {
#       pnorm1 <- pnorm
#       dnorm1 <- dnorm
#     }
#     if (posdrift)
#       denom <- pmax(pnorm1(v/sv), 1e-10) else
#         denom <- 1
#     A_small <- A < 1e-10
#     if (any(A_small)) {
#       out <- numeric(length(dt))
#       out[A_small] <- pmin(1, pmax(0, (pnorm1(b[A_small]/dt[A_small],
#                                               mean = v[A_small], sd = sv[A_small],
#                                               lower.tail = FALSE))/denom[A_small]))
#       zs <- dt[!A_small] * sv[!A_small]
#       zu <- dt[!A_small] * v[!A_small]
#       chiminuszu <- b[!A_small] - zu
#       xx <- chiminuszu - A[!A_small]
#       chizu <- chiminuszu/zs
#       chizumax <- xx/zs
#       tmp1 <- zs * (dnorm1(chizumax) - dnorm1(chizu))
#       tmp2 <- xx * pnorm1(chizumax) - chiminuszu * pnorm1(chizu)
#       out[!A_small] <- pmin(pmax(0, (1 + (tmp1 + tmp2)/A[!A_small])/denom[!A_small]),1)
#       return(out)
#     } else {
#       zs <- dt * sv
#       zu <- dt * v
#       chiminuszu <- b - zu
#       xx <- chiminuszu - A
#       chizu <- chiminuszu/zs
#       chizumax <- xx/zs
#       tmp1 <- zs * (dnorm1(chizumax) - dnorm1(chizu))
#       tmp2 <- xx * pnorm1(chizumax) - chiminuszu * pnorm1(chizu)
#       return(pmin(pmax(0, (1 + (tmp1 + tmp2)/A)/denom), 1))
#     }
# }
#
#

dLBA <- function (rt, pars, posdrift = TRUE, robust = FALSE)
  # posdrift = truncated positive normal rates
  # robust slower, deals with extreme rate values
{
  dt <- rt - pars[,"t0"]
  ok <- (dt>0) & (pars[,"b"] >= pars[,"A"])
  ok[is.na(ok) | !is.finite(dt)] <- FALSE
  out <- numeric(length(dt))
  out[ok] <- dlba(t = dt[ok], A = pars[ok,"A"], b = pars[ok,"b"],
                         v = pars[ok,"v"], sv = pars[ok,"sv"],
                         posdrift = posdrift, robust = robust)
  out
}

pLBA <- function (rt, pars, posdrift = TRUE, robust = FALSE)
  # posdrift = truncated positive normal rates
  # robust slower, deals with extreme rate values
{
  dt <- rt - pars[,"t0"]
  ok <- (dt>0) & (pars[,"b"] >= pars[,"A"])
  ok[is.na(ok) | !is.finite(dt)] <- FALSE
  out <- numeric(length(dt))
  out[ok] <- plba(t = dt[ok], A = pars[ok,"A"], b = pars[ok,"b"],
                         v = pars[ok,"v"], sv = pars[ok,"sv"],
                         posdrift = posdrift, robust = robust)
  out
}


rLBA <- function(lR,pars,p_types=c("v","sv","b","A","t0"),posdrift = TRUE,
                 ok=rep(TRUE,length(lR)))
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars is a matrix of corresponding parameter values named as in p_types
  # pars must be sorted so accumulators and parameter for each trial are in
  # contiguous rows.
{
  bad <- rep(NA, length(lR)/length(levels(lR)))
  out <- data.frame(R = bad, rt = bad)
  nr <- length(levels(lR))
  dt <- matrix(Inf,nrow=nr,ncol=nrow(pars)/nr)
  t0 <- pars[,"t0"]
  pars <- pars[ok,]
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ",paste(p_types,collapse = " "))
  dt[ok] <- (pars[,"b"]-pars[,"A"]*runif(dim(pars)[1]))/
    msm::rtnorm(dim(pars)[1],pars[,"v"],pars[,"sv"],ifelse(posdrift,0,-Inf))
  dt[dt<0] <- Inf
  bad <- apply(dt,2,function(x){all(is.infinite(x))})
  R <- apply(dt,2,which.min)
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  # Any t0 difference with lR due to response production time (no effect on race)
  rt <- matrix(t0,nrow=nr)[pick] + dt[pick]
  R <- factor(levels(lR)[R],levels=levels(lR))
  R[bad] <- NA
  rt[bad] <- Inf
  ok <- matrix(ok,nrow=length(levels(lR)))[1,]
  out$R[ok] <- levels(lR)[R][ok]
  out$R <- factor(out$R,levels=levels(lR))
  out$rt[ok] <- rt[ok]
  out
}
#### Model functions ----

#' The Linear Ballistic Accumulator model
#'
#' Model file to estimate the Linear Ballistic Accumulator (LBA) in EMC2.
#'
#' Model files are almost exclusively used in `design()`.
#'
#' @details
#'
#' Default values are used for all parameters that are not explicitly listed in the `formula`
#' argument of `design()`.They can also be accessed with `LBA()$p_types`.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default**   | **Mapping**                    | **Interpretation**                                            |
#' |-----------|-----------|---------------|-----------|----------------------------|-----------------------------------------------------------|
#' | *v*       | -         | \[-Inf, Inf\] | 1         |                            | Mean evidence-accumulation rate                                              |
#' | *A*       | log       | \[0, Inf\]    | log(0)    |                            | Between-trial variation (range) in start point                     |
#' | *B*       | log       | \[0, Inf\]    | log(1)    | *b* = *B*+*A*              | Distance from *A* to *b* (response threshold)                                       |
#' | *t0*      | log       | \[0, Inf\]    | log(0)    |                            | Non-decision time                                         |
#' | *sv*      | log       | \[0, Inf\]    | log(1)    |                            | Between-trial variation in evidence-accumulation rate                      |
#'
#'
#' All parameters are estimated on the log scale, except for the drift rate which is estimated on the real line.
#'
#' Conventionally, `sv` is fixed to 1 to satisfy scaling constraints.
#'
#' The *b* = *B* + *A* parameterization ensures that the response threshold is always higher than the between trial variation in start point of the drift rate.
#'
#' Because the LBA is a race model, it has one accumulator per response option.
#' EMC2 automatically constructs a factor representing the accumulators `lR` (i.e., the
#' latent response) with level names taken from the `R` column in the data.
#'
#' The `lR` factor is mainly used to allow for response bias, analogous to `Z` in the
#' DDM. For example, in the LBA, response thresholds are determined by the *B*
#' parameters, so `B~lR` allows for different thresholds for the accumulator
#' corresponding to left and right stimuli (e.g., a bias to respond left occurs
#' if the left threshold is less than the right threshold).
#' For race models, the `design()` argument `matchfun` can be provided, a
#' function that takes the `lR` factor (defined in the augmented data (d)
#' in the following function) and returns a logical defining the correct response.
#' In the example below, the match is simply such that the `S` factor equals the
#' latent response factor: `matchfun=function(d)d$S==d$lR`. Then `matchfun` is
#' used to automatically create a latent match (`lM`) factor with
#' levels `FALSE` (i.e., the stimulus does not match the accumulator) and `TRUE`
#' (i.e., the stimulus does match the accumulator). This is added internally
#' and can also be used in model formula, typically for parameters related to
#' the rate of accumulation.
#'
#' Brown, S. D., & Heathcote, A. (2008). The simplest complete model of choice response time: Linear ballistic accumulation.
#' *Cognitive Psychology, 57*(3), 153-178. https://doi.org/10.1016/j.cogpsych.2007.12.002
#'
#' @return A model list with all the necessary functions for EMC2 to sample
#' @examples
#' # When working with lM it is useful to design  an "average and difference"
#' # contrast matrix, which for binary responses has a simple canonical from:

#' ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
#' # We also define a match function for lM
#' matchfun=function(d)d$S==d$lR
#' # We now construct our design, with v ~ lM and the contrast for lM the ADmat.
#' design_LBABE <- design(data = forstmann,model=LBA,matchfun=matchfun,
#'                        formula=list(v~lM,sv~lM,B~E+lR,A~1,t0~1),
#'                        contrasts=list(v=list(lM=ADmat)),constants=c(sv=log(1)))
#' # For all parameters that are not defined in the formula, default values are assumed
#' # (see Table above).
#' @export
#'

LBA <- function(){
  list(
    type="RACE",
    c_name = "LBA",
    # p_vector transform, sets sv as a scaling parameter
    p_types=c("v" = 1,"sv" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0)),
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
    rfun=function(lR=NULL,pars) {
      ok <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0)
      if (is.null(lR)) ok else rLBA(lR,pars,posdrift=TRUE,ok=ok)
    },
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



