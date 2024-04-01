
dLNR <- function(rt,pars){
  rt <- rt - pars[,"t0"]
  out <- numeric(length(rt))
  ok <- rt > 0
  ok[is.na(ok) | is.infinite(rt)] <- FALSE
  out[ok] <- stats::dlnorm(rt[ok],meanlog=pars[ok,"m"],sdlog=pars[ok,"s"])
  out
}

pLNR <- function(rt,pars){
  rt <- rt - pars[,"t0"]
  out <- numeric(length(rt))
  ok <- rt > 0
  ok[is.na(ok) | is.infinite(rt)] <- FALSE
  out[ok] <- stats::plnorm(rt[ok],meanlog=pars[ok,"m"],sdlog=pars[ok,"s"])
  out

}

rLNR <- function(lR,pars,p_types=c("m","s","t0")){
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ",paste(p_types,collapse = " "))
  dt <- matrix(stats::rlnorm(dim(pars)[1],meanlog=pars[,"m"],sdlog=pars[,"s"]),
               nrow=length(levels(lR)))
  R <- apply(dt,2,which.min)
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  # Any t0 difference with lR due to response production time (no effect on race)
  rt <- matrix(pars[,"t0"],nrow=length(levels(lR)))[pick] + dt[pick]
  R <- factor(levels(lR)[R],levels=levels(lR))
  cbind.data.frame(R=R,rt=rt)
}


#' The log-normal race (LNR) model
#'
#' The LNR proposes a race between accumulators associated with each choice alternative.
#' These accumulators have bounds, and the first accumulator to reach its bound determines the response time and response choice.
#' The times at which accumulator reaches its bound is assumed to be lognormally distributed.
#' For details see `Rouder et al., 2015`
#'
#' The parameters of the LNR are the scale `m`,  shape `s` and non-decision time `t0`.
#' Unlike the other race models, no parameters need to be constrained for scaling purposes.
#'
#' Parameters `s` and `t0` are estimated on the log-scale.
#'
#' @return A model list with all the necessary functions to sample
#' @examples
#' # As the LNR is a race model it has one accumulator representing each possible
#' # response. EMC2 uses the Rlevels specification given to make_design to
#' # automatically construct a factor representing the accumulators "lR" (the
#' # latent response) with level names taken from Rlevels.

#' # matchfun is used to automatically create a latent match (lM) factor with
#' # levels "FALSE" (i.e., the stimulus does not match the accumulator) and "TRUE"
#' # (i.e., the stimulus does match the accumulator). This is added internally
#' # and can also be used in model formula, typically for parameters related to
#' # the rate of accumulation.

#' # When working with lM it is useful to design  an "average and difference"
#' # contrast matrix, which for binary responses has a simple canonical from:

#' ADmat <- matrix(c(-1/2,1/2),ncol=1,dimnames=list(NULL,"d"))
#' # We now construct our design, with v ~ lM and the contrast for lM the ADmat.
#' design_LNRmE <- make_design(data = forstmann,model=LNR,matchfun=matchfun,
#' formula=list(m~lM + E,s~1,t0~1),
#' contrasts=list(m=list(lM=ADmat)))
#' # For all parameters that aren't defined in the formula, default values are assumed.
#' # These default values can be found in Appendix A of the EMC paper, or accessed using:
#' LNR()$p_types
#' @export
#'


LNR <- function() {
  list(
    type="RACE",
    c_name = "LNR",
    p_types=c("m" = 1,"s" = log(1),"t0" = log(0)),
    Ntransform=function(x) {
      # Transform to natural scale
      x[,dimnames(x)[[2]] != "m"] <- exp(x[,dimnames(x)[[2]] != "m"])
      x
    },
    # p_vector transform scaling parameter by s=1 assumed in lnr.R
    transform = function(x) x,
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) pars,
    # Random function for racing accumulators
    rfun=function(lR=NULL,pars) {
      if (is.null(lR)) return(rep(TRUE,dim(pars)[1]))
      rLNR(lR,pars)
    },
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dLNR(rt,pars),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pLNR(rt,pars),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
      log_likelihood_race(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
  )
}






