#' RDM_B parameterization with missing values
#'
#' @return A list defining the cognitive model
#' @export
#'
MrdmB <- function(){
  list(
    type="RACE",
    p_types=c("v","B","A","t0","s","pContaminant"),
    # Transform to natural scale
    Ntransform=function(x) {
      # transform parameters back to real line
      doprobit <- dimnames(x)[[2]] == "pContaminant"
      x[,doprobit] <- pnorm(x[,doprobit])
      x[,!doprobit] <- exp(x[,!doprobit])
      x
    },
    # p_vector transform
    transform = function(x) x,
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      attr(pars,"ok") <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0)
      pars
    },
    # Random function for racing accumulators
    rfun=function(lR=NULL,pars) {
      ok <- (pars[,"t0"] > .05) & ((pars[,"A"] > 1e-6) | pars[,"A"] == 0)
      if (is.null(lR)) ok else rRDM(lR,pars,ok=ok)
    },
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dRDM(rt,pars),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pRDM(rt,pars),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
      log_likelihood_race_missing(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
  )
}

