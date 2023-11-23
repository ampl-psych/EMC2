#' lnrMS Model accommodating missing values (truncation and censoring)
#'
#' @return A model list with all the necessary functions to sample
#' @export
MlnrMS <- function(){
  list(
    type="RACE",
    p_types=c("m","s","t0","pContaminant"),
    Ntransform=function(x) {
      # Transform to natural scale
      doexp <- !(dimnames(x)[[2]] %in% c("m","pContaminant"))
      x[,doexp] <- exp(x[,doexp])
      doprobit <- dimnames(x)[[2]] == "pContaminant"
      x[,doprobit] <- pnorm(x[,doprobit])
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
      log_likelihood_race_missing(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
  )
}


#' lnr model where s proportional to m
#'
#' @return
#' @export
#'
#' @examples
lnrMP <- function() {
  list(
    type="RACE",
    c_name = "lnrMS",
    p_types=c("m","s","t0","sp"),
    Ntransform=function(x) {
      # Transform to natural scale
      x[,dimnames(x)[[2]] != "m"] <- exp(x[,dimnames(x)[[2]] != "m"])
      x[,"s"] <- x[,"s"] + exp(x[,"m"])*x[,"sp"]
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
