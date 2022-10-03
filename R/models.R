#' dLNR()
#'
#' density for single accumulator
#'
#' @param rt
#' @param pars
#'
#' @return
#' @export
#'
#' @examples
dLNR <- function(rt,pars){
  rt <- rt - pars[,"t0"]
  out <- numeric(length(rt))
  ok <- rt > 0
  out[ok] <- dlnorm(rt[ok],meanlog=pars[ok,"m"],sdlog=pars[ok,"s"])
  out
}


#' pLNR()
#'
#' Cumulative density for single accumulator
#'
#' @param rt
#' @param pars
#'
#' @return
#' @export
#'
#' @examples
pLNR <- function(rt,pars){
  rt <- rt - pars[,"t0"]
  out <- numeric(length(rt))
  ok <- rt > 0
  out[ok] <- plnorm(rt[ok],meanlog=pars[ok,"m"],sdlog=pars[ok,"s"])
  out

}


#' rLNR()
#'
#' @param lR empty latent response factor lR with one level for each accumulator.
#' @param pars matrix of corresponding parameter values named as in p_types
#' @param p_types
#'
#' @return
#' @export
#'
#' @examples
rLNR <- function(lR,pars,p_types=c("m","s","t0")){
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ",paste(p_types,collapse = " "))
  dt <- matrix(rlnorm(dim(pars)[1],meanlog=pars[,"m"],sdlog=pars[,"s"]),
               nrow=length(levels(lR)))
  R <- apply(dt,2,which.min)
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  # Any t0 difference with lR due to response production time (no effect on race)
  rt <- matrix(pars[,"t0"],nrow=length(levels(lR)))[pick] + dt[pick]
  R <- factor(levels(lR)[R],levels=levels(lR))
  cbind.data.frame(R=R,rt=rt)
}


#' lnrMS()
#'
#' LNR mu and sigma parameterization
#'
#' @return
#' @export
#'
#' @examples
lnrMS <- function() {
list(
    type="RACE",
    p_types=c("m","s","t0"),
    Ntransform=function(x) {
      # Transform to natural scale

      get_p_types <- function(nams)
        unlist(lapply(strsplit(nams,"_"),function(x){x[[1]]}))

      if (!is.matrix(x)) {
        nams <- get_p_types(names(x))
        x[nams != "m"] <- exp(x[nams != "m"])
      } else {
        nams <- get_p_types(dimnames(x)[[2]])
        x[,nams != "m"] <- exp(x[,nams != "m"])
      }
      x
    },
    # p_vector transform scaling parameter by s=1 assumed in lnr.R
    transform = function(x) x,
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) pars,
    # Random function for racing accumulators
    rfun=function(lR,pars) rLNR(lR,pars),
    # Density function (PDF) for single accumulator
    dfun=function(rt,pars) dLNR(rt,pars),
    # Probability function (CDF) for single accumulator
    pfun=function(rt,pars) pLNR(rt,pars),
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10))
      log_likelihood_race(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
  )
  }



