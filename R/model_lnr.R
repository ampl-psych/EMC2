
dLNR <- function(rt,pars){
  rt <- rt - pars[,"t0"]
  out <- numeric(length(rt))
  ok <- rt > 0
  out[ok] <- stats::dlnorm(rt[ok],meanlog=pars[ok,"m"],sdlog=pars[ok,"s"])
  out
}

pLNR <- function(rt,pars){
  rt <- rt - pars[,"t0"]
  out <- numeric(length(rt))
  ok <- rt > 0
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
      x[,dimnames(x)[[2]] != "m"] <- exp(x[,dimnames(x)[[2]] != "m"])
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




