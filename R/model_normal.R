

dNORMAL <- function(rt,pars)
  # density for single accumulator
{
  dnorm(rt,mean=pars[,"mean"],sd=pars[,"sd"])
}
## RNORMAL?

#' Normal density function for use with fMRI GLM
#'
#' @return
#' @export
#'
#' @examples
glm_normal <- function() {list(
  type="MRI",
  p_types=c("sd"), #This is a bit hacky for now
  # Transform to natural scale
  Ntransform=function(x) {

    x[,dimnames(x)[[2]] == "sd"] <- exp(x[,dimnames(x)[[2]] == "sd"])
    x
  },
  # Trial dependent parameter transform
  Ttransform = function(pars,dadm=NULL)
  {
    pars
  },
  # p_vector transform
  transform = function(x) x,
  # Random function for racing accumulators
  rfun=function(lR,pars) rNORMAL(lR,pars),
  # Density function (PDF) for single accumulator
  dfun=function(rt,pars) dNORMAL(rt,pars),
  # Probability function (CDF) for single accumulator
  pfun=function(rt,pars) pNORMAL(rt,pars),
  # Race likelihood combining pfun and dfun
  log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)){
    y <- dadm[,!colnames(dadm) %in% c("subjects", "trials")]
    X <- attr(dadm, "model")$design_matrix[[dadm$subjects[1]]]
    is_sd <- grepl("sd", names(p_vector))
    betas <- p_vector[!is_sd]
    sigma <- exp(p_vector[is_sd])
    return(max(sum(dnorm(y, mean=X %*% betas, sd = sigma, log = T)), min_ll*length(y)))
  }
)
}
