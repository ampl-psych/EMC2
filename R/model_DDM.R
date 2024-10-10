find_duplicate_indices <- function(df) {
  df <- as.data.frame(df)
  # Convert dataframe rows to character strings for fast comparison
  row_strings <- do.call(paste, c(df, sep = "_"))
  first_occurrence <- !duplicated(row_strings)
  index_map <- match(row_strings, row_strings[first_occurrence])
  return(index_map)
}

rDDM <- function(lR,pars,precision=5e-3,ok=rep(TRUE,length(lR)))
  # lR is an empty latent response factor lR with one level for each boundary
  # pars is a matrix of parameter values named as in p_types
  # lower is mapped to first level of lR and upper to second
  # test
  # pars=cbind.data.frame(a=c(1,1),v=c(-1,1),t0=c(.2,.2),z=c(.5,.5),d=c(0,0),
  #                       sz=c(0,0),sv=c(0,0),st0=c(0,0),s=c(1,1))
  # lR <- factor(c("left","right"))

{
  bad <- rep(NA,length(lR))
  out <- data.frame(response=bad,rt=bad)
  out_ok <- out[ok,]
  pars <- pars[ok,]
  lR <- lR[ok]
  pars <- as.matrix(pars);
  idx <- find_duplicate_indices(pars)
  for(id in unique(idx)){
    is_id <- which(idx == id)
    first <- is_id[1]
    tmp <- rWDM(N = length(is_id), a = pars[first,"a"]/pars[first, "s"], v = pars[first,"v"]/pars[first, "s"], t0 = pars[first,"t0"],
                       w = pars[first,"Z"], sw = pars[first,"SZ"], sv = pars[first,"sv"]/pars[first, "s"],
                       st0 = pars[first,"st0"], precision = precision)
    tmp <- data.frame(response = tmp$response, rt = tmp$q)
    out_ok[is_id,] <- tmp[sample(nrow(tmp)),]
  }
  out[ok,] <- out_ok
  # cbind.data.frame(R=factor(out[,"response"],levels=c("lower","upper"),labels=levels(lR)),rt=out[,"rt"])
  cbind.data.frame(R=factor(out[,"response"], labels = levels(lR), levels = c("lower", "upper")),rt=out[,"rt"])
}


dDDM <- function(rt,R,pars,precision=5e-3)
  # DDM density for response factor R with rt
  # lower is mapped to first level of R and upper to second
  # test
  # pars=cbind.data.frame(a=c(1,1),v=c(-1,1),t0=c(.2,.2),z=c(.5,.5),d=c(0,0),
  #                       sz=c(0,0),sv=c(0,0),st0=c(0,0),s=c(1,1))
  # R <- factor(c("left","right")); rt=c(1,1)
{
  levels(R) <- c("lower","upper")
  res <- dWDM(rt, response=as.character(R), a = pars[,"a"]/pars[, "s"], v = pars[,"v"]/pars[, "s"], t0 = pars[,"t0"], w = pars[,"Z"],
       sw = pars[,"SZ"], sv = pars[,"sv"]/pars[, "s"],st0 = pars[,"st0"], precision = precision)$value
  return(res)
}

pDDM <- function(rt,R,pars,precision=5e-3)
  # DDM cdf for response factor R with rt
  # lower is mapped to first level of R and upper to second
{
  levels(R) <- c("lower","upper")
  pWDM(rt, response=as.character(R), a = pars[,"a"]/pars[, "s"], v = pars[,"v"]/pars[, "s"], t0 = pars[,"t0"], w = pars[,"Z"],
              sw = pars[,"SZ"], sv = pars[,"sv"]/pars[, "s"],st0 = pars[,"st0"], precision = precision)$value
}


#' The Diffusion Decision Model
#'
#' Model file to estimate the Diffusion Decision Model (DDM) in EMC2.
#'
#' Model files are almost exclusively used in `design()`.
#'
#' @details
#'
#' Default values are used for all parameters that are not explicitly listed in the `formula`
#' argument of `design()`.They can also be accessed with `DDM()$p_types`.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default**   | **Mapping**                    | **Interpretation**                                            |
#' |-----------|-----------|---------------|-----------|----------------------------|-----------------------------------------------------------|
#' | *v*       | -         | \[-Inf, Inf\]     | 1         |                            | Mean evidence-accumulation rate (drift rate)              |
#' | *a*       | log       | \[0, Inf\]        | log(1)    |                            | Boundary separation                                       |
#' | *t0*      | log       | \[0, Inf\]        | log(0)    |                            | Non-decision time                                         |
#' | *s*       | log       | \[0, Inf\]        | log(1)    |                            | Within-trial standard deviation of drift rate            |
#' | *Z*       | probit    | \[0, 1\]        | qnorm(0.5)| *z* = *Z* x *a*                  | Relative start point (bias)                              |
#' | *SZ*      | probit    | \[0, 1\]        | qnorm(0)  | *sz* = 2 x *SZ* x min(*a* x *Z*, *a* x (1-*Z*)) | Relative between-trial variation in start point       |
#' | *sv*      | log       | \[0, Inf\]        | log(0)    |                            | Between-trial standard deviation of drift rate           |
#' | *st0*     | log       | \[0, Inf\]        | log(0)    |                            | Between-trial variation (range) in non-decision time    |
#'
#' `a`, `t0`, `sv`, `st0`, `s` are sampled on the log scale because these parameters are strictly positive,
#' `Z`, `SZ` and `DP` are sampled on the probit scale because they should be strictly between 0 and 1.
#'
#' `Z` is estimated as the ratio of bias to one boundary where 0.5 means no bias.
#' `DP` comprises the difference in non-decision time for each response option.
#'
#' Conventionally, `sv` is fixed to 1 to satisfy scaling constraints.
#'
#' See Ratcliff, R., & McKoon, G. (2008).
#' The diffusion decision model: theory and data for two-choice decision tasks.
#' *Neural computation, 20*(4), 873-922. doi:10.1162/neco.2008.12-06-420.
#'
#' @return A model list with all the necessary functions for EMC2 to sample
#' @examples
#' design_DDMaE <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#' # For all parameters that are not defined in the formula, default values are assumed
#' # (see Table above).
#'
#' @export
#'

DDM <- function(){
  list(
    type="DDM",
    c_name = "DDM",
    p_types=c("v" = 1,"a" = log(1),"sv" = log(0),"t0" = log(0),"st0" = log(0),"s" = log(1),"Z" = qnorm(0.5),"SZ" = qnorm(0)),
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
      pars[,"SZ"] <- 2*pars[,"SZ"]*apply(cbind(pars[,"Z"],1-pars[,"Z"]),1,min)
      pars <- cbind(pars,z=pars[,"Z"]*pars[,"a"],
                    sz = pars[,"SZ"]*pars[,"a"])
      attr(pars,"ok") <- !( abs(pars[,"v"])> 20 | pars[,"a"]> 10 | pars[,"sv"]> 10 | pars[,"SZ"]> .999 |
             pars[,"t0"] < .05 | pars[,"st0"]>5)
      if (pars[1,"sv"] !=0) attr(pars,"ok") <- attr(pars,"ok") & pars[,"sv"] > .001
      if (pars[1,"SZ"] !=0) attr(pars,"ok") <- attr(pars,"ok") & pars[,"SZ"] > .001
      attr(pars, "ok") <- rep(TRUE, nrow(pars))
      pars
    },
    # Random function
    rfun=function(lR=NULL,pars) {
      ok <- !(abs(pars[,"v"])> 20 | pars[,"a"]> 10 | pars[,"sv"]> 10 | pars[,"SZ"]> .999 |
                 pars[,"st0"]>5 | pars[,"t0"] < .05)
      if (pars[1,"sv"] !=0) attr(pars,"ok") <- attr(pars,"ok") & pars[,"sv"] > .001
      if (pars[1,"SZ"] !=0) attr(pars,"ok") <- attr(pars,"ok") & pars[,"SZ"] > .001
      if (is.null(lR)) ok else rDDM(lR,pars,precision=1e-3,ok)
    },
    # Density function (PDF)
    dfun=function(rt,R,pars) dDDM(rt,R,pars,precision=1e-3),
    # Probability function (CDF)
    pfun=function(rt,R,pars) pDDM(rt,R,pars,precision=1e-3),
    log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)){
      log_likelihood_ddm(p_vector=p_vector, dadm = dadm, min_ll = min_ll)
    }
  )
}
