find_duplicate_indices <- function(df) {
  df <- as.data.frame(df)
  # Convert dataframe rows to character strings for fast comparison
  row_strings <- do.call(paste, c(df, sep = "_"))
  first_occurrence <- !duplicated(row_strings)
  index_map <- match(row_strings, row_strings[first_occurrence])
  return(index_map)
}

# Unfortunately there's some unwanted print statements in the code of rWDM
# From the original package
suppress_output <- function(expr) {
  sink(tempfile())  # Redirect output to a temporary file
  on.exit(sink())   # Ensure sink is reset afterward
  invisible(force(expr))  # Run the expression
}

rDDM <- function(R,pars,ok=rep(TRUE,length(R)), precision=5e-3)
{
  bad <- rep(NA,nrow(pars))
  out <- data.frame(response=bad,rt=bad)
  out_ok <- out[ok,]
  pars <- pars[ok,]
  R <- R[ok]
  pars <- as.matrix(pars);
  idx <- find_duplicate_indices(pars)
  for(id in unique(idx)){
    is_id <- which(idx == id)
    cur_pars <- pars[is_id[1],]
    tmp <- suppress_output(rWDM(N = length(is_id), a = cur_pars["a"]/cur_pars[ "s"], v = cur_pars["v"]/cur_pars[ "s"], t0 = cur_pars["t0"],
                       w = cur_pars["Z"], sw = cur_pars["SZ"], sv = cur_pars["sv"]/cur_pars[ "s"],
                       st0 = cur_pars["st0"], precision = precision))
    tmp <- data.frame(response = tmp$response, rt = tmp$q)
    out_ok[is_id,] <- tmp[sample(nrow(tmp)),]
  }
  out[ok,] <- out_ok
  cbind.data.frame(R=factor(out[,"response"], labels = levels(R), levels = c("lower", "upper")),rt=out[,"rt"])
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
#' Conventionally, `s` is fixed to 1 to satisfy scaling constraints.
#'
#' The DDM was introduced by Ratcliff (1978) with parameters `v`, `a`, `t0`, `s`, and `Z`.
#' Between-trial variability parameters `SZ` and `sv` were introduced in Ratcliff & Rouder (1998). Between trial-variability parameter `st0` was introduced in Ratcliff (2002).
#' Several review papers provide further background on the model (e.g., Forstmann et al., 2016; Ratcliff & McKoon, 2008; Ratcliff, 2013; Voss et al., 2004).
#'
#' The C++ implementation of the model's log likelihood function was adapted from code released with the `WienR` package (Hartmann & Klauer, 2021), available under a GPL-2 | GPL-3 license.
#'
#' @references
#' Forstmann, B. U., Ratcliff, R., & Wagenmakers, E. J. (2016). Sequential sampling models in cognitive neuroscience: Advantages, applications, and extensions. *Annual Review of Psychology*, *67*(1), 641-666. \doi{10.1146/annurev-psych-122414-033645}
#'
#' Hartmann, R., & Klauer, K. C. (2021). Partial derivatives for the first-passage time distribution in Wiener diffusion models. *Journal of Mathematical Psychology*, *103*, 102550. \doi{10.1016/j.jmp.2021.102550}
#'
#' Ratcliff, R. (1978). A theory of memory retrieval. *Psychological Review*, *85*(2), 59. \doi{10.1037/0033-295X.85.2.59}
#'
#' Ratcliff, R. (2002). A diffusion model account of response time and accuracy in a brightness discrimination task: Fitting real data and failing to fit fake but plausible data. *Psychonomic Bulletin & Review*, *9*(2), 278-291. \doi{10.3758/BF03196302}
#'
#' Ratcliff, R. (2013). Parameter variability and distributional assumptions in the diffusion model. *Psychological Review*, *120*(1), 281. \doi{10.1037/a0030775}
#'
#' Ratcliff, R., & Rouder, J. N. (1998). Modeling response times for two-choice decisions. *Psychological Science*, *9*(5), 347-356. \doi{10.1111/1467-9280.00067}
#'
#' Ratcliff, R., & McKoon, G. (2008). The diffusion decision model: theory and data for two-choice decision tasks. *Neural computation*, *20*(4), 873-922. \doi{10.1162/neco.2008.12-06-420}
#'
#' Voss, A., Rothermund, K., & Voss, J. (2004). Interpreting the parameters of the diffusion model: An empirical validation. *Memory & Cognition*, *32*, 1206-1220. \doi{10.3758/BF03196893}
#'
#' @return A model list with all the necessary functions for EMC2 to sample
#'
#' @examples
#' design_DDMaE <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#' # For all parameters that are not defined in the formula, default values are assumed
#' # (see Table above).
#'
#' @export
DDM <- function(){
  list(
    c_name = "DDM",
    type="DDM",
    p_types=c("v" = 1,"a" = log(1),"sv" = log(0),"t0" = log(0),"st0" = log(0),"s" = log(1),"Z" = qnorm(0.5),"SZ" = qnorm(0)),
    # Trial dependent parameter transform
    transform=list(func=c(v = "identity",a = "exp",sv = "exp",t0 = "exp",
                          st0 = "exp",s = "exp",Z = "pnorm",SZ = "pnorm")),
    bound=list(minmax=cbind(v=c(-20,20),a=c(0,10),Z=c(.01,.99),t0=c(0.05,Inf),
                            sv=c(.01,10),s=c(0,Inf),SZ=c(.01,.99),st0=c(0,.5)),
               exception=c(sv=0,SZ=0,st0=0)),
    Ttransform = function(pars,dadm) {
      pars[,"SZ"] <- 2*pars[,"SZ"]*apply(cbind(pars[,"Z"],1-pars[,"Z"]),1,min)
      pars <- cbind(pars,z=pars[,"Z"]*pars[,"a"], sz = pars[,"SZ"]*pars[,"a"])
      pars
    },
    # Random function
    rfun=function(data=NULL,pars) rDDM(data$R,pars, attr(pars, "ok")),
    # Density function (PDF)
    dfun=function(rt,R,pars) dDDM(rt,R,pars),
    # Probability function (CDF)
    pfun=function(rt,R,pars) pDDM(rt,R,pars),
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10)){
      log_likelihood_ddm(pars=pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}


#' The DDM with response omissions
#'
#' @details
#' See ?DDM for most details
#'
#' pContaminant is a probit scaled parameter for non-process (contaminant) miss.
#' For background see :
#'
#' Damaso, K. A. M., Castro, S. C., Todd, J., Strayer, D. L., Provost, A.,
#' Matzke, D., & Heathcote, A. (2021). A cognitive model of response omissions
#' in distraction paradigms. Memory & Cognition, 1–17.
#' https://doi.org/10.3758/s13421-021-01265-z
#'
#' Note: This model does not have an Rcpp version.
#'
#' @return A list defining the cognitive model
#' @export

MDDM <- function(){
  list(
    type="DDM",
    p_types=c("v" = 1,"a" = log(1),"sv" = log(0),"t0" = log(0),"st0" = log(0),
      "s" = log(1),"Z" = qnorm(0.5),"SZ" = qnorm(0),pContaminant=qnorm(0)),
    # Trial dependent parameter transform
    transform=list(func=c(v = "identity",a = "exp",sv = "exp",t0 = "exp",
      st0 = "exp",s = "exp",Z = "pnorm",SZ = "pnorm",pContaminant="pnorm")),
    bound=list(minmax=cbind(v=c(-20,20),a=c(0,10),Z=c(.01,.99),t0=c(0.05,Inf),
      sv=c(.01,10),s=c(0,Inf),SZ=c(.01,.99),st0=c(0,.5),pContaminant=c(0,1)),
               exception=c(sv=0,SZ=0,st0=0, pContaminant=0)),
    Ttransform = function(pars,dadm) {
      pars[,"SZ"] <- 2*pars[,"SZ"]*apply(cbind(pars[,"Z"],1-pars[,"Z"]),1,min)
      pars <- cbind(pars,z=pars[,"Z"]*pars[,"a"], sz = pars[,"SZ"]*pars[,"a"])
      pars
    },
    # Random function
    rfun=function(data=NULL,pars) rDDM(data$R,pars, attr(pars, "ok")),
    # Density function (PDF)
    dfun=function(rt,R,pars) dDDM(rt,R,pars),
    # Probability function (CDF)
    pfun=function(rt,R,pars) pDDM(rt,R,pars),
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10)){
      log_likelihood_ddm_missing(p_vector=pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}


#### GNG ----

#' The GNG (go/nogo) Diffusion Decision Model (DDMGNGnoC)
#'
#' In the GNG paradigm one of the two possible choices results in a response
#' being withheld (a non-response), which is indicated in the data by an NA for
#' the rt, with the corresponding level of the R (response) factor still being
#' specified. For example, suppose the go response is coded as "yes" and nogo is
#' coded as "no", then for a non-response (R,rt) = ("no",NA) and for a response
#' e.g., (R,rt) = ("yes",1.36). The GNG paradigm must also have a response
#' # window (i.e., a length of time, TIMEOUT period, after which withholding is
#' assumed).
#'
#' The model used is described in the following paper, with the addition of
#' modeling the TIMEOUT (which is considered but not used in this paper).
#'
#' Gomez, P., Ratcliff, R., & Perea, M. (2007). A Model of the Go/No-Go Task.
#' Journal of Experimental Psychology: General, 136(3), 389–413.
#' https://doi.org/10.1037/0096-3445.136.3.389
#'
#' The likelihood of non-responses requires and evaluation of the DDM cdf,
#' specifically 1 - p(hitting the yes boundary before TIMEOUT).
#'
#' To use these models three functions must be supplied in the design's function
#' argument with the names TIMEOUT, Rnogo and Rgo. For example, assuming a
#' 2.5 second timeout, and R factor with levels c("no","yes") and "no" mapping
#' to a non-response.
#'
#' TIMEOUT=function(d)rep(2.5,nrow(d))
#' Rnogo=function(d)factor(rep("no",nrow(d)),levels=c("no","yes"))
#' Rgo=function(d)factor(rep("yes",nrow(d)),levels=c("no","yes")))
#'
#' See the help for DDM for further details. At present this model is not fully
#' implemented in C, so is a little slower to use than the DDM, but not greatly.
#'
#' @return A model list with all the necessary functions to sample
#' @examples
#' dGNG <- design(Rlevels = c("left","right"),
#'                factors=list(subjects=1,S=c("left","right")),
#'                functions=list(
#'                TIMEOUT=function(d)rep(2.5,nrow(d)),
#'                # no go response level
#'                Rnogo=function(d)factor(rep("left",nrow(d)),levels=c("left","right")),
#'                # go response level
#'                Rgo=function(d)factor(rep("right",nrow(d)),levels=c("left","right"))),
#'                formula=list(v~S,a~1, Z~1, t0~1),
#'                model=DDMGNG)
#'
#' p_vector <- sampled_pars(dGNG)
#' @export
DDMGNG <- function(){
  list(
    type="DDM",
    p_types=c("v" = 1,"a" = log(1),"sv" = log(0),"t0" = log(0),"st0" = log(0),
              "s" = log(1),"Z" = qnorm(0.5),"SZ" = qnorm(0)),
    # Trial dependent parameter transform
    transform=list(func=c(v = "identity",a = "exp",sv = "exp",t0 = "exp",
                          st0 = "exp",s = "exp",Z = "pnorm",SZ = "pnorm")),
    bound=list(minmax=cbind(v=c(-20,20),a=c(0,10),Z=c(.001,.999),t0=c(0.05,Inf),
                            sv=c(.01,10),s=c(0,Inf),SZ=c(.001,.999),st0=c(0,.5)),
               exception=c(sv=0,SZ=0,st0=0)),
    Ttransform = function(pars,dadm) {
      pars[,"SZ"] <- 2*pars[,"SZ"]*apply(cbind(pars[,"Z"],1-pars[,"Z"]),1,min)
      pars <- cbind(pars,z=pars[,"Z"]*pars[,"a"], sz = pars[,"SZ"]*pars[,"a"],
                    TIMEOUT=dadm$TIMEOUT,Rnogo=as.numeric(dadm$Rnogo))
      pars
    },
    # Random function
    rfun=function(data,pars) {
      out <- rDDM(data$R,pars, attr(pars, "ok"))
      out$rt[out$rt>pars[,"TIMEOUT"]] <- NA
      out$rt[as.numeric(out$R)==pars[,"Rnogo"]] <- NA
      out
    },
    # Density function (PDF)
    dfun=function(rt,R,pars) dDDM(rt,R,pars),
    # Probability function (CDF)
    pfun=function(rt,R,pars) pDDM(rt,R,pars),
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10)){
      log_likelihood_ddmgng(pars=pars, dadm = dadm, model = model, min_ll = min_ll)
    }
  )
}

