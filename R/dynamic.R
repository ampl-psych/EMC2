#### Dynfuns ----

#### Learning functions ----

run_delta_i <- function(q0,alpha,target,winner,return_extras=FALSE) {

  q <- pe <- numeric(length(target))
  q[1] <- q0
  alpha <- pnorm(alpha)
  for (i in 2:length(q)) {
    if (is.na(target[i-1]) | winner[i-1] != 1) {
      pe[i-1] <- NA
      q[i] <- q[i-1]
    } else {
      pe[i-1] <- target[i-1]-q[i-1]
      q[i] <- q[i-1] + alpha[i-1]*pe[i-1]
    }
  }
  if (return_extras) {
    pe <- target - q
    cbind(q=q,pe=pe)
  } else q
}

run_delta <- function(q0,p,target,return_extras=FALSE) {
  if (any(dimnames(target)[[2]]=="winner")) {
    winner <- as.numeric(target[,1]) # R always first
    target <- target[,-1,drop=FALSE]
  } else winner <- rep(1,nrow(target))
  if (!return_extras) out <- matrix(ncol=ncol(target),nrow=nrow(target)) else
    out <- array(dim=c(ncol(target),nrow(target),2),
                 dimnames=list(dimnames(target)[[2]],NULL,c("q","pe")))
  for (i in 1:ncol(target)) {
    if (!return_extras)
      out[,i] <- run_delta_i(q0,p,target[,i],winner,return_extras) else
        out[i,,] <- run_delta_i(q0,p,target[,i],winner,return_extras)
  }
  out
}



run_delta2_i <- function(q0,alphaSlow,dFast,dSwitch,target,
                         winner,return_extras=FALSE) {

  qOut <- qFast <- qSlow <- peFast <- peSlow <- numeric(length(target))
  qOut[1] <- qFast[1] <- qSlow[1] <- q0
  alphaFast <- pnorm(alphaSlow + dFast)
  alphaSlow <- pnorm(alphaSlow)
  dSwitch <- pnorm(dSwitch)
  for (i in 2:length(target)) {
    if (is.na(target[i-1]) | winner[i-1] != 1) {
      peSlow[i-1] <- peFast[i-1] <- NA
      qFast[i] <- qFast[i-1]
      qSlow[i] <- qSlow[i-1]
    } else {
      peFast[i-1] <- target[i-1]-qFast[i-1]
      peSlow[i-1] <- target[i-1]-qSlow[i-1]
      qFast[i] <- qFast[i-1] + alphaFast[i-1]*peFast[i-1]
      qSlow[i] <- qSlow[i-1] + alphaSlow[i-1]*peSlow[i-1]
      if (abs(qFast[i]-qSlow[i])>dSwitch[i])
        qOut[i] <- qFast[i] else qOut[i] <- qSlow[i]
    }
  }
  if (return_extras) {
    peFast <- target-qFast
    peSlow <- target-qSlow
    pe <- qOut-target
    # Didnt return pe for Fast and Slow, can reconstruct
    cbind(q=qOut,pe=pe,fast=qFast,slow=qSlow)
  } else qOut
}


run_delta2 <- function(q0,p,target,return_extras=FALSE) {

  if (any(dimnames(target)[[2]]=="winner")) {
    winner <- as.numeric(target[,1]) # winner always first
    target <- target[,-1,drop=FALSE]
  } else winner <- rep(1,nrow(target))
  if (!return_extras) out <- matrix(ncol=ncol(target),nrow=nrow(target)) else
    out <- array(dim=c(ncol(target),nrow(target),4),
                 dimnames=list(dimnames(target)[[2]],NULL,c("q","pe","fast","slow")))

  for (i in 1:ncol(target)) {
    if (!return_extras)
      out[,i] <- run_delta2_i(q0,p[,1],p[,2],p[,3],target[,i],winner,return_extras) else
        out[i,,] <- run_delta2_i(q0,p[,1],p[,2],p[,3],target[,i],winner,return_extras)
  }
  out
}


run_d1b <- function(base,scale,q0,alpha,target,return_extras=FALSE) {

  q <- pe <- numeric(length(target))
  q[1] <- q0
  alpha <- pnorm(alpha)
  for (i in 2:length(q)) {
    # if (is.na(target[i-1])) { # remove? not needed for RL
    #   pe[i-1] <- NA
    #   q[i] <- q[i-1]
    #   base[i] <- base[i-1]
    # } else {
    target[i-1] <- target[i-1]*base[i-1]
    pe[i-1] <- target[i-1]-q[i-1]
    q[i] <- q[i-1] + alpha[i-1]*pe[i-1]
    base[i] <- base[i] + scale[i]*q[i]
    # }
  }
  if (return_extras) {
    pe <- target - q
    cbind(q=q,pe=pe)
  } else cbind(q=q,b=base)
}



rep_dadm <- function(out1,cs1,Rlevs,lR)
  # Expand out1 to fit with dadm if required.
{
  if (is.array(out1)) {
    if (length(dim(out1)==3) && dim(out1)[1]==1) out1 <- out1[1,,]
    if (is.matrix(out1)) {
      out <- matrix(nrow=length(cs1),ncol=dim(out1)[2],
                    dimnames=list(NULL,dimnames(out1)[[2]]))
      if (is.na(lR) || lR) {
        for (i in Rlevs) out[cs1==i,] <- out1
      } else return(out1)
    } else {
      # This turns it into a matrix ...
      #   out <- matrix(nrow = length(cs1), ncol = prod(dim(out1)[c(1,3)]),
      #     dimnames = list(NULL, as.vector(outer(dimnames(out1)[[3]],
      #                     dimnames(out1)[[1]], paste, sep = "_"))))
      #   for (j in 1:dim(out1)[1])
      #     out[,(1+(j-1)*dim(out1)[3]):(j*dim(out1)[3])] <- out1[j,,]
      out <- out1
    }
  } else {
    out <- numeric(length(cs1))
    if (is.na(lR) || lR) {
      for (i in Rlevs) out[cs1==i] <- out1
    } else return(out1)
  }
  out
}

# tmp <- array(dim=c(2,10000,2),dimnames=list(1:2,NULL,dimnames(out1)[[3]]))
# tmp[1,,] <- out1[1,,]
# tmp[2,,] <- out1[1,,]
# out1=tmp

dynfuns <- function(base,dp,dyntype,cs,lR,maptype="lin",
                    return_extras=FALSE) {
  Rlevs <- levels(cs[,1])
  isRL <- any(names(cs)=="winner")
  if (isRL) {
    remove_cs <- 2
  } else remove_cs <- 1
  if (is.na(lR) || lR) {
    is1 <- cs[,1]==Rlevs[1]
  } else is1 <- !logical(nrow(cs))
  if (is.na(lR)) { # learning and adaptive
    is1map <- !logical(nrow(cs))
  } else is1map <- is1
  if (is.matrix(dp)) {
    out1 <- switch(dyntype, # kernel
                   ld = -cs[is1,2],
                   li = cs[is1,2],
                   ed = exp(-exp(dp[is1,2])*cs[is1,2]),
                   ei = (1-exp(-exp(dp[is1,2])*cs[is1,2])),
                   pd = (1+cs[is1,2])^(-exp(dp[is1,2])),
                   pi = (1-(1+cs[is1,2])^(-exp(dp[is1,2]))),
                   p2 = dp[is1,1]*cs[is1,2] + dp[is1,2]*cs[is1,2]^2,
                   p3 = dp[is1,1]*cs[is1,2] + dp[is1,2]*cs[is1,2]^2 + dp[is1,3]*cs[is1,2]^3,
                   p4 = dp[is1,1]*cs[is1,2] + dp[is1,2]*cs[is1,2]^2 + dp[is1,3]*cs[is1,2]^3 +
                     dp[is1,4]*cs[is1,2]^4,
                   d1U = run_delta(pnorm(dp[is1,2][1]),dp[is1,3],
                                   cs[is1,-c(1:remove_cs),drop=FALSE],return_extras=return_extras),
                   d2U = run_delta2(pnorm(dp[is1,2][1]),dp[is1,3:5],
                                    cs[is1,-c(1:remove_cs),drop=FALSE],return_extras=return_extras),
                   d1P = run_delta(exp(dp[is1,2][1]),dp[is1,3],
                                   cs[is1,-c(1:remove_cs),drop=FALSE],return_extras=return_extras),
                   d2P = run_delta2(exp(dp[is1,2][1]),dp[is1,3:5],
                                    cs[is1,-c(1:remove_cs),drop=FALSE],return_extras=return_extras),
                   d1 = run_delta(dp[is1,2][1],dp[is1,3],
                                  cs[is1,-c(1:remove_cs),drop=FALSE],return_extras=return_extras),
                   d2 = run_delta2(dp[is1,2][1],dp[is1,3:5],
                                   cs[is1,-c(1:remove_cs),drop=FALSE],return_extras=return_extras),
                   d1b = run_d1b(base=base[is1],scale=dp[is1,1],q0=dp[is1,2][1],alpha=dp[is1,3],
                                 cs[is1,-c(1:remove_cs)],return_extras=return_extras)
    )
    if (!return_extras) {
      if (isRL) out1 <- out1[cbind(1:nrow(out1),cs[,2])]
      if (length(base)>1) base <- base[is1map]
      if (is.na(lR)) {
        out1 <- rep_dadm(out1,cs[,1],Rlevs,lR)
        lR <- FALSE
      }
      if (is.null(maptype)) {
        out1 <- out1[,"b"]
      } else out1 <- switch(maptype, # map to parameter
                            lin = base + dp[is1map,1]*out1,
                            plin = exp(base) + exp(dp[is1map,1])*out1,
                            add = base + out1,
                            ucent = base + dp[is1map,1]*(out1-.5),
                            identity = out1
      )
    }
  } else {
    n <- sum(is1)
    out1 <- switch(dyntype, # kernel
                   ld = -cs[is1,2],
                   li = cs[is1,2],
                   ed = exp(-exp(dp[2])*cs[is1,2]),
                   ei = (1-exp(-exp(dp[2])*cs[is1,2])),
                   pd = (1+cs[is1,2])^(-exp(dp[2])),
                   pi = (1-(1+cs[is1,2])^(-exp(dp[2]))),
                   p2 = dp[1]*cs[is1,2] + dp[2]*cs[is1,2]^2,
                   p3 = dp[1]*cs[is1,2] + dp[2]*cs[is1,2]^2 + dp[3]*cs[is1,2]^3,
                   p4 = dp[1]*cs[is1,2] + dp[2]*cs[is1,2]^2 + dp[3]*cs[is1,2]^3 +
                     dp[4]*cs[is1,2]^4,
                   d1U = run_delta(pnorm(dp[2]),rep(dp[3],n),
                                   cs[is1,-c(1:remove_cs),drop=FALSE],return_extras=return_extras),
                   d2U = run_delta2(pnorm(dp[2]),matrix(rep(dp[3:5],each=n*3),ncol=3),
                                    cs[is1,-c(1:remove_cs),drop=FALSE],return_extras=return_extras),
                   d1P = run_delta(exp(dp[2]),rep(dp[3],n),
                                   cs[is1,-c(1:remove_cs),drop=FALSE],return_extras=return_extras),
                   d2P = run_delta2(exp(dp[2]),matrix(rep(dp[3:5],each=n*3),ncol=3),
                                    cs[is1,-c(1:remove_cs),drop=FALSE],return_extras=return_extras),
                   d1 = run_delta(dp[2],rep(dp[3],n),
                                  cs[is1,-c(1:remove_cs),drop=FALSE],return_extras=return_extras),
                   d2 = run_delta2(dp[2],matrix(rep(dp[3:5],each=n*3),ncol=3),
                                   cs[is1,-c(1:remove_cs),drop=FALSE],return_extras=return_extras),
                   d1b = run_d1b(base=base[is1],scale=dp[is1,1],q0=dp[2],alpha=dp[3],
                                 cs[is1,-c(1:remove_cs)],return_extras=return_extras)
    )
    if (!return_extras) { # map to parameter
      if (isRL) out1 <- out1[cbind(1:nrow(out1),cs[,2])]
      if (length(base)>1) base <- base[is1map]
      if (is.na(lR)) {
        out1 <- rep_dadm(out1,cs[,1],Rlevs,lR)
        lR <- FALSE
      }
      if (is.null(maptype)) {
        out1 <- out1[,"b"]
      } else out1 <- switch(maptype, # map to parameter
                            lin = base + dp[1]*out1,
                            plin = exp(base) + exp(dp[1])*out1,
                            add = base + out1,
                            ucent = base + dp[1]*(out1-.5),
                            identity <- out1
      )
    }
  }
  rep_dadm(out1,cs[,1],Rlevs,lR)
}

#### Utility functions ----

do_adaptive <- function(pars,dadm)
  # Used by adaptive to map parameters
{
  for (j in names(attr(dadm,"adaptive"))) {
    dyntype <- attr(dadm,"adaptive")[[j]]$dyntype
    aptypes <- attr(dadm,"adaptive")[[j]]$aptypes
    covnames <- attr(dadm,"adaptive")[[j]]$covnames
    if (any(covnames=="winner")) { # For RL
      if (!any(names(dadm)=="winner")) dadm$winner <- dadm$R==dadm$lR
      lR <- FALSE
    } else lR <- attr(dadm,"adaptive")[[j]]$lR
    if (!is.null(attr(dadm,"adaptive")[[j]]$pcovnames) &&
        j!=attr(dadm,"adaptive")[[j]]$pcovnames) { # Works with standard learning rules
      Rlevs <- levels(dadm$lR)
      # Average pcovnames parameter over accumulators and update covariate
      dadm[dadm$lR==Rlevs[1],covnames] <- dadm[dadm$lR==Rlevs[1],covnames]*
        apply(matrix(pars[,attr(dadm,"adaptive")[[j]]$pcovnames],nrow=length(Rlevs)),2,mean)
    }
    maptype <- attr(dadm,"adaptive")[[j]]$maptype
    pars[,j] <- dynfuns(pars[,j],pars[,aptypes],dyntype,dadm[,c("lR",covnames),drop=FALSE],
                        lR,maptype)
  }
  pars
}


#' Get adaptive parameter names
#'
#' @param adaptive A list containing a dyntype specification
#'
#' @return List of names generated from an adaptive "dyntype" specification
#' @export
adaptive_pars <- function(adaptive) {
  pnams <- check_adaptive(adaptive)
  lapply(pnams,\(x) x$aptypes)
}

dynamic_names <- function(nam=NULL) {
  # Return names associated with currently supported dynamic types
  if (is.null(nam)) {
    list(mtypes=c("lin","plin","add","ucent"),
         dtypes=c("ld","li","ei","ed","pi","pd","p2","p3","p4",
                  "d1U","d2U","d1P","d2P","d1","d2","d1b")
    )
  } else if (nam=="learn") { # number of q values to track
    setNames(c(1,3,1,3,1,3,1),c("d1U","d2U","d1P","d2P","d1","d2","d1b"))
  } else if (nam=="poly") {
    c("p2","p3","p4")
  } else switch(nam,
                li=c("1li"),
                ld=c("1ld"),
                ei=c("1ei","2ei"),
                ed=c("1ed","2ed"),
                pi=c("1pi","2pi"),
                pd=c("1pd","2pd"),
                p2=c("1p2","2p2"),
                p3=c("1p3","2p3","3p3"),
                p4=c("1p4","2p4","3p4","4p4"),
                d1U=c("1d1U","2d1U","3d1U"),
                d2U=c("1d2U","2d2U","3d2U","4d2U","5d2U"),
                d1P=c("1d1P","2d1P","3d1P"),
                d2P=c("1d2P","2d2P","3d2P","4d2P","5d2P"),
                d1=c("1d1","2d1","3d1"),
                d2=c("1d2","2d2","3d2","4d2","5d2"),
                d1b=c("1d1b","2d1b","3d1b")
  )
}

set_dynamic_bounds <- function(di) {
  switch(di$dyntype,
         ei=list(index=2,lower=c(.001),upper=c(NA),exception=list(list=NA)),
         ed=list(index=2,lower=c(.001),upper=c(NA),exception=list(list=NA)),
         pi=list(index=2,lower=c(.001),upper=c(NA),exception=list(list=NA)),
         pd=list(index=2,lower=c(.001),upper=c(NA),exception=list(list=NA)),
         d1U=list(index=3,lower=c(.001),upper=c(NA),exception=list(list=NA)),
         d2U=list(index=4:5,lower=c(.001,.001),upper=c(NA,NA),exception=list(list=NA,list=NA)),
         d1P=list(index=3,lower=c(.001),upper=c(NA),exception=list(list=NA)),
         d2P=list(index=4:5,lower=c(.001,.001),upper=c(NA,NA),exception=list(list=NA,list=NA)),
         d1=list(index=3,lower=c(.001),upper=c(NA),exception=list(list=NA)),
         d2=list(index=4:5,lower=c(.001,.001),upper=c(NA,NA),exception=list(list=NA,list=NA)),
         d1b=list(index=3,lower=c(.001),upper=c(NA),exception=list(list=NA)),
  )
}

check_dynamic <- function(dynamic, covariates) {
  if (is.null(covariates)) stop("must specify covariates when using dynamic")
  covnames <- unlist(lapply(dynamic,function(x)x$covnames))
  if (!all(covnames %in% names(covariates))) stop("dynamic argument has covnames not in covariates")
  if (any(duplicated(names(dynamic)))) stop("Duplicate names in dynamic")
  for (i in names(dynamic)) {
    dnams <- names(dynamic[[i]])
    nams <- dnams[!(dnams %in% c("dnames","maptype","transform","equal_accumulators","shared","S"))]
    if (length(nams)!=2 || (!all(nams %in% c("covnames","dyntype"))))
      stop("dynamic argument ",i," must have names: covnames, dyntype")
    if (!(dynamic[[i]]$dyntype %in% dynamic_names()$dtypes))
      stop("dynamic argument ",i," has an unsupported entry in dyntype, use one of:\n",
           paste(dynamic_names()$dtypes,collapse=","))
    if (dynamic[[i]]$dyntype=="d1b") stop("Learning rule d1b only supported for adaptive")
    if (is.null(dynamic[[i]]$dnames))
      dynamic[[i]]$dnames <- dynamic_names(dynamic[[i]]$dyntype)
    if (length(dynamic[[i]]$dnames) != length(dynamic_names(dynamic[[i]]$dyntype)))
      stop("dnames argument must be length ",length(dynamic_names(dynamic[[i]]$dyntype)))
    if (is.null(dynamic[[i]]$shared))
      dynamic[[i]]$shared <- rep(FALSE,length(dynamic[[i]]$dnames))
    if (length(dynamic[[i]]$shared)==1)
      dynamic[[i]]$shared <- rep(dynamic[[i]]$shared,length(dynamic[[i]]$dnames))
    if (!all(is.logical(dynamic[[i]]$shared)) |
        length(dynamic[[i]]$shared) != length(dynamic[[i]]$dnames))
      stop("shared argument must be a logical length 1 or ",length(dynamic[[i]]$dnames))
    if (is.null(dynamic[[i]]$maptype)) {
      if (dynamic[[i]]$dyntype %in% dynamic_names("poly"))
        dynamic[[i]]$maptype <- "add" else
          dynamic[[i]]$maptype <- "lin"
    }
    if (is.null(dynamic[[i]]$transform)) dynamic[[i]]$transform <- TRUE
    if (!(dynamic[[i]]$maptype %in% dynamic_names()$mtypes))
      stop("dynamic argument ",i," has an unsupported entry in maptype, use one of:\n",
           paste(dynamic_names()$mtypes,collapse=","))
    dynamic[[i]]$dpnames <- dynamic[[i]]$dnames
    dynamic[[i]]$dpnames[!dynamic[[i]]$shared] <- paste0(i,dynamic[[i]]$dnames[!dynamic[[i]]$shared])
    if (!is.null(dynamic[[i]]$S)) {
      if (!(dynamic[[i]]$dyntype %in% names(dynamic_names("learn"))))
        stop("Must specify a learning dyntype when S is supplied.")
      dynamic[[i]]$lR1 <- FALSE
    } else if (!is.null(dynamic[[i]]$equal_accumulators)) {
      if (length(dynamic[[i]]$equal_accumulators)!=1 || !is.logical(dynamic[[i]]$equal_accumulators))
        stop("equal_accumulators must be a logical length 1")
      dynamic[[i]]$lR1 <- dynamic[[i]]$equal_accumulators
      dynamic[[i]]$equal_accumulators <- NULL
    } else {
      dynamic[[i]]$lR1 <- TRUE
    }
    if (is.null(dynamic[[i]]$bounds)) dynamic[[i]] <- set_dynamic_bounds(dynamic[[i]])
  }
  dynamic
}

check_pars_dynamic <- function(dynamic, p_vector, design){
  if (!all(names(dynamic) %in% names(add_constants(p_vector,design$constants))))
    stop("dynamic argument has parameter names not in the model")
  pt <- lapply(attr(attr(design, "p_vector"), "map"),
               function(x) unlist(dimnames(x)[[2]]))
  isd <- setNames(!logical(length(pt)),names(pt))
  dp <- names(lapply(dynamic, function(x) names(x$dpnames)))
  dp <- dp[unlist(lapply(dynamic,function(x){!x$transform}))]
  isd <- unlist(lapply(pt, function(x) any(dp %in% x)))
  # If name included do transformation
  attr(dynamic,"transform_names") <- names(isd[!isd])
  for (i in names(dynamic)) {
    dynamic[[i]]$ptype <- names(pt)[unlist(lapply(pt,function(x){i %in% x}))]
  }
  return(dynamic)
}


check_adaptive <- function(adaptive,model, covariates,formula=NULL) {
  if (!all(names(adaptive) %in% names(model()$p_types)))
    stop("adaptive argument has a parameter type names not in the model")
  if (is.null(covariates))
    stop("must specify covariates when using adaptive")
  covnames <- unlist(lapply(adaptive,function(x)x$covnames))
  if (!all(covnames %in% names(covariates)))
    stop("adaptive argument has covnames not in covariates")
  if (any(duplicated(names(adaptive)))) stop("Duplicate names in adaptive")
  for (i in names(adaptive)) {
    anams <- names(adaptive[[i]])
    nams <- anams[!(anams %in% c("anames","maptype","transform","equal_accumulators","shared","S","pcovnames"))]
    if (length(nams)!=2 || !all(nams %in% c("covnames","dyntype"))) {
      stop("adaptive argument ",i," must have names: covnames, dyntype")
    }
    if (!(adaptive[[i]]$dyntype %in% dynamic_names()$dtypes))
      stop("adaptive argument ",i," has an unsupported entry in dyntype, use one of:\n",
           paste(dynamic_names()$dtypes,collapse=","))
    if (is.null(adaptive[[i]]$maptype) & adaptive[[i]]$dyntype != "d1b") {
      if (adaptive[[i]]$dyntype %in% dynamic_names("poly"))
        adaptive[[i]]$maptype <- "add" else
          adaptive[[i]]$maptype <- "lin"
    }
    if (!is.null(adaptive[[i]]$maptype) && !(adaptive[[i]]$maptype %in% dynamic_names()$mtypes))
      stop("adaptive argument ",i," has an unsupported entry in maptype, use one of:\n",
           paste(dynamic_names()$mtypes,collapse=","))
    if (is.null(adaptive[[i]]$transform)) adaptive[[i]]$transform <- FALSE
    if (is.null(adaptive[[i]]$anames))
      adaptive[[i]]$anames <- dynamic_names(adaptive[[i]]$dyntype)
    if (length(adaptive[[i]]$anames) != length(dynamic_names(adaptive[[i]]$dyntype)))
      stop("anames argument must be length ",length(dynamic_names(adaptive[[i]]$dyntype)))
    if (is.null(adaptive[[i]]$shared))
      adaptive[[i]]$shared <- rep(FALSE,length(adaptive[[i]]$anames))
    if (length(adaptive[[i]]$shared)==1)
      adaptive[[i]]$shared <- rep(adaptive[[i]]$shared,length(adaptive[[i]]$anames))
    if (!all(is.logical(adaptive[[i]]$shared)) |
        length(adaptive[[i]]$shared) != length(adaptive[[i]]$anames))
      stop("shared argument must be a logical length 1 or ",length(adaptive[[i]]$anames))
    adaptive[[i]]$aptypes <- adaptive[[i]]$anames
    adaptive[[i]]$aptypes[!adaptive[[i]]$shared] <-
      paste(i,adaptive[[i]]$anames[!adaptive[[i]]$shared],sep="")
    if (!is.null(formula)) {
      isin <-  adaptive[[i]]$aptypes %in% unlist(lapply(formula,function(x)all.vars(x)[1]))
      if (!all(isin)) stop("Missing in formula: ",paste(adaptive[[i]]$aptypes[!isin],collapse=","))
    }
    if (!is.null(adaptive[[i]]$S)) {
      if (!(adaptive[[i]]$dyntype %in% names(dynamic_names("learn"))))
        stop("Must specify a learning dyntype when S is supplied.")
      adaptive[[i]]$lR1 <- FALSE
    } else if (!is.null(adaptive[[i]]$equal_accumulators)) {
      if (length(adaptive[[i]]$equal_accumulators)!=1 || !is.logical(adaptive[[i]]$equal_accumulators))
        stop("equal_accumulators must be a logical length 1")
      adaptive[[i]]$lR1 <- adaptive[[i]]$equal_accumulators
      adaptive[[i]]$equal_accumulators <- NULL
    } else {
      adaptive[[i]]$lR1 <- FALSE
    }
    if (adaptive[[i]]$dyntype %in% names(dynamic_names("learn"))) adaptive[[i]]$lR1 <- NA
    if (is.null(adaptive[[i]]$bounds)) adaptive[[i]] <- set_dynamic_bounds(adaptive[[i]])
  }
  attr(adaptive,"aptypes") <- unique(unlist(lapply(adaptive,function(x)x$aptypes),use.names = FALSE))
  adaptive
}

check_pars_adaptive <- function(adaptive, design){
  pt <- lapply(attr(attr(design, "p_vector"), "map"),
               function(x) unlist(dimnames(x)[[2]]))
  istrans <- setNames(!logical(length(pt)),names(pt))
  # Parameters not to transform
  ap <- names(adaptive)[unlist(lapply(adaptive,function(x)!x$transform))]
  # ap <- c(ap,unlist(lapply(adaptive, function(x)
  #   if (!x$transform) x$apnames),use.names=FALSE))
  ap <- c(ap,unlist(lapply(adaptive, function(x) {
    if (!x$transform) x$aptypes}),use.names = FALSE))
  istrans[ap] <- FALSE
  # If name included do transformation
  attr(adaptive,"transform_names") <- names(istrans[istrans])
  return(adaptive)
}




#' Get dynamic parameter, or underlying non-linear components
#'
#' Non-linear components are q values and prediction errors for RL models and
#' kernel for trend models
#'
#' @param pname parameter to provide results for
#' @param p_vector parameter vector
#' @param data data frame providing covariate(s)
#' @param design design with dynamic RL component
#' @param return_dadm return full dadm (accumulators) otherwise just trials
#' @param compress compress dadm? Default FALSE appropriate for RL models
#' @param return_p return dynamic parameter, default FALSE, return non-linear
#'
#' @return dadm/data q values, prediction errors etc. (etc. depends on RL model dyntype)
#' @export
get_dynamic <- function(pname,p_vector,data,design,return_dadm=FALSE,return_p=FALSE,compress=FALSE)
{
  if (!pname %in% names(p_vector)) stop("pname must be a name in p")
  dadm <- design_model(data,design,compress=compress,verbose=FALSE)
  if (!all(sort(names(add_constants(p_vector,attributes(dadm)$constants)))==sort(attr(dadm,"p_names"))))
    stop("p_vector names must be: ",paste(attr(dadm,"p_names"),collapse=", "))
  if (is.null(attr(dadm,"dynamic"))) stop("design must be dynamic")
  if (!(pname %in% names(attr(dadm,"dynamic")))) stop(pname," is not dynamic")
  p_vector <- add_constants(p_vector,attributes(dadm)$constants)
  dyntype <- attr(dadm,"dynamic")[[pname]]$dyntype
  dpnames <- attr(dadm,"dynamic")[[pname]]$dpnames
  covnames <- attr(dadm,"dynamic")[[pname]]$covnames
  if (!return_p & any(covnames=="winner")) return_dadm <- FALSE
  lR <- attr(dadm,"dynamic")[[pname]]$lR
  maptype <- attr(dadm,"dynamic")[[pname]]$maptype
  pm <- p_vector[pname]
  dp <- p_vector[dpnames]
  kernel <- dynfuns(pm,dp,dyntype,dadm[,c("lR",covnames),drop=FALSE],
                    lR,maptype,return_extras=!return_p)
  if (length(dim(kernel))!=3) {
    if (return_dadm) {
      out <- cbind(dadm,kernel)
    } else {
      is1 <- dadm$lR==levels(dadm$lR)[1]
      if (is.matrix(kernel)) kernel <- kernel[is1,] else kernel <- kernel[is1]
      out <- cbind(dadm[is1,],kernel)
    }
    if (return_p) dimnames(out)[[2]][length(dimnames(out)[[2]])] <- pname
    return(out)
  } else {
    return(kernel)
  }
}


#' Get adaptive parameter, or underlying non-linear components
#'
#' Non-linear components are q values and prediction errors for RL models and
#' kernel for trend models
#'
#' @param pname parameter to provide results for
#' @param p_vector parameter vector
#' @param data data frame providing covariate(s)
#' @param design design with dynamic RL component
#' @param return_dadm return full dadm (accumulators) otherwise just trials
#' @param compress compress dadm? Default FALSE appropriate for RL models
#' @param return_p return dynamic parameter, default FALSE, return non-linear
#'
#' @return dadm/data q values, prediction errors etc. (etc. depends on RL model dyntype)
#' @export
get_adaptive <- function(pname,p_vector,data,design,return_dadm=TRUE,return_p=TRUE,compress=FALSE)
{
  if (!pname %in% names(design$model()$p_types)) stop("pname must be a parameter type")
  dadm <- design_model(data,design,compress=compress,verbose=FALSE)
  if (is.null(attr(dadm,"adaptive"))) stop("design must be adaptive")
  if (!(pname %in% names(attr(dadm,"adaptive")))) stop(pname," is not adaptive")
  model <- attributes(dadm)$model
  p <- model()$Ttransform(model()$Ntransform(map_p(
    design$model()$transform(add_constants(p_vector,design$constants)),dadm
  ),attr(design,"transform_names")),dadm)
  covnames <- attr(dadm,"adaptive")[[pname]]$covnames
  isRL <- any(covnames=="winner")
  if (return_p) out <- p[,pname] else {
    dyntype <- attr(dadm,"adaptive")[[pname]]$dyntype
    aptypes <- attr(dadm,"adaptive")[[pname]]$aptypes
    if (!is.null(attr(dadm,"adaptive")[[pname]]$pcovnames) &&
        attr(dadm,"adaptive")[[pname]]$pcovnames!=pname) {
      Rlevs <- levels(dadm$lR)
      # Average pcovnames parameter over accumulators and update covariate
      dadm[dadm$lR==Rlevs[1],covnames] <- dadm[dadm$lR==Rlevs[1],covnames]*
        apply(matrix(p[,attr(dadm,"adaptive")[[pname]]$pcovnames],nrow=length(Rlevs)),2,mean)
    }
    if (isRL) {
      if (!any(names(dadm)=="winner")) dadm$winner <- dadm$R==dadm$lR
      lR <- FALSE
    } else lR <- attr(dadm,"adaptive")[[pname]]$lR
    maptype <- attr(dadm,"adaptive")[[pname]]$maptype
    out <- dynfuns(p[,pname],p[,aptypes],dyntype,dadm[,c("lR",covnames),drop=FALSE],
                   lR,maptype,return_extras=!return_p)
  }
  if (isRL & !return_p) return(out) else if (return_dadm) {
    out <- cbind(dadm,out)
  } else {
    is1 <- dadm$lR==levels(dadm$lR)[1]
    if (is.matrix(out)) out <- out[is1,] else out <- out[is1]
    out <- cbind(dadm[is1,],out)
    names(out)[ncol(out)] <- "kernel"
  }
  if (return_p) dimnames(out)[[2]][length(dimnames(out)[[2]])] <- pname
  out
}

#### Dynamic covariates ----

# One step versions for data simulation,

one_delta <- function(q0,alpha,target) {
  q0 + exp(alpha)*(target-q0)
}


one_delta2 <- function(q0,alphas,target) {
  qSlow <- q0[1] + pnorm(alphas[1])*(target-q0[1])
  qFast <- q0[2] + pnorm(sum(alphas[1:2]))*(target-q0[2])
  if (abs(qFast-qSlow)>pnorm(alphas[3]))
    q <- qFast else q <- qSlow
  c(qSlow=qSlow,qFast=qFast,q=q)
}

one_d <- function(q0s,alphas,target,dyntype) {
  switch(substr(dyntype,1,2),
         d1=one_delta(q0s,alphas,target),
         d2=one_delta2(q0s,alphas,target)
  )
}


dynamic_rfun <- function(dynamic_cv,p_vector,dadm,design)
{

  if (!is.null(design$adaptive) & !is.null(design$dynamic))
    stop("Cant mix dynamic and adaptive")

  nacc <- length(levels(dadm$lR))
  ndadadm <- dadm # no dynamic or adaptive, use in Ttransform
  attr(ndadadm,"dynamic") <- NULL
  attr(ndadadm,"adaptive") <- NULL
  pars <- map_p(design$model()$transform(add_constants(p_vector,design$constants)), dadm)
  data_out <- data.frame()
  unq_subs <- unique(dadm$subjects)
  k <- 1
  for(sub in unique(dadm$subjects)){
    idx <- dadm$subjects == sub
    tmp_p_vector <- p_vector[k,,drop = F]
    k <- k + 1
    tmp_dadm <- dadm[idx,]
    ntrials <- nrow(tmp_dadm)/nacc
    tmp_dadm$subjects <- factor(tmp_dadm$subjects)
    tmp_pars <- pars[idx,]
    if (!is.null(design$dynamic)) {
      # Now make dynamic pars
      tmp_p_vector <- add_constants(tmp_p_vector,design$constants)
      dynamic <- attr(tmp_dadm,"dynamic")
      # Get number of q values for each dynamic parameter
      nq <- lapply(dynamic,function(x){dynamic_names("learn")[x$dyntype]})
      nq[is.na(nq)] <- NULL # remove non-dynamic
      if (length(nq)==0) stop("No learning functions!")
      # List to contain q values
      dqlist <- setNames(vector(mode="list",length=length(nq)),names(nq))
      # Collect together bits to do sequential updating
      dpars <- setNames(vector(mode="list",length=length(dqlist)),names(dqlist))
      for (p in names(dpars)) {
        # Copy relevant bits of dynamic
        dpars[[p]]$dyntype <- attr(tmp_dadm, "dynamic")[[p]]$dyntype
        dpars[[p]]$maptype <- attr(tmp_dadm, "dynamic")[[p]]$maptype
        dpars[[p]]$dpnames <- attr(tmp_dadm, "dynamic")[[p]]$dpnames
        dpars[[p]]$covnames <- attr(tmp_dadm, "dynamic")[[p]]$covnames
        dpars[[p]]$lR <- attr(tmp_dadm, "dynamic")[[p]]$lR
        # As in map_p make matrix of base parameters
        pm <- t(matrix(tmp_p_vector[,dimnames(attr(tmp_dadm, "designs")[[dynamic[[p]]$ptype]])[[2]]],
                       dimnames=list(dimnames(attr(tmp_dadm, "designs")[[dynamic[[p]]$ptype]])[[2]],NULL))
        )
        dpars[[p]]$pm <- pm[rep(1, dim(tmp_pars)[1]), , drop = FALSE]
        # Linear design mapping matrix
        dpars[[p]]$dmat <- attr(tmp_dadm, "designs")[[dynamic[[p]]$ptype]][attr(attr(tmp_dadm,
                                                                                     "designs")[[dynamic[[p]]$ptype]], "expand"), , drop = FALSE]
        dpars[[p]]$dmat <- dpars[[p]]$dmat[idx,, drop = F]
        # Do appropriate transform on q0 parameter
        dqlist[[p]] <- switch(dpars[[p]]$dyntype,
                              d1U = pnorm(tmp_p_vector[,dpars[[p]]$dpnames[2]]),
                              d2U = rep(pnorm(tmp_p_vector[,dpars[[p]]$dpnames[2]]),2),
                              d1P = exp(tmp_p_vector[,dpars[[p]]$dpnames[2]]),
                              d2P = rep(exp(tmp_p_vector[,dpars[[p]]$dpnames[2]]),2),
                              d1 =  tmp_p_vector[,dpars[[p]]$dpnames[2]],
                              d2 =  rep(tmp_p_vector[,dpars[[p]]$dpnames[2]],2)
        )
        covnames <- attributes(tmp_dadm)$dynamic[[p]]$covnames
        if (length(covnames)>1 && covnames[2]=="winner") {
          dqlist[[p]] <- setNames(rep(dqlist[[p]],length(covnames)-2),covnames[-c(1:2)])
        }
      }
    }
    if (!is.null(design$adaptive)) {
      adaptive <- design$adaptive
      # Get number of q values for each adaptive parameter
      nq <- lapply(adaptive,function(x){dynamic_names("learn")[x$dyntype]})
      nq[is.na(nq)] <- NULL # remove non-dynamic
      if (length(nq)==0) stop("No learning functions!")
      # List to contain q values
      aqlist <- setNames(vector(mode="list",length=length(nq)),names(nq))
      # Collect together bits to do sequential updating
      apars <- setNames(vector(mode="list",length=length(aqlist)),names(aqlist))
      for (p in names(apars)) {
        # Copy relevant bits of adaptive
        apars[[p]]$dyntype <- attr(tmp_dadm, "adaptive")[[p]]$dyntype
        apars[[p]]$maptype <- attr(tmp_dadm, "adaptive")[[p]]$maptype
        if (apars[[p]]$dyntype=="d1b") apars[[p]]$maptype <- "lin"
        apars[[p]]$apnames <- attr(tmp_dadm, "adaptive")[[p]]$apnames
        apars[[p]]$aptypes <- attr(tmp_dadm, "adaptive")[[p]]$aptypes
        apars[[p]]$covnames <- attr(tmp_dadm, "adaptive")[[p]]$covnames
        apars[[p]]$pcovnames <- attr(tmp_dadm, "adaptive")[[p]]$pcovnames
        apars[[p]]$lR <- attr(tmp_dadm, "adaptive")[[p]]$lR
        # Do appropriate transform on q0 parameter
        aqlist[[p]] <- switch(apars[[p]]$dyntype,
                              d1U = pnorm(tmp_pars[1,apars[[p]]$aptypes[2]]),
                              d2U = rep(pnorm(tmp_pars[1,apars[[p]]$aptypes[2]]),2),
                              d1P = exp(tmp_pars[1,apars[[p]]$aptypes[2]]),
                              d2P = rep(exp(tmp_pars[1,apars[[p]]$aptypes[2]]),2),
                              d1 =  tmp_pars[1,apars[[p]]$aptypes[2]],
                              d2 =  rep(tmp_pars[1,apars[[p]]$aptypes[2]],2),
                              d1b = 0
        )
        covnames <- attributes(tmp_dadm)$adaptive[[p]]$covnames
        if (length(covnames) > 1 && covnames[2]=="winner")
          aqlist[[p]] <- setNames(rep(aqlist[[p]],length(covnames)-2),covnames[-c(1:2)])

      }
    }
    # Make container for simulated data and fill in first line
    data <- tmp_dadm[tmp_dadm$lR==levels(tmp_dadm$lR)[1],]
    acci <- 1:nacc
    tmp_pars[acci,] <- design$model()$Ttransform(design$model()$Ntransform(
      tmp_pars[acci,,drop=FALSE],attr(design, "transform_names")),ndadadm)
    if (!is.null(design$adaptive)) {
      for (p in names(adaptive)) {
        tmp_pars[acci,p] <- switch(apars[[p]]$maptype, # map q to parameter
                                   lin = tmp_pars[acci,p] + tmp_p_vector[,apars[[p]]$aptypes[1]]*aqlist[[p]][1],
                                   plin = exp(tmp_pars[acci,p]) + exp(tmp_p_vector[,apars[[p]]$aptypes[1]])*aqlist[[p]][1],
                                   ucent = tmp_pars[acci,p] + tmp_p_vector[,apars[[p]]$aptypes[1]]*(aqlist[[p]][1]-.5))
      }
    }
    data[1,c("R","rt")] <- design$model()$rfun(tmp_dadm$lR[acci],tmp_pars[acci,,drop=FALSE])
    for (cv in names(dynamic_cv)) {
      if (cv=="winner") {
        value <- t(tmp_dadm[acci,][,attr(tmp_dadm, "isRL")])
      } else if (!is.null(design$adaptive) && !is.null(apars[[p]]$pcovnames)) {
        data[1,cv] <- dynamic_cv[[cv]](data[1,,drop=FALSE],tmp_pars[acci,apars[[p]]$pcovnames])
      } else data[1,cv] <- dynamic_cv[[cv]](data[1,,drop=FALSE])
    }
    # Main loop
    for (i in 2:ntrials) {
      j <- (1+(i-1)*nacc):(i*nacc) # pick out rows of dadm for current trial
      if (!is.null(design$dynamic)) {
        for (p in names(dynamic)) {  # loop over dynamic parameters, updating q and pars
          if (any(dpars[[p]]$covnames=="winner"))  {
            resp <- as.numeric(data[i-1,"R"]) # only update last winning response
            learn <- names(value[,resp])[!is.na(value[,resp])]
            useq <- tmp_dadm[[dpars[[p]]$covnames[1]]][j] # q values to use in computing  current v
            dqlist[[p]][learn] <- one_d(dqlist[[p]][learn],tmp_p_vector[,dpars[[p]]$dpnames][-c(1:2)],
                                        value[!is.na(value)][resp],dpars[[p]]$dyntype)
          } else {
            target <- data[i - 1, dpars[[p]]$covnames]
            dqlist[[p]] <- one_d(dqlist[[p]],tmp_p_vector[,dpars[[p]]$dpnames][-c(1:2)],
                                 target,dpars[[p]]$dyntype)
            useq <- length(dqlist[[p]])
          }
          pm <- dpars[[p]]$pm[j,,drop=FALSE] # Current base parameter
          pm[,p] <- switch(dpars[[p]]$maptype, # map q to parameter
                           lin = pm[,p] + tmp_p_vector[,dpars[[p]]$dpnames[1]]*dqlist[[p]][useq],
                           plin = exp(pm[,p]) + exp(tmp_p_vector[,dpars[[p]]$dpnames[1]])*dqlist[[p]][useq],
                           ucent = pm[,p] + tmp_p_vector[,dpars[[p]]$dpnames[1]]*(dqlist[[p]][useq]-.5)
          )
          tmp_pars[j, dynamic[[p]]$ptype] <- apply(pm*dpars[[p]]$dmat[j,], 1, sum) # Linear design mapping
        }
      }
      if (!is.null(design$adaptive)) {
        for (p in names(adaptive)) {
          if (any(apars[[p]]$covnames=="winner"))  {
            resp <- as.numeric(data[i-1,"R"]) # only update last winning response
            learn <- names(value[,resp])[!is.na(value[,resp])]
            aqlist[[p]][learn] <- one_d(aqlist[[p]][learn],tmp_pars[j[1],apars[[p]]$aptypes[-c(1:2)]],
                                        value[!is.na(value)][resp],apars[[p]]$dyntype)
            useq <- tmp_dadm[[apars[[p]]$covnames[1]]][j] # q values to use in computing  current v
          } else {
            target <- data[i-1,apars[[p]]$covnames]
            aqlist[[p]] <- one_d(aqlist[[p]],tmp_pars[j[1],apars[[p]]$aptypes[-c(1:2)]],
                                 target,apars[[p]]$dyntype)
            useq <- length(aqlist[[p]])
          }
          tmp_pars[j,p] <- switch(apars[[p]]$maptype, # map q to parameter
                                  lin = tmp_pars[j,p] + tmp_pars[j,apars[[p]]$aptypes[1]]*aqlist[[p]][useq],
                                  plin = exp(tmp_pars[,p]) + exp(tmp_pars[j,apars[[p]]$aptypes[1]])*aqlist[[p]][useq],
                                  ucent = tmp_pars[j,p] + tmp_pars[j,apars[[p]]$aptypes[1]]*(aqlist[[p]][useq]-.5))
        }
      }
      # Having done all dynamic and adaptive, so updates simulate data
      tmp_pars[j,] <- design$model()$Ttransform(design$model()$Ntransform(
        tmp_pars[j,,drop=FALSE],attr(design, "transform_names")),ndadadm)
      # if (!all(design$model()$rfun(NULL,tmp_pars[j,,drop=FALSE])))
      #   stop("Parameter generation for ",p," not ok: ",tmp_pars[j,p])
      data[i,c("R","rt")] <- design$model()$rfun(tmp_dadm$lR[j],tmp_pars[j,,drop=FALSE])
      # Finally update cvs
      for (cv in names(dynamic_cv)) {
        if (cv=="winner") value <- t(tmp_dadm[j,][,attr(tmp_dadm, "isRL")]) else {
          if (!is.null(design$adaptive)) {
            if (!is.null(apars[[p]]$pcovnames)) {
              data[i,cv] <- dynamic_cv[[cv]](data[i,,drop=FALSE],tmp_pars[j,apars[[p]]$pcovnames])
            } else if (apars[[p]]$dyntype=="d1b") {
              data[i,cv] <- dynamic_cv[[cv]](data[i,,drop=FALSE],tmp_pars[j,p])
            }
          } else data[i,cv] <- dynamic_cv[[cv]](data[i,,drop=FALSE])
        }
      }
    }
    nams <- c("R","rt",names(dynamic_cv))
    nams <- nams[nams!="winner"]
    # Revert cv back to original form based only on data
    if (!is.null(design$adaptive)) for (p in names(adaptive)) for (cv in names(dynamic_cv)) {
      if (!is.null(apars[[p]]$pcovnames) | (apars[[p]]$dyntype=="d1b"))
        data[,cv] <- dynamic_cv[[cv]](data[,,drop=FALSE])
    }
    data_out <- rbind(data_out, data[,nams])
  }
  return(data_out)
}
