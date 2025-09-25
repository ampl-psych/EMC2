#### Staircase ----
#' Assign Stop-Signal Delays (SSDs) to trials
#'
#' @description
#' This function assigns stop-signal delays (SSDs) to trials based on specified probabilities.
#'
#' @param d A data frame containing trial data with a response factor column 'lR'
#' @param SSD A vector of stop-signal delays to be assigned to trials
#' @param pSSD A vector of probabilities for each SSD value. If length is one less than SSD,
#'             the remaining probability is calculated automatically. Default is 0.25.
#'
#' @return A vector of SSDs with the same length as the number of rows in 'd'.
#'         Trials without a stop signal are assigned Inf.
SSD_function <- function(d,SSD=NA,pSSD=.25) {
  if (sum(pSSD)>1) stop("pSSD sum cannot exceed 1.")
  if (length(pSSD)==length(SSD)-1) pSSD <- c(pSSD,1-sum(pSSD))
  if (length(pSSD)!=length(SSD))
    stop("pSSD must be the same length or 1 less than the length of SSD")
  n_acc <- length(levels(d$lR))
  n_trial <- nrow(d)
  out <- rep(Inf,n_trial)
  trials <- c(1:n_trial)
  for (i in 1:length(pSSD)) {
    pick <- sample(trials,floor(pSSD[i]*n_trial))
    out[pick] <- SSD[i]
    trials <- trials[!(trials %in% pick)]
  }
  return(out)
}

staircase_function <- function(dts,staircase) {
  ns <- ncol(dts)
  SSD <- sR <- srt <- numeric()
  SSD[1] <- staircase$SSD0
  rules <- staircase$rules
  if (is.null(rules)) rules <- list(up = NULL, down = NULL)
  labels <- staircase$labels
  accST <- staircase$accST
  iSSD <- 1
  if (!is.null(accST)) iSSD <- c(iSSD, accST)
  match_rule <- function(label, rule) {
    if (is.null(rule) || !length(rule)) return(FALSE)
    if (is.na(label)) {
      any(is.na(rule))
    } else {
      label %in% rule[!is.na(rule)]
    }
  }
  for (i in 1:ns) {
    if (SSD[i]<staircase$stairmin) SSD[i] <- staircase$stairmin
    if (SSD[i]>staircase$stairmax) SSD[i] <- staircase$stairmax
    trial <- dts[,i]
    trial[iSSD] <- trial[iSSD] + SSD[i]
    if (all(is.infinite(trial[-1]))) {
      Ri <- 1
    } else {
      Ri <- which.min(trial)
    }
    if (Ri==1) {
      if (!is.null(accST) && length(accST) > 0) {
        st_cols <- accST
        st_finish <- trial[st_cols]
        st_idx <- st_cols[which.min(st_finish)]
        sR[i] <- st_idx - 1
        srt[i] <- st_finish[which.min(st_finish)]
        label <- if (!is.null(labels) && (st_idx-1) <= length(labels)) labels[st_idx-1] else NA_character_
      } else {
        sR[i] <- srt[i] <- NA
        label <- NA_character_
      }
    } else {
      sR[i] <- Ri-1
      if (Ri==1) {
        srt[i] <- NA
      } else {
        srt[i] <- trial[Ri]
      }
      label <- if (!is.null(labels) && (Ri-1) <= length(labels)) labels[Ri-1] else NA_character_
    }
    step_dir <- NULL
    if (!is.null(rules$up) || !is.null(rules$down)) {
      success <- match_rule(label, rules$up)
      failure <- match_rule(label, rules$down)
      if (!is.null(rules$down) && !is.null(rules$up) && success && failure) {
        stop("`staircase_up` and `staircase_down` overlap for label ", label)
      }
      if (is.null(rules$down) && !is.null(rules$up)) {
        failure <- !success
      }
      if (success) {
        step_dir <- "up"
      } else if (failure) {
        step_dir <- "down"
      }
    }
    if (is.null(step_dir)) {
      if (Ri==1) step_dir <- "up" else step_dir <- "down"
    }
    if (i<ns) {
      if (identical(step_dir, "up")) {
        SSD[i+1] <- round(SSD[i] + staircase$stairstep,3)
      } else if (identical(step_dir, "down")) {
        SSD[i+1] <- round(SSD[i] - staircase$stairstep,3)
      }
    }
  }
  list(sR=sR,srt=srt,SSD=SSD)
}


apply_staircase_trials <- function(dts, staircase, accST = NULL) {
  if (inherits(staircase, "emc_staircase") && !is.null(staircase$specs)) {
    return(apply_grouped_staircase(dts, staircase, accST))
  }

  stair_fun <- attr(staircase, "staircase_function")
  if (length(accST) > 0) {
    staircase$accST <- 1 + accST
  }
  if (is.null(stair_fun)) {
    stair_fun <- staircase_function
  }
  stair_fun(dts, staircase)
}


apply_grouped_staircase <- function(dts, staircase, accST = NULL) {
  specs <- staircase$specs
  group_id <- staircase$group_id
  data_meta <- staircase$data
  if (is.null(staircase$rules)) staircase$rules <- attr(specs, "rules")
  if (is.null(specs) || is.null(group_id)) {
    stop("Grouped staircase specifications are incomplete.")
  }
  if (length(group_id) != ncol(dts)) {
    stop("Grouped staircase information does not match the number of staircase trials.")
  }

  res <- list(
    sR = rep(NA_real_, ncol(dts)),
    srt = rep(NA_real_, ncol(dts)),
    SSD = rep(NA_real_, ncol(dts))
  )

  base_fun <- attr(staircase, "staircase_function")
  specs_fun <- attr(specs, "staircase_function")

  for (lvl in names(specs)) {
    idx <- which(group_id == lvl)
    if (!length(idx)) next
    spec <- specs[[lvl]]
    spec$rules <- spec$rules %||% staircase$rules
    spec$labels <- spec$labels %||% staircase$labels
    stair_fun <- spec$staircase_function %||% attr(spec, "staircase_function") %||%
      specs_fun %||% base_fun
    if (length(accST) > 0) {
      spec$accST <- 1 + accST
    }
    if (is.null(stair_fun)) {
      stair_fun <- staircase_function
    }
    if (!is.null(data_meta)) {
      spec$data <- data_meta[idx, , drop = FALSE]
    }
    res_group <- stair_fun(dts[, idx, drop = FALSE], spec)
    if (!is.null(res_group$SSD) && length(res_group$SSD)) {
      SSD0 <- spec$SSD0
      if (is.null(SSD0)) {
        SSD0 <- attr(specs, "base_spec")$SSD0
      }
      if (!is.null(SSD0)) {
        res_group$SSD[1] <- SSD0
      }
    }
    if (!is.null(res_group$sR)) res$sR[idx] <- res_group$sR
    if (!is.null(res_group$srt)) res$srt[idx] <- res_group$srt
    if (!is.null(res_group$SSD)) res$SSD[idx] <- res_group$SSD
  }

  res
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}





#### Single exGaussian functions ----


dexGaussian <- function(rt,pars)
  # exGaussian pdf (returns normal or exponential for small tau/sigma)
{
  isexp <- pars[,"sigma"] < 1e-4 # shifted exponential
  rt[isexp] <- dexp(rt[isexp]-pars[isexp,"mu"],1/pars[isexp,"tau"])
  isnorm <- !isexp & pars[,"tau"] < 0.05 * pars[,"sigma"] # normal
  rt[isnorm] <- dnorm(rt[isnorm], mean = pars[isnorm,"mu"],
                      sd = pars[isnorm,"sigma"])
  isexg <- !(isexp | isnorm)
  if (any(isexg)) {
    s2 <- pars[isexg,"sigma"]^2
    z <- rt[isexg] - pars[isexg,"mu"] - (s2/pars[isexg,"tau"])
    rt[isexg] <- exp(
      log(pnorm(z/pars[isexg,"sigma"])) -
        log(pars[isexg,"tau"]) -
        (z + (s2/(2 *  pars[isexg,"tau"])))/pars[isexg,"tau"]
    )
  }
  rt
}

pexGaussian <- function(rt,pars)
  # exGaussian cdf (returns normal or exponential for small tau/sigma)
{
  isexp <- pars[,"sigma"] < 1e-4 # shifted exponential
  rt[isexp] <- pexp(rt[isexp]-pars[isexp,"mu"],1/pars[isexp,"tau"])
  isnorm <- !isexp & pars[,"tau"] < 0.05 * pars[,"sigma"] # normal
  rt[isnorm] <- pnorm(rt[isnorm], mean = pars[isnorm,"mu"],
                      sd = pars[isnorm,"sigma"])
  isexg <- !(isexp | isnorm)
  if (any(isexg)) {
    s2 <- pars[isexg,"sigma"]^2
    z <- rt[isexg] - pars[isexg,"mu"] - (s2/pars[isexg,"tau"])
    rt[isexg] <-
      pnorm((rt[isexg] - pars[isexg,"mu"])/pars[isexg,"sigma"]) -
      exp(log(pnorm(z/pars[isexg,"sigma"])) +
            ((pars[isexg,"mu"] + (s2/pars[isexg,"tau"]))^2 - (pars[isexg,"mu"]^2) -
               2 * rt[isexg] * (s2/pars[isexg,"tau"]))/(2 * s2))
  }
  rt
}

# following adapted from code released with Tanis et al. (2024)
# https://osf.io/u3k5f/
# DynamicCognitivePsychometrics/dmc/models/WALD-SSEXG/dists.R: functions dEXG and pEXG
# for simplicity assuming no upper truncation

dtexGaussian <- function(rt,pars) {
  out <- Fa <- numeric(length(rt))
  a <- pars[,"exg_lb"]
  ok <- (rt>a) & (rt<Inf) & (a<Inf)
  a_real <- (a>-Inf)
  Fa[ok & a_real] <- pexGaussian(a[ok & a_real], pars[ok & a_real,,drop=FALSE])
  normaliser <- 1 - Fa
  out[ok] <- dexGaussian(rt[ok], pars[ok,,drop=FALSE]) / normaliser[ok]
  return(out)
}

ptexGaussian <- function(rt,pars) {
  out <- Fa <- numeric(length(rt))
  a <- pars[,"exg_lb"]
  ok <- (rt>a) & (rt<Inf) & (a<Inf)
  out[rt==Inf] <- 1
  a_real <- (a>-Inf)
  Fa[ok & a_real] <- pexGaussian(a[ok & a_real], pars[ok & a_real,,drop=FALSE])
  normaliser <- 1 - Fa
  out[ok] <- (pexGaussian(rt[ok], pars[ok,,drop=FALSE]) - Fa[ok]) / normaliser[ok]
  return(out)
}


#### Go Single ExGaussian ----

# Go cdf/pdf (strips out NAs and calls d/pexGaussian)
# Is stripping out rt NAs really necessary?

dtexGaussianG <- function(rt,pars)
{
  out <- numeric(length(rt))
  ok <- !is.na(rt)
  out[ok] <- dtexGaussian(rt[ok],pars[ok,,drop=FALSE])
  out
}

ptexGaussianG <- function(rt,pars)
{
  out <- numeric(length(rt))
  ok <- !is.na(rt)
  out[ok] <- ptexGaussian(rt[ok],pars[ok,,drop=FALSE])
  out
}



#### Stop Single ExGaussian ----
dtexGaussianS <- function(rt,pars)
{
  rt <- rt - pars[,"SSD"]
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="muS"] <- "mu"
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="sigmaS"] <- "sigma"
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="tauS"] <- "tau"
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="exgS_lb"] <- "exg_lb"
  dtexGaussian(rt,pars)
}


ptexGaussianS <- function(rt,pars)
{
  rt <- rt - pars[,"SSD"]
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="muS"] <- "mu"
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="sigmaS"] <- "sigma"
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="tauS"] <- "tau"
  dimnames(pars)[[2]][dimnames(pars)[[2]]=="exgS_lb"] <- "exg_lb"
  ptexGaussian(rt,pars)
}


#### ExG Race function ----

# Following functions moved to C++ model_SS_EXG.cpp

# dEXGrace <- function(dt,mu,sigma,tau)
#   # Generates defective PDF for win by first runner, dt (decison time) is
#   # a matrix with length(mu) rows, one row for each runner, and one column
#   # for each decision time for which a defective density value will be
#   # returned.
# {
#   dt[1,] <- dEXG(dt[1,],mu[1],sigma[1],tau[1])
#   if (length(mu)>1) for (i in 2:length(mu))
#     dt[1,] <- dt[1,]*pEXG(dt[i,],mu[i],sigma[i],tau[i],lower_tail=FALSE)
#   dt[1,]
# }

#
#
# stopfn_exg <- function(t,mu,sigma,tau,SSD)
#   # Used by my.integrate, t = vector of times, SSD is a scalar stop-signal delay.
# {
#   dt <- matrix(rep(t+SSD,each=length(mu)),nrow=length(mu))
#   dt[1,] <- dt[1,]-SSD
#   dEXGrace(dt,mu,sigma,tau)
# }

#### ExGaussian random ----

rexG <- function(n,mu,sigma,tau) rnorm(n,mean=mu,sd=sigma) + rexp(n,rate=1/tau)

# Truncated (lower) ex-Gaussian sampler matching likelihood's lower bound handling
rtexG <- function(n, mu, sigma, tau, lb) {
  # Vectorized over parameters; draws from exG truncated at lb
  # n should equal length(mu)==length(sigma)==length(tau)==length(lb)
  out <- numeric(n)
  need <- rep(TRUE, n)
  if (length(mu) != n || length(sigma) != n || length(tau) != n || length(lb) != n)
    stop("rtexG parameter lengths must equal n")
  while (any(need)) {
    k <- sum(need)
    x <- rnorm(k, mean = mu[need], sd = sigma[need]) + rexp(k, rate = 1/tau[need])
    ok <- x > lb[need]
    if (any(ok)) {
      idx <- which(need)[ok]
      out[idx] <- x[ok]
      need[idx] <- FALSE
    }
  }
  out
}



rexGaussian <- function(lR,pars,p_types=c("mu","sigma","tau"),
                        ok=rep(TRUE,dim(pars)[1]))
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars is a matrix of corresponding parameter values named as in p_types
  # pars must be sorted so accumulators and parameter for each trial are in
  # contiguous rows.
  #
  # test
  # pars=cbind(mu=c(.5,.6),sigma=c(.1,.1),tau=c(.2,.2)); lR=factor(c(1))
{
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ",paste(p_types,collapse = " "))
  dt <- matrix(rexG(dim(pars)[1],pars[,"mu"],pars[,"sigma"],pars[,"tau"]),
               nrow=length(levels(lR)))
  R <- apply(dt,2,which.min)
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  rt <- dt[pick]
  R <- factor(levels(lR)[R],levels=levels(lR))
  cbind.data.frame(R=R,rt=rt)
}


#### EXG Stop signal random -----

rSSexGaussian <- function(data,pars,ok=rep(TRUE,dim(pars)[1]))
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars must contain an SSD column and an lI column indicating if an
  # accumulator is triggered by the stop signal (ST = stop-triggered).

  # NB1: Go failures will only apply to accumulators where lI = TRUE
  #      and can still have a stop-triggered response on a go-failure trial.
{
  lR <- data$lR
  nacc <- length(levels(lR))   # Does not include stop runner
  ntrials <- dim(pars)[1]/nacc # Number of trials to simulate
  is1 <- lR==levels(lR)[1]     # First go accumulator
  acc <- 1:nacc                # choice (go and ST) accumulator index
  allSSD <- data$SSD[is1]
  # stop-triggered racers
  isST <- pars[,"lI"]==1              # Boolean for all pars
  accST <- acc[pars[1:nacc,"lI"]==1]  # Index of ST accumulator, from 1st trial

  # Setup race among nacc+1 (includes stop accumulator)
  # Default Inf finishing time so if not changed always looses
  dt <- matrix(Inf,nrow=nacc+1,ncol=ntrials)

  # Go failure trials (rep for each accumulator)
  isgf <- rep(pars[is1,"gf"] > runif(ntrials),each=nacc)
  # Go accumulators that don't fail (removing ST)
  isGO <- !isgf & !isST
  ngo <- sum(isGO)

  # Fill in go accumulators (apply lower bound exg_lb)
  if (any(isGO)) dt[-1,][isGO] <- rtexG(
    ngo,
    mu = pars[isGO, "mu"], sigma = pars[isGO, "sigma"], tau = pars[isGO, "tau"],
    lb = pars[isGO, "exg_lb"]
  )

  # pick out stop trials and races with SSD that is not Inf (i.e., finite or
  # NA, the latter so a staircase can be filled in)
  isStrial <- !is.infinite(pars[is1,"SSD"])
  isSrace <- rep(isStrial,each=nacc)

  # pick out stop trials that are triggered
  isT <- pars[is1,"tf"][isStrial] < runif(sum(isStrial))

  # Logical to pick stop-triggered accumulators that are triggered
  isSTT <- logical(ntrials*nacc) # Default false

  # Pick out stop-triggered accumulators that are triggered
  # NB: can occur on gf trials
  isSTT[isSrace][rep(isT,each=nacc) & isST[isSrace]] <- TRUE
  nst <- sum(isSTT) # Number of triggered ST accumulators

  # Fill in stop-triggered accumulators (apply lower bound exg_lb)
  if (any(isSTT)) dt[-1,][isSTT] <- rtexG(
    nst,
    mu = pars[isSTT, "mu"], sigma = pars[isSTT, "sigma"], tau = pars[isSTT, "tau"],
    lb = pars[isSTT, "exg_lb"]
  )

  # pick out triggered stop racers
  isTS <- logical(ntrials)
  isTS[isStrial][isT] <- TRUE
  ns <- sum(isTS)

  # Fill in stop accumulators (apply lower bound exgS_lb)
  if (any(isTS)) dt[1, isTS] <- rtexG(
    ns,
    mu = pars[is1, "muS"][isTS], sigma = pars[is1, "sigmaS"][isTS], tau = pars[is1, "tauS"][isTS],
    lb = pars[is1, "exgS_lb"][isTS]
  )

  # staircase algorithm
  pstair <- is.na(pars[,"SSD"])
  stair <- pstair[is1]
  if (any(stair)) {
    if (is.null(attr(data,"staircase")))
      stop("When SSD has NAs a staircase list must be supplied!")
    staircase <- attr(data,"staircase")

    allR <- allrt <- allSSD <- numeric(ncol(dt))  # to store unified results
    allSSD[] <- Inf
    dts <- dt[,stair,drop=F]

    # Non-staircase trials
    dt <- dt[,!stair,drop=F]
    spars <- pars[pstair,,drop=F]
    pars <- pars[!pstair,,drop=F]
    isTS <- isTS[!stair]
    isSTT <- isSTT[!pstair]
    isST <- isST[!pstair]
    ntrials <- sum(!stair)
    is1 <- is1[!pstair]
    stair_res <- apply_staircase_trials(dts, staircase, accST)
    allR[stair] <- stair_res$sR
    allrt[stair] <- stair_res$srt
    if (inherits(staircase, "emc_staircase") && !is.null(staircase$specs)) {
      stair_idx <- which(stair)
      gid <- as.character(staircase$group_id)
      base_spec <- attr(staircase$specs, "base_spec")
      for (lvl in names(staircase$specs)) {
        cols <- which(gid == lvl)
        if (!length(cols)) next
        pos <- stair_idx[cols[1]]
        spec <- staircase$specs[[lvl]]
        ssd0 <- spec$SSD0
        if (is.null(ssd0) && !is.null(base_spec)) ssd0 <- base_spec$SSD0
        if (!is.null(ssd0)) {
          stair_res$SSD[cols[1]] <- ssd0
          allSSD[pos] <- ssd0
        }
      }
    }
    allSSD[stair] <- stair_res$SSD
  }


  # All SSD already filled in so return R and rt

  # Add SSD to triggered stop accumulators
  if (any(isTS)) dt[1,isTS] <- dt[1,isTS] + pars[is1,"SSD"][isTS]

  # Add SSD to stop-triggered accumulators that are triggered
  if (any(isSTT)) dt[-1,][isSTT] <- dt[-1,][isSTT] + pars[isSTT,"SSD"]

  # R <- factor(rep(NA,ntrials),levels=levels(lR))
  R <- rt <- rep(NA,ntrials)

  # All accumulators Inf (usually when both gf and tf)
  allinf <- apply(dt,2,function(x)all(is.infinite(x)))

  # get winner of stop and go where there is a race
  r <- c(0, acc)[apply(dt[,!allinf,drop=FALSE],2,which.min)]

  # stop wins
  stopwins <- r==0

  # First fill in cases where stop looses, can include wins by ST
  if (any(!stopwins)) {
    rgo <- r[!stopwins]
    R[!allinf][!stopwins] <- rgo
    pick <- cbind(rgo,c(1:sum(!stopwins))) # Matrix to pick winner
    rt[!allinf][!stopwins] <-
      dt[-1,!allinf,drop=FALSE][,!stopwins,drop=FALSE][pick]
  }

  # then if stop wins and kills go, find ST winner
  if (any(isST) & any(stopwins)) {
    # stop triggered accumulators that are racing
    rst <- dt[-1,!allinf,drop=FALSE][accST,stopwins,drop=FALSE]
    # stop-triggered winners
    rtw <- apply(rst,2,which.min)
    # index for stop-triggered
    R[!allinf][stopwins] <- accST[rtw]
    pick <- cbind(rtw,1:ncol(rst))
    rt[!allinf][stopwins] <- rst[pick]
  }

  rt[is.na(R)] <- NA

  if (any(stair)) {
    allrt[!stair] <- rt
    allR[!stair] <- R
    allSSD[stair] <- stair_res$SSD
    out <- cbind.data.frame(R=factor(allR,levels=1:nacc,labels=levels(lR)),
                            rt=allrt, SSD = allSSD)
    return(out)
  }
  cbind.data.frame(R=factor(R,levels=1:nacc,labels=levels(lR)),rt=rt)
}


#### ExG stop probability (no stop triggered) ----

my.integrate <- function(..., upper = Inf, big = 10) {
  out <- try(
    integrate(..., upper = upper, rel.tol = 1e-6, abs.tol = 1e-8),
    silent = TRUE
  )
  if (inherits(out, "try-error")) {
    return(0)
  }
  if (upper == Inf && out$subdivisions == 1) {
    out <- try(
      integrate(..., upper = big, rel.tol = 1e-6, abs.tol = 1e-8),
      silent = TRUE
    )
    if (inherits(out, "try-error") || out$subdivisions == 1) {
      return(0)
    }
  }
  return(out$value)
}

pstopTEXG <- function(
    parstop, n_acc, upper=Inf,
    gpars=c("mu","sigma","tau","exg_lb"), spars=c("muS","sigmaS","tauS","exgS_lb")
) {
  sindex <- seq(1,nrow(parstop),by=n_acc)  # Stop accumulator index
  ps <- parstop[sindex,spars,drop=FALSE]   # Stop accumulator parameters
  SSDs <- parstop[sindex,"SSD",drop=FALSE] # SSDs
  ntrials <- length(SSDs)
  if (length(upper)==1) upper <- rep(upper,length.out=ntrials)
  pgo <- array(parstop[,gpars],dim=c(n_acc,ntrials,length(gpars)),
               dimnames=list(NULL,NULL,gpars))
  cells <- apply(
    cbind(SSDs,ps,upper,matrix(as.vector(aperm(pgo,c(2,1,3))),nrow=ntrials))
    ,1,paste,collapse="")
  # cells <- character(ntrials)
  # for (i in 1:ntrials)
  #   cells[i] <- paste(SSDs[i],ps[i,],pgo[,i,],upper[i],collapse="")
  uniq <- !duplicated(cells)
  ups <- sapply(which(uniq),function(i){
    my.integrate(f=stopfn_texg,lower=ps[i,"exgS_lb"],SSD=SSDs[i],upper=upper[i],
                 mu=c(ps[i,"muS"],pgo[,i,"mu"]),
                 sigma=c(ps[i,"sigmaS"],pgo[,i,"sigma"]),
                 tau=c(ps[i,"tauS"],pgo[,i,"tau"]),
                 lb=c(ps[i,"exgS_lb"],pgo[,i,"exg_lb"]))
  })
  ups[as.numeric(factor(cells,levels=cells[uniq]))]
}


#' The ex-Gaussian race model of the stop signal task
#'
#' Model file to estimate the ex-Gaussian race model of stop signal task data in EMC2.
#'
#' Model files are almost exclusively used in `design()`.
#'
#' @details
#'
#' Default values are used for all parameters that are not explicitly listed in the `formula`
#' argument of `design()`.They can also be accessed with `SSexG()$p_types`.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default**   | **Mapping**                    | **Interpretation**                                            |
#' |-----------|-----------|---------------|-----------|----------------------------|-----------------------------------------------------------|
#' | *mu*       | log         | \[0, Inf\]     | log(.4)         |                            | Mean of Gaussian component of ex-Gaussian go finish time distribution              |
#' | *sigma*       | log       | \[0, Inf\]        | log(.05)    |                            | Standard deviation of Gaussian component of ex-Gaussian go finish time distribution                                        |
#' | *tau*      | log       | \[0, Inf\]        | log(.1)    |                            | Mean (inverse rate) of exponential component of ex-Gaussian go finish time distribution                                          |
#' | *muS*       | log       | \[0, Inf\]        | log(.3)    |                            | Mean of Gaussian component of ex-Gaussian stop finish time distribution           |
#' | *sigmaS*       | log    | \[0, Inf\]        | log(.025)|                   | Standard deviation of Gaussian component of ex-Gaussian stop finish time distribution                              |
#' | *tauS*      | log    | \[0, Inf\]        | log(.05)  |  | Mean (inverse rate) of exponential component of ex-Gaussian stop finish time distribution       |
#' | *tf*      | probit       | \[0, 1\]        | qnorm(0)    |                            | Attentional lapse rate for stop process ("trigger failure")           |
#' | *gf*     | probit       | \[0, 1\]        | qnorm(0)    |                            | Attentional lapse rate for go process ("go failure")    |
#' | *exg_lb*      | -       | \[-Inf, Inf\]        | .05    |                            | Lower bound of ex-Gaussian go finish time distribution           |
#' | *exgS_lb*     | -       | \[-Inf, Inf\]        | .05    |                            | Lower bound of ex-Gaussian stop finish time distribution    |
#'
#' Because the ex-Gaussian stop signal model is a race model, it has one accumulator per response option.
#' EMC2 automatically constructs a factor representing the accumulators `lR` (i.e., the
#' latent response) with level names taken from the `R` column in the data.
#'
#' For race models, the `design()` argument `matchfun` can be provided, a
#' function that takes the `lR` factor (defined in the augmented data (d)
#' in the following function) and returns a logical defining the correct response.
#' In the example below, the match is simply such that the `S` factor equals the
#' latent response factor: `matchfun=function(d)d$S==d$lR`. Then `matchfun` is
#' used to automatically create a latent match (`lM`) factor with
#' levels `FALSE` (i.e., the stimulus does not match the accumulator) and `TRUE`
#' (i.e., the stimulus does match the accumulator). This is added internally
#' and can also be used in model formula.
#'
#' The race model of the stop signal task was first introduced by Logan & Cowan (1984), who modeled task performance as a race between two stochastic processes: a go process for response execution (triggered by the choice stimulus), and a stop process for response inhibition (triggered by the stop signal).
#' Several review papers provide further background on the task paradigm and modelling framework (e.g., Matzke et al., 2018; Verbruggen et al., 2019; Colonius et al., 2023).
#'
#' Early work on the ex-Gaussian distribution as a descriptive model of RT data includes Ratcliff & Murdock (1976), Ratcliff (1979), and Heathcote et al. (1991).
#' Matzke et al. (2013) introduced the ex-Gaussian distribution as a parametric model of the finish time distributions of the go and stop processes.
#' This work was extended in Matzke et al. (2019) to allow for an arbitrary number of go processes.
#' Lastly, attentional lapse rate parameters for the stop process and go process, termed "trigger failure" (`tf`) and "go failure" (`gf`) respectively, were introduced in Matzke et al. (2017) and Matzke et al. (2019).
#'
#' Note that the ex-Gaussian parameters `mu`, `sigma`, and `tau` do not have clear psychological interpretations (Matzke & Wagenmakers, 2009).
#' Inference is typically based on the mean of the ex-Gaussian distribution, which is given by `mu + tau`.
#' The mean of the ex-Gaussian stop process finish time distribution (`muS + tauS`) is taken as the stop signal reaction time (SSRT), which is typically the primary modelling outcome of interest.
#'
#' The ex-Gaussian distribution has support on the real line \eqn{\left(-\infty, \infty\right)}.
#' To prevent evaluation of impossible (i.e., negative) or implausibly fast finish times, lower truncation is applied to both the go and stop finish time distributions, using the parameters `exg_lb` and `exgS_lb`, respectively.
#' The default values for these lower bounds are `.05`, based on empirical estimates of the onset latency of early sensory processing (Schmolesky et al., 1998) and in line with Tanis et al. (2024).
#' If strict replication of the "original" ex-Gaussian race models (e.g., Matzke et al., 2013; 2019) is desired, both `exg_lb` and `exgS_lb` should be set to `-Inf`, using the `constants` argument of the [design()] function.
#'
#' @references
#' Colonius, H., & Diederich, A. (2023). Modeling response inhibition in the stop signal task. *F. G. Ashby, H. Colonius, & E. N. Dzhafarov (Eds.), New Handbook of Mathematical Psychology*, *3*, 311-356. \doi{10.1017/9781108902724.008}
#'
#' Heathcote, A., Popiel, S. J., & Mewhort, D. J. (1991). Analysis of response time distributions: an example using the Stroop task. *Psychological Bulletin*, *109*(2), 340. \doi{10.1037/0033-2909.109.2.340}
#'
#' Logan, G. D., & Cowan, W. B. (1984). On the ability to inhibit thought and action: A theory of an act of control. *Psychological Review*, *91*(3), 295. \doi{10.1037/0033-295X.91.3.295}
#'
#' Matzke, D., & Wagenmakers, E. J. (2009). Psychological interpretation of the ex-Gaussian and shifted Wald parameters: A diffusion model analysis. *Psychonomic Bulletin & Review*, *16*, 798-817. \doi{10.3758/PBR.16.5.798}
#'
#' Matzke, D., Dolan, C. V., Logan, G. D., Brown, S. D., & Wagenmakers, E. J. (2013). Bayesian parametric estimation of stop-signal reaction time distributions. *Journal of Experimental Psychology: General*, *142*(4), 1047. \doi{10.1037/a0030543}
#'
#' Matzke, D., Love, J., & Heathcote, A. (2017). A Bayesian approach for estimating the probability of trigger failures in the stop-signal paradigm. *Behavior Research Methods*, *49*, 267-281. \doi{10.3758/s13428-015-0695-8}
#'
#' Matzke, D., Verbruggen, F., & Logan, G. (2018). The stop-signal paradigm. *Stevens’ handbook of experimental psychology and cognitive neuroscience*, *5*, 383-427. \doi{10.1002/9781119170174.epcn510}
#'
#' Matzke, D., Curley, S., Gong, C. Q., & Heathcote, A. (2019). Inhibiting responses to difficult choices. *Journal of Experimental Psychology: General*, *148*(1), 124. \doi{10.1037/xge0000525}
#'
#' Ratcliff, R., & Murdock, B. B. (1976). Retrieval processes in recognition memory. *Psychological Review*, *83*(3), 190. \doi{10.1037/0033-295X.83.3.190}
#'
#' Ratcliff, R. (1979). Group reaction time distributions and an analysis of distribution statistics. *Psychological Bulletin*, *86*(3), 446. \doi{10.1037/0033-2909.86.3.446}
#'
#' Schmolesky, M. T., Wang, Y., Hanes, D. P., Thompson, K. G., Leutgeb, S., Schall, J. D., & Leventhal, A. G. (1998). Signal timing across the macaque visual system. *Journal of Neurophysiology*, *79*(6), 3272-3278. \doi{10.1152/jn.1998.79.6.3272}
#'
#' Tanis, C. C., Heathcote, A., Zrubka, M., & Matzke, D. (2024). A hybrid approach to dynamic cognitive psychometrics: Dynamic cognitive psychometrics. *Behavior Research Methods*, *56*(6), 5647-5666. \doi{10.3758/s13428-023-02295-y}
#'
#' Verbruggen, F., Aron, A. R., Band, G. P., Beste, C., Bissett, P. G., Brockett, A. T., ... & Boehler, C. N. (2019). A consensus guide to capturing the ability to inhibit actions and impulsive behaviors in the stop-signal task. *elife*, *8*, e46323. \doi{10.7554/eLife.46323}
#'
#' @return A model list with all the necessary functions to sample
#' @export
SSEXG <- function() {
  list(
    type = "RACE",
    c_name = "SSEXG",
    p_types = c(
      mu = log(.4), sigma = log(.05), tau = log(.1),
      muS = log(.3), sigmaS = log(.025), tauS = log(.05),
      tf = qnorm(0), gf = qnorm(0),
      exg_lb = .05, exgS_lb = .05
    ),
    transform = list(
      func = c(
        mu = "exp", sigma = "exp", tau = "exp",
        muS = "exp", sigmaS = "exp", tauS = "exp",
        tf = "pnorm", gf = "pnorm",
        exg_lb = "identity", exgS_lb = "identity"
      )
    ),
    bound = list(
      minmax = cbind(
        mu = c(0, Inf), sigma = c(1e-4, Inf), tau = c(1e-4, Inf),
        muS = c(0, Inf), sigmaS = c(1e-4, Inf), tauS = c(1e-4, Inf),
        tf = c(.001, .999), gf = c(.001, .999),
        exg_lb = c(-Inf, Inf), exgS_lb = c(-Inf, Inf)
      ),
      exception = c(
        tf = 0, gf = 0,
        exg_lb = -Inf, exgS_lb = -Inf
      )
    ),
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      # if (any(names(dadm)=="SSD")) pars <- cbind(pars,SSD=dadm$SSD) else
      #   pars <- cbind(pars,SSD=rep(Inf,dim(pars)[1]))
      pars <- cbind(pars, SSD = dadm$SSD)
      pars <- cbind(pars, lI = as.numeric(dadm$lI))  # Only necessary for data generation.
      return(pars)
    },
    # Density function (PDF) for single go racer
    dfunG = function(rt, pars) return(dtexGaussianG(rt, pars)),
    # Probability function (CDF) for single go racer
    pfunG = function(rt, pars) return(ptexGaussianG(rt, pars)),
    # Density function (PDF) for single stop racer
    dfunS = function(rt, pars) {
      parsS <- pars[ , c("muS", "sigmaS", "tauS", "SSD", "exgS_lb"), drop=FALSE]
      return(dtexGaussianS(rt, parsS))
    },
    # Probability function (CDF) for single stop racer
    pfunS = function(rt, pars) {
      parsS <- pars[ , c("muS", "sigmaS", "tauS", "SSD", "exgS_lb"), drop=FALSE]
      return(ptexGaussianS(rt, parsS))
    },
    # Stop probability integral
    sfun = function(pars, n_acc, upper = Inf) {
      return(pstopTEXG(pars, n_acc, upper = upper))
    },
    # Random function for SS race
    # TODO
    rfun = function(data = NULL, pars) {
      return(rSSexGaussian(data, pars, ok = attr(pars, "ok")))
    },
    # Race likelihood combining pfun and dfun
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      return(log_likelihood_race_ss(pars, dadm, model, min_ll = min_ll))
    }
  )
}


#####################  RDEX ----

#### RDEX SS random ----
rSShybrid <- function(data,pars,ok=rep(TRUE,dim(pars)[1]))
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars must contain an SSD column and an lI column indicating if an
  # accumulator is triggered by the stop signal (ST = stop-triggered).

  # NB1: Go failures will only apply to accumulators where lI = TRUE
  #      and can still have a stop-triggered response on a go-failure trial.
{
  lR <- data$lR # For Michelle as an example
  pars[,c("A","B","v")] <- pars[,c("A","B","v")]/pars[ok,"s"]

  nacc <- length(levels(lR))   # Does not include stop runner
  ntrials <- dim(pars)[1]/nacc # Number of trials to simulate
  is1 <- lR==levels(lR)[1]     # First go accumulator
  acc <- 1:nacc                # choice (go and ST) accumulator index
  allSSD <- data$SSD[is1]

  # stop-triggered racers
  isST <- pars[,"lI"]==1              # Boolean for all pars
  accST <- acc[pars[1:nacc,"lI"]==1]  # Index of ST accumulator, from 1st trial

  # Setup race among nacc+1 (includes stop accumulator)
  # Default Inf finishing time so if not changed always looses
  dt <- matrix(Inf,nrow=nacc+1,ncol=ntrials)

  # Go failure trials (rep for each accumulator)
  isgf <- rep(pars[is1,"gf"] > runif(ntrials),each=nacc)
  # Go accumulators that don't fail (removing ST)
  isGO <- !isgf & !isST
  ngo <- sum(isGO)

  # Fill in go accumulators
  if (any(isGO)) dt[-1,][isGO] <-
    pars[isGO,"t0"] + rWald(ngo,pars[isGO,"B"],pars[isGO,"v"],pars[isGO,"A"])

  # pick out stop trials and races with SSD that is not Inf (i.e., finite or
  # NA, the latter so a staircase can be filled in)
  isStrial <- !is.infinite(pars[is1,"SSD"])
  isSrace <- rep(isStrial,each=nacc)

  # pick out stop trials that are triggered
  isT <- pars[is1,"tf"][isStrial] < runif(sum(isStrial))

  # Logical to pick stop-triggered accumulators that are triggered
  isSTT <- logical(ntrials*nacc) # Default false

  # Pick out stop-triggered accumulators that are triggered
  # NB: can occur on gf trials
  isSTT[isSrace][rep(isT,each=nacc) & isST[isSrace]] <- TRUE
  nst <- sum(isSTT)

  # Fill in stop-triggered accumulators (ST; Wald go process)
  if (any(isSTT)) dt[-1,][isSTT] <-
    pars[isSTT, "t0"] + rWald(nst, pars[isSTT, "B"], pars[isSTT, "v"], pars[isSTT, "A"])

  # pick out triggered stop racers
  isTS <- logical(ntrials)
  isTS[isStrial][isT] <- TRUE
  ns <- sum(isTS)

  # Fill in stop accumulators (apply lower bound exgS_lb)
  if (any(isTS)) dt[1, isTS] <- rtexG(
    ns,
    mu = pars[is1, "muS"][isTS], sigma = pars[is1, "sigmaS"][isTS], tau = pars[is1, "tauS"][isTS],
    lb = pars[is1, "exgS_lb"][isTS]
  )

  # staircase algorithm
  pstair <- is.na(pars[,"SSD"])
  stair <- pstair[is1]
  if (any(stair)) {
    if (is.null(attr(pars,"staircase")))
      stop("When SSD has NAs a staircase list must be supplied!")

    staircase <- attr(pars,"staircase")

    allR <- allrt <- numeric(ncol(dt))  # to store unified results
    dts <- dt[,stair,drop=F]

    # Non-staircase trials
    dt <- dt[,!stair,drop=F]
    spars <- pars[pstair,,drop=F]
    pars <- pars[!pstair,,drop=F]
    isTS <- isTS[!stair]
    isSTT <- isSTT[!pstair]
    isST <- isST[!pstair]
    ntrials <- sum(!stair)
    is1 <- is1[!pstair]
    stair_res <- apply_staircase_trials(dts, staircase, accST)
    if (inherits(staircase, "emc_staircase") && !is.null(staircase$specs)) {
      stair_idx <- which(stair)
      gid <- as.character(staircase$group_id)
      base_spec <- attr(staircase$specs, "base_spec")
      for (lvl in names(staircase$specs)) {
        cols <- which(gid == lvl)
        if (!length(cols)) next
        pos <- stair_idx[cols[1]]
        spec <- staircase$specs[[lvl]]
        ssd0 <- spec$SSD0
        if (is.null(ssd0) && !is.null(base_spec)) ssd0 <- base_spec$SSD0
        if (!is.null(ssd0)) {
          stair_res$SSD[cols[1]] <- ssd0
          allSSD[pos] <- ssd0
        }
      }
    }
    allR[stair] <- stair_res$sR
    allrt[stair] <- stair_res$srt
    allSSD[stair] <- stair_res$SSD
  }


  # All SSD already filled in so return R and rt

  # Add SSD to triggered stop accumulators
  if (any(isTS)) dt[1,isTS] <- dt[1,isTS] + pars[is1,"SSD"][isTS]

  # Add SSD to stop-triggered accumulators that are triggered
  if (any(isSTT)) dt[-1,][isSTT] <- dt[-1,][isSTT] + pars[isSTT,"SSD"]

  # R <- factor(rep(NA,ntrials),levels=levels(lR))
  R <- rt <- rep(NA,ntrials)

  # All accumulators Inf (usually when both gf and tf)
  allinf <- apply(dt,2,function(x)all(is.infinite(x)))

  # get winner of stop and go where there is a race
  r <- c(0, acc)[apply(dt[,!allinf,drop=FALSE],2,which.min)]

  # stop wins
  stopwins <- r==0

  # First fill in cases where stop looses, can include wins by ST
  if (any(!stopwins)) {
    rgo <- r[!stopwins]
    R[!allinf][!stopwins] <- rgo
    pick <- cbind(rgo,c(1:sum(!stopwins))) # Matrix to pick winner
    rt[!allinf][!stopwins] <-
      dt[-1,!allinf,drop=FALSE][,!stopwins,drop=FALSE][pick]
  }

  # then if stop wins and kills go, find ST winner
  if (any(isST) & any(stopwins)) {
    # stop triggered accumulators that are racing
    rst <- dt[-1,!allinf,drop=FALSE][accST,stopwins,drop=FALSE]
    # stop-triggered winners
    rtw <- apply(rst,2,which.min)
    # index for stop-triggered
    R[!allinf][stopwins] <- accST[rtw]
    pick <- cbind(rtw,1:ncol(rst))
    rt[!allinf][stopwins] <- rst[pick]
  }

  rt[is.na(R)] <- NA

  if (any(stair)) {
    allrt[!stair] <- rt
    allR[!stair] <- R
    out <- cbind.data.frame(R=factor(allR,levels=1:nacc,labels=levels(lR)),
                            rt=allrt, SSD = allSSD)
    return(out)
  }
  cbind.data.frame(R=factor(R,levels=1:nacc,labels=levels(lR)),rt=rt)
}

#### RDEX stop probability ----

# # NB: these functions are in Rcpp
# dWald_RDEX
# pWald_RDEX
# stopfn_rdex


pstopHybrid <- function(
    parstop, n_acc, upper = Inf,
    gpars = c("v", "B", "A", "t0", "s"), spars = c("muS", "sigmaS", "tauS", "exgS_lb")
) {
  sindex <- seq(1,nrow(parstop),by=n_acc)  # Stop accumulator index
  ps <- parstop[sindex,spars,drop=FALSE]   # Stop accumulator parameters
  SSDs <- parstop[sindex,"SSD",drop=FALSE] # SSDs
  ntrials <- length(SSDs)
  if (length(upper)==1) upper <- rep(upper,length.out=ntrials)
  pgo <- array(parstop[,gpars],dim=c(n_acc,ntrials,length(gpars)),
               dimnames=list(NULL,NULL,gpars))
  cells <- apply(
    cbind(SSDs,ps,upper,matrix(as.vector(aperm(pgo,c(2,1,3))),nrow=ntrials))
    ,1,paste,collapse="")
  uniq <- !duplicated(cells)
  ups <- sapply(which(uniq),function(i){
    my.integrate(
      # args passed to `my.integrate`
      f = stopfn_rdex, lower = ps[i, "exgS_lb"], upper = upper[i],
      # args passed to `stopfn_rdex`
      n_acc = n_acc,
      mu = ps[i, "muS"], sigma = ps[i, "sigmaS"], tau = ps[i, "tauS"], lb = ps[i, "exgS_lb"],
      v = pgo[ , i, "v"], B = pgo[ , i, "B"], A = pgo[ , i, "A"], t0 = pgo[ , i, "t0"], s = pgo[ , i, "s"],
      SSD = SSDs[i]
    )
  })
  ups[as.numeric(factor(cells,levels=cells[uniq]))]
}


#### RDEX model list ----

#' The hybrid Wald / ex-Gaussian race model of the stop signal task
#'
#' Model file to estimate the hybrid Wald / ex-Gaussian race model of stop signal task data in EMC2.
#'
#' Model files are almost exclusively used in `design()`.
#'
#' @details
#'
#' Default values are used for all parameters that are not explicitly listed in the `formula`
#' argument of `design()`.They can also be accessed with `SShybrid()$p_types`.
#'
#' | **Parameter** | **Transform** | **Natural scale** | **Default**   | **Mapping**                    | **Interpretation**                                            |
#' |-----------|-----------|---------------|-----------|----------------------------|-----------------------------------------------------------|
#' | *v*       | log       | \[0, Inf\]      | log(1)    |                  | Evidence-accumulation rate (drift rate) of the go process                       |
#' | *A*       | log       | \[0, Inf\]      | log(0)    |                  | Between-trial variation (range) in start point of the go process                |
#' | *B*       | log       | \[0, Inf\]      | log(1)    | *b* = *B* + *A*      | Distance from *A* to *b* (response threshold) of the go process                  |
#' | *t0*      | log       | \[0, Inf\]      | log(0)    |                  | Non-decision time of the go process                                            |
#' | *s*       | log       | \[0, Inf\]      | log(1)    |                  | Within-trial standard deviation of drift rate of the go process                |
#' | *muS*       | log       | \[0, Inf\]        | log(.3)    |                            | Mean of Gaussian component of ex-Gaussian stop finish time distribution           |
#' | *sigmaS*       | log    | \[0, Inf\]        | log(.025)|                   | Standard deviation of Gaussian component of ex-Gaussian stop finish time distribution                              |
#' | *tauS*      | log    | \[0, Inf\]        | log(.05)  |  | Mean (inverse rate) of exponential component of ex-Gaussian stop finish time distribution       |
#' | *tf*      | probit       | \[0, 1\]        | qnorm(0)    |                            | Attentional lapse rate for stop process ("trigger failure")           |
#' | *gf*     | probit       | \[0, 1\]        | qnorm(0)    |                            | Attentional lapse rate for go process ("go failure")    |
#' | *exgS_lb*     | -       | \[-Inf, Inf\]        | .05    |                            | Lower bound of ex-Gaussian stop finish time distribution    |
#'
#' All parameters are estimated on the log scale, with the exception of the attentional failure parameters (`tf` and `gf`), which are estimated on the probit scale.
#'
#' The parameterization *b* = *B* + *A* ensures that the response threshold is
#' always higher than the between trial variation in start point.
#'
#' Conventionally, `s` is fixed to 1 to satisfy scaling constraints.
#'
#' Because the ex-Gaussian stop signal model is a race model, it has one accumulator per response option.
#' EMC2 automatically constructs a factor representing the accumulators `lR` (i.e., the
#' latent response) with level names taken from the `R` column in the data.
#'
#' For race models, the `design()` argument `matchfun` can be provided, a
#' function that takes the `lR` factor (defined in the augmented data (d)
#' in the following function) and returns a logical defining the correct response.
#' In the example below, the match is simply such that the `S` factor equals the
#' latent response factor: `matchfun=function(d)d$S==d$lR`. Then `matchfun` is
#' used to automatically create a latent match (`lM`) factor with
#' levels `FALSE` (i.e., the stimulus does not match the accumulator) and `TRUE`
#' (i.e., the stimulus does match the accumulator). This is added internally
#' and can also be used in model formula.
#'
#' This hybrid race model of the stop signal task was introduced in Tanis et al. (2024).
#'
#' Note that, in contrast to the parameters of the go process, the ex-Gaussian stop process parameters `muS`, `sigmaS`, and `tauS` do not have clear psychological interpretations (Matzke & Wagenmakers, 2009).
#' Matzke et al. (2020) showed that evidence accumulation models of the stop process have poor psychometric properties.
#' Thus, this hybrid race model - with an evidence accumulation account of _going_ and a descriptive account of _stopping_ - represents a compromise between process realism and practically useful measurement properties (Tanis et al., 2024).
#'
#' The mean of the ex-Gaussian stop distribution (`muS + tauS`) is taken as the stop signal reaction time (SSRT).
#'
#' The ex-Gaussian distribution has support on the real line \eqn{\left(-\infty, \infty\right)}.
#' To prevent evaluation of impossible (i.e., negative) or implausibly fast stop finish times, lower truncation is applied to the stop finish time distribution using the parameter `exgS_lb`.
#' The default value for this lower bound is `.05`, based on empirical estimates of the onset latency of early sensory processing (Schmolesky et al., 1998) and in line with Tanis et al. (2024).
#'
#' @references
#'
#' Matzke, D., & Wagenmakers, E. J. (2009). Psychological interpretation of the ex-Gaussian and shifted Wald parameters: A diffusion model analysis. *Psychonomic Bulletin & Review*, *16*, 798-817. \doi{10.3758/PBR.16.5.798}
#'
#' Matzke, D., Logan, G. D., & Heathcote, A. (2020). A cautionary note on evidence-accumulation models of response inhibition in the stop-signal paradigm. *Computational Brain & Behavior*, *3*(3), 269-288. \doi{10.1007/s42113-020-00075-x}
#'
#' Tanis, C. C., Heathcote, A., Zrubka, M., & Matzke, D. (2024). A hybrid approach to dynamic cognitive psychometrics: Dynamic cognitive psychometrics. *Behavior Research Methods*, *56*(6), 5647-5666. \doi{10.3758/s13428-023-02295-y}
#'
#' @return A model list with all the necessary functions to sample
#' @export
SSRDEX <- function() {
  list(
    type = "RACE",
    c_name = "SSRDEX",
    p_types = c(
      v = log(1), B = log(1), A = log(0), t0 = log(0), s = log(1),
      muS = log(.3), sigmaS = log(.025), tauS = log(.05),
      tf = qnorm(0), gf = qnorm(0),
      exgS_lb = .05
    ),
    transform = list(
      func = c(
        v = "exp", B = "exp", A = "exp", t0 = "exp", s = "exp",
        muS = "exp", sigmaS = "exp", tauS = "exp",
        tf = "pnorm", gf = "pnorm",
        exgS_lb = "identity"
      )
    ),
    bound = list(
      minmax = cbind(
        v = c(1e-3, Inf), B = c(0, Inf), A = c(1e-4, Inf), t0 = c(0.05, Inf), s = c(0, Inf),
        muS = c(0, Inf), sigmaS = c(1e-4, Inf), tauS = c(1e-4,Inf),
        tf = c(.001, .999), gf = c(.001, .999),
        exgS_lb = c(-Inf, Inf)
      ),
      exception = c(
        v = 0, A = 0,
        tf = 0, gf = 0,
        exgS_lb = -Inf
      )
    ),
    # Trial dependent parameter transform
    Ttransform = function(pars, dadm) {
      pars <- cbind(pars, b = pars[,"B"] + pars[,"A"])
      pars <- cbind(pars, SSD = dadm$SSD)
      pars <- cbind(pars, lI = as.numeric(dadm$lI))  # Only necessary for data generation.
      return(pars)
    },
    # Density function (PDF) for single go racer
    dfunG = function(rt, pars) {
      return(dRDM(rt, pars))
    },
    # Probability function (CDF) for single go racer
    pfunG = function(rt, pars) {
      return(pRDM(rt, pars))
    },
    # Density function (PDF) for single stop racer
    dfunS = function(rt, pars) {
      parsS <- pars[ , c("muS", "sigmaS", "tauS", "SSD", "exgS_lb"), drop=FALSE]
      return(dtexGaussianS(rt, parsS))
    },
    # Probability function (CDF) for single stop racer
    pfunS = function(rt, pars) {
      parsS <- pars[ , c("muS", "sigmaS", "tauS", "SSD", "exgS_lb"), drop=FALSE]
      return(ptexGaussianS(rt, parsS))
    },
    # Stop probability integral
    sfun = function(pars, n_acc, upper = Inf) {
      return(pstopHybrid(pars, n_acc, upper = upper))
    },
    # Random function for SS race
    rfun = function(data = NULL, pars) {
      return(rSShybrid(data, pars, ok = attr(pars, "ok")))
    },
    # Race likelihood combining pfun and dfun
    log_likelihood = function(pars, dadm, model, min_ll = log(1e-10)) {
      return(log_likelihood_race_ss(pars, dadm, model, min_ll = min_ll))
    }
  )
}


log_likelihood_race_ss <- function(pars,dadm,model,min_ll=log(1e-10))
{
  # All bad?
  if (is.null(attr(pars,"ok")))
    ok <- !logical(dim(pars)[1]) else ok <- attr(pars,"ok")
    if (!any(ok)) return(min_ll*length(attr(dadm, "expand")))

    # Nomenclature for indices:
    # "is" = logical, "isp" pars/dadm index,
    # "t" trials index, "ist" logical on trial index
    # "n_" number of integer

    # Counts
    n_acc <- length(levels(dadm$lR))                   # total number of accumulators
    n_trials <- nrow(dadm)/n_acc                       # number of trials
    n_accG <- sum(as.numeric(dadm[1:n_acc,"lI"])==2)   # go accumulators
    n_accST <- sum(as.numeric(dadm[1:n_acc,"lI"])==1)

    # Likelihood for all trials and for ok trials
    allLL <- rep(min_ll,n_trials)
    # used to put results into allLL[allok]
    allok <- ok[dadm$lR==levels(dadm$lR)[1]]

    # Remove trials not ok
    pars <- pars[ok,,drop=FALSE]
    dadm <- dadm[ok,,drop=FALSE]

    # Only keep OK trials for stop trial computations
    n_trials <- nrow(dadm)/n_acc # number of trials
    trials <- 1:n_trials         # trial number

    # Booleans for ok pars/dadm
    isp1 <- dadm$lR==levels(dadm$lR)[1]      # 1st accumulator rows
    ispGOacc <- dadm$lI==levels(dadm$lI)[2]  # Go accumulator rows
    ispStop <- is.finite(dadm$SSD)           # Stop-trial rows

    # Failure parameters for each trial
    gf <- pars[isp1,"gf"]
    tf <- pars[isp1,"tf"]

    # Go response names
    GoR <- as.character(dadm[1:n_acc,"lR"][dadm[1:n_acc,"lI"]==2])

    # No response
    ispNR <- is.na(dadm$R)
    if ( any(ispNR) ) {  # Note by definition no ST present

      # Go failures
      ispgoNR <- ispNR & !ispStop            # No response and no stop signal
      tgoNR <- c(1:sum(isp1))[ispgoNR[isp1]] # trial number
      if (any(ispgoNR))
        allLL[allok][tgoNR] <- log(gf[tgoNR])

      # Stop trial with no response and no ST accumulator
      ispstopNR <- ispNR & ispStop & (n_accST == 0)
      if ( any(ispstopNR) ) { # Stop trial probability
        # Non-response and stop trial & go/stop accumulator parameters
        pStop <- pmin(1,pmax(0,  # protection to keep in 0-1
                             model$sfun(
                               pars[ispNR & ispStop & ispGOacc,,drop=FALSE],n_acc=n_accG)
        ))
        # Fill in stop-trial non-response probabilities, either 1) go failure
        # 2) not go failure and not trigger failure and stop wins
        tstopNR <- trials[ispstopNR[isp1]]  # trial number
        allLL[allok][tstopNR] <- log(gf[tstopNR] + (1-gf[tstopNR])*(1-tf[tstopNR])*pStop)
      }
    }

    # Response made
    if (any(!ispNR)) {
      # Only keep response trials for further computation
      allr <- !ispNR[isp1] # used to put back into allLL[allok]
      pars <- pars[!ispNR,,drop=FALSE]
      dadm <- dadm[!ispNR,,drop=FALSE]
      n_trials <- nrow(dadm)/n_acc # number of trials
      trials <- 1:n_trials
      ptrials <- rep(trials,each=n_acc) # trial number for pars/dadm

      isp1 <- dadm$lR==levels(dadm$lR)[1]      # 1st accumulator rows
      ispGOacc <- dadm$lI==levels(dadm$lI)[2]  # Go accumulator rows
      ispStop <- is.finite(dadm$SSD)           # stop-trial rows with a response
      gf <- pars[isp1,"gf"]
      tf <- pars[isp1,"tf"]

      # likelihoods for trials with a response (eventually, intermediately probabilities)
      like <- numeric(n_trials)
      lds <- numeric(nrow(dadm)) # log density and survivor, used for both go and stop trials

      # Go trials with response
      if (any(!ispStop)) {
        ispGOwin <-  !ispStop & dadm$winner # Winner go accumulator rows
        tGO <- ptrials[ispGOwin]  # Go trials
        # Winner density
        like[tGO] <- log(model$dfunG(
          rt=dadm$rt[ispGOwin],pars=pars[ispGOwin,,drop=FALSE]))
        if (n_accG >1) {  # Looser survivor go accumulator(s)
          ispGOloss <- !ispStop & !dadm$winner & ispGOacc # Looser go accumulator rows
          like[tGO] <- like[tGO] + apply(matrix(log(1-model$pfunG(
            rt=dadm$rt[ispGOloss],pars=pars[ispGOloss,,drop=FALSE])),nrow=n_accG-1),2,sum)
        }
        # Transform back to densities to include go failure
        like[tGO] <- (1-gf[tGO])*exp(like[tGO])
      }

      # Stop trials with a response, occurs if
      # 1) All triggered
      #   a) Stop does not beat go before rt, produces go and ST
      #   b) Stop beats go before rt, produces ST only
      # 2) go failure, stop triggered, produces ST only
      # 3) go triggered, stop failure, produces go only
      # NB: cant have both go and stop failure as then no response

      if (any(ispStop)) {                       # Stop trials
        ispSGO <- ispStop & (dadm$R %in% GoR)   # Go responses
        ispSST <- ispStop & !(dadm$R %in% GoR)  # ST responses
        # Go beats stop and ST (if any)
        if (any(ispSGO)) {
          ispGOwin <-  ispSGO & dadm$winner # Winner go accumulator rows
          tGO <- ptrials[ispGOwin]          # Go trials
          like[tGO] <- log(model$dfunG(
            rt=dadm$rt[ispGOwin],pars=pars[ispGOwin,,drop=FALSE]))
          if (n_accG > 1) {  # Looser survivor gp accumulators
            ispGOloss <- ispSGO & !dadm$winner & ispGOacc
            like[tGO] <- like[tGO] + apply(matrix(log(1-model$pfunG(
              rt=dadm$rt[ispGOloss],pars=pars[ispGOloss,,drop=FALSE])),
              nrow=n_accG-1),2,sum)
          }
          # trigger stop, add in stop survivor
          ts <- like[tGO] + log(1-model$pfunS(
            rt=dadm$rt[ispGOwin],pars=pars[ispGOwin,,drop=FALSE]))
          # ST loosers
          if (n_accST == 0) stl <- 0 else {
            ispSTloss <- ispSGO & !ispGOacc
            stl <- apply(matrix(log(1-model$pfunG(
              rt=dadm$rt[ispSTloss]-pars[ispSTloss,"SSD"], # correct for SSD
              pars=pars[ispSTloss,,drop=FALSE])),
              nrow=n_accST),2,sum)
          }

          # Transform back to densities to include failures
          like[tGO] <- (1-gf[tGO])*(tf[tGO]*exp(like[tGO]) +
                                      (1-tf[tGO])*exp(ts+stl))
        }

        # ST WINS (never tf)
        if (any(ispSST)) { # Go triggered 1a) with (1-pStop), all race,
          #              1b) with pStop, only ST races
          # Go failure    2) Only ST races
          # ST winner rows
          ispSSTwin <-  dadm$winner &  ispSST
          # Stop probability on ST win trials, only needs go accumulators
          pStop <- model$sfun(pars[ispSST & ispGOacc,,drop=FALSE],n_acc=n_accG,
                              upper=dadm$rt[ispSSTwin]) # pStop before observed RT
          # ST win ll
          tST <- ptrials[ispSSTwin]
          like[tST] <- log(model$dfunG(
            rt=dadm$rt[ispSSTwin]-dadm$SSD[ispSSTwin], # correct ST racers for SSD delay
            pars=pars[ispSSTwin,,drop=FALSE]))
          # ST looser survivors
          if (n_accST > 1) {  # Survivor for looser for ST accumulator(s)
            ispSSTloss <-  !dadm$winner &  ispSST & !ispGOacc
            llST <-  log(1-model$pfunG(
              rt=dadm$rt[ispSSTloss]-dadm$SSD[ispSSTloss],
              pars=pars[ispSSTloss,,drop=FALSE]))
            if (n_accST == 2) # Could remove branch, maybe faster as no matrix sum?
              like[tST] <- like[tST] + llST else
                like[tST] <- like[tST] + apply(matrix(llST,nrow=n_accST-1),2,sum)
          }
          # Go looser survivor
          ispSGloss <- ispSST & ispGOacc
          llG <- apply(matrix(log(1-model$pfunG(
            rt=dadm$rt[ispSGloss],pars=pars[ispSGloss,,drop=FALSE])),
            nrow=n_accG),2,sum)
          like[tST] <- (1-tf[tST])*(                 # Never trigger failure
            gf[tST]*exp(like[tST]) +             # Case 2, gf only ST race
              (1-gf[tST])*(pStop*exp(like[tST]) +      # Case 1b, no gf, stop beats go, ST race
                             (1-pStop)*exp(like[tST]+llG))) # Case 1a, no  gf, all race (no stop win)
        }
      }
      allLL[allok][allr] <- log(like)
    }

    allLL[is.na(allLL)|is.nan(allLL)] <- min_ll
    allLL <- pmax(min_ll,allLL)

    # # TODO remove when we're happy with testing
    # llR <<- allLL

    sum(allLL[attr(dadm,"expand")])
}
