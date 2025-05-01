#### Multiple threshold models ----

# Called in design_model to add "dL" attribute to dadm
get_dL <- function(NR,type="BE") {
  # BE = balance of evidence
  # TC = threshold count

  even <- function(n) n%%2==0

  getRR <- function(RR,nr,type) {

    opposite <- function(n) c(2:1)[n]

    TR <- (NR+1)/2           # dividing rating response
    R <- ceiling(abs(RR-TR)) # rating strength
    if (type=="BE") { # Balance of evidence
      sa <- as.numeric(RR > TR)+1 # stopping accumulator
      nsa <- opposite(sa)  # not stopping accumulator
      tt <- nr             # triggering threshold
      lt <- nr - R         # lower threshold
    } else { # Threshold counting
      if (even(RR)) {
        sa <- 2
        tt <- RR/2
        lt <- nr - (RR %/% 2)
      } else {
        sa <- 1
        tt <- (nr+1) - (RR+1)/2
        lt <- RR %/% 2
      }
      nsa <- opposite(sa)
    }
    ut <- lt + 1         # upper threshold
    c(sa=sa,tt=tt+1,nsa=nsa,lt=lt+1,ut=ut+1) # move from 0 to 1 base for thresholds
  }


  if (!even(NR)) {
    # Use next NR up then adjust
    NR <- NR+1
    nr = NR/2
    out <- t(sapply(c(1:(2*nr)),getRR,nr=nr,type=type))
    return(cbind(RR=c(1:nr,nr:(NR-1)),out)) # collapse two middle
  } else {
    nr = NR/2   # rating strengths least 1..nr most
    out <- t(sapply(c(1:(2*nr)),getRR,nr=nr,type=type))
    return(cbind(RR=1:NR,out))
  }
}


# get_dL(4)
# get_dL(4,"TC")
# get_dL(5)
# get_dL(5,"TC")
# get_dL(6)
# get_dL(6,"TC")
# get_dL(7)
# get_dL(7,"TC")
# get_dL(8)
# get_dL(8,"TC")


#### Random (BE = balance of evidence, TC = threshold count) ----

# Probabilities of guessing rating 1..n from pg1, pg2 ... parameters
get_pgm <- function(pm) {
  ps <- 0
  for (i in 2:ncol(pm)) {
    ps <- ps + pm[,i-1]
    pm[,i] <- (1-ps)*pm[,i]
  }
  pm
}
guess <- function(out,pars) {

  pg <- function(pm) {
    pm <- get_pgm(pm)
    apply(apply(pm,1,cumsum),2,\(x) .bincode(runif(1),c(0,x,1)))
  }

  isguess <- runif(nrow(out)) < pars[is2,"gp"]
  if (any(isguess)) {
    nt <- sum(isguess)
    p <- pars[isguess,c("v", "sv", "gb", "A", "t0"),drop=FALSE]
    names(p)[3] <- "b"
    Rrt <- rLBA(factor(rep(1:2,nt)),p)
    Rrt$R <- as.numeric(Rrt$R)
    isp2 <- 1:nrow(pars) %% 2 == 0
    isR2 <- Rrt$R==2
    nc <- sum(ispg)+1
    pgnams <- sort(colnames(pars)[substr(colnames(pars),1,2)=="pg"])
    Rrt$R[!isR2] <- (nc+1) - pg(pars[!isp2 & isguess,pgnams,drop=FALSE][!isR2,,drop=FALSE])
    Rrt$R[ isR2] <-     nc + pg(pars[ isp2 & isguess,pgnams,drop=FALSE][ isR2,,drop=FALSE])
    out[isguess,c("R","rt")] <- Rrt
  }
  out
}

rLBA_BE <- function(lR,pars,ok=rep(TRUE,dim(pars)[1]),return_activation = FALSE)
  # Balance of evidence
  # Rating from 1=high confidence left to length(d1)+length(d2)+2 for
  # high confidence right, assumes posdrift=TRUE.
  # pars is a matrix with columns A, b, t0, v, sv and r sorted so parameters for
  # each trial are in contiguous rows. Returns nrow(pars)/2 trials
  # d1 and d2 are matrices of lower thresholds for first and second accumulator
  # with one row per trial.
{

  getq <- function(pmat)
    # sample from mvtnorm
  {
    sv2 <- pmat[,"sv"]^2
    sigma <- matrix(c(sv2[1], rep(prod(pmat[,"sv"]) * pmat[1,"r"], 2), sv2[2]),2, 2, byrow = TRUE)
    rtmvnorm(n = 1, mean = pmat[,"v"], sigma = sigma, lower = rep(0, 2))
  }

  get_rate <- function(x){.bincode(x[1],x[-1])}

  bad <- rep(NA, dim(pars)[1]/2)
  out <- data.frame(R = bad, rt = bad)
  pars <- pars[ok,]

  pmats <- array(pars,dim=c(2,nrow(pars)/2,ncol(pars)),
                 dimnames=list(NULL,NULL,colnames(pars)))
  n <- dim(pmats)[2] # number of trials
  q <- apply(pmats,2,getq)
  A <- matrix(runif(length(pmats[,,"A"]),min=0,max=pmats[,,"A"]),nrow=2)
  rts <- (pmats[,,"b"]-A) / q
  rts <- rts + pmats[,,"t0"]
  response <- apply(rts, 2, which.min)
  rt <- rts[cbind(response,seq_len(n))]
  # rating
  looser2 <- response==1
  lmat <- cbind(1 + as.numeric(looser2),seq_len(n))
  ylooser <- A[lmat] + q[lmat]*(rt-pmats[1,,"t0"])
  rating <- numeric(length(response))
  # 2=nacc x n x #DT
  d <- pmats[,,substr(dimnames(pmats)[[3]],1,2)=="DT",drop=FALSE]
  nr <- dim(d)[3]+1 # number of ratings = DT1, DT2 ..., b
  if (any(looser2)) {
    if (sum(looser2)==1)
      rating[looser2] <- get_rate(
        c(ylooser[looser2],rep(0,sum(looser2)),d[2,looser2,],pmats[2,looser2,"b"])) else
      rating[looser2] <- apply(
        cbind(ylooser[looser2],rep(0,sum(looser2)),d[2,looser2,],pmats[2,looser2,"b"]),1,get_rate)
  }
  if (any(!looser2)) {
    if (sum(!looser2)==1)
      rating[!looser2] <- (2*nr+1)-get_rate(
        c(ylooser[!looser2],rep(0,sum(!looser2)),d[1,!looser2,],pmats[2,!looser2,"b"])) else
      rating[!looser2] <- (2*nr+1)-apply(
        cbind(ylooser[!looser2],rep(0,sum(!looser2)),d[1,!looser2,],pmats[2,!looser2,"b"]),1,get_rate)
  }
  ok <- matrix(ok,nrow=2)[1,]
  out[ok,] <- cbind(response = rating,rt = rt)
  if (return_activation) {
    tmp <- cbind(out[ok,],y = t(A) + t(q) * (rt - pmats[1,,"t0"]))
    out <- cbind(out,y1=bad,y2=bad)
    out[ok,] <- tmp
  }
  out[ok,] <- guess(out[ok,],pars[ok,])
  out$R <- factor(out$R,levels=1:(2*nr))
  if (!(length(lR)%%2==0)) {
    levels(out$R) <- c(1:nr,nr:(2*nr-1))
  }
  out$R <- factor(out$R,levels=lR)
  out
}

rLBA_TC <- function(lR,pars,ok=rep(TRUE,dim(pars)[1]),return_activation = FALSE)
  # Threshold count
  # Rating from 1=high confidence left to length(d1)+length(d2)+2 for
  # high confidence right, assumes posdrift=TRUE.
  # pars is a matrix with columns A, b, t0, v, sv and r sorted so parameters for
  # each trial are in contiguous rows. Returns nrow(pars)/2 trials
  # d1 and d2 are matrices of lower thresholds for first and second accumulator
  # with one row per trial.
{

  getq <- function(pmat)
    # sample from mvtnorm
  {
    sv2 <- pmat[,"sv"]^2
    sigma <- matrix(c(sv2[1], rep(prod(pmat[,"sv"]) * pmat[1,"r"], 2), sv2[2]),2, 2, byrow = TRUE)
    rtmvnorm(n = 1, mean = pmat[,"v"], sigma = sigma, lower = rep(0, 2))
  }

  bad <- rep(NA, dim(pars)[1]/2)
  out <- data.frame(R = bad, rt = bad)
  pars <- pars[ok,]

  Dnams <- dimnames(pars)[[2]][substr(dimnames(pars)[[2]],1,2)=="DT"]
  tmats <- array(cbind(pars[,c(Dnams,"b")]),dim=c(2,nrow(pars)/2,1+length(Dnams)),
                 dimnames=list(NULL,NULL,c(Dnams,"b")))
  pnams <- dimnames(pars)[[2]][!(dimnames(pars)[[2]] %in% c(Dnams,"b"))]
  pmats <- array(pars[,pnams],dim=c(2,nrow(pars)/2,length(pnams)),
                 dimnames=list(NULL,NULL,pnams))

  nr <- dim(tmats)[3] # number of ratings = triggering count
  nt=nr*2

  # rating lookup
  dL <- matrix(nrow=nt,ncol=2)
  even <- array(c(1:nt),dim=c(2,nr))
  odd <- as.vector(even[1,])
  even <- as.vector(even[2,])
  dL[odd,1] <- 1; dL[even,1] <- 2
  dL[odd,2] <- nr:1; dL[even,2] <- 1:nr
  rate <- dL[,1] + (dL[,2]-1)*2
  n <- dim(pmats)[2]        # number of trials
  q <- apply(pmats,2,getq)  # rates
  A <- matrix(runif(length(pmats[,,"A"]),min=0,max=pmats[,,"A"]),nrow=2)
  rts=array(apply(tmats,3,function(x){(x-A)/q}),dim=dim(tmats))
  rts[rts<0] <- 0
  ok <- matrix(ok,nrow=2)[1,]
  out[ok,] <- t(apply(rts,2,function(x){
    rt <- sort(x)[nr]
    c(rate[x==rt],rt)
  }))
  out[ok,"rt"] <- out[ok,"rt"] + pmats[1,,"t0"]
  if (return_activation) {
    tmp <- cbind(out[ok,],y = t(A) + t(q) * (out[ok,"rt"] - pmats[1,,"t0"]))
    out <- cbind(out,y1=bad,y2=bad)
    out[ok,] <- tmp
  }
  out[ok,] <- guess(out[ok,],pars[ok,])
  out$R <- factor(out$R,levels=1:(2*nr))
  if (!(length(lR)%%2==0)) {
    levels(out$R) <- c(1:nr,nr:(2*nr-1))
  }
  out$R <- factor(out$R,levels=lR)
  out
}

#### Density (cdf implict as multi-variate) ----

# Taken from MomTrunc
recintab0 <- function (kappa, a, b, mu, Sigma)
{
  p = length(kappa)
  seqq = seq_len(p)
  if (p == 1) {
    M = matrix(data = 0, nrow = kappa + 1, ncol = 1)
    s1 = sqrt(Sigma)
    aa = (a - mu)/s1
    bb = (b - mu)/s1
    M[1] = pnorm(bb) - pnorm(aa)
    if (kappa > 0) {
      pdfa = s1 * dnorm(aa)
      pdfb = s1 * dnorm(bb)
      M[2] = mu * M[1] + pdfa - pdfb
      if (a == -Inf) {
        a = 0
      }
      if (b == Inf) {
        b = 0
      }
      for (i in seq_len(kappa - 1) + 1) {
        pdfa = pdfa * a
        pdfb = pdfb * b
        M[i + 1] = mu * M[i] + (i - 1) * Sigma * M[i -
                                                     1] + pdfa - pdfb
      }
    }
  }
  else {
    pk = prod(kappa + 1)
    M = rep(NA, pk)
    nn = round(pk/(kappa + 1), 0)
    begind = cumsum(c(0, nn))
    pk1 = begind[p + 1]
    cp = matrix(0, p, p)
    for (i in seqq) {
      kk = kappa
      kk[i] = 0
      cp[i, ] = c(1, cumprod(kk[1:p - 1] + 1))
    }
    G = rep(0, pk1)
    H = rep(0, pk1)
    s = sqrt(diag(Sigma))
    pdfa = dnorm(x = a, mean = mu, sd = s)
    pdfb = dnorm(x = b, mean = mu, sd = s)
    for (i in seqq) {
      ind2 = seqq[-i]
      kappai = kappa[ind2]
      ai = a[ind2]
      bi = b[ind2]
      mui = mu[ind2]
      Si = Sigma[ind2, i]
      SSi = Sigma[ind2, ind2] - Si %*% t(Si)/Sigma[i,
                                                   i]
      ind = (begind[i] + 1):begind[i + 1]
      if (a[i] != -Inf) {
        mai = mui + Si/Sigma[i, i] * (a[i] - mu[i])
        G[ind] = pdfa[i] * recintab0(kappai, ai, bi,
                                     mai, SSi)
      }
      if (b[i] != Inf) {
        mbi = mui + Si/Sigma[i, i] * (b[i] - mu[i])
        H[ind] = pdfb[i] * recintab0(kappai, ai, bi,
                                     mbi, SSi)
      }
    }
    M[1] = pmvnorm(lower = a, upper = b, mean = mu, sigma = Sigma)
    a[a == -Inf] = 0
    b[b == Inf] = 0
    cp1 = cp[p, ]
    dims = lapply(X = kappa + 1, FUN = seq_len)
    grid = do.call(expand.grid, dims)
    for (i in seq_len(pk - 1) + 1) {
      kk = as.numeric(grid[i, ])
      ii = sum((kk - 1) * cp1) + 1
      i1 = min(seqq[kk > 1])
      kk1 = kk
      kk1[i1] = kk1[i1] - 1
      ind3 = ii - cp1[i1]
      M[ii] = mu[i1] * M[ind3]
      for (j in seqq) {
        kk2 = kk1[j] - 1
        if (kk2 > 0) {
          M[ii] = M[ii] + Sigma[i1, j] * kk2 * M[ind3 -
                                                   cp1[j]]
        }
        ind4 = as.numeric(begind[j] + sum(cp[j, ] *
                                            (kk1 - 1)) - cp[j, j] * kk2 + 1)
        M[ii] = M[ii] + Sigma[i1, j] * (a[j]^kk2 * G[ind4] -
                                          b[j]^kk2 * H[ind4])
      }
    }
  }
  return(M)
}

n1PDF_MTR_1 <- function(rt, pars,dl,du,b)
  # Works for Balance of Evidence and Threshold count.
  # Assumes single rt value, posdrift = TRUE
  # pars is a 2 row matrix of parameters with columns A, t0, v, sv, r
  # with first row corresponding to the accumulator that stops the race
  # b is stopping threshold, dl and du are the upper and lower thresholds for
  # the accumulator that does not stop the race.
  # r and t0 same for both (second ignored)

{

  if (rt <= pars[1,"t0"]) {
    return(0)
  }
  n <- 2
  t <- rt - pars[1,"t0"]
  sv2 <- pars[,"sv"]^2
  sigma <- matrix(c(sv2[1], rep(prod(pars[,"sv"]) * pars[1,"r"], 2), sv2[2]),2, 2, byrow = TRUE)

  Pqpos <- recintab0(kappa = rep(0, n),
                                mu = pars[,"v"],
                                Sigma = sigma,
                                a = rep(0, n),
                                b = rep(Inf, n))

  # Integral bounds
  lb <- ub <- numeric(4)
  lb[1] <- 0
  lb[2] <- ifelse((du - pars[2,"A"]) / t <= 0, - du / pars[2,"A"], -1)
  lb[3] <- 0
  lb[4] <- ifelse((dl - pars[2,"A"]) / t <= 0, - dl / pars[2,"A"], -1)

  ub[1] <- ifelse((du - pars[2,"A"]) / t <= 0, 0, (du - pars[2,"A"]) / t)
  ub[2] <- ifelse(du / t          <= 0, - du / pars[2,"A"], 0)
  ub[3] <- ifelse((dl - pars[2,"A"]) / t <= 0, 0, (dl - pars[2,"A"]) / t)
  ub[4] <- ifelse(dl / t       <= 0, - dl / pars[2,"A"], 0)

  # signs
  signs <- c(1, -1, -1, 1)

  # Cs
  cs <- cbind(NA, c(1, t / pars[2,"A"]), NA, c(1, t / pars[2,"A"]))

  # ms
  ms <- cbind(NA, c(0, - du / pars[2,"A"]), NA, c(0, - dl / pars[2,"A"]))

  # index for which integrals relevant
  index <- which(lb != ub)

  ints <- sapply(index, function(x) {

    if (x == 1 || x == 3) {

      mom <- recintab0(kappa = c(1, 0),
                                  mu = pars[,"v"],
                                  Sigma = sigma,
                                  a = c((b - pars[1,"A"]) / t, lb[x]),
                                  b = c( b / t,         ub[x]))
      return(signs[x] * mom[length(mom)])

    } else if (x == 2 || x == 4) {

      C_j <- diag(cs[,x])
      m_j <- ms[,x]

      m_x <- as.vector(m_j + C_j %*% pars[,"v"])
      sigma_x <- C_j %*% sigma %*% t(C_j)

      mom <- recintab0(kappa = rep(1, 2),
                                  mu = m_x,
                                  Sigma = sigma_x,
                                  a = c((b - pars[1,"A"]) / t, lb[x]),
                                  b = c( b / t,         ub[x]))
      return(signs[x] * mom[length(mom)])

    }

  })

  out <- sum(ints) / (Pqpos * pars[1,"A"])

  return(pmax(out,0))

}



#' Balance of Evidence 2 Threshold LBA model
#'
#' All thresholds greater than A, DT1 <- A+DT1, b <- DT1 + B
#'
#' @return A model list with all the necessary functions to sample
#' @export
BE2LBA <- function(){
  list(
    type="BE",
    p_types=c("v" = 1,"sv" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0),
              r=qnorm(.5),DT1=log(.5),
              gB=log(1),gp=qnorm(1),pg1=qnorm(1)),
    transform=list(func=c(v = "identity",sv = "exp", B = "exp", A = "exp",t0 = "exp",
                          r="pnorm",DT1="exp",
                          gB="exp",gp="pnorm",pg1="pnorm")),
    bound=list(minmax=cbind(v=c(-Inf,Inf),sv = c(0, Inf), A=c(1e-4,Inf),b=c(0,Inf),
      t0=c(0.05,Inf),r=c(-1,1),DT1=c(0,Inf),
      gb=c(0,Inf),gp=c(0,1),pg1=c(0,1)),exception=c(A=0)),
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars[,"r"] <- 2*pars[,"r"]-1
      pars <- cbind(pars,gb=pars[,"gB"] + pars[,"A"])
      pnams <- c("A",dimnames(pars)[[2]][substr(dimnames(pars)[[2]],1,2)=="DT"],"B")
      pars[,pnams] <- t(apply(pars[,pnams],1,cumsum))
      dimnames(pars)[[2]][dimnames(pars)[[2]]=="B"] <- "b"
      pars
    },
    # Random function for racing accumulator
    rfun=function(lR=NULL,pars) {
      rLBA_BE(lR,pars,ok=attr(pars, "ok"))
    },
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,model,min_ll=log(1e-10)){
      log_likelihood_mt(pars=p_vector, dadm = dadm, model=model, min_ll = min_ll,n_cores=3)
    }
  )
}



#' Balance of Evidence 2 Threshold LBA model
#'
#' Intermediate thresholds proportional to b, can be less than A
#'
#' @return A model list with all the necessary functions to sample
#' @export
BE2PLBA <- function(){
  list(
    type="BE",
    p_types=c("v" = 1,"sv" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0),
              r=qnorm(.5),DT1=pnorm(.5),
              gB=log(1),gp=qnorm(0),pg1=qnorm(1)),
    transform=list(func=c(v = "identity",sv = "exp", B = "exp", A = "exp",t0 = "exp",
                          r="pnorm",DT1="pnorm",
                          gB="exp",gp="pnorm",pg1="pnorm")),
    bound=list(minmax=cbind(v=c(-Inf,Inf),sv = c(0, Inf), A=c(1e-4,Inf),b=c(0,Inf),
      t0=c(0.05,Inf),r=c(-1,1),DT1=c(0,1),
      gb=c(0,Inf),gp=c(0,1),pg1=c(0,1)),exception=c(A=0)),
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars[,"r"] <- 2*pars[,"r"]-1
      pars <- cbind(pars,b=pars[,"B"] + pars[,"A"],gb=pars[,"gB"] + pars[,"A"])
      isDT <- substr(dimnames(pars)[[2]],1,2)=="DT"
      pars[,isDT] <- pars[,isDT]*pars[,"b"]
      dimnames(pars)[[2]][dimnames(pars)[[2]]=="B"] <- "b"
      pars
    },
    # Random function for racing accumulator
    rfun=function(lR=NULL,pars) {
      rLBA_BE(lR,pars,ok=attr(pars, "ok"))
    },
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,model,min_ll=log(1e-10)){
      log_likelihood_mt(pars=p_vector, dadm = dadm, model=model, min_ll = min_ll,n_cores=3)
    }
  )
}


#' Threshold Count 2 LBA model
#'
#' All thresholds greater than A, DT1 <- A+DT1, b <- DT1 + B
#'
#' @return A model list with all the necessary functions to sample
#' @export
TC2LBA <- function(){
  list(
    type="TC",
    p_types=c("v" = 1,"sv" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0),
              r=qnorm(.5),DT1=log(.5),
              gB=log(1),gp=qnorm(0),pg1=qnorm(1)),
    transform=list(func=c(v = "identity",sv = "exp", B = "exp", A = "exp",t0 = "exp",
                          r="pnorm",DT1="exp",
                          gB="exp",gp="pnorm",pg1="pnorm")),
    bound=list(minmax=cbind(v=c(-Inf,Inf),sv = c(0, Inf), A=c(1e-4,Inf),b=c(0,Inf),
      t0=c(0.05,Inf),r=c(-1,1),DT1=c(0,Inf),
      gb=c(0,Inf),gp=c(0,1),pg1=c(0,1)),exception=c(A=0)),
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars[,"r"] <- 2*pars[,"r"]-1
      pars <- cbind(pars,gb=pars[,"gB"] + pars[,"A"])
      pnams <- c("A",dimnames(pars)[[2]][substr(dimnames(pars)[[2]],1,2)=="DT"],"B")
      pars[,pnams] <- t(apply(pars[,pnams],1,cumsum))
      dimnames(pars)[[2]][dimnames(pars)[[2]]=="B"] <- "b"
      pars
    },
    # Random function for racing accumulator
    rfun=function(lR=NULL,pars) {
      rLBA_TC(lR,pars,ok=attr(pars, "ok"))
    },
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,model,min_ll=log(1e-10)){
      log_likelihood_mt(pars=p_vector, dadm = dadm, model=model, min_ll = min_ll,n_cores=3)
    }
  )
}


#' Balance of Evidence 3 Threshold LBA model
#'
#' All thresholds greater than A, DT1 <- A+DT1, DT2 <- DT1 + DT2, b <- DT2 + B
#'
#' @return A model list with all the necessary functions to sample
#' @export
BE3LBA <- function(){
  list(
    type="BE",
    p_types=c("v" = 1,"sv" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0),
              r=qnorm(.5),DT1=log(1/3),DT2=log(1/3),
              gB=log(1),gp=qnorm(0),pg1=qnorm(1),pg2=qnorm(0)),
    transform=list(func=c(v = "identity",sv = "exp", B = "exp", A = "exp",t0 = "exp",
                          r="pnorm",DT1="exp",DT2="exp",
                          gB="exp",gp="pnorm",g1="pnorm",g2="pnorm")),
    bound=list(minmax=cbind(v=c(-Inf,Inf),sv = c(0, Inf), A=c(1e-4,Inf),b=c(0,Inf),
      t0=c(0.05,Inf),r=c(-1,1),DT1=c(0,Inf),DT2=c(0,Inf),
      gb=c(0,Inf),gp=c(0,1),pg1=c(0,1)),exception=c(A=0)),
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars[,"r"] <- 2*pars[,"r"]-1
      pars <- cbind(pars,gb=pars[,"gB"] + pars[,"A"])
      pnams <- c("A",dimnames(pars)[[2]][substr(dimnames(pars)[[2]],1,2)=="DT"],"B")
      pars[,pnams] <- t(apply(pars[,pnams],1,cumsum))
      dimnames(pars)[[2]][dimnames(pars)[[2]]=="B"] <- "b"
      pars
    },
    # Random function for racing accumulator
    rfun=function(lR=NULL,pars) {
      rLBA_BE(lR,pars,ok=attr(pars, "ok"))
    },
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,model,min_ll=log(1e-10)){
      log_likelihood_mt(pars=p_vector, dadm = dadm, model=model, min_ll = min_ll,n_cores=3)
    }
  )
}

#' Balance of Evidence 3 Threshold LBA model
#'
#' Intermediate thresholds proportional to b, can be less than A
#'
#' @return A model list with all the necessary functions to sample
#' @export
BE3PLBA <- function(){
  list(
    type="BE",
    p_types=c("v" = 1,"sv" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0),
              r=qnorm(.5),DT1=pnorm(1/3),DT2=pnorm(2/3),
              gB=log(1),gp=qnorm(0),pg1=qnorm(1),pg2=qnorm(0)),
    transform=list(func=c(v = "identity",sv = "exp", B = "exp", A = "exp",t0 = "exp",
                          r="pnorm",DT1="pnorm",DT2="pnorm",
                          gB="exp",gp="pnorm",pg1="pnorm",pg2="pnorm")),
    bound=list(minmax=cbind(v=c(-Inf,Inf),sv = c(0, Inf), A=c(1e-4,Inf),b=c(0,Inf),
      t0=c(0.05,Inf),r=c(-1,1),DT1=c(0,1),DT2=c(0,1),
      gb=c(0,Inf),gp=c(0,1),pg1=c(0,1),pg2=c(0,1)),exception=c(A=0)),
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars[,"r"] <- 2*pars[,"r"]-1
      pars <- cbind(pars,b=pars[,"B"] + pars[,"A"],gb=pars[,"gB"] + pars[,"A"])
      isDT <- substr(dimnames(pars)[[2]],1,2)=="DT"
      pars[,isDT] <- pars[,isDT]*pars[,"b"]
      dimnames(pars)[[2]][dimnames(pars)[[2]]=="B"] <- "b"
      pars
    },
    # Random function for racing accumulator
    rfun=function(lR=NULL,pars) {
      rLBA_BE(lR,pars,ok=attr(pars, "ok"))
    },
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,model,min_ll=log(1e-10)){
      log_likelihood_mt(pars=p_vector, dadm = dadm, model=model, min_ll = min_ll,n_cores=3)
    }
  )
}

#' Threshold Count 3 LBA model
#'
#' All thresholds greater than A, DT1 <- A+DT1, DT2 <- DT1 + DT2, b <- DT2 + B
#'
#' @return A model list with all the necessary functions to sample
#' @export
TC3LBA <- function(){
  list(
    type="TC",
    p_types=c("v" = 1,"sv" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0),
              r=qnorm(.5),DT1=log(1/3),DT2=log(1/3),
              gB=log(1),gp=qnorm(0),pg1=qnorm(1),pg2=qnorm(0)),
    transform=list(func=c(v = "identity",sv = "exp", B = "exp", A = "exp",t0 = "exp",
                          r="pnorm",DT1="exp",DT2="exp",
                          gB="exp",gp="pnorm",pg1="pnorm",pg2="pnorm")),
    bound=list(minmax=cbind(v=c(-Inf,Inf),sv = c(0, Inf), A=c(1e-4,Inf),b=c(0,Inf),
      t0=c(0.05,Inf),r=c(-1,1),DT1=c(0,Inf),DT2=c(0,Inf),
      gb=c(0,Inf),gp=c(0,1),pg1=c(0,1),pg2=c(0,1)),exception=c(A=0)),
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars[,"r"] <- 2*pars[,"r"]-1
      pars <- cbind(pars,gb=pars[,"gB"] + pars[,"A"])
      pnams <- c("A",dimnames(pars)[[2]][substr(dimnames(pars)[[2]],1,2)=="DT"],"B")
      pars[,pnams] <- t(apply(pars[,pnams],1,cumsum))
      dimnames(pars)[[2]][dimnames(pars)[[2]]=="B"] <- "b"
      pars
    },
    # Random function for racing accumulator
    rfun=function(lR=NULL,pars) {
      rLBA_TC(lR,pars,ok=attr(pars, "ok"))
    },
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,model,min_ll=log(1e-10)){
      log_likelihood_mt(pars=p_vector, dadm = dadm, model=model, min_ll = min_ll,n_cores=3)
    }
  )
}


#' Balance of Evidence 4 Threshold LBA model
#'
#' All thresholds greater than A, DT1 <- A+DT1, DT2 <- DT1 + DT2, DT3 <- DT2 + DT3, b <- DT3 + B
#'
#' @return A model list with all the necessary functions to sample
#' @export
BE4LBA <- function(){
  list(
    type="BE",
    p_types=c("v" = 1,"sv" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0),
              r=qnorm(.5),DT1=log(1/4),DT2=log(1/4),DT3=log(1/4),
              gB=log(1),gp=qnorm(0),pg1=qnorm(1),pg2=qnorm(0),pg3=qnorm(0)),
    transform=list(func=c(v = "identity",sv = "exp", B = "exp", A = "exp",t0 = "exp",
                          r="pnorm",DT1="exp",DT2="exp",DT3="exp",
                          gB="exp",gp="pnorm",pg1="pnorm",pg2="pnorm",pg3="pnorm")),
    bound=list(minmax=cbind(v=c(-Inf,Inf),sv = c(0, Inf), A=c(1e-4,Inf),b=c(0,Inf),
      t0=c(0.05,Inf),r=c(-1,1),DT1=c(0,Inf),DT2=c(0,Inf),DT3=c(0,Inf),
      gb=c(0,Inf),gp=c(0,1),pg1=c(0,1),pg2=c(0,1),pg3=c(0,1)),exception=c(A=0)),
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars[,"r"] <- 2*pars[,"r"]-1
      pars <- cbind(pars,gb=pars[,"gB"] + pars[,"A"])
      pnams <- c("A",dimnames(pars)[[2]][substr(dimnames(pars)[[2]],1,2)=="DT"],"B")
      pars[,pnams] <- t(apply(pars[,pnams],1,cumsum))
      dimnames(pars)[[2]][dimnames(pars)[[2]]=="B"] <- "b"
      pars
    },
    # Random function for racing accumulator
    rfun=function(lR=NULL,pars) {
      rLBA_BE(lR,pars,ok=attr(pars, "ok"))
    },
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,model,min_ll=log(1e-10)){
      log_likelihood_mt(pars=p_vector, dadm = dadm, model=model, min_ll = min_ll,n_cores=3)
    }
  )
}

#' Balance of Evidence 4 Threshold LBA model
#'
#' Intermediate thresholds proportional to b, can be less than A
#'
#' @return A model list with all the necessary functions to sample
#' @export
BE4PLBA <- function(){
  list(
    type="BE",
    p_types=c("v" = 1,"sv" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0),
              r=qnorm(.5),DT1=pnorm(1/4),DT2=pnorm(1/2),DT3=pnorm(3/4),
              gB=log(1),gp=qnorm(0),pg1=qnorm(1),pg2=qnorm(0),pg3=qnorm(0)),
    transform=list(func=c(v = "identity",sv = "exp", B = "exp", A = "exp",t0 = "exp",
                          r="pnorm",DT1="pnorm",DT2="pnorm",DT3="pnorm",
                          gB="exp",gp="pnorm",pg1="pnorm",pg2="pnorm",pg3="pnorm")),
    bound=list(minmax=cbind(v=c(-Inf,Inf),sv = c(0, Inf), A=c(1e-4,Inf),b=c(0,Inf),
      t0=c(0.05,Inf),r=c(-1,1),DT1=c(0,1),DT2=c(0,1),DT3=c(0,1),
      gb=c(0,Inf),gp=c(0,1),pg1=c(0,1),pg2=c(0,1),pg3=c(0,1)),exception=c(A=0)),
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars[,"r"] <- 2*pars[,"r"]-1
      pars <- cbind(pars,b=pars[,"B"] + pars[,"A"],gb=pars[,"gB"] + pars[,"A"])
      isDT <- substr(dimnames(pars)[[2]],1,2)=="DT"
      pars[,isDT] <- pars[,isDT]*pars[,"b"]
      dimnames(pars)[[2]][dimnames(pars)[[2]]=="B"] <- "b"
      pars
    },
    # Random function for racing accumulator
    rfun=function(lR=NULL,pars) {
      rLBA_BE(lR,pars,ok=attr(pars, "ok"))
    },
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,model,min_ll=log(1e-10)){
      log_likelihood_mt(pars=p_vector, dadm = dadm, model=model, min_ll = min_ll,n_cores=3)
    }
  )
}


#' Threshold Count 4 LBA model
#'
#' All thresholds greater than A, DT1 <- A+DT1, DT2 <- DT1 + DT2, DT3 <- DT2 + DT3, b <- DT3 + B
#'
#' @return A model list with all the necessary functions to sample
#' @export
TC4LBA <- function(){
  list(
    type="TC",
    p_types=c("v" = 1,"sv" = log(1),"B" = log(1),"A" = log(0),"t0" = log(0),
              r=qnorm(.5),DT1=log(1/4),DT2=log(1/4),DT3=log(1/4),
              gB=log(1),gp=qnorm(0),pg1=qnorm(1),pg2=qnorm(0),pg3=qnorm(0)),
    transform=list(func=c(v = "identity",sv = "exp", B = "exp", A = "exp",t0 = "exp",
                          r="pnorm",DT1="exp",DT2="exp",DT3="exp",
                          gB="exp",gp="pnorm",pg1="pnorm",pg2="pnorm",pg3="pnorm")),
    bound=list(minmax=cbind(v=c(-Inf,Inf),sv = c(0, Inf), A=c(1e-4,Inf),b=c(0,Inf),
      t0=c(0.05,Inf),r=c(-1,1),DT1=c(0,Inf),DT2=c(0,Inf),DT3=c(0,Inf),
      gb=c(0,Inf),gp=c(0,1),pg1=c(0,1),pg2=c(0,1),pg3=c(0,1)),exception=c(A=0)),
    # Trial dependent parameter transform
    Ttransform = function(pars,dadm) {
      pars[,"r"] <- 2*pars[,"r"]-1
      pars <- cbind(pars,gb=pars[,"gB"] + pars[,"A"])
      pnams <- c("A",dimnames(pars)[[2]][substr(dimnames(pars)[[2]],1,2)=="DT"],"B")
      pars[,pnams] <- t(apply(pars[,pnams],1,cumsum))
      dimnames(pars)[[2]][dimnames(pars)[[2]]=="B"] <- "b"
      pars
    },
    # Random function for racing accumulator
    rfun=function(lR=NULL,pars) {
      rLBA_TC(lR,pars,ok=attr(pars, "ok"))
    },
    # Race likelihood combining pfun and dfun
    log_likelihood=function(p_vector,dadm,model,min_ll=log(1e-10)){
      log_likelihood_mt(pars=p_vector, dadm = dadm, model=model, min_ll = min_ll,n_cores=3)
    }
  )
}
