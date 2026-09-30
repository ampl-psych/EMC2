# VERBATIM extract (extract_refs.R) from the don't-know paper's DMC code:
#   DMC_DK_OSF/rod/dmc/models/MTR-LBA/dists.R  (https://osf.io/download/akzxw/)
# Do not edit.
# md5 of source file: f8d1312767e59ae3716485f121d6303c

# ---- n1PDFfixedt0.DK (source lines 190-252) ----
n1PDFfixedt0.DK <- function(dt,A,b,d,mean_v,sd_v, 
                            r2="",posdrift=TRUE,robust = FALSE) 
  # Generates defective PDF for responses on node= 1
  # dt (decison time) is a matrix with length(mean_v) rows, one row for
  # each accumulator to allow for different start times
{
  
  n_acc <- 2
  if (is.null(dim(dt))) dt <- matrix(rep(dt,each=n_acc),nrow=n_acc)
  if (!is.matrix(mean_v)) mean_v <- matrix(rep(mean_v,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(sd_v)) sd_v <- matrix(rep(sd_v,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(b)) b <- matrix(rep(b,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(d)) d <- matrix(rep(d,dim(dt)[2]),nrow=n_acc)
  if (!is.matrix(A)) A <- matrix(rep(A,dim(dt)[2]),nrow=n_acc)
  
  r2 <- rep(r2,length.out=dim(dt)[2])
  is.dk <- r2 =="DK"
  
  #   D1b <- dlba.norm(dt[1,!is.dk],A=A[1,!is.dk],b=b[1,!is.dk],mean_v=mean_v[1,!is.dk],
  #       sd_v=sd_v[1,!is.dk],posdrift=posdrift,robust=robust)
  #   S2d <- (1-plba.norm(dt[2,!is.dk],A=A[2,!is.dk],b=d[2,!is.dk],
  #       mean_v=mean_v[2,!is.dk],sd_v=sd_v[2,!is.dk],posdrift=posdrift,robust=robust))
  # 
  #   if ( any(!is.dk) ) dt[1,!is.dk] <- D1b*S2d
  # 
  #   D1b <- dlba.norm(dt[1,is.dk],A=A[1,is.dk],b=b[1,is.dk],mean_v=mean_v[1,is.dk],
  #       sd_v=sd_v[1,is.dk],posdrift=posdrift,robust=robust)
  #   D2b <- dlba.norm(dt[2,is.dk],A=A[2,is.dk],b=b[2,is.dk],mean_v=mean_v[2,is.dk],
  #       sd_v=sd_v[2,is.dk],posdrift=posdrift,robust=robust)
  # 
  #   S1d <- (1-plba.norm(dt[1,is.dk],A=A[1,is.dk],b=d[1,is.dk],
  #       mean_v=mean_v[1,is.dk],sd_v=sd_v[1,is.dk],posdrift=posdrift,robust=robust))
  #   S1b <- (1-plba.norm(dt[1,is.dk],A=A[1,is.dk],b=b[1,is.dk],
  #       mean_v=mean_v[1,is.dk],sd_v=sd_v[1,is.dk],posdrift=posdrift,robust=robust))
  #   S2b <- (1-plba.norm(dt[2,is.dk],A=A[2,is.dk],b=b[2,is.dk],
  #       mean_v=mean_v[2,is.dk],sd_v=sd_v[2,is.dk],posdrift=posdrift,robust=robust))
  #   S2d <- (1-plba.norm(dt[2,is.dk],A=A[2,is.dk],b=d[2,is.dk],
  #       mean_v=mean_v[2,is.dk],sd_v=sd_v[2,is.dk],posdrift=posdrift,robust=robust))
  #   
  #   if ( any( is.dk) ) dt[1, is.dk] <- D1b*(S2b-S2d) + D2b*(S1b-S1d)
  
  
  if ( any(!is.dk) ) dt[1,!is.dk] <- 
    dlba.norm(dt[1,!is.dk],A=A[1,!is.dk],b=b[1,!is.dk],mean_v=mean_v[1,!is.dk],
              sd_v=sd_v[1,!is.dk],posdrift=posdrift,robust=robust)*
    (1-plba.norm(dt[2,!is.dk],A=A[2,!is.dk],b=d[2,!is.dk],
                 mean_v=mean_v[2,!is.dk],sd_v=sd_v[2,!is.dk],posdrift=posdrift,robust=robust))
  if ( any(is.dk) ) dt[1,is.dk] <- 
    dlba.norm(dt[1,is.dk],A=A[1,is.dk],b=b[1,is.dk],
              mean_v=mean_v[1,is.dk],sd_v=sd_v[1,is.dk],posdrift=posdrift,robust=robust)*
    ((1-plba.norm(dt[2,is.dk],A=A[2,is.dk],b=b[2,is.dk],
                  mean_v=mean_v[2,is.dk],sd_v=sd_v[2,is.dk],posdrift=posdrift,robust=robust))-
       (1-plba.norm(dt[2,is.dk],A=A[2,is.dk],b=d[2,is.dk],
                    mean_v=mean_v[2,is.dk],sd_v=sd_v[2,is.dk],posdrift=posdrift,robust=robust))) +
    dlba.norm(dt[2,is.dk],A=A[2,is.dk],b=b[2,is.dk],
              mean_v=mean_v[2,is.dk],sd_v=sd_v[2,is.dk],posdrift=posdrift,robust=robust)*
    ((1-plba.norm(dt[1,is.dk],A=A[1,is.dk],b=b[1,is.dk],
                  mean_v=mean_v[1,is.dk],sd_v=sd_v[1,is.dk],posdrift=posdrift,robust=robust))-
       (1-plba.norm(dt[1,is.dk],A=A[1,is.dk],b=d[1,is.dk],
                    mean_v=mean_v[1,is.dk],sd_v=sd_v[1,is.dk],posdrift=posdrift,robust=robust)))
  
  dt[1,]
}

# ---- n1PDF.DK (source lines 255-276) ----
n1PDF.DK <- function (t,r2,A,b,d,mean_v,sd_v,t0,st0=0,posdrift=TRUE,robust = FALSE)
  # r2 is second respond, d is dont know crierion
{
  dt <- pmax(t-t0[1], 0)
  if (st0[1] == 0) 
    return(n1PDFfixedt0.DK(dt,r2=r2,A=A,b=b,d=d,mean_v=mean_v,sd_v=sd_v, 
                           posdrift=posdrift,robust = robust))
  else {
    tmpf <- function(tau,ti0,r2,A,b,d,mean_v,sd_v,scale,shape) 
      n1PDFfixedt0.DK(tau,r2,A,b,mean_v,sd_v)/st0
    outs <- numeric(length(dt))
    for (i in c(1:length(outs))) {
      tmp <- try(integrate(f=tmpf,
                           lower=0,upper=dt[i],ti0=dt[i],r2=r2[i],
                           A=A,b=b,d=d,mean_v=mean_v,sd_v=sd_v)$value,silent=T)
      if (is.numeric(tmp)) 
        outs[i] <- tmp else
          outs[i] <- 0
    }
    return(outs)
  }
}

