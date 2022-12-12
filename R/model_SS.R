# Used by different stop-signal models

update_ssd <- function(isstop,idx,idx1,ssd,stairstep,stairmin,stairmax)
  # Used in random function
{
  if (isstop) {
    if (ssd[idx]+ stairstep < stairmax)
      ssd[idx1] <- ssd[idx] + stairstep else ssd[idx1] <- ssd[idx]
    } else {
      if (ssd[idx] - stairstep > stairmin)
          ssd[idx1] <- ssd[idx] - stairstep else ssd[idx1] <- ssd[idx]
    }
  ssd
}


# p stop functions

my.integrate <- function(...,big=10)
  # Avoids bug in integrate upper=Inf that uses only 1  subdivision
  # Use of  big=10 is arbitrary ...
{
  out <- try(integrate(...,upper=Inf),silent=TRUE)
  if (is(out,"try-error")) 0 else
  {
    if (out$subdivisions==1)
    {
      out <- try(integrate(...,upper=big),silent=TRUE)
      if (is(out,"try-error")) 0 else
      {
        if (out$subdivisions==1) 0 else out$value
      }
    } else out$value
  }
}
