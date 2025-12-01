# Last observation carried forward
# Replaces NA values with the last non-NA value
# Mimics zoo::na.locf behavior (vectorized for speed)
na_locf <- function(x, na.rm = FALSE) {
  if (length(x) == 0) return(x)
  
  # Check if all values are NA
  na_mask <- is.na(x)
  if (all(na_mask)) {
    if (na.rm) {
      return(x[0])  # Return empty vector
    } else {
      return(x)  # Return as is
    }
  }
  
  # Vectorized approach: create an index that tracks the last non-NA position
  # For each position, we need the index of the last non-NA value up to that point
  idx <- seq_along(x)
  idx[na_mask] <- NA  # Set NA positions to NA in index
  idx <- cummax(ifelse(na_mask, 0, idx))  # Cumulative max gives us last non-NA index
  
  # Replace values: use the index to look up the last non-NA value
  # For positions where idx is 0 (leading NAs), keep as NA
  result <- x
  non_zero <- idx > 0
  result[non_zero] <- x[idx[non_zero]]
  
  # If na.rm=TRUE, remove leading NAs
  if (na.rm) {
    first_non_na <- which(!na_mask)[1]
    if (!is.na(first_non_na)) {
      result <- result[first_non_na:length(result)]
    }
  }
  
  return(result)
}

#
# augment <- function(s,da,design)
#   # Adds attributes to augmented data
#   # learn: empty array for Q values with dim = choice alternative (low,high) x
#   #   stimulus x trials (max across stimuli)
#   # index: look up (row number) for stimuli in da, a matrix dim  = max trials x
#   #   stimulus matrix (rows for each choice alternative contiguous)
# {
#   if (!is.null(design$adapt$stimulus)) {
#     targets <- design$adapt$stimulus$targets
#     par <- design$adapt$stimulus$output_name
#     maxn <- max(sapply(dimnames(targets)[[1]],function(x){table(da[da$subjects==s,x])}))
#     # da index x stimulus
#     out <- sapply(targets[1,],getIndex,cname=dimnames(targets)[[1]][1],
#                   da=da[da$subjects==s,],maxn=maxn)
#     stimulus <- list(index=out)
#     # accumulator x stimulus x trials
#     stimulus$learn <- array(NA,dim=c(dim(targets),maxn/dim(targets)[1]),
#                             dimnames=list(rownames(targets),targets[1,],NULL))
#     stimulus$targets <- targets
#     stimulus$par <- par
#   } # add other types here
#   list(stimulus=stimulus)
# }
#
# getIndex <- function(typei,cname,da,maxn) {
#   out <- which(da[,cname]==typei)
#   c(out,rep(NA,maxn-length(out)))
# }
