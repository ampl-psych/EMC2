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
