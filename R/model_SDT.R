
pPROBIT <- function(lt,ut,pars)
  # probability between lt and ut
{
  pnorm(ut,mean=pars[,"mean"],sd=pars[,"sd"]) - pnorm(lt,mean=pars[,"mean"],sd=pars[,"sd"])
}


rPROBIT <- function(lR,pars,p_types=c("mean","sd","threshold"),lt=-Inf)
  # lR is an empty latent response factor lR with one level for response.
  # pars is a matrix of corresponding parameter values named as in p_types
  # pars must be sorted so accumulators and parameter for each trial are in
  # contiguous rows.

{
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ",paste(p_types,collapse = " "))
  nr <- length(levels(lR)) # Number of responses
  n <- dim(pars)[1]/nr     # Number of simulated trials
  first <- seq(1,dim(pars)[1]-nr+1,length.out=n)  # pick out mean and sd
  threshold <- matrix(pars[,"threshold"],nrow=nr) # format thresholds
  threshold[dim(threshold)[1],] <- Inf
  pmat <- rbind(rnorm(n,pars[first,"mean"],pars[first,"sd"]), # sample normal
                rep(lt,dim(threshold)[2]),threshold)          # lt ... ut
  pmat[dim(pmat)[1],] <- Inf
  R <- factor(apply(pmat,2,function(x){.bincode(x[1],x[-1])}),
              levels=1:length(levels(lR)),labels=levels(lR))
  cbind.data.frame(R=R,rt=NA)
}

# #' Gaussian Signal Detection Theory Model
# #'
# #' Discrete choice based on continuous Gaussian latent, with no rt. Model
# #' parameters are mean (unbounded) sd (log scale) and threshold, with a first
# #' value is  on the natural scale, and others for designs with with more than
# #' two responses are threshold increases on a log scale to enforce monotonic
# #' increase on the natural scale.
# #'
# #' @return A model list with all the necessary functions to sample
#' @export
#'

probit <- function(){
  list(
  type="SDT",
  p_types=c("mean" = 0,"sd" = log(1),"threshold" = 0),
  # Trial dependent parameter transform
  transform=list(func=c(mean = "identity",sd = "exp",threshold="identity")),
  bound=list(minmax=cbind(mean=c(-Inf,Inf),sd = c(0, Inf), threshold=c(-Inf,Inf))),
  Ttransform = function(pars,dadm) {
    pars
  },
  # Random function for discrete choices
  rfun=function(lR=NULL,pars) {
    rPROBIT(lR,pars)
  },
  # probability of choice between lower and upper thresholds (lt & ut)
  pfun=function(lt,ut,pars) pPROBIT(lt,ut,pars),
  # quantile function, p = probability, used in making linear ROCs
  qfun=function(p) qnorm(p),
  # Likelihood, lb is lower bound threshold for first response
  log_likelihood=function(pars,dadm,model,min_ll=log(1e-10)){
    log_likelihood_sdt(pars=pars, dadm = dadm, model = model, min_ll = min_ll, lb=-Inf)
  })
}
