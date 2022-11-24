
pPROBIT <- function(lt,ut,pars)
  # probability between lt and ut
{
  pnorm(ut,mean=pars[,"mean"],sd=pars[,"sd"]) - pnorm(lt,mean=pars[,"mean"],sd=pars[,"sd"])
}

#### random


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

# Normal, natural mu and log(sigma), increasing threshold (first natural,
# others on log scale) parameterization
probit <- function(){
  list(
  type="SDT", # Discrete choice based on continuous latent, no RT
  p_types=c("mean","sd","threshold"),
  # Transform to natural scale
  Ntransform=function(x) {
    is_sd <- grepl("sd",dimnames(x)[[2]])
    x[,is_sd] <- exp(x[,is_sd])
    x
  },
  # p_vector transform
  transform = function(x) {
    if (!is.matrix(x)) {
      increasing <- grepl("threshold",names(x)) & grepl(":lR",names(x)) | grepl("threshold_lR",names(x))
      x[increasing] <- exp(x[increasing])
      x
    } else {
      increasing <- grepl("threshold",dimnames(x)[[2]]) & grepl(":lR",dimnames(x)[[2]]) |
        grepl("threshold_lR",dimnames(x)[[2]])
      x[,increasing] <- exp(x[,increasing])
      x
    }
  },
  # Trial dependent parameter transform
  Ttransform = function(pars,dadm) {
    pars
  },
  # Random function for discrete choices
  rfun=function(lR,pars) rPROBIT(lR,pars),
  # probability of choice between lower and upper thresholds (lt & ut)
  pfun=function(lt,ut,pars) pPROBIT(lt,ut,pars),
  # quantile function, p = probability, used in making linear ROCs
  qfun=function(p) qnorm(p),
  # Likelihood, lb is lower bound threshold for first response
  log_likelihood=function(p_vector,dadm,min_ll=log(1e-10)){
    log_likelihood_sdt(p_vector=p_vector, dadm = dadm, min_ll = min_ll, lb=-Inf)
  })
}
