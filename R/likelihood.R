log_likelihood_race <- function(p_vector,dadm,min_ll=log(1e-10))
  # Race model summed log likelihood
{

  pars <- get_pars(p_vector,dadm)
  if (is.null(attr(pars,"ok")))
    ok <- !logical(dim(pars)[1]) else ok <- attr(pars,"ok")

    lds <- numeric(dim(dadm)[1]) # log pdf (winner) or survivor (losers)
    lds[dadm$winner] <- log(attr(dadm,"model")()$dfun(rt=dadm$rt[dadm$winner],
                                                    pars=pars[dadm$winner,]))
    n_acc <- length(levels(dadm$R))
    if (n_acc>1) lds[!dadm$winner] <- log(1-attr(dadm,"model")()$pfun(rt=dadm$rt[!dadm$winner],pars=pars[!dadm$winner,]))
    lds[is.na(lds) | !ok] <- min_ll
    lds <- lds[attr(dadm,"expand")] # decompress
    if (n_acc>1) {
      winner <- dadm$winner[attr(dadm,"expand")]
      ll <- lds[winner]
      if (n_acc==2) {
        ll <- ll + lds[!winner]
      } else {
        ll <- ll + apply(matrix(lds[!winner],nrow=n_acc-1),2,sum)
      }
      ll[is.na(ll)] <- min_ll
      return(sum(pmax(min_ll,ll)))
    } else return(sum(pmax(min_ll,lds)))
}



log_likelihood_ddm <- function(p_vector,dadm,min_ll=log(1e-10))
  # DDM summed log likelihood, with protection against numerical issues
{
  pars <- get_pars(p_vector,dadm)
  like <- numeric(dim(dadm)[1])
  if (any(attr(pars,"ok"))) like[attr(pars,"ok")] <-
    attr(dadm,"model")()$dfun(dadm$rt[attr(pars,"ok")],dadm$R[attr(pars,"ok")],pars[attr(pars,"ok"),,drop=FALSE])
  like[attr(pars,"ok")][is.na(like[attr(pars,"ok")])] <- 0
  sum(pmax(min_ll,log(like[attr(dadm,"expand")])))
}

#### sdt choice likelihoods ----

log_likelihood_sdt <- function(p_vector,dadm,lb=-Inf,min_ll=log(1e-10))
  # probability of ordered discrete choices based on integrals of a continuous
  # distribution between thresholds, with fixed lower bound for first response
  # lb. Upper bound for last response is a fixed value in threshold vector
{

  pars <- get_pars(p_vector,dadm)
  first <- dadm$lR==levels(dadm$lR)[1]
  last <- dadm$lR==levels(dadm$lR)[length(levels(dadm$lR))]
  pars[last,"threshold"] <- Inf
  # upper threshold
  ut <- pars[dadm$winner,"threshold"]
  # lower threshold fixed at lb for first response
  pars[first &  dadm$winner,"threshold"] <- lb
  # otherwise threshold of response before one made
  notfirst <- !first &  dadm$winner
  pars[notfirst,"threshold"] <- pars[which(notfirst)-1,"threshold"]
  lt <- pars[dadm$winner,"threshold"]
  # fix race expand to suit SDT
  nr <- length(levels(dadm$lR))      # number of responses
  ne <- length(attr(dadm,"expand"))  # length of expand
  # Shorten expand to only one per lR set
  expand <- (attr(dadm,"expand")[(c(1:ne) %% nr)==0] + 1 ) %/% nr
  # log probability
  ll <- numeric(sum(dadm$winner))
  if (!is.null(attr(pars,"ok"))) { # Bad parameter region
    ok <- attr(pars,"ok")
    okw <- ok[dadm$winner]
    ll[ok] <- log(attr(dadm,"model")()$pfun(lt=lt[okw],ut=ut[okw],pars=pars[dadm$winner & ok,]))
  } else ll <- log(attr(dadm,"model")()$pfun(lt=lt,ut=ut,pars=pars[dadm$winner,]))
  ll <- ll[expand]
  ll[is.na(ll)] <- 0
  sum(pmax(min_ll,ll))
}

log_likelihood_joint <- function(pars, dadms, component = NULL)
{
  parPreFixs <- unique(gsub("[|].*", "", names(pars)))
  i <- 0
  total_ll <- 0
  if(!is.null(component)) dadms <- dadms[component]
  for(dadm in dadms){
    if(is.data.frame(dadm)){
      i <- i + 1
      parPrefix <- parPreFixs[i]
      currentPars <- pars[grep(paste0(parPrefix, "|"), names(pars), fixed = T)]
      names(currentPars) <- gsub(".*[|]", "", names(currentPars))
      total_ll <- total_ll +  attr(dadm, "model")$log_likelihood(currentPars, dadm)
    }
  }
  return(total_ll)
}


