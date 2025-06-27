calc_ll_R <- function(p_vector, model, dadm){
  if(!is.null(model$transform)){
    pars <- get_pars_matrix(p_vector, dadm, model)
  } else{
    pars <- p_vector
  }
  ll <- model$log_likelihood(pars, dadm, model)
  return(ll)
}


log_likelihood_race <- function(pars,dadm,model,min_ll=log(1e-10))
  # Race model summed log likelihood
{
  if (any(names(dadm)=="RACE")){# Some accumulators not present
    pars[as.numeric(dadm$lR)>as.numeric(as.character(dadm$RACE)),] <- NA
  }

  if (is.null(attr(pars,"ok"))){
    ok <- !logical(dim(pars)[1])
  } else ok <- attr(pars,"ok")
  lds <- numeric(dim(dadm)[1]) # log pdf (winner) or survivor (losers)
  lds[dadm$winner] <- log(model$dfun(rt=dadm$rt[dadm$winner],
                                                    pars=pars[dadm$winner,]))
  n_acc <- length(levels(dadm$R))
  if (n_acc>1) lds[!dadm$winner] <- log(1-model$pfun(rt=dadm$rt[!dadm$winner],pars=pars[!dadm$winner,]))
  lds[is.na(lds) | !ok] <- min_ll
  if (n_acc>1) {
    ll <- lds[dadm$winner]
    if (n_acc==2) {
      ll <- ll + lds[!dadm$winner]
    } else {
      ll <- ll + apply(matrix(lds[!dadm$winner],nrow=n_acc-1),2,sum)
    }
    ll[is.na(ll)] <- min_ll
    return(sum(pmax(min_ll,ll[attr(dadm,"expand")])))
  } else return(sum(pmax(min_ll,lds[attr(dadm,"expand")])))
}



log_likelihood_ddm <- function(pars,dadm,model,min_ll=log(1e-10))
  # DDM summed log likelihood, with protection against numerical issues
{
  like <- numeric(dim(dadm)[1])
  if (any(attr(pars,"ok")))
    like[attr(pars,"ok")] <- model$dfun(dadm$rt[attr(pars,"ok")],dadm$R[attr(pars,"ok")],
                                                       pars[attr(pars,"ok"),,drop=FALSE])
  like[attr(pars,"ok")][is.na(like[attr(pars,"ok")])] <- 0
  sum(pmax(min_ll,log(like[attr(dadm,"expand")])))
}


log_likelihood_ddmgng <- function(pars,dadm,model,min_ll=log(1e-10))
  # DDM summed log likelihood for go/nogo model
{
  like <- numeric(dim(dadm)[1])
  if (any(attr(pars,"ok"))) {
    isna <- is.na(dadm$rt)
    ok <- attr(pars,"ok") & !isna
    like[ok] <- model$dfun(dadm$rt[ok],dadm$R[ok],pars[ok,,drop=FALSE])
    ok <- attr(pars,"ok") & isna
    like[ok] <- # dont terminate on go boundary before timeout
      pmax(0,pmin(1,(1-model$pfun(dadm$TIMEOUT[ok],dadm$Rgo[ok],pars[ok,,drop=FALSE]))))

  }
  like[attr(pars,"ok")][is.na(like[attr(pars,"ok")])] <- 0
  sum(pmax(min_ll,log(like[attr(dadm,"expand")])))
}



#### sdt choice likelihoods ----

log_likelihood_sdt <- function(pars,dadm, model,lb=-Inf, min_ll=log(1e-10))
  # probability of ordered discrete choices based on integrals of a continuous
  # distribution between thresholds, with fixed lower bound for first response
  # lb. Upper bound for last response is a fixed value in threshold vector
{
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
  # log probability
  ll <- numeric(sum(dadm$winner))
  if (!is.null(attr(pars,"ok"))) { # Bad parameter region
    ok <- attr(pars,"ok")
    okw <- ok[dadm$winner]
    ll[ok] <- log(model$pfun(lt=lt[okw],ut=ut[okw],pars=pars[dadm$winner & ok,,drop=FALSE]))
  } else ll <- log(model$pfun(lt=lt,ut=ut,pars=pars[dadm$winner,,drop=FALSE]))
  ll <- ll[attr(dadm,"expand")]
  ll[is.na(ll)] <- 0
  sum(pmax(min_ll,ll))
}

# Two options:
# 1) component = NULL in which case we do all likelihoods in one block
# 2) component = integer, in which case we are blocking the ll and only want that one
log_likelihood_joint <- function(proposals, dadms, model_list, component = NULL){
  parPreFixs <- unique(gsub("[|].*", "", colnames(proposals)))
  i <- 0
  k <- 0
  total_ll <- 0
  for(dadm in dadms){
    i <- i + 1
    if(is.null(component) || component == i){
      if(is.data.frame(dadm)){
        k <- k + 1
        # Sometimes designs are the same across models (typically in fMRI)
        # Instead of storing multiple, we just store a pointer to the first original design
        if(is.numeric(attr(dadm, "designs"))){
          ref_idx <- attr(dadm, "designs")
          attr(dadm, "designs") <- attr(dadms[[ref_idx]], "designs")
        }
        parPrefix <- parPreFixs[k]
        # Unfortunately indexing can't get quicker than this as far as I can tell.
        columns_to_use <- sapply(strsplit(colnames(proposals), "|", fixed = TRUE), function(x) x == parPrefix)[1,]
        currentPars <- proposals[,columns_to_use, drop = F]
        colnames(currentPars) <- gsub(".*[|]", "", colnames(currentPars))
        total_ll <- total_ll +  calc_ll_manager(currentPars, dadm, model_list[[i]])
      }
    }
  }
  return(total_ll)
}
