log_likelihood_race <- function(p_vector,dadm,min_ll=log(1e-10))
  # Race model summed log likelihood
{

  pars <- get_pars(p_vector,dadm)
  if (is.null(attr(pars,"ok")))
    ok <- !logical(dim(pars)[1]) else ok <- attr(pars,"ok")

    lds <- numeric(dim(dadm)[1]) # log pdf (winner) or survivor (losers)
    lds[dadm$winner] <- log(attr(dadm,"model")$dfun(rt=dadm$rt[dadm$winner],
                                                    pars=pars[dadm$winner,]))
    n_acc <- length(levels(dadm$R))
    if (n_acc>1) lds[!dadm$winner] <-
      log(1-attr(dadm,"model")$pfun(rt=dadm$rt[!dadm$winner],pars=pars[!dadm$winner,]))
    lds[is.na(lds) | !ok] <- 0
    lds <- lds[attr(dadm,"expand")] # decompress
    if (n_acc>1) {
      winner <- dadm$winner[attr(dadm,"expand")]
      ll <- lds[winner]
      if (n_acc==2)
        ll <- ll + lds[!winner] else
          ll <- ll + apply(matrix(lds[!winner],nrow=n_acc-1),2,sum)
      ll[is.na(ll)] <- 0
      sum(pmax(min_ll,ll))
    } else sum(pmax(min_ll,lds))
}
