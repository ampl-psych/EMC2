## ALL FUNCTIONS OVERWRITTEN BY SM

augment = function(s,da,design)
  # Adds attributes to augmented data
  # learn: empty array for Q values with dim = choice alternative (low,high) x
  #   stimulus x trials (max across stimuli)
  # index: look up (row number) for stimuli in da, a matrix dim  = max trials x
  #   stimulus matrix (rows for each choice alternative contiguous)
{
  getIndex <- function(typei,cname,da,maxn) {
    out <- which(da[,cname]==typei)
    c(out,rep(NA,maxn-length(out)))
  }

  getIndexThis <- function(da, outcomes) {
    ## gets index of Q-value corresponding to the accumulator
    # NB: outcomes is sorted by trial order, da is *not*
    # so we need to return an index that takes this into account
    #da$lS <- da[cbind(1:nrow(da), match(as.character(da$lR), colnames(da)))]
    da$colNumbers <- match(as.character(da$lS), colnames(outcomes))

    return(cbind(da$trials, da$colNumbers))
  }

  getIndexOther <- function(da, outcomes) {
    ## ONLY WORKS FOR 2AFC!!
    ## gets index of Q-value corresponding to the *OTHER* accumulator (in an advantage framework)
    lROptions <- unique(as.character(da$lR))
    da$lRother <- ifelse(da$lR==lROptions[1], lROptions[2], lROptions[1])
    da$lSOther <- da[cbind(1:nrow(da), match(paste0('s_', as.character(da$lRother)), colnames(da)))]
    da$colNumbers <- match(as.character(da$lSOther), colnames(outcomes))
    return(cbind(da$trials, da$colNumbers))
    # return(cbind(da$trials, match(da[cbind(da$trials, match(da$lRother, colnames(da)))], colnames(outcomes))))
  }

  makepArray <- function(x) {
    stim <- design$adapt$stimulus$targets
    x <- x[x$R==x$lR,]          # NB: this ONLY allows the winning accumulator/stimulus to receive feedback!!
    x <- x[order(x$trials),]    # *must* be ordered by trial n for this matrix

    pArray <- matrix(NA, nrow=max(x$trials), ncol=length(stim))
    colnames(pArray) <- stim
    for(i in 1:nrow(x)) {
      trial <- x[i,'trials']
      pArray[trial,as.character(x[i,'s_left'])] <- x[i,'p_left']
      pArray[trial,as.character(x[i,'s_right'])] <- x[i,'p_right']
    }
    return(pArray)
  }

  makeDCT <- function(frame_times, high_pass=1/128) {
      n_frames = length(frame_times)
      n_times <- 1:n_frames

      dt = (frame_times[length(frame_times)]-frame_times[1]) / (n_frames-1)  # should be 1 with trials
      order = pmin(n_frames-1, floor(2*n_frames*high_pass*dt))
      cosine_drift = matrix(0, nrow=n_frames, ncol=order+1)
      normalizer = sqrt(2/n_frames)

      for(k in seq(2, order+2)) {
        cosine_drift[,k-1] = normalizer*cos((pi/n_frames)*(n_times+0.5)*k)
      }
      return(cosine_drift)
  }

  makeOutcomes <- function(x) {
    stim <- design$adapt$stimulus$targets
    x <- x[x$R==x$lR,]          # NB: this ONLY allows the winning accumulator/stimulus to receive feedback!!
    # x <- x[x$winner,]
    x <- x[order(x$trials),]    # *must* be ordered by trial n for this matrix

    outcomes <- data.frame(matrix(NA, nrow=max(x$trials), ncol=length(stim)))
    colnames(outcomes) <- stim
    for(trial in unique(x$trials)) {
      outcomes[trial,as.character(x[trial,'lS'])] <- x[trial,'reward']
    }
    return(outcomes)
  }

  if (!is.null(design$adapt$stimulus)) {
    targets <- design$adapt$stimulus$targets
    par <- design$adapt$stimulus$output_name
    #    maxn <- max(sapply(dimnames(targets)[[1]], function(x){table(da[da$subjects==s,x])}))

    # da index x stimulus
    # out <- sapply(targets[1,],getIndex,cname=dimnames(targets)[[1]][1],
    #               da=da[da$subjects==s,],maxn=maxn)
    stimulus <- list() #index=out)
    # accumulator x stimulus x trials
    ## SM: why not stimulus x trials, and match to accumulators at a later point (via getIndex)? would make using c much easier
    # stimulus$learn <- array(NA,dim=c(dim(targets),maxn/dim(targets)[1]),
    #                          dimnames=list(rownames(targets),targets[1,],NULL))
    stimulus$targets <- targets
    stimulus$par <- par

    if('trials' %in% colnames(da)) {
      outcomes <- makeOutcomes(da[da$subjects==s,])
      pArray <- makepArray(da[da$subjects==s,])
      return(list(stimulus=stimulus, outcomes=outcomes, pArray=pArray,
                  index=getIndexThis(da[da$subjects==s,], outcomes),
                  indexOther=getIndexOther(da[da$subjects==s,], outcomes)))
    }
  } else if(!is.null(design$adapt$dynamic)) { # add other types here
    # outcomes becomes a matrix of nTrials x [accuracy, choice]
    # no pArray
    columns_to_include <- design$adapt$dynamic$columns_to_include
    if('trials' %in% colnames(da)) {

      tmp <- da[da$subjects==s,]
      tmp <- tmp[order(tmp$trials),]
      outcomes <- tmp[, columns_to_include]
      if(length(columns_to_include) == 1) outcomes <- matrix(outcomes, ncol=1)
      index <- match(da[da$subjects==s,'trials'], tmp$trials)

      if(is.null(design$adapt$includeDCT)) {
        return(list(outcomes=outcomes, index=index))
      } else {
        cdt <- makeDCT(high_pass=1/64, frame_times=1:nrow(tmp))
        return(list(outcomes=outcomes, index=index, cdt=cdt))  ##
      }

    }
  } else {
    return(list(stimulus=stimulus))
  }
}

adapt.c.emc <- function(...){
  return(NA)
}

# ARCHITECTURE OF FOLLOWING FUNCTION WILL NEED UPDATING FOR MULTIPLE ADAPT TYPES
# s="10";npars=pars;da=data; rfun=model$rfun;return_learning=FALSE;mapped_p=FALSE
update_pars = function(s,npars,da,rfun=NULL,return_learning=FALSE,mapped_p=FALSE,
                       return_all=FALSE)
  # for subject s
  # Either return da filling in responses and RT if da has no responses (in
  # in which case rfun must be supplied), or return npars filling in adapted
  # parameters or if return_learning returns learning (e.g., Q values)
{
  adapt <- attr(da,"adapt")$design

  if(adapt$useC & !mapped_p & !return_all) {
    outcomes <- attr(da, 'adapt')[[s]]$outcomes
    index <- attr(da, 'adapt')[[s]]$index
    npars <- npars[da$subjects==s,]
    da <- da[da$subjects==s,]

    # update
    if(!is.null(adapt$stimulus) & !any(is.na(da$R))) {
      if(adapt$stimulus$adapt_fun_name=='delta') {
        learningRates <- matrix(npars[1,'alpha'], nrow=nrow(outcomes), ncol=ncol(outcomes))  # TODO make this more flexible
        startValues <- rep(npars[1,'q0'], ncol(outcomes))  # TODO make this more flexible
        updated <- adapt.c.emc(feedback=as.matrix(outcomes),
                               arguments=list(startValues = startValues,
                                              learningRates = learningRates),
                               learningRule='delta')
        # updated <- adapt.c.dmc(startValues = startValues,
        #                        learningRates = learningRates,
        #                        feedback = as.matrix(outcomes),
        #                        learningRule='SARSA')
        colnames(updated$adaptedValues) <- colnames(outcomes)
        allQs <- updated$adaptedValues
      } else if(adapt$stimulus$adapt_fun_name=='vkf') {
        volatilityLearningRates <- matrix(npars[1,'alpha'], nrow=nrow(outcomes), ncol=ncol(outcomes))  # TODO make this more flexible
        predictionsStartValues <- rep(npars[1,'q0'], ncol(outcomes))  # TODO make this more flexible
        volatilitiesStartValues <- rep(npars[1,'volatility0'], ncol(outcomes))  # TODO make this more flexible
        uncertaintiesStartValues <- rep(npars[1,'w0'], ncol(outcomes))
        updated <- adapt.c.emc(feedback=as.matrix(outcomes),
                               arguments=list(volatilityLearningRates = volatilityLearningRates,
                                              predictionsStartValues = predictionsStartValues,
                                              volatilitiesStartValues = volatilitiesStartValues,
                                              uncertaintiesStartValues = uncertaintiesStartValues),
                               learningRule='vkf')
        allQs <- updated$adaptedPredictions
      } else if(adapt$stimulus$adapt_fun_name=='vkfbinary') {
        volatilityLearningRates <- matrix(npars[1,'alpha'], nrow=nrow(outcomes), ncol=ncol(outcomes))  # TODO make this more flexible
        predictionsStartValues <- rep(npars[1,'q0'], ncol(outcomes))  # TODO make this more flexible
        volatilitiesStartValues <- rep(npars[1,'volatility0'], ncol(outcomes))  # TODO make this more flexible
        uncertaintiesStartValues <- rep(npars[1,'w0'], ncol(outcomes))
        updated <- adapt.c.emc(feedback=as.matrix(outcomes),
                               arguments=list(volatilityLearningRates = volatilityLearningRates,
                                              predictionsStartValues = predictionsStartValues,
                                              volatilitiesStartValues = volatilitiesStartValues,
                                              uncertaintiesStartValues = uncertaintiesStartValues),
                               learningRule='vkf')
        allQs <- updated$adaptedPredictions
      }
      if(return_learning) {
        return(updated)
      }

      Q <- allQs[index]

      # Advantage framework (2AFC ONLY!!)
      if(('ws' %in% adapt$stimulus$output_par_names) & ('wd' %in% adapt$stimulus$output_par_names)) {
        indexOther <- attr(da, 'adapt')[[s]]$indexOther
        Q <- cbind(Q, allQs[indexOther])
      }

      ## function
      npars[,adapt$stimulus$output_name] <- adapt$stimulus$output_fun(npars[,adapt$stimulus$output_par_names], Q)
      ## hacky way of preventing errors due to extreme values
      npars[,'v'] <- pmin(npars[,'v'], 1e3)
      npars[,'v'] <- pmax(npars[,'v'], 0)

      ## add prediction errors and other latent learning variables!
      attr(npars, 'learn') <- updated
      return(npars)

    } else if(!is.null(adapt$dynamic)) {

      ###### DYNAMIC
      if(adapt$dynamic$adapt_fun_name=='delta') {

        learningRates <- as.matrix(cbind(npars[,'alpha1'], npars[,'alpha2'], npars[,'alpha3']))  #matrix(npars[1,'alpha'], nrow=nrow(outcomes), ncol=ncol(outcomes))  # TODO make this more flexible
        startValues <- npars[1, c('q01', 'q02', 'q03')]

        # remove not-updated columns
        alphaIsZero <- round(learningRates[1,],6)==0
        learningRates <- learningRates[,!alphaIsZero]
        startValues <- startValues[!alphaIsZero]

        if(!all(alphaIsZero)) {
          updated <- adapt.c.emc(feedback=as.matrix(outcomes),
                                 arguments=list(startValues = startValues,
                                                learningRates = learningRates),
                                 learningRule='delta')
        } else {
          #### Nothing to update, return npars!
          return(npars)
        }

        allQs <- matrix(NA, nrow=nrow(outcomes), ncol=3)
        if(sum(alphaIsZero)>0) allQs[,alphaIsZero] <- npars[, c('q01', 'q02', 'q03')][,alphaIsZero]    # re-add not-updated columns
        if(!all(alphaIsZero)) allQs[,!alphaIsZero] <- updated$adaptedValues
      }

      # if(!is.null(adapt$dynamic$output_name)) {
        ##### FOR RDM
        npars[,adapt$dynamic$output_name] <- adapt$dynamic$output_fun(npars, allQs[index,], da)
      # } else {
        ###### FOR DDM STILL HARDCODED!
        # npars[,'a'] <- npars[,'a'] + npars[,'weight1']*(1-allQs[index,1])
        # npars[,'a'] <- pmin(npars[,'a'], 10)  # hacky hacky
        # npars[,'a'] <- pmax(npars[,'a'], 0.1)

        # npars[,'v'] <- npars[,'v'] + npars[,'weight2']*(allQs[index,2])

        # npars[,'Z'] <- npars[,'Z'] + (npars[,'weight3']-.5)*(allQs[index,3])
        # npars[,'Z'] <- pmin(npars[,'Z'], 0.9)  # hacky hacky
        # npars[,'Z'] <- pmax(npars[,'Z'], 0.1)

        # ## SM: this shouldn't be here, move this to adapt_data later!!
        # if(any(is.na(da$R))) {
        #   pars <- npars
        #   pars <- cbind(pars,z=pars[,"a"]*pars[,"Z"],
        #                 sz = 2*pars[,"SZ"]*pars[,"a"]*apply(cbind(pars[,"Z"],1-pars[,"Z"]),1,min))
        #   pars <- cbind(pars, d = pars[,"t0"]*(2*pars[,"DP"]-1))
        #
        #   attr(pars,"ok") <-
        #     !( abs(pars[,"v"])> 20 | pars[,"a"]> 10 | pars[,"sv"]> 10 | pars[,"SZ"]> .999 | pars[,"st0"]>.2)
        #   if (pars[1,"sv"] !=0) attr(pars,"ok") <- attr(pars,"ok") & pars[,"sv"] > .001
        #   if (pars[1,"SZ"] !=0) attr(pars,"ok") <- attr(pars,"ok") & pars[,"SZ"] > .001
        #   npars <- pars
        # }
       # }

#       if(any(is.na(da$R))) {
#         if(is.null(rfun)) rfun <- attr(da, 'model')()$rfun
#         Rrt <- rfun(da$lR, npars)
# #        da <- da[seq(1, nrow(da), 2),]  # unique trials
#         da[,c('R', 'rt')] <- Rrt[rep(c(1:nrow(Rrt)), each=nrow(da)/nrow(Rrt)),]
#         return(da)
#       }
      return(npars)
    }
  } else if(adapt$useSMsApproach) {
    outcomes <- attr(da, 'adapt')[[s]]$outcomes
    index <- attr(da, 'adapt')[[s]]$index
    npars <- npars[da$subjects==s,]
    da <- da[da$subjects==s,]

    add_response <- any(is.na(da$R))
    if(add_response) {
      if(is.null(rfun)) rfun <- attr(da, 'model')()$rfun
      ## extract probability of winning per trial & stimulus
      pArray <- attr(da, 'adapt')[[s]]$pArray

      # reset
      outcomes <- matrix(NA, ncol=ncol(outcomes), nrow=nrow(outcomes), dimnames=dimnames(outcomes))
    }

    # todo: make this a 3dimensional array, of nTrials x nChoices x learningThing
    # with 'learningThing' being predictions(/Q), uncertainty, volatility; possibly also PEs?
    learn <- array(NA, dim=c(nrow(outcomes), ncol(outcomes), length(adapt$stimulus$init_par)))
    dimnames(learn) <- list(NULL, colnames(outcomes), adapt$stimulus$init_par)

    # learn <- matrix(NA, nrow=nrow(outcomes), ncol=ncol(outcomes), dimnames = dimnames(outcomes))
    if(length(npars[1,adapt$stimulus$init_par]) > 1) {
      learn[1,,] <- matrix(npars[1,adapt$stimulus$init_par],
                           ncol=length(npars[1,adapt$stimulus$init_par]),
                           nrow=nrow(learn[1,,]), byrow=TRUE) #  # initialize
    } else {
      learn[1,,] <- npars[1,adapt$stimulus$init_par]
    }

    # update
    for(trial in 1:nrow(learn)) {
      Ri <- da$trials==trial

      # Update parameters
      Q <- learn[cbind(index,1)][Ri]
      if(('ws' %in% adapt$stimulus$output_par_names) & ('wd' %in% adapt$stimulus$output_par_names)) {
        indexOther <- attr(da, 'adapt')[[s]]$indexOther
        Q <- cbind(Q, learn[cbind(indexOther,1)][Ri])
      }

      npars[Ri,adapt$stimulus$output_name] <- adapt$stimulus$output_fun(npars[Ri,adapt$stimulus$output_par_names], Q)

      ## check if response needs to be generated
      if(add_response) {
        # simulate response
        Rrt <- rfun(factor(levels=da[Ri,'lS']),npars[Ri,])
        da[Ri, 'rt'] <- Rrt[,'rt']
        da[Ri, 'lRS'] <- as.character(Rrt[,'R'])
        da[Ri, 'R'] <- ifelse(all(as.character(da[Ri,'s_left']) == as.character(Rrt[,'R'])), 'left', 'right')
        da[Ri, 'winner'] <- da[Ri, 'lRS']==da[Ri,'lS']

        ## AGAIN, here we assume *ONLY* the winner receives feedback!
        reward <- setNames(rbinom(length(da[Ri,'lS']), 1, pArray[trial,as.character(da[Ri,'lS'])]), da[Ri,'lS'])
        outcomes[trial,as.character(da[Ri&da$winner,'lRS'])] <- reward[da[Ri,'winner']]
      }

      ## update
      if(trial < nrow(learn)) {
        for(column in 1:ncol(learn)) {
          if(!is.na(outcomes[trial,column])) {
            # feedback was given, so update
            learn[trial+1,column,] <- adapt$stimulus$adapt_fun(lastValues=learn[trial,column,],
                                                               parameters=npars[trial,adapt$stimulus$adapt_par],
                                                               reward=outcomes[trial,column])
          } else {
            learn[trial+1,column,] <- learn[trial,column,]
          }
        }
      }
    }
    if (return_all) return(list(learn=learn,pars=npars,data=da))
    if (return_learning) return(learn)
    if (mapped_p) return(list(data=da,pars=npars))
    if (add_response) return(da)
    return(npars)
  } else {
    index <- attr(da,"adapt")[[s]]$stimulus$index
    learn <- attr(da,"adapt")[[s]]$stimulus$learn
    npars <- npars[da$subjects==s,]
    da <- da[da$subjects==s,]
    add_response <- any(is.na(da$R))
    nAcc <- dim(adapt$stimulus$targets)[1]
    namAcc <- dimnames(adapt$stimulus$targets)[[1]]
    nStim <- dim(adapt$stimulus$targets)[2]
    namStim <- colnames(adapt$stimulus$targets)
    # Maximum number of Q updates
    maxQupdates <- dim(index)[1]/nAcc
    # fill an array of parameters, dim = accumulator x trial x
    # stimulus x parameter type ("v0"    "B"     "t0"    "alpha" "w"     "q0"    "A" )
    parArr <- aperm(
      array(apply(index,2,function(x){npars[x,]}),
            dim=c(nAcc,maxQupdates,dim(npars)[2],nStim),
            dimnames=list(namAcc,NULL,dimnames(npars)[[2]],namStim)),
      c(1,2,4,3)) # reorder to make look up quick in loop
    # fill prob reward array: trials x stimulus x stimulus component
    # pReward <- array(as.vector(sapply(dimnames(adapt$stimulus$targets)[[1]],function(x){
    #   da[index,paste("p",x,sep="_")]})),dim=c(nAcc,maxQupdates,nStim,nAcc),
    #   dimnames=list(namAcc,NULL,namStim,namAcc))[1,,,]
    pReward <- array(array(as.vector(sapply(dimnames(adapt$stimulus$targets)[[1]],function(x){
      da[index,paste("p",x,sep="_")]})),dim=c(nAcc,maxQupdates,nStim,nAcc))[1,,,],
      dim=c(maxQupdates,nStim,nAcc),dimnames=list(NULL,namStim,namAcc))

    # Extract Q values and update
    if (add_response) {
      da_reward <- da_rt <- rep(NA,dim(da)[1])         # Rewards and rts
      da_R <- factor(da_reward,levels=namAcc)          # Response factor
      Ri <- array(index,dim=c(nAcc,maxQupdates,nStim)) # Row indices
    } else {
      # maxQupdates rows, nStim columns
      Rmat <- array(da[index,"R"], dim=c(nAcc,maxQupdates,nStim))[1,,]
      reward_mat <- array(da[index,"reward"], dim=c(nAcc,maxQupdates,nStim))[1,,]
    }
    for (i in 1:maxQupdates)  {
      ok <- !is.na(parArr[1,i,,1]) # stimuli that need updating
      if (i==1) learn[,ok,i] <- parArr[,1,,adapt$stimulus$init_par] # Initialize
      # pick out fixed pars
      pars <- setNames(data.frame(matrix(parArr[,i,,][,ok,adapt$fixed_pars,drop=FALSE],
                                         ncol=length(adapt$fixed_pars))),adapt$fixed_pars)

      # SM: advantage framework, ONLY 2AFC SO FAR!
      Q <- as.vector(learn[,,i][,ok])
      if(('ws' %in% adapt$stimulus$output_par_names) & ('wd' %in% adapt$stimulus$output_par_names)) {
        Q <- cbind(Q, as.vector(learn[c(2,1),,i][,ok]))
      }

      # calculate and add output_par
      pars[[adapt$stimulus$output_name]] <- adapt$stimulus$output_fun(
        output_pars=matrix(parArr[,i,ok,adapt$stimulus$output_par_names],
                           ncol=length(adapt$stimulus$output_par_names)),
        Q=Q) #as.vector(learn[,,i][,ok]))


      if (add_response) { # Simulate trial
        Rrt <- rfun(factor(levels=dimnames(adapt$stimulus$targets)[[1]]),pars)
        Rfac <- Rrt[,"R"]
        reward <- rbinom(length(Rfac),1,pReward[i,ok,][cbind(1:length(Rfac),as.numeric(Rfac))])
        # harvest new trial info
        da_rt[Ri[,i,][,ok,drop=FALSE]] <- rep(Rrt[,"rt"],each=nAcc)
        da_reward[Ri[,i,][,ok,drop=FALSE]] <- rep(reward,each=nAcc)
        da_R[Ri[,i,][,ok,drop=FALSE]] <- rep(Rfac,each=nAcc)
      } else { # Extract trial information
        Rfac <- factor(Rmat[i,ok],levels=namAcc)
        reward <- reward_mat[i,ok]
      }
      if (i<maxQupdates) { # Update Q value
        learn[,ok,i+1] <- learn[,ok,i] # Copy last for all components
        imat <- cbind(as.numeric(Rfac),1:sum(ok)) # Components to update
        learn[,ok,i+1][imat] <- adapt$stimulus$adapt_fun(Qlast=learn[,,i][,ok,drop=FALSE][imat],
                                                         adapt_par=parArr[,i,,adapt$stimulus$adapt_par][,ok,drop=FALSE][imat],reward)
      }
      if (mapped_p | !add_response | return_all) # Output will be new so update parrArr
        parArr[,i,,adapt$stimulus$output_name][!is.na(parArr[,i,,adapt$stimulus$output_name])] <-
        pars[,adapt$stimulus$output_name]
    }
    if (add_response) {
      da$R <- da_R
      da$reward <- da_reward
      da$rt <- da_rt
      attr(da,"adapt")[[s]]$stimulus$learn <- learn
    }
    if (mapped_p | !add_response | return_all) {
      stim_output <- parArr[,,,adapt$stimulus$output_name][!is.na(index)]
      # Protect against numerical problems, may screw up dfun for some models
      stim_output[is.na(stim_output)|is.nan(stim_output)] <- -Inf
      npars[index[!is.na(index)],adapt$stimulus$output_name] <- stim_output
    }

    if (return_all) return(list(learn=learn,pars=npars,data=da))
    if (return_learning) return(learn)
    if (mapped_p) return(list(data=da,pars=npars))
    if (add_response) return(da)
    npars
  }
}




# data=dadm
adapt_data <- function(data,design,model,pars,
                       add_response=FALSE,return_learning=FALSE,mapped_p=FALSE,return_all=FALSE)
  # runs learning, and either returns learning, or data with responses and rt
  # added (either if data has NA in R or if add_Response), or adapted parameters
  # or data + parameters (mapped_p)
{
  if (return_all & mapped_p) {
    warning("return_all and mapped_p incompatible, going with the latter")
    return_all <- FALSE
  }
  if (add_response) data$R <- NA # Force new data generation ignoring old data
  if ( all(is.na(data$R)) ) add_response <- TRUE # make new if none supplied
  # add augmentation
  if (is.null(attr(data,"adapt"))) {
    attr(data,"adapt") <- setNames(
      lapply(levels(data$subjects),augment,da=data,design=design),
      levels(data$subjects))
    attr(data,"adapt")$design <- design$adapt
  }
  daList <- setNames(
    lapply(levels(data$subjects),update_pars,npars=pars,da=data,return_all=return_all,
           rfun=model()$rfun,return_learning=return_learning,mapped_p=mapped_p),
    levels(data$subjects))
  if (return_learning) return(daList)
  adapt <- attr(data,"adapt")
  if (mapped_p | return_all) {
    data <- do.call(rbind,lapply(daList,function(x)x$data))
    pars <- do.call(rbind,lapply(daList,function(x)x$pars))
    if (mapped_p) return(cbind(data[,!(names(data) %in% c("R","rt"))],pars))
    for (i in names(daList)) adapt[[i]] <- attr(daList[[i]]$data,"adapt")[[i]]
    data <- cbind(data,pars)
    attr(data,"adapt") <- adapt
    return(list(learn=lapply(daList,function(x)x$learn),data=data))
  }
  for (i in names(daList)) adapt[[i]] <- attr(daList[[i]],"adapt")[[i]]
  data <- do.call(rbind,daList)
  attr(data,"adapt") <- adapt
  data
}


# augment = function(s,da,design)
#   # Adds attributes to augmented data
#   # learn: empty array for Q values with dim = choice alternative (low,high) x
#   #   stimulus x trials (max across stimuli)
#   # index: look up (row number) for stimuli in da, a matrix dim  = max trials x
#   #   stimulus matrix (rows for each choice alternative contiguous)
# {
#   getIndex <- function(typei,cname,da,maxn) {
#     out <- which(da[,cname]==typei)
#     c(out,rep(NA,maxn-length(out)))
#   }
#
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
# # ARCHITECTURE OF FOLLOWING FUNCTION WILL NEED UPDATING FOR MULTIPLE ADAPT TYPES
# # s="10";npars=pars;da=data; rfun=model()$rfun;return_learning=FALSE;mapped_p=FALSE
# update_pars = function(s,npars,da,rfun=NULL,return_learning=FALSE,mapped_p=FALSE,
#                        return_all=FALSE)
#   # for subject s
#   # Either return da filling in responses and RT if da has no responses (in
#   # in which case rfun must be supplied), or return npars filling in adapted
#   # parameters or if return_learning returns learning (e.g., Q values)
# {
#
#   adapt <- attr(da,"adapt")$design
#   index <- attr(da,"adapt")[[s]]$stimulus$index
#   learn <- attr(da,"adapt")[[s]]$stimulus$learn
#   npars <- npars[da$subjects==s,]
#   da <- da[da$subjects==s,]
#   add_response <- any(is.na(da$R))
#   nAcc <- dim(adapt$stimulus$targets)[1]
#   namAcc <- dimnames(adapt$stimulus$targets)[[1]]
#   nStim <- dim(adapt$stimulus$targets)[2]
#   namStim <- colnames(adapt$stimulus$targets)
#   # Maximum number of Q updates
#   maxQupdates <- dim(index)[1]/nAcc
#   # fill an array of parameters, dim = accumulator x trial x
#   # stimulus x parameter type ("v0"    "B"     "t0"    "alpha" "w"     "q0"    "A" )
#   parArr <- aperm(
#     array(apply(index,2,function(x){npars[x,]}),
#           dim=c(nAcc,maxQupdates,dim(npars)[2],nStim),
#           dimnames=list(namAcc,NULL,dimnames(npars)[[2]],namStim)), c(1,2,4,3)) # reorder to make look up quick in loop
#   # fill prob reward array: trials x stimulus x stimulus component
#   # pReward <- array(as.vector(sapply(dimnames(adapt$stimulus$targets)[[1]],function(x){
#   #   da[index,paste("p",x,sep="_")]})),dim=c(nAcc,maxQupdates,nStim,nAcc),
#   #   dimnames=list(namAcc,NULL,namStim,namAcc))[1,,,]
#   pReward <- array(array(as.vector(sapply(dimnames(adapt$stimulus$targets)[[1]],function(x){
#     da[index,paste("p",x,sep="_")]})),dim=c(nAcc,maxQupdates,nStim,nAcc))[1,,,],
#     dim=c(maxQupdates,nStim,nAcc),dimnames=list(NULL,namStim,namAcc))
#
#   # Extract Q values and update
#   if (add_response) {
#     da_reward <- da_rt <- rep(NA,dim(da)[1])         # Rewards and rts
#     da_R <- factor(da_reward,levels=namAcc)          # Response factor
#     Ri <- array(index,dim=c(nAcc,maxQupdates,nStim)) # Row indices
#   } else {
#     # maxQupdates rows, nStim columns
#     Rmat <- array(da[index,"R"], dim=c(nAcc,maxQupdates,nStim))[1,,]
#     reward_mat <- array(da[index,"reward"], dim=c(nAcc,maxQupdates,nStim))[1,,]
#   }
#   for (i in 1:maxQupdates)  {
#     ok <- !is.na(parArr[1,i,,1]) # stimuli that need updating
#     if (i==1) learn[,ok,i] <- parArr[,1,,adapt$stimulus$init_par] # Initialize
#     # pick out fixed pars
#     pars <- setNames(data.frame(matrix(parArr[,i,,][,ok,adapt$fixed_pars,drop=FALSE],
#                                        ncol=length(adapt$fixed_pars))),adapt$fixed_pars)
#     # calculate and add output_par
#     pars[[adapt$stimulus$output_name]] <- adapt$stimulus$output_fun(
#       output_pars=matrix(parArr[,i,ok,adapt$stimulus$output_par_names],
#                          ncol=length(adapt$stimulus$output_par_names)),
#       Q=as.vector(learn[,,i][,ok]))
#     if (add_response) { # Simulate trial
#       Rrt <- rfun(factor(levels=dimnames(adapt$stimulus$targets)[[1]]),pars)
#       Rfac <- Rrt[,"R"]
#       reward <- rbinom(length(Rfac),1,pReward[i,ok,][cbind(1:length(Rfac),as.numeric(Rfac))])
#       # harvest new trial info
#       da_rt[Ri[,i,][,ok,drop=FALSE]] <- rep(Rrt[,"rt"],each=nAcc)
#       da_reward[Ri[,i,][,ok,drop=FALSE]] <- rep(reward,each=nAcc)
#       da_R[Ri[,i,][,ok,drop=FALSE]] <- rep(Rfac,each=nAcc)
#     } else { # Extract trial information
#       Rfac <- factor(Rmat[i,ok],levels=namAcc)
#       reward <- reward_mat[i,ok]
#     }
#     if (i<maxQupdates) { # Update Q value
#       learn[,ok,i+1] <- learn[,ok,i] # Copy last for all components
#       imat <- cbind(as.numeric(Rfac),1:sum(ok)) # Components to update
#       learn[,ok,i+1][imat] <- adapt$stimulus$adapt_fun(Qlast=learn[,,i][,ok,drop=FALSE][imat],
#                                                        adapt_par=parArr[,i,,adapt$stimulus$adapt_par][,ok,drop=FALSE][imat],reward)
#     }
#     if (mapped_p | !add_response | return_all){
#       # Output will be new so update parrArr
#       parArr[,i,,adapt$stimulus$output_name][!is.na(parArr[,i,,adapt$stimulus$output_name])] <-pars[,adapt$stimulus$output_name]
#     }
#   }
#   if (add_response) {
#     da$R <- da_R
#     da$reward <- da_reward
#     da$rt <- da_rt
#     attr(da,"adapt")[[s]]$stimulus$learn <- learn
#   }
#   if (mapped_p | !add_response | return_all) {
#     stim_output <- parArr[,,,adapt$stimulus$output_name][!is.na(index)]
#     # Protect against numerical problems, may screw up dfun for some models
#     stim_output[is.na(stim_output)|is.nan(stim_output)] <- -Inf
#     npars[index[!is.na(index)],adapt$stimulus$output_name] <- stim_output
#   }
#   if (return_all) return(list(learn=learn,pars=npars,data=da))
#   if (return_learning) return(learn)
#   if (mapped_p) return(list(data=da,pars=npars))
#   if (add_response) return(da)
#   npars
# }
#
# # data=dadm
# adapt_data <- function(data,design,model,pars,
#                        add_response=FALSE,return_learning=FALSE,mapped_p=FALSE,return_all=FALSE)
#   # runs learning, and either returns learning, or data with responses and rt
#   # added (either if data has NA in R or if add_Response), or adapted parameters
#   # or data + parameters (mapped_p)
# {
#   if (return_all & mapped_p) {
#     warning("return_all and mapped_p incompatible, going with the latter")
#     return_all <- FALSE
#   }
#   if (add_response) data$R <- NA # Force new data generation ignoring old data
#   if ( all(is.na(data$R)) ) add_response <- TRUE # make new if none supplied
#   # add augmentation
#   if (is.null(attr(data,"adapt"))) {
#     attr(data,"adapt") <- setNames(
#       lapply(levels(data$subjects),augment,da=data,design=design),
#       levels(data$subjects))
#     attr(data,"adapt")$design <- design$adapt
#   }
#   daList <- setNames(
#     lapply(levels(data$subjects),update_pars,npars=pars,da=data,return_all=return_all,
#            rfun=model()$rfun,return_learning=return_learning,mapped_p=mapped_p),
#     levels(data$subjects))
#   if (return_learning) {
#     return(daList)
#   }
#   adapt <- attr(data,"adapt")
#   if (mapped_p | return_all) {
#     data <- do.call(rbind,lapply(daList,function(x)x$data))
#     pars <- do.call(rbind,lapply(daList,function(x)x$pars))
#     if (mapped_p) {
#       return(cbind(data[,!(names(data) %in% c("R","rt"))],pars))
#     }
#     for (i in names(daList)){
#       adapt[[i]] <- attr(daList[[i]]$data,"adapt")[[i]]
#     }
#     data <- cbind(data,pars)
#     attr(data,"adapt") <- adapt
#     return(list(learn=lapply(daList,function(x)x$learn),data=data))
#   }
#   for (i in names(daList)){
#     adapt[[i]] <- attr(daList[[i]],"adapt")[[i]]
#   }
#   data <- do.call(rbind,daList)
#   attr(data,"adapt") <- adapt
#   data
# }
