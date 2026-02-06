minimal_design <- function(design, covariates = NULL, drop_subjects = TRUE,
                           n_trials = 1, add_acc = TRUE, drop_R = TRUE,
                           drop_R_levels = TRUE, do_functions = TRUE,
                           group_design = NULL, ...) {
  dots <- add_defaults(list(...), verbose = TRUE)
  if(!is.null(design$Ffactors)) design <- list(design)
  out <- list()
  if(!is.null(group_design)){
    group_factors <- unique(unlist(lapply(attr(group_design, "Flist"), function(x) all.vars(x[[3]]))))
  }
  for(i in 1:length(design)){
    cur_des <- design[[i]]
    ## 1.  Check and construct the factorial backbone
    if (is.null(cur_des$Ffactors) || !is.list(cur_des$Ffactors)) {
      stop("`cur_des` must contain a list element called `Ffactors`.")
    }
    # Add Rlevels
    cur_des$Ffactors$R <- cur_des$Rlevels
    # Remove subjects
    if(drop_subjects){
      cur_des$Ffactors <- cur_des$Ffactors[names(cur_des$Ffactors) != "subjects"]
    } else{
      if(is.null(cur_des$Ffactors$subjects))cur_des$Ffactors$subjects <- factor(1)
    }
    # if(!is.null(group_design)){
    #   cur_des$Ffactors <- cur_des$Ffactors[!names(cur_des$Ffactors) %in% group_factors]
    #   cur_des$Fcovariates <- cur_des$Fcovariates[!names(cur_des$Fcovariates) %in% group_factors]
    # }


    if(n_trials > 1){
      cur_des$Ffactors <- c("trials" = list(1:n_trials), cur_des$Ffactors)
    }
    fac_df <- do.call(expand.grid, c(cur_des$Ffactors))
    n <- nrow(fac_df)
    fac_df[] <- lapply(fac_df, as.factor)
    if(n_trials > 1 & !drop_subjects){
      fac_df$trials <- as.numeric(ave(as.character(fac_df$subjects), as.character(fac_df$subjects), FUN = seq_along))
    } else if(n_trials > 1){
      fac_df$trials <- 1:nrow(fac_df)
    }

    ## 2.  Add covariates (if requested)
    if (!is.null(cur_des$Fcovariates)) {
      for (cv in cur_des$Fcovariates) {
        if(!is.null(list(...)$emc)){
          dat <- list(...)$emc
          if(is.data.frame(dat)){
            covariates[[cv]] <- dat[,cv]
          } else if(!(is(dat, "emc") & is.null(dat[[1]]$subjects))){
            dat <- get_data(dat)
            covariates[[cv]] <- dat[,cv]
          }
        }
        if(!cv %in% names(covariates)){
          if(dots$verbose) message(paste0("Imputing ", cv, " with random values"))
          vec <- rnorm(n)
        } else{
          vec <- covariates[[cv]]
          if(length(vec) < nrow(fac_df)){
            vec <- rep(vec, ceiling(nrow(fac_df)/length(vec)))
          }
        }
        fac_df[[cv]] <- vec[seq_len(n)]
      }
    }

    ## 3. Add accumulators
    if(add_acc){
      fac_df <- add_accumulators(fac_df, matchfun = cur_des$matchfun, type = cur_des$model()$type)
    }
    if(!is.null(fac_df$R) & drop_R_levels){
      fac_df <- fac_df[fac_df$R == unique(fac_df$R)[1],]
    }
    if(!is.null(fac_df$R) & drop_R){
      fac_df <- fac_df[,!colnames(fac_df) %in% c("R", "winner")]
    }

    ## 4.  Derive additional columns via Ffunctions (if any)
    if (!is.null(cur_des$Ffunctions) & do_functions) {
      funs <- cur_des$Ffunctions
      if (!is.list(funs)) funs <- list(funs)

      for (j in seq_along(funs)) {
        res <- funs[[j]](fac_df)
        nm <- names(funs)[j]
        fac_df[[nm]] <- res
      }
    }
    if(!is.null(group_design)){
      fac_df <- minimal_group_design(fac_df, group_factors, attr(group_design, "data"))
    }

    out[[i]] <- fac_df
  }
  if(length(out) == 1) out <- out[[1]]
  return(out)
}

minimal_group_design <- function(data, group_factors, group_data){
  if(nrow(group_data) != length(unique(data$subjects))){
    if(length(unique(data$subjects)) != 1) stop("Group design doesn't match subject inputs")
    input_data <- data[!duplicated(data[,colnames(data) != "trials"]),]
    input_group <- unique(group_data)
    out <- vector("list", nrow(input_group))
    for(i in 1:length(out)){
      tmp <- input_data
      for(j in 1:length(group_factors)){
        tmp[,group_factors[j]] <- input_group[i,j]
      }
      tmp$subjects <- i
      out[[i]] <- tmp
    }

  } else{
    unq_subjects <- unique(data$subjects)
    out <- setNames(vector("list", length(unq_subjects)), unq_subjects)
    k <- 0
    for(sub in unq_subjects){
      k <- k + 1
      idx <- data$subjects == sub
      group_tmp <- group_data[rep(k, each = sum(idx)), , drop = F]
      out[[sub]] <- cbind(data[idx,], group_tmp)
    }
  }

  data <- do.call(rbind, out)
  rownames(data) <- NULL
  data <- add_trials(data[order(data$subjects),])
  return(data)
}


do_map <- function(draws, map, by_subject, design,
                   n_trials = 1, data = NULL, functions = NULL,
                   add_recalculated = FALSE, selection = "alpha",
                   group_design = NULL,...){
  par_idx <- 0
  all_out <- list()
  if(!is.null(data$subjects)) data <- list(data)
  for(i in 1:length(design)){
    cur_des <- design[[i]]
    # See which design we have to get out
    cur_idx <- par_idx + seq_len(length(sampled_pars(cur_des)))
    par_idx <- max(cur_idx)
    cur_draws <- draws[cur_idx,,,drop = F]
    if(grepl("MRI", cur_des$model()$type)){
      # MRI mapping is rudimentary for now
      all_out[[i]] <- map_MRI(cur_draws)
      next
    }
    joint <- FALSE
    if(any(grepl("|", rownames(cur_draws), fixed =  T))){
      # Get out joint
      joint <- TRUE
      prefix <- unique(gsub("[|].*", "", rownames(cur_draws)))
      rownames(cur_draws) <- sub("^[^|]*\\|", "", rownames(cur_draws))
    }

    out <- mapper_wrapper(map = map, by_subject, cur_draws, cur_des, n_trials = n_trials,
                          data = data[[i]], functions = functions,
                          add_recalculated = add_recalculated,
                          group_design = group_design, ...)

    if(joint){
      rownames(out) <- paste0(prefix, "|", rownames(out))
    }
    all_out[[i]] <- out
  }
  out <- abind(all_out, along = 1)
  return(out)
}

# To fix, apply add_recalculated
mapper_wrapper <- function(map, by_subject = FALSE, par_mcmc, design, n_trials = NULL, data = NULL,
                           functions = NULL, add_recalculated = FALSE,
                           group_design = NULL, ...){
  res <- par_data_map(par_mcmc, design, n_trials = n_trials,
                      data = data, functions = functions,
                      add_recalculated = add_recalculated,
                      group_design = group_design, ...)

  pars <- res$pars
  data <- res$data
  if(!by_subject){
    data$subjects <- 1
    factor(data$subjects)
  }
  n_pars <- dim(pars)[3]

  fmls <- setNames(
    lapply(design$Flist, function(f) as.formula(paste("~", deparse(f[[3]])))),
    vapply(design$Flist, function(f) deparse(f[[2]]), character(1))
  )

  fml_names <- names(fmls)
  par_names <- dimnames(pars)[[3]]

  # Some map wrangling
  if(!is.logical(map)){ # Not just a simple boolean
    # Map specified as a list
    if(is.list(map)){
      if(is.null(names(map))){
        # Unnamed map
        if(length(map) != 1 & length(map) != n_pars) stop("map list length not equal to number of parameters")
        if(length(map) == 1){
          # Map of length 1, assume that they should all be replicated
          map <- replicate(n_pars, map)
        }
      } else{ # Named map
        if(any(!par_names %in% names(map))){
          if(is.character(map[[1]])){ # Factor input
            map[setdiff(par_names, names(map))] <- ""
          } else{ # Otherwise assume it's a formula
            map[setdiff(par_names, names(map))] <- as.formula("~ 1")
          }
          map <- map[par_names]
        }

      }
    } else{
      map <- replicate(n_pars, list(map))
    }
  } else{
    map <- rep(map, n_pars)
  }

  # Before this we should compress covariates
  colnams <- colnames(data)
  for(i in 1:ncol(data)){
    if(colnams[i] %in% c("subjects", "rt", "R", "trials")) next
    cur_col <- data[,i]
    if(!is.factor(cur_col) && length(unique(cur_col)) > 6){
      new_col <- rep("high", length(cur_col))
      new_col[cur_col < median(cur_col)] <- "low"
      data[,i] <- new_col
    }
  }

  df <- data
  df[colnames(df)] <- lapply(colnames(df), function(v) {
    f <- df[[v]]
    factor(paste0(v, f))
  })
  subjects <- unique(data$subjects)
  all_pars <- setNames(vector("list", length(subjects)), subjects)
  for(sub in subjects){
    out <- list()
    idx <- data$subjects == sub
    for(i in 1:n_pars){
      # First case, map = TRUE
      if(isTRUE(map[i]) || is.character(map[[i]])){
        if(isTRUE(map[i])){
          vars <- all.vars(fmls[[par_names[i]]])
        } else{
          vars <- map[[i]]
          if(length(vars) == 1 && vars == "") vars <- character(0)
        }

        if(length(vars) > 0){
          cells <- interaction(df[idx,vars], sep = "_", drop = TRUE, lex.order = TRUE)
          g <- nlevels(cells)
          ng <- as.vector(table(cells))
          res <- rowsum(pars[idx,,i], cells) / ng            # (g x k)

        } else{
          if (sum(idx)==1) res <- pars[idx,,i] else
            res <- colMeans(pars[idx,,i])
        }
      } else{ # Third case map is a formula
        S  <- model.matrix(map[[i]], data = data[idx,])
        ng <- colSums(S)
        # g x k means:
        res <- (t(S) %*% pars[idx,,i]) / ng
      }
      if(is.vector(res)){
        res <- t(res)
        rownames(res) <- par_names[i]
      } else{
        rownames(res) <- paste0(par_names[i], "_", rownames(res))
      }
      out[[i]] <- round(res, 8)
    }
    all_pars[[sub]] <- do.call(rbind, out)
  }
  # This is to ensure no-one is missing a cell
  unq_rows <- unique(unlist(lapply(all_pars, rownames)))
  res <- lapply(all_pars, function(x){
    missing <- setdiff(unq_rows, rownames(x))

    # Create NA matrix for missing rows
    add <- matrix(NA, nrow = length(missing), ncol = ncol(x),
                  dimnames = list(missing, colnames(x)))

    # Bind and reorder
    mat_full <- rbind(x, add)[unq_rows,]
    return(mat_full)
  })

  arr <- abind(res, along = 3)    # rows x cols x list
  arr <- aperm(arr, c(1, 3, 2))   # rows x list x cols
  return(arr)
}


par_data_map <- function(par_mcmc, design, n_trials = NULL, data = NULL,
                         functions = NULL, add_recalculated = FALSE,
                         group_design = NULL,...){
  design <- design
  model <- design$model
  if ( is.null(data) ) {
    design$Ffactors$subjects <- colnames(par_mcmc)
    if ( is.null(n_trials) )
      stop("If data is not provided need to specify number of trials")
    design$Fcovariates <- design$Fcovariates[!design$Fcovariates %in% names(functions)]
    data <- minimal_design(design, covariates = list(...)$covariates,
                           drop_subjects = F, n_trials = n_trials, add_acc=F,
                           drop_R = F, group_design = group_design)
  }
  if(!is.null(functions)){
    for(i in 1:length(functions)){
      data[[names(functions)[i]]] <- functions[[i]](data)
    }
  }
  if (!is.factor(data$subjects)) data$subjects <- factor(data$subjects)
  if (!is.null(model)) {
    if (!is.function(model)) stop("model argument must  be a function")
    if ( is.null(model()$p_types) ) stop("model()$p_types must be specified")
    if ( is.null(model()$Ttransform) ) stop("model()$Ttransform must be specified")
  }
  data <- do.call(rbind, lapply(split(data, data$subjects), function(x){
    x <- x[!duplicated(x[,colnames(x) != "rt"]),]
    return(x)
  }))
  data <- add_trials(data[order(data$subjects),])


  model <- design$model
  data <- design_model(
    add_accumulators(data,design$matchfun,simulate=TRUE,type=model()$type,Fcovariates=design$Fcovariates),
    design,model,add_acc=FALSE,compress=FALSE,verbose=FALSE,
    rt_check=FALSE)

  n_mcmc <- dim(par_mcmc)[3]
  n_pars <- length(model()$p_types)
  n_subs <- ncol(par_mcmc)
  for(i in 1:n_mcmc){
    parameters <- t(as.matrix(par_mcmc[,,i], nrow = n_pars, ncol = n_subs))
    pars <- do_transform(parameters, model()$pre_transform) #t(apply(parameters, 1, do_pre_transform, model()$pre_transform))
    if(nrow(parameters) == length(unique(data$subjects))){
      design$Ffactors$subjects <- unique(data$subjects)
    }

    rownames(pars) <- design$Ffactors$subjects
    pars <- map_p(add_constants(pars,design$constants),data, model(), return_trend_pars = TRUE)
    exclude_transform <- rep(FALSE, ncol(pars))
    if (!is.null(model()$trend)) {
      phases <- vapply(model()$trend, function(x) x$phase, character(1))
      if (any(phases == "pretransform")) pars <- prep_trend_phase(data, model()$trend, pars, "pretransform",
                                                                  return_trend_pars = TRUE)
      exclude_transform <- unlist(lapply(model()$trend, function(x){
        if(x$phase != "posttransform"){
          return(x$trend_pnames)
        } else{
          return(NULL)
        }
      }))
      exclude_transform <- colnames(pars) %in% exclude_transform
    }
    pars[,!exclude_transform] <- do_transform(pars[,!exclude_transform], model()$transform)
    if (!is.null(model()$trend)) {
      if (any(phases == "posttransform")) pars <- prep_trend_phase(data, model()$trend, pars, "posttransform",
                                                                   return_trend_pars = TRUE)
    }
    if(add_recalculated) pars <- model()$Ttransform(pars, data)
    if(i == 1){
      # Ttransform could add unwanted friends, so safest to just
      # figure out the dimensions here
      out <- array(NA, dim = c(nrow(pars), n_mcmc, ncol(pars)),
                   dimnames = list(NULL, NULL, colnames(pars)))
    }
    out[,i,] <- pars
  }
  return(list(data = data, pars = out))
}

group_mapping <- function(samples, selection) {
  if(ncol(samples) < 2) stop("Please use type = 'single' for designs with 1 subject")
  n_params <- dim(samples)[1]
  n_iters  <- dim(samples)[3]
  if(selection == "mu"){
    return(list(samples = list(theta_mu_mean = apply(samples, c(1,3), mean, na.rm = TRUE))))
  } # Else we need theta_var
  covs  <- array(NA_real_, dim = c(n_params, n_params, n_iters),
                 dimnames = list(dimnames(samples)[[1]],
                                 dimnames(samples)[[1]], NULL))
  for (iter in seq_len(n_iters)) {
    tmp <- t(samples[,,iter])
    covs[,,iter] <- cov(tmp, use = "pairwise.complete.obs")
  }
  return(list(samples = list(theta_var = covs)))
}

map_selecter <- function(map, selection){
  if(!isFALSE(map)){
    if(selection %in% c("mu", "sigma2", "Sigma", "correlation", "covariance")){
      return("alpha")
    }
  }
  return(selection)
}

