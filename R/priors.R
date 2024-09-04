#' Prior specification
#'
#' Specify priors for the chosen model. These values are entered manually by default but can be
#' recycled from another prior (given in the `update` argument).
#'
#' Where a value is not supplied, the user is prompted to enter
#' numeric values (or functions that evaluate to numbers).
#'
#' To get the default prior for a type, run: `get_prior_{type}(design = design, sample = F)`
#'
#' E.g.: `get_prior_diagonal(design = design, sample = F)`
#'
#' @param design Design list for which a prior is constructed, typically the output of `design()`
#' @param update Prior list from which to copy values
#' @param type Character. What type of group-level model you plan on using i.e. `diagonal`
#' @param ask Character. For which parameter types to ask for prior specification, i.e. `Sigma`, `mu` or `loadings` for factor models
#' @param fill_default Boolean, If `TRUE` will fill all non-specified parameters, and parameters outside of `ask`, to default values
#' @param ... Either values to prefill, i.e. `theta_mu_mean = c(1:6)`, or additional arguments such as `n_factors = 2`
#' @return A prior list object
#' @examples
#' # First define a design for the model
#' design_DDMaE <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#' # Then set up a prior using prior
#' p_vector=c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
#'                      t0=log(.2),Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))
#' psd <- c(v_Sleft=1,v_Sright=1,a=.3,a_Eneutral=.3,a_Eaccuracy=.3,
#'                      t0=.4,Z=1,sv=.4,SZ=1)
#' # Here we left the variance prior at default
#' prior_DDMaE <- prior(design_DDMaE,mu_mean=p_vector,mu_sd=psd)
#' # Also add a group-level variance prior:
#' pscale <- c(v_Sleft=.6,v_Sright=.6,a=.3,a_Eneutral=.3,a_Eaccuracy=.3,
#'                              t0=.2,Z=.5,sv=.4,SZ=.3)
#' df <- .4
#' prior_DDMaE <- prior(design_DDMaE,mu_mean=p_vector,mu_sd=psd, A = pscale, df = df)
#' # If we specify a new design
#' design_DDMat0E <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~E, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#' # We can easily update the prior
#' prior_DDMat0E <- prior(design_DDMat0E, update = prior_DDMaE)
#' @export
prior <- function(design, type = "standard", update = NULL,
                      ask = NULL, fill_default = TRUE, ...){
  if(!is.null(update)){
    type <- attr(update, "type")
  }
  prior <- get_objects(design = design, type = type, ...)
  args <- list(...)
  if(!is.null(args$mu_mean)){
    args$theta_mu_mean <- args$mu_mean
    if(is.null(names(args$theta_mu_mean))) names(args$theta_mu_mean) <- names(prior$prior$theta_mu_mean)
  }
  if(!is.null(args$mu_sd)){
    args$theta_mu_var <- args$mu_sd^2
    if(is.null(names(args$theta_mu_var))) names(args$theta_mu_var) <- names(prior$prior$theta_mu_mean)
  }
  if(!is.null(args$pmean)){
    args$theta_mu_mean <- args$pmean
    if(is.null(names(args$theta_mu_mean))) names(args$theta_mu_mean) <- names(prior$prior$theta_mu_mean)
  }
  if(!is.null(args$psd)){
    args$theta_mu_var <- args$psd^2
    if(is.null(names(args$theta_mu_var))) names(args$theta_mu_var) <- names(prior$prior$theta_mu_mean)
  }

  if(!is.null(update)){
    for(name in names(update)){
      if(is.null(args[[name]])){
        if(!is.null(dim(update[[name]]))){
          update[[name]] <- diag(update[[name]])
        } else{
          update[[name]] <- update[[name]]
        }
        args[[name]] <- update[[name]]
      } else{
        # First make sure we only take diagonals
        if(!is.null(dim(update[[name]]))){
          prior_par_names <- names(diag(prior$prior[[name]]))
          upd <- diag(update[[name]])
        } else{
          prior_par_names <- names(prior$prior[[name]])
          upd <- update[[name]]
        }
        # now fill in
        upd_names <- names(upd)
        if(length(upd) == 1){
          if(!name %in% names(args)) args[[name]] <- upd
        } else{
          for(i in 1:length(upd)){
            if(upd_names[i] %in% prior_par_names & (!upd_names[i] %in% names(args[[name]]))){
              args[[name]][upd_names[i]] <- upd[i]
            }
          }
        }
      }
    }
  }
  for(group in names(prior$groups)){
    if(is.null(ask) & !fill_default){
      group_to_do <- TRUE
    } else{
      group_to_do <- group %in% ask
    }
    if(group_to_do) cat(paste0(prior$group_descriptions[[group]], "\n \n"))
    for(pri in prior$groups[[group]]){
      if(pri %in% names(prior$descriptions)){ # This excluded prior$theta_mu_invar
        if(pri %in% names(args)){ # Already specified in ellipsis, so fill in
          input <- args[[pri]]
          if(!is.null(dim(prior$prior[[pri]]))){
            to_check <- diag(prior$prior[[pri]])
            if(length(input) == length(to_check) & is.null(names(input))){
              prior$prior[[pri]] <- input
              to_do <- rep(F, length(to_check))
            } else if (length(input) == 1 & is.null(names(input))){
              prior$prior[[pri]] <- rep(input, length(to_check))
              to_do <- rep(F, length(to_check))
            } else{
              to_do <- !(names(to_check) %in% names(input))
              to_check[!to_do] <- input[names(input) %in% names(to_check)][names(to_check)[!to_do]]
              prior$prior[[pri]][,] <- diag(to_check)
            }
          } else{
            to_check <- prior$prior[[pri]]
            if(length(input) == length(to_check) & is.null(names(input))){
              prior$prior[[pri]] <- input
              to_do <- rep(F, length(to_check))
            } else if (length(input) == 1 & is.null(names(input))){
              prior$prior[[pri]] <- rep(input, length(to_check))
              to_do <- rep(F, length(to_check))
            } else{
              to_do <- !(names(to_check) %in% names(input))
              to_check[!to_do] <- input[names(input) %in% names(to_check)][names(to_check)[!to_do]]
              prior$prior[[pri]] <- to_check
            }
          }
        } else{
          if(!is.null(dim(prior$prior[[pri]]))){
            to_do <- rep(T, nrow(prior$prior[[pri]]))
          } else{
            to_do <- rep(T, length(prior$prior[[pri]]))
          }
        }
        # Ask the user to manually specify
        if(any(to_do)){
          tmp <- ask_user_prior(prior, pri, to_do, fill_default, group_to_do)
          if(!is.null(dim(prior$prior[[pri]]))){
            input <- diag(prior$prior[[pri]][,])
            input[to_do] <- tmp
            prior$prior[[pri]][,] <- diag(input)
          } else{
            prior$prior[[pri]][to_do] <- tmp
          }
        }
      }
    }
  }
  prior <- get_objects(design = design, type = type, prior = prior$prior, ...)
  return(prior$prior)
}

check_prior <- function(prior, sampled_p_names){
  if (is.null(names(prior$theta_mu_mean)))
    names(prior$theta_mu_mean) <- sampled_p_names
  if(length(prior$theta_mu_mean) != length(sampled_p_names))
    stop("prior mu should be same length as estimated parameters (p_vector)")
  if (!all(sort(names(prior$theta_mu_mean)) == sort(sampled_p_names)))
    stop("theta_mu_mean names not the same as sampled paramter names")
  # Make sure theta_mu_mean has same order as sampled parameters
  pnams <- names(prior$theta_mu_mean)
  prior$theta_mu_mean <- prior$theta_mu_mean[sampled_p_names]
  if(!is.matrix(prior$theta_mu_var)) {
    if(length(prior$theta_mu_var) != length(sampled_p_names))
      stop("prior theta should be same length as estimated parameters (p_vector)")
    # Make sure theta_mu_var has same order as sampled parameters
    names(prior$theta_mu_var) <- pnams
    prior$theta_mu_var <- prior$theta_mu_var[sampled_p_names]
  } else {
    if(nrow(prior$theta_mu_var) != length(sampled_p_names))
      stop("prior theta should have same number of rows as estimated parameters (p_vector)")
    if(ncol(prior$theta_mu_var) != length(sampled_p_names))
      stop("prior theta should have same number of columns as estimated parameters (p_vector)")
    # Make sure theta_mu_var has same order as sampled parameters
    dimnames(prior$theta_mu_var) <- list(pnams,pnams)
    prior$theta_mu_var <- prior$theta_mu_var[sampled_p_names,sampled_p_names]
  }
  return(prior)
}

merge_priors <- function(prior_list){
  out <- prior_list[[1]]
  if(length(prior_list) > 1){
    out_names <- names(out)
    for(j in 2:length(prior_list)){
      for(hyper in out_names){
        if(length(out[[hyper]]) > 1){
          if(length(dim(out[[hyper]])) == 2){
            out[[hyper]] <- adiag(out[[hyper]], prior_list[[j]][[hyper]])
          } else{
            out[[hyper]] <- c(out[[hyper]], prior_list[[j]][[hyper]])
          }
        }
      }
    }
  }
  return(out)
}


ask_user_prior <- function(prior, cur_idx, to_do, fill_default, group_to_do){
  if(!is.null(dim(prior$prior[[cur_idx]]))){
    tmp <- diag(prior$prior[[cur_idx]])[to_do]
  } else{
    tmp <- prior$prior[[cur_idx]][to_do]
  }
  if(!fill_default || group_to_do){
    cat(paste0("Specify ", prior$descriptions[[cur_idx]]," \n"))
    if(length(tmp) == 1){
      cat("Press enter to use default value (", tmp[1], ")")
    } else{
      cat("Press enter to fill with default value (", tmp[1], ")")
    }
    for(i in 1:length(tmp)){ # Loop over the length of the prior object
      name <- ifelse(length(to_do) == 1, "", names(tmp)[i])
      repeat {
        ans <- try(eval(parse(text=readline(paste0(name,": ")))),silent=TRUE)
        if(!is.null(ans)){
          if (any(class(ans) %in% c("warning", "error", "try-error")) || is.na(ans) || !is.numeric(ans)) {
            cat("Must provide a numeric value\n")
          } else {
            tmp[i] <- ans
            break
          }
        } else{
          break
        }
      }
    }
  }
  return(tmp)
}



#' Title
#'
#' @param prior A prior list created with `prior`
#' @param design A design list created with `design`
#' @param do_plot Boolean. If `FALSE` will only return prior samples and omit plotting.
#' @param covariates dataframe/functions as specified by the design
#' @inheritParams plot_pars
#' @param ... Optional arguments that can be passed to get_pars, histogram, plot.default (see par()),
#' or arguments required for the types of models e.g. n_factors for type = "factor"
#' @return An mcmc.list object with prior samples of the selected type
#' @export
#'
#' @examples \donttest{
#' # First define a design for the model
#' design_DDMaE <- design(data = forstmann,model=DDM,
#'                            formula =list(v~0+S,a~E, t0~1, s~1, Z~1, sv~1, SZ~1),
#'                            constants=c(s=log(1)))
#' # Then set up a prior using make_prior
#' p_vector=c(v_Sleft=-2,v_Sright=2,a=log(1),a_Eneutral=log(1.5),a_Eaccuracy=log(2),
#'           t0=log(.2),Z=qnorm(.5),sv=log(.5),SZ=qnorm(.5))
#' psd <- c(v_Sleft=1,v_Sright=1,a=.3,a_Eneutral=.3,a_Eaccuracy=.3,
#'           t0=.4,Z=1,sv=.4,SZ=1)
#' # Here we left the variance prior at default
#' prior_DDMaE <- prior(design_DDMaE,mu_mean=p_vector,mu_sd=psd)
#' # Now we can plot all sorts of (implied) priors
#' plot_prior(prior_DDMaE, design_DDMaE, selection = "mu", N = 1e3)
#' plot_prior(prior_DDMaE, design_DDMaE, selection = "mu", mapped = FALSE, N=1e3)
#' # We can also plot the implied prior on the participant level effects.
#' plot_prior(prior_DDMaE, design_DDMaE, selection = "alpha", col = "green", N = 1e3)
#' }

plot_prior <- function(prior, design, selection = "mu", do_plot = TRUE, covariates = NULL,
                           layout = NA, N = 5e4, ...){
  dots <- add_defaults(list(...), breaks = 30, cut_off = 0.0015, prob = TRUE, by_subject = TRUE, map = TRUE)
  oldpar <- par(no.readonly = TRUE) # code line i
  on.exit(par(oldpar)) # code line i + 1
  if(is.null(design$Ffactors)){
    dots$map <- FALSE
    warning("For joint models, map = TRUE is not yet implemented")
  }
  type <- attr(prior, "type")
  if(type == "single") selection <- "alpha"
  samples <-  get_objects(design = design, prior = prior, type = type, sample_prior = T,
                          selection = selection, N = N, ...)
  MCMC_samples <- do.call(get_pars, c(list(samples, selection = selection, type = type, covariates = covariates), fix_dots(dots, get_pars)))
  if(do_plot){
    for(i in 1:length(MCMC_samples)){
      xlab <- ifelse(is.null(names(MCMC_samples)[i]), selection, names(MCMC_samples)[i])
      if(any(is.na(layout))){
        par(mfrow = coda_setmfrow(Nchains = length(MCMC_samples[[1]]),
                                  Nparms = ncol(MCMC_samples[[1]][[1]]), nplots = 1))
      } else{
        par(mfrow=layout)
      }
      for(j in 1:ncol(MCMC_samples[[i]][[1]])){
        do.call(robust_hist, c(list(MCMC_samples[[i]][[1]][,j], dots$breaks, dots$cut_off, dots$prob),
                               fix_dots_plot(add_defaults(dots, ylab = "Density", xlab = xlab,
                                                          main = colnames(MCMC_samples[[i]][[1]])[j],
                                                          cex.lab = 1.25, cex.main = 1.5))))
      }
    }
  }
  return(invisible(MCMC_samples))
}


