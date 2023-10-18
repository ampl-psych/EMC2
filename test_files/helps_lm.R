# rm(list = ls())
# devtools::load_all()
# source("test_files/utils_lm.R")
# library(lme4)
# library(BayesFactor)
# n_subjects <- 6
# n_cond <- 3
# n_trials <- 50
# n_groups <- 2
#
# subject <- rep(1:n_subjects, each = n_trials*n_cond)
# cond <- rep(rep(1:n_cond, n_subjects), each = n_trials)
#
# df <- data.frame(subjects = subject, cond = cond)
# for(i in 2:n_groups){
#   new_df <- df
#   new_df$subjects <- new_df$subjects +n_subjects
#   df <- rbind(df, new_df)
# }
# df$group <- rep(1:n_groups, each = n_trials*n_cond*n_subjects)
# df$subjects <- as.factor(df$subjects)
# df$cond <- as.factor(df$cond)
# df$group <- as.factor(df$group)
# df$R <- as.factor(c(1,2))
# df$S <- as.factor(c(1,2,1,2,1,2,2,1))
# df$rt <- .1 + exp(rnorm(nrow(df), .2, .5))
# hist(df$rt, breaks = 30)


make_design_lm <- function(Flist,matchfun=NULL,constants=NULL,report_p_vector=TRUE, ddata = NULL, model){
  # This will reorder the data for computational purposes! The design matrices will be on the reordered data.
  if(is.null(ddata)){
    stop("For now, ddata has to be specified in type lm")
  }
  all_vars <- unique(unlist(lapply(form, FUN = function(x) all.vars(delete.response(terms(x))))))
  df_covs <- sapply(df, is.factor)
  df_covs <- df_covs[!df_covs]
  Fcovariates <- all_vars[!all_vars %in% c("subjects", "RT", "R") & all_vars %in% df_covs]
  if(length(Fcovariates == 0)) Fcovariates <- NULL
  ddata <- EMC2:::add_accumulators(ddata, matchfun = matchfun, type = model()$type, Fcovariates = Fcovariates)

  # Set up placeholders
  DM_fixed <- vector("list", length(Flist))
  DM_random <- vector("list", length(Flist))
  g_fixed <- vector("list", length(Flist))
  g_random <- vector("list", length(Flist))
  pars_random <- character()
  pars_fixed <- character()
  names(g_random) <- lapply(Flist, FUN = function(x) return(x[[2]]))
  names(g_fixed) <- lapply(Flist, FUN = function(x) return(x[[2]]))
  uniques_in_group <- character()
  names(DM_fixed) <- lapply(Flist, FUN = function(x) return(x[[2]]))
  names(DM_random) <- lapply(Flist, FUN = function(x) return(x[[2]]))
  for(form in Flist){
    name <- form[[2]]
    use <- types_and_baseForm(form)
    trms <- use$terms
    dataTypes <- use$types
    if(length(trms) > 0){ # More than just an intercept
      Xs = lapply(trms, function(trm, data, dataTypes) {
        oneDM(trm, data = ddata, dataTypes = dataTypes)
      }, data = ddata, dataTypes = dataTypes)
      DM <- sepFixRandom(Xs, use$types)
      par_names_curr <- get_parNames(DM$fixed,DM$random, as.character(name))
      colnames(DM$fixed) <- par_names_curr$fixed
      colnames(DM$random) <- par_names_curr$random
      g <- getGMap(DM, use, ddata, constants)
      DM_random[[name]] <- DM$random
      DM_fixed[[name]] <- DM$fixed
      g_random[[name]] <- g$random
      g_fixed[[name]] <-  g$fixed
      uniques_in_group <- c(uniques_in_group, par_names_curr$uniques_in_group)
      pars_fixed <- c(pars_fixed, par_names_curr$fixed)
      pars_random <- c(pars_random, par_names_curr$random)
    } else{
      pars_fixed <- c(pars_fixed, as.character(name)) # Intercept only
      intercept <- as.matrix(rep(1, nrow(ddata)))
      colnames(intercept) <- as.character(name)
      DM_fixed[[name]] <- intercept
      if(!as.character(name) %in% names(constants)){
        g_curr <- 1
        names(g_curr) <- as.character(name)
        g_fixed[[name]] <- g_curr
      }
    }
  }
  if (report_p_vector) {
    cat("fixed : \n")
    cat(pars_fixed, "\n\n")
    cat("random per subject : \n")
    cat(uniques_in_group)
  }
  design <- list(
    DM_fixed = DM_fixed,
    DM_random = DM_random,
    constants = constants,
    model = model,
    Fcovariates = Fcovariates,
    per_subject = uniques_in_group
  )
  attr(design, "p_vector_fixed") <- pars_fixed[!pars_fixed %in% names(constants)]
  attr(design, "p_vector_random") <- pars_random[!pars_random %in% names(constants)]
  attr(design, "g_fixed") <- g_fixed
  attr(design, "g_random") <- g_random
  return(design)
}


# form <- list(v~cond+lM + (1|subjects), B~cond+ (cond|subjects), A ~ 1,
#             t0~ 1 + (1|group), sv ~ 1)
#
# design <- make_design_lm(form, ddata = df, model = lbaB, matchfun=function(d)d$S==d$lR,
#                        constants = c(sv = 1))

# names_fixed <- attr(test, "p_vector_fixed")
# names_random <- attr(test, "p_vector_random")
#
# REs <- rnorm(length(names_random), sd = .2)
# names(REs) <- names_random
# fixed <- c(1, .5, .3, 1, log(1.5), .2 , .1, log(.4), log(.2))
# names(fixed) <- names_fixed
# p_vector <- c(fixed, REs, test$constants)
#
# DM <- mapply(cbind, test$DM_fixed, test$DM_random)
# pars <- matrix(0, nrow(DM[[1]]), length(DM))
# colnames(pars) <- names(DM)
# for(i in 1:length(DM)){
#   pars[,i] <- as.vector(DM[[i]] %*% p_vector[names(p_vector) %in% colnames(DM[[i]])])
# }
filtered_samples_lm <- function(sampler, filter){
  out <- list(
    random = sampler$samples$random[,filter, drop = F],
    fixed = sampler$samples$fixed[,filter, drop = F],
    g_map_random = sampler$g_map_random,
    g_fixed = sampler$samples$g_fixed[,filter, drop = F],
    g_random = sampler$samples$g_random[,filter, drop = F],
    iteration = length(filter)
  )
}

get_conditionals_lm_between <- function(samples, iteration = NULL){
  iteration <- ifelse(is.null(iteration), samples$iteration, iteration)
  between_idx <- !grepl("subjects", rownames(samples$random))
  g_random_idx <- unique(rownames(samples$g_random)[samples$g_map_random][between_idx])
  all_samples <- rbind(samples$fixed, samples$random[between_idx,],samples$g_fixed, samples$g_random[g_random_idx,])
  mu_tilde <- rowMeans(all_samples)
  var_tilde <- stats::cov(t(all_samples))
  n_pars <- nrow(samples$fixed) + sum(between_idx)
  condmvn <- condMVN(mean = mu_tilde, sigma = var_tilde,
                     dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                     X.given = c(samples$random[between_idx,iteration], samples$g_fixed[,iteration], samples$g_random[g_random_idx,iteration]))
  return(list(eff_mu = condmvn$condMean, eff_var = condmvn$condVar))
}

get_conditionals_lm <- function(s, samples, iteration = NULL){
  iteration <- ifelse(is.null(iteration), samples$iteration, iteration)
  sub_idx <- get_sub_idx(s, rownames(samples$random))
  within_idx <- grepl("subjects", rownames(samples$random))
  g_random_idx <- unique(rownames(samples$g_random)[samples$g_map_random][within_idx])
  all_samples <- rbind(samples$random[sub_idx,],samples$g_random[g_random_idx,])
  mu_tilde <- rowMeans(all_samples)
  var_tilde <- stats::cov(t(all_samples))
  n_pars <- sum(sub_idx)
  condmvn <- condMVN(mean = mu_tilde, sigma = var_tilde,
                     dependent.ind = 1:n_pars, given.ind = (n_pars + 1):length(mu_tilde),
                     X.given = samples$g_random[g_random_idx,iteration])
  return(list(eff_mu = condmvn$condMean, eff_var = condmvn$condVar))
}



make_samplers_lm <- function(data_list,design_list, model_list =NULL,
                            n_chains=3,compress=TRUE,rt_resolution=0.02,
                            prior = NULL, cov_type = "standard")

{
  if (is(data_list, "data.frame")) data_list <- list(data_list)
  if (!is.null(names(design_list)[1]) && names(design_list)[1]=="DM_fixed")
    design_list <- list(design_list)
  if (length(design_list)!=length(data_list))
    design_list <- rep(design_list,length(data_list))
  if (is.null(model_list)) model_list <- lapply(design_list,function(x){x$model})
  if (any(unlist(lapply(model_list,is.null))))
    stop("Must supply model_list if model is not in all design_list components")
  if (!is.null(names(model_list)[1]) && names(model_list)[1]=="type")
    model_list <- list(model_list)
  if (length(model_list)!=length(data_list))
    model_list <- rep(model_list,length(data_list))
  data_list <- lapply(data_list,function(d){
    if (!is.factor(d$subjects)) d$subjects <- factor(d$subjects)
    d <- d[order(d$subjects),]
    EMC2:::add_trials(d)
  })
  dadm_list <- vector(mode="list",length=length(data_list))
  rt_resolution <- rep(rt_resolution,length.out=length(data_list))
  for (i in 1:length(dadm_list)) {
    message("Processing data set ",i)
    dadm_list[[i]] <- design_model(data=data_list[[i]],design=design_list[[i]],
                                   compress=compress,model=model_list[[i]],rt_resolution=rt_resolution[i])
  }

  attr(dadm_list[[1]], "prior") <- prior

  # if(!is.null(subject_covariates)) attr(dadm_list, "subject_covariates") <- subject_covariates

  out <- pmwgs_lm(dadm_list, type = "standard")
  # replicate chains
  dadm_lists <- rep(list(out),n_chains)
  # For post predict
  attr(dadm_lists,"data_list") <- data_list
  attr(dadm_lists,"design_list") <- design_list
  attr(dadm_lists,"model_list") <- model_list
  return(dadm_lists)
}
# devtools::load_all()
#
#
# samplers <- make_samplers_lm(df, design)
# samplers <- run_samplers(samplers, stage = "preburn", cores_for_chains = 1, cores_per_chain = 14,iter = 100,
#                          verbose = T)
#
# test <- auto_burn(samplers, cores_for_chains = 1, cores_per_chain = 14,
#                          verbose = T)

#
# m1 <- lmer(form, data)
#
# getME(m1, "Z")


