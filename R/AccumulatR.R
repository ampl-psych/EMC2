AccumulatR_model <- function(model_spec){
  model_list <- list(
    c_name = "AccumulatR",
    type = "AccumulatR",
    p_types = c("q" = qnorm(0),  "w" = qnorm(1),
                "t0" = log(0.05), "p1" = log(.5), "p2" = log(1), "p3" = log(1)),
    transform=list(func=c(q = "pnorm", w = "pnorm", t0 = "exp",
                          p1 = "identity",p2 = "exp", p3 = "exp")),
    Ttransform = function(pars,dadm) pars,
    bound=list(minmax=cbind(q=c(0,1), w = c(0, 1), t0=c(0.05,Inf),
                            p1=c(-Inf, Inf), p2 = c(0, Inf), p3 = c(0, Inf)),
               exception=c(q=0,w=1,t0=0.05, p2 = 0, p3 = 0)),
    spec = model_spec
  )
  return(function() {return(model_list)})
}

AccumulatR_expand_data <- function(model_spec, data){
  if(is.null(model_spec)) stop("model_specification needs to be made for AccumulatR models")
  if(nrow(model_spec$component) > 1){
    if(is.null(data$component)){
      warning("The AccumulatR model contains multiple components. Specify the components
              column in your data if you do not want to treat them as mixtures
              marginalized in the likelihood")
      cat("expected components :", model_spec$components$component_id)
    }
  }
  datar <- nest_accumulators(model_spec, data)
  return(datar)
}




AccumulatR_add_context <- function(dadm){
  model <- attr(dadm, "model")
  if(is.null(model)) return(dadm)
  if(model()$type != "AccumulatR") return(dadm)
  nacc <- unique
  data <- dadm[!duplicated(dadm$trial),]
  model_list <- model()
  model_spec <- model_list$spec
  ctx <- build_likelihood_context(model_spec, data)

  context <- list(
    native_ctx = ctx$native_ctx,  # externalptr (the “real” native context)
    structure  = ctx$structure,   # model structure (R list)
    param_layout = ctx$param_layout,
    data_df    = ctx$data_df,     # observed data
    rel_tol    = ctx$rel_tol,
    abs_tol    = ctx$abs_tol,
    max_depth  = ctx$max_depth
  )


  attr(dadm, "AccumulatR_context") <- context
  return(dadm)
}

AccumulatR_check_context <- function(emc){
  if(!emc[[1]]$model()$type == "AccumulatR") return(emc)
  for(i in 1:length(emc)){
    for(j in 1:length(emc[[i]]$data)){
      dat <- emc[[i]]$data[[j]] # Looping over chains, and subjects
      if(is.data.frame(dat)){
        ctx <- ensure_native_ctx(attr(dat, "AccumulatR_context"))
        context <- list(
          native_ctx = ctx$native_ctx,  # externalptr (the “real” native context)
          structure  = ctx$structure,   # model structure (R list)
          param_layout = ctx$param_layout,
          data_df    = ctx$data_df,     # observed data
          rel_tol    = ctx$rel_tol,
          abs_tol    = ctx$abs_tol,
          max_depth  = ctx$max_depth
        )
        attr(dat, "AccumulatR_context") <- context
        emc[[i]]$data[[j]] <- dat
      }
    }
  }
  return(emc)
}


