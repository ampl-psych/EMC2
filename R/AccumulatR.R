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


AccumulatR_add_context <- function(dadm){
  model <- attr(dadm, "model")
  if(is.null(model)) return(dadm)
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
