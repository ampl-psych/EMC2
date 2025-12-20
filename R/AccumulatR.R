new_accumulatR_model <- function(model_spec, data){
  structure <- finalize_model(model_spec)
  ctx <- build_likelihood_context(structure, data)

  context <- list(
    native_ctx = ctx$native_ctx,  # externalptr (the “real” native context)
    structure  = ctx$structure,   # model structure (R list)
    data_df    = ctx$data_df,     # observed data
    rel_tol    = ctx$rel_tol,
    abs_tol    = ctx$abs_tol,
    max_depth  = ctx$max_depth
  )

  specs <- list(
    c_name = "AccumulatR",
    type = "AccumulatR",
    p_types = c("q" = qnorm(0),  "w" = qnorm(1),
                "t0" = log(0.05), "p1" = log(.5), "p2" = log(1), "p3" = log(1)),
    transform=list(func=c(q = "pnorm", w = "pnorm", t0 = "exp",
                          p1 = "identity",p2 = "exp", p3 = "exp")),
    bound=list(minmax=cbind(q=c(0,1), w = c(0, 1), t0=c(0.05,Inf),
                            p1=c(-Inf, Inf), p2 = c(0, Inf), p3 = c(0, Inf))),
    context = context
  )
  return(specs)
}
