# Package-level store for loo warnings, cleared at the start of each compare()
# call and accessible afterwards via waic_warnings().
.loo_warn_store <- new.env(parent=emptyenv())
.loo_warn_store$msgs <- character(0)

#' Retrieve loo WAIC Warnings
#'
#' Returns any p_waic warnings issued by the \pkg{loo} package during the most
#' recent call to \code{\link{compare}} or \code{\link{compare_subject}}.
#' Warnings are cleared at the start of each such call.
#'
#' @return A character vector of warning messages, or \code{character(0)} if
#'   none occurred.
#' @export
waic_warnings <- function() .loo_warn_store$msgs

# Evaluates thunk(), muffling all loo warnings during execution and storing
# them in .loo_warn_store so they can be retrieved with waic_warnings().
# Returns structure(result, warned=<logical>).
.waic_call <- function(thunk){
  result <- withCallingHandlers(thunk(), warning = function(w){
    .loo_warn_store$msgs <- c(.loo_warn_store$msgs, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  structure(result, warned=length(.loo_warn_store$msgs) > 0)
}

waic_subject <- function(emc, stage="sample", filter=0, subject){
  # Returns scalar WAIC estimate for a single subject.
  alpha <- get_pars(emc, selection="alpha", stage=stage, filter=filter,
                    by_subject=TRUE, merge_chains=TRUE)
  proposals <- do.call(rbind, alpha[[subject]])   # [n_iter x n_pars]
  ll_mat <- calc_ll_manager_pw(proposals, emc[[1]]$data[[subject]], emc[[1]]$model)
  .waic_call(function() loo::waic(ll_mat)$estimates["waic","Estimate"])
}

waic_pooled <- function(emc, stage="sample", filter=0){
  # Returns scalar WAIC estimate pooled across all subjects.
  alpha <- get_pars(emc, selection="alpha", stage=stage, filter=filter,
                    by_subject=TRUE, merge_chains=TRUE)
  data  <- emc[[1]]$data
  model <- emc[[1]]$model
  ll_list <- lapply(names(alpha), function(sub){
    proposals <- do.call(rbind, alpha[[sub]])     # [n_iter x n_pars]
    calc_ll_manager_pw(proposals, data[[sub]], model)  # [n_iter x n_trials_sub]
  })
  ll_all <- do.call(cbind, ll_list)              # [n_iter x total_trials]
  .waic_call(function() loo::waic(ll_all)$estimates["waic","Estimate"])
}
