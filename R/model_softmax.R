
pSOFTMAX <- function(trials, pars)
  # probability between lt and ut
{
  V <- pars[,'x']*pars[,'beta']
  # Compute max(V) within each trial (for numerical stability)
  Vmax <- ave(V, trials, FUN = max)
  expV <- exp(V - Vmax)

  # Trial-wise denominator
  denom <- ave(expV, trials, FUN = sum)

  # Probabilities
  expV / denom
}


rSOFTMAX <- function(lR,pars,p_types=c("x","beta"))
  # lR is an empty latent response factor lR with one level for response.
  # pars is a matrix of corresponding parameter values named as in p_types
  # pars must be sorted so accumulators and parameter for each trial are in
  # contiguous rows.

{
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ",paste(p_types,collapse = " "))
  nr <- length(levels(lR)) # Number of responses
  n <- dim(pars)[1]/nr     # Number of simulated trials

  # trial index for each accumulator row
  trial_id <- rep(seq_len(n), each=nr)
  p <- pSOFTMAX(trial_id, pars)

  # --- sample a response for each trial (vectorized using cumulative sums) -----
  # generate uniform random variable for each trial
  u <- runif(n)
  u_expanded <- u[trial_id]   # match accumulator rows

  # 3. compute cumulative probability per trial
  #    Using cumsum grouped with diff() trick (way faster than ave)
  cp <- cumsum(p)

  # reset cumulative sums at group boundaries:
  # Find where new trials begin
  idx <- c(TRUE, diff(trial_id) != 0)
  cp[idx] <- p[idx]          # restart cumsum at trial start

  # 4. Determine first accumulator where cp >= u
  is_hit <- cp >= u_expanded   # logical matrix "hit" per accumulator row

  # For each trial: we want the FIRST hit
  # Convert trial structure to a 2D matrix: nr rows × n trials
  hit_mat <- matrix(is_hit, nrow = nr, ncol = n)

  # First TRUE per column (trial); returns NA if none
  chosen <- apply(hit_mat, 2, function(x) {
    w <- which(x)
    if (length(w) == 0) return(nr)   # fallback (rare)
    w[1]
  })


  # # cumulative probability per trial
  # cp <- ave(p, trial_id, FUN=cumsum)
  #
  # # choose one per trial
  # chosen <- tapply(cp >= u_expanded, trial_id, function(x) {
  #   w <- which(x)
  #   if (length(w) == 0) return(nr)  # fallback
  #   w[1]
  # })

  # chosen <- ave(cp >= u_expanded, trial_id, FUN=function(x) {
  #   idx <- which(x)[1]
  #   if (is.na(idx)) idx <- nr  # safety fallback
  #   rep(idx, length(x))
  # })
  #
  # chosen <- chosen_per_trial[seq(1, n, nr)]  # pick one per trial
  # --- build response factor ---------------------------------------------------
  R_levels <- levels(lR)
  R <- factor(R_levels[chosen], levels=R_levels)

  # return same structure as rPROBIT
  data.frame(R = R, rt = NA_real_)
}

#' Softmax Discrete Choice Model
#'
#' A multinomial (or binary) discrete-choice model using the softmax
#' (Luce/Logit) rule. The probability of choosing option \eqn{i} in trial
#' \eqn{t} is
#'
#' \deqn{
#'   P(R_t = i \mid x_{ti}, \beta_t)
#'      = \frac{\exp(\beta_t \, x_{ti})}
#'             {\sum_{j} \exp(\beta_t \, x_{tj})}.
#' }
#'
#' The model operates on *long-format* accumulator data, where each trial
#' consists of a block of rows—one per response option (accumulator). No
#' response time is modeled; \code{rt} must be \code{NA} in the input data.
#'
#' @section Parameters:
#' The softmax model uses the following parameters:
#'
#' \describe{
#'   \item{\code{x}}{A continuous latent Q-value or score for each response
#'        option. Unbounded.}
#'
#'   \item{\code{beta}}{Inverse temperature (precision) parameter. Must be
#'        non-negative. Larger values produce more deterministic choice.}
#' }
#'
#' Returned sampled data consist of a response factor \code{R} and an
#' \code{rt} column always set to \code{NA}.
#'
#' @return A model list with all the necessary functions to sample
#'
#' @examples
#' soft <- design(
#'   Rlevels = c("left","right"),
#'   factors = list(subjects = 1, S = c("left","right")),
#'   formula = list(x ~ 0 + S, beta ~ 1),
#'   matchfun = function(d) d$S == d$lR,
#'   constants = c(beta = log(1)),
#'   model = softmax
#' )
#'
#' # Sample parameter vector
#' p_vector <- sampled_pars(soft)
#'
#' @export
softmax <- function(){
  list(
    type="SOFTMAX",
    c_name = "SOFTMAX",
    p_types=c("x" = 0,"beta" = log(1)),
    # Trial dependent parameter transform
    transform=list(func=c(x = "identity", beta = "exp")),
    bound=list(minmax=cbind(x=c(-Inf,Inf),beta = c(0, Inf))),
    Ttransform = function(pars,dadm) {
      pars
    },
    # Random function for discrete choices
    rfun=function(data=NULL,pars) {
      rSOFTMAX(data$lR,pars)
    },
    # probability of choice between lower and upper thresholds (lt & ut)
    pfun=function(trials,pars) pSOFTMAX(trials,pars),
    # Likelihood, lb is lower bound threshold for first response
    log_likelihood=function(pars,dadm,model,min_ll=log(1e-10)){
      log_likelihood_softmax(pars=pars, dadm = dadm, model = model, min_ll = min_ll)
    })
}
