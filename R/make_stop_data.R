#' Simulate Latent Stop and Go Finishing Times
#'
#' For the stop-signal race models [SSEXG()] and [SSRDEX()], simulates data as
#' [make_data()] (parameter vector or matrix) or [predict()] (an `emc` object
#' with posterior samples) do, and adds the latent finishing times that
#' produced each simulated trial: the stop racer's finishing time and the
#' finishing time and identity of the winner of the go race.
#'
#' The latent times come from the same random draws as the simulated `R` and
#' `rt`, so on a go trial with a response `rt == goRT`, on a signal-respond
#' trial `goRT < SSD + SSRT`, and on a successful stop `goRT > SSD + SSRT`.
#' They follow the model as fitted:
#'
#' * `SSRT` is drawn from the ex-Gaussian stop distribution with lower bound
#'   `exgS_lb`, measured from stop-signal onset (`SSRT > exgS_lb`), and is `Inf`
#'   on a trigger failure. For [SSEXG()], go finishing times have lower bound
#'   `exg_lb`.
#' * `goRT` is the finishing time of the fastest go accumulator, measured from
#'   go onset. Stop-triggered accumulators (`lI`) are not part of the go race.
#'   `goRT` is `Inf` (and `goR` `NA`) on a go failure.
#' * Trends on go and/or stop parameters are applied trial by trial; when a
#'   trend covariate depends on behaviour the data are simulated trial by trial
#'   (see [predict()]).
#' * Censoring and truncation (`TC` in [design()], [make_data()] or columns of
#'   the data, see [make_missing()]) apply as for the observed data: trials
#'   removed by truncation are removed with their latent times, and the go
#'   race is censored at the same `LC`/`UC` bounds as `rt`.
#'
#' @param parameters An `emc` object fitted with [SSEXG()] or [SSRDEX()], in
#'   which case simulations use `n_post` posterior draws (as [predict()]), or a
#'   parameter vector or matrix (as [make_data()]).
#' @param design Design list from [design()]; only used when `parameters` is
#'   not an `emc` object.
#' @param n_trials Integer, trials per design cell when `parameters` is not an
#'   `emc` object and `data` is not supplied.
#' @param data Data frame whose design (trials, SSDs) is used for the
#'   simulation. For an `emc` object the default is the fitted data.
#' @param n_post Integer, number of posterior draws (`emc` objects only).
#' @param hyper Logical, simulate from group-level rather than subject-level
#'   draws (`emc` objects only, see [predict()]).
#' @param stat Character, `"random"` (default) uses random posterior draws,
#'   `"mean"` or `"median"` the posterior mean or median (`emc` objects only).
#' @param censor_go Logical, default `TRUE`, censor `goRT` at the `LC`/`UC`
#'   bounds used for `rt`. `FALSE` returns uncensored go-race times.
#' @param stop_on_go_trials Logical, default `FALSE`. If `TRUE`, also draw a
#'   stop finishing time (with trigger failure) on go trials, where no stop
#'   signal was presented; otherwise `SSRT` is `NA` on go trials. The extra
#'   draws do not change the simulated `R` and `rt`.
#' @param n_cores Integer, cores used across posterior draws (`emc` objects
#'   only).
#' @param ... Further arguments passed to [predict()] or [make_data()].
#'
#' @return The simulated data frame (with a `postn` column for `emc` objects)
#'   with added columns:
#'   \describe{
#'     \item{goR}{Winning go accumulator (factor with the response levels);
#'       `NA` on a go failure or, unless `LCresponse`/`UCresponse`, when
#'       `goRT` is censored.}
#'     \item{goRT}{Finishing time of the go-race winner; `Inf` on a go
#'       failure when there is no finite upper censoring bound; `NA` when
#'       censored.}
#'     \item{SSRT}{Stop finishing time from stop-signal onset; `Inf` on a
#'       trigger failure, `NA` on go trials unless `stop_on_go_trials = TRUE`.}
#'     \item{goMissingness}{Censoring of `goRT`, coded as `missingness`:
#'       `NA` observed, 1 below `LC`, 2 above `UC` (only when
#'       `censor_go = TRUE`).}
#'   }
#' @examples
#' design_ss <- design(model = SSEXG,
#'   factors = list(subjects = 1, S = c("left", "right")),
#'   Rlevels = c("left", "right"),
#'   matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR),
#'   formula = list(mu ~ lM, sigma ~ 1, tau ~ 1, muS ~ 1, sigmaS ~ 1, tauS ~ 1,
#'                  gf ~ 1, tf ~ 1))
#' p_vector <- c(mu = log(.6), mu_lMTRUE = log(.8), sigma = log(.05),
#'               tau = log(.2), muS = log(.2), sigmaS = log(.03),
#'               tauS = log(.05), gf = qnorm(.1), tf = qnorm(.1))
#' dat <- make_stop_data(p_vector, design_ss, n_trials = 100,
#'                       functions = list(SSD = make_ssd()))
#' # Stop finishing time distribution (stop trials with a triggered stop racer)
#' summary(dat$SSRT[is.finite(dat$SSRT)])
#' # For a fitted model: make_stop_data(emc, n_post = 100)
#' @export
make_stop_data <- function(parameters, design = NULL, n_trials = NULL, data = NULL,
                           n_post = 50, hyper = FALSE, stat = c("random", "mean", "median")[1],
                           censor_go = TRUE, stop_on_go_trials = FALSE, n_cores = 1, ...)
{
  is_ss <- function(des) isTRUE(des$model()$c_name %in% c("SSEXG", "SSRDEX"))
  if (is(parameters, "emc")) {
    des <- get_design(parameters)
    if (length(des) > 1)
      stop("make_stop_data() does not yet support joint models")
    if (!is_ss(des[[1]]))
      stop("make_stop_data() requires a stop-signal model (SSEXG or SSRDEX)")
    args <- list(parameters, hyper = hyper, n_post = n_post, n_cores = n_cores,
                 stat = stat, latent = TRUE, censor_go = censor_go,
                 stop_on_go_trials = stop_on_go_trials, ...)
    if (!is.null(data)) args$data <- data
    return(do.call(predict, args))
  }
  if (is.null(design)) stop("design must be supplied when parameters is not an emc object")
  if (!is_ss(design))
    stop("make_stop_data() requires a stop-signal model (SSEXG or SSRDEX)")
  make_data(parameters, design = design, n_trials = n_trials, data = data,
            latent = TRUE, censor_go = censor_go,
            stop_on_go_trials = stop_on_go_trials, ...)
}
