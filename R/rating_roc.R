# ROC goodness-of-fit displays for choice + rating (+ RT) data, after
# Weidemann & Kahana (2016) as used by Heathcote et al. (MTR confidence paper;
# its FirstSubmitPrep/functions.R plotCL_ROC, plotCL_RROC, plotL_CRROC): the
# confidence ROC (C-ROC) and zROC, the latency ROC (L-ROC) from RT quantile bins,
# the response-contingent ROCs (C-RROC, L-RROC, L-RROC within a confidence
# level), their areas and the random ROC through the classification point.
# Ported from rating-work/mtlba/stageM8g/R/roc_functions.R; the arithmetic is
# kept operation for operation so the values equal that code's (and the
# paper's stored values) exactly.

# ---- curves and areas ------------------------------------------------------

# Cumulative proportions over ordered levels 1 .. L for signal and noise
# trials; L + 1 points starting at (0, 0). No trials gives a flat 0 curve.
roc_cum <- function(level, sig, L) {
  cnt_s <- tabulate(level[sig], L); cnt_n <- tabulate(level[!sig], L)
  div <- function(v) if (sum(v) > 0) cumsum(v) / sum(v) else rep(0, L)
  data.frame(k = 0:L, FAR = c(0, div(cnt_n)), HR = c(0, div(cnt_s)))
}

# Trapezoid area under a curve.
roc_auc <- function(FAR, HR) sum(diff(FAR) * (HR[-length(HR)] + HR[-1]) / 2)

# Area under the random ROC: (0, 0), the point after k_class levels, (1, 1).
roc_random_auc <- function(FAR, HR, k_class) {
  f <- FAR[k_class + 1]; h <- HR[k_class + 1]
  f * h / 2 + (1 - f) * (h + 1) / 2
}

# RT quantile bins within one response: 1 = fastest .. n_quantiles - 1 = slowest,
# cut at the n_quantiles quantiles (type 7) of these RTs.
roc_rt_bins <- function(rt, n_quantiles) {
  if (!length(rt)) return(integer(0))
  q <- stats::quantile(rt, probs = seq(0, 1, length.out = n_quantiles), names = FALSE)
  findInterval(rt, q[-c(1, length(q))]) + 1L
}

# ---- preparation -----------------------------------------------------------

# Rows kept by bounds (strict rule as make_missing: LC <= rt <= UC). bounds is
# c(LC, UC) or a data frame with LC and UC (and subjects when it has more than
# one row).
roc_in_bounds <- function(x, bounds) {
  if (is.numeric(bounds) && !is.data.frame(bounds)) {
    if (length(bounds) != 2) stop("bounds must be c(LC, UC) or a data frame with columns LC and UC")
    return(x$rt >= bounds[1] & x$rt <= bounds[2])
  }
  if (!is.data.frame(bounds) || !all(c("LC", "UC") %in% names(bounds)))
    stop("bounds must be c(LC, UC) or a data frame with columns LC and UC")
  if (nrow(bounds) == 1 && is.null(bounds$subjects))
    return(x$rt >= bounds$LC & x$rt <= bounds$UC)
  if (is.null(bounds$subjects) || is.null(x$subjects))
    stop("bounds with more than one row need a subjects column, matched to the data's subjects")
  i <- match(as.character(x$subjects), as.character(bounds$subjects))
  if (anyNA(i)) stop("bounds has no row for subjects ",
                     paste(unique(x$subjects[is.na(i)]), collapse = ", "))
  x$rt >= bounds$LC[i] & x$rt <= bounds$UC[i]
}

# Trials in a common form. ref = the set-up of the data, so predictions get the
# same K, signal response and cells. Rows with R or RR NA (censored) drop out;
# with need_rt or bounds also rows without a finite rt; the missingness column
# is not used (predict() can return a stale one).
roc_setup <- function(input, signal_factor, noise, signal_response = NULL, by = NULL,
                      bounds = NULL, n_ratings = NULL, dont_know = FALSE, need_rt = FALSE,
                      ref = NULL) {
  x <- rating_input(input)
  if (length(signal_factor) != 1 || !signal_factor %in% names(x))
    stop("signal_factor must name one column of the data")
  if (nlevels(x$R) != 2) stop("ROCs need exactly two response levels in R")
  if (!is.null(by) && !all(by %in% names(x)))
    stop("by factor(s) not in the data: ", paste(setdiff(by, names(x)), collapse = ", "))
  keep <- !is.na(x$R) & !is.na(x$RR)
  if (need_rt || !is.null(bounds)) keep <- keep & is.finite(x$rt)
  if (!is.null(bounds)) keep <- keep & roc_in_bounds(x, bounds)
  x <- x[keep, , drop = FALSE]
  stim <- as.character(x[[signal_factor]])
  if (is.null(ref)) {
    if (!any(stim == noise)) stop("no trials with ", signal_factor, " == ", noise)
    sig_lev <- setdiff(levels(factor(x[[signal_factor]])), noise)
    if (!length(sig_lev)) stop("no signal trials (", signal_factor, " other than ", noise, ")")
    if (is.null(signal_response))
      signal_response <- if (length(sig_lev) == 1 && sig_lev %in% levels(x$R)) sig_lev else levels(x$R)[2]
    if (!signal_response %in% levels(x$R)) stop("signal_response must be a level of R")
    K <- if (is.null(n_ratings)) max(x$RR) else n_ratings
    cells <- if (is.null(by)) "all" else
      levels(interaction(x[, by, drop = FALSE], drop = TRUE, sep = " ", lex.order = TRUE))
  } else {
    sig_lev <- ref$signal_levels; signal_response <- ref$signal_response
    K <- ref$K; cells <- ref$cells
  }
  fold <- rating_fold(x$R, x$RR, K, dont_know)
  idx <- as.integer(fold)
  if (anyNA(idx)) stop("RR outside 1 .. ", K, " (set n_ratings)")
  L <- nlevels(fold)
  sig_second <- signal_response == levels(x$R)[2]
  cell <- if (is.null(by)) rep("all", nrow(x)) else
    as.character(interaction(x[, by, drop = FALSE], drop = TRUE, sep = " ", lex.order = TRUE))
  by_values <- NULL
  if (!is.null(by)) {
    first <- !duplicated(cell)
    by_values <- x[first, by, drop = FALSE]
    by_values <- by_values[match(cells, cell[first]), , drop = FALSE]
    rownames(by_values) <- NULL
  }
  list(tab = data.frame(cell = factor(cell, levels = cells), stim = stim,
                        resp_sig = as.character(x$R) == signal_response,
                        conf = as.integer(x$RR),
                        # C-ROC order: most confident signal response first
                        lev = if (sig_second) L + 1L - idx else idx,
                        rt = x$rt, postn = if (is.null(x$postn)) 1L else x$postn,
                        stringsAsFactors = FALSE),
       K = K, L = L, dont_know = dont_know, noise = noise, signal_factor = signal_factor,
       signal_levels = sig_lev, signal_response = signal_response,
       noise_response = setdiff(levels(x$R), signal_response),
       response_labels = if (sig_second) rev(levels(fold)) else levels(fold),
       cells = cells, by = by, by_values = by_values)
}

# Apply fun to every postn x cell x signal level (each signal level against the
# noise trials). fun(tab) returns list(stats = named numeric, <table> = data
# frame, ...). Returns list(stats = wide data frame, <table> = long data frame).
roc_apply <- function(s, fun) {
  tab <- s$tab
  groups <- split(seq_len(nrow(tab)), list(tab$postn, tab$cell), drop = TRUE, sep = "\r")
  out <- list()
  for (g in groups) {
    tg <- tab[g, , drop = FALSE]
    for (sl in s$signal_levels) {
      w <- tg$stim == s$noise | tg$stim == sl
      if (!any(w)) next
      tw <- tg[w, , drop = FALSE]
      tw$sig <- tw$stim != s$noise
      r <- fun(tw)
      key <- data.frame(postn = tw$postn[1], cell = as.character(tw$cell[1]), signal = sl,
                        stringsAsFactors = FALSE)
      out[[length(out) + 1]] <- c(list(stats = cbind(key, n = nrow(tw), t(r$stats))),
                                  lapply(r[names(r) != "stats"], function(d) cbind(key, d, row.names = NULL)))
    }
  }
  nms <- names(out[[1]])
  stats::setNames(lapply(nms, function(nm) {
    d <- do.call(rbind, lapply(out, `[[`, nm)); rownames(d) <- NULL; d
  }), nms)
}

# The result object: data tables (by columns instead of the internal cell key),
# the long statistics table with the predictive summary, and the per-draw tables.
roc_result <- function(type, s, dat, pred, probs, ...) {
  add_by <- function(d, keep_postn = FALSE) {
    if (is.null(d)) return(NULL)
    lead <- if (is.null(s$by)) NULL else s$by_values[match(d$cell, s$cells), , drop = FALSE]
    rest <- d[, setdiff(names(d), c("cell", if (!keep_postn) "postn")), drop = FALSE]
    out <- if (is.null(lead)) rest else cbind(lead, rest)
    rownames(out) <- NULL
    out
  }
  stat_names <- setdiff(names(dat$stats), c("postn", "cell", "signal", "n"))
  long <- do.call(rbind, lapply(seq_len(nrow(dat$stats)), function(i)
    data.frame(cell = dat$stats$cell[i], signal = dat$stats$signal[i], stat = stat_names,
               data = unlist(dat$stats[i, stat_names]), n = dat$stats$n[i],
               stringsAsFactors = FALSE)))
  if (!is.null(pred)) {
    ps <- pred$stats
    pkey <- paste(ps$cell, ps$signal, sep = "\r")
    summ <- lapply(seq_len(nrow(long)), function(i) {
      v <- ps[[long$stat[i]]][pkey == paste(long$cell[i], long$signal[i], sep = "\r")]
      c(pred_median = stats::median(v, na.rm = TRUE),
        pred_lo = stats::quantile(v, probs[1], na.rm = TRUE, names = FALSE),
        pred_hi = stats::quantile(v, probs[2], na.rm = TRUE, names = FALSE),
        n_draws = length(v))
    })
    summ <- as.data.frame(do.call(rbind, summ))
    long$pred_median <- summ$pred_median; long$pred_lo <- summ$pred_lo; long$pred_hi <- summ$pred_hi
    long$inside <- long$data >= long$pred_lo & long$data <= long$pred_hi
    long$n_draws <- summ$n_draws
  }
  rownames(long) <- NULL
  tables <- setdiff(names(dat), "stats")
  res <- c(stats::setNames(lapply(dat[tables], add_by), tables),
           list(stats = add_by(long),
                pred = if (!is.null(pred))
                  stats::setNames(lapply(pred, add_by, keep_postn = TRUE), names(pred))))
  res <- c(res, list(type = type, K = s$K, dont_know = s$dont_know, signal_factor = s$signal_factor,
                     noise = s$noise, signal_levels = s$signal_levels,
                     signal_response = s$signal_response, noise_response = s$noise_response,
                     by = s$by, probs = probs), list(...))
  class(res) <- "emc_roc"
  res
}

roc_run <- function(type, data, signal_factor, noise, signal_response, post_predict, by,
                    bounds, n_ratings, dont_know, probs, need_rt, fun, ...) {
  dk <- rating_dont_know(data, dont_know)
  s <- roc_setup(data, signal_factor, noise, signal_response, by, bounds, n_ratings, dk, need_rt)
  if (dk && type != "croc") stop(type, "() does not handle a shared don't-know response")
  f <- function(tw) fun(tw, s)
  dat <- roc_apply(s, f)
  pred <- NULL
  if (!is.null(post_predict)) {
    if (is.null(post_predict$postn)) stop("post_predict needs a postn column (as from predict())")
    sp <- roc_setup(post_predict, signal_factor, noise, by = by, bounds = bounds, dont_know = dk,
                    need_rt = need_rt, ref = s)
    pred <- roc_apply(sp, f)
  }
  roc_result(type, s, dat, pred, probs, ...)
}

# ---- statistics per data set -------------------------------------------------

croc_one <- function(tw, s) {
  K <- s$K; L <- s$L
  cv <- roc_cum(tw$lev, tw$sig, L)
  cv$response <- c(NA, s$response_labels)
  cAUC <- roc_auc(cv$FAR, cv$HR)
  rAUC <- if (s$dont_know) NA_real_ else roc_random_auc(cv$FAR, cv$HR, K)
  # zROC on the L - 1 interior points; least-squares slope and intercept,
  # curvature = the quadratic coefficient
  zF <- stats::qnorm(cv$FAR[-c(1, L + 1)]); zH <- stats::qnorm(cv$HR[-c(1, L + 1)])
  ok <- is.finite(zF) & is.finite(zH)
  zfit <- if (sum(ok) >= 2) stats::coef(stats::lm(zH[ok] ~ zF[ok])) else c(NA, NA)
  zcurv <- if (sum(ok) >= 3) stats::coef(stats::lm(zH[ok] ~ zF[ok] + I(zF[ok]^2)))[3] else NA
  list(stats = c(cAUC = cAUC, rAUC = rAUC, cRes = cAUC - rAUC, zSlope = unname(zfit[2]),
                 zInt = unname(zfit[1]), zCurv = unname(zcurv)),
       curve = cv,
       zroc = data.frame(k = seq_len(L - 1), zF = ifelse(is.finite(zF), zF, NA),
                         zH = ifelse(is.finite(zH), zH, NA)))
}

lroc_one <- function(tw, s, nq) {
  rs <- tw$resp_sig
  b <- integer(nrow(tw)); b[rs] <- roc_rt_bins(tw$rt[rs], nq); b[!rs] <- roc_rt_bins(tw$rt[!rs], nq)
  # signal responses fastest -> slowest, then noise responses slowest -> fastest
  lev <- ifelse(rs, b, (nq - 1L) + (nq - b))
  cv <- roc_cum(lev, tw$sig, 2 * (nq - 1))
  lAUC <- roc_auc(cv$FAR, cv$HR); lrAUC <- roc_random_auc(cv$FAR, cv$HR, nq - 1)
  list(stats = c(lAUC = lAUC, lrAUC = lrAUC, lRes = lAUC - lrAUC), curve = cv)
}

rroc_one <- function(tw, s, type, conf, nq) {
  rs <- tw$resp_sig; K <- s$K
  pair <- function(cs, cn, suffix) {
    curve <- rbind(cbind(response = "signal", cs), cbind(response = "noise", cn))
    if (!is.null(suffix)) curve <- cbind(conf = suffix, curve)
    list(curve = curve, auc = c(roc_auc(cs$FAR, cs$HR), roc_auc(cn$FAR, cn$HR)))
  }
  if (type == "confidence") {
    # signal responses high -> low confidence, noise responses low -> high
    p <- pair(roc_cum((K - tw$conf + 1L)[rs], tw$sig[rs], K), roc_cum(tw$conf[!rs], tw$sig[!rs], K), NULL)
    return(list(stats = c(cAUC_sig = p$auc[1], cAUC_noise = p$auc[2]), curve = p$curve))
  }
  # latency: signal responses fastest first, noise responses slowest first;
  # within a confidence level the bins are recomputed within that level
  one <- function(ws, wn, suffix) {
    bs <- roc_rt_bins(tw$rt[ws], nq); bn <- roc_rt_bins(tw$rt[wn], nq)
    pair(roc_cum(bs, tw$sig[ws], nq - 1), roc_cum(nq - bn, tw$sig[wn], nq - 1), suffix)
  }
  if (is.null(conf)) {
    p <- one(rs, !rs, NULL)
    return(list(stats = c(lAUC_sig = p$auc[1], lAUC_noise = p$auc[2]), curve = p$curve))
  }
  ps <- lapply(conf, function(cc) one(rs & tw$conf == cc, !rs & tw$conf == cc, cc))
  st <- unlist(lapply(seq_along(conf), function(i)
    stats::setNames(ps[[i]]$auc, paste0(c("lAUC_sig_c", "lAUC_noise_c"), conf[i]))))
  st <- st[c(paste0("lAUC_sig_c", conf), paste0("lAUC_noise_c", conf))]
  list(stats = st, curve = do.call(rbind, lapply(ps, `[[`, "curve")))
}

# ---- exported functions --------------------------------------------------------

#' Confidence ROC, Latency ROC and Response-Contingent ROCs for Rating Data
#'
#' ROC analyses of two-choice data with a confidence rating, after Weidemann and
#' Kahana (2016), with their areas (AUC) and, for posterior predictives, the
#' same per posterior draw with a median and interval. `croc()` gives the
#' confidence ROC (C-ROC) and the zROC, `lroc()` the latency ROC (L-ROC), from
#' response-time quantiles, and `rroc()` the response-contingent ROCs (C-RROC,
#' L-RROC, optionally within one confidence level).
#'
#' The conventions are those of Weidemann and Kahana (2016) as used by
#' Reynolds, Osth, Kvam and Heathcote for the multiple-threshold race (MTR;
#' their Figure 2 and analysis code); AUCs are comparable with those papers
#' only under them.
#'
#' * **Signal and noise.** Trials with `signal_factor == noise` are noise
#'   trials and every other level is a signal level. Each signal level is
#'   analysed against the noise trials (one curve per level). The signal
#'   response is `signal_response`: by default the level of `R` named like the
#'   (single) signal level, otherwise the second level of `R`.
#' * **C-ROC.** Hit rate (HR, signal trials) against false-alarm rate (FAR,
#'   noise trials), cumulated from the most confident signal response down to
#'   the least confident, then from the least confident noise response up to
#'   the most confident: `2K + 1` points from (0, 0) to (1, 1) for `K` ratings
#'   (column `response` names the response added at each point). AUC by
#'   trapezoids. The random ROC
#'   joins (0, 0), the classification point (point `K`, after all signal
#'   responses) and (1, 1); its area `rAUC` depends on accuracy alone and the
#'   confidence resolution is `cRes = cAUC - rAUC`.
#' * **zROC.** The probits of the `2K - 1` interior points (NA where a
#'   proportion is 0 or 1), with the least-squares slope `zSlope` and intercept
#'   `zInt` of z(HR) on z(FAR), and `zCurv`, the quadratic coefficient of a
#'   least-squares quadratic (negative = convex, i.e. bent towards the lower
#'   right). [zroc()] gives the same points in the opposite order.
#' * **L-ROC.** Rating is replaced by RT: the RTs of each response are cut at
#'   their own `n_quantiles` quantiles (type 7, pooled over signal and noise
#'   trials), giving `n_quantiles - 1` bins per response. The curve runs
#'   through the signal responses from the fastest bin to the slowest, then
#'   through the noise responses from the slowest to the fastest, `2
#'   n_quantiles - 1` points; `lAUC`, the random ROC through point
#'   `n_quantiles - 1` (`lrAUC`) and the latency resolution `lRes`.
#' * **Response-contingent ROCs.** One response at a time, HR and FAR
#'   cumulated within that response: confidence (`cAUC_sig`, `cAUC_noise`;
#'   signal responses from high to low confidence, noise responses from low to
#'   high, `K + 1` points) or latency (`lAUC_sig`, `lAUC_noise`; signal
#'   responses fastest first, noise responses slowest first, `n_quantiles`
#'   points). With `conf`, the latency version within each given confidence
#'   level (`lAUC_sig_c<k>`, `lAUC_noise_c<k>`), with the RT bins recomputed
#'   within the level. A response's AUC is below .5 when its slow (signal) or
#'   fast (noise) trials are the more accurate.
#' * **Pooling and bounds.** Without `by`, trials are pooled over everything
#'   else, the subject factor included (as [zroc()]); `by = "subjects"` (or
#'   any factors) gives one set per level. Censored trials (`R` or `RR` NA, as
#'   coded by [make_missing()]) are dropped, and for the latency ROCs also
#'   trials without a finite RT. `bounds` keeps only trials with `LC <= rt <=
#'   UC` (the rule of [make_missing()]) in the data and the predictions alike,
#'   for a comparison within a censoring window. Quantile bins are taken from
#'   the trials pooled within a `by` level (and within a draw for the
#'   predictions).
#'
#' A shared don't-know response (`dont_know`) is placed between the signal and
#' noise responses of the C-ROC (`2K` points), which then has no classification
#' point (`rAUC` and `cRes` NA); `lroc()` and `rroc()` do not handle it.
#'
#' @param data A data frame with columns `R` (two levels), `RR` (rating, 1 =
#'   lowest confidence, `K` = highest), `rt` and `signal_factor`, or an emc
#'   object (its data).
#' @param signal_factor Name of the stimulus factor that separates noise from
#'   signal trials.
#' @param noise Level of `signal_factor` that is noise (all others are signal).
#' @param signal_response Level of `R` that is the signal response (see
#'   Details for the default).
#' @param post_predict Optional posterior predictives (a data frame with a
#'   `postn` column, as from [predict()]): the same curves and statistics per
#'   draw, summarised by the median and the `probs` interval.
#' @param by Optional factor name(s), e.g. `"subjects"`: one analysis per
#'   level. Default: pooled.
#' @param bounds Optional censoring window: `c(LC, UC)`, or a data frame with
#'   columns `LC`, `UC` and (for several rows) `subjects`.
#' @param n_ratings Number of ratings `K` (default: the largest `RR` in `data`).
#' @param dont_know Whether rating 1 is a shared don't-know response (see
#'   [rating_summary()]).
#' @param probs Probabilities of the predictive interval.
#' @return An object of class `emc_roc`, a list with
#'   * `curve`: the data's curve(s): the `by` columns, `signal` (signal level),
#'     `k`, `FAR`, `HR` (and `response` for `croc()`; `response`
#'     = `"signal"` / `"noise"` response, and `conf`, for `rroc()`);
#'   * `zroc` (`croc()` only): `k`, `zF`, `zH`;
#'   * `stats`: one row per statistic (and `by` level and signal level): the
#'     data value `data`, the number of trials `n` and, with `post_predict`,
#'     `pred_median`, `pred_lo`, `pred_hi`, whether the data value lies inside
#'     the interval (`inside`) and `n_draws`;
#'   * `pred`: with `post_predict`, the same tables per draw (`postn`), the
#'     statistics in wide form;
#'   * the settings (`type`, `K`, `signal_response`, `by`, ...).
#' @references Weidemann, C. T., & Kahana, M. J. (2016). Assessing
#'   recognition memory using confidence ratings and response times. *Royal
#'   Society Open Science, 3*, 150670.
#'
#'   Reynolds, A., Osth, A. F., Kvam, P. D., & Heathcote, A. Measuring and
#'   modeling the relationships among speed, confidence, and accuracy.
#'   Manuscript.
#' @seealso [plot_croc()], [plot_roc_auc()], [zroc()].
#' @examples
#' # simulated old/new recognition data with three confidence levels
#' sim <- function(n) {
#'   S <- factor(rep(c("new", "old"), each = n / 2))
#'   e <- rnorm(n, ifelse(S == "old", 1, 0), ifelse(S == "old", 1.25, 1))
#'   R <- factor(ifelse(e > .5, "old", "new"), levels = c("new", "old"))
#'   RR <- as.integer(cut(abs(e - .5), c(0, .4, 1, Inf)))
#'   data.frame(S = S, R = R, RR = RR, rt = .3 + rexp(n, 1 + 2 * abs(e - .5)))
#' }
#' set.seed(1)
#' dat <- sim(1000)
#' pp <- do.call(rbind, lapply(1:20, function(i) cbind(postn = i, sim(1000))))
#' x <- croc(dat, "S", "new", post_predict = pp)
#' x$stats
#' @export
croc <- function(data, signal_factor, noise, signal_response = NULL, post_predict = NULL,
                 by = NULL, bounds = NULL, n_ratings = NULL, dont_know = NULL,
                 probs = c(.025, .975)) {
  roc_run("croc", data, signal_factor, noise, signal_response, post_predict, by, bounds,
          n_ratings, dont_know, probs, need_rt = FALSE, fun = croc_one)
}

#' @rdname croc
#' @param n_quantiles Number of RT quantiles per response, giving
#'   `n_quantiles - 1` bins (28, as in Reynolds et al.).
#' @examples
#' plot_lroc(lroc(dat, "S", "new", post_predict = pp))
#' @export
lroc <- function(data, signal_factor, noise, signal_response = NULL, post_predict = NULL,
                 by = NULL, bounds = NULL, n_quantiles = 28, n_ratings = NULL,
                 probs = c(.025, .975)) {
  roc_run("lroc", data, signal_factor, noise, signal_response, post_predict, by, bounds,
          n_ratings, NULL, probs, need_rt = TRUE,
          fun = function(tw, s) lroc_one(tw, s, n_quantiles), n_quantiles = n_quantiles)
}

#' @rdname croc
#' @param type `"confidence"` (C-RROC) or `"latency"` (L-RROC).
#' @param conf For `type = "latency"`, optional confidence level(s): the
#'   L-RROCs within each level.
#' @examples
#' rroc(dat, "S", "new", type = "latency", conf = 3)$stats
#' @export
rroc <- function(data, signal_factor, noise, signal_response = NULL, post_predict = NULL,
                 by = NULL, bounds = NULL, type = c("confidence", "latency"), conf = NULL,
                 n_quantiles = 28, n_ratings = NULL, probs = c(.025, .975)) {
  type <- match.arg(type)
  if (!is.null(conf) && type == "confidence") stop("conf applies to type = \"latency\" only")
  if (!is.null(conf)) conf <- as.integer(conf)
  roc_run("rroc", data, signal_factor, noise, signal_response, post_predict, by, bounds,
          n_ratings, NULL, probs, need_rt = type == "latency",
          fun = function(tw, s) {
            if (!is.null(conf) && any(conf < 1 | conf > s$K)) stop("conf must be in 1 .. ", s$K)
            rroc_one(tw, s, type, conf, n_quantiles)
          },
          rroc_type = type, conf = conf, n_quantiles = n_quantiles)
}

#' @export
print.emc_roc <- function(x, ...) {
  lab <- c(croc = "C-ROC", lroc = "L-ROC", rroc = "Response-contingent ROC")[x$type]
  cat(lab, if (x$type == "rroc") paste0(" (", x$rroc_type, ")"), ": signal response ",
      x$signal_response, ", ", length(x$signal_levels), " signal level(s)",
      if (!is.null(x$by)) paste0(", by ", paste(x$by, collapse = " x ")),
      if (!is.null(x$pred)) paste0(", ", length(unique(x$pred$stats$postn)), " predictive draws"),
      "\n", sep = "")
  print(x$stats, ...)
  invisible(x)
}

# ---- plotting ------------------------------------------------------------------

# Pointwise median curve and interval over draws of a fixed-length curve.
roc_band <- function(curves, probs = c(.025, .975)) {
  FAR <- sapply(curves, `[[`, "FAR"); HR <- sapply(curves, `[[`, "HR")
  if (is.null(dim(FAR))) { FAR <- matrix(FAR, ncol = length(curves)); HR <- matrix(HR, ncol = length(curves)) }
  q <- function(m, p) apply(m, 1, stats::quantile, p, na.rm = TRUE, names = FALSE)
  data.frame(FAR = apply(FAR, 1, stats::median, na.rm = TRUE), HR = apply(HR, 1, stats::median, na.rm = TRUE),
             FAR_lo = q(FAR, probs[1]), FAR_hi = q(FAR, probs[2]),
             HR_lo = q(HR, probs[1]), HR_hi = q(HR, probs[2]))
}

# Every draw thin, the pointwise median thick, crosshairs at the points.
roc_draw_pred <- function(curves, col, probs, crosshairs = TRUE, lty = 1, alpha = .08) {
  if (!length(curves)) return(invisible(NULL))
  for (cv in curves) graphics::lines(cv$FAR, cv$HR, col = grDevices::adjustcolor(col, alpha), lty = lty)
  b <- roc_band(curves, probs)
  graphics::lines(b$FAR, b$HR, col = col, lwd = 2, lty = lty)
  if (crosshairs) {
    graphics::segments(b$FAR_lo, b$HR, b$FAR_hi, b$HR, col = col)
    graphics::segments(b$FAR, b$HR_lo, b$FAR, b$HR_hi, col = col)
  }
  invisible(b)
}

roc_frame <- function(main, xlab = "False-alarm rate", ylab = "Hit rate", lim = c(0, 1), ...) {
  graphics::plot(NA, xlim = lim, ylim = lim, xlab = xlab, ylab = ylab, main = main, ...)
  graphics::abline(0, 1, lty = 3, col = "grey50")
}

# Panel keys (by levels in order) and their labels; sets the panel layout.
roc_panels <- function(x, d, mfrow, n_extra = 1) {
  key <- if (is.null(x$by)) rep("all", nrow(d)) else
    do.call(paste, c(unname(as.list(d[, x$by, drop = FALSE])), sep = " "))
  cells <- unique(key)
  n <- length(cells) * n_extra
  op <- NULL
  if (n > 1 || !is.null(mfrow)) {
    if (is.null(mfrow)) mfrow <- if (n <= 4) c(1, n) else rev(grDevices::n2mfrow(n))
    op <- graphics::par(mfrow = mfrow)
  }
  list(key = key, cells = cells, op = op)
}

roc_key <- function(x, d) {
  if (is.null(d)) return(NULL)
  if (is.null(x$by)) rep("all", nrow(d)) else do.call(paste, c(unname(as.list(d[, x$by, drop = FALSE])), sep = " "))
}

roc_title <- function(main, i, default, cell) {
  if (!is.null(main)) return(rep_len(main, i)[i])
  if (cell == "all") default else paste0(default, ": ", cell)
}

roc_cols <- function(n, col_data, col_pred) {
  if (is.null(col_data)) col_data <- if (n == 1) "black" else c("#0072B2", "#D55E00", "#009E73", "#CC79A7")
  col_data <- rep_len(col_data, n)
  if (is.null(col_pred)) col_pred <- if (n == 1) "#0072B2" else col_data
  list(data = col_data, pred = rep_len(col_pred, n))
}

roc_draws <- function(p, xcol = "FAR", ycol = "HR") {
  if (is.null(p) || !nrow(p)) return(list())
  lapply(split(p, p$postn), function(d) data.frame(FAR = d[[xcol]], HR = d[[ycol]]))
}

#' Plot Confidence, Latency and Response-Contingent ROCs
#'
#' Plots of the results of [croc()], [lroc()] and [rroc()]: the data as a
#' line (with points at the discrete criteria), posterior predictives as every
#' draw's curve (thin), the pointwise median curve (thick) and, at the points
#' of the C-ROC, zROC and C-RROC, crosshairs spanning the pointwise
#' predictive interval. The random ROC through the classification point is
#' dashed (the classification point is the large dot). One panel per `by`
#' level; several signal levels share a panel in different colours.
#' `plot_roc_auc()` shows statistics as bars (data) with points and intervals
#' (predictive median and interval), for one analysis or, with `by` (e.g.
#' subjects), grouped by level.
#'
#' @param x For `plot_croc()`, `plot_lroc()` and `plot_rroc()` the result of
#'   [croc()], [lroc()] or [rroc()] (for `plot_rroc()` also a list of [rroc()]
#'   results, from which `type` and `conf` select); for `plot_roc_auc()` one
#'   such result or a list of them.
#' @param z For `plot_croc()`: the zROC instead of the C-ROC.
#' @param main Panel title(s) (default: the curve's name and the `by` level).
#' @param col_data,col_pred Colours of the data and the predictions, one per
#'   signal level (`plot_rroc()`: `col_pred` per signal and noise response).
#' @param random Whether to draw the random ROC.
#' @param legend Whether to draw a legend.
#' @param mfrow Panel layout for several panels (default: one row up to four
#'   panels).
#' @param ... Further arguments passed to `plot()` (`barplot()` for
#'   `plot_roc_auc()`).
#' @return Invisibly, `x`.
#' @examples
#' sim <- function(n) {
#'   S <- factor(rep(c("new", "old"), each = n / 2))
#'   e <- rnorm(n, ifelse(S == "old", 1, 0), ifelse(S == "old", 1.25, 1))
#'   R <- factor(ifelse(e > .5, "old", "new"), levels = c("new", "old"))
#'   RR <- as.integer(cut(abs(e - .5), c(0, .4, 1, Inf)))
#'   data.frame(S = S, R = R, RR = RR, rt = .3 + rexp(n, 1 + 2 * abs(e - .5)))
#' }
#' set.seed(1)
#' dat <- sim(1000)
#' pp <- do.call(rbind, lapply(1:20, function(i) cbind(postn = i, sim(1000))))
#' x <- croc(dat, "S", "new", post_predict = pp)
#' oldpar <- par(mfrow = c(1, 2))
#' plot_croc(x)
#' plot_croc(x, z = TRUE)
#' par(oldpar)
#' @export
plot_croc <- function(x, z = FALSE, main = NULL, col_data = NULL, col_pred = NULL,
                      random = TRUE, legend = TRUE, mfrow = NULL, ...) {
  if (!inherits(x, "emc_roc") || x$type != "croc") stop("x must be the result of croc()")
  d <- if (z) x$zroc else x$curve
  p <- if (z) x$pred$zroc else x$pred$curve
  xc <- if (z) "zF" else "FAR"; yc <- if (z) "zH" else "HR"
  pan <- roc_panels(x, d, mfrow)
  if (!is.null(pan$op)) on.exit(graphics::par(pan$op))
  pkey <- roc_key(x, p)
  cols <- roc_cols(length(x$signal_levels), col_data, col_pred)
  lim <- c(0, 1)
  if (z) {
    v <- c(d$zF, d$zH, p$zF, p$zH); v <- v[is.finite(v)]
    lim <- if (length(v)) range(v) else c(-2, 2)
    lim <- lim + c(-.05, .05) * max(diff(lim), .1)
  }
  for (i in seq_along(pan$cells)) {
    ce <- pan$cells[i]
    roc_frame(roc_title(main, i, if (z) "zROC" else "C-ROC", ce),
              if (z) "z(false-alarm rate)" else "False-alarm rate",
              if (z) "z(hit rate)" else "Hit rate", lim, ...)
    for (j in seq_along(x$signal_levels)) {
      sl <- x$signal_levels[j]
      dd <- d[pan$key == ce & d$signal == sl, ]
      if (!is.null(p)) roc_draw_pred(roc_draws(p[pkey == ce & p$signal == sl, ], xc, yc), cols$pred[j], x$probs)
      graphics::lines(dd[[xc]], dd[[yc]], col = cols$data[j], lwd = 2)
      if (z) {
        graphics::points(dd$zF, dd$zH, pch = 16, col = cols$data[j])
      } else {
        cls <- !x$dont_know & dd$k == x$K
        graphics::points(dd$FAR, dd$HR, pch = ifelse(cls, 16, 1), cex = ifelse(cls, 1.6, 1), col = cols$data[j])
        if (random && !x$dont_know)
          graphics::lines(c(0, dd$FAR[x$K + 1], 1), c(0, dd$HR[x$K + 1], 1), lty = 2, lwd = 2, col = cols$data[j])
      }
    }
    if (legend && i == 1) roc_legend(x, cols, if (z) "topleft" else "bottomright",
                                     if (z) "zROC" else "C-ROC", random && !z && !x$dont_know, TRUE)
  }
  invisible(x)
}

roc_legend <- function(x, cols, pos, what, random, crosshairs) {
  model <- if (crosshairs) "model (draws, median, 95%)" else "model (draws, median)"
  if (length(x$signal_levels) == 1) {
    lab <- c(paste(what, "data"), if (random) "random ROC", if (!is.null(x$pred)) model)
    graphics::legend(pos, lab, col = c(cols$data[1], if (random) cols$data[1], if (!is.null(x$pred)) cols$pred[1]),
                     lty = c(1, if (random) 2, if (!is.null(x$pred)) 1), lwd = 2, bty = "n")
  } else {
    graphics::legend(pos, x$signal_levels, col = cols$data, lwd = 2, bty = "n")
  }
}

#' @rdname plot_croc
#' @examples
#' plot_lroc(lroc(dat, "S", "new", post_predict = pp))
#' @export
plot_lroc <- function(x, main = NULL, col_data = NULL, col_pred = NULL, random = TRUE,
                      legend = TRUE, mfrow = NULL, ...) {
  if (!inherits(x, "emc_roc") || x$type != "lroc") stop("x must be the result of lroc()")
  d <- x$curve; p <- x$pred$curve; nq <- x$n_quantiles
  pan <- roc_panels(x, d, mfrow)
  if (!is.null(pan$op)) on.exit(graphics::par(pan$op))
  pkey <- roc_key(x, p)
  cols <- roc_cols(length(x$signal_levels), col_data, col_pred)
  for (i in seq_along(pan$cells)) {
    ce <- pan$cells[i]
    roc_frame(roc_title(main, i, "L-ROC", ce), ...)
    for (j in seq_along(x$signal_levels)) {
      sl <- x$signal_levels[j]
      dd <- d[pan$key == ce & d$signal == sl, ]
      if (!is.null(p)) roc_draw_pred(roc_draws(p[pkey == ce & p$signal == sl, ]), cols$pred[j], x$probs,
                                     crosshairs = FALSE)
      graphics::lines(dd$FAR, dd$HR, col = cols$data[j], lwd = 2)
      graphics::points(dd$FAR[nq], dd$HR[nq], pch = 16, cex = 1.6, col = cols$data[j])
      if (random) graphics::lines(c(0, dd$FAR[nq], 1), c(0, dd$HR[nq], 1), lty = 2, lwd = 2, col = cols$data[j])
    }
    if (legend && i == 1) roc_legend(x, cols, "bottomright", "L-ROC", random, FALSE)
  }
  invisible(x)
}

#' @rdname plot_croc
#' @param type,conf For `plot_rroc()`: which response-contingent ROC to plot
#'   (default: that of `x`); `conf` selects confidence level(s) of a latency
#'   result made with `conf`, one panel each.
#' @param labels For `plot_rroc()`: legend labels of the signal and noise
#'   responses (default: their `R` levels); for `plot_roc_auc()`: labels of
#'   `stats`.
#' @examples
#' plot_rroc(rroc(dat, "S", "new", post_predict = pp))
#' @export
plot_rroc <- function(x, type = NULL, conf = NULL, main = NULL, col_data = "black",
                      col_pred = c("#0072B2", "#D55E00"), labels = NULL, legend = TRUE,
                      mfrow = NULL, ...) {
  xs <- if (inherits(x, "emc_roc")) list(x) else x
  if (!all(vapply(xs, function(r) inherits(r, "emc_roc") && r$type == "rroc", TRUE)))
    stop("x must be the result of rroc() or a list of them")
  if (!is.null(type)) {
    type <- match.arg(type, c("confidence", "latency"))
    xs <- xs[vapply(xs, function(r) r$rroc_type == type, TRUE)]
  }
  if (!is.null(conf)) xs <- xs[vapply(xs, function(r) all(conf %in% r$conf), TRUE)]
  else if (length(xs) > 1) xs <- xs[vapply(xs, function(r) is.null(r$conf), TRUE)]
  if (!length(xs)) stop("no rroc() result in x has this type and conf")
  x <- xs[[1]]
  if (is.null(conf)) conf <- x$conf
  confs <- if (is.null(conf)) list(NULL) else as.list(conf)
  d <- x$curve; p <- x$pred$curve
  pan <- roc_panels(x, d, mfrow, length(x$signal_levels) * length(confs))
  if (!is.null(pan$op)) on.exit(graphics::par(pan$op))
  pkey <- roc_key(x, p)
  if (is.null(labels)) labels <- paste(c(x$signal_response, x$noise_response), "responses")
  col_pred <- rep_len(col_pred, 2)
  is_conf <- x$rroc_type == "confidence"
  i <- 0
  for (ce in pan$cells) for (sl in x$signal_levels) for (cc in confs) {
    i <- i + 1
    sel <- function(t, key) {
      w <- key == ce & t$signal == sl
      if (!is.null(cc)) w <- w & t$conf == cc
      t[w, ]
    }
    dd <- sel(d, pan$key)
    default <- paste0(if (is_conf) "C-RROC" else "L-RROC", if (!is.null(cc)) paste0(", confidence ", cc),
                      if (length(x$signal_levels) > 1) paste0(", ", sl))
    roc_frame(roc_title(main, i, default, ce), ...)
    if (!is.null(p)) {
      pp <- sel(p, pkey)
      roc_draw_pred(roc_draws(pp[pp$response == "signal", ]), col_pred[1], x$probs, crosshairs = is_conf)
      roc_draw_pred(roc_draws(pp[pp$response == "noise", ]), col_pred[2], x$probs, crosshairs = is_conf, lty = 2)
    }
    ds <- dd[dd$response == "signal", ]; dn <- dd[dd$response == "noise", ]
    graphics::lines(ds$FAR, ds$HR, col = col_data, lwd = 2)
    graphics::lines(dn$FAR, dn$HR, col = col_data, lwd = 2, lty = 2)
    if (is_conf) {
      graphics::points(ds$FAR, ds$HR, pch = 16, col = col_data)
      graphics::points(dn$FAR, dn$HR, pch = 1, col = col_data)
    }
    if (legend && i == 1)
      graphics::legend("bottomright", c(paste(labels, "data"), if (!is.null(p)) paste(labels, "model")),
                       col = c(col_data, col_data, if (!is.null(p)) col_pred),
                       lty = c(1, 2, if (!is.null(p)) c(1, 2)), lwd = 2, bty = "n")
  }
  invisible(x)
}

#' @rdname plot_croc
#' @param stats For `plot_roc_auc()`: names of the statistics to show (default:
#'   every AUC in `x`), e.g. `c("rAUC", "lAUC", "cAUC")`.
#' @param ylim For `plot_roc_auc()`: y-axis limits (default: .4 to 1, widened
#'   to cover the values).
#' @param h For `plot_roc_auc()`: height of a dotted reference line (chance).
#' @examples
#' plot_roc_auc(list(croc(dat, "S", "new", post_predict = pp), lroc(dat, "S", "new", post_predict = pp)),
#'              stats = c("rAUC", "lAUC", "cAUC"), labels = c("random", "latency", "confidence"))
#' @export
plot_roc_auc <- function(x, stats = NULL, labels = NULL, main = "AUC", ylim = NULL, h = .5,
                         col_pred = "#0072B2", ...) {
  xs <- if (inherits(x, "emc_roc")) list(x) else x
  if (!length(xs) || !all(vapply(xs, inherits, TRUE, "emc_roc")))
    stop("x must be the result of croc(), lroc() or rroc(), or a list of them")
  by <- xs[[1]]$by
  if (!all(vapply(xs, function(r) identical(r$by, by), TRUE))) stop("the results in x must have the same by")
  has_pred <- all(vapply(xs, function(r) !is.null(r$pred), TRUE))
  cols <- c(by, "signal", "stat", "data", if (has_pred) c("pred_median", "pred_lo", "pred_hi"))
  tab <- do.call(rbind, lapply(xs, function(r) r$stats[, cols, drop = FALSE]))
  if (is.null(stats)) stats <- unique(tab$stat[grepl("AUC", tab$stat)])
  if (!all(stats %in% tab$stat)) stop("not in x: ", paste(setdiff(stats, tab$stat), collapse = ", "))
  if (is.null(labels)) labels <- stats
  gcols <- c(by, if (length(unique(tab$signal)) > 1) "signal")
  tab$group <- if (length(gcols)) do.call(paste, c(unname(as.list(tab[, gcols, drop = FALSE])), sep = " ")) else "all"
  groups <- unique(tab$group)
  get <- function(what) matrix(sapply(groups, function(g) tab[[what]][tab$group == g][match(stats, tab$stat[tab$group == g])]),
                               nrow = length(stats), dimnames = list(labels, groups))
  dat <- get("data")
  vals <- c(dat, if (has_pred) c(get("pred_lo"), get("pred_hi")))
  if (is.null(ylim)) {
    r <- range(c(.4, 1, vals), na.rm = TRUE)
    ylim <- r + c(-.02, .02) * diff(r) * (r != c(.4, 1))
  }
  one <- length(groups) == 1
  mp <- if (one) {
    graphics::barplot(dat[, 1], names.arg = labels, ylim = ylim, xpd = FALSE, col = "grey85", main = main,
                      ylab = "value", ...)
  } else {
    graphics::barplot(dat, beside = TRUE, ylim = ylim, xpd = FALSE,
                      col = grDevices::gray.colors(length(stats), .95, .6), main = main,
                      xlab = paste(gcols, collapse = " "),
                      ylab = "value", legend.text = if (length(stats) > 1) labels,
                      args.legend = list(x = "top", bty = "n", horiz = TRUE, cex = .8), ...)
  }
  graphics::box()
  if (!is.null(h)) graphics::abline(h = h, lty = 3)
  if (has_pred) {
    mp <- matrix(mp, nrow = length(stats))
    graphics::points(mp, get("pred_median"), pch = 16, col = col_pred)
    ok <- is.finite(get("pred_lo")) & is.finite(get("pred_hi")) & get("pred_hi") > get("pred_lo")
    graphics::arrows(mp[ok], get("pred_lo")[ok], mp[ok], get("pred_hi")[ok], angle = 90, code = 3,
                     length = .03, col = col_pred)
  }
  invisible(x)
}
