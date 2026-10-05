# Plotting helpers for choice + rating data (rating layer; nothing here is
# model specific). The response is folded into one ordered factor running from
# the highest rating of the first response level to the highest rating of the
# second (rating_fold), e.g. NEW.3, NEW.2, NEW.1, OLD.1, OLD.2, OLD.3.

# Data frame from a data frame or an emc object; RR and R are required.
rating_input <- function(x) {
  if (inherits(x, "emc")) x <- get_data(x)
  if (!is.data.frame(x) || !all(c("R", "rt", "RR") %in% names(x)))
    stop("need a data frame (or emc object) with columns R, rt and RR")
  if (!is.factor(x$R)) x$R <- factor(x$R)
  x
}

rating_cells <- function(x, factors) {
  if (is.null(factors) || length(factors) == 0) return(factor(rep("all", nrow(x))))
  interaction(x[, factors, drop = FALSE], drop = TRUE, sep = " ", lex.order = TRUE)
}

#' Summarise Choice and Rating Data
#'
#' Response proportions and response-time quantiles for each folded response
#' (choice by rating), per design cell. Statistics are computed per subject and
#' then averaged over subjects. For posterior predictives (a `postn` column, as
#' returned by `predict()`), they are computed per posterior draw.
#'
#' @param input A data frame with columns `subjects`, `R`, `rt` and `RR`
#'   (`RR` = 1 the lowest rating), or an emc object.
#' @param factors Character vector of design factors defining the cells.
#' @param probs RT quantile probabilities.
#' @param n_ratings Number of rating categories K (default: the largest RR).
#' @return A data frame with columns for the cell, `postn` (if present), the
#'   folded response `resp`, proportion `p` and one column per quantile.
#' @examples
#' dat <- data.frame(subjects = factor(1), R = factor(c("a", "b", "a")),
#'                   rt = c(.5, .7, .6), RR = c(3, 1, 2))
#' rating_summary(dat)
#' @export
rating_summary <- function(input, factors = NULL, probs = c(.1, .5, .9), n_ratings = NULL) {
  x <- rating_input(input)
  K <- if (is.null(n_ratings)) max(x$RR, na.rm = TRUE) else n_ratings
  if (is.null(x$subjects)) x$subjects <- factor(1)
  if (is.null(x$postn)) x$postn <- 1
  x <- x[!is.na(x$R) & !is.na(x$RR), ]
  x$resp <- rating_fold(x$R, x$RR, K)
  x$cell <- rating_cells(x, factors)
  qn <- paste0("q", round(100 * probs))
  one <- function(y) {
    n <- tabulate(as.integer(y$resp), nlevels(y$resp))
    q <- t(sapply(split(y$rt, y$resp), function(r)
      if (length(r)) stats::quantile(r, probs, names = FALSE) else rep(NA, length(probs))))
    if (length(probs) == 1) q <- t(q)
    colnames(q) <- qn
    data.frame(resp = factor(levels(y$resp), levels = levels(y$resp)), p = n / sum(n), q,
               row.names = NULL)
  }
  keys <- interaction(x$postn, x$cell, x$subjects, drop = TRUE, sep = "\r")
  per <- lapply(split(x, keys), function(y)
    cbind(postn = y$postn[1], cell = y$cell[1], subjects = y$subjects[1], one(y)))
  per <- do.call(rbind, per)
  agg <- stats::aggregate(per[, c("p", qn)], by = per[, c("postn", "cell", "resp")],
                          FUN = mean, na.rm = TRUE, na.action = stats::na.pass)
  agg[order(agg$postn, agg$cell, agg$resp), ]
}

# Draw data (points) and predictions (line + band) for one statistic across the
# folded responses.
rating_draw <- function(xs, dat, pp, col_data, col_pred, band = c(.025, .975)) {
  if (!is.null(pp)) {
    lo <- apply(pp, 2, stats::quantile, band[1], na.rm = TRUE)
    hi <- apply(pp, 2, stats::quantile, band[2], na.rm = TRUE)
    ok <- is.finite(lo) & is.finite(hi)
    graphics::polygon(c(xs[ok], rev(xs[ok])), c(lo[ok], rev(hi[ok])),
                      col = grDevices::adjustcolor(col_pred, .25), border = NA)
    graphics::lines(xs, colMeans(pp, na.rm = TRUE), col = col_pred, lwd = 2)
  }
  graphics::points(xs, dat, pch = 16, col = col_data)
}

#' Plot Choice and Rating Data by Folded Response
#'
#' For each design cell, response proportions (top row) and response-time
#' quantiles (bottom row) against the folded response (choice by rating), with
#' posterior predictives as lines and 95% bands.
#'
#' @inheritParams rating_summary
#' @param post_predict Optional posterior predictives (from `predict()`).
#' @param col_data,col_pred Colours for data and predictions.
#' @param main Optional overall title.
#' @param ... Further arguments passed to `plot()`.
#' @return Invisibly, a list with the data and prediction summaries.
#' @examples
#' dat <- data.frame(subjects = factor(1), R = factor(c("a", "b", "a")),
#'                   rt = c(.5, .7, .6), RR = c(3, 1, 2))
#' plot_ratings(dat)
#' @export
plot_ratings <- function(input, post_predict = NULL, factors = NULL,
                         probs = c(.1, .5, .9), n_ratings = NULL,
                         col_data = "black", col_pred = "#0072B2", main = NULL, ...) {
  x <- rating_input(input)
  K <- if (is.null(n_ratings)) max(x$RR, na.rm = TRUE) else n_ratings
  sd <- rating_summary(x, factors, probs, K)
  sp <- if (!is.null(post_predict)) rating_summary(rating_input(post_predict), factors, probs, K)
  qn <- paste0("q", round(100 * probs))
  cells <- levels(factor(sd$cell))
  labs <- levels(sd$resp)
  xs <- seq_along(labs)
  op <- graphics::par(mfcol = c(2, length(cells)), mar = c(4, 4, 2, .5), oma = c(0, 0, if (is.null(main)) 0 else 2, 0))
  on.exit(graphics::par(op))
  qlim <- range(c(unlist(sd[, qn]), if (!is.null(sp)) unlist(sp[, qn])), na.rm = TRUE)
  plim <- c(0, max(c(sd$p, if (!is.null(sp)) sp$p), na.rm = TRUE))
  wide <- function(s, v) {
    if (is.null(s)) return(NULL)
    w <- stats::reshape(s[, c("postn", "resp", v)], idvar = "postn", timevar = "resp", direction = "wide")
    as.matrix(w[, paste(v, labs, sep = "."), drop = FALSE])
  }
  for (ce in cells) {
    d1 <- sd[sd$cell == ce, ]
    p1 <- if (!is.null(sp)) sp[sp$cell == ce, ]
    graphics::plot(xs, d1$p, type = "n", ylim = plim, xaxt = "n", xlab = "", ylab = "Proportion",
                   main = ce, ...)
    graphics::axis(1, xs, labs, cex.axis = .8)
    rating_draw(xs, d1$p, wide(p1, "p"), col_data, col_pred)
    graphics::plot(xs, d1[[qn[1]]], type = "n", ylim = qlim, xaxt = "n", xlab = "Response",
                   ylab = "RT quantiles (s)", ...)
    graphics::axis(1, xs, labs, cex.axis = .8)
    for (q in qn) rating_draw(xs, d1[[q]], wide(p1, q), col_data, col_pred)
  }
  if (!is.null(main)) graphics::mtext(main, outer = TRUE, font = 2)
  invisible(list(data = sd, predicted = sp))
}

#' z-transformed Receiver Operating Characteristics for Rating Data
#'
#' Cumulative response probabilities from the most confident response of the
#' second response level (e.g. "old") down to the most confident of the first,
#' z-transformed, for each signal cell against the noise cell. Probabilities are
#' pooled over subjects.
#'
#' @inheritParams rating_summary
#' @param signal_factor Factor distinguishing noise from signal trials.
#' @param noise Level of `signal_factor` that is noise (all others are signal).
#' @param by Further factors: one ROC per cell.
#' @return `zroc()`: a data frame with columns for the `by` cell, the signal
#'   level, the criterion `k` and `zF`, `zH` (NA where a proportion is 0 or 1).
#' @examples
#' dat <- data.frame(S = factor(rep(c("n", "s"), each = 6)),
#'                   R = factor(rep(c("a", "a", "a", "b", "b", "b"), 2)),
#'                   rt = 1, RR = rep(c(3, 2, 1, 1, 2, 3), 2))
#' zroc(dat, "S", "n")
#' @export
zroc <- function(input, signal_factor, noise, by = NULL, n_ratings = NULL) {
  x <- rating_input(input)
  K <- if (is.null(n_ratings)) max(x$RR, na.rm = TRUE) else n_ratings
  x <- x[!is.na(x$R) & !is.na(x$RR), ]
  x$resp <- rating_fold(x$R, x$RR, K)
  cum <- function(y) {
    p <- tabulate(as.integer(y$resp), nlevels(x$resp)) / nrow(y)
    rev(cumsum(rev(p)))[-1]  # P(resp >= k), k = 2 .. 2K
  }
  x$cell <- rating_cells(x, by)
  lev <- setdiff(levels(factor(x[[signal_factor]])), noise)
  out <- list()
  for (ce in levels(x$cell)) {
    xc <- x[x$cell == ce, ]
    Fp <- cum(xc[xc[[signal_factor]] == noise, ])
    for (s in lev) {
      Hp <- cum(xc[xc[[signal_factor]] == s, ])
      z <- function(p) ifelse(p > 0 & p < 1, stats::qnorm(p), NA)
      out[[length(out) + 1]] <- data.frame(cell = ce, signal = s, k = seq_along(Fp),
                                           F = Fp, H = Hp, zF = z(Fp), zH = z(Hp))
    }
  }
  do.call(rbind, out)
}

#' @rdname zroc
#' @param post_predict Optional posterior predictives (from `predict()`): the
#'   zROC of each posterior draw is drawn as a thin line.
#' @param col Colours for the signal levels.
#' @param ... Further arguments passed to `plot()`.
#' @export
plot_zroc <- function(input, signal_factor, noise, by = NULL, post_predict = NULL,
                      n_ratings = NULL, col = NULL, ...) {
  x <- rating_input(input)
  zd <- zroc(x, signal_factor, noise, by, n_ratings)
  zp <- NULL
  if (!is.null(post_predict)) {
    pp <- rating_input(post_predict)
    zp <- lapply(split(pp, pp$postn), zroc, signal_factor = signal_factor,
                 noise = noise, by = by, n_ratings = n_ratings)
  }
  sig <- unique(zd$signal)
  if (is.null(col)) col <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")[seq_along(sig)]
  cells <- unique(zd$cell)
  op <- graphics::par(mfrow = c(1, length(cells)), mar = c(4, 4, 2, .5))
  on.exit(graphics::par(op))
  lim <- range(c(zd$zF, zd$zH), na.rm = TRUE)
  for (ce in cells) {
    graphics::plot(NA, xlim = lim, ylim = lim, xlab = "z(false alarm)", ylab = "z(hit)",
                   main = if (ce == "all") "" else ce, ...)
    graphics::abline(0, 1, lty = 3, col = "grey")
    for (i in seq_along(sig)) {
      if (!is.null(zp)) for (z in zp) {
        zz <- z[z$cell == ce & z$signal == sig[i], ]
        graphics::lines(zz$zF, zz$zH, col = grDevices::adjustcolor(col[i], .15))
      }
      zz <- zd[zd$cell == ce & zd$signal == sig[i], ]
      graphics::points(zz$zF, zz$zH, pch = 16, col = col[i], type = "b")
    }
    graphics::legend("topleft", legend = sig, col = col, pch = 16, bty = "n")
  }
  invisible(list(data = zd, predicted = zp))
}
