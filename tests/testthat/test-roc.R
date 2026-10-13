# croc(), lroc(), rroc() and their plots. Reference values (golden/roc/roc_functions_ref.rds) are
# from the Stage M8g scripts (rating-work/mtlba/stageM8g/R/roc_functions.R, validated against the
# MTR paper's stored values), made by rating-work/mtlba/stageM8g/m8g2/make_ref.R.

roc_sim_design <- function() {
  design(factors = list(subjects = 1:2, S = c("new", "old")), Rlevels = c("new", "old"),
         model = function() MTLBA(n_ratings = 3, rule = "BE", drifts = "positive",
                                  lowest_threshold = "above_A"),
         matchfun = function(d) d$S == d$lR,
         formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, c1 ~ lR, c2 ~ lR), report_p_vector = FALSE)
}
roc_sim_p <- function(shift = 0) {
  p <- c(v = 1, v_lMTRUE = 1.2, B = log(1), A = log(.5), t0 = log(.25),
         c1 = log(.6), c1_lRold = .15, c2 = log(.6), c2_lRold = -.1)
  p <- rbind(p, p + c(0, .4, .1, 0, 0, .1, 0, -.1, 0))
  p[, "v_lMTRUE"] <- p[, "v_lMTRUE"] + shift
  rownames(p) <- 1:2
  p
}
# 2 subjects x 300 trials and 5 "posterior predictive" draws (other parameters)
roc_sim_data <- function() {
  des <- suppressMessages(roc_sim_design())
  # R's default generators whatever an earlier test file set (several set L'Ecuyer-CMRG)
  kind <- RNGkind()
  on.exit(RNGkind(kind[1], kind[2], kind[3]))
  set.seed(2026, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  dat <- make_data(roc_sim_p(), des, n_trials = 150)
  pp <- do.call(rbind, lapply(1:5, function(i)
    cbind(postn = i, make_data(roc_sim_p((i - 3) / 10), des, n_trials = 150))))
  list(data = dat, pp = pp)
}
roc_ref <- function() readRDS(test_path("golden", "roc", "roc_functions_ref.rds"))
roc_fp <- function(x) c(n = nrow(x), rt = sum(x$rt), RR = sum(x$RR * as.integer(x$R)))
roc_maxdiff <- function(a, b) {
  a <- as.numeric(unlist(a)); b <- as.numeric(unlist(b))
  stopifnot(length(a) == length(b), identical(is.na(a), is.na(b)))
  if (all(is.na(a))) 0 else max(abs(a - b), na.rm = TRUE)
}
fin <- function(v) ifelse(is.finite(v), v, NA)

# Every statistic and curve of one analysis (no by, one signal level) against the reference
roc_check_against <- function(dat, ref, pp = NULL) {
  x <- croc(dat, "S", "new", post_predict = pp)
  l <- lroc(dat, "S", "new", post_predict = pp)
  rc <- rroc(dat, "S", "new", post_predict = pp)
  rl <- rroc(dat, "S", "new", post_predict = pp, type = "latency")
  rk <- rroc(dat, "S", "new", post_predict = pp, type = "latency", conf = 1:3)
  st <- rbind(x$stats, l$stats, rc$stats, rl$stats, rk$stats)
  r <- ref$summary[match(st$stat, ref$summary$stat), ]
  expect_equal(st$stat, r$stat)
  expect_lt(roc_maxdiff(st$data, r$data), 1e-12)
  expect_true(all(st$n == ref$summary$data[ref$summary$stat == "n"]))
  cur <- function(t, resp = NULL, cc = NULL) {
    if (!is.null(resp)) t <- t[t$response == resp, ]
    if (!is.null(cc)) t <- t[t$conf == cc, ]
    t[, c("k", "FAR", "HR")]
  }
  expect_lt(roc_maxdiff(cur(x$curve), ref$croc), 1e-12)
  expect_lt(roc_maxdiff(x$zroc[, c("k", "zF", "zH")], data.frame(ref$zroc$k, fin(ref$zroc$zF), fin(ref$zroc$zH))), 1e-12)
  expect_lt(roc_maxdiff(cur(l$curve), ref$lroc), 1e-12)
  for (resp in c("signal", "noise")) {
    rr <- if (resp == "signal") "sig" else "noise"
    expect_lt(roc_maxdiff(cur(rc$curve, resp), ref$crroc[[rr]]), 1e-12)
    expect_lt(roc_maxdiff(cur(rl$curve, resp), ref$lrroc[[rr]]), 1e-12)
    for (cc in 1:3) expect_lt(roc_maxdiff(cur(rk$curve, resp, cc), ref$lrroc_conf[[cc]][[rr]]), 1e-12)
  }
  invisible(st)
}

test_that("croc, lroc, rroc equal the M8g reference on simulated MTLBA data", {
  s <- roc_sim_data()
  ref <- roc_ref()
  # the simulation must be the one the reference was made from
  expect_equal(roc_fp(s$data), ref$sim_fingerprint$data)
  expect_equal(roc_fp(s$pp), ref$sim_fingerprint$pp)
  roc_check_against(s$data, ref$sim)
})

test_that("croc, lroc, rroc equal the M8g reference on r94 subject 1 (all raw trials)", {
  f <- "~/Documents/Papers/SUBMITTED/MTR/EMC2/r94/r94_emc.rds"
  skip_if_not(file.exists(f), "r94 data not available")
  r94 <- readRDS(f)
  d <- r94$data[r94$data$subjects == 1, ]
  d$S <- factor(d$NWS != "new", labels = c("new", "old"))
  roc_check_against(d, roc_ref()$r94_s1)
})

test_that("data statistics of all 11 r94 subjects equal the paper's stored values", {
  f <- "~/Documents/Papers/SUBMITTED/MTR/EMC2/r94/r94_emc.rds"
  fa <- "~/Documents/Papers/SUBMITTED/MTR/FirstSubmitPrep/Fits94/AUC.RData"
  skip_if_not(file.exists(f) && file.exists(fa), "r94 data or the paper's AUC.RData not available")
  d <- readRDS(f)$data
  d$S <- factor(d$NWS != "new", labels = c("new", "old"))
  e <- new.env(); load(fa, e)
  st <- rbind(croc(d, "S", "new", by = "subjects")$stats, lroc(d, "S", "new", by = "subjects")$stats,
              rroc(d, "S", "new", by = "subjects")$stats,
              rroc(d, "S", "new", by = "subjects", type = "latency")$stats,
              rroc(d, "S", "new", by = "subjects", type = "latency", conf = 1:3)$stats)
  # the mapping of rating-work/mtlba/stageM8g/R/run_roc.R (results/dmc_paper_check.csv)
  paper <- function(si) {
    A <- e$AUCs[[si]]; R <- e$AUCRs[[si]]
    H <- e$AUCCRhighs[[si]]; M <- e$AUCCRmeds[[si]]; L <- e$AUCCRlos[[si]]
    c(cAUC = A$cAUC, rAUC = A$rAUC, lAUC = A$rtAUC, zSlope = unname(e$av.zroc[si, 2]),
      cAUC_sig = R$cAUCO, cAUC_noise = R$cAUCN, lAUC_sig = R$rtAUCO, lAUC_noise = R$rtAUCN,
      lAUC_sig_c3 = H$rtAUCO, lAUC_noise_c3 = H$rtAUCN, lAUC_sig_c2 = M$rtAUCO,
      lAUC_noise_c2 = M$rtAUCN, lAUC_sig_c1 = L$rtAUCO, lAUC_noise_c1 = L$rtAUCN)
  }
  n <- 0
  for (si in 1:11) {
    p <- paper(si)
    ours <- st$data[st$subjects == si][match(names(p), st$stat[st$subjects == si])]
    expect_lt(roc_maxdiff(ours, p), 1e-12)
    n <- n + length(p)
  }
  expect_equal(n, 154)
})

test_that("post_predict: per-draw statistics, intervals and inside", {
  s <- roc_sim_data()
  ref <- roc_ref()$sim
  x <- croc(s$data, "S", "new", post_predict = s$pp)
  l <- lroc(s$data, "S", "new", post_predict = s$pp)
  rk <- rroc(s$data, "S", "new", post_predict = s$pp, type = "latency", conf = 1:3)
  st <- rbind(x$stats, l$stats, rk$stats)
  r <- ref$summary[match(st$stat, ref$summary$stat), ]
  for (v in c("pred_median", "pred_lo", "pred_hi")) expect_lt(roc_maxdiff(st[[v]], r[[v]]), 1e-12)
  expect_identical(st$inside, r$inside)
  expect_true(all(st$n_draws == 5))
  # per draw, against the reference and by hand
  expect_equal(x$pred$stats$postn, 1:5)
  for (nm in c("cAUC", "zSlope", "zCurv"))
    expect_lt(roc_maxdiff(x$pred$stats[[nm]], ref$pred_stats[, nm]), 1e-12)
  expect_lt(roc_maxdiff(l$pred$stats$lAUC, ref$pred_stats[, "lAUC"]), 1e-12)
  v <- x$pred$stats$cAUC
  row <- x$stats[x$stats$stat == "cAUC", ]
  expect_equal(row$pred_median, median(v))
  expect_equal(c(row$pred_lo, row$pred_hi), unname(quantile(v, c(.025, .975))))
  expect_equal(row$inside, row$data >= row$pred_lo && row$data <= row$pred_hi)
  # other interval probabilities
  x50 <- croc(s$data, "S", "new", post_predict = s$pp, probs = c(.25, .75))
  expect_equal(x50$stats$pred_lo[1], unname(quantile(x50$pred$stats$cAUC, .25)))
  # per-draw curves: draw 2 alone gives the same curve as the second draw
  d2 <- croc(s$pp[s$pp$postn == 2, ], "S", "new")
  expect_equal(x$pred$curve[x$pred$curve$postn == 2, c("k", "FAR", "HR")], d2$curve[, c("k", "FAR", "HR")],
               ignore_attr = TRUE)
})

test_that("by = 'subjects' equals the per-subject calls", {
  s <- roc_sim_data()
  calls <- list(function(d, p, ...) croc(d, "S", "new", post_predict = p, ...),
                function(d, p, ...) lroc(d, "S", "new", post_predict = p, n_quantiles = 10, ...),
                function(d, p, ...) rroc(d, "S", "new", post_predict = p, ...),
                function(d, p, ...) rroc(d, "S", "new", post_predict = p, type = "latency", conf = 2:3, ...))
  for (f in calls) {
    b <- f(s$data, s$pp, by = "subjects")
    expect_equal(sort(unique(as.character(b$stats$subjects))), c("1", "2"))
    for (sub in 1:2) {
      one <- f(s$data[s$data$subjects == sub, ], s$pp[s$pp$subjects == sub, ])
      bs <- b$stats[b$stats$subjects == sub, names(one$stats)]
      expect_equal(bs, one$stats, ignore_attr = TRUE)
      bc <- b$curve[b$curve$subjects == sub, names(one$curve)]
      expect_equal(bc, one$curve, ignore_attr = TRUE)
      bp <- b$pred$stats[b$pred$stats$subjects == sub, names(one$pred$stats)]
      expect_equal(bp, one$pred$stats, ignore_attr = TRUE)
    }
  }
})

test_that("bounds, zroc() order, several signal levels, don't know, errors", {
  s <- roc_sim_data()
  dat <- s$data
  # bounds restrict data and predictions alike (LC <= rt <= UC)
  bnd <- data.frame(subjects = 1:2, LC = c(.4, .45), UC = c(1.2, 1.5))
  lc <- bnd$LC[dat$subjects]; uc <- bnd$UC[dat$subjects]
  dw <- dat[dat$rt >= lc & dat$rt <= uc, ]
  pw <- s$pp[s$pp$rt >= bnd$LC[s$pp$subjects] & s$pp$rt <= bnd$UC[s$pp$subjects], ]
  expect_equal(lroc(dat, "S", "new", post_predict = s$pp, bounds = bnd)$stats,
               lroc(dw, "S", "new", post_predict = pw)$stats)
  expect_equal(croc(dat, "S", "new", bounds = c(.4, 1.2))$stats,
               croc(dat[dat$rt >= .4 & dat$rt <= 1.2, ], "S", "new")$stats)
  expect_error(croc(dat, "S", "new", bounds = bnd[1, c("subjects", "LC", "UC")][0, ]), "bounds")
  # censored trials (R, RR NA) drop out; for the latency ROCs also a missing rt
  dc <- dat; dc$R[1:5] <- NA; dc$RR[1:5] <- NA; dc$rt[6:8] <- NA
  expect_equal(croc(dc, "S", "new")$stats$n[1], nrow(dat) - 5)
  expect_equal(lroc(dc, "S", "new")$stats$n[1], nrow(dat) - 8)
  # zroc() has the same points in the opposite order
  x <- croc(dat, "S", "new", by = "subjects")
  z <- zroc(dat, "S", "new", by = "subjects")
  for (sub in 1:2) {
    expect_equal(x$zroc$zF[x$zroc$subjects == sub], rev(z$zF[z$cell == sub]))
    expect_equal(x$zroc$zH[x$zroc$subjects == sub], rev(z$zH[z$cell == sub]))
  }
  # the classification point: hit and false-alarm rates of the signal response; rAUC from them
  h <- mean(dat$R[dat$S == "old"] == "old"); f <- mean(dat$R[dat$S == "new"] == "old")
  x <- croc(dat, "S", "new")
  expect_equal(c(x$curve$FAR[4], x$curve$HR[4]), c(f, h))
  expect_equal(x$stats$data[x$stats$stat == "rAUC"], f * h / 2 + (1 - f) * (h + 1) / 2)
  expect_equal(x$curve$response, c(NA, "old.3", "old.2", "old.1", "new.1", "new.2", "new.3"))
  # signal_response = the first level reverses the roles; areas of a mirrored curve are equal
  y <- croc(dat, "S", "old", signal_response = "new")
  expect_equal(y$stats$data[y$stats$stat == "cAUC"], x$stats$data[x$stats$stat == "cAUC"])
  # several signal levels: each against the noise trials, as zroc()
  d3 <- dat; d3$S3 <- factor(ifelse(d3$S == "new", "new", ifelse(seq_len(nrow(d3)) %% 2, "weak", "strong")))
  x3 <- croc(d3, "S3", "new")
  expect_equal(x3$signal_levels, c("strong", "weak"))
  expect_equal(x3$signal_response, "old")
  expect_equal(x3$stats$data[x3$stats$signal == "weak"],
               croc(d3[d3$S3 != "strong", ], "S3", "new", signal_response = "old")$stats$data)
  z3 <- zroc(d3, "S3", "new")
  expect_equal(x3$zroc$zH[x3$zroc$signal == "weak"], rev(z3$zH[z3$signal == "weak"]))
  # don't know: 2K - 1 responses, no classification point, zROC as zroc()
  dk <- dat; dk$R[dk$RR == 1] <- "new"
  xd <- croc(dk, "S", "new", dont_know = TRUE)
  expect_equal(nrow(xd$curve), 6)
  expect_true(is.na(xd$stats$data[xd$stats$stat == "rAUC"]))
  expect_equal(xd$zroc$zF, rev(zroc(dk, "S", "new", dont_know = TRUE)$zF))
  expect_error(rroc(dat, "S", "new", conf = 2), "latency")
  expect_error(rroc(dat, "S", "new", type = "latency", conf = 4), "conf")
  expect_error(croc(dat, "S", "nothing"), "no trials")
  expect_error(croc(dat, "S", "new", post_predict = dat), "postn")
  expect_output(print(x), "C-ROC")
})

test_that("plots run on a pdf device", {
  s <- roc_sim_data()
  x <- croc(s$data, "S", "new", post_predict = s$pp)
  xb <- croc(s$data, "S", "new", post_predict = s$pp, by = "subjects")
  l <- lroc(s$data, "S", "new", post_predict = s$pp)
  lb <- lroc(s$data, "S", "new", by = "subjects")
  rc <- rroc(s$data, "S", "new", post_predict = s$pp)
  rl <- rroc(s$data, "S", "new", post_predict = s$pp, type = "latency")
  rk <- rroc(s$data, "S", "new", post_predict = s$pp, type = "latency", conf = 1:3)
  d3 <- s$data; d3$S3 <- factor(ifelse(d3$S == "new", "new", ifelse(seq_len(nrow(d3)) %% 2, "weak", "strong")))
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot_croc(x))
  expect_silent(plot_croc(x, z = TRUE))
  expect_silent(plot_croc(xb, z = TRUE))
  expect_silent(plot_croc(croc(s$data, "S", "new")))
  p3 <- s$pp; p3$S3 <- factor(ifelse(p3$S == "new", "new", ifelse(seq_len(nrow(p3)) %% 2, "weak", "strong")))
  x3 <- croc(d3, "S3", "new", post_predict = p3)
  expect_true(all(x3$stats$n_draws == 5))
  expect_silent(plot_croc(x3))
  expect_silent(plot_croc(x3, z = TRUE))
  expect_silent(plot_lroc(lroc(d3, "S3", "new", post_predict = p3)))
  expect_silent(plot_rroc(rroc(d3, "S3", "new", post_predict = p3)))
  expect_silent(plot_roc_auc(x3))
  expect_silent(plot_lroc(l))
  expect_silent(plot_lroc(lb))
  expect_silent(plot_rroc(rc))
  expect_silent(plot_rroc(rl))
  expect_silent(plot_rroc(rk, conf = 3))
  expect_silent(plot_rroc(list(rc, rl, rk), type = "latency"))
  expect_silent(plot_rroc(list(rc, rl, rk), type = "latency", conf = 1))
  expect_silent(plot_roc_auc(list(x, l, rc, rl)))
  expect_silent(plot_roc_auc(list(x, l), stats = c("rAUC", "lAUC", "cAUC"), labels = c("random", "L", "C")))
  expect_silent(plot_roc_auc(xb, stats = "zSlope", h = 1))
  expect_silent(plot_roc_auc(list(xb, croc(s$data, "S", "new", by = "subjects", post_predict = s$pp))))
  expect_error(plot_roc_auc(list(x, xb)), "same by")
  expect_error(plot_croc(l), "croc")
  expect_silent(plot_zroc(s$data, "S", "new", by = "subjects", post_predict = s$pp))
  out <- plot_zroc(s$data, "S", "new", post_predict = s$pp)
  expect_equal(out$data, zroc(s$data, "S", "new"))
  expect_length(out$predicted, 5)
})

test_that("by levels keep the order of the factor's levels", {
  s <- roc_sim_data()
  d <- s$data
  d$subjects <- factor(ifelse(d$subjects == 1, "10", "2"), levels = c("2", "10"))
  x <- croc(d, "S", "new", by = "subjects")
  expect_equal(as.character(unique(x$stats$subjects)), c("2", "10"))
  expect_equal(as.character(unique(x$curve$subjects)), c("2", "10"))
})
