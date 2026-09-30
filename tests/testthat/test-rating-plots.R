rating_toy <- function() {
  set.seed(1)
  n <- 400
  data.frame(subjects = factor(rep(1:2, each = n / 2)),
             S = factor(rep(c("n", "s"), n / 2)),
             R = factor(sample(c("a", "b"), n, TRUE)),
             rt = runif(n, .3, 1.5), RR = sample(1:3, n, TRUE))
}

test_that("rating_summary: proportions sum to 1 per subject-averaged cell", {
  dat <- rating_toy()
  s <- rating_summary(dat, factors = "S")
  expect_equal(levels(s$resp), c("a.3", "a.2", "a.1", "b.1", "b.2", "b.3"))
  expect_equal(as.vector(tapply(s$p, s$cell, sum)), c(1, 1))
  # hand check: subject 1, cell n, response a.3
  x <- dat[dat$S == "n", ]
  p1 <- mean(x$R[x$subjects == 1] == "a" & x$RR[x$subjects == 1] == 3)
  p2 <- mean(x$R[x$subjects == 2] == "a" & x$RR[x$subjects == 2] == 3)
  expect_equal(s$p[s$cell == "n" & s$resp == "a.3"], (p1 + p2) / 2)
  q <- sapply(1:2, function(i) quantile(x$rt[x$subjects == i & x$R == "a" & x$RR == 3], .5, names = FALSE))
  expect_equal(s$q50[s$cell == "n" & s$resp == "a.3"], mean(q))
})

test_that("zroc: cumulative proportions from the top of the second response", {
  dat <- rating_toy()
  z <- zroc(dat, "S", "n")
  nz <- dat[dat$S == "n", ]
  f <- EMC2:::rating_fold(nz$R, nz$RR, 3)
  expect_equal(z$F[1], mean(as.integer(f) >= 2))
  expect_equal(z$F[5], mean(as.integer(f) == 6))
  expect_equal(z$zH, qnorm(z$H))
  expect_true(all(diff(z$F) <= 0))
})

test_that("rating plots run with and without predictions", {
  dat <- rating_toy()
  pp <- do.call(rbind, lapply(1:3, function(i) cbind(postn = i, dat[sample(nrow(dat)), ])))
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot_ratings(dat, pp, factors = "S"))
  expect_silent(plot_zroc(dat, "S", "n", post_predict = pp))
})
