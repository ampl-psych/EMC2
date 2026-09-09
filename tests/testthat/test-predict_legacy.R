# predict() on objects fitted with older EMC2 versions, and the row-preservation
# guarantees it relies on. See update2version() / refresh_model() in design.R.

# Emulate a fit whose model closure was frozen before May 2025: rfun(lR, pars)
make_stale <- function(emc) {
  pr <- get_prior(emc)
  des <- attr(pr, "design")
  ml <- des[[1]]$model()
  ml$rfun <- function(lR, pars) rLNR(lR, pars, ok = attr(pars, "ok"))
  ml$bound$minmax[, "t0"] <- c(0.1, Inf)          # a fit-time custom bound to preserve
  des[[1]]$model <- function() ml
  attr(pr, "design") <- des
  emc <- lapply(emc, function(x) { x$prior <- pr; x })
  emc[[1]]$model <- function() ml
  class(emc) <- "emc"
  emc
}

test_that("stale model closures are detected and refreshed for predict()", {
  stale <- make_stale(samples_LNR)
  expect_false(emc_has_stale_model(samples_LNR))
  expect_true(emc_has_stale_model(stale))
  # get_design() hands back a usable design without touching the object
  d <- get_design(stale)
  expect_identical(names(formals(d[[1]]$model()$rfun))[1], "data")
  expect_equal(d[[1]]$model()$bound$minmax[, "t0"], c(0.1, Inf))
  expect_true(emc_has_stale_model(stale))
  set.seed(1)
  expect_message(p <- predict(stale, n_post = 2), "update2version")
  expect_equal(nrow(p), 2 * nrow(get_data(stale)))
})

test_that("update2version() keeps one predicted row per trial and the stored bound", {
  for (obj in list(samples_LNR, make_stale(samples_LNR))) {
    up <- update2version(obj)
    expect_false(emc_has_stale_model(up))
    expect_equal(nrow(get_data(up)), nrow(get_data(samples_LNR)))
    expect_equal(attr(up[[1]]$data[[1]], "expand"), attr(samples_LNR[[1]]$data[[1]], "expand"))
    set.seed(2)
    expect_equal(nrow(predict(up, n_post = 1)), nrow(get_data(up)))
  }
  up <- update2version(make_stale(samples_LNR))
  expect_equal(get_design(up)[[1]]$model()$bound$minmax[, "t0"], c(0.1, Inf))
})

test_that("update2version() converts the pre-2025 all-rows expand convention", {
  x <- samples_LNR[[1]]$data[[1]]
  new_exp <- attr(x, "expand")
  # Rebuild the old index: every row of the compressed dadm, trial by trial
  win_rows <- which(x$winner)
  old_exp <- as.vector(vapply(new_exp, function(w) {
    tr <- x$trials[win_rows[w]]; s <- x$subjects[win_rows[w]]
    which(x$trials == tr & x$subjects == s)
  }, numeric(length(unique(x$lR)))))
  expect_gt(max(old_exp), sum(x$winner))
  old <- samples_LNR
  attr(old[[1]]$data[[1]], "expand") <- old_exp
  up <- update2version(old)
  expect_equal(attr(up[[1]]$data[[1]], "expand"), new_exp)
})

test_that("rfuns return one row per trial with NA for out-of-bound trials", {
  lR <- factor(rep(c("left", "right"), 5), levels = c("left", "right"))
  # LNR / RDM: 2 accumulators x 5 trials; trial 3 out of bounds
  ok <- rep(TRUE, 10); ok[5:6] <- FALSE
  p_lnr <- cbind(m = rep(0, 10), s = 1, t0 = .2)
  r <- rLNR(lR, p_lnr, ok = ok)
  expect_equal(nrow(r), 5); expect_true(is.na(r$rt[3])); expect_true(is.na(r$R[3]))
  expect_true(all(is.finite(r$rt[-3])))
  p_rdm <- cbind(v = rep(2, 10), B = 1, A = 0, t0 = .2, s = 1)
  r <- rRDM(lR, p_rdm, ok = ok)
  expect_equal(nrow(r), 5); expect_true(is.na(r$rt[3])); expect_true(is.na(r$R[3]))
  # DDM: one row per trial, trials 2 and 4 out of bounds
  R <- factor(rep("left", 5), levels = c("left", "right"))
  p_ddm <- cbind(v = rep(1, 5), a = 1, sv = 0, t0 = .2, st0 = 0, s = 1, Z = .5, SZ = 0)
  r <- rDDM(R, p_ddm, ok = c(TRUE, FALSE, TRUE, FALSE, TRUE))
  expect_equal(nrow(r), 5)
  expect_equal(is.na(r$rt), c(FALSE, TRUE, FALSE, TRUE, FALSE))
  expect_equal(levels(r$R), c("left", "right"))
  expect_error(rDDM(factor(1:3), p_ddm[1:3, ]), "exactly two response levels")
})

test_that("design() refuses a DDM with more than two response levels", {
  dat <- forstmann
  levels(dat$R) <- c("left", "right")
  dat$R <- factor(as.character(dat$R), levels = c("left", "right", "none"))
  expect_error(design(data = dat, model = DDM, formula = list(v ~ 1, a ~ 1, t0 ~ 1)),
               "exactly two response levels")
  dat$R <- droplevels(dat$R)
  expect_s3_class(design(data = dat, model = DDM, formula = list(v ~ 1, a ~ 1, t0 ~ 1),
                         report_p_vector = FALSE), "emc.design")
})

test_that("make_data() checks bounds before Ttransform, as the sampler does", {
  des <- design(data = forstmann, model = DDM, formula = list(v ~ 1, a ~ 1, t0 ~ 1, Z ~ 1, SZ ~ 1),
                report_p_vector = FALSE)
  # SZ = .015 is inside the (.01, .99) bound the sampler checks, but the derived
  # 2 * SZ * min(Z, 1 - Z) = .009 fails it if checked after Ttransform
  p <- c(v = 1, a = log(1), t0 = log(.2), Z = qnorm(.3), SZ = qnorm(.015))
  d <- make_data(p, des, data = forstmann)
  expect_s3_class(d, "data.frame")
  expect_equal(nrow(d), nrow(forstmann))
})

test_that("rDDM() refuses a/s beyond what rWDM can simulate instead of hanging", {
  R <- factor(rep("left", 6), levels = c("left", "right"))
  p <- cbind(v = rep(1, 6), a = 1, sv = 0.5, t0 = .2, st0 = 0, s = c(1, 1, 1, 1e-4, 1e-10, 1), Z = .5, SZ = 0.1)
  setTimeLimit(elapsed = 20, transient = TRUE)
  expect_warning(r <- rDDM(R, p), "a/s > 1000")
  setTimeLimit()
  expect_equal(nrow(r), 6)
  expect_equal(is.na(r$rt), c(FALSE, FALSE, FALSE, TRUE, TRUE, FALSE))
})
