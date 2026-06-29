# Tests for the full RDM go/nogo likelihood (G-3): the C++ path (RDMGNG) scores
# observed go trials with the ordinary race density and withheld (inhibited) nogo
# trials with the integral P(no-go finishes first), wired through CensorSpec
# (missingness == 4). Validated against an independent R reference built from the
# model's own dRDM/pRDM.

gng_rdm_design <- function() {
  design(factors = list(subjects = 1, S = c("left", "right")), Rlevels = c("go", "nogo"),
         formula = list(v ~ lR, B ~ lR, A ~ 1, t0 ~ 1), constants = c(s = log(1)),
         model = RDM)
}

# Independent R likelihood: go trials -> race density (winner d x loser survivors);
# withheld trials -> integral of f_nogo(t) * prod_go S_go(t) over [0, cap], using
# the model's own dRDM/pRDM. cap matches the C++ GNG_UPPER_CAP (30).
ref_rdmgng_ll <- function(p, dadm, model, cap = 30, min_ll = log(1e-10)) {
  pars <- EMC2:::get_pars_matrix_oo(p, dadm, model)
  mk <- function(r) matrix(pars[r, ], nrow = 1, dimnames = list(NULL, colnames(pars)))
  ll <- 0
  for (tr in unique(dadm$trials)) {
    idx  <- which(dadm$trials == tr)
    lev  <- as.character(dadm$lR[idx])
    win  <- dadm$winner[idx]
    rt   <- dadm$rt[idx]
    miss <- dadm$missingness[idx][1]
    if (!is.na(miss) && miss == 4) {
      ng <- which(lev == "nogo")
      f <- function(tt) vapply(tt, function(t) {
        d <- dRDM(t, mk(idx[ng]))
        for (g in setdiff(seq_along(idx), ng)) d <- d * (1 - pRDM(t, mk(idx[g])))
        d
      }, numeric(1))
      pr <- stats::integrate(f, 0, cap, rel.tol = 1e-8)$value
      ll <- ll + max(min_ll, log(max(pr, 0)))
    } else {
      wi <- which(win)
      L <- dRDM(rt[wi], mk(idx[wi]))
      for (g in which(!win)) L <- L * (1 - pRDM(rt[g], mk(idx[g])))
      ll <- ll + max(min_ll, log(max(L, 0)))
    }
  }
  ll
}

cpp_ll <- function(emc, p) {
  dadm  <- emc[[1]]$data[[1]]
  model <- emc[[1]]$model()
  p_types <- names(model$p_types)
  designs <- list()
  for (pp in p_types)
    designs[[pp]] <- attr(dadm, "designs")[[pp]][attr(attr(dadm, "designs")[[pp]], "expand"), , drop = FALSE]
  constants <- attr(dadm, "constants"); if (is.null(constants)) constants <- NA
  p_mat <- matrix(p, nrow = 1, dimnames = list(NULL, names(p)))
  EMC2:::calc_ll(p_mat, dadm, constants = constants, designs = designs, type = model$c_name,
                 model$bound, model$transform, model$pre_transform, p_types = p_types,
                 min_ll = log(1e-10), model$trend)
}

test_that("RDMGNG C++ likelihood matches the dRDM/pRDM reference", {
  for (sd in c(1, 7, 13)) {
    set.seed(sd)
    des <- gng_rdm_design()
    p <- sampled_pars(des); p[] <- log(1); p["v_lRnogo"] <- -0.3; p["v"] <- 0.5
    dat <- make_data(p, des, n_trials = 60)
    emc <- make_emc(dat, des, type = "single", compress = FALSE, n_chains = 1)
    expect_gt(sum(emc[[1]]$data[[1]]$missingness == 4, na.rm = TRUE), 0)  # has withheld
    expect_equal(cpp_ll(emc, p), ref_rdmgng_ll(p, emc[[1]]$data[[1]], emc[[1]]$model()),
                 tolerance = 1e-5)
  }
})

test_that("RDMGNG likelihood is finite and varies with parameters", {
  set.seed(2)
  des <- gng_rdm_design()
  p <- sampled_pars(des); p[] <- log(1); p["v_lRnogo"] <- -0.3; p["v"] <- 0.5
  dat <- make_data(p, des, n_trials = 60)
  emc <- make_emc(dat, des, type = "single", compress = FALSE, n_chains = 1)
  ll0 <- cpp_ll(emc, p)
  p2 <- p; p2["v_lRnogo"] <- 0.5     # different no-go drift -> different likelihood
  expect_true(is.finite(ll0))
  expect_false(isTRUE(all.equal(ll0, cpp_ll(emc, p2))))
})

test_that("run_sbc runs end-to-end on an RDM go/nogo model", {
  skip_on_cran()
  set.seed(3)
  des <- gng_rdm_design()
  pars <- sampled_pars(des)
  pri  <- prior(des, type = "single", mu_mean = pars, mu_sd = rep(0.3, length(pars)))
  out <- run_sbc(des, pri, replicates = 3, trials = 80, verbose = FALSE,
                 fileName = tempfile(fileext = ".RData"),
                 cores_per_chain = 1, cores_for_chains = 1,
                 stop_criteria = list(iter = 100, selection = "alpha"))
  expect_true(all(c("rank", "med", "bias", "coverage") %in% names(out)))
  ranks <- out$rank[[1]]
  expect_equal(nrow(ranks), 3L)
  expect_true(all(ranks >= 0 & ranks <= 1, na.rm = TRUE))
})
