make_ss_forstmann_data <- function(seed = 123, slope = 8, slowdown_factor = 0.50) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)

  dat <- forstmann
  dat$SSD <- Inf

  # SSD grid in seconds, starting at 0.25
  ssd_grid <- seq(0.25, 0.60, by = 0.05)
  ssd_min <- min(ssd_grid)
  ssd_mid <- mean(range(ssd_grid))

  subject_ids <- unique(dat$subjects)

  for (s in subject_ids) {
    idx_s <- which(dat$subjects == s)
    n_s <- length(idx_s)
    if (n_s == 0) next

    # 25% stop trials (at least 1)
    n_stop_s <- max(1, floor(0.25 * n_s))
    stop_idx_s <- sort(sample(idx_s, n_stop_s, replace = FALSE))

    # Assign SSDs cycling so each subject includes 0.25 first
    ssd_assign <- rep(ssd_grid, length.out = n_stop_s)
    dat$SSD[stop_idx_s] <- ssd_assign

    # Response probability increases with SSD (around 50% overall)
    a <- -slope * ssd_mid
    p_resp <- plogis(a + slope * dat$SSD[stop_idx_s])
    responded <- as.logical(rbinom(n_stop_s, 1, p_resp))

    # Successful stops: set R and rt to NA
    succ_s_idx <- stop_idx_s[!responded]
    if (length(succ_s_idx) > 0) {
      dat$R[succ_s_idx]  <- NA
      dat$rt[succ_s_idx] <- NA
    }

    # Failed stops: RT slower for longer SSDs
    fail_s_idx <- stop_idx_s[responded]
    if (length(fail_s_idx) > 0) {
      ssd_vals <- dat$SSD[fail_s_idx]
      slowdown <- slowdown_factor * (ssd_vals - ssd_min)
      dat$rt[fail_s_idx] <- dat$rt[fail_s_idx] + slowdown
    }
  }

  dat
}

# Fixed trials cover go responses, failed stops, stop nonresponses and go failure.
ss_test_data <- function() {
  data.frame(
    subjects = factor(rep(1, 8)), trials = 1:8,
    S = factor(rep(c("left", "right"), 4)),
    R = factor(c("left", "left", "left", "left", NA, NA, NA, "right"),
               levels = c("left", "right")),
    rt = c(.65, .8, .65, .8, NA, NA, NA, .9),
    SSD = c(Inf, Inf, .3, .3, .2, .4, Inf, Inf)
  )
}

ss_test_design <- function(go, model = stop_signal(go = go)) {
  go_formula <- if (go == "exgaussian") {
    list(mu ~ lM, sigma ~ 1, tau ~ 1)
  } else {
    list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1)
  }
  design(model = model,
         factors = list(subjects = 1, S = c("left", "right")),
         Rlevels = c("left", "right"), matchfun = function(d) d$S == d$lR,
         formula = c(go_formula, list(muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)),
         report_p_vector = FALSE)
}

make_stop_cpp_p <- function(design, values = c(
  mu = log(.55), mu_lMTRUE = log(.85), sigma = log(.07), tau = log(.15),
  v = log(1.8), v_lMTRUE = log(.35), B = log(.8), A = log(.1), t0 = log(.2), s = 0,
  muS = log(.24), sigmaS = log(.06), tauS = log(.06),
  gf = qnorm(.05), tf = qnorm(.05)
)) {
  p <- sampled_pars(design, doMap = FALSE)
  missing <- setdiff(names(p), names(values))
  if (length(missing)) stop("Missing values for: ", paste(missing, collapse = ", "))
  p[] <- values[names(p)]
  p
}

expect_stop_cpp_equal_r <- function(design, p, tolerance = 1e-5) {
  withr::local_options(emc.ll_backend = "multiprocess")
  dat <- ss_test_data()
  emc <- make_emc(dat, design, type = "single", n_chains = 1,
                  compress = FALSE, verbose = FALSE)
  model_c <- emc[[1]]$model
  model_r_list <- model_c()
  model_r_list$c_name <- NULL
  model_r <- function() model_r_list

  proposals <- rbind(baseline = p, go_changed = p, stop_changed = p)
  go_parameter <- if ("mu" %in% names(p)) "mu" else "v"
  stop_parameter <- intersect(c("muS", "meanlogS", "shapeS"), names(p))
  proposals[2, go_parameter] <- p[go_parameter] + .1
  proposals[3, stop_parameter] <- p[stop_parameter] + .1

  ll_c <- calc_ll_manager(proposals, emc[[1]]$data[[1]], model_c,
                          return_trialwise = TRUE)
  # The R route returns totals, so evaluate each fixed trial separately.
  ll_r <- t(vapply(seq_len(nrow(dat)), function(i) {
    one <- make_emc(dat[i, , drop = FALSE], design, type = "single",
                    n_chains = 1, compress = FALSE, verbose = FALSE)
    as.numeric(calc_ll_manager(proposals, one[[1]]$data[[1]], model_r))
  }, numeric(nrow(proposals))))

  expect_equal(dim(ll_c), c(nrow(dat), nrow(proposals)))
  expect_true(all(is.finite(ll_c)) && all(is.finite(ll_r)))
  expect_gt(min(ll_c), log(1e-10))
  expect_lt(max(abs(ll_c - ll_r)), tolerance)
  total <- as.numeric(calc_ll_manager(proposals, emc[[1]]$data[[1]], model_c))
  expect_equal(unname(colSums(ll_c)), total, tolerance = 1e-10)

  go <- is.infinite(dat$SSD)
  stop <- !go
  expect_gt(max(abs(ll_c[go, 2] - ll_c[go, 1])), 1e-8)
  expect_gt(max(abs(ll_c[stop, 3] - ll_c[stop, 1])), 1e-8)
  expect_equal(ll_c[go, 3], ll_c[go, 1], tolerance = 1e-12)
  # Go-trial nonresponse is gf; stop-trial nonresponse also includes inhibition.
  expect_equal(as.numeric(ll_c[7, ]), rep(pnorm(p[["gf"]], log.p = TRUE), 3),
               tolerance = 1e-10)
  expect_true(all(ll_c[c(5, 6), 1] > ll_c[7, 1]))
}
