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

