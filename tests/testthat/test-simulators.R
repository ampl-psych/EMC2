simulation_case <- function(model, formula, centre, levels = c("left", "right"), functions = NULL) {
  des <- design(model = model, factors = list(subjects = 1),
    Rlevels = levels, formula = formula, functions = functions)
  list(design = des, centre = centre, levels = levels)
}
cases <- list(
  LNR = simulation_case(LNR, list(m ~ 0 + lR, s ~ 1, t0 ~ 1),
    c(m_lRleft = -1, m_lRright = -.5, s = log(.4), t0 = log(.2))),
  LBA = simulation_case(LBA, list(v ~ 0 + lR, sv ~ 1, B ~ 1, A ~ 1, t0 ~ 1),
    c(v_lRleft = 3, v_lRright = 2, sv = 0, B = log(.5), A = log(.2), t0 = log(.2))),
  RDM = simulation_case(RDM, list(v ~ 0 + lR, s ~ 1, B ~ 1, A ~ 1, t0 ~ 1),
    c(v_lRleft = log(3), v_lRright = log(2), s = 0, B = log(.5), A = log(.2), t0 = log(.2))),
  WDM = simulation_case(DDM, list(v ~ 1, a ~ 1, Z ~ 1, t0 ~ 1),
    c(v = 1, a = log(1.5), Z = qnorm(.4), t0 = log(.2))),
  DDM = simulation_case(DDM, list(v ~ 1, a ~ 1, Z ~ 1, t0 ~ 1, sv ~ 1, SZ ~ 1, st0 ~ 1),
    c(v = 1, a = log(1.5), Z = qnorm(.4), t0 = log(.2), sv = log(.5), SZ = qnorm(.2), st0 = log(.05))),
  SSEXG = simulation_case(EMC2:::SSEXG,
    list(mu ~ 0 + lR, sigma ~ 1, tau ~ 1, muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1),
    c(mu_lRleft = log(.5), mu_lRright = log(.65), sigma = log(.05), tau = log(.15),
      muS = log(.2), sigmaS = log(.03), tauS = log(.05), gf = qnorm(.1), tf = qnorm(.1)),
    functions = list(SSD = function(d) rep(.25, nrow(d)))),
  SSRDEX = simulation_case(EMC2:::SSRDEX,
    list(v ~ 0 + lR, s ~ 1, B ~ 1, A ~ 1, t0 ~ 1, muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1),
    c(v_lRleft = log(3), v_lRright = log(2), s = 0, B = log(.5), A = log(.2), t0 = log(.2),
      muS = log(.2), sigmaS = log(.03), tauS = log(.05), gf = qnorm(.1), tf = qnorm(.1)),
    functions = list(SSD = function(d) rep(.25, nrow(d)))),
  ordered_probit = simulation_case(ordered_probit, list(location ~ 1, scale ~ 1, cut ~ 1),
    c(location = .3, scale = 0, cut_lRleft = -.5, cut_lRmiddle = log(1.5)),
    levels = c("left", "middle", "right")),
  ordered_logit = simulation_case(ordered_logit, list(location ~ 1, scale ~ 1, cut ~ 1),
    c(location = .3, scale = 0, cut_lRleft = -.5, cut_lRmiddle = log(1.5)),
    levels = c("left", "middle", "right")),
  multinomial_logit = simulation_case(multinomial_logit, list(utility ~ 0 + lR),
    c(utility_lRleft = .5, utility_lRmiddle = 0, utility_lRright = -.5),
    levels = c("left", "middle", "right")),
  multinomial_probit = simulation_case(multinomial_probit, list(utility ~ 0 + lR),
    c(utility_lRleft = .5, utility_lRmiddle = 0, utility_lRright = -.5),
    levels = c("left", "middle", "right"))
)
# Direct model density/CDF evaluation at the same mapped coefficients.
reference <- function(case, name) {
  template <- data.frame(subjects = factor(1), R = factor(case$levels, levels = case$levels), rt = .6)
  emc <- make_emc(template, case$design, type = "single",
    n_chains = 1, compress = FALSE, rt_resolution = 1e-12)
  dadm <- emc[[1]]$data[[1]]
  model <- emc[[1]]$model()
  pars <- EMC2:::get_pars_matrix_oo(case$centre, dadm, model)
  if (name %in% c("ordered_probit", "ordered_logit", "multinomial_logit", "multinomial_probit")) {
    probs <- vapply(case$levels, function(r) {
      d <- dadm[dadm$R == r, , drop = FALSE]
      p <- pars[dadm$R == r, , drop = FALSE]
      attr(p, "ok") <- rep(TRUE, nrow(p))
      attr(d, "expand") <- 1L
      exp(model$log_likelihood(p, d, model))
    }, numeric(1))
    return(data.frame(response = case$levels, cutoff = Inf, expected = probs))
  }
  cuts <- c(.4, .6, 1, 2, Inf)
  levels <- case$levels
  grid <- expand.grid(response = levels, cutoff = cuts, stringsAsFactors = FALSE)
  p <- if (name %in% c("WDM", "DDM")) pars[1, , drop = FALSE] else pars[seq_along(case$levels), , drop = FALSE]
  grid$expected <- mapply(function(r, cutoff) {
    if (name %in% c("WDM", "DDM")) {
      return(model$pfun(cutoff, factor(r, levels = case$levels), p))
    }
    k <- match(r, case$levels)
    density <- function(t) {
      winner <- p[rep(k, length(t)), , drop = FALSE]
      loser <- p[rep(3 - k, length(t)), , drop = FALSE]
      if (name %in% c("SSEXG", "SSRDEX")) {
        return(model$dfunG(t, winner) * (1 - model$pfunG(t, loser)) *
          (1 - winner[, "gf"]) * (winner[, "tf"] +
          (1 - winner[, "tf"]) * (1 - model$pfunS(t, winner))))
      }
      model$dfun(t, winner) * (1 - model$pfun(t, loser))
    }
    integrate(density, 0, cutoff, rel.tol = 1e-8, abs.tol = 1e-9)$value
  }, grid$response, grid$cutoff)
  if (name %in% c("SSEXG", "SSRDEX")) {
    grid <- rbind(grid, data.frame(response = "none", cutoff = Inf,
      expected = p[1, "gf"] + (1 - p[1, "gf"]) * (1 - p[1, "tf"]) * model$sfun(p, 2)))
  }
  grid
}
# Six worst-case binomial standard errors; no dependence on exact RNG draws.
n <- if (nzchar(Sys.getenv("CI")) || !identical(Sys.getenv("NOT_CRAN"), "true")) 10000L else 25000L
for (name in names(cases)) {
  test_that(paste(name, "simulator agrees with response and RT probabilities"), {
    withr::local_seed(123, .rng_kind = "Mersenne-Twister", .rng_normal_kind = "Inversion")
    case <- cases[[name]]
    ref <- reference(case, name)
    dat <- make_data(case$centre, case$design, n_trials = n)
    observed <- mapply(function(r, t) {
      if (r == "none") return(mean(is.na(dat$rt)))
      if (!"rt" %in% names(dat)) return(mean(dat$R == r))
      mean(!is.na(dat$rt) & dat$R == r & dat$rt <= t)
    }, ref$response, ref$cutoff)
    expect_lt(max(abs(observed - ref$expected)), 3 / sqrt(n))
  })
}
