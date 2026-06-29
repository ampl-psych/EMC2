RNGkind("L'Ecuyer-CMRG")
set.seed(123)

make_ssrdex_cpp_design <- function() {
  design(
    model = function() stop_signal(go = "racing_diffusion", stop = "exgaussian"),
    factors = list(subjects = 1, S = c("left", "right")),
    Rlevels = c("left", "right"),
    matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR),
    formula = list(v ~ lM, B ~ 1, A ~ 1, t0 ~ 1, s ~ 1,
                   muS ~ 1, sigmaS ~ 1, tauS ~ 1, gf ~ 1, tf ~ 1)
  )
}

make_stop_cpp_p <- function(design, values) {
  p <- sampled_pars(design, doMap = FALSE)
  missing <- setdiff(names(p), names(values))
  if (length(missing)) stop("Missing values for: ", paste(missing, collapse = ", "))
  p[] <- values[names(p)]
  p
}

expect_stop_cpp_equal_r <- function(design, p, tolerance = 1e-5) {
  dat <- make_data(p, design, n_trials = 30, functions = list(SSD = make_ssd()))
  emc <- make_emc(dat, design, type = "single", n_chains = 1, compress = FALSE)

  proposals <- rbind(
    p,
    p + rep(c(-0.003, 0.002, 0.001), length.out = length(p)),
    p + rep(c(0.002, -0.001, 0.003), length.out = length(p))
  )
  colnames(proposals) <- names(p)

  model_c <- emc[[1]]$model
  model_r_list <- model_c()
  model_r_list$c_name <- NULL
  model_r <- function() model_r_list

  ll_c <- calc_ll_manager(proposals, emc[[1]]$data[[1]], model_c)
  ll_r <- calc_ll_manager(proposals, emc[[1]]$data[[1]], model_r)

  expect_equal(unname(ll_c), unname(ll_r), tolerance = tolerance)
}

test_that("SSRDEX C++ likelihood matches R reference path", {
  design_ssrdex <- make_ssrdex_cpp_design()
  p_ssrdex <- make_stop_cpp_p(
    design_ssrdex,
    c(
      v = log(1.8),
      v_lMTRUE = log(.35),
      B = log(.8),
      A = log(.10),
      t0 = log(.20),
      s = log(1),
      muS = log(.17),
      sigmaS = log(.025),
      tauS = log(.05),
      gf = qnorm(.1),
      tf = qnorm(.1)
    )
  )

  expect_identical(design_ssrdex$model()$c_name, "SSRDEX")
  expect_stop_cpp_equal_r(design_ssrdex, p_ssrdex)
})
