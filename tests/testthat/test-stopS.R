RNGkind("L'Ecuyer-CMRG")
set.seed(123)

# Compare the optimized C/C++ SSEXG likelihood to the R reference path by
# removing `c_name` from the model while keeping the same data and proposals.
expect_ssexg_cr_equal <- function(emc, p_vector, n = 25, proposal_var = 0.0025,
                                  tolerance = 1e-5) {
  proposals <- mvtnorm::rmvnorm(
    n,
    mean = p_vector,
    sigma = diag(length(p_vector)) * proposal_var
  )
  proposals <- rbind(p_vector, proposals)
  colnames(proposals) <- names(p_vector)

  model_c <- emc[[1]]$model
  model_r_list <- model_c()
  model_r_list$c_name <- NULL
  model_r <- function() model_r_list

  ll_r <- calc_ll_manager(proposals, emc[[1]]$data[[1]], model_r)
  ll_c <- calc_ll_manager(proposals, emc[[1]]$data[[1]], model_c)

  expect_equal(ll_c, ll_r, tolerance = tolerance)
}

# Baseline single-subject SSEXG fixture using staircase SSD generation.
mySSD_function <- make_ssd()

designSSexG <- design(model=SSEXG,
                      factors=list(subjects=1,S=c("left","right")),Rlevels=c("left","right"),
                      matchfun=function(d) as.numeric(d$S)==as.numeric(d$lR),
                      formula=list(mu~lM,sigma~1,tau~1,muS~1,sigmaS~1,tauS~1, gf~1,tf~1)
)

p_vector <- sampled_pars(designSSexG,doMap = FALSE)
p_vector[1:length(p_vector)] <- c(log(.6), log(.8), log(0.05), log(0.2),
                                  log(0.2),log(0.03), log(0.05),
                                  qnorm(.1),qnorm(.1))

dat <- make_data(p_vector, designSSexG, n_trials = 10,
                 functions = list(SSD = mySSD_function))
emc <- make_emc(dat, designSSexG, type = "single")


test_that("exG", {
  # Smoke-test that SSEXG can initialize and returns structurally valid samples.
  samples <- init_chains(emc, particles = 10, cores_for_chains = 1)[[1]]$samples
  expect_equal(dim(samples$alpha), c(length(p_vector), 1, 1))
  expect_named(samples$alpha[, 1, 1], names(p_vector))
  expect_equal(as.character(samples$stage), "init")
  expect_equal(samples$idx, 1)
  expect_true(is.finite(samples$subj_ll[1, 1]))
})

test_that("SSEXG mapped_pars supplies internal SSD and hides bookkeeping columns", {
  mapped <- mapped_pars(designSSexG, p_vector)

  expect_s3_class(mapped, "data.frame")
  expect_false("SSD" %in% names(mapped))
  expect_false("lI" %in% names(mapped))
  expect_true(all(c("S", "lR", "lM", "mu", "sigma", "tau", "muS", "sigmaS", "tauS", "gf", "tf") %in% names(mapped)))
})

test_that("SSEXG mapped_pars supports the R reference model path", {
  SSEXG_R <- function() {
    model <- SSEXG()
    model$c_name <- NULL
    model
  }
  design_R <- design(model = SSEXG_R,
                     factors = list(subjects = 1, S = c("left", "right")),
                     Rlevels = c("left", "right"),
                     matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR),
                     formula = list(mu ~ lM, sigma ~ 1, tau ~ 1,
                                    muS ~ 1, sigmaS ~ 1, tauS ~ 1,
                                    gf ~ 1, tf ~ 1))

  mapped <- mapped_pars(design_R, p_vector)

  expect_s3_class(mapped, "data.frame")
  expect_false("SSD" %in% names(mapped))
  expect_false("lI" %in% names(mapped))
})

test_that("SSEXG mapped_pars keeps lI when it varies by accumulator", {
  design_lI <- design(model = SSEXG,
                      factors = list(subjects = 1, S = c("left", "right"), lI = 1:2),
                      Rlevels = c("left", "right"),
                      matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR),
                      formula = list(mu ~ lM, sigma ~ 1, tau ~ 1,
                                     muS ~ 1, sigmaS ~ 1, tauS ~ 1,
                                     gf ~ 1, tf ~ 1))
  p_lI <- p_vector
  names(p_lI) <- names(sampled_pars(design_lI, doMap = FALSE))

  mapped <- mapped_pars(design_lI, p_lI)

  expect_false("SSD" %in% names(mapped))
  expect_true("lI" %in% names(mapped))
  expect_equal(sort(unique(as.numeric(as.character(mapped$lI)))), c(1, 2))
})

test_that("SSEXG mapped prior plotting supplies internal SSD", {
  prior_SS <- prior(designSSexG, type = "single")
  samples <- get_objects(design = designSSexG, prior = prior_SS, type = "single",
                         sample_prior = TRUE, selection = "alpha", N = 10)

  expect_error(get_pars(samples, selection = "alpha", type = "single",
                        stage = "sample", map = TRUE), NA)
})



test_that("stop-signal models disable compression and RT binning", {
  tiny_dat <- data.frame(
    subjects = factor(c(1, 1)),
    S = factor(c("left", "left"), levels = c("left", "right")),
    R = factor(c("left", "left"), levels = c("left", "right")),
    rt = c(0.123, 0.178),
    SSD = c(Inf, Inf)
  )

  expect_error(
    EMC2:::design_model(
      tiny_dat, designSSexG,
      compress = TRUE, verbose = FALSE, rt_check = FALSE
    ),
    "Compression is not supported for stop-signal models"
  )

  expect_message(
    emc_stop <- make_emc(
      tiny_dat, designSSexG, type = "single",
      n_chains = 1, compress = TRUE, rt_resolution = 0.05
    ),
    "Stop-signal models do not support compression"
  )

  winner_rows <- emc_stop[[1]]$data[[1]]$winner
  expect_equal(emc_stop[[1]]$data[[1]]$rt[winner_rows], tiny_dat$rt)
})

test_that("staircase resets per subject", {
  # Each subject-specific staircase should start from the same initial SSD.
  design_multi <- design(model = SSEXG,
                         factors = list(subjects = 1:3, S = c("left", "right")),
                         Rlevels = c("left", "right"),
                         matchfun = function(d) as.numeric(d$S) == as.numeric(d$lR),
                         formula = list(mu ~ lM, sigma ~ 1, tau ~ 1,
                                         muS ~ 1, sigmaS ~ 1, tauS ~ 1,
                                         gf ~ 1, tf ~ 1))

  p_multi <- sampled_pars(design_multi, doMap = FALSE)
  p_multi[1:length(p_multi)] <- c(log(.6), log(.8), log(0.05), log(0.2),
                                  log(0.2), log(0.03), log(0.05),
                                  qnorm(.1), qnorm(.1))

  dat_multi <- make_data(p_multi, design_multi, n_trials = 400,
                         functions = list(SSD = make_ssd(p_stop = 1)))

  staircase <- attr(dat_multi, "staircase")
  expect_false(is.null(staircase))
  expect_equal(
    unname(vapply(staircase$specs, `[[`, numeric(1), "SSD0")),
    rep(0.25, length(staircase$specs))
  )
})

test_that("posterior-only stop-signal plot summaries are stable", {
  tiny_ss_data <- data.frame(
    subjects = factor(rep(1:2, each = 8)),
    E = factor(rep(rep(c("speed", "accuracy"), each = 4), 2),
               levels = c("speed", "accuracy")),
    S = factor(rep(c("left", "right"), 8), levels = c("left", "right")),
    SSD = rep(c(0.10, 0.20, 0.30, 0.40), 4),
    R = factor(c(
      "left", NA, "left", NA,
      "right", "right", NA, "right",
      "left", "left", NA, "left",
      NA, "right", "right", NA
    ), levels = c("left", "right")),
    rt = c(
      0.31, NA, 0.35, NA,
      0.33, 0.34, NA, 0.38,
      0.32, 0.33, NA, 0.36,
      NA, 0.35, 0.37, NA
    )
  )

  tiny_post <- do.call(rbind, lapply(1:3, function(i) {
    out <- tiny_ss_data
    out$postn <- i
    out$rt <- out$rt + c(0, 0.01, -0.01)[i]
    if (i == 2) {
      miss <- is.finite(out$SSD) & out$SSD == 0.20
      out$R[miss] <- NA
      out$rt[miss] <- NA
    }
    out[, c("postn", setdiff(names(out), "postn"))]
  }))
  rownames(tiny_post) <- NULL

  draw1 <- data.frame(
    x_plot = c(0.25, 0.50, 0.75, 1.00),
    ssd = c("[0.1,0.2]", "(0.2,0.3]", "(0.3,0.4]", "(0.4,0.5]"),
    srrt = c(0.42, 0.45, 0.48, 0.51)
  )
  draw2 <- data.frame(
    x_plot = c(0.25, 0.75, 1.00),
    ssd = c("[0.1,0.2]", "(0.3,0.4]", "(0.4,0.5]"),
    srrt = c(0.40, 0.47, 0.50)
  )
  postn_list <- list(`1` = list(left = draw1), `2` = list(left = draw2))

  aligned <- get_stop_signal_postn_quantiles(
    postn_list = postn_list,
    level = "left",
    value_col = "srrt",
    quants = c(0.025, 0.975)
  )
  check <- prep_data_plot(
    input = tiny_ss_data,
    post_predict = tiny_post,
    prior_predict = NULL,
    to_plot = "posterior",
    limits = "posterior",
    factors = "E",
    defective_factor = "S",
    subject = NULL,
    n_cores = 1,
    n_post = 3,
    functions = NULL,
    remove_na = FALSE
  )

  summary_lines <- c(
    paste("sources:", paste(unname(check$sources), collapse = ", ")),
    paste("datasets:", paste(names(check$datasets), collapse = ", ")),
    paste("labels:", paste(colnames(aligned), collapse = " | ")),
    paste("x_plot:", paste(sprintf("%.2f", unname(aligned["x_plot", ])), collapse = ", ")),
    paste("median:", paste(sprintf("%.3f", unname(aligned["50%", ])), collapse = ", "))
  )
  expect_snapshot(writeLines(summary_lines))

  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  expect_no_error(plot_ss_if(
    tiny_ss_data,
    post_predict = tiny_post,
    probs = c(0, 0.5, 1),
    use_global_quantiles = TRUE,
    to_plot = "posterior",
    use_lim = "posterior",
    factors = "E",
    within_plot = "S",
    layout = c(1, 2),
    legendpos = c(NA, NA)
  ))
  expect_error(
    plot_ss_if(
      tiny_ss_data,
      probs = c(0, 0.5, 0.5, 1),
      factors = "E",
      within_plot = "S"
    ),
    "`probs` must not contain duplicate values",
    fixed = TRUE
  )
  expect_error(
    plot_ss_srrt(
      tiny_ss_data,
      probs = c(0.2, 0.6, 1),
      factors = "E",
      within_plot = "S"
    ),
    "`probs` must start at 0 and end at 1",
    fixed = TRUE
  )
  expect_message(
    expect_no_error(plot_ss_if(
      tiny_ss_data,
      ssd_binning = "value",
      probs = c(0.5, 0.5),
      factors = "E",
      within_plot = "S",
      layout = c(1, 2),
      legendpos = c(NA, NA)
    )),
    "ignored when `ssd_binning = \"value\"`",
    fixed = TRUE
  )
  expect_error(
    plot_ss_if(
      tiny_ss_data,
      post_predict = tiny_post,
      quants = c(0.025, 0.5, 0.975),
      to_plot = "posterior",
      use_lim = "posterior",
      factors = "E",
      within_plot = "S"
    ),
    "`quants` must be a numeric vector of exactly two finite values",
    fixed = TRUE
  )
  expect_error(
    plot_ss_srrt(
      tiny_ss_data,
      post_predict = tiny_post,
      quants = c(0.6, 0.9),
      to_plot = "posterior",
      use_lim = "posterior",
      factors = "E",
      within_plot = "S"
    ),
    "`quants` must bracket 0.5",
    fixed = TRUE
  )
  expect_no_error(plot_ss_srrt(
    tiny_ss_data,
    post_predict = tiny_post,
    probs = c(0, 0.5, 1),
    use_global_quantiles = TRUE,
    to_plot = "posterior",
    use_lim = "posterior",
    factors = "E",
    within_plot = "S",
    layout = c(1, 2),
    legendpos = c(NA, NA)
  ))
  expect_no_error(plot_ss_if(
    tiny_ss_data,
    post_predict = tiny_post,
    probs = seq(0, 1, length.out = 6),
    use_global_quantiles = FALSE,
    on_duplicate_quantiles = "value",
    to_plot = c("data", "posterior"),
    use_lim = c("data", "posterior"),
    factors = "E",
    within_plot = "S",
    layout = c(1, 2),
    legendpos = c(NA, NA)
  ))
  expect_no_error(plot_ss_if(
    tiny_ss_data,
    post_predict = tiny_post,
    probs = seq(0, 1, length.out = 7),
    use_global_quantiles = TRUE,
    on_duplicate_quantiles = "reduce",
    to_plot = "posterior",
    use_lim = "posterior",
    factors = "E",
    within_plot = "S",
    layout = c(1, 2),
    legendpos = c(NA, NA)
  ))
  tiny_post_mismatched_ssd <- tiny_post
  tiny_post_mismatched_ssd$SSD <- tiny_post_mismatched_ssd$SSD + 0.50
  expect_error(
    plot_ss_if(
      tiny_ss_data,
      post_predict = tiny_post_mismatched_ssd,
      probs = seq(0, 1, length.out = 7),
      use_global_quantiles = TRUE,
      on_duplicate_quantiles = "error",
      to_plot = c("data", "posterior"),
      use_lim = c("data", "posterior"),
      factors = "E",
      within_plot = "S",
      layout = c(1, 2),
      legendpos = c(NA, NA)
    ),
    "different SSD support"
  )
  expect_no_error(plot_ss_if(
    tiny_ss_data,
    ssd_binning = "value",
    ssd_round = 0.10,
    factors = "E",
    within_plot = "S",
    layout = c(1, 2),
    legendpos = c(NA, NA)
  ))

  observed_if <- plot_ss_if(
    tiny_ss_data,
    probs = c(0, 0.5, 1),
    to_plot = "data",
    use_lim = "data",
    factors = "E",
    within_plot = "S",
    layout = c(1, 2),
    legendpos = c(NA, NA)
  )
  expect_s3_class(observed_if, "data.frame")
  expect_named(observed_if, c("source", "dataset", "binning", "group_key",
                              "within_level", "x_plot", "bin_label", "ssd",
                              "value", "se", "n", "n_obs"))
  expect_true(all(observed_if$source == "data"))
  expect_true(all(observed_if$binning == "individual_quantile"))
  expect_true(all(observed_if$x_plot %in% c(0.5, 1)))
  expect_true(all(!is.na(observed_if$bin_label)))
  expect_true(all(observed_if$n_obs >= observed_if$n))
  expect_true(all(is.finite(observed_if$value)))

  observed_srrt <- plot_ss_srrt(
    tiny_ss_data,
    probs = c(0, 0.5, 1),
    use_global_quantiles = FALSE,
    to_plot = "data",
    use_lim = "data",
    factors = "E",
    within_plot = "S",
    layout = c(1, 2),
    legendpos = c(NA, NA)
  )
  expect_s3_class(observed_srrt, "data.frame")
  expect_true(all(observed_srrt$source == "data"))
  expect_true(all(observed_srrt$binning == "individual_quantile"))
  expect_true(all(observed_srrt$x_plot %in% c(0.5, 1)))
  expect_true(all(!is.na(observed_srrt$bin_label)))
  expect_true(all(is.finite(observed_srrt$value)))

  posterior_if <- plot_ss_if(
    tiny_ss_data,
    post_predict = tiny_post,
    probs = c(0, 0.5, 1),
    use_global_quantiles = TRUE,
    to_plot = "posterior",
    use_lim = "posterior",
    factors = "E",
    within_plot = "S",
    layout = c(1, 2),
    legendpos = c(NA, NA)
  )
  expect_true(all(posterior_if$source == "posterior"))
  expect_true(all(posterior_if$binning == "global_quantile"))
  expect_true(all(c("lower", "median", "upper", "n_draws") %in% names(posterior_if)))
  expect_true(all(is.finite(posterior_if$median)))
  expect_true(all(posterior_if$n_draws > 0))

  mixed_value <- NULL
  expect_message(
    mixed_value <- plot_ss_if(
      tiny_ss_data,
      post_predict = tiny_post,
      ssd_binning = "value",
      to_plot = c("data", "posterior"),
      use_lim = c("data", "posterior"),
      factors = "E",
      within_plot = "S",
      layout = c(1, 2),
      legendpos = c(NA, NA)
    ),
    "ignored when `ssd_binning = \"value\"`",
    fixed = TRUE
  )
  expect_setequal(unique(mixed_value$source), c("data", "posterior"))
  expect_true(all(mixed_value$binning == "value"))
  expect_true(all(!is.na(mixed_value$ssd)))
  expect_true(all(!is.na(mixed_value$bin_label)))

  reduced_srrt <- plot_ss_srrt(
    tiny_ss_data,
    post_predict = tiny_post,
    probs = seq(0, 1, length.out = 7),
    use_global_quantiles = TRUE,
    on_duplicate_quantiles = "reduce",
    to_plot = "posterior",
    use_lim = "posterior",
    factors = "E",
    within_plot = "S",
    layout = c(1, 2),
    legendpos = c(NA, NA)
  )
  expect_true(all(reduced_srrt$source == "posterior"))
  expect_true(all(reduced_srrt$binning == "global_quantile"))
  expect_true(length(unique(reduced_srrt$x_plot)) < length(seq(0, 1, length.out = 7)) - 1)
  expect_true(all(reduced_srrt$n_draws > 0))

  printed_if <- capture.output(plot_ss_if(
    tiny_ss_data,
    post_predict = tiny_post,
    probs = c(0, 0.5, 1),
    use_global_quantiles = TRUE,
    to_plot = c("data", "posterior"),
    use_lim = c("data", "posterior"),
    factors = "E",
    within_plot = "S",
    layout = c(1, 2),
    legendpos = c(NA, NA),
    print_plot_data = TRUE
  ))
  expect_true(any(grepl("Predictive plotted summaries:", printed_if, fixed = TRUE)))
  expect_true(any(grepl("Observed plotted summaries:", printed_if, fixed = TRUE)))
  expect_no_error(plot_ss_srrt(
    tiny_ss_data,
    post_predict = tiny_post,
    probs = seq(0, 1, length.out = 7),
    use_global_quantiles = TRUE,
    on_duplicate_quantiles = "reduce",
    to_plot = "posterior",
    use_lim = "posterior",
    factors = "E",
    within_plot = "S",
    layout = c(1, 2),
    legendpos = c(NA, NA)
  ))
  expect_no_error(plot_ss_srrt(
    tiny_ss_data,
    post_predict = tiny_post,
    probs = seq(0, 1, length.out = 7),
    use_global_quantiles = TRUE,
    on_duplicate_quantiles = "value",
    to_plot = "posterior",
    use_lim = "posterior",
    factors = "E",
    within_plot = "S",
    layout = c(1, 2),
    legendpos = c(NA, NA)
  ))
})
