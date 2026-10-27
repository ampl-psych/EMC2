for (name in names(circular_parameters)) {
  test_that(paste(name, "fixed-input likelihood agrees with the R density"), {
    dat <- data.frame(subjects = factor(1), rt = c(.25, .4, .8, 1.6),
                      R = c(.2, .8, 1.6, 2.8))
    if (name == "CDM") dat$R <- dat$R - pi
    if (name %in% c("SDM", "HSDM", "PHSDM")) dat$R2 <- c(.3, .9, 1.7, 2.8)
    if (name == "SDM") dat$R2 <- dat$R2 - pi
    if (name == "HSDM") dat$R3 <- c(-2, -.3, .6, 2.4)

    emc <- make_emc(dat, circular_design(name), type = "single",
                    compress = FALSE, n_chains = 1, verbose = FALSE)
    parameters <- likelihood_parameters(circular_parameters[[name]])
    parameters[1:8, "sv"] <- -Inf
    parameters <- rbind(parameters, invalid = circular_parameters[[name]])
    parameters[nrow(parameters), "t0"] <- log(10)

    model_r <- emc[[1]]$model()
    model_r$c_name <- NULL
    expected <- calc_ll_manager(parameters, emc[[1]]$data[[1]], function() model_r)
    actual <- model_likelihood(emc, parameters)
    expect_lt(max(abs(as.numeric(actual) - expected)), 1e-7)
  })
}
