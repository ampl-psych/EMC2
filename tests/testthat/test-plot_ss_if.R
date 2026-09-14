dat <- make_ss_forstmann_data(seed = 123)

test_that("plot_ss_if produces stable snapshot", {
  vdiffr::expect_doppelganger(
    "IF_plot_ss",
    plot_ss_if(
      dat,
      probs = seq(0, 1, length.out = 5),
      use_global_quantiles = TRUE,
      to_plot = c("data"),
      use_lim = c("data"),factors = "E", within_plot = "S"
      )
  )
})

