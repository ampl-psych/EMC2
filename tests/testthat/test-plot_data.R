# Per-panel `main` support in plot_cdf / plot_delta / plot_caf.

test_that("panel_main resolves the title for each panel", {
  expect_equal(panel_main(NULL, "gk", 1), "gk")          # default is the group key
  expect_equal(panel_main("Pre ", "gk", 1), "Pre gk")    # scalar is prefixed
  expect_equal(panel_main("", "gk", 1), "")              # "" blanks the title
  expect_equal(panel_main("Pre", "All Data", 1), "Pre")  # "All Data" drops the key
  expect_equal(panel_main(c("X", "Y", "Z"), "gk", 2), "Y")  # vector used verbatim, by panel index
})

test_that("check_panel_main validates the vector length", {
  expect_null(check_panel_main(NULL, 3))
  expect_silent(check_panel_main("one", 3))            # scalar always allowed
  expect_silent(check_panel_main(c("a", "b", "c"), 3)) # one per panel
  expect_error(check_panel_main(c("a", "b"), 3), "one per panel")
  expect_error(check_panel_main(1:3, 3), "character")
})

test_that("plot_cdf/plot_delta/plot_caf accept a per-panel main vector", {
  skip_on_os("windows")
  tmp <- tempfile(fileext = ".pdf")
  pdf(tmp)
  on.exit({ dev.off(); unlink(tmp) }, add = TRUE)

  # factors = "S" produces 2 panels for each of these functions.
  expect_error(plot_cdf(forstmann, factors = "S"), NA)
  expect_error(plot_cdf(forstmann, factors = "S", main = c("Left", "Right")), NA)
  expect_error(plot_delta(forstmann, factors = "S", main = c("Left", "Right")), NA)
  expect_error(plot_caf(forstmann, factors = "S", caf_factor = "S",
                        main = c("Left", "Right")), NA)

  # Wrong-length vectors are rejected with the panel-count message.
  expect_error(plot_cdf(forstmann, factors = "S", main = c("a", "b", "c")),
               "one per panel")
  expect_error(plot_delta(forstmann, factors = "S", main = c("a", "b", "c")),
               "one per panel")
  expect_error(plot_caf(forstmann, factors = "S", caf_factor = "S",
                        main = c("a", "b", "c")), "one per panel")
})
