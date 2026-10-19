# nn_golden_check(): the shipped artefacts against their golden vectors, and a
# perturbed golden row fails loudly.

golden_dir <- system.file("extdata", "flownn", "golden", package = "EMC2")

test_that("the shipped artefacts pass their golden checks", {
  r <- nn_golden_check("rdm_small")
  expect_true(r$pass); expect_identical(r$kind, "flow_race"); expect_equal(r$n, 575)
  expect_output(print(r), "^PASS golden check of 'rdm_small'")
  z <- nn_golden_check("ddm_st0zero")
  expect_true(z$pass); expect_lt(max(z$max_abs), 1e-9)
  # the cap256w golden was written with an earlier classifier: the flow part
  expect_true(nn_golden_check("ddm_cap256w_c4", flow_only = TRUE)$pass)
  expect_false(nn_golden_check("ddm_cap256w_c4", error = FALSE)$pass)
})

test_that("a perturbed golden row fails loudly, with the output and size named", {
  g <- read.csv(file.path(golden_dir, "rdm_small_golden.csv"))
  g$log_pdf[42] <- g$log_pdf[42] + 1e-6
  f <- tempfile(fileext = ".csv"); write.csv(g, f, row.names = FALSE)
  expect_error(nn_golden_check("rdm_small", golden = f), "FAIL golden check.*max \\|diff\\| 1.00e-06 \\(log_pdf\\)")
  r <- nn_golden_check("rdm_small", golden = f, error = FALSE)
  expect_false(r$pass)
  expect_equal(unname(r$max_abs["log_pdf"]), 1e-6, tolerance = 1e-6)
  # a golden file for another network fails by a wide margin
  expect_error(nn_golden_check("ddm_st0zero", golden = file.path(golden_dir, "ddm_cap256w_c4_golden.csv")),
               "lacks the input column|FAIL")
})

test_that("a card file is checked against the golden file beside it; a missing one is an error", {
  card <- readRDS(system.file("extdata", "flownn", "rdm_small.rds", package = "EMC2"))
  d <- tempfile("golden"); dir.create(d)
  saveRDS(card, file.path(d, "my_flow_weights.rds"))
  expect_error(nn_golden_check(file.path(d, "my_flow_weights.rds")), "No golden file")
  file.copy(file.path(golden_dir, "rdm_small_golden.csv"), file.path(d, "my_flow_golden.csv"))
  expect_true(nn_golden_check(file.path(d, "my_flow_weights.rds"))$pass)
})

test_that("a LAN card is checked against onnxruntime's values at float32 precision", {
  lan <- test_path("golden", "lan_ddm_uniform_st.rds")
  skip_if_not(file.exists(lan))
  lan_twin <- function() list(
    type = "DDM", p_types = c(v = 0, a = log(.75), z = 0, t = log(.55), st = log(.075)),
    transform = list(func = c(v = "identity", a = "exp", z = "pnorm", t = "exp", st = "exp")))
  r <- nn_golden_check(lan, golden = test_path("golden", "lan_ddm_uniform_st_ort.csv"), twin = lan_twin, tol = 1e-4)
  expect_true(r$pass)
  expect_error(nn_golden_check(lan, golden = test_path("golden", "lan_ddm_uniform_st_ort.csv"),
                               twin = lan_twin, tol = 1e-12), "FAIL")
})
