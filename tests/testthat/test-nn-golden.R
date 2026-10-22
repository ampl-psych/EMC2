# nn_golden_check(): the shipped artefacts against their golden vectors, and a
# perturbed golden row fails loudly.

golden_dir <- system.file("extdata", "flownn", "golden", package = "EMC2")
rdm_card <- "rdm_ss_base_s801_e271_weights.json"

test_that("the shipped artefacts pass their golden checks", {
  r <- nn_golden_check(rdm_card)
  expect_true(r$pass); expect_identical(r$kind, "flow_race"); expect_equal(r$n, 575)
  expect_output(print(r), paste0("^PASS golden check of '", rdm_card, "'"))
  z <- nn_golden_check("ddm_sc_s1106_e245_weights.json")
  expect_true(z$pass); expect_lt(max(z$max_abs), 1e-9)
  # RTCON2 (flow_cat): the model-level rows (mix, convolution, censoring) and the
  # card's embedded network-level vectors; 1e-8 as nle_verify.R (the
  # upper-censored masses 1 - F lose digits to cancellation)
  k <- nn_golden_check(RTCONnn, tol = 1e-8)
  expect_true(k$pass); expect_identical(k$kind, "flow_cat"); expect_equal(k$n, 528)
  expect_setequal(names(k$max_abs), c("log_pdf", "cdf", "net_log_pdf", "net_cdf"))
})

test_that("a perturbed golden row fails loudly, with the output and size named", {
  g <- read.csv(file.path(golden_dir, "rdm_ss_base_s801_e271_golden.csv"))
  g$log_pdf[42] <- g$log_pdf[42] + 1e-6
  f <- tempfile(fileext = ".csv"); write.csv(g, f, row.names = FALSE)
  expect_error(nn_golden_check(rdm_card, golden = f), "FAIL golden check.*max \\|diff\\| 1.00e-06 \\(log_pdf\\)")
  r <- nn_golden_check(rdm_card, golden = f, error = FALSE)
  expect_false(r$pass)
  expect_equal(unname(r$max_abs["log_pdf"]), 1e-6, tolerance = 1e-6)
})

test_that("a card file is checked against the golden file beside it; a missing one is an error", {
  src <- system.file("extdata", "flownn", rdm_card, package = "EMC2")
  d <- tempfile("golden"); dir.create(d)
  file.copy(src, file.path(d, "my_flow_weights.json"))
  expect_error(nn_golden_check(file.path(d, "my_flow_weights.json")), "No golden file")
  file.copy(file.path(golden_dir, "rdm_ss_base_s801_e271_golden.csv"),
            file.path(d, "my_flow_golden.csv"))
  expect_true(nn_golden_check(file.path(d, "my_flow_weights.json"))$pass)
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
