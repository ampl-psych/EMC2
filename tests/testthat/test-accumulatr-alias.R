make_accumulatr_alias_spec <- function(parameters) {
  AccumulatR::race_spec() |>
    AccumulatR::add_accumulator("A", "lognormal") |>
    AccumulatR::add_accumulator("B", "lognormal") |>
    AccumulatR::add_outcome("A", "A") |>
    AccumulatR::add_outcome("B", "B") |>
    AccumulatR::set_parameters(parameters) |>
    AccumulatR::finalize_model()
}

make_accumulatr_alias_data <- function() {
  data.frame(
    subjects = factor(c("1", "1", "1")),
    condition = factor(c("low", "high", "low")),
    cov1 = c(0.1, 0.2, 0.3),
    R = factor(c("A", "B", "A")),
    rt = c(0.3, 0.5, 0.4)
  )
}

test_that("AccumulatR alias bridge exposes public sampled names and compiles runtime designs", {
  spec <- make_accumulatr_alias_spec(list(
    v = c("A.m", "B.m"),
    s = c("A.s", "B.s"),
    onset = c("A.t0", "B.t0"),
    gate = c("A.q", "B.q")
  ))
  mod <- AccumulatR_model(spec)
  dat <- make_accumulatr_alias_data()

  des <- design(
    data = dat,
    model = mod,
    formula = list(v ~ condition, s ~ 1, onset ~ 1, gate ~ 1),
    report_p_vector = FALSE
  )

  expect_true(all(!grepl("^p[0-9]+$", names(mod()$p_types))))
  expect_true(all(!grepl("^p[0-9]", names(sampled_pars(des)))))

  dadm <- EMC2:::design_model(dat, des, model = des$model, compress = FALSE, verbose = FALSE, rt_check = FALSE)

  expect_setequal(names(attr(dadm, "designs")), c("p1", "p2", "q", "t0"))
  expect_false(any(grepl("^p[0-9]", attr(dadm, "sampled_p_names"))))
  expect_false(is.null(attr(dadm, "accumulatr_bridge")))
})

test_that("AccumulatR alias bridge maps back to public parameter names", {
  spec <- make_accumulatr_alias_spec(list(
    v = c("A.m", "B.m"),
    s = c("A.s", "B.s"),
    onset = c("A.t0", "B.t0"),
    gate = c("A.q", "B.q")
  ))
  mod <- AccumulatR_model(spec)
  dat <- make_accumulatr_alias_data()

  des <- design(
    data = dat,
    model = mod,
    formula = list(v ~ condition, s ~ 1, onset ~ 1, gate ~ 1),
    report_p_vector = FALSE
  )
  dadm <- EMC2:::design_model(dat, des, model = des$model, compress = FALSE, verbose = FALSE, rt_check = FALSE)

  extra <- setdiff(attr(dadm, "sampled_p_names"), c("v", "s", "onset", "gate"))
  p_vector <- c(v = 2, s = 0, onset = log(0.2), gate = qnorm(0.25))
  if (length(extra) > 0) {
    p_vector[extra[[1]]] <- -1
  }
  p_vector <- p_vector[attr(dadm, "sampled_p_names")]

  pars <- EMC2:::get_pars_matrix_oo(p_vector, dadm, des$model())
  expect_equal(colnames(pars), c("v", "s", "gate", "onset"))
  expect_equal(unique(pars[, "gate"]), 0.25)
  expect_equal(unique(pars[, "onset"]), 0.2)
  expect_equal(sort(unique(pars[, "v"])), c(1, 2))

  mapped <- mapped_pars(des, p_vector)
  expect_true(all(c("v", "s", "gate", "onset") %in% names(mapped)))
  expect_false(any(duplicated(names(mapped))))
  expect_false(any(grepl("^p[0-9]", names(mapped))))
})

test_that("AccumulatR alias bridge rejects mixed runtime-family mappings", {
  spec <- make_accumulatr_alias_spec(list(
    bad = c("A.m", "A.s")
  ))

  expect_error(
    AccumulatR_model(spec),
    "must map to exactly one runtime slot family"
  )
})

test_that("AccumulatR alias bridge enforces post-map trend safety", {
  spec <- make_accumulatr_alias_spec(list(
    v = "A.m",
    b = "B.m",
    s = c("A.s", "B.s"),
    onset = c("A.t0", "B.t0"),
    gate = c("A.q", "B.q")
  ))
  mod <- AccumulatR_model(spec)
  dat <- make_accumulatr_alias_data()

  trend_premap <- make_trend(
    par_names = "v",
    cov_names = list("cov1"),
    kernels = "exp_incr",
    phase = "premap",
    at = NULL
  )
  expect_no_error(
    design(
      data = dat,
      model = mod,
      formula = list(v ~ 1, b ~ 1, s ~ 1, onset ~ 1, gate ~ 1),
      trend = trend_premap,
      report_p_vector = FALSE
    )
  )

  trend_post <- make_trend(
    par_names = "v",
    cov_names = list("cov1"),
    kernels = "exp_incr",
    phase = "posttransform",
    at = NULL
  )
  expect_error(
    design(
      data = dat,
      model = mod,
      formula = list(v ~ 1, b ~ 1, s ~ 1, onset ~ 1, gate ~ 1),
      trend = trend_post,
      report_p_vector = FALSE
    ),
    "posttransform trend on 'v' is ambiguous"
  )
})
