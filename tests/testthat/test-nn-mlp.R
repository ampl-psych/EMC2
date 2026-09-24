# Plain-MLP likelihood networks (kinds "regression_joint" and "mlp_joint",
# src/mlp_lik.cpp): the LAN converted from ONNX by inst/scripts/onnx_to_card.py
# (HSSM's ddm_uniform_st, huggingface.co/Eitanm/ddm-st-lans) against
# onnxruntime, the network's meaning against EMC2's own DDM density, the
# out-of-box floor, and the C++ forward pass against plain R matrix algebra.
# Golden files: tests/testthat/golden/lan_ddm_uniform_st{.rds,_ort.csv}
# (card = the converter's output; csv = onnxruntime, float32, on 300 random rows).

nle <- function(f) getFromNamespace(f, "EMC2")
lan_card_path <- test_path("golden", "lan_ddm_uniform_st.rds")
lan_p_types <- c(v = 0, a = 1, z = .5, t = .5, st = .05)
lan_model <- function(...) register_nn_model(lan_card_path, p_types = lan_p_types, ...)

# A random tanh / GELU network, and its forward pass in plain R
mk_net <- function(dims, act) {
  layers <- lapply(seq_len(length(dims) - 1), function(i)
    list(W = matrix(rnorm(dims[i] * dims[i + 1], sd = 1 / sqrt(dims[i])), dims[i]),
         b = rnorm(dims[i + 1], sd = .1)))
  list(activation = act, layers = layers, use_norm = FALSE)
}
r_forward <- function(mlp, X) {
  H <- X
  for (i in seq_along(mlp$layers)) {
    H <- sweep(H %*% mlp$layers[[i]]$W, 2L, mlp$layers[[i]]$b, "+")
    if (i < length(mlp$layers))
      H <- if (mlp$activation == "tanh") tanh(H) else
        0.5 * H * (1 + tanh(sqrt(2 / pi) * (H + 0.044715 * H^3)))
  }
  H[, 1]
}

# --- a. the converted LAN against onnxruntime ------------------------------------
test_that("the converted LAN reproduces onnxruntime (float32) on 300 random rows", {
  g <- read.csv(test_path("golden", "lan_ddm_uniform_st_ort.csv"))
  m <- lan_model()()
  expect_identical(m$nn$kind, "mlp_joint")
  expect_identical(m$c_name, "NN")
  pars <- as.matrix(g[, c("v", "a", "z", "t", "st")])
  R <- ifelse(g$choice > 0, 2L, 1L)        # EMC2: response 1 = lower (-1), 2 = upper (+1)
  lp <- log(m$dfun(g$rt, R, pars))
  expect_lt(max(abs(lp - g$logp_ort)), 1e-4)      # observed 5e-6: float32 vs float64
  # the same numbers straight from the evaluator, without the exp/log round trip
  ptr <- nle("nn_ptr")(m$nn)
  expect_lt(max(abs(nle("mlp_lik_eval_cpp")(ptr, pars, g$rt, R) - g$logp_ort)), 1e-4)
})

test_that("the card-level facts of the LAN are read from the card, not assumed", {
  m <- lan_model()()
  expect_identical(m$nn$context_names, c("v", "a", "z", "t", "st"))
  expect_equal(m$nn$ll_floor, -16.11809565095832)
  expect_equal(m$nn$oob, m$nn$ll_floor)
  expect_equal(unname(m$nn$lower), c(-3, .3, .3, .25, .001))
  expect_equal(unname(m$nn$upper), c(3, 2.5, .7, 2.25, .25))
  expect_false(m$nn$cdf)
  expect_error(m$pfun(1, 1L, cbind(v = 0, a = 1, z = .5, t = .5, st = .05)), "without a CDF")
  # no analytic twin: the defaults must be supplied
  expect_error(register_nn_model(lan_card_path), "p_types")
})

# --- b. what the network means: EMC2's DDM density ----------------------------------
# The LAN's parameterisation (HSSM): a = half the boundary separation, z relative
# start, t = centre of the uniform non-decision time, st = its half-width;
# input choice +1 = upper. In EMC2's DDM: a = 2 a_lan, Z = z, t0 = t - st,
# st0 = 2 st, sv = SZ = 0, s = 1; R = 2 = upper. A LAN is an approximation, so
# the agreement is that of a trained network (median 0.04 log units where the
# density is not tiny), and the point of the test is that every wrong mapping
# is far worse.
test_that("the LAN's parameterisation maps onto EMC2's DDM as documented", {
  set.seed(31)
  n <- 3000
  th <- cbind(v = runif(n, -2, 2), a = runif(n, .5, 2), z = runif(n, .4, .6),
              t = runif(n, .35, 1), st = runif(n, .02, .2))
  rt <- th[, "t"] + th[, "st"] + rexp(n, 1.2) + .05          # after the non-decision time
  R <- sample(1:2, n, TRUE)
  lan <- log(lan_model()()$dfun(rt, R, th))
  d <- DDM()
  ddm_ll <- function(a = 2 * th[, "a"], t0 = th[, "t"] - th[, "st"], R2 = R, sv = 0)
    log(d$dfun(rt, R2, d$Ttransform(cbind(v = th[, "v"], a = a, sv = sv, t0 = t0,
                                          st0 = 2 * th[, "st"], s = 1, Z = th[, "z"], SZ = 0), NULL)))
  ld <- ddm_ll()
  dense <- is.finite(ld) & ld > -5
  expect_gt(sum(dense), 1000)
  err <- abs(lan - ld)[dense]
  expect_lt(median(err), .1)
  expect_lt(unname(quantile(err, .95)), .35)
  for (bad in list(ddm_ll(a = th[, "a"]),                 # a not doubled
                   ddm_ll(t0 = th[, "t"]),                # t taken as the lower end
                   ddm_ll(R2 = 3L - R))) {                # response coding swapped
    e2 <- abs(lan - bad)[dense]
    expect_gt(median(e2, na.rm = TRUE), 3 * median(err))
  }
})

# --- c. out-of-box rows, times at or below zero ---------------------------------------
test_that("out-of-box rows return the card's log floor; time <= 0 has zero density", {
  m <- lan_model()()
  th <- cbind(v = c(0, 4, 0, 0, 0), a = c(1, 1, 3, 1, 1), z = .5, t = .5, st = c(.05, .05, .05, .05, .3))
  R <- rep(1L, 5)
  rt <- rep(1, 5)
  ptr <- nle("nn_ptr")(m$nn)
  lp <- nle("mlp_lik_eval_cpp")(ptr, th, rt, R)
  expect_true(is.finite(lp[1]))
  expect_identical(lp[4], lp[1])                                  # same row, in the box
  expect_identical(lp[c(2, 3, 5)], rep(m$nn$ll_floor, 3))         # v, a, st outside
  expect_equal(m$dfun(c(1, 0, -1), c(1L, 1L, 2L), th[c(1, 1, 1), ]), c(exp(lp[1]), 0, 0),
               tolerance = 1e-13)
  th[1, "z"] <- NA
  expect_identical(nle("mlp_lik_eval_cpp")(ptr, th[1, , drop = FALSE], 1, 1L), m$nn$ll_floor)
})

# --- d. the evaluator against plain R matrix algebra -------------------------------------
test_that("mlp_lik matches an R forward pass: tanh and GELU, layouts, scalers, blocks", {
  set.seed(4)
  ctx <- c("p1", "p2", "p3")
  for (act in c("tanh", "gelu_tanh")) {
    lay <- c("p2", "log_rt", "p1", "R", "p3", "rt")
    card <- list(kind = "regression_joint", context_names = ctx,
                 context_transforms = list(p1 = "identity", p2 = "log", p3 = "probit"),
                 bounds_natural = list(lower = c(-1, .1, .05), upper = c(1, 3, .95)),
                 model = NULL, input_layout = lay, response_values = c(-1, 1),
                 scaler = list(mean = rnorm(6, sd = .2), scale = runif(6, .5, 2)),
                 output_scaler = list(mean = -3, scale = 2.5),
                 mlp = mk_net(c(6, 40, 24, 1), act))
    tmp <- tempfile(fileext = ".rds"); saveRDS(card, tmp)
    m <- register_nn_model(tmp, p_types = c(p1 = 0, p2 = 0, p3 = 0))()
    n <- 700                                                  # spans three 256-row blocks
    th <- cbind(p1 = runif(n, -1, 1), p2 = runif(n, .1, 3), p3 = runif(n, .05, .95))
    rt <- exp(runif(n, log(.1), log(3))); R <- sample(1:2, n, TRUE)
    X <- cbind(p2 = log(th[, "p2"]), log_rt = log(rt), p1 = th[, "p1"],
               R = c(-1, 1)[R], p3 = qnorm(th[, "p3"]), rt = rt)
    X <- sweep(sweep(X, 2L, card$scaler$mean), 2L, card$scaler$scale, "/")
    ref <- r_forward(card$mlp, X) * 2.5 - 3
    ptr <- nle("nn_ptr")(m$nn)
    got <- nle("mlp_lik_eval_cpp")(ptr, cbind(th[, 1], log(th[, 2]), qnorm(th[, 3])), rt, R)
    expect_lt(max(abs(got - ref)), 1e-12)
    expect_lt(max(abs(log(m$dfun(rt, R, th)) - ref)), 1e-12)
    expect_identical(nle("mlp_lik_valid")(ptr), TRUE)
  }
})

# --- e. card checks ---------------------------------------------------------------------------
test_that("mlp_joint card checks: floor, layout, activation, single output, cdf", {
  card <- readRDS(lan_card_path)
  reg <- function(x, ...) {
    tmp <- tempfile(fileext = ".rds"); saveRDS(x, tmp)
    register_nn_model(tmp, p_types = lan_p_types, ...)
  }
  expect_silent(reg(card))
  c1 <- card; c1$ll_floor_log <- NULL
  expect_error(reg(c1), "ll_floor_log")
  c1 <- card; c1$input_layout <- NULL
  expect_error(reg(c1), "input_layout")
  c1 <- card; c1$mlp$activation <- "relu"
  expect_error(reg(c1), "unsupported activation")
  c1 <- card; c1$input_layout[2] <- "zz"                             # an input nobody supplies
  expect_error(reg(c1), "neither context_names")
  c1 <- card; c1$input_layout <- card$input_layout[c(1, 1, 3:7)]      # an input listed twice
  expect_error(reg(c1), "exactly once")
  expect_error(reg(card, cdf = TRUE), "no CDF")
  expect_error(reg(card, kind = "flow_joint"), "not flow_joint")
  # a kind declared in the card must agree with the network it holds
  c1 <- card; c1$kind <- "flow_race"
  expect_error(reg(c1), "declares kind")
})

# --- f. JSON cards (the converter's output format) ---------------------------------------------------
test_that("the converter's JSON card registers the same network as the .rds", {
  skip_if_not_installed("jsonlite")
  card <- readRDS(lan_card_path)
  js <- tempfile(fileext = ".json")
  writeLines(jsonlite::toJSON(card, auto_unbox = TRUE, digits = NA), js)
  m <- register_nn_model(js, p_types = lan_p_types)()
  g <- read.csv(test_path("golden", "lan_ddm_uniform_st_ort.csv"))
  R <- ifelse(g$choice > 0, 2L, 1L)
  pars <- as.matrix(g[, c("v", "a", "z", "t", "st")])
  expect_lt(max(abs(log(m$dfun(g$rt, R, pars)) - g$logp_ort)), 1e-4)
})

# The documented mapping is where the LAN and the DDM density agree best: for
# each mapped quantity, scan it around its documented value; the best value
# must sit within a small band of it (a LAN is an approximation, so the minimum
# is not exactly at the documented value: measured minima are a = 2.05, t0 shift
# +0.01, st0 factor 2.1; the bands below are about twice that) and the
# documented value must beat values one band further out on both sides.
test_that("the mapping's parameters are recovered as the best agreement (profile scan)", {
  set.seed(31)
  n <- 4000
  th <- cbind(v = runif(n, -2, 2), a = runif(n, .5, 2), z = runif(n, .4, .6),
              t = runif(n, .35, 1), st = runif(n, .02, .2))
  rt <- th[, "t"] + th[, "st"] + rexp(n, 1.2) + .05
  R <- sample(1:2, n, TRUE)
  lan <- log(lan_model()()$dfun(rt, R, th))
  d <- DDM()
  ll <- function(a_scale = 2, dt = 0, st_scale = 2, dz = 0, v_scale = 1)
    log(d$dfun(rt, R, d$Ttransform(cbind(
      v = v_scale * th[, "v"], a = a_scale * th[, "a"], sv = 0, t0 = th[, "t"] - th[, "st"] + dt,
      st0 = st_scale * th[, "st"], s = 1, Z = th[, "z"] + dz, SZ = 0), NULL)))
  dense <- is.finite(ll()) & ll() > -5
  score <- function(x) median(abs(lan - x)[dense & is.finite(x)])
  scan <- function(arg, doc, band) {
    grid <- doc + band * c(-2, -1, -.5, 0, .5, 1, 2)
    sc <- vapply(grid, function(g) score(do.call(ll, stats::setNames(list(g), arg))), 0)
    best <- grid[which.min(sc)]
    expect_lte(abs(best - doc), band, label = paste0("best ", arg, " (", format(best), ")"))
    expect_lt(sc[4], sc[1]); expect_lt(sc[4], sc[7])          # beats the values two bands away on both sides
  }
  scan("a_scale", 2, .1)
  scan("dt", 0, .02)
  scan("st_scale", 2, .3)
  scan("dz", 0, .01)
  scan("v_scale", 1, .05)
})

# --- g. the converter itself (inst/scripts/onnx_to_card.py) ---------------------------------------
# Needs python with onnx and onnxruntime; point EMC2_NLE_PYTHON at it (default python3).
test_that("onnx_to_card.py turns the original ONNX file into the golden card, --verify passes", {
  skip_if_not_installed("jsonlite")
  py <- Sys.getenv("EMC2_NLE_PYTHON", "python3")
  has <- suppressWarnings(system2(py, c("-c", shQuote("import onnx, onnxruntime, numpy")),
                                  stdout = FALSE, stderr = FALSE)) == 0
  skip_if_not(has, "python with onnx and onnxruntime not available (set EMC2_NLE_PYTHON)")
  script <- system.file("scripts", "onnx_to_card.py", package = "EMC2")
  if (!nzchar(script)) script <- test_path("..", "..", "inst", "scripts", "onnx_to_card.py")
  skip_if_not(file.exists(script))
  out <- tempfile(fileext = ".json")
  msg <- system2(py, c(shQuote(script), shQuote(test_path("golden", "ddm_uniform_st.onnx")), shQuote(out),
                       "--preset", "ddm_uniform_st", "--verify", "500"), stdout = TRUE, stderr = TRUE)
  expect_null(attr(msg, "status"))
  diff <- as.numeric(sub(".*max \\|diff\\| = ([0-9.e+-]+).*", "\\1", grep("--verify", msg, value = TRUE)))
  expect_lt(diff, 1e-4)
  new <- jsonlite::fromJSON(out, simplifyDataFrame = FALSE)
  old <- readRDS(lan_card_path)
  for (f in c("kind", "context_names", "context_transforms", "bounds_natural", "input_layout",
              "response_values", "ll_floor_log", "mlp", "source_sha256"))
    expect_identical(new[[f]], old[[f]], label = f)
  # refuses what it cannot convert or was not told
  bad <- system2(py, c(shQuote(script), shQuote(test_path("golden", "ddm_uniform_st.onnx")), shQuote(out)),
                 stdout = TRUE, stderr = TRUE)
  expect_true(!is.null(attr(bad, "status")) && any(grepl("--params is required", bad)))
})
