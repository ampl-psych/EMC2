# register_nn_model(): the registration contract (R/nn_register.R).
#
# Complements test-nn-port.R (which checks the compiled evaluators against
# golden vectors and the R port); this file checks that register_nn_model()
# reads a card correctly, reproduces the evaluator's numbers through the
# ordinary dfun/pfun interface, refuses malformed cards, and that the
# registry (cache, save/load, sha256 pin) behaves.

nle <- function(f) getFromNamespace(f, "EMC2")

port <- new.env()
sys.source(test_path("port", "flow_race.R"), envir = port)
sys.source(test_path("port", "flow_ddm.R"), envir = port)

box_pars <- function(reg, n, lo_q = .1, hi_q = .9) {
  # random in-box rows on the natural scale, named reg$pars order
  m <- t(replicate(n, reg$lower + (reg$upper - reg$lower) * runif(length(reg$lower), lo_q, hi_q)))
  colnames(m) <- reg$pars
  m
}
hand_theta <- function(reg, pars) {
  # natural-scale pars -> sampled-scale theta, built from the card's
  # transforms directly (not via nn_context)
  theta <- matrix(0, nrow(pars), length(reg$context_names))
  for (j in seq_along(reg$context_names)) {
    tf <- reg$transforms[j]; p <- reg$context_names[j]
    theta[, j] <- switch(tf, identity = pars[, p], log = log(pars[, p]), probit = qnorm(pars[, p]))
  }
  theta
}

# --- a. shipped artefacts register by name and by path ----------------------
test_that("shipped artefacts register by name and by full path", {
  by_name <- register_nn_model("rdm_small")()$nn
  by_path <- register_nn_model(file.path(nle("nle_dir")(), "rdm_small.rds"))()$nn
  expect_identical(by_name$kind, "flow_race")
  expect_identical(by_name$pars, by_name$context_names)
  expect_identical(by_name$sha256, by_path$sha256)
  expect_identical(by_path$artefact, "rdm_small")   # shipped, however it was addressed

  codes <- c(identity = 0, log = 1, probit = 2)[by_name$transforms]
  expect_equal(unname(codes), by_name$transform_codes)

  man <- nle("nle_manifest")()
  expect_identical(by_name$sha256, man[["rdm_small.rds"]])

  expect_identical(DDMnn()$p_types, DDM()$p_types)
  expect_identical(RDMnn()$p_types, RDM()$p_types)
  expect_identical(DDMnn()$transform, DDM()$transform)
  expect_identical(RDMnn()$transform, RDM()$transform)
  expect_identical(DDMnn("ddm_st0zero")$p_types, DDM()$p_types)
  expect_identical(DDMnn("ddm_st0zero")$transform, DDM()$transform)

  d0 <- DDMnn("ddm_st0zero")
  expect_equal(d0$bound$exception, c(st0 = 0))
  expect_equal(unname(d0$bound$minmax[, "st0"]), c(0, 0))

  reg <- register_nn_model("ddm_cap256w_c4")()$nn
  bn <- reg$card$bounds_natural; nm <- names(reg$card$context_transforms)
  minmax <- DDMnn()$bound$minmax
  for (p in nm) expect_equal(unname(minmax[, p]), c(bn$lower[match(p, nm)], bn$upper[match(p, nm)]))
})

# --- b. registered dfun/pfun reproduce the compiled evaluator ---------------
test_that("DDMnn dfun/pfun reproduce ddm_ens_eval_trials_cpp on hand-built theta", {
  set.seed(101)
  for (name in c("ddm_cap256w_c4", "ddm_st0zero")) {
    reg <- register_nn_model(name)()$nn
    m <- register_nn_model(name)()
    ptr <- nle("nn_ptr")(reg)
    rows <- box_pars(reg, 5)                     # a handful of distinct rows
    idx <- sample(5, 50, TRUE)                   # interleaved, ~50 trials
    pars <- rows[idx, , drop = FALSE]
    rt <- exp(runif(50, log(.06), log(2)))
    R <- factor(sample(1:2, 50, TRUE))
    theta <- hand_theta(reg, pars)
    ev <- nle("ddm_ens_eval_trials_cpp")(ptr, theta, rt, as.integer(as.character(R)))
    expect_identical(m$dfun(rt, R, pars), ev$pdf)
    expect_identical(m$pfun(rt, R, pars), ev$cdf)
  }
})

test_that("RDMnn dfun/pfun reproduce flow_eval_trials_cpp on hand-built theta", {
  set.seed(102)
  reg <- register_nn_model("rdm_small")()$nn
  m <- register_nn_model("rdm_small")()
  ptr <- nle("nn_ptr")(reg)
  rows <- box_pars(reg, 5)
  idx <- sample(5, 50, TRUE)
  pars <- rows[idx, , drop = FALSE]
  rt <- exp(runif(50, log(.06), log(2)))
  theta <- hand_theta(reg, pars)
  ev <- nle("flow_eval_trials_cpp")(ptr, theta, rt)
  expect_identical(m$dfun(rt, pars), ev$pdf)
  expect_identical(m$pfun(rt, pars), ev$cdf)
})

test_that("DDMnn('ddm_st0zero') reproduces the golden joint density (1e-9 relative)", {
  golden <- read.csv(test_path("golden", "ddm_z0_s503_e282_golden.csv"))
  m <- register_nn_model("ddm_st0zero")()
  pars <- cbind(v = golden$v, a = exp(golden$a), t0 = exp(golden$t0), s = exp(golden$s),
                Z = pnorm(golden$Z), SZ = pnorm(golden$SZ), sv = exp(golden$sv), st0 = 0)
  d1 <- m$dfun(golden$rt, factor(golden$R), pars)
  ref <- exp(golden$log_pdf_joint)
  expect_lt(max(abs(d1 - ref) / pmax(abs(ref), 1e-300)), 1e-9)
})

# --- c. pipeline sum: calc_ll_manager against the model's own dfun/pfun -----
test_that("DDMnn: calc_ll_manager's sum equals a hand-summed dfun over dadm", {
  set.seed(103)
  des <- suppressMessages(design(
    factors = list(subjects = 1, S = c("a", "b")), Rlevels = c("a", "b"), model = DDMnn,
    formula = list(v ~ S, a ~ 1, t0 ~ 1),
    constants = c(s = log(1), sv = log(.5), SZ = qnorm(.3), st0 = log(.1), Z = qnorm(.5)),
    report_p_vector = FALSE))
  p_vector <- sampled_pars(des, doMap = FALSE)
  p_vector[] <- c(0.5, -.5, log(1.2), log(.25))[seq_along(p_vector)]
  dat <- make_data(p_vector, des, n_trials = 30)
  emc <- make_emc(dat, des, type = "single", n_chains = 1, compress = FALSE,
                  verbose = FALSE, rt_resolution = .001)
  dadm <- emc[[1]]$data[[1]]
  pm <- matrix(p_vector, nrow = 1, dimnames = list(NULL, names(p_vector)))
  ll1 <- calc_ll_manager(pm, dadm, des$model)
  m <- des$model()
  pars <- get_pars_matrix_oo(p_vector, dadm, des$model)
  ll2 <- sum(pmax(log(1e-10), log(m$dfun(dadm$rt, dadm$R, pars))))
  expect_equal(ll1, ll2)
})

test_that("RDMnn: calc_ll_manager's sum equals a hand-summed race likelihood over dadm", {
  set.seed(104)
  des <- suppressMessages(design(
    factors = list(subjects = 1, S = c("a", "b")), Rlevels = c("a", "b"), model = RDMnn,
    formula = list(v ~ S, B ~ 1, t0 ~ 1, A ~ 1), constants = c(s = log(1)),
    report_p_vector = FALSE))
  p_vector <- sampled_pars(des, doMap = FALSE)
  p_vector[] <- c(log(1.5), log(.8), log(1.1), log(.25), log(.3))[seq_along(p_vector)]
  dat <- make_data(p_vector, des, n_trials = 30)
  emc <- make_emc(dat, des, type = "single", n_chains = 1, compress = FALSE,
                  verbose = FALSE, rt_resolution = .001)
  dadm <- emc[[1]]$data[[1]]
  pm <- matrix(p_vector, nrow = 1, dimnames = list(NULL, names(p_vector)))
  ll1 <- calc_ll_manager(pm, dadm, des$model)
  m <- des$model()
  pars <- get_pars_matrix_oo(p_vector, dadm, des$model)
  lds <- numeric(nrow(dadm))
  lds[dadm$winner] <- log(m$dfun(dadm$rt[dadm$winner], pars[dadm$winner, , drop = FALSE]))
  lds[!dadm$winner] <- log(1 - m$pfun(dadm$rt[!dadm$winner], pars[!dadm$winner, , drop = FALSE]))
  ll <- (lds[dadm$winner] + lds[!dadm$winner])[attr(dadm, "expand")]
  ll2 <- sum(pmax(log(1e-10), ll))
  expect_equal(ll1, ll2)
})

# --- d. JSON cards ------------------------------------------------------------
test_that("a JSON RDM card round-trips the weights and agrees with the .rds registration", {
  skip_if_not_installed("jsonlite")
  x <- readRDS(file.path(nle("nle_dir")(), "rdm_small.rds"))
  tmp <- tempfile(fileext = ".json")
  writeLines(jsonlite::toJSON(x, auto_unbox = TRUE, digits = NA), tmp)
  raw <- jsonlite::fromJSON(tmp, simplifyDataFrame = FALSE)
  W <- as.matrix(raw$mlp$layers[[1]]$W)
  expect_equal(dim(W), dim(x$mlp$layers[[1]]$W))   # n_in x n_out, as in the .rds

  m_rds <- register_nn_model("rdm_small")()
  m_json <- register_nn_model(tmp, twin = RDM)()
  rt <- exp(runif(30, log(.1), log(2)))
  pars <- cbind(v = 1.5, B = .8, t0 = .25, s = 1, A = .3)[rep(1, 30), ]
  expect_lt(max(abs(m_rds$dfun(rt, pars) - m_json$dfun(rt, pars))), 1e-9)
  expect_lt(max(abs(m_rds$pfun(rt, pars) - m_json$pfun(rt, pars))), 1e-9)
})

test_that("a JSON DDM (ddm_st0zero) member card agrees with the .rds registration", {
  skip_if_not_installed("jsonlite")
  x <- readRDS(file.path(nle("nle_dir")(), "ddm_st0zero.rds"))
  tmp <- tempfile(fileext = ".json")
  writeLines(jsonlite::toJSON(x$members[[1]], auto_unbox = TRUE, digits = NA), tmp)

  m_rds <- register_nn_model("ddm_st0zero")()
  m_json <- register_nn_model(tmp, twin = DDM)()
  rt <- exp(runif(30, log(.1), log(2)))
  R <- factor(sample(1:2, 30, TRUE))
  pars <- cbind(v = 1, a = 1, t0 = .3, s = 1, Z = .5, SZ = .3, sv = .5, st0 = 0)[rep(1, 30), ]
  expect_lt(max(abs(m_rds$dfun(rt, R, pars) - m_json$dfun(rt, R, pars))), 1e-9)
  expect_lt(max(abs(m_rds$pfun(rt, R, pars) - m_json$pfun(rt, R, pars))), 1e-9)
})

# --- e. refusals --------------------------------------------------------------
test_that("a card whose bounds_sampled disagrees with transformed bounds_natural is refused", {
  x <- readRDS(file.path(nle("nle_dir")(), "rdm_small.rds"))
  ord <- match(c("v", "t0", "B", "s", "A"), x$context_names)   # B/t0 permuted
  x$context_names <- x$context_names[ord]
  tmp <- tempfile(fileext = ".rds"); saveRDS(x, tmp)
  expect_error(register_nn_model(tmp, twin = RDM)(), "permuted or mislabelled")
})

test_that("pars that would permute the card's inputs is refused", {
  expect_error(register_nn_model("rdm_small", pars = c("B", "v", "t0", "s", "A")), "would permute")
})

test_that("kind mismatch, disagreeing transforms are refused", {
  expect_error(register_nn_model("rdm_small", kind = "flow_joint"), "not flow_joint")
  expect_error(register_nn_model("rdm_small",
    transforms = c(v = "identity", B = "log", t0 = "log", s = "log", A = "log")), "disagrees")
})

test_that("st0zero without context_encoding needs exceptions, and registers correctly with it", {
  x <- readRDS(file.path(nle("nle_dir")(), "ddm_st0zero.rds"))
  x$members[[1]]$context_encoding <- NULL
  tmp <- tempfile(fileext = ".rds"); saveRDS(x, tmp)
  expect_error(register_nn_model(tmp, twin = DDM)(), "exceptions = c\\(st0")

  m <- register_nn_model(tmp, twin = DDM, exceptions = c(st0 = 0))()
  m_ref <- register_nn_model("ddm_st0zero")()
  rt <- exp(runif(10, log(.1), log(2))); R <- factor(sample(1:2, 10, TRUE))
  pars <- cbind(v = 1, a = 1, t0 = .3, s = 1, Z = .5, SZ = .3, sv = .5, st0 = 0)[rep(1, 10), ]
  expect_identical(m$dfun(rt, R, pars), m_ref$dfun(rt, R, pars))
})

test_that("exceptions cannot fix a network input", {
  expect_error(register_nn_model("ddm_cap256w_c4", twin = DDM, exceptions = c(sv = 0)),
               "is an input of the network")
})

test_that("a twin parameter the network would ignore is refused", {
  x <- readRDS(file.path(nle("nle_dir")(), "ddm_st0zero.rds"))
  x$members[[1]]$context_transforms$st0 <- NULL
  x$members[[1]]$bounds_natural$lower <- x$members[[1]]$bounds_natural$lower[1:7]
  x$members[[1]]$bounds_natural$upper <- x$members[[1]]$bounds_natural$upper[1:7]
  tmp <- tempfile(fileext = ".rds"); saveRDS(x, tmp)
  expect_error(register_nn_model(tmp, twin = DDM)(), "would ignore")
})

test_that("cdf = FALSE on a flow_race is refused; on a flow_joint pfun alone refuses", {
  expect_error(register_nn_model("rdm_small", cdf = FALSE), "survivor")
  m <- register_nn_model("ddm_st0zero", cdf = FALSE)()
  rt <- exp(runif(10, log(.1), log(2))); R <- factor(sample(1:2, 10, TRUE))
  pars <- cbind(v = 1, a = 1, t0 = .3, s = 1, Z = .5, SZ = .3, sv = .5, st0 = 0)[rep(1, 10), ]
  expect_true(all(is.finite(m$dfun(rt, R, pars))))
  expect_error(m$pfun(rt, R, pars), "without a CDF")
})

test_that("a card naming no model and given no twin is refused", {
  x <- readRDS(file.path(nle("nle_dir")(), "rdm_small.rds"))
  x$model <- NULL
  tmp <- tempfile(fileext = ".rds"); saveRDS(x, tmp)
  expect_error(register_nn_model(tmp), "no analytic twin")
})

# --- f. renaming is allowed ---------------------------------------------------
test_that("renaming a context name (in context_names and context_transforms) is allowed", {
  x <- readRDS(file.path(nle("nle_dir")(), "rdm_small.rds"))
  x$context_names[x$context_names == "v"] <- "drift"
  names(x$context_transforms)[names(x$context_transforms) == "v"] <- "drift"
  tmp <- tempfile(fileext = ".rds"); saveRDS(x, tmp)
  m <- register_nn_model(tmp, pars = c("v", "B", "t0", "s", "A"))()
  rt <- exp(runif(20, log(.1), log(2)))
  pars <- cbind(v = 1.5, B = .8, t0 = .25, s = 1, A = .3)[rep(1, 20), ]
  expect_identical(m$dfun(rt, pars), RDMnn()$dfun(rt, pars))
})

# --- g. pre: shifting rt before the network sees it --------------------------
test_that("pre (a parameter name, or a shift function) shifts the network's time", {
  reg <- register_nn_model("rdm_small")()$nn
  m_name <- register_nn_model("rdm_small", pre = "t0")()
  ptr <- nle("nn_ptr")(m_name$nn)
  t0 <- .25
  rt <- exp(runif(40, log(t0 + .01), log(2 + t0)))
  pars <- cbind(v = 1.5, B = .8, t0 = t0, s = 1, A = .3)[rep(1, 40), ]
  d1 <- m_name$dfun(rt, pars)
  theta <- nle("nn_context")(pars, m_name$nn)
  ref <- nle("flow_eval_trials_cpp")(ptr, theta, rt - t0)$pdf
  expect_identical(d1, ref)

  rt2 <- c(seq(.01, t0, length.out = 5), rt)      # at/below t0: exactly 0
  pars2 <- cbind(v = 1.5, B = .8, t0 = t0, s = 1, A = .3)[rep(1, length(rt2)), ]
  expect_true(all(m_name$dfun(rt2, pars2)[1:5] == 0))

  m_fun <- register_nn_model("rdm_small", pre = function(rt, pars) rt - pars[, "t0"])()
  expect_identical(m_fun$dfun(rt, pars), d1)

  expect_error(register_nn_model("rdm_small", pre = function(rt, pars) 2 * rt), "must shift")
})

test_that("nn_check_pre reports which parameter(s) a shift function consumes", {
  reg <- register_nn_model("rdm_small")()$nn
  pars <- c("v", "B", "s", "A")                   # t0 held out of the network inputs
  box <- list(lower = reg$lower[pars], upper = reg$upper[pars])
  used <- nle("nn_check_pre")(function(rt, pars) rt - pars[, "t0"], pars, reg$defaults, box, "x: ")
  expect_identical(used, "t0")
})

# --- h. regression_joint (synthetic card) ------------------------------------
mk_mlp <- function(dims, norm) {
  layers <- lapply(seq_len(length(dims) - 1), function(i)
    list(W = matrix(rnorm(dims[i] * dims[i + 1], sd = 1 / sqrt(dims[i])), dims[i]),
         b = rnorm(dims[i + 1], sd = .1)))
  mlp <- list(activation = "gelu_tanh", layers = layers, use_norm = norm)
  if (norm) mlp$norms <- lapply(dims[2:(length(dims) - 1)], function(d)
    list(scale = runif(d, .5, 1.5), bias = rnorm(d, sd = .1), eps = 1e-5))
  mlp
}
mk_regression_card <- function() {
  m1 <- readRDS(file.path(nle("nle_dir")(), "ddm_cap256w_c4.rds"))$members[[1]]
  list(context_names = m1$context_names, context_transforms = m1$context_transforms,
       bounds_natural = m1$bounds_natural, bounds_sampled = m1$bounds_sampled,
       model = "DDM", mlp = mk_mlp(c(10, 32, 16, 1), FALSE),
       input_layout = c("v", "a", "t0", "log_rt", "R", "s", "Z", "SZ", "sv", "st0"),
       scaler = list(mean = rnorm(10, sd = .1), scale = runif(10, .5, 2)),
       response_values = c(-1, 1), output_scaler = list(mean = -1, scale = 2))
}

test_that("a regression_joint card registers and its dfun matches a hand-written forward pass", {
  set.seed(105)
  card <- mk_regression_card()
  tmp <- tempfile(fileext = ".rds"); saveRDS(card, tmp)
  m <- register_nn_model(tmp)()
  expect_identical(m$nn$kind, "regression_joint")
  expect_false(m$nn$cdf)

  reg <- m$nn; ctx <- reg$context_names
  n <- 30
  pars <- box_pars(reg, n)
  rt <- exp(runif(n, log(.1), log(2)))
  R <- sample(1:2, n, TRUE)
  d1 <- m$dfun(rt, factor(R), pars)

  theta <- hand_theta(reg, pars)
  lay <- card$input_layout; rv <- card$response_values
  X <- matrix(0, n, length(lay))
  for (i in seq_along(lay))
    X[, i] <- switch(lay[i], rt = rt, log_rt = log(rt), R = rv[R], theta[, match(lay[i], ctx)])
  Xs <- sweep(sweep(X, 2L, card$scaler$mean), 2L, card$scaler$scale, "/")
  y <- port$.mlp_fwd(card$mlp, Xs, use_norm = FALSE)[, 1L]
  ref <- exp(y * card$output_scaler$scale + card$output_scaler$mean)
  expect_lt(max(abs(d1 - ref)), 1e-10)

  pars_oob <- pars[1, , drop = FALSE]; pars_oob[, "v"] <- reg$upper[["v"]] + 1
  expect_identical(m$dfun(rt[1], factor(R[1]), pars_oob), 0)
  expect_error(m$pfun(rt, factor(R), pars), "without a CDF")
})

test_that("regression_joint card checks: input_layout, single output, cdf", {
  card <- mk_regression_card()

  card_no_layout <- card; card_no_layout$input_layout <- NULL
  tmp <- tempfile(fileext = ".rds"); saveRDS(card_no_layout, tmp)
  expect_error(register_nn_model(tmp)(), "input_layout")

  card_short <- card; card_short$input_layout <- card_short$input_layout[1:9]
  tmp <- tempfile(fileext = ".rds"); saveRDS(card_short, tmp)
  expect_error(register_nn_model(tmp)(), "input_layout")

  card_2out <- card; card_2out$mlp <- mk_mlp(c(10, 32, 16, 2), FALSE)
  tmp <- tempfile(fileext = ".rds"); saveRDS(card_2out, tmp)
  expect_error(register_nn_model(tmp)(), "single output")

  tmp_ok <- tempfile(fileext = ".rds"); saveRDS(card, tmp_ok)
  expect_error(register_nn_model(tmp_ok, cdf = TRUE)(), "no CDF")
})

test_that("a regression_joint model runs through design() + make_emc() + calc_ll_manager", {
  set.seed(106)
  card <- mk_regression_card()
  tmp <- tempfile(fileext = ".rds"); saveRDS(card, tmp)
  regf <- register_nn_model(tmp)
  des <- suppressMessages(design(
    factors = list(subjects = 1, S = c("a", "b")), Rlevels = c("a", "b"), model = regf,
    formula = list(v ~ S, a ~ 1, t0 ~ 1),
    constants = c(s = log(1), sv = log(.5), SZ = qnorm(.3), st0 = log(.1), Z = qnorm(.5)),
    report_p_vector = FALSE))
  p_vector <- sampled_pars(des, doMap = FALSE)
  p_vector[] <- c(0.5, -.5, log(1.2), log(.25))[seq_along(p_vector)]
  dat <- make_data(p_vector, des, n_trials = 20)
  emc <- make_emc(dat, des, type = "single", n_chains = 1, compress = FALSE,
                  verbose = FALSE, rt_resolution = .001)
  dadm <- emc[[1]]$data[[1]]
  pm <- matrix(p_vector, nrow = 1, dimnames = list(NULL, names(p_vector)))
  ll <- calc_ll_manager(pm, dadm, des$model)
  expect_true(is.finite(ll))
})

# --- i. registry: lazy rebuild, save/load, sha256 pin, no weights ------------
test_that("the evaluator is rebuilt lazily after the cache is cleared, with identical values", {
  m <- register_nn_model("rdm_small")()
  rt <- exp(runif(10, log(.1), log(2)))
  pars <- cbind(v = 1.5, B = .8, t0 = .25, s = 1, A = .3)[rep(1, 10), ]
  d1 <- m$dfun(rt, pars)
  rm(list = ls(nle("nle_cache")), envir = nle("nle_cache"))
  expect_identical(m$dfun(rt, pars), d1)
})

test_that("saveRDS/readRDS of a registered model list round-trips and evaluates", {
  m <- register_nn_model("rdm_small")()
  rt <- exp(runif(10, log(.1), log(2)))
  pars <- cbind(v = 1.5, B = .8, t0 = .25, s = 1, A = .3)[rep(1, 10), ]
  d1 <- m$dfun(rt, pars)
  tmp <- tempfile(fileext = ".rds"); saveRDS(m, tmp)
  m2 <- readRDS(tmp)
  expect_identical(m2$dfun(rt, pars), d1)
})

test_that("a card modified after registration is refused (sha256 changed)", {
  usercard <- readRDS(file.path(nle("nle_dir")(), "rdm_small.rds"))
  tmp <- tempfile(fileext = ".rds"); saveRDS(usercard, tmp)
  m <- register_nn_model(tmp)()
  rm(list = ls(nle("nle_cache")), envir = nle("nle_cache"))
  usercard$ll_floor_log <- -30
  Sys.sleep(.05)
  saveRDS(usercard, tmp)
  rt <- exp(runif(10, log(.1), log(2)))
  pars <- cbind(v = 1.5, B = .8, t0 = .25, s = 1, A = .3)[rep(1, 10), ]
  expect_error(m$dfun(rt, pars), "changed since the model was registered")
})

test_that("a registration carries no weights and is small", {
  m <- DDMnn()
  expect_null(m$nn$card$mlp)
  expect_null(m$nn$card$flow_mlp)
  expect_lt(length(serialize(m$nn, NULL)), 2e4)
})

# --- j. design()-time refusals -------------------------------------------------
test_that("design()-time refusals: defaults trap, out-of-box constants, fixed parameters", {
  f <- function(...) suppressMessages(design(...))

  expect_error(f(factors = list(subjects = 1, S = c("a", "b")), Rlevels = c("a", "b"), model = DDMnn,
                 formula = list(v ~ S, a ~ 1, t0 ~ 1),
                 constants = c(s = log(1), SZ = qnorm(.3), st0 = log(.1), Z = qnorm(.5)),
                 report_p_vector = FALSE),
               "Parameter 'sv' is at the model default")

  expect_error(f(factors = list(subjects = 1, S = c("a", "b")), Rlevels = c("a", "b"), model = DDMnn,
                 formula = list(v ~ S, t0 ~ 1),
                 constants = c(s = log(1), sv = log(.5), SZ = qnorm(.3), st0 = log(.1),
                               Z = qnorm(.5), a = log(10)),
                 report_p_vector = FALSE),
               "outside")

  ddm_st0 <- function() DDMnn("ddm_st0zero")
  expect_error(f(factors = list(subjects = 1, S = c("a", "b")), Rlevels = c("a", "b"), model = ddm_st0,
                 formula = list(v ~ S, a ~ 1, t0 ~ 1, st0 ~ 1),
                 constants = c(s = log(1), sv = log(.5), SZ = qnorm(.3), Z = qnorm(.5)),
                 report_p_vector = FALSE),
               "fixed at 0")

  expect_error(f(factors = list(subjects = 1, S = c("a", "b")), Rlevels = c("a", "b"), model = ddm_st0,
                 formula = list(v ~ S, a ~ 1, t0 ~ 1),
                 constants = c(s = log(1), sv = log(.5), SZ = qnorm(.3), Z = qnorm(.5), st0 = log(.1)),
                 report_p_vector = FALSE),
               "fixed at 0")

  expect_error(f(factors = list(subjects = 1, S = c("a", "b", "c")), Rlevels = c("a", "b", "c"),
                 model = DDMnn, formula = list(v ~ S, a ~ 1, t0 ~ 1),
                 constants = c(s = log(1), sv = log(.5), SZ = qnorm(.3), st0 = log(.1), Z = qnorm(.5)),
                 report_p_vector = FALSE),
               "two-response")

  expect_error(f(factors = list(subjects = 1, S = c("a", "b")), Rlevels = c("a", "b"), model = RDMnn,
                 formula = list(v ~ S, B ~ 1, t0 ~ 1), constants = c(s = log(1)),
                 report_p_vector = FALSE),
               "'A' is at the model default")

  des_ok <- f(factors = list(subjects = 1, S = c("a", "b")), Rlevels = c("a", "b"), model = ddm_st0,
              formula = list(v ~ S, a ~ 1, t0 ~ 1),
              constants = c(s = log(1), sv = log(.5), SZ = qnorm(.3), Z = qnorm(.5)),
              report_p_vector = FALSE)
  expect_s3_class(des_ok, "emc.design")
})

# --- k. the simulator is the twin's -------------------------------------------
test_that("make_data from a DDMnn design equals make_data from the twin DDM design", {
  formula <- list(v ~ S, a ~ 1, t0 ~ 1)
  constants <- c(s = log(1), sv = log(.5), SZ = qnorm(.3), st0 = log(.1), Z = qnorm(.5))
  des1 <- suppressMessages(design(factors = list(subjects = 1:2, S = c("a", "b")),
    Rlevels = c("a", "b"), model = DDMnn, formula = formula, constants = constants,
    report_p_vector = FALSE))
  des2 <- suppressMessages(design(factors = list(subjects = 1:2, S = c("a", "b")),
    Rlevels = c("a", "b"), model = DDM, formula = formula, constants = constants,
    report_p_vector = FALSE))
  p <- sampled_pars(des1, doMap = FALSE)
  p[] <- c(.5, -.5, log(1.2), log(.25))[seq_along(p)]
  set.seed(201); d1 <- make_data(p, des1, n_trials = 20)
  set.seed(201); d2 <- make_data(p, des2, n_trials = 20)
  expect_equal(d1, d2)
})

test_that("make_data from an RDMnn design equals make_data from the twin RDM design", {
  formula <- list(v ~ S, B ~ 1, t0 ~ 1, A ~ 1)
  constants <- c(s = log(1))
  des1 <- suppressMessages(design(factors = list(subjects = 1:2, S = c("a", "b")),
    Rlevels = c("a", "b"), model = RDMnn, formula = formula, constants = constants,
    report_p_vector = FALSE))
  des2 <- suppressMessages(design(factors = list(subjects = 1:2, S = c("a", "b")),
    Rlevels = c("a", "b"), model = RDM, formula = formula, constants = constants,
    report_p_vector = FALSE))
  p <- sampled_pars(des1, doMap = FALSE)
  p[] <- c(log(1.5), log(.8), log(1.1), log(.25), log(.3))[seq_along(p)]
  set.seed(202); d1 <- make_data(p, des1, n_trials = 20)
  set.seed(202); d2 <- make_data(p, des2, n_trials = 20)
  expect_equal(d1, d2)
})
