# Tests for recover_sbc() and the shared SBC assembly helpers.
# These fabricate the intermediate rep_<i>.rds / prior_samples.rds files that
# run_sbc() writes, so no model fitting is needed: the assembly logic is
# exercised directly and deterministically.

# --- helpers to fabricate intermediate files --------------------------------

make_single_rep <- function(par_names, seed) {
  set.seed(seed)
  p <- length(par_names)
  list(
    rank     = setNames(runif(p), par_names),
    med      = setNames(rnorm(p), par_names),
    bias     = setNames(rnorm(p), par_names),
    coverage = setNames(runif(p) > 0.1, par_names)
  )
}

make_hier_rep <- function(par_names, var_names, seed) {
  set.seed(seed)
  list(
    rank_mu_row   = setNames(runif(length(par_names)), NULL),
    rank_var_row  = setNames(runif(length(var_names)), NULL),
    var_col_names = var_names,
    rand_effects  = matrix(rnorm(length(par_names) * 2), nrow = 2,
                           dimnames = list(NULL, par_names))
  )
}

# Writes rep_<i>.rds (+ prior_samples.rds) into a temp dir and returns that dir,
# which is what the new recover_sbc() takes as its first argument (`tempdir`).
write_single_run <- function(dir, par_names, reps_idx, replicates,
                             write_prior = TRUE) {
  temp_dir <- file.path(dir, "run_temp")
  dir.create(temp_dir, showWarnings = FALSE)
  for (i in reps_idx)
    saveRDS(make_single_rep(par_names, i), file.path(temp_dir, paste0("rep_", i, ".rds")))
  if (write_prior) {
    prior_alpha <- matrix(rnorm(replicates * length(par_names)),
                          nrow = replicates, dimnames = list(NULL, par_names))
    saveRDS(prior_alpha, file.path(temp_dir, "prior_samples.rds"))
  }
  temp_dir
}

write_hier_run <- function(dir, par_names, var_names, reps_idx, replicates,
                           write_prior = TRUE) {
  temp_dir <- file.path(dir, "run_temp")
  dir.create(temp_dir, showWarnings = FALSE)
  for (i in reps_idx)
    saveRDS(make_hier_rep(par_names, var_names, i), file.path(temp_dir, paste0("rep_", i, ".rds")))
  if (write_prior) {
    prior_mu  <- matrix(rnorm(length(par_names) * replicates), nrow = length(par_names),
                        dimnames = list(par_names, NULL))
    prior_var <- array(rnorm(length(par_names)^2 * replicates),
                       dim = c(length(par_names), length(par_names), replicates))
    saveRDS(list(prior_mu = prior_mu, prior_var = prior_var),
            file.path(temp_dir, "prior_samples.rds"))
  }
  temp_dir
}

# --- single ------------------------------------------------------------------

test_that("recover_sbc reassembles a complete single run", {
  dir <- tempfile("sbc_"); dir.create(dir)
  par_names <- c("m", "m_lMTRUE", "s", "t0")
  tdir <- write_single_run(dir, par_names, reps_idx = 1:5, replicates = 5)

  save_to <- file.path(dir, "out.RData")
  SBC <- recover_sbc(tdir, par_names, fileName = save_to, type = "single", verbose = FALSE)
  expect_named(SBC, c("rank", "med", "bias", "coverage"))
  expect_equal(nrow(SBC$rank$alpha), 5)
  expect_equal(colnames(SBC$rank$alpha), par_names)
  expect_equal(attr(SBC, "recovered_reps"), 1:5)

  # rows match what split_list_to_dfs would produce directly from the reps
  ref <- split_list_to_dfs(lapply(1:5, function(i) make_single_rep(par_names, i)))
  expect_equal(unname(SBC$rank$alpha), unname(ref$rank$alpha))

  # saved file round-trips in run_sbc's single format (SBC + prior_alpha)
  expect_true(file.exists(save_to))
  env <- new.env(); load(save_to, envir = env)
  expect_true(all(c("SBC", "prior_alpha") %in% ls(env)))
  expect_equal(env$SBC$rank$alpha, SBC$rank$alpha)
})

test_that("recover_sbc returns without saving when fileName is NULL", {
  dir <- tempfile("sbc_"); dir.create(dir)
  par_names <- c("m", "s", "t0")
  tdir <- write_single_run(dir, par_names, reps_idx = 1:4, replicates = 4)

  SBC <- recover_sbc(tdir, par_names, verbose = FALSE)   # no fileName
  expect_equal(nrow(SBC$rank$alpha), 4)
  expect_false(file.exists("recover_sbc.RData"))          # nothing written anywhere
  expect_length(list.files(tdir, pattern = "\\.RData$"), 0)
})

test_that("recover_sbc auto-detects the base of persistent <base>_rep_<i>.rds files", {
  dir <- tempfile("sbc_"); dir.create(dir)
  par_names <- c("m", "s", "t0")
  for (i in 1:4)
    saveRDS(make_single_rep(par_names, i),
            file.path(dir, sprintf("SBC_DDM_singleC_rep_%d.rds", i)))

  SBC <- recover_sbc(dir, par_names, type = "single", verbose = FALSE)
  expect_equal(nrow(SBC$rank$alpha), 4)
  expect_equal(attr(SBC, "recovered_reps"), 1:4)
})

test_that("recover_sbc requires tempdir and design", {
  dir <- tempfile("sbc_"); dir.create(dir)
  write_single_run(dir, c("m", "s"), reps_idx = 1:2, replicates = 2)
  expect_error(recover_sbc(design = c("m", "s")), "tempdir")
  expect_error(recover_sbc(file.path(dir, "run_temp")), "design")
})

test_that("recover_sbc reassembles a partial single run and reports missing", {
  dir <- tempfile("sbc_"); dir.create(dir)
  par_names <- c("m", "s", "t0")
  tdir <- write_single_run(dir, par_names, reps_idx = c(1, 2, 5), replicates = 10)

  expect_message(SBC <- recover_sbc(tdir, par_names, type = "single"),
                 "Recovered 3 single replicate")
  expect_equal(nrow(SBC$rank$alpha), 3)
  expect_equal(attr(SBC, "recovered_reps"), c(1, 2, 5))
})

test_that("recover_sbc auto-detects single and works without a prior file", {
  dir <- tempfile("sbc_"); dir.create(dir)
  par_names <- c("m", "s")
  tdir <- write_single_run(dir, par_names, reps_idx = 1:4, replicates = 4,
                           write_prior = FALSE)

  # inspect-only message when prior draws are absent
  expect_message(SBC <- recover_sbc(tdir, par_names),
                 "Inspect-only")
  expect_equal(nrow(SBC$rank$alpha), 4)
  expect_equal(attr(SBC, "recovered_reps"), 1:4)
})

test_that("recover_sbc skips failed single replicates", {
  dir <- tempfile("sbc_"); dir.create(dir)
  par_names <- c("m", "s")
  temp_dir <- file.path(dir, "run_temp")
  dir.create(temp_dir)
  saveRDS(make_single_rep(par_names, 1), file.path(temp_dir, "rep_1.rds"))
  saveRDS(list(rank = NULL, med = NULL, bias = NULL, coverage = NULL, failed = TRUE),
          file.path(temp_dir, "rep_2.rds"))
  saveRDS(make_single_rep(par_names, 3), file.path(temp_dir, "rep_3.rds"))

  SBC <- recover_sbc(temp_dir, par_names, type = "single", verbose = FALSE)
  expect_equal(nrow(SBC$rank$alpha), 2)
  expect_equal(attr(SBC, "recovered_reps"), c(1, 3))
})

test_that("recover_sbc reads reps from `tempdir` and saves the result elsewhere", {
  in_dir  <- tempfile("sbc_in_");  dir.create(in_dir)
  out_dir <- tempfile("sbc_out_"); dir.create(out_dir)
  par_names <- c("m", "s", "t0")
  tdir <- write_single_run(in_dir, par_names, reps_idx = 1:4, replicates = 4)
  save_to <- file.path(out_dir, "run.RData")

  SBC <- recover_sbc(tdir, par_names, fileName = save_to, type = "single", verbose = FALSE)
  expect_equal(nrow(SBC$rank$alpha), 4)
  expect_true(file.exists(save_to))            # saved to fileName, not in tempdir
  expect_length(list.files(tdir, pattern = "\\.RData$"), 0)
})

# --- hierarchical ------------------------------------------------------------

test_that("recover_sbc reassembles a complete hierarchical run", {
  dir <- tempfile("sbc_"); dir.create(dir)
  par_names <- c("v", "a", "t0")
  var_names <- c("v.v", "a.a", "t0.t0")
  tdir <- write_hier_run(dir, par_names, var_names, reps_idx = 1:6, replicates = 6)

  save_to <- file.path(dir, "out.RData")
  out <- recover_sbc(tdir, par_names, fileName = save_to, type = "hierarchical", verbose = FALSE)
  expect_named(out, c("rank", "prior", "rand_effects"))
  expect_equal(nrow(out$rank$mu), 6)
  expect_equal(colnames(out$rank$mu), par_names)         # from design par names
  expect_equal(colnames(out$rank$var), var_names)        # from rep var_col_names
  expect_false(is.null(out$prior$mu))
  expect_length(out$rand_effects, 6)
  expect_equal(attr(out, "recovered_reps"), 1:6)

  # saved file round-trips in run_sbc's hierarchical format (SBC_temp)
  expect_true(file.exists(save_to))
  env <- new.env(); load(save_to, envir = env)
  expect_true("SBC_temp" %in% ls(env))
  expect_equal(env$SBC_temp$rank$mu, out$rank$mu)
})

test_that("recover_sbc recovers hierarchical ranks without prior (prior NULL + warning)", {
  dir <- tempfile("sbc_"); dir.create(dir)
  par_names <- c("v", "a")
  var_names <- c("v.v", "a.a")
  tdir <- write_hier_run(dir, par_names, var_names, reps_idx = 1:3, replicates = 3,
                         write_prior = FALSE)

  expect_warning(out <- recover_sbc(tdir, par_names, type = "hierarchical", verbose = FALSE),
                 "Prior draws not found")
  expect_equal(nrow(out$rank$mu), 3)
  expect_equal(colnames(out$rank$var), var_names)
  expect_null(out$prior$mu)
})

# --- errors ------------------------------------------------------------------

test_that("recover_sbc errors when no rep files exist", {
  dir <- tempfile("sbc_"); dir.create(dir)
  expect_error(recover_sbc(dir, c("m", "s"), verbose = FALSE),
               "No replicate files found")
})
