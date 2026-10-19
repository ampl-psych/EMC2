# Regenerates vignettes/data/using-a-trained-flow.RData, the shown output of the
# vignette "Using a trained flow", from the card it uses (ddm_svst0_e299, not
# shipped with EMC2; in the NLE repository, handover/nle-methods/port/export/).
# The RData holds summaries only: the card's fields, the printed output of each
# step, and the posterior summaries and timings of the two fits.
#   R_LIBS=<library with EMC2 >= the card's `fixed` support> R_PROFILE_USER=/dev/null \
#     Rscript inst/scripts/using_a_trained_flow_data.R <ddm_svst0_e299_weights.json> \
#       [vignettes/data/using-a-trained-flow.RData]
# Six cores: three chains, two cores each.
suppressMessages(library(EMC2))
a <- commandArgs(TRUE)
if (!length(a)) stop("usage: Rscript using_a_trained_flow_data.R <card.json> [out.RData]")
card <- normalizePath(a[1], mustWork = TRUE)
out <- if (length(a) > 1) a[2] else file.path("vignettes", "data", "using-a-trained-flow.RData")
SHA256 <- "67739e2dca9da692e0b241ea16389eadc5933d68e6873f535f47bef314a04374"
CORES_PER_CHAIN <- 2; N_CHAINS <- 3
RNGkind("L'Ecuyer-CMRG")    # the chains run in forked processes: reproducible streams need this generator

# Runs expr and keeps what it printed, its messages and its warnings.
grab <- function(expr) {
  msgs <- character(); warns <- character()
  printed <- utils::capture.output(value <- withCallingHandlers(expr,
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") },
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") }))
  list(value = value, printed = printed, messages = msgs, warnings = warns)
}

# --- 1. the card -----------------------------------------------------------------------
j <- jsonlite::fromJSON(card)
nat <- c(j$context_names, names(j$fixed))
box <- data.frame(input = j$context_names, transform = unlist(j$context_transforms[j$context_names]),
                  lower = j$bounds_natural$lower[seq_along(j$context_names)],
                  upper = j$bounds_natural$upper[seq_along(j$context_names)], row.names = NULL)
box$zerobox_c <- unlist(j$zerobox_c)[box$input]
card_info <- list(file = basename(card), sha256 = SHA256, card_version = j$card_version, kind = j$kind,
                  model = j$model, label = j$label, context_names = j$context_names,
                  context_encoding = j$context_encoding, fixed = unlist(j$fixed), box = box,
                  spline = j$spline[c("num_bins", "range_min", "range_max")], ll_floor_log = j$ll_floor_log,
                  training = list(epochs = j$training$epochs, val_loss = j$training$val_loss,
                                  checkpoint = j$checkpoint_step_flow), export_gate = j$export_gate,
                  exporter = j$exporter_version, fields = names(j))

reg <- grab(register_nn_model(card, twin = DDM))
m <- reg$value
if (!identical(m()$nn$sha256, SHA256)) stop("the card's sha256 is ", m()$nn$sha256, ", not ", SHA256)
reg_out <- list(messages = reg$messages, warnings = reg$warnings, p_types = m()$p_types,
                exception = m()$nn$exception, fixed = m()$nn$fixed)
golden <- grab(print(nn_golden_check(m)))
verify <- system2(file.path(R.home("bin"), "Rscript"),
                  c("-e", shQuote('source(system.file("scripts", "nle_verify.R", package = "EMC2"))'), shQuote(card)),
                  stdout = TRUE, stderr = TRUE)
verify <- sub(" @ .*$", "", verify)                        # not the local library path

# --- 2. a design and a prior -------------------------------------------------------------
FORM <- list(v ~ S, a ~ 1, t0 ~ 1, Z ~ 1, sv ~ 1, st0 ~ 1)
CONST <- c(s = log(1), SZ = qnorm(0))                      # the card's fixed: s = 1, SZ = 0
mk_design <- function(model) design(factors = list(subjects = 1, S = c("left", "right")),
                                    Rlevels = c("left", "right"), formula = FORM, constants = CONST,
                                    model = model)
des_nn <- grab(mk_design(m))
des_an <- grab(mk_design(DDM))
PMEAN <- c(v = 0, v_Sright = 1, a = log(1.2), t0 = log(.25), Z = 0, sv = log(.8), st0 = log(.1))
WIDE <- c(v = 1, v_Sright = 1, a = .3, t0 = .2, Z = .2, sv = .3, st0 = .3)
PSD <- c(v = .6, v_Sright = .6, a = .15, t0 = .2, Z = .2, sv = .3, st0 = .3)
pri_wide <- prior(des_nn$value, type = "single", pmean = PMEAN, psd = WIDE)
pri_nn <- prior(des_nn$value, type = "single", pmean = PMEAN, psd = PSD)
pri_an <- prior(des_an$value, type = "single", pmean = PMEAN, psd = PSD)
box_wide <- grab(nn_in_box(des_nn$value, pri_wide))
box_ok <- grab(nn_in_box(des_nn$value, pri_nn))

# --- 3. one subject, two fits ------------------------------------------------------------
TRUTH <- c(v = -1, v_Sright = 2, a = log(1.3), t0 = log(.3), Z = qnorm(.45), sv = log(1), st0 = log(.15))
set.seed(20261020)
dat <- make_data(TRUTH, des_an$value, n_trials = 200)      # 400 trials from the analytic DDM
data_summary <- list(n = nrow(dat), table = table(S = dat$S, R = dat$R),
                     rt = aggregate(rt ~ S + R, dat, function(x) round(c(mean = mean(x), min = min(x)), 3)))
coarse <- grab(make_emc(dat, des_nn$value, prior_list = pri_nn, n_chains = N_CHAINS))   # rt_resolution 1/60

posterior <- function(emc) {
  x <- sapply(get_pars(emc, selection = "alpha", merge_chains = TRUE), function(p) as.numeric(p[[1]]))
  data.frame(mean = colMeans(x), sd = apply(x, 2, sd), lower = apply(x, 2, quantile, .025),
             upper = apply(x, 2, quantile, .975))
}
fit_one <- function(design, pri, seed) {
  emc <- make_emc(dat, design, prior_list = pri, n_chains = N_CHAINS, rt_resolution = 1 / 1000)
  set.seed(seed)
  t <- system.time(emc <- suppressMessages(fit(emc, cores_per_chain = CORES_PER_CHAIN,
                                               cores_for_chains = N_CHAINS, verbose = FALSE)))
  list(post = posterior(emc), seconds = unname(t[["elapsed"]]),
       gd_max = gd_summary(emc, selection = "alpha")[1, "max"],
       ess_min = ess_summary(emc, selection = "alpha")[1, "min"],
       iter = chain_n(emc)[, "sample"])
}
fits <- list(card = fit_one(des_nn$value, pri_nn, 1), analytic = fit_one(des_an$value, pri_an, 1))
pc <- fits$card$post; pa <- fits$analytic$post
overlap <- pmax(0, pmin(pc$upper, pa$upper) - pmax(pc$lower, pa$lower)) / (pa$upper - pa$lower)
comparison <- data.frame(parameter = rownames(pa), true = TRUTH[rownames(pa)],
                         card_mean = pc$mean, card_lower = pc$lower, card_upper = pc$upper,
                         analytic_mean = pa$mean, analytic_lower = pa$lower, analytic_upper = pa$upper,
                         shift_sd = (pc$mean - pa$mean) / pa$sd, sd_ratio = pc$sd / pa$sd,
                         overlap = overlap, row.names = NULL)
for (f in names(fits)) fits[[f]]$post <- NULL

env <- list(R = R.version.string, EMC2 = as.character(utils::packageVersion("EMC2")),
            machine = paste(Sys.info()[c("sysname", "machine")], collapse = " "),
            cores = CORES_PER_CHAIN * N_CHAINS, made = format(Sys.Date()))
ut <- list(card_info = card_info, reg = reg_out, golden = golden$printed, golden_msg = golden$messages,
           verify = verify, design_nn = des_nn$printed, design_messages = des_nn$messages,
           box_wide = box_wide$value, box_wide_warning = box_wide$warnings, box_ok = box_ok$value,
           truth = TRUTH, pmean = PMEAN, psd = PSD, data_summary = data_summary,
           coarse_warning = coarse$warnings, comparison = comparison, fits = fits, env = env)
save(ut, file = out, compress = "xz")
cat(sprintf("wrote %s (%.1f kB)\n", out, file.size(out) / 1024))
print(comparison, digits = 3)
str(fits)
