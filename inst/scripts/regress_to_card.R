# Convert a direct-regression net from the NLE training pipeline into an EMC2
# model card for register_nn_model() (kind "regression_joint").
#
# The training output (results/model_states/<tag>_budget/regress.json, shipped
# in the NLE handover as artefacts/<tag>_regress.json) holds `dims`, `layers`
# (W as n_in x n_out, b), `scaler_mean`/`scaler_scale` over (context on the
# sampled scale, log rt), `context` and `epoch`. The network's input is that
# standardised vector with the raw response code (1 = lower, 2 = upper)
# appended unscaled; its single output is the joint log density
# log p(rt, R | theta), hidden layers GELU (tanh approximation), linear output.
# This is the evaluator of the NLE project's evaluation/DDMreg.R.
#
# The training output records neither the parameter scales nor the training
# box, so both come from a reference card trained on the same simulations
# (reg_s402 and ddm_cap256w_c4 both train on data/train/tpt40p400k/DDM; the
# sampled-scale ranges of its parameters.parquet fill cap256w's box). The
# reference card's context_names must equal the net's `context`.
#
# Usage (Rscript, jsonlite needed):
#   Rscript regress_to_card.R <regress.json> <reference card .rds> <out.rds> [budget]
# e.g. the shipped ddm_reg_s402.rds:
#   Rscript regress_to_card.R reg_s402_regress.json \
#     $(Rscript -e 'cat(system.file("extdata","flownn","ddm_cap256w_c4.rds",package="EMC2"))') \
#     ddm_reg_s402.rds reg_s402

regress_to_card <- function(json, reference, budget = NULL) {
  net <- jsonlite::fromJSON(json, simplifyVector = TRUE)
  ref <- readRDS(reference)
  if (!is.null(ref$members)) ref <- ref$members[[1]]
  ctx <- as.character(net$context)
  k <- length(ctx)
  if (!identical(ctx, as.character(ref$context_names)))
    stop("the net's context (", paste(ctx, collapse = ", "), ") differs from the reference card's ",
         "context_names (", paste(ref$context_names, collapse = ", "), ")")
  if (length(net$scaler_mean) != k + 1L || length(net$scaler_scale) != k + 1L)
    stop("expected a scaler over the ", k, " contexts plus log rt")
  layers <- lapply(seq_along(net$layers$W), function(i)
    list(W = as.matrix(net$layers$W[[i]]), b = as.numeric(net$layers$b[[i]])))
  if (nrow(layers[[1]]$W) != k + 2L || ncol(layers[[length(layers)]]$W) != 1L)
    stop("expected ", k + 2L, " inputs (contexts, log rt, R) and one output")
  if (!is.null(net$dims) && !identical(as.integer(net$dims),
                                       c(nrow(layers[[1]]$W), vapply(layers, function(l) ncol(l$W), 0L))))
    stop("layer shapes disagree with dims")
  hash <- if (exists("sha256sum", envir = asNamespace("tools"), inherits = FALSE))
    unname(tools::sha256sum(json)) else NA_character_
  list(
    model = "DDM", kind = "regression_joint",
    budget = if (is.null(budget)) sub("_regress\\.json$|\\.json$", "", basename(json)) else budget,
    checkpoint_step = as.integer(net$epoch),
    exported = format(Sys.Date()),
    source = basename(json), source_sha256 = hash,
    context_names = ctx,
    context_transforms = ref$context_transforms,
    parameter_scale = "sampled",
    bounds_natural = ref$bounds_natural,
    bounds_sampled = ref$bounds_sampled,
    response_coding = "raw R (1 = lower, 2 = upper) appended unscaled after the standardised (context, log rt)",
    input_layout = c(ctx, "log_rt", "R"),
    # the whole input vector: R passes through as (R - 0) / 1
    scaler = list(mean = c(as.numeric(net$scaler_mean), 0), scale = c(as.numeric(net$scaler_scale), 1)),
    response_values = c(1, 2),
    mlp = list(activation = "gelu_tanh", use_norm = FALSE, layers = layers))
}

if (!interactive() && sys.nframe() == 0L) {
  args <- commandArgs(TRUE)
  if (length(args) < 3L) stop("usage: regress_to_card.R <regress.json> <reference card .rds> <out.rds> [budget]")
  card <- regress_to_card(args[1], args[2], if (length(args) >= 4L) args[4] else NULL)
  saveRDS(card, args[3])
  cat("wrote", args[3], "\n")
}
