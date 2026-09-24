# Reference R implementation of the exported DDM flow + choice classifier.
#
# Evaluates the two-network DDM likelihood exported by
# flow/scripts/export_flow_ddm.py (see flow/src/flow/config.py
# `joint_log_prob`): a spline flow over rt conditioned on
# [standardized params, raw response code R in {1, 2}], and a classifier
# MLP mapping standardized params to the logit of P(R = 2). With
# u = log(rt), z = g_R^{-1}(u) (the spline for response R):
#
#   joint pdf(rt, R)      = dnorm(z) * |dz/du| / rt * P(R)
#   defective cdf(rt, R)  = pnorm(z) * P(R)        (-> P(R) as rt -> Inf)
#
# Bounding box exactly as in flow_race.R: outside the training region both
# pdf and cdf are 0 (log_pdf = -Inf).
#
# Requires flow_race.R to be sourced first (shares .gelu_tanh, .softplus,
# .layernorm, .rqs_inverse).

ddm_load <- function(path) {
  fl <- jsonlite::fromJSON(path, simplifyDataFrame = FALSE)
  fix_layers <- function(mlp) {
    mlp$layers <- lapply(mlp$layers, function(l)
      list(W = as.matrix(l$W), b = as.numeric(l$b)))
    mlp
  }
  fl$flow_mlp <- fix_layers(fl$flow_mlp)
  fl$classifier_mlp <- fix_layers(fl$classifier_mlp)
  stopifnot(fl$spline$num_splines == 1,
            fl$spline$output_transform == "exp",
            fl$spline$base_distribution == "standard_normal",
            !isTRUE(fl$flow_mlp$use_norm) || !is.null(fl$flow_mlp$norms))
  fl
}

.mlp_fwd <- function(mlp, x, use_norm = FALSE, norms = NULL) {
  n_layers <- length(mlp$layers)
  h <- x
  for (i in seq_len(n_layers - 1L)) {
    l <- mlp$layers[[i]]
    h <- h %*% l$W + matrix(l$b, nrow(h), length(l$b), byrow = TRUE)
    if (use_norm) {
      nrm <- norms[[i]]
      h <- .layernorm(h, nrm$scale, nrm$bias, nrm$eps)
    }
    h <- .gelu_tanh(h)
  }
  l <- mlp$layers[[n_layers]]
  h %*% l$W + matrix(l$b, nrow(h), length(l$b), byrow = TRUE)
}

.build_knots_row <- function(raw, sp) {
  K <- sp$num_bins
  range_size <- sp$range_max - sp$range_min
  norm_bins <- function(u) {
    e <- exp(u - max(u))
    (e / sum(e)) * (range_size - K * sp$min_bin_size) + sp$min_bin_size
  }
  widths <- norm_bins(raw[1:K])
  heights <- norm_bins(raw[(K + 1):(2 * K)])
  offset <- log(exp(1 - sp$min_knot_slope) - 1)
  slopes <- .softplus(raw[(2 * K + 1):(3 * K + 1)] + offset) + sp$min_knot_slope
  slopes[1L] <- 1
  slopes[K + 1L] <- 1
  list(x_pos = c(sp$range_min, sp$range_min + cumsum(widths[-K]), sp$range_max),
       y_pos = c(sp$range_min, sp$range_min + cumsum(heights[-K]), sp$range_max),
       slopes = slopes)
}

ddm_in_box <- function(fl, theta) {
  all(theta >= fl$bounds_sampled$lower & theta <= fl$bounds_sampled$upper)
}

# Amortizable conditioning: one parameter vector -> spline knots for BOTH
# responses plus log P(R = 1) and log P(R = 2).
# EMC2 (v, log a, log t0, log s, qnorm Z, qnorm SZ, log sv, log st0) ->
# (v, log SZ0, log t0, log s, log SZ, log SZ1, log sv, log st0). Mirrors
# emc2_to_triple() in flow/scripts/export_flow_ddm.py. SZ0 and SZ1 are
# non-negative by construction because sw <= 2*min(w, 1-w) always holds; the
# pmax only guards float cancellation.
.emc2_to_triple <- function(theta) {
  a  <- exp(theta[2L])
  w  <- pnorm(theta[5L])
  sw <- 2 * pnorm(theta[6L]) * min(w, 1 - w)
  SZ  <- sw * a
  SZ0 <- a * w - SZ / 2
  SZ1 <- a - SZ0 - SZ
  eps <- 1e-12
  theta[2L] <- log(max(SZ0, eps))
  theta[5L] <- log(max(SZ,  eps))
  theta[6L] <- log(max(SZ1, eps))
  theta
}

# theta is ALWAYS supplied on EMC2's scale - that is the port's public
# interface. Under context_encoding == "triple" the flow is conditioned on the
# re-encoded parameters while the frozen classifier keeps the EMC2 ones, so the
# two heads use different vectors AND different scalers.
ddm_condition <- function(fl, theta) {
  triple <- identical(fl$context_encoding, "triple")
  theta_flow <- if (triple) .emc2_to_triple(theta) else theta
  ctx_std <- (theta_flow - fl$scaler$mean) / fl$scaler$scale
  clf_std <- if (triple)
    (theta - fl$classifier_scaler$mean) / fl$classifier_scaler$scale else ctx_std
  kn <- lapply(c(1, 2), function(R) {
    raw <- .mlp_fwd(fl$flow_mlp, matrix(c(ctx_std, R), nrow = 1L),
                    use_norm = isTRUE(fl$flow_mlp$use_norm),
                    norms = fl$flow_mlp$norms)
    .build_knots_row(as.numeric(raw), fl$spline)
  })
  logit <- as.numeric(.mlp_fwd(fl$classifier_mlp, matrix(clf_std, nrow = 1L)))
  # log P(R=1) = -softplus(logit); log P(R=2) = -softplus(-logit)
  list(knots = kn, log_p_R = c(-.softplus(logit), -.softplus(-logit)),
       ctx_std = ctx_std, logit = logit)
}

# Main evaluator: one parameter vector (sampled scale, order
# fl$context_names), rt vector, R vector (1/2, recycled if scalar).
ddm_eval <- function(fl, theta, rt, R) {
  n <- length(rt)
  R <- rep_len(as.integer(R), n)
  if (!ddm_in_box(fl, theta)) {
    return(list(pdf = numeric(n), cdf = numeric(n),
                log_pdf = rep(-Inf, n), p_R = rep(NA_real_, n),
                z = rep(NA_real_, n), in_box = FALSE))
  }
  cond <- ddm_condition(fl, theta)
  u <- log(rt)
  z <- logdet <- numeric(n)
  for (r in 1:2) {
    sel <- R == r
    if (any(sel)) {
      inv <- .rqs_inverse(u[sel], cond$knots[[r]]$x_pos,
                          cond$knots[[r]]$y_pos, cond$knots[[r]]$slopes)
      z[sel] <- inv$x
      logdet[sel] <- inv$logdet
    }
  }
  log_p_R <- cond$log_p_R[R]
  log_pdf <- dnorm(z, log = TRUE) + logdet - u + log_p_R
  list(pdf = exp(log_pdf),
       cdf = pnorm(z) * exp(log_p_R),   # defective CDF
       log_pdf = log_pdf,
       p_R = exp(log_p_R),
       z = z, in_box = TRUE)
}
