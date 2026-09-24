# Reference R implementation of exported race-model flows (LNR, RDM, ...).
#
# Evaluates any flow exported by flow/scripts/export_flow.py: a single
# rational-quadratic spline flow over rt, conditioned via an MLP on the model
# parameters. This is the readable spec for the C++ (Rcpp) port
# (port/src/flow_race.cpp); it is validated against float64 golden vectors
# from Python (port/tests/validate_r.R).
#
# Math (see flow/src/flow/config.py and distrax rational_quadratic_spline.py):
#   rt = exp(g(z)),  z ~ N(0,1),  g = rational-quadratic spline on
#   [range_min, range_max] with identity tails. The MLP maps standardized
#   parameters to the 3K+1 unconstrained spline parameters.
# Hence with u = log(rt), z = g^{-1}(u):
#   pdf(rt) = dnorm(z) * |dg^{-1}/du| / rt
#   cdf(rt) = pnorm(z)          (bijector is monotone increasing)
#   sf(rt)  = pnorm(z, lower = FALSE)
#
# Bounding box: parameter vectors outside the training region (recorded in
# the JSON, sampled scale) return pdf = 0 and cdf = 0 immediately
# (log_pdf = -Inf, log_sf = 0, consistent with cdf = 0). Rejection sentinel,
# not a statement about the flow's extrapolation.

flow_load <- function(path) {
  fl <- jsonlite::fromJSON(path, simplifyDataFrame = FALSE)
  fl$mlp$layers <- lapply(fl$mlp$layers, function(l)
    list(W = as.matrix(l$W), b = as.numeric(l$b)))
  if (isTRUE(fl$mlp$use_norm))
    fl$mlp$norms <- lapply(fl$mlp$norms, function(n)
      list(scale = as.numeric(n$scale), bias = as.numeric(n$bias),
           eps = as.numeric(n$eps)))
  stopifnot(fl$spline$num_splines == 1,
            fl$spline$output_transform == "exp",
            fl$spline$base_distribution == "standard_normal")
  fl
}

.gelu_tanh <- function(x) {
  # jax.nn.gelu(approximate=TRUE), the variant used in training
  0.5 * x * (1 + tanh(sqrt(2 / pi) * (x + 0.044715 * x^3)))
}

.softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

.layernorm <- function(x, scale, bias, eps) {
  mu <- rowMeans(x)
  v <- rowMeans((x - mu)^2)
  (x - mu) / sqrt(v + eps) * matrix(scale, nrow(x), length(scale), byrow = TRUE) +
    matrix(bias, nrow(x), length(bias), byrow = TRUE)
}

.row_cumsum <- function(m) {
  if (nrow(m) == 1L) matrix(cumsum(m[1L, ]), nrow = 1L) else t(apply(m, 1L, cumsum))
}

# Conditioner: parameter matrix (n x K_ctx, sampled scale) -> spline knots.
# This is the amortizable part: call once per unique parameter vector, then
# evaluate any number of rts against the returned knots.
flow_knots <- function(fl, theta) {
  theta <- matrix(theta, ncol = length(fl$context_names))
  sp <- fl$spline
  K <- sp$num_bins
  range_size <- sp$range_max - sp$range_min

  h <- sweep(sweep(theta, 2L, fl$scaler$mean, "-"), 2L, fl$scaler$scale, "/")
  ctx_std <- h
  n_layers <- length(fl$mlp$layers)
  for (i in seq_len(n_layers - 1L)) {
    l <- fl$mlp$layers[[i]]
    h <- h %*% l$W + matrix(l$b, nrow(h), length(l$b), byrow = TRUE)
    if (isTRUE(fl$mlp$use_norm)) {
      nrm <- fl$mlp$norms[[i]]
      h <- .layernorm(h, nrm$scale, nrm$bias, nrm$eps)
    }
    h <- .gelu_tanh(h)
  }
  l <- fl$mlp$layers[[n_layers]]
  raw <- h %*% l$W + matrix(l$b, nrow(h), length(l$b), byrow = TRUE)

  norm_bins <- function(u) {
    e <- exp(u - apply(u, 1L, max))
    (e / rowSums(e)) * (range_size - K * sp$min_bin_size) + sp$min_bin_size
  }
  widths  <- norm_bins(raw[, 1:K, drop = FALSE])
  heights <- norm_bins(raw[, (K + 1):(2 * K), drop = FALSE])

  offset <- log(exp(1 - sp$min_knot_slope) - 1)
  slopes <- .softplus(raw[, (2 * K + 1):(3 * K + 1), drop = FALSE] + offset) +
    sp$min_knot_slope
  slopes[, 1L] <- 1        # boundary_slopes = "identity"
  slopes[, K + 1L] <- 1

  x_pos <- cbind(sp$range_min,
                 sp$range_min + .row_cumsum(widths[, -K, drop = FALSE]),
                 sp$range_max)
  y_pos <- cbind(sp$range_min,
                 sp$range_min + .row_cumsum(heights[, -K, drop = FALSE]),
                 sp$range_max)

  list(x_pos = x_pos, y_pos = y_pos, slopes = slopes,
       raw = raw, ctx_std = ctx_std, num_bins = K)
}

# Inverse spline + log|d inverse/dy| for one knot row, vectorized over y.
# Mirrors distrax _rational_quadratic_spline_inv (incl. stable quadratic root).
.rqs_inverse <- function(y, x_pos, y_pos, slopes) {
  K <- length(x_pos) - 1L
  below <- y <= y_pos[1L]
  above <- y >= y_pos[K + 1L]
  idx <- pmin(pmax(findInterval(y, y_pos), 1L), K)

  x_lo <- x_pos[idx]; x_hi <- x_pos[idx + 1L]
  y_lo <- y_pos[idx]; y_hi <- y_pos[idx + 1L]
  d_lo <- slopes[idx]; d_hi <- slopes[idx + 1L]

  bin_width  <- x_hi - x_lo
  bin_height <- y_hi - y_lo
  bin_slope  <- bin_height / bin_width

  w <- pmin(pmax((y - y_lo) / bin_height, 0), 1)
  slopes_term <- d_hi + d_lo - 2 * bin_slope
  cc <- -bin_slope * w
  bb <- d_lo - slopes_term * w
  aa <- bin_slope - bb

  sqrt_diff <- bb^2 - 4 * aa * cc
  safe_sqrt <- sqrt(pmax(sqrt_diff, .Machine$double.xmin))
  safe_sqrt[sqrt_diff <= 0] <- 0
  num <- ifelse(bb >= 0, 2 * cc, -bb + safe_sqrt)
  den <- ifelse(bb >= 0, -bb - safe_sqrt, 2 * aa)
  z <- pmin(pmax(num / den, 0), 1)

  x <- x_lo + bin_width * z

  sq_z <- z * z
  z1mz <- z - sq_z
  sq_1mz <- (1 - z)^2
  denominator <- bin_slope + slopes_term * z1mz
  logdet <- -2 * log(bin_slope) -
    log(d_hi * sq_z + 2 * bin_slope * z1mz + d_lo * sq_1mz) +
    2 * log(denominator)

  # identity tails (boundary slopes are 1)
  x[below] <- y[below] - y_pos[1L] + x_pos[1L]
  x[above] <- y[above] - y_pos[K + 1L] + x_pos[K + 1L]
  logdet[below | above] <- 0

  list(x = x, logdet = logdet)
}

in_box <- function(fl, theta) {
  all(theta >= fl$bounds_sampled$lower & theta <= fl$bounds_sampled$upper)
}

# Main evaluator: one parameter vector (sampled scale), many rts.
# Returns pdf, cdf, and their log/survivor companions.
flow_eval <- function(fl, theta, rt) {
  n <- length(rt)
  if (!in_box(fl, theta)) {
    return(list(pdf = numeric(n), cdf = numeric(n),
                log_pdf = rep(-Inf, n), log_sf = numeric(n),
                z = rep(NA_real_, n), in_box = FALSE))
  }
  kn <- flow_knots(fl, theta)
  u <- log(rt)
  inv <- .rqs_inverse(u, kn$x_pos[1L, ], kn$y_pos[1L, ], kn$slopes[1L, ])
  log_pdf <- dnorm(inv$x, log = TRUE) + inv$logdet - u
  list(pdf = exp(log_pdf),
       cdf = pnorm(inv$x),
       log_pdf = log_pdf,
       log_sf = pnorm(inv$x, lower.tail = FALSE, log.p = TRUE),
       z = inv$x, in_box = TRUE)
}
