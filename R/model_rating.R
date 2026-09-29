# Rating layer: model-independent support for choice + rating data.
#
# A rating model observes, on each trial, a choice R, its response time rt and a
# rating RR, an integer 1..K with 1 the LOWEST and K the HIGHEST rating. The
# rating is read from K - 1 ordered thresholds, expressed as proportions of the
# choice threshold, 0 = d_0 < d_1 < ... < d_{K-1} < d_K = 1, on the balance of
# evidence at the decision: a rating model reports RR = c when that evidence
# (as a proportion of the threshold) lies in (d_{K-c}, d_{K-c+1}), so a small
# proportion means high confidence.
#
# Sampling scale: criteria on the log balance-of-evidence scale. Criterion k
# separates rating k from rating k + 1; the first, -log(d_{K-1}), is exp(c1)
# and each further one adds a positive increment exp(ck), so the thresholds are
# ordered by construction and an effect on c1 moves all of them together:
#
#     d_{K-k} = exp(-(c1 + ... + ck)),  k = 1 .. K - 1  (natural-scale c's)
#
# The sampled p_types are c1 .. c{K-1}; the model's Ttransform adds the
# natural-scale thresholds d1 .. d{K-1} (d1 the smallest, the paper's and DMC's
# numbering), so they can be reported with add_recalculated = TRUE.
#
# Nothing here depends on the process model that produces the evidence.

# Number of ratings K: a single integer >= 1.
rating_check_n <- function(n_ratings) {
  if (length(n_ratings) != 1L)
    stop("n_ratings must be a single integer: unequal rating counts per choice are not yet supported")
  if (!is.numeric(n_ratings) || is.na(n_ratings) || n_ratings < 1 ||
      n_ratings != round(n_ratings))
    stop("n_ratings must be an integer >= 1")
  as.integer(n_ratings)
}

# Sampled criterion p_types c1 .. c{K-1} (character(0) when K = 1).
rating_criterion_names <- function(K) {
  if (K < 2) character(0) else paste0("c", seq_len(K - 1))
}

# Natural-scale threshold names d1 .. d{K-1} (d1 the smallest).
rating_threshold_names <- function(K) {
  if (K < 2) character(0) else paste0("d", seq_len(K - 1))
}

# Model-list pieces for the criteria: default values (sampled scale), transforms
# and bounds. Defaults give d_{K-1} = exp(-0.5) and each further threshold a
# factor exp(-0.5) lower.
rating_p_types <- function(K) {
  nams <- rating_criterion_names(K)
  stats::setNames(rep(log(0.5), length(nams)), nams)
}

rating_transform <- function(K) {
  nams <- rating_criterion_names(K)
  stats::setNames(rep("exp", length(nams)), nams)
}

rating_bound <- function(K) {
  nams <- rating_criterion_names(K)
  matrix(rep(c(0, Inf), length(nams)), nrow = 2,
         dimnames = list(NULL, nams))
}

# Ttransform piece: add the natural-scale thresholds d1 .. d{K-1} to a
# parameter matrix holding the natural-scale criteria c1 .. c{K-1} (row-wise).
rating_add_thresholds <- function(pars, K) {
  if (K < 2) return(pars)
  cn <- rating_criterion_names(K)
  C <- pars[, cn, drop = FALSE]
  if (K > 2) for (k in 2:(K - 1)) C[, k] <- C[, k - 1] + C[, k]
  # criterion k gives d_{K-k}: reverse so column j is d_j
  d <- exp(-C[, rev(seq_len(K - 1)), drop = FALSE])
  colnames(d) <- rating_threshold_names(K)
  keep <- attributes(pars)[c("ok")]
  pars <- cbind(pars[, !(colnames(pars) %in% colnames(d)), drop = FALSE], d)
  if (!is.null(keep$ok)) attr(pars, "ok") <- keep$ok
  pars
}

# Lower and upper threshold proportions bounding rating RR, for thresholds
# d (a matrix with columns d1 .. d{K-1}, one row per trial):
# RR = c  <=>  d_{K-c} < evidence proportion < d_{K-c+1}, d_0 = 0, d_K = 1.
rating_interval <- function(d, RR, K) {
  n <- length(RR)
  full <- cbind(rep(0, n), if (K > 1) d[, rating_threshold_names(K), drop = FALSE], rep(1, n))
  # column j + 1 of full holds d_j
  lo <- full[cbind(seq_len(n), K - RR + 1)]
  hi <- full[cbind(seq_len(n), K - RR + 2)]
  list(lower = lo, upper = hi)
}

# Rating from an evidence proportion E in [0, 1] (row-wise thresholds d):
# the number m of thresholds below E puts E in (d_m, d_{m+1}), so RR = K - m.
rating_from_evidence <- function(E, d, K) {
  if (K < 2) return(rep(1, length(E)))
  d <- d[, rating_threshold_names(K), drop = FALSE]
  K - rowSums(d < E)
}

# Fold choice and rating into one ordered factor running from the highest
# rating of the first response level to the highest rating of the second, e.g.
# for K = 3: left.3, left.2, left.1, right.1, right.2, right.3.
rating_fold <- function(R, RR, K) {
  if (!is.factor(R)) stop("R must be a factor")
  Rl <- levels(R)
  labs <- unlist(lapply(seq_along(Rl), function(i) {
    ks <- if (i == 1) K:1 else 1:K
    paste(Rl[i], ks, sep = ".")
  }))
  factor(paste(as.character(R), RR, sep = "."), levels = labs, ordered = TRUE)
}

rating_unfold <- function(x) {
  s <- as.character(x)
  cut <- regexpr("\\.[^.]*$", s)
  R <- substr(s, 1, cut - 1)
  RR <- as.integer(substr(s, cut + 1, nchar(s)))
  lev <- unique(sub("\\.[^.]*$", "", levels(x)))
  data.frame(R = factor(R, levels = lev), RR = RR)
}

# Data checks for the rating column (called from design_model through the
# model's check_data hook).
rating_check_data <- function(data, K) {
  if (!"RR" %in% names(data)) stop("Rating models need a rating column RR in the data")
  RR <- data$RR
  if (is.factor(RR) || !is.numeric(RR))
    stop("RR must be numeric (integer-valued, 1 = lowest rating); convert a factor with as.numeric(as.character(RR))")
  obs <- !is.na(data$R)
  if (any(is.na(RR[obs]))) stop("RR is NA on a trial with an observed response R")
  rr <- RR[obs & !is.na(RR)]
  if (any(rr != round(rr))) stop("RR must be integer-valued")
  if (any(rr < 1 | rr > K))
    stop("RR must lie in 1 .. ", K, " (n_ratings = ", K, ")")
  unused <- setdiff(seq_len(K), rr)
  if (length(unused))
    warning("Rating(s) ", paste(unused, collapse = ", "), " never occur in the data")
  invisible(TRUE)
}

# prepare_design hook factory: refuse other than n_choices responses and keep
# trial-level parameters (e.g. a correlation) off the accumulator factors.
rating_prepare_design <- function(n_choices = 2, trial_level = character(0)) {
  function(formula, constants, Rlevels = NULL, ...) {
    if (!is.null(n_choices) && length(Rlevels) != n_choices)
      stop("This rating model needs exactly ", n_choices, " response levels (Rlevels has ",
           length(Rlevels), ")")
    for (f in formula) {
      lhs <- as.character(stats::terms(f)[[2]])
      rhs <- all.vars(f[[3]])
      if ("RR" %in% rhs) stop("RR is a response and cannot be a predictor (", lhs, ")")
      if (lhs %in% trial_level && any(c("lR", "lM") %in% rhs))
        stop(lhs, " is a trial-level parameter and cannot depend on lR or lM")
    }
    list(formula = formula, constants = constants)
  }
}

# Stable log(exp(la) - exp(lb)) for la >= lb.
log_diff_exp <- function(la, lb) {
  x <- lb - la
  out <- la + ifelse(x > -log(2), log(-expm1(x)), log1p(-exp(x)))
  out[la == lb] <- -Inf
  out
}

# log P(z_hi < Z < z_lo) for a standard normal Z and z_hi <= z_lo, accurate in
# both tails and across 0 (Phi(z) - 1/2 = sign(z) * pchisq(z^2, 1) / 2).
log_pnorm_interval <- function(z_hi, z_lo) {
  out <- numeric(length(z_hi))
  up <- z_hi >= 0
  lo <- z_lo <= 0 & !up
  mid <- !up & !lo
  if (any(up)) out[up] <- log_diff_exp(stats::pnorm(z_hi[up], lower.tail = FALSE, log.p = TRUE),
                                       stats::pnorm(z_lo[up], lower.tail = FALSE, log.p = TRUE))
  if (any(lo)) out[lo] <- log_diff_exp(stats::pnorm(z_lo[lo], log.p = TRUE),
                                       stats::pnorm(z_hi[lo], log.p = TRUE))
  if (any(mid)) out[mid] <- log((stats::pchisq(z_hi[mid]^2, 1) + stats::pchisq(z_lo[mid]^2, 1)) / 2)
  out
}
