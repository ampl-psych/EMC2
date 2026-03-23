.flat_loc_stats_per_col <- function(mat, p1 = 1/3, p2 = 1/3) {
  if (!is.matrix(mat) || ncol(mat) == 0) return(numeric(0))
  if (nrow(mat) < 10) return(rep(Inf, ncol(mat)))

  n <- nrow(mat)
  xlen <- max(1L, round(n * p1))
  ylen <- max(1L, round(n * p2))

  apply(mat, 2, function(x) {
    denom <- stats::IQR(x)
    if (!is.finite(denom) || denom == 0) return(0)
    m1 <- stats::median(x[1:xlen])
    m2 <- stats::median(x[(n - ylen + 1):n])
    abs(m1 - m2) / denom
  })
}

.flat_stage_indices <- function(ch, stage = "sample") {
  which(ch$samples$stage %in% stage & seq_along(ch$samples$stage) <= ch$samples$idx)
}

.flat_mat_subj_ll <- function(emc, stage = "sample") {
  c_n <- length(emc)
  idx_list <- lapply(emc, .flat_stage_indices, stage = stage)
  lens <- vapply(idx_list, length, integer(1))
  iter_n <- min(lens)
  if (!is.finite(iter_n) || iter_n <= 0) return(matrix(numeric(0), nrow = 0, ncol = 0))

  sll1 <- emc[[1]]$samples$subj_ll
  s_n <- nrow(sll1)
  out <- array(NA_real_, dim = c(c_n, iter_n, s_n))
  for (i in seq_len(c_n)) {
    idx <- tail(idx_list[[i]], iter_n)
    out[i, , ] <- t(emc[[i]]$samples$subj_ll[, idx, drop = FALSE])
  }
  mat <- matrix(out, ncol = s_n)
  colnames(mat) <- rownames(sll1)
  mat
}

.flat_mat_alpha <- function(emc, stage = "sample") {
  c_n <- length(emc)
  idx_list <- lapply(emc, .flat_stage_indices, stage = stage)
  lens <- vapply(idx_list, length, integer(1))
  iter_n <- min(lens)
  if (!is.finite(iter_n) || iter_n <= 0) return(matrix(numeric(0), nrow = 0, ncol = 0))

  a1 <- emc[[1]]$samples$alpha
  p_n <- dim(a1)[1]
  s_n <- dim(a1)[2]
  out <- array(NA_real_, dim = c(c_n, iter_n, p_n * s_n))
  for (i in seq_len(c_n)) {
    idx <- tail(idx_list[[i]], iter_n)
    a <- emc[[i]]$samples$alpha[, , idx, drop = FALSE]
    a <- aperm(a, c(3, 1, 2))
    dim(a) <- c(iter_n, p_n * s_n)
    out[i, , ] <- a
  }
  mat <- matrix(out, ncol = p_n * s_n)
  p_names <- dimnames(a1)[[1]]
  s_names <- dimnames(a1)[[2]]
  if (!is.null(p_names) && !is.null(s_names)) {
    colnames(mat) <- as.vector(outer(p_names, s_names, paste, sep = "_"))
  }
  mat
}

.flat_mat_theta_mu <- function(emc, stage = "sample") {
  if (is.null(emc[[1]]$samples$theta_mu)) return(NULL)

  c_n <- length(emc)
  idx_list <- lapply(emc, .flat_stage_indices, stage = stage)
  lens <- vapply(idx_list, length, integer(1))
  iter_n <- min(lens)
  if (!is.finite(iter_n) || iter_n <= 0) return(matrix(numeric(0), nrow = 0, ncol = 0))

  mu1 <- emc[[1]]$samples$theta_mu
  p_n <- nrow(mu1)
  out <- array(NA_real_, dim = c(c_n, iter_n, p_n))
  for (i in seq_len(c_n)) {
    idx <- tail(idx_list[[i]], iter_n)
    out[i, , ] <- t(emc[[i]]$samples$theta_mu[, idx, drop = FALSE])
  }
  mat <- matrix(out, ncol = p_n)
  colnames(mat) <- rownames(mu1)
  mat
}

.flat_mat_theta_var <- function(emc, stage = "sample") {
  if (is.null(emc[[1]]$samples$theta_var)) return(NULL)

  c_n <- length(emc)
  idx_list <- lapply(emc, .flat_stage_indices, stage = stage)
  lens <- vapply(idx_list, length, integer(1))
  iter_n <- min(lens)
  if (!is.finite(iter_n) || iter_n <= 0) return(matrix(numeric(0), nrow = 0, ncol = 0))

  var1 <- emc[[1]]$samples$theta_var
  p_n <- dim(var1)[1]
  idx_ut <- which(upper.tri(matrix(NA_real_, p_n, p_n), diag = TRUE), arr.ind = TRUE)
  k_n <- nrow(idx_ut)
  out <- array(NA_real_, dim = c(c_n, iter_n, k_n))
  for (i in seq_len(c_n)) {
    idx <- tail(idx_list[[i]], iter_n)
    var_i <- emc[[i]]$samples$theta_var[, , idx, drop = FALSE]
    for (k in seq_len(k_n)) {
      out[i, , k] <- var_i[idx_ut[k, 1], idx_ut[k, 2], ]
    }
  }
  mat <- matrix(out, ncol = k_n)
  p_names <- dimnames(var1)[[1]]
  if (!is.null(p_names)) {
    colnames(mat) <- paste0("cov_", p_names[idx_ut[, 1]], "_", p_names[idx_ut[, 2]])
  }
  mat
}

.get_flat_matrix <- function(emc, selection, stage = "sample") {
  if (selection == "alpha") return(.flat_mat_alpha(emc, stage = stage))
  if (selection == "subj_ll") return(.flat_mat_subj_ll(emc, stage = stage))
  if (selection == "theta_mu") return(.flat_mat_theta_mu(emc, stage = stage))
  if (selection == "theta_var") return(.flat_mat_theta_var(emc, stage = stage))
  stop("Unknown flatness selection: ", selection)
}

.check_flatness_stop <- function(emc, max_flat_loc, flat_selection, flat_p1 = 1/3, flat_p2 = 1/3,
                                 stage = "sample") {
  stats_by_sel <- list()
  max_by_sel <- numeric(0)

  for (sel in flat_selection) {
    mat <- .get_flat_matrix(emc, sel, stage = stage)
    if (is.null(mat)) next
    cur <- .flat_loc_stats_per_col(mat, p1 = flat_p1, p2 = flat_p2)
    if (length(cur) == 0) next
    stats_by_sel[[sel]] <- cur
    max_by_sel[sel] <- max(cur, na.rm = TRUE)
  }

  if (length(max_by_sel) == 0) {
    return(list(
      flat_done = FALSE,
      max_flat = Inf,
      max_by_selection = max_by_sel,
      stats_by_selection = stats_by_sel
    ))
  }

  max_flat <- max(max_by_sel, na.rm = TRUE)
  list(
    flat_done = is.finite(max_flat) && (max_flat <= max_flat_loc),
    max_flat = max_flat,
    max_by_selection = max_by_sel,
    stats_by_selection = stats_by_sel
  )
}

.trim_sample_stage <- function(emc, max_sample_iter) {
  cur_sample <- chain_n(emc)[, "sample"][1]
  if (cur_sample <= max_sample_iter) return(emc)
  n_drop <- cur_sample - max_sample_iter
  subset(emc, stage = "sample", filter = n_drop, keep_stages = TRUE)
}
