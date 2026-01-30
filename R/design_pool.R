default_memory_keep_cols <- function() {
  c("subjects", "trials", "rt", "R", "lR", "winner", "lM", "RACE", "TIMEOUT", "Rgo")
}

resolve_memory_keep_cols <- function(dadm, model, keep_cols = NULL) {
  keep <- default_memory_keep_cols()
  if (!is.null(model) && !is.function(model) && !is.null(model$trend)) {
    trend <- model$trend
    covs <- unique(unlist(lapply(trend, function(x) x$covariate)))
    ats <- unique(unlist(lapply(trend, function(x) x$at)))
    covs <- covs[!is.na(covs)]
    ats <- ats[!is.na(ats)]
    keep <- c(keep, covs, ats)
  }
  if (!is.null(keep_cols)) {
    keep <- c(keep, keep_cols)
  }
  keep <- unique(keep)
  intersect(names(dadm), keep)
}

has_trend_map <- function(model) {
  if (is.null(model) || is.function(model) || is.null(model$trend)) return(FALSE)
  any(vapply(model$trend, function(x) !is.null(x$map) && length(x$map) > 0, logical(1)))
}

build_design_pool <- function(dadm, designs) {
  if (is.null(designs)) return(NULL)
  n_rows <- nrow(dadm)
  pool_cols <- list()
  pool_names <- character()
  pool_map <- vector("list", length(designs))
  names(pool_map) <- names(designs)
  sig_env <- new.env(hash = TRUE, parent = emptyenv())

  fmt_num <- function(x) formatC(x, digits = 17, format = "g")

  for (i in seq_along(designs)) {
    dm <- designs[[i]]
    if (is.null(dm)) {
      pool_map[[i]] <- list(idx = integer(), colnames = character())
      next
    }
    exp <- attr(dm, "expand")
    n_cols <- ncol(dm)
    col_names <- colnames(dm)
    idx <- integer(n_cols)
    for (j in seq_len(n_cols)) {
      col_vec <- dm[, j]
      if (!is.null(exp)) col_vec <- col_vec[exp]
      na_count <- sum(is.na(col_vec))
      if (na_count == length(col_vec)) {
        sig <- paste0("allna|", length(col_vec))
      } else {
        s <- sum(col_vec, na.rm = TRUE)
        ss <- sum(col_vec * col_vec, na.rm = TRUE)
        mn <- min(col_vec, na.rm = TRUE)
        mx <- max(col_vec, na.rm = TRUE)
        sig <- paste(length(col_vec), na_count, fmt_num(s), fmt_num(ss),
                     fmt_num(mn), fmt_num(mx), sep = "|")
      }
      cand <- sig_env[[sig]]
      if (is.null(cand)) {
        pool_cols[[length(pool_cols) + 1]] <- col_vec
        pool_names[length(pool_cols)] <- col_names[j]
        idx[j] <- length(pool_cols)
        sig_env[[sig]] <- idx[j]
      } else {
        cand_idx <- cand
        found <- FALSE
        for (ci in cand_idx) {
          if (identical(col_vec, pool_cols[[ci]])) {
            idx[j] <- ci
            found <- TRUE
            break
          }
        }
        if (!found) {
          pool_cols[[length(pool_cols) + 1]] <- col_vec
          pool_names[length(pool_cols)] <- col_names[j]
          idx[j] <- length(pool_cols)
          sig_env[[sig]] <- c(cand_idx, idx[j])
        }
      }
    }
    pool_map[[i]] <- list(idx = idx, colnames = col_names)
  }
  if (length(pool_cols) == 0) {
    pool <- matrix(numeric(0), nrow = n_rows, ncol = 0)
  } else {
    pool <- do.call(cbind, pool_cols)
    colnames(pool) <- pool_names
  }
  return(list(pool = pool, map = pool_map))
}

get_designs_expanded <- function(dadm, model = NULL) {
  pool <- attr(dadm, "design_pool")
  pool_map <- attr(dadm, "design_pool_map")
  if (!is.null(pool) && !is.null(pool_map)) {
    if (is.null(model)) model <- attr(dadm, "model")
    if (is.function(model)) model <- model()
    p_types <- names(model$p_types)
    designs <- vector("list", length(p_types))
    names(designs) <- p_types
    for (i in seq_along(p_types)) {
      map_i <- pool_map[[p_types[i]]]
      if (is.null(map_i)) {
        designs[[i]] <- matrix(nrow = nrow(dadm), ncol = 0)
        next
      }
      idx <- map_i$idx
      if (length(idx) == 0) {
        mat <- matrix(nrow = nrow(dadm), ncol = 0)
      } else {
        mat <- pool[, idx, drop = FALSE]
        colnames(mat) <- map_i$colnames
      }
      designs[[i]] <- mat
    }
    return(designs)
  }
  designs <- attr(dadm, "designs")
  if (is.null(designs)) return(NULL)
  if (is.numeric(designs)) stop("designs reference not resolved for pooled design lookup")
  if (!is.null(model)) {
    if (is.function(model)) model <- model()
    p_types <- names(model$p_types)
    if (!is.null(names(designs)) && all(p_types %in% names(designs))) {
      out <- lapply(p_types, function(p) {
        x <- designs[[p]]
        if (is.null(x)) return(NULL)
        exp <- attr(x, "expand")
        if (is.null(exp)) return(x)
        x[exp, , drop = FALSE]
      })
      names(out) <- p_types
      return(out)
    }
  }
  out <- lapply(designs, function(x) {
    if (is.null(x)) return(NULL)
    exp <- attr(x, "expand")
    if (is.null(exp)) return(x)
    x[exp, , drop = FALSE]
  })
  names(out) <- names(designs)
  return(out)
}
