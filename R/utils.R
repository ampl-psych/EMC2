#' Print EMC2 build configuration
#'
#' Prints the compiler flags and options that were active when the package
#' was installed. Useful for debugging performance or numerical issues.
#'
#' @export
emc2_build_info <- function() invisible(.emc2_build_info())


# From Zach's branch
.floor_to_rt_resolution <- function(x, rt_resolution=1/60) {
  if (is.null(rt_resolution)) return(x)
  finite <- is.finite(x)
  if (!any(finite)) return(x)

  # Avoid dropping exact-bin values one bin lower due to binary rounding, e.g.
  # floor(0.7 / 0.05) == 13 instead of 14 on typical platforms.
  q <- x[finite] / rt_resolution
  tol <- 16 * .Machine$double.eps * pmax(1, abs(q))
  x[finite] <- floor(q + tol) * rt_resolution
  x
}

.emc2_ll_cache_version <- 2L

.is_valid_ll_cache <- function(dadm, n_trials, n_lR, has_RACE_col) {
  if (!is.data.frame(dadm)) return(FALSE)
  if (!identical(attr(dadm, "emc2_ll_cache_version"), .emc2_ll_cache_version)) return(FALSE)
  if (is.null(attr(dadm, "emc2_all_finite_trials"))) return(FALSE)
  all_finite <- attr(dadm, "emc2_all_finite_trials")
  if (!is.logical(all_finite) || length(all_finite) != 1L || is.na(all_finite)) return(FALSE)

  if (has_RACE_col) {
    race_nacc <- attr(dadm, "RACE_nacc_by_row")
    race_mask <- attr(dadm, "RACE_mask")
    if (is.null(race_nacc) || is.null(race_mask)) return(FALSE)
    if (length(race_nacc) != n_trials || length(race_mask) != n_trials) return(FALSE)
  }

  if (isTRUE(all_finite)) return(TRUE)
  if (n_lR <= 0L || (n_trials %% n_lR) != 0L) return(FALSE)
  n_unique_trials <- n_trials %/% n_lR

  finite_rt_mask <- attr(dadm, "finite_rt_mask")
  finite_rt_unique <- attr(dadm, "finite_rt_unique_trial_indices")
  other_unique <- attr(dadm, "other_unique_trial_indices")
  if (is.null(finite_rt_mask) || is.null(finite_rt_unique) || is.null(other_unique)) return(FALSE)
  if (length(finite_rt_mask) != n_trials) return(FALSE)
  if ((length(finite_rt_unique) + length(other_unique)) != n_unique_trials) return(FALSE)

  valid_idx <- function(x) {
    if (!is.numeric(x) && !is.integer(x)) return(FALSE)
    if (length(x) == 0L) return(TRUE)
    all(is.finite(x)) && all(x == as.integer(x)) &&
      all(x >= 0L) && all(x <= (n_unique_trials - 1L))
  }
  if (!valid_idx(finite_rt_unique) || !valid_idx(other_unique)) return(FALSE)
  finite_rt_unique <- as.integer(finite_rt_unique)
  other_unique <- as.integer(other_unique)
  if (length(unique(finite_rt_unique)) != length(finite_rt_unique)) return(FALSE)
  if (length(unique(other_unique)) != length(other_unique)) return(FALSE)
  if (length(intersect(finite_rt_unique, other_unique)) > 0L) return(FALSE)
  if (!setequal(c(finite_rt_unique, other_unique), 0:(n_unique_trials - 1L))) return(FALSE)

  # Validate cache partition against current dadm content to guard against stale attributes.
  rts <- dadm[["rt"]]
  R_idx <- dadm[["R"]]
  if (is.null(rts) || is.null(R_idx)) return(FALSE)
  race_nacc <- if (has_RACE_col) attr(dadm, "RACE_nacc_by_row") else NULL
  if (has_RACE_col && (is.null(race_nacc) || length(race_nacc) != n_trials)) return(FALSE)
  finite_set <- finite_rt_unique
  other_set <- other_unique
  for (j0 in 0:(n_unique_trials - 1L)) {
    start_row <- 1L + j0 * n_lR
    n_lR_j <- if (has_RACE_col) race_nacc[start_row] else n_lR
    if (!is.finite(n_lR_j) || n_lR_j < 1L || n_lR_j > n_lR) return(FALSE)
    finite_trial <- is.finite(rts[start_row]) && rts[start_row] > 0 && !is.na(R_idx[start_row])
    if (finite_trial) {
      if (!(j0 %in% finite_set)) return(FALSE)
      rows <- start_row:(start_row + n_lR_j - 1L)
      if (!all(finite_rt_mask[rows])) return(FALSE)
    } else {
      if (!(j0 %in% other_set)) return(FALSE)
    }
  }
  TRUE
}

.cache_ll_data_attrs <- function(dadm, force_rebuild = FALSE) {
  if (!is.data.frame(dadm)) return(dadm)

  cols <- names(dadm)
  if (!all(c("lR", "rt", "R") %in% cols)) {
    attr(dadm, "emc2_ll_cache_version") <- .emc2_ll_cache_version
    return(dadm)
  }

  n_trials <- nrow(dadm)
  lR <- dadm[["lR"]]
  lR_codes <- as.integer(lR)
  n_lR <- length(unique(lR_codes))

  has_RACE_col <- "RACE" %in% cols && is.factor(dadm[["RACE"]])
  if (!isTRUE(force_rebuild) &&
      .is_valid_ll_cache(dadm, n_trials = n_trials, n_lR = n_lR, has_RACE_col = has_RACE_col)) {
    return(dadm)
  }

  if (isTRUE(force_rebuild)) {
    attr(dadm, "emc2_all_finite_trials") <- NULL
    attr(dadm, "finite_rt_mask") <- NULL
    attr(dadm, "finite_rt_unique_trial_indices") <- NULL
    attr(dadm, "other_unique_trial_indices") <- NULL
    attr(dadm, "RACE_nacc_by_row") <- NULL
    attr(dadm, "RACE_mask") <- NULL
  }
  if (has_RACE_col) {
    race_idx <- dadm[["RACE"]]
    race_levels <- levels(race_idx)
    nacc_by_level <- suppressWarnings(as.integer(race_levels))
    if (anyNA(nacc_by_level)) stop("RACE column levels must be integer-valued (e.g., '2', '3').")

    race_nacc_by_row <- rep.int(as.integer(n_lR), n_trials)
    race_mask <- rep.int(TRUE, n_trials)
    race_codes <- as.integer(race_idx)
    for (i in seq_len(n_trials)) {
      code <- race_codes[i]
      if (is.na(code)) next
      nacc <- nacc_by_level[code]
      race_nacc_by_row[i] <- nacc
      lR_i <- lR_codes[i]
      if (!is.na(lR_i) && lR_i > nacc) race_mask[i] <- FALSE
    }
    attr(dadm, "RACE_nacc_by_row") <- race_nacc_by_row
    attr(dadm, "RACE_mask") <- race_mask
  }

  all_finite_trials <- TRUE
  if (n_trials > 0L) {
    if (n_lR <= 0L || (n_trials %% n_lR) != 0L) {
      all_finite_trials <- FALSE
    } else {
      start_idx <- seq.int(1L, n_trials, by = n_lR)
      rts <- dadm[["rt"]][start_idx]
      R_idx <- dadm[["R"]][start_idx]
      ok <- is.finite(rts) & rts > 0 & !is.na(R_idx)
      LT <- if ("LT" %in% cols) dadm[["LT"]][start_idx] else rep(0, length(start_idx))
      UT <- if ("UT" %in% cols) dadm[["UT"]][start_idx] else rep(Inf, length(start_idx))
      ok <- ok & (LT == 0) & is.infinite(UT)
      if (!all(ok)) all_finite_trials <- FALSE
    }
  }
  attr(dadm, "emc2_all_finite_trials") <- all_finite_trials

  if (!all_finite_trials && n_trials > 0L && n_lR > 0L && (n_trials %% n_lR) == 0L) {
    n_unique_trials <- n_trials %/% n_lR
    finite_rt_mask <- rep.int(FALSE, n_trials)
    finite_rt_unique_trial_indices <- integer(0)
    other_unique_trial_indices <- integer(0)

    start_idx <- seq.int(1L, n_trials, by = n_lR)
    rts_dadm <- dadm[["rt"]]
    R_idxs_dadm <- dadm[["R"]]
    race_nacc_by_row <- if (has_RACE_col) attr(dadm, "RACE_nacc_by_row") else NULL

    for (j0 in 0:(n_unique_trials - 1L)) {
      start_row_idx <- start_idx[j0 + 1L]
      rt_j <- rts_dadm[start_row_idx]
      R_j <- R_idxs_dadm[start_row_idx]
      n_lR_j <- if (has_RACE_col && length(race_nacc_by_row) == n_trials) {
        race_nacc_by_row[start_row_idx]
      } else {
        n_lR
      }
      if (is.finite(rt_j) && rt_j > 0 && !is.na(R_j)) {
        finite_rt_unique_trial_indices <- c(finite_rt_unique_trial_indices, j0)
        finite_rt_mask[start_row_idx + 0:(n_lR_j - 1L)] <- TRUE
      } else {
        other_unique_trial_indices <- c(other_unique_trial_indices, j0)
      }
    }

    attr(dadm, "finite_rt_mask") <- finite_rt_mask
    attr(dadm, "finite_rt_unique_trial_indices") <- finite_rt_unique_trial_indices
    attr(dadm, "other_unique_trial_indices") <- other_unique_trial_indices
  }

  attr(dadm, "emc2_ll_cache_version") <- .emc2_ll_cache_version
  dadm
}

# # Last observation carried forward
# # Replaces NA values with the last non-NA value
# # Mimics zoo::na.locf behavior (vectorized for speed)
# na_locf <- function(x, na.rm = FALSE) {
#   if (length(x) == 0) return(x)
#
#   # Check if all values are NA
#   na_mask <- is.na(x)
#   if (all(na_mask)) {
#     if (na.rm) {
#       return(x[0])  # Return empty vector
#     } else {
#       return(x)  # Return as is
#     }
#   }
#
#   # Vectorized approach: create an index that tracks the last non-NA position
#   # For each position, we need the index of the last non-NA value up to that point
#   idx <- seq_along(x)
#   idx[na_mask] <- NA  # Set NA positions to NA in index
#   idx <- cummax(ifelse(na_mask, 0, idx))  # Cumulative max gives us last non-NA index
#
#   # Replace values: use the index to look up the last non-NA value
#   # For positions where idx is 0 (leading NAs), keep as NA
#   result <- x
#   non_zero <- idx > 0
#   result[non_zero] <- x[idx[non_zero]]
#
#   # If na.rm=TRUE, remove leading NAs
#   if (na.rm) {
#     first_non_na <- which(!na_mask)[1]
#     if (!is.na(first_non_na)) {
#       result <- result[first_non_na:length(result)]
#     }
#   }
#
#   return(result)
# }

#
# augment <- function(s,da,design)
#   # Adds attributes to augmented data
#   # learn: empty array for Q values with dim = choice alternative (low,high) x
#   #   stimulus x trials (max across stimuli)
#   # index: look up (row number) for stimuli in da, a matrix dim  = max trials x
#   #   stimulus matrix (rows for each choice alternative contiguous)
# {
#   if (!is.null(design$adapt$stimulus)) {
#     targets <- design$adapt$stimulus$targets
#     par <- design$adapt$stimulus$output_name
#     maxn <- max(sapply(dimnames(targets)[[1]],function(x){table(da[da$subjects==s,x])}))
#     # da index x stimulus
#     out <- sapply(targets[1,],getIndex,cname=dimnames(targets)[[1]][1],
#                   da=da[da$subjects==s,],maxn=maxn)
#     stimulus <- list(index=out)
#     # accumulator x stimulus x trials
#     stimulus$learn <- array(NA,dim=c(dim(targets),maxn/dim(targets)[1]),
#                             dimnames=list(rownames(targets),targets[1,],NULL))
#     stimulus$targets <- targets
#     stimulus$par <- par
#   } # add other types here
#   list(stimulus=stimulus)
# }
#
# getIndex <- function(typei,cname,da,maxn) {
#   out <- which(da[,cname]==typei)
#   c(out,rep(NA,maxn-length(out)))
# }
