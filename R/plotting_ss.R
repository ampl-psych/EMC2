get_stop_signal_postn_quantiles <- function(postn_list, level, value_col, quants) {
  # Posterior/prior draws can contain different SSD bins after simulation because
  # response presence can vary by draw. Align draws by x position before taking
  # quantiles so missing bins do not shift summaries into the wrong column.
  draw_summaries <- lapply(postn_list, function(draw) draw[[level]])
  valid_summaries <- Filter(function(x) {
    !is.null(x) && nrow(x) > 0 && all(c("x_plot", value_col) %in% names(x))
  }, draw_summaries)

  if (!length(valid_summaries)) return(NULL)

  x_plot <- sort(unique(unlist(lapply(valid_summaries, function(x) x$x_plot))))
  ssd_labels <- NULL
  if (all(vapply(valid_summaries, function(x) "ssd" %in% names(x), logical(1)))) {
    ssd_lookup <- do.call(rbind, lapply(valid_summaries, function(x) {
      unique(x[, c("x_plot", "ssd"), drop = FALSE])
    }))
    ssd_lookup <- ssd_lookup[!duplicated(ssd_lookup$x_plot), , drop = FALSE]
    ssd_labels <- ssd_lookup$ssd[match(x_plot, ssd_lookup$x_plot)]
  }
  probs <- sort(c(quants, 0.5))

  y_mat <- do.call(cbind, lapply(draw_summaries, function(x) {
    y <- rep(NA_real_, length(x_plot))
    if (!is.null(x) && nrow(x) > 0 && all(c("x_plot", value_col) %in% names(x))) {
      y[match(x$x_plot, x_plot)] <- x[[value_col]]
    }
    y
  }))

  qy <- apply(y_mat, 1, function(y) {
    if (all(is.na(y))) {
      rep(NA_real_, length(probs))
    } else {
      quantile(y, probs = probs, na.rm = TRUE)
    }
  })

  out <- rbind(qy, x_plot)
  if (!is.null(ssd_labels)) colnames(out) <- ssd_labels
  out
}

round_ssd_values <- function(SSD, ssd_round = NULL) {
  if (is.null(ssd_round)) return(SSD)
  round(SSD / ssd_round) * ssd_round
}

get_stop_signal_bin_ssd <- function(df, ssd_round = NULL) {
  round_ssd_values(df$SSD, ssd_round)
}

format_ssd_values <- function(x) {
  format(round(x, 4), trim = TRUE, scientific = FALSE)
}

format_stop_signal_quantile_intervals <- function(probs) {
  up <- round(probs[-1] * 100)
  low <- round(probs[-length(probs)] * 100)
  intervals <- paste0("(", paste(low, up, sep = ","), "]")
  intervals[1] <- sub("^\\(", "[", intervals[1])
  intervals
}

message_ignored_value_binning_args <- function() {
  message("`probs`, `use_global_quantiles`, and `on_duplicate_quantiles` ",
          "are ignored when `ssd_binning = \"value\"`.")
}

validate_stop_signal_probs <- function(probs) {
  if (!is.numeric(probs) || length(probs) < 2 ||
      anyNA(probs) || any(!is.finite(probs))) {
    stop("`probs` must be a numeric vector with at least two finite values.")
  }
  if (any(probs < 0 | probs > 1)) {
    stop("`probs` must contain values between 0 and 1.")
  }
  if (any(duplicated(probs))) {
    stop("`probs` must not contain duplicate values.")
  }
  if (is.unsorted(probs, strictly = TRUE)) {
    stop("`probs` must be sorted in strictly increasing order.")
  }
  if (!isTRUE(all.equal(probs[1], 0)) || !isTRUE(all.equal(probs[length(probs)], 1))) {
    stop("`probs` must start at 0 and end at 1.")
  }

  invisible(NULL)
}

validate_stop_signal_quants <- function(quants) {
  if (!is.numeric(quants) || length(quants) != 2 ||
      anyNA(quants) || any(!is.finite(quants))) {
    stop("`quants` must be a numeric vector of exactly two finite values.")
  }
  if (any(quants < 0 | quants > 1)) {
    stop("`quants` must contain values between 0 and 1.")
  }
  if (quants[1] >= quants[2]) {
    stop("`quants` must be in strictly increasing order.")
  }
  if (!(quants[1] < 0.5 && quants[2] > 0.5)) {
    stop("`quants` must bracket 0.5.")
  }

  invisible(NULL)
}

validate_stop_signal_sources <- function(to_plot, use_lim) {
  valid_sources <- c("data", "posterior", "prior")

  if (!is.character(to_plot) || !length(to_plot) || anyNA(to_plot)) {
    stop("`to_plot` must be a non-empty character vector.")
  }
  if (!all(to_plot %in% valid_sources)) {
    stop("`to_plot` must contain only \"data\", \"posterior\", and/or \"prior\".")
  }
  if (any(duplicated(to_plot))) {
    stop("`to_plot` must not contain duplicate values.")
  }

  if (!is.character(use_lim) || !length(use_lim) || anyNA(use_lim)) {
    stop("`use_lim` must be a non-empty character vector.")
  }
  if (!all(use_lim %in% valid_sources)) {
    stop("`use_lim` must contain only \"data\", \"posterior\", and/or \"prior\".")
  }
  if (any(duplicated(use_lim))) {
    stop("`use_lim` must not contain duplicate values.")
  }
  if (!any(use_lim %in% to_plot)) {
    stop("`use_lim` must include at least one source that is also in `to_plot`.")
  }

  invisible(NULL)
}

can_generate_stop_signal_predictives <- function(input) {
  if (inherits(input, "emc")) return(TRUE)
  if (is.data.frame(input) || !is.list(input) || !length(input)) return(FALSE)
  all(vapply(input, inherits, logical(1), "emc"))
}

validate_stop_signal_predictive_request <- function(input, post_predict, prior_predict,
                                                    to_plot, to_plot_missing) {
  if (to_plot_missing) return(invisible(NULL))

  can_generate <- can_generate_stop_signal_predictives(input)

  if ("posterior" %in% to_plot && is.null(post_predict) && !can_generate) {
    stop("`post_predict` must be supplied when `to_plot` includes \"posterior\" ",
         "and `input` is not an emc object.")
  }
  if ("prior" %in% to_plot && is.null(prior_predict) && !can_generate) {
    stop("`prior_predict` must be supplied when `to_plot` includes \"prior\" ",
         "and `input` is not an emc object.")
  }

  invisible(NULL)
}

normalize_stop_signal_layout <- function(layout) {
  if (is.null(layout)) return(NA)
  if (length(layout) == 1 && is.na(layout)) return(layout)
  if (!is.numeric(layout) || length(layout) != 2 ||
      anyNA(layout) || any(!is.finite(layout)) ||
      any(layout < 1) || any(layout != as.integer(layout))) {
    stop("`layout` must be NA, NULL, or a length-2 positive integer vector.")
  }

  as.integer(layout)
}

normalize_stop_signal_legendpos <- function(legendpos) {
  if (is.null(legendpos)) return(c(NA_character_, NA_character_))
  if (!is.atomic(legendpos) || length(legendpos) < 1 || length(legendpos) > 2) {
    stop("`legendpos` must be NULL, NA, or a character vector of length 1 or 2.")
  }
  if (all(is.na(legendpos))) return(c(NA_character_, NA_character_))
  if (!is.character(legendpos)) {
    stop("`legendpos` must be NULL, NA, or a character vector of length 1 or 2.")
  }
  if (length(legendpos) == 1) return(c(legendpos, NA_character_))

  legendpos
}

draw_stop_signal_se <- function(x, y, se) {
  # Base graphics warns on zero-height arrows; skip those while keeping all
  # finite, nonzero SE intervals.
  x <- as.numeric(x)
  y <- as.numeric(y)
  se <- as.numeric(se)
  n <- min(length(x), length(y), length(se))
  if (!n) return(invisible(NULL))

  x <- x[seq_len(n)]
  y <- y[seq_len(n)]
  se <- se[seq_len(n)]
  y0 <- y - se
  y1 <- y + se
  keep <- is.finite(x) & is.finite(y0) & is.finite(y1) &
    abs(y1 - y0) > sqrt(.Machine$double.eps)
  if (!any(keep)) return(invisible(NULL))

  arrows(x0 = x[keep], y0 = y0[keep],
         x1 = x[keep], y1 = y1[keep],
         angle = 90, code = 3, length = 0.05)
  invisible(NULL)
}

get_stop_signal_x_range_values <- function(x) {
  x <- as.numeric(x)
  x[is.finite(x)]
}

get_stop_signal_y_range_values <- function(df, value_col) {
  if (is.null(df[[value_col]])) return(numeric(0))

  y <- df[[value_col]]
  if (!is.null(df$se)) {
    y <- c(y, df[[value_col]] - df$se, df[[value_col]] + df$se)
  }

  y[is.finite(y)]
}

draw_stop_signal_predictive_points <- function(postn_list, level, value_col,
                                               src_args, point_args = list()) {
  if (is.null(postn_list) || !length(postn_list)) return(invisible(NULL))

  # Draw-level points make the empirical distribution behind the predictive
  # ribbon visible. They are deliberately opt-in because dense posterior draws
  # can otherwise dominate the plot.
  point_args <- add_defaults(point_args,
                             pch = 16,
                             cex = 0.35,
                             alpha.f = 0.35,
                             jitter = 0)
  jitter <- point_args$jitter
  if (!is.numeric(jitter) || length(jitter) != 1 || is.na(jitter) || jitter < 0) {
    stop("`predictive_point_args$jitter` must be a single non-negative numeric value.")
  }

  pts <- do.call(rbind, lapply(postn_list, function(draw) {
    df <- draw[[level]]
    if (is.null(df) || !all(c("x_plot", value_col) %in% names(df))) return(NULL)
    data.frame(x = as.numeric(df$x_plot),
               y = as.numeric(df[[value_col]]))
  }))
  if (is.null(pts) || !nrow(pts)) return(invisible(NULL))

  pts <- pts[is.finite(pts$x) & is.finite(pts$y), , drop = FALSE]
  if (!nrow(pts)) return(invisible(NULL))

  if (jitter > 0) {
    # Deterministic jitter spreads points that share an SSD/bin x-position
    # without adding random variation to reproduced plots.
    point_idx <- split(seq_len(nrow(pts)), pts$x)
    for (idx in point_idx) {
      if (length(idx) > 1) {
        pts$x[idx] <- pts$x[idx] + seq(-jitter, jitter, length.out = length(idx))
      }
    }
  }

  point_args$col <- do.call(adjustcolor,
                            fix_dots(add_defaults(point_args,
                                                  col = src_args$col[1]),
                                     adjustcolor))
  point_args$alpha.f <- NULL
  point_args$jitter <- NULL
  point_args <- fix_dots_plot(point_args)
  do.call(points, c(list(x = pts$x, y = pts$y), point_args))

  invisible(NULL)
}

get_stop_signal_predictive_count_summary <- function(postn_list, level, value_col, x_plot) {
  # Count how many posterior/prior draws contribute a finite plotted summary at
  # this x-position. For SRRT, this is useful because high SSD bins can have few
  # or no signal responses in some predictive draws.
  y_vals <- numeric(0)
  n_vals <- numeric(0)
  n_obs_vals <- numeric(0)

  for (draw in postn_list) {
    df <- draw[[level]]
    if (is.null(df) || !all(c("x_plot", value_col) %in% names(df))) next
    match_idx <- match(x_plot, df$x_plot)
    if (is.na(match_idx)) next

    y <- as.numeric(df[[value_col]][match_idx])
    if (!is.finite(y)) next

    y_vals <- c(y_vals, y)
    if ("n" %in% names(df)) {
      n <- as.numeric(df$n[match_idx])
      if (is.finite(n)) n_vals <- c(n_vals, n)
    }
    if ("n_obs" %in% names(df)) {
      n_obs <- as.numeric(df$n_obs[match_idx])
      if (is.finite(n_obs)) n_obs_vals <- c(n_obs_vals, n_obs)
    }
  }

  out <- data.frame(n_draws = length(y_vals),
                    n_min = NA_real_,
                    n_median = NA_real_,
                    n_max = NA_real_,
                    n_obs_min = NA_real_,
                    n_obs_median = NA_real_,
                    n_obs_max = NA_real_)
  if (length(n_vals)) {
    out$n_min <- min(n_vals, na.rm = TRUE)
    out$n_median <- median(n_vals, na.rm = TRUE)
    out$n_max <- max(n_vals, na.rm = TRUE)
  }
  if (length(n_obs_vals)) {
    out$n_obs_min <- min(n_obs_vals, na.rm = TRUE)
    out$n_obs_median <- median(n_obs_vals, na.rm = TRUE)
    out$n_obs_max <- max(n_obs_vals, na.rm = TRUE)
  }

  out
}

stop_signal_rbind_fill <- function(x) {
  # Combine observed and predictive summaries even though they have different
  # columns: observed rows have value/se, predictive rows have lower/median/upper.
  all_names <- unique(unlist(lapply(x, names)))
  x <- lapply(x, function(df) {
    missing_names <- setdiff(all_names, names(df))
    for (nm in missing_names) df[[nm]] <- NA
    df[all_names]
  })

  do.call(rbind, x)
}

drop_all_na_columns <- function(x) {
  # Printing diagnostics should be compact, but the returned object keeps the
  # full schema for programmatic consistency across plot types.
  if (!nrow(x)) return(x)
  x[, !vapply(x, function(col) all(is.na(col)), logical(1)), drop = FALSE]
}

print_stop_signal_plot_data <- function(x) {
  # Print observed and predictive rows separately so columns that are meaningful
  # for only one row type do not appear as long runs of NA in the other.
  predictive_rows <- "median" %in% names(x) & !is.na(x$median)
  observed_rows <- "value" %in% names(x) & !is.na(x$value)

  if (any(predictive_rows)) {
    cat("Predictive plotted summaries:\n")
    print(drop_all_na_columns(x[predictive_rows, , drop = FALSE]), row.names = FALSE)
  }
  if (any(observed_rows)) {
    if (any(predictive_rows)) cat("\n")
    cat("Observed plotted summaries:\n")
    print(drop_all_na_columns(x[observed_rows, , drop = FALSE]), row.names = FALSE)
  }
  if (!any(predictive_rows) && !any(observed_rows)) {
    print(drop_all_na_columns(x), row.names = FALSE)
  }

  invisible(NULL)
}

get_stop_signal_bin_label <- function(x_plot, binning, ssd = NULL) {
  # Participant-quantile plots do not have a single shared absolute SSD interval,
  # so the diagnostic table reports the quantile-bin label separately from SSD.
  if (identical(unname(binning), "individual_quantile")) {
    return(format_stop_signal_quantile_intervals(c(0, x_plot)))
  }
  if (!is.null(ssd)) return(as.character(ssd))

  NA_character_
}

get_stop_signal_plot_data <- function(summary_by_source, draw_quantiles_by_source, sources,
                                      source_bin_modes, value_col) {
  # Build the table that mirrors what is actually drawn. This is returned
  # invisibly and optionally printed for checking sparse bins or draw counts.
  out <- list()

  for (sname in names(sources)) {
    group_names <- names(summary_by_source[[sname]])
    if (is.null(group_names)) next

    for (group_key in group_names) {
      source_group <- summary_by_source[[sname]][[group_key]]
      if (is.null(source_group)) next

      if (!is.null(draw_quantiles_by_source[[sname]])) {
        quant_group <- draw_quantiles_by_source[[sname]][[group_key]]
        if (is.null(quant_group)) next

        for (level in names(quant_group)) {
          mat <- quant_group[[level]]
          if (is.null(mat)) next

          for (i in seq_len(ncol(mat))) {
            ssd <- if (!is.null(colnames(mat))) colnames(mat)[i] else NA_character_
            counts <- get_stop_signal_predictive_count_summary(
              source_group, level, value_col, unname(mat[nrow(mat), i])
            )
            out[[length(out) + 1]] <- cbind(
              data.frame(source = unname(sources[sname]),
                         dataset = sname,
                         binning = unname(source_bin_modes[sname]),
                         group_key = group_key,
                         within_level = level,
                         x_plot = unname(mat[nrow(mat), i]),
                         bin_label = get_stop_signal_bin_label(unname(mat[nrow(mat), i]),
                                                               source_bin_modes[sname],
                                                               ssd),
                         ssd = ssd,
                         lower = unname(mat[1, i]),
                         median = unname(mat[2, i]),
                         upper = unname(mat[3, i])),
              counts
            )
          }
        }
      } else {
        for (level in names(source_group)) {
          df <- source_group[[level]]
          if (is.null(df) || !all(c("x_plot", value_col) %in% names(df))) next

          out[[length(out) + 1]] <- data.frame(
            source = unname(sources[sname]),
            dataset = sname,
            binning = unname(source_bin_modes[sname]),
            group_key = group_key,
            within_level = level,
            x_plot = df$x_plot,
            bin_label = get_stop_signal_bin_label(df$x_plot,
                                                  source_bin_modes[sname],
                                                  if ("ssd" %in% names(df)) df$ssd else NULL),
            ssd = if ("ssd" %in% names(df)) as.character(df$ssd) else NA_character_,
            value = df[[value_col]],
            se = if ("se" %in% names(df)) df$se else NA_real_,
            n = if ("n" %in% names(df)) df$n else NA_real_,
            n_obs = if ("n_obs" %in% names(df)) df$n_obs else NA_real_
          )
        }
      }
    }
  }

  if (!length(out)) return(data.frame())
  stop_signal_rbind_fill(out)
}

has_duplicate_individual_ssd_quantiles <- function(df, probs, within_plot = NULL,
                                                   ssd_round = NULL) {
  # Individual quantile bins are computed separately within each panel,
  # within-plot level, and subject. If any subject has collapsed quantile
  # breaks, cut() cannot form the requested bins.
  groups <- split(df, df$group_key)

  any(vapply(groups, function(group_df) {
    within_groups <- if (is.null(within_plot)) {
      list(group_df)
    } else {
      split(group_df, group_df[[within_plot]])
    }

    any(vapply(within_groups, function(within_df) {
      # Posterior/prior predictives are summarized per draw. Check duplicate
      # quantile breaks at the same postn level where binning is later applied,
      # otherwise sparse draw-level SRRT data can be hidden by pooled draws.
      postn_groups <- if ("postn" %in% names(within_df)) {
        split(within_df, within_df$postn)
      } else {
        list(within_df)
      }

      any(vapply(postn_groups, function(postn_df) {
        subjects <- unique(postn_df$subjects)
        SSD_all <- get_stop_signal_bin_ssd(postn_df, ssd_round)
        any(vapply(subjects, function(s) {
          SSD <- SSD_all[postn_df$subjects == s & is.finite(SSD_all)]
          if (!length(SSD)) return(FALSE)
          any(duplicated(quantile(SSD, probs = probs, na.rm = TRUE)))
        }, logical(1)))
      }, logical(1)))
    }, logical(1)))
  }, logical(1)))
}

has_duplicate_global_ssd_quantiles <- function(df, probs, ssd_round = NULL) {
  # Global quantile bins pool SSDs within a plotted source. Duplicate breaks
  # indicate that the requested quantile grid is too fine for the SSD design.
  SSD <- get_stop_signal_bin_ssd(df, ssd_round)
  SSD <- SSD[is.finite(SSD)]
  if (!length(SSD)) return(FALSE)
  any(duplicated(quantile(SSD, probs = probs, na.rm = TRUE)))
}

get_reference_stop_signal_source_index <- function(data_sources, sources) {
  # Global quantile plots need one absolute SSD grid for all overlaid sources.
  # Prefer the observed data when it is plotted; otherwise use the first source
  # that contains finite SSDs (posterior-only and prior-only plots).
  candidate_idx <- seq_along(data_sources)
  if (!is.null(sources) && "data" %in% sources) {
    candidate_idx <- unique(c(which(sources == "data"), candidate_idx))
  }

  for (i in candidate_idx) {
    df <- data_sources[[i]]
    if (!is.null(df) && nrow(df) && !is.null(df$SSD)) {
      SSD <- df$SSD
      if (any(is.finite(SSD))) return(i)
    }
  }

  stop("No finite SSD values available for global quantile bins.")
}

get_reference_global_ssd_breaks <- function(data_sources, sources, probs, ssd_round = NULL) {
  reference_idx <- get_reference_stop_signal_source_index(data_sources, sources)
  SSD <- get_stop_signal_bin_ssd(data_sources[[reference_idx]], ssd_round)
  SSD <- SSD[is.finite(SSD)]

  quantile(SSD, probs = probs, na.rm = TRUE)
}

get_reference_global_ssd_values <- function(data_sources, sources, ssd_round = NULL) {
  reference_idx <- get_reference_stop_signal_source_index(data_sources, sources)
  SSD <- get_stop_signal_bin_ssd(data_sources[[reference_idx]], ssd_round)
  SSD <- SSD[is.finite(SSD)]

  sort(unique(SSD))
}

validate_global_ssd_breaks <- function(data_sources, sources, global_ssd_breaks,
                                       ssd_round = NULL) {
  if (any(duplicated(global_ssd_breaks))) return(FALSE)

  reference_ssd <- get_reference_global_ssd_values(data_sources, sources, ssd_round)
  all(vapply(data_sources, function(df) {
    if (is.null(df) || !nrow(df) || is.null(df$SSD)) return(TRUE)
    SSD <- get_stop_signal_bin_ssd(df, ssd_round)
    SSD <- SSD[is.finite(SSD)]
    if (!length(SSD)) return(TRUE)
    isTRUE(all.equal(sort(unique(SSD)), reference_ssd, tolerance = 1e-8))
  }, logical(1)))
}

has_duplicate_stop_signal_quantiles <- function(data_sources, sources, probs,
                                                use_global_quantiles, within_plot,
                                                ssd_round = NULL) {
  if (use_global_quantiles) {
    global_ssd_breaks <- get_reference_global_ssd_breaks(data_sources, sources, probs, ssd_round)
    return(!validate_global_ssd_breaks(data_sources, sources, global_ssd_breaks, ssd_round))
  }

  any(vapply(data_sources, function(df) {
    if (is.null(df) || !nrow(df) || is.null(df$SSD)) return(FALSE)
    has_duplicate_individual_ssd_quantiles(df, probs, within_plot, ssd_round)
  }, logical(1)))
}

get_reduced_stop_signal_probs <- function(data_sources, sources, probs, use_global_quantiles,
                                          within_plot, ssd_round = NULL) {
  # The reduce fallback preserves quantile binning but uses fewer evenly spaced
  # breaks. Search downward from the requested number of breaks and keep the
  # largest grid that is valid for every plotted source.
  for (n_breaks in seq(length(probs), 2)) {
    candidate_probs <- seq(0, 1, length.out = n_breaks)
    if (!has_duplicate_stop_signal_quantiles(data_sources, sources, candidate_probs,
                                             use_global_quantiles, within_plot,
                                             ssd_round)) {
      return(candidate_probs)
    }
  }

  stop("Duplicate quantile values detected even with the smallest quantile grid. ",
       "Please set `on_duplicate_quantiles = \"value\"`.")
}

draw_stop_signal_x_axis <- function(tick_data, bin_mode, global_title, participant_title,
                                    value_title = "SSD Value (sec.)",
                                    probs = NULL) {
  if (is.null(tick_data) || !"x_plot" %in% names(tick_data)) return(invisible(NULL))

  bin_mode <- unname(bin_mode)
  # Unknown modes are treated as individual-quantile axes first. The check below
  # then detects whether the supplied x values are actually SSD values.
  if (length(bin_mode) != 1 || is.na(bin_mode) ||
      !bin_mode %in% c("value", "global_quantile", "individual_quantile")) {
    bin_mode <- "individual_quantile"
  }

  if (identical(bin_mode, "individual_quantile") && !is.null(probs)) {
    # Individual quantile axes should be at probs[-1]. If the data reaching the
    # axis has absolute SSD positions instead, draw a value axis to avoid
    # misleading quantile-bin labels.
    expected_x <- probs[-1]
    if (length(tick_data$x_plot) != length(expected_x) ||
        !isTRUE(all.equal(unname(tick_data$x_plot), unname(expected_x), tolerance = 1e-8))) {
      bin_mode <- "value"
    }
  }

  if (identical(bin_mode, "value")) {
    labels <- if ("ssd" %in% names(tick_data)) tick_data$ssd else format_ssd_values(tick_data$x_plot)
    axis(1,
         at = tick_data$x_plot,
         labels = labels,
         cex.axis = 0.7)
    title(xlab = value_title, line = 2)
  } else if (identical(bin_mode, "global_quantile")) {
    if (!"ssd" %in% names(tick_data)) return(invisible(NULL))
    ssd_labels <- gsub("\\s+", "", tick_data$ssd)
    if (!is.null(probs) && length(probs[-1]) == nrow(tick_data)) {
      quantile_labels <- paste(format_stop_signal_quantile_intervals(probs),
                               ssd_labels,
                               sep = "\n")
    } else {
      quantile_labels <- ssd_labels
    }
    axis(1,
         at = tick_data$x_plot,
         labels = quantile_labels,
         cex.axis = 0.65)
    title(xlab = global_title, line = 3)
  } else {
    quantile_labels <- unlist(lapply(strsplit(format_stop_signal_quantile_intervals(c(0, tick_data$x_plot)), ","),
                                     function(x) paste0(x[1],",\n",x[2])))
    axis(1,
         at = tick_data$x_plot,
         labels = quantile_labels,
         cex.axis = 0.7)
    title(xlab = participant_title, line = 2)
  }

  invisible(NULL)
}

plot_stop_signal_summary <- function(input,
                                     post_predict = NULL,
                                     prior_predict = NULL,
                                     probs = seq(0, 1, length.out = 5),
                                     factors = NULL,
                                     within_plot = NULL,
                                     use_global_quantiles = FALSE,
                                     ssd_binning = c("quantile", "value"),
                                     ssd_round = NULL,
                                     on_duplicate_quantiles = c("error", "value", "reduce"),
                                     subject = NULL,
                                     quants = c(0.025, 0.975),
                                     functions = NULL,
                                     n_cores = 1,
                                     n_post = 50,
                                     layout = NA,
                                     to_plot = c("data", "posterior", "prior")[1:2],
                                     use_lim = c("data", "posterior", "prior")[1:2],
                                     legendpos = c("topleft", "bottomright"),
                                     posterior_args = list(),
                                     prior_args = list(),
                                     predictive_points = FALSE,
                                     predictive_point_args = list(),
                                     print_plot_data = FALSE,
                                     dots,
                                     value_col,
                                     y_label,
                                     summary_funs,
                                     default_ylim = NULL,
                                     initial_y_min = Inf,
                                     initial_y_max = -Inf,
                                     to_plot_missing = FALSE) {
  # Shared plotting engine for plot_ss_if() and plot_ss_srrt(). The wrappers
  # supply the metric-specific summary functions and y-axis label; this helper
  # owns validation, predictive handling, bin fallback selection, axis limits,
  # panel drawing, and the returned plot-data table.
  ssd_binning <- match.arg(ssd_binning)
  on_duplicate_quantiles <- match.arg(on_duplicate_quantiles)
  if (!is.null(ssd_round) &&
      (!is.numeric(ssd_round) || length(ssd_round) != 1 || is.na(ssd_round) || ssd_round <= 0)) {
    stop("`ssd_round` must be NULL or a single positive numeric value.")
  }
  if (ssd_binning == "value") message_ignored_value_binning_args()
  if (ssd_binning == "quantile") validate_stop_signal_probs(probs)
  validate_stop_signal_quants(quants)
  validate_stop_signal_sources(to_plot, use_lim)
  if (!is.logical(predictive_points) || length(predictive_points) != 1 ||
      is.na(predictive_points)) {
    stop("`predictive_points` must be TRUE or FALSE.")
  }
  if (!is.logical(print_plot_data) || length(print_plot_data) != 1 ||
      is.na(print_plot_data)) {
    stop("`print_plot_data` must be TRUE or FALSE.")
  }
  validate_stop_signal_predictive_request(input, post_predict, prior_predict,
                                          to_plot, to_plot_missing)
  layout <- normalize_stop_signal_layout(layout)
  legendpos <- normalize_stop_signal_legendpos(legendpos)

  check <- prep_data_plot(input, post_predict, prior_predict, to_plot, use_lim,
                          factors, within_plot, subject, n_cores, n_post,
                          remove_na = FALSE, functions)
  data_sources <- check$datasets
  sources <- check$sources

  dots <- add_defaults(dots, col = c("black", "#A9A9A9", "#666666"))
  posterior_args <- add_defaults(posterior_args, col = c("darkgreen", "#0000FF", "#008B8B"))
  prior_args <- add_defaults(prior_args, col = c("red", "#800080", "#CC00FF"))

  unique_group_keys <- levels(factor(data_sources[[1]]$group_key))
  if (is.null(within_plot)) {
    within_levels <- "All"
  } else {
    within_levels <- levels(factor(data_sources[[1]][[within_plot]]))
  }
  line_types <- seq_along(within_levels)

  summary_by_source <- list()
  draw_quantiles_by_source <- list()
  source_bin_modes <- character()

  # Choose one SSD binning rule for all plotted sources. If quantile breaks are
  # duplicated anywhere, the requested fallback is applied globally so observed
  # and predictive summaries are drawn on the same x scale.
  effective_ssd_binning <- ssd_binning
  effective_probs <- probs
  plot_global_ssd_breaks <- NULL
  if (ssd_binning == "quantile") {
    has_duplicate_breaks <- has_duplicate_stop_signal_quantiles(
      data_sources, sources, probs, use_global_quantiles, within_plot, ssd_round
    )

    if (has_duplicate_breaks) {
      if (on_duplicate_quantiles == "value") {
        effective_ssd_binning <- "value"
      } else if (on_duplicate_quantiles == "reduce") {
        effective_probs <- get_reduced_stop_signal_probs(
          data_sources, sources, probs, use_global_quantiles, within_plot, ssd_round
        )
      } else {
        stop("Duplicate quantile values detected, or plotted sources have different ",
             "SSD support for a shared global quantile grid. Please use fewer bins, ",
             "match the SSD design across sources, or set `on_duplicate_quantiles` ",
             "to \"value\" or \"reduce\".")
      }
    }
  }
  if (effective_ssd_binning == "quantile" && use_global_quantiles) {
    plot_global_ssd_breaks <- get_reference_global_ssd_breaks(data_sources, sources,
                                                              effective_probs, ssd_round)
  }

  y_min <- initial_y_min
  y_max <- initial_y_max
  x_min <- Inf
  x_max <- -Inf

  # First pass: summarize every source. Observed data produce one summary per
  # group/within level; predictive sources first summarize each postn draw and
  # then take quantiles across those draw-level summaries.
  for (k in seq_along(data_sources)) {
    df <- data_sources[[k]]
    styp <- sources[k]
    sname <- names(sources)[k]

    if (is.null(df) || !nrow(df)) next
    if (is.null(df$SSD)) stop("No SSD column in data.")

    splitted <- split(df, df$group_key)
    dots$ssd_round <- ssd_round

    if (effective_ssd_binning == "value") {
      quantile_fun <- summary_funs$value
      source_bin_modes[sname] <- "value"
    } else if (use_global_quantiles) {
      quantile_fun <- summary_funs$global
      dots$global_ssd_breaks <- plot_global_ssd_breaks
      source_bin_modes[sname] <- "global_quantile"
    } else {
      quantile_fun <- summary_funs$individual
      source_bin_modes[sname] <- "individual_quantile"
    }

    if ("postn" %in% names(df)) {
      summary_by_source[[sname]] <- lapply(splitted, function(sub_grp) {
        postn_splits <- split(sub_grp, sub_grp$postn)
        lapply(postn_splits, quantile_fun,
               group_factor = within_plot, probs = effective_probs, dots = dots)
      })

      draw_quantiles_by_source[[sname]] <- lapply(summary_by_source[[sname]], function(postn_list) {
        out <- list()
        for (lev in within_levels) {
          out[[lev]] <- get_stop_signal_postn_quantiles(
            postn_list = postn_list,
            level = lev,
            value_col = value_col,
            quants = quants
          )
        }
        out
      })

      if (styp %in% use_lim) {
        x_vals <- unlist(lapply(draw_quantiles_by_source[[sname]], function(group_val) {
          sapply(group_val, function(mat4) {
            if (is.null(mat4)) return(NULL)
            get_stop_signal_x_range_values(mat4[nrow(mat4), ])
          })
        }))
        x_min <- min(x_min, x_vals, na.rm = TRUE)
        x_max <- max(x_max, x_vals, na.rm = TRUE)

        y_vals <- unlist(lapply(draw_quantiles_by_source[[sname]], function(group_val) {
          lapply(group_val, function(mat4) {
            if (is.null(mat4)) return(NULL)
            range(c(mat4[1, ], mat4[3, ]), na.rm = TRUE)
          })
        }))
        y_min <- min(y_min, y_vals, na.rm = TRUE)
        y_max <- max(y_max, y_vals, na.rm = TRUE)
      }
    } else {
      summary_by_source[[sname]] <- lapply(splitted, quantile_fun,
                                     group_factor = within_plot,
                                     probs = effective_probs,
                                     dots = dots)

      if (styp %in% use_lim) {
        all_x_vals <- c()
        all_y_vals <- c()
        for (grp_name in names(summary_by_source[[sname]])) {
          summary_group <- summary_by_source[[sname]][[grp_name]]
          if (!is.null(summary_group)) {
            for (draw in summary_group) {
              if (!is.null(draw[["x_plot"]]) && !is.null(draw[[value_col]])) {
                all_x_vals <- c(all_x_vals, get_stop_signal_x_range_values(draw$x_plot))
                all_y_vals <- c(all_y_vals, get_stop_signal_y_range_values(draw, value_col))
              }
            }
          }
        }

        x_min <- min(x_min, all_x_vals, na.rm = TRUE)
        x_max <- max(x_max, all_x_vals, na.rm = TRUE)
        y_min <- min(y_min, all_y_vals, na.rm = TRUE)
        y_max <- max(y_max, all_y_vals, na.rm = TRUE)
      }
    }
  }

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  if (any(is.na(layout))) {
    par(mfrow = coda_setmfrow(Nchains = 1, Nparms = length(unique_group_keys), nplots = 1))
  } else {
    par(mfrow = layout)
  }

  if (!is.finite(y_min) || !is.finite(y_max) || y_min >= y_max) {
    ylim <- default_ylim
  } else {
    y_buffer <- 0.05 * (y_max - y_min)
    ylim <- c(y_min - y_buffer, y_max + y_buffer)
  }

  if (!is.finite(x_min) || !is.finite(x_max) || x_min >= x_max) {
    xlim <- NULL
  } else {
    x_buffer <- 0.05 * (x_max - x_min)
    xlim <- c(x_min - x_buffer, x_max + x_buffer)
  }

  # Second pass: draw one panel per factor combination, reusing the summaries
  # from the first pass. The same loop handles IF and SRRT via value_col.
  for (group_key in unique_group_keys) {
    tmp_dots <- dots
    tmp_posterior_args <- posterior_args
    tmp_prior_args <- prior_args

    plot_args <- add_defaults(dots, ylim = ylim, xlim = xlim,
                              main = group_key, xlab = "", ylab = "")
    plot_args <- fix_dots_plot(plot_args)
    plot_args$axes <- FALSE
    do.call(plot, c(list(NA), plot_args))
    axis(2)
    mtext(y_label, side = 2, line = 3)

    legend_map <- character(0)
    lwd_map <- numeric()
    for (k in seq_along(data_sources)) {
      styp <- sources[k]
      sname <- names(sources)[k]
      if (styp == "data") {
        src_args <- tmp_dots
        tmp_dots$col <- tmp_dots$col[-1]
      } else if (styp == "posterior") {
        src_args <- tmp_posterior_args
        tmp_posterior_args$col <- tmp_posterior_args$col[-1]
      } else if (styp == "prior") {
        src_args <- tmp_prior_args
        tmp_prior_args$col <- tmp_prior_args$col[-1]
      }
      legend_map[sname] <- src_args$col[1]
      lwd_map[sname] <- ifelse(is.null(src_args$lwd), 1, src_args$lwd)

      if (!is.null(summary_by_source[[sname]])) {
        if (!("postn" %in% names(data_sources[[k]]))) {
          summary_df <- summary_by_source[[sname]][[group_key]]
          if (!is.null(summary_df)) {
            ilev <- 1
            for (df in summary_df) {
              lines_args <- add_defaults(src_args, lty = line_types[ilev])
              lines_args <- fix_dots_plot(lines_args)
              do.call(lines, c(list(x = df$x_plot, y = df[[value_col]]), lines_args))
              do.call(points, c(list(x = df$x_plot, y = df[[value_col]]), lines_args))
              if (!any(sources == "posterior")) {
                draw_stop_signal_se(df$x_plot, df[[value_col]], df$se)
              }
              ilev <- ilev + 1
            }
          }
        } else {
          if (!is.null(draw_quantiles_by_source[[sname]])) {
            draw_quantiles_for_group <- draw_quantiles_by_source[[sname]][[group_key]]
            if (!is.null(draw_quantiles_for_group)) {
              ilev <- 1
              for (lev in within_levels) {
                mat4 <- draw_quantiles_for_group[[lev]]
                if (!is.null(mat4)) {
                  y_lower <- mat4[1, ]
                  y_med <- mat4[2, ]
                  y_upper <- mat4[3, ]
                  x_plot <- mat4[4, ]

                  lines_args <- add_defaults(src_args, lty = line_types[ilev])
                  lines_args <- fix_dots_plot(lines_args)
                  adj_color <- do.call(adjustcolor,
                                       fix_dots(add_defaults(src_args, alpha.f = 0.2),
                                                adjustcolor))
                  poly_args <- src_args
                  poly_args$col <- adj_color
                  poly_args <- fix_dots_plot(poly_args)
                  do.call(polygon, c(list(
                    y = c(y_lower, rev(y_upper)),
                    x = c(x_plot, rev(x_plot)),
                    border = NA
                  ), poly_args))
                  if (predictive_points) {
                    draw_stop_signal_predictive_points(
                      summary_by_source[[sname]][[group_key]], lev, value_col,
                      src_args, predictive_point_args
                    )
                  }
                  do.call(lines, c(list(x = x_plot, y = y_med), lines_args))
                }
                ilev <- ilev + 1
              }
            }
          }
        }
      }
    }

    for (k in seq_along(data_sources)) {
      sname <- names(sources)[k]
      if (!("postn" %in% names(data_sources[[k]]))) {
        tick_group <- summary_by_source[[sname]][[group_key]]
        if (!is.null(tick_group)) {
          tick_data <- tick_group[[names(tick_group)[1]]]
          draw_stop_signal_x_axis(tick_data, source_bin_modes[sname],
                                  "Global SSD Quantile Bin (%, sec.)",
                                  "Participant SSD Quantile Bin (%)",
                                  probs = effective_probs)
          break
        }
      } else {
        tick_group <- draw_quantiles_by_source[[sname]][[group_key]]
        if (!is.null(tick_group)) {
          tick_mat <- tick_group[[names(tick_group)[1]]]
          if (!is.null(tick_mat)) {
            tick_data <- data.frame(x_plot = unname(tick_mat[nrow(tick_mat), ]))
            tick_labels <- names(tick_mat[nrow(tick_mat), ])
            if (source_bin_modes[sname] %in% c("value", "global_quantile") &&
                !is.null(tick_labels)) {
              tick_data$ssd <- tick_labels
            }
            draw_stop_signal_x_axis(tick_data, source_bin_modes[sname],
                                    "Global SSD Quantile Bin (%, sec.)",
                                    "Participant SSD Quantile Bin (%)",
                                    probs = effective_probs)
            break
          }
        }
      }
    }

    if (!is.na(legendpos[1]) && !(length(within_levels) == 1 && within_levels == "All")) {
      legend(legendpos[1], legend = within_levels, lty = line_types, col = "black",
             title = within_plot, bty = "n")
    }

    if (length(data_sources) > 1) {
      if (!is.na(legendpos[2])) {
        legend(legendpos[2], legend = names(legend_map), lty = 1, col = legend_map,
               title = "Source", bty = "n", lwd = lwd_map)
      }
    }
  }

  plot_data <- get_stop_signal_plot_data(summary_by_source, draw_quantiles_by_source, sources,
                                         source_bin_modes, value_col)
  if (print_plot_data) print_stop_signal_plot_data(plot_data)

  invisible(plot_data)
}

#' Plot Inhibition Functions
#'

#' Plots panels of the inhibition functions (probability of responding; Pr(R)) for each
#' level of specified factors of stop-signal data as a function of user-defined
#' SSD bins/categories. Optionally, posterior and/or prior predictive
#' inhibition functions can be overlaid.
#'
#' Per default, the SSD-categories are defined in terms of the quantiles of the
#' SSD distribution for each participant, and then averaged over participants (see `use_global_quantiles`).
#'
#' Posterior/prior predictive intervals are quantile intervals of the plotted
#' summary across predictive draws. Observed data are plotted with error bars
#' (plus/minus the standard error per SSD bin/category).
#'
#' @param input Either an emc object or a stop-signal data frame, or a *list* of such objects. SSD column in data required.
#' @param post_predict Optional posterior predictive data (matching columns) or *list* thereof.
#' @param prior_predict Optional prior predictive data (matching columns) or list thereof.
#' @param probs Strictly increasing numeric vector from 0 to 1 that defines SSD quantile bins/categories.
#' @param factors Character vector of factor names to aggregate over; defaults to plotting full data set ungrouped by factors if NULL.
#' @param within_plot Character indicating factor for which inhibition functions are plotted in the same panel
#' @param use_global_quantiles If set to `TRUE`, SSDs are pooled over participants before calculating quantiles, so
#' the same absolute SSD range is used to get Pr(R) for each participant,
#' and then these probabilities are averaged over participants.
#' @param ssd_binning Character. `"quantile"` uses SSD quantile bins; `"value"` groups by SSD values.
#' @param ssd_round Optional numeric bin width for rounding SSD values before stop-signal binning.
#' @param on_duplicate_quantiles Character. What to do when quantile breaks are duplicated:
#'   `"error"` stops with a message; `"value"` falls back to value-based SSD bins
#'   for all plotted sources; `"reduce"` uses the largest evenly spaced quantile
#'   grid without duplicated breaks.
#' @param subject Subset the data to a single subject (by index or name).
#' @param quants Numeric vector of exactly two predictive interval quantile bounds
#'   that bracket 0.5 (e.g. c(0.025, 0.975)).
#' @param functions A function (or list of functions) that create new columns in the datasets or predictives
#' @param n_cores Number of CPU cores to use if generating predictives from an emc object.
#' @param n_post Number of posterior draws to simulate if needed for predictives.
#' @param layout Numeric vector used in par(mfrow=...); use NA for auto-layout.
#' @param to_plot Character vector: any of "data", "posterior", "prior".
#' @param use_lim Character vector controlling which source(s) define axis limits
#' @param legendpos Character vector controlling the positions of the legends
#' @param posterior_args Optional list of graphical parameters for posterior lines/ribbons.
#' @param prior_args Optional list of graphical parameters for prior lines/ribbons.
#' @param predictive_points Logical. If `TRUE`, plot the draw-level posterior/prior
#'   predictive summaries as points.
#' @param predictive_point_args Optional list of graphical parameters for
#'   predictive points. Supports `jitter` for deterministic horizontal jitter.
#' @param print_plot_data Logical. If `TRUE`, print the plotted summary data and
#'   count columns used to form each plotted point/bin.
#' @param ... Other graphical parameters for the real data lines.
#'
#' @return Invisibly returns the plotted summary data.
#' @export
plot_ss_if <- function(input,
                    post_predict = NULL,
                    prior_predict = NULL,
                    probs = seq(0,1, length.out=5),
                    factors = NULL,
                    within_plot = NULL,
                    use_global_quantiles = FALSE,
                    ssd_binning = c("quantile", "value"),
                    ssd_round = NULL,
                    on_duplicate_quantiles = c("error", "value", "reduce"),
                    subject = NULL,
                    quants = c(0.025, 0.975), functions = NULL,
                    n_cores = 1,
                    n_post = 50,
                    layout = NA,
                    to_plot = c('data','posterior','prior')[1:2],
                    use_lim = c('data','posterior','prior')[1:2],
                    legendpos = c('topleft', 'bottomright'),
                    posterior_args = list(),
                    prior_args = list(),
                    predictive_points = FALSE,
                    predictive_point_args = list(),
                    print_plot_data = FALSE,
                    ...) {

  plot_stop_signal_summary(
    input = input,
    post_predict = post_predict,
    prior_predict = prior_predict,
    probs = probs,
    factors = factors,
    within_plot = within_plot,
    use_global_quantiles = use_global_quantiles,
    ssd_binning = ssd_binning,
    ssd_round = ssd_round,
    on_duplicate_quantiles = on_duplicate_quantiles,
    subject = subject,
    quants = quants,
    functions = functions,
    n_cores = n_cores,
    n_post = n_post,
    layout = layout,
    to_plot = to_plot,
    use_lim = use_lim,
    legendpos = legendpos,
    posterior_args = posterior_args,
    prior_args = prior_args,
    predictive_points = predictive_points,
    predictive_point_args = predictive_point_args,
    print_plot_data = print_plot_data,
    dots = list(...),
    value_col = "p_response",
    y_label = "Pr(R)",
    summary_funs = list(
      individual = get_response_probability_by_individual_ssd_quantile,
      global = get_response_probability_by_global_ssd_quantile,
      value = get_response_probability_by_ssd_value
    ),
    default_ylim = c(0, 1),
    initial_y_min = 0,
    initial_y_max = 0,
    to_plot_missing = missing(to_plot)
  )
}

get_response_probability_by_individual_ssd_quantile <- function(x, group_factor, probs, dots) {
  # Filter: only rows with finite SSD and subject info
  x <- x[is.finite(x[["SSD"]]) & !is.na(x[["subjects"]]), ]
  x$SSD_bin_value <- get_stop_signal_bin_ssd(x, dots$ssd_round)
  subjects <- unique(x$subjects)

  # Get grouping levels if group_factor is specified
  if (!is.null(group_factor)) {
    group_vals <- unique(x[[group_factor]])
    group_names <- as.character(group_vals)
  } else {
    group_names <- "All"
  }

  # Function to compute stats for a subset of data
  compute_group_stats <- function(data_subset) {
    subj_stats <- lapply(unique(data_subset$subjects), function(s) {
      df <- data_subset[data_subset$subjects == s, ]

      ssd_quants <- quantile(df$SSD_bin_value, probs = probs, na.rm = TRUE)
      if (anyDuplicated(ssd_quants)) {
        stop("Duplicate quantile values detected. Please use fewer bins.")
      }
      ssd_quants <- unique(ssd_quants)

      df$ssd_bin <- cut(df$SSD_bin_value, breaks = ssd_quants, include.lowest = TRUE, right = TRUE)

      bin_stats <- aggregate(I(!is.na(R)) ~ ssd_bin, data = df, FUN = function(v) {
        n <- sum(!is.na(v))
        mean_p <- mean(v, na.rm = TRUE)
        se_p <- sqrt(mean_p * (1 - mean_p) / n)
        c(mean = mean_p, se = se_p, n = n)
      }, na.action = na.pass)

      bin_stats <- do.call(data.frame, bin_stats)

      names(bin_stats) <- c("ssd_bin", "p_response", "se", "n")


      bin_stats$ssd <- as.numeric(gsub(".*,", "", gsub("\\[|\\]|\\(|\\)", "", bin_stats$ssd_bin)))
      bin_stats$subject <- s
      bin_stats$x_plot <- probs[-1][match(as.character(bin_stats$ssd_bin),
                                          levels(df$ssd_bin))]

      return(bin_stats)
    })

    if(length(subj_stats) > 1){
      all_stats <- do.call(rbind, subj_stats)
      summary_stats <- aggregate(p_response ~ x_plot, data = all_stats,
                                 FUN = function(v) {
                                   n <- length(v)
                                   mean_p <- mean(v, na.rm = TRUE)
                                   se_p <- sd(v, na.rm = TRUE)/sqrt(n)
                                   c(mean = mean_p, se = se_p, n = n)
                                 },
                                 na.action = na.pass)
      summary_stats <- do.call(data.frame, summary_stats)

      summary_stats <- summary_stats[order(summary_stats$x_plot), ]
      names(summary_stats) <- c("x_plot", "p_response", "se", "n")
      n_obs_lookup <- aggregate(n ~ x_plot, data = all_stats, sum, na.rm = TRUE)
      names(n_obs_lookup) <- c("x_plot", "n_obs")
      summary_stats <- merge(summary_stats, n_obs_lookup, by = "x_plot", all.x = TRUE)
      summary_stats <- summary_stats[order(summary_stats$x_plot), ]
    } else {
      summary_stats <- subj_stats[[1]]
      summary_stats$n_obs <- summary_stats$n
    }

    return(summary_stats)

  }

  # Loop over group levels (or just "All")
  out <- lapply(group_names, function(g) {
    if (g == "All") {
      data_subset <- x
    } else {
      data_subset <- x[x[[group_factor]] == g, ]
    }
    compute_group_stats(data_subset)
  })

  names(out) <- group_names
  return(out)
}

get_response_probability_by_global_ssd_quantile <- function(x, group_factor, probs, dots) {
  # Check global SSD breaks
  if (!"global_ssd_breaks" %in% names(dots)) {
    stop("Missing 'global_ssd_breaks' in dots")
  }
  global_ssd_breaks <- dots$global_ssd_breaks

  # Filter: only rows with finite SSD and subject info
  x <- x[is.finite(x[["SSD"]]) & !is.na(x[["subjects"]]), ]
  x$SSD_bin_value <- get_stop_signal_bin_ssd(x, dots$ssd_round)
  subjects <- unique(x$subjects)

  # Get grouping levels if group_factor is specified
  if (!is.null(group_factor)) {
    group_vals <- unique(x[[group_factor]])
    group_names <- as.character(group_vals)
  } else {
    group_names <- "All"
  }

  # Create factor labels for breaks
  bin_labels <- cut(global_ssd_breaks[-1], breaks = global_ssd_breaks, include.lowest = TRUE, right = TRUE)
  bin_labels <- levels(cut(x$SSD_bin_value, breaks = global_ssd_breaks, include.lowest = TRUE, right = TRUE))

  # Function to compute stats for a subset of data
  compute_group_stats <- function(data_subset) {
    subj_stats <- lapply(unique(data_subset$subjects), function(s) {
      df <- data_subset[data_subset$subjects == s, ]

      df$ssd_bin <- cut(df$SSD_bin_value, breaks = global_ssd_breaks, include.lowest = TRUE, right = TRUE)

      bin_stats <- aggregate(I(!is.na(R)) ~ ssd_bin, data = df, FUN = function(v) {
        n <- sum(!is.na(v))
        mean_p <- mean(v, na.rm = TRUE)
        se_p <- sqrt(mean_p * (1 - mean_p) / n)
        c(mean = mean_p, se = se_p, n = n)
      }, na.action = na.pass)

      bin_stats <- do.call(data.frame, bin_stats)

      names(bin_stats) <- c("ssd_bin", "p_response", "se", "n")

      bin_stats$ssd <- as.character(bin_stats$ssd_bin)
      bin_stats$x_plot <- as.numeric(gsub(".*,", "", gsub("\\[|\\]|\\(|\\)", "", bin_stats$ssd_bin)))
      bin_stats$subject <- s

      return(bin_stats)
    })

    if(length(subj_stats) > 1){
      all_stats <- do.call(rbind, subj_stats)
      summary_stats <- aggregate(p_response ~ ssd, data = all_stats,
                                 FUN = function(v) {
                                   n <- sum(!is.na(v))
                                   mean_p <- mean(v, na.rm = TRUE)
                                   se_p <- sd(v, na.rm = TRUE)/sqrt(n)
                                   c(mean = mean_p, se = se_p, n = n)
                                 },
                                 na.action = na.pass)
      summary_stats <- do.call(data.frame, summary_stats)

      # Add x_plot by matching back from all_stats
      x_plot_lookup <- aggregate(x_plot ~ ssd, data = all_stats, FUN = function(x) unique(x)[1])
      summary_stats <- merge(summary_stats, x_plot_lookup, by = "ssd")

      summary_stats <- summary_stats[order(summary_stats$x_plot), ]
      names(summary_stats) <- c("ssd", "p_response", "se", "n", "x_plot")
      n_obs_lookup <- aggregate(n ~ ssd, data = all_stats, sum, na.rm = TRUE)
      names(n_obs_lookup) <- c("ssd", "n_obs")
      summary_stats <- merge(summary_stats, n_obs_lookup, by = "ssd", all.x = TRUE)
      summary_stats <- summary_stats[order(summary_stats$x_plot), ]
    } else {
      summary_stats <- subj_stats[[1]]
      summary_stats$n_obs <- summary_stats$n
    }
    return(summary_stats)
  }

  # Loop over group levels (or just "All")
  out <- lapply(group_names, function(g) {
    if (g == "All") {
      data_subset <- x
    } else {
      data_subset <- x[x[[group_factor]] == g, ]
    }
    compute_group_stats(data_subset)
  })

  names(out) <- group_names
  return(out)
}

get_response_probability_by_ssd_value <- function(x, group_factor, probs, dots) {
  # Value binning is for discrete/staircase SSD designs. SSDs have already been
  # routed through the same rounding helper used by the quantile paths.
  x <- x[is.finite(x[["SSD"]]) & !is.na(x[["subjects"]]), ]
  x$SSD_value <- get_stop_signal_bin_ssd(x, dots$ssd_round)

  if (!is.null(group_factor)) {
    group_vals <- unique(x[[group_factor]])
    group_names <- as.character(group_vals)
  } else {
    group_names <- "All"
  }

  compute_group_stats <- function(data_subset) {
    subj_stats <- lapply(unique(data_subset$subjects), function(s) {
      df <- data_subset[data_subset$subjects == s, ]

      bin_stats <- aggregate(I(!is.na(R)) ~ SSD_value, data = df, FUN = function(v) {
        n <- sum(!is.na(v))
        mean_p <- mean(v, na.rm = TRUE)
        se_p <- sqrt(mean_p * (1 - mean_p) / n)
        c(mean = mean_p, se = se_p, n = n)
      }, na.action = na.pass)

      bin_stats <- do.call(data.frame, bin_stats)
      names(bin_stats) <- c("x_plot", "p_response", "se", "n")
      bin_stats$ssd <- format_ssd_values(bin_stats$x_plot)
      bin_stats$subject <- s
      bin_stats <- bin_stats[order(bin_stats$x_plot), ]

      bin_stats
    })

    if(length(subj_stats) > 1){
      all_stats <- do.call(rbind, subj_stats)
      summary_stats <- aggregate(p_response ~ x_plot, data = all_stats,
                                 FUN = function(v) {
                                   n <- sum(!is.na(v))
                                   mean_p <- mean(v, na.rm = TRUE)
                                   se_p <- sd(v, na.rm = TRUE)/sqrt(n)
                                   c(mean = mean_p, se = se_p, n = n)
                                 },
                                 na.action = na.pass)
      summary_stats <- do.call(data.frame, summary_stats)
      names(summary_stats) <- c("x_plot", "p_response", "se", "n")
      n_obs_lookup <- aggregate(n ~ x_plot, data = all_stats, sum, na.rm = TRUE)
      names(n_obs_lookup) <- c("x_plot", "n_obs")
      summary_stats <- merge(summary_stats, n_obs_lookup, by = "x_plot", all.x = TRUE)
      summary_stats$ssd <- format_ssd_values(summary_stats$x_plot)
      summary_stats <- summary_stats[order(summary_stats$x_plot), ]
    } else {
      summary_stats <- subj_stats[[1]]
      summary_stats$n_obs <- summary_stats$n
    }

    summary_stats
  }

  out <- lapply(group_names, function(g) {
    if (g == "All") {
      data_subset <- x
    } else {
      data_subset <- x[x[[group_factor]] == g, ]
    }
    compute_group_stats(data_subset)
  })

  names(out) <- group_names
  return(out)
}

#' Plot Mean SRRT
#'
#' Plots panels of the mean signal-respond response time (SRRT) as a function of user-defined
#' SSD bins/categories for each level of specified factors of stop-signal data.
#' Optionally, posterior and/or prior predictive inhibition functions can be overlaid.
#'
#' Per default, SSDs are pooled over participants before calculating quantiles, so
#' the same absolute SSD range is used to get mean SRRT for each participant,
#' and then these probabilities are averaged over participants (see `use_global_quantiles`).
#'
#' Posterior/prior predictive intervals are quantile intervals of the plotted
#' summary across predictive draws. Observed data are plotted with error bars
#' (plus/minus the standard error per SSD bin/category).
#'
#' For SRRT summaries, the number of signal-respond RTs contributing to a
#' predictive draw can vary across SSD bins because response occurrence is itself
#' simulated. Sparse high-SSD bins can therefore produce predictive intervals based
#' on fewer contributing draw-level summaries.
#'
#' @param input Either an emc object or a stop-signal data frame, or a list of such objects. SSD column in data required.
#' @param post_predict Optional posterior predictive data (matching columns) or list thereof.
#' @param prior_predict Optional prior predictive data (matching columns) or list thereof.
#' @param probs Strictly increasing numeric vector from 0 to 1 that defines SSD quantile bins/categories.
#' @param factors Character vector of factor names to aggregate over; defaults to plotting full data set ungrouped by factors if NULL.
#' @param within_plot Character indicating factor for which inhibition functions are plotted in the same panel
#' @param use_global_quantiles If set to FALSE, the SSD-categories are defined in terms of the quantiles of the
#' SSD distribution for each participant, and then averaged over participants.
#' @param ssd_binning Character. `"quantile"` uses SSD quantile bins; `"value"` groups by SSD values.
#' @param ssd_round Optional numeric bin width for rounding SSD values before stop-signal binning.
#' @param on_duplicate_quantiles Character. What to do when quantile breaks are duplicated:
#'   `"error"` stops with a message; `"value"` falls back to value-based SSD bins
#'   for all plotted sources; `"reduce"` uses the largest evenly spaced quantile
#'   grid without duplicated breaks.
#' @param subject Subset the data to a single subject (by index or name).
#' @param quants Numeric vector of exactly two predictive interval quantile bounds
#'   that bracket 0.5 (e.g. c(0.025, 0.975)).
#' @param functions A function (or list of functions) that create new columns in the datasets or predictives
#' @param n_cores Number of CPU cores to use if generating predictives from an emc object.
#' @param n_post Number of posterior draws to simulate if needed for predictives.
#' @param layout Numeric vector used in par(mfrow=...); use NA for auto-layout.
#' @param to_plot Character vector: any of "data", "posterior", "prior".
#' @param use_lim Character vector controlling which source(s) define axis limits
#' @param legendpos Character vector controlling the positions of the legends
#' @param posterior_args Optional list of graphical parameters for posterior lines/ribbons.
#' @param prior_args Optional list of graphical parameters for prior lines/ribbons.
#' @param predictive_points Logical. If `TRUE`, plot the draw-level posterior/prior
#'   predictive summaries as points.
#' @param predictive_point_args Optional list of graphical parameters for
#'   predictive points. Supports `jitter` for deterministic horizontal jitter.
#' @param print_plot_data Logical. If `TRUE`, print the plotted summary data and
#'   count columns used to form each plotted point/bin.
#' @param ... Other graphical parameters for the real data lines.
#'
#' @return Invisibly returns the plotted summary data.
#' @export
plot_ss_srrt <- function(input,
                      post_predict = NULL,
                      prior_predict = NULL,
                      probs = seq(0,1, length.out=5),
                      factors = NULL,
                      within_plot = NULL,
                      use_global_quantiles = TRUE,
                      ssd_binning = c("quantile", "value"),
                      ssd_round = NULL,
                      on_duplicate_quantiles = c("error", "value", "reduce"),
                      subject = NULL,
                      quants = c(0.025, 0.975), functions = NULL,
                      n_cores = 1,
                      n_post = 50,
                      layout = NA,
                      to_plot = c('data','posterior','prior')[1:2],
                      use_lim = c('data','posterior','prior')[1:2],
                      legendpos = c('topleft', 'bottomright'),
                      posterior_args = list(),
                      prior_args = list(),
                      predictive_points = FALSE,
                      predictive_point_args = list(),
                      print_plot_data = FALSE,
                      ...) {

  plot_stop_signal_summary(
    input = input,
    post_predict = post_predict,
    prior_predict = prior_predict,
    probs = probs,
    factors = factors,
    within_plot = within_plot,
    use_global_quantiles = use_global_quantiles,
    ssd_binning = ssd_binning,
    ssd_round = ssd_round,
    on_duplicate_quantiles = on_duplicate_quantiles,
    subject = subject,
    quants = quants,
    functions = functions,
    n_cores = n_cores,
    n_post = n_post,
    layout = layout,
    to_plot = to_plot,
    use_lim = use_lim,
    legendpos = legendpos,
    posterior_args = posterior_args,
    prior_args = prior_args,
    predictive_points = predictive_points,
    predictive_point_args = predictive_point_args,
    print_plot_data = print_plot_data,
    dots = list(...),
    value_col = "srrt",
    y_label = "Mean SRRT (sec.)",
    summary_funs = list(
      individual = get_srrt_by_individual_ssd_quantile,
      global = get_srrt_by_global_ssd_quantile,
      value = get_srrt_by_ssd_value
    ),
    default_ylim = NULL,
    initial_y_min = Inf,
    initial_y_max = -Inf,
    to_plot_missing = missing(to_plot)
  )
}

get_srrt_by_individual_ssd_quantile <- function(x, group_factor, probs, dots) {
  # Filter: only rows with finite SSD and subject info
  x <- x[is.finite(x[["SSD"]]) & !is.na(x[["subjects"]]), ]
  x$SSD_bin_value <- get_stop_signal_bin_ssd(x, dots$ssd_round)
  subjects <- unique(x$subjects)

  # Get grouping levels if group_factor is specified
  if (!is.null(group_factor)) {
    group_vals <- unique(x[[group_factor]])
    group_names <- as.character(group_vals)
  } else {
    group_names <- "All"
  }

  # Function to compute stats for a subset of data
  compute_group_stats <- function(data_subset) {
    subj_stats <- lapply(unique(data_subset$subjects), function(s) {
      df <- data_subset[data_subset$subjects == s, ]
      ssd_quants <- quantile(df$SSD_bin_value, probs = probs, na.rm = TRUE)
      if (anyDuplicated(ssd_quants)) {
        stop("Duplicate quantile values detected. Please use fewer bins.")
      }
      ssd_quants <- unique(ssd_quants)

      df$ssd_bin <- cut(df$SSD_bin_value, breaks = ssd_quants, include.lowest = TRUE, right = TRUE)

      # Compute mean RT (srrt), count (n), and SD (rt_sd)
      bin_stats <- aggregate(rt ~ ssd_bin, data = df,
                             FUN = function(v) {
                               n <- sum(!is.na(v))
                               mean_rt <- mean(v, na.rm = TRUE)
                               se_rt <- sd(v, na.rm = TRUE) / sqrt(n)
                               c(mean = mean_rt, se = se_rt, n = n)
                             },
                             na.action = na.pass)
      # Convert matrix columns to individual columns
      bin_stats <- do.call(data.frame, bin_stats)
      names(bin_stats)[2:4] <- c("srrt", "se", "n")

      bin_stats$ssd <- as.numeric(gsub(".*,", "", gsub("\\[|\\]|\\(|\\)", "", bin_stats$ssd_bin)))
      bin_stats$subject <- s
      bin_stats$x_plot <- probs[-1][match(as.character(bin_stats$ssd_bin),
                                          levels(df$ssd_bin))]

      return(bin_stats)
    })

    if(length(subj_stats) > 1){
      all_stats <- do.call(rbind, subj_stats)
      summary_stats <- aggregate(srrt ~ x_plot, data = all_stats, FUN = function(v) {
        n <- sum(!is.na(v))
        mean_rt <- mean(v, na.rm = TRUE)
        se_rt <- sd(v, na.rm = TRUE) / sqrt(n)
        c(mean = mean_rt, se = se_rt, n = n)
      },
      na.action = na.pass)
      summary_stats <- do.call(data.frame, summary_stats)
      names(summary_stats) <- c("x_plot", "srrt", "se", "n")
      n_obs_lookup <- aggregate(n ~ x_plot, data = all_stats, sum, na.rm = TRUE)
      names(n_obs_lookup) <- c("x_plot", "n_obs")
      summary_stats <- merge(summary_stats, n_obs_lookup, by = "x_plot", all.x = TRUE)
      summary_stats <- summary_stats[order(summary_stats$x_plot), ]
    } else {
      summary_stats <- subj_stats[[1]][, c("x_plot", "srrt", "ssd", "se", "n")]
      summary_stats$n_obs <- summary_stats$n
    }

    return(summary_stats)
  }

  # Loop over group levels (or just "All")
  out <- lapply(group_names, function(g) {
    if (g == "All") {
      data_subset <- x
    } else {
      data_subset <- x[x[[group_factor]] == g, ]
    }
    compute_group_stats(data_subset)
  })

  names(out) <- group_names
  return(out)
}

get_srrt_by_global_ssd_quantile <- function(x, group_factor, probs, dots) {
  # Check global SSD breaks
  if (!"global_ssd_breaks" %in% names(dots)) {
    stop("Missing 'global_ssd_breaks' in dots")
  }
  global_ssd_breaks <- dots$global_ssd_breaks

  # Filter: only rows with finite SSD and subject info
  x <- x[is.finite(x[["SSD"]]) & !is.na(x[["subjects"]]), ]
  x$SSD_bin_value <- get_stop_signal_bin_ssd(x, dots$ssd_round)
  subjects <- unique(x$subjects)

  # Get grouping levels if group_factor is specified
  if (!is.null(group_factor)) {
    group_vals <- unique(x[[group_factor]])
    group_names <- as.character(group_vals)
  } else {
    group_names <- "All"
  }

  # Create factor labels for breaks (used for consistent bin labeling)
  bin_labels <- cut(global_ssd_breaks[-1], breaks = global_ssd_breaks, include.lowest = TRUE, right = TRUE)
  bin_labels <- levels(cut(x$SSD_bin_value, breaks = global_ssd_breaks, include.lowest = TRUE, right = TRUE))

  # Function to compute stats for a subset of data
  compute_group_stats <- function(data_subset) {
    subj_stats <- lapply(unique(data_subset$subjects), function(s) {
      df <- data_subset[data_subset$subjects == s, ]

      df$ssd_bin <- cut(df$SSD_bin_value, breaks = global_ssd_breaks, include.lowest = TRUE, right = TRUE)


      # Compute mean RT (srrt), count (n), and SD (rt_sd)
      bin_stats <- aggregate(rt ~ ssd_bin, data = df,
                             FUN = function(v) {
                               n <- sum(!is.na(v))
                               mean_rt <- mean(v, na.rm = TRUE)
                               se_rt <- sd(v, na.rm = TRUE) / sqrt(n)
                               c(mean = mean_rt, se = se_rt, n = n)
                             },
                             na.action = na.pass)
      # Convert matrix columns to individual columns
      bin_stats <- do.call(data.frame, bin_stats)
      names(bin_stats)[2:4] <- c("srrt", "se", "n")

      bin_stats$ssd <- as.character(bin_stats$ssd_bin)
      bin_stats$x_plot <- as.numeric(gsub(".*,", "", gsub("\\[|\\]|\\(|\\)", "", bin_stats$ssd_bin)))

      bin_stats$subject <- s

      return(bin_stats)
    })

    if(length(subj_stats) > 1){
      all_stats <- do.call(rbind, subj_stats)
      summary_stats <- aggregate(srrt ~ ssd, data = all_stats, FUN = function(v) {
        n <- sum(!is.na(v))
        mean_rt <- mean(v, na.rm = TRUE)
        se_rt <- sd(v, na.rm = TRUE) / sqrt(n)
        c(mean = mean_rt, se = se_rt, n = n)
      },
      na.action = na.pass)
      summary_stats <- do.call(data.frame, summary_stats)
      names(summary_stats) <- c("ssd", "srrt", "se", "n")

      x_plot_lookup <- aggregate(x_plot ~ ssd, data = all_stats, FUN = function(x) unique(x)[1])

      summary_stats <- merge(summary_stats, x_plot_lookup, by = "ssd")

      summary_stats <- summary_stats[order(summary_stats$x_plot), ]
      n_obs_lookup <- aggregate(n ~ ssd, data = all_stats, sum, na.rm = TRUE)
      names(n_obs_lookup) <- c("ssd", "n_obs")
      summary_stats <- merge(summary_stats, n_obs_lookup, by = "ssd", all.x = TRUE)
      summary_stats <- summary_stats[order(summary_stats$x_plot), ]
    } else {
      summary_stats <- subj_stats[[1]][, c("x_plot", "srrt", "ssd", "se", "n")]
      summary_stats$n_obs <- summary_stats$n
      summary_stats <- summary_stats[order(summary_stats$x_plot), ]
    }

    return(summary_stats)
  }

  # Loop over group levels (or just "All")
  out <- lapply(group_names, function(g) {
    if (g == "All") {
      data_subset <- x
    } else {
      data_subset <- x[x[[group_factor]] == g, ]
    }
    compute_group_stats(data_subset)
  })

  names(out) <- group_names
  return(out)
}

get_srrt_by_ssd_value <- function(x, group_factor, probs, dots) {
  # Value binning is for discrete/staircase SSD designs. SSDs have already been
  # routed through the same rounding helper used by the quantile paths.
  x <- x[is.finite(x[["SSD"]]) & !is.na(x[["subjects"]]), ]
  x$SSD_value <- get_stop_signal_bin_ssd(x, dots$ssd_round)

  if (!is.null(group_factor)) {
    group_vals <- unique(x[[group_factor]])
    group_names <- as.character(group_vals)
  } else {
    group_names <- "All"
  }

  compute_group_stats <- function(data_subset) {
    subj_stats <- lapply(unique(data_subset$subjects), function(s) {
      df <- data_subset[data_subset$subjects == s, ]

      bin_stats <- aggregate(rt ~ SSD_value, data = df,
                             FUN = function(v) {
                               n <- sum(!is.na(v))
                               mean_rt <- mean(v, na.rm = TRUE)
                               se_rt <- sd(v, na.rm = TRUE) / sqrt(n)
                               c(mean = mean_rt, se = se_rt, n = n)
                             },
                             na.action = na.pass)

      bin_stats <- do.call(data.frame, bin_stats)
      names(bin_stats) <- c("x_plot", "srrt", "se", "n")
      bin_stats$ssd <- format_ssd_values(bin_stats$x_plot)
      bin_stats$subject <- s
      bin_stats <- bin_stats[order(bin_stats$x_plot), ]

      bin_stats
    })

    if(length(subj_stats) > 1){
      all_stats <- do.call(rbind, subj_stats)
      summary_stats <- aggregate(srrt ~ x_plot, data = all_stats, FUN = function(v) {
        n <- sum(!is.na(v))
        mean_rt <- mean(v, na.rm = TRUE)
        se_rt <- sd(v, na.rm = TRUE) / sqrt(n)
        c(mean = mean_rt, se = se_rt, n = n)
      },
      na.action = na.pass)
      summary_stats <- do.call(data.frame, summary_stats)
      names(summary_stats) <- c("x_plot", "srrt", "se", "n")
      n_obs_lookup <- aggregate(n ~ x_plot, data = all_stats, sum, na.rm = TRUE)
      names(n_obs_lookup) <- c("x_plot", "n_obs")
      summary_stats <- merge(summary_stats, n_obs_lookup, by = "x_plot", all.x = TRUE)
      summary_stats$ssd <- format_ssd_values(summary_stats$x_plot)
      summary_stats <- summary_stats[order(summary_stats$x_plot), ]
    } else {
      summary_stats <- subj_stats[[1]]
      summary_stats$n_obs <- summary_stats$n
    }

    summary_stats
  }

  out <- lapply(group_names, function(g) {
    if (g == "All") {
      data_subset <- x
    } else {
      data_subset <- x[x[[group_factor]] == g, ]
    }
    compute_group_stats(data_subset)
  })

  names(out) <- group_names
  return(out)
}

compare_obs_vs_postpred <- function(obs_df, pred_df) {

  # Only stop trials
  obs_df <- obs_df[is.finite(obs_df$SSD), ]
  pred_df <- pred_df[is.finite(pred_df$SSD), ]

  # Compute observed response rates (using proportion of missing R as inhibition)
  obs_rates <- aggregate(is.na(R) ~ subjects + SSD, data = obs_df, FUN = mean)
  names(obs_rates)[names(obs_rates) == "is.na(R)"] <- "obs_resp_rate"

  # Compute observed mean RTs
  mean_rts <- aggregate(rt ~ subjects + SSD, data = obs_df, FUN = function(x) mean(x, na.rm = TRUE),
                        na.action = na.pass)
  names(mean_rts)[names(mean_rts) == "rt"] <- "obs_mean_rt"

  # Compute posterior predictive response rates per sample
  pred_df$sample_id <- pred_df$postn
  pred_rates <- aggregate(is.na(R) ~ subjects + SSD + sample_id, data = pred_df, FUN = mean,
                          na.action = na.pass)
  names(pred_rates)[names(pred_rates) == "is.na(R)"] <- "pred_resp_rate"

  # Compute posterior predictive mean RTs per sample
  pred_rts <- aggregate(rt ~ subjects + SSD + sample_id, data = pred_df, FUN = function(x) mean(x, na.rm = TRUE),
                        na.action = na.pass)
  names(pred_rts)[names(pred_rts) == "rt"] <- "pred_mean_rt"

  # Ensure consistent types
  obs_rates$subjects <- as.character(obs_rates$subjects)
  obs_rates$SSD <- as.character(obs_rates$SSD)
  mean_rts$subjects <- as.character(mean_rts$subjects)
  mean_rts$SSD <- as.character(mean_rts$SSD)
  pred_rates$subjects <- as.character(pred_rates$subjects)
  pred_rates$SSD <- as.character(pred_rates$SSD)
  pred_rts$subjects <- as.character(pred_rts$subjects)
  pred_rts$SSD <- as.character(pred_rts$SSD)

  # Unique combinations
  unique_combos <- unique(obs_rates[, c("subjects", "SSD")])

  results <- do.call(rbind, lapply(1:nrow(unique_combos), function(i) {
    subj <- unique_combos$subjects[i]
    ssd  <- unique_combos$SSD[i]

    # Observed values
    obs_resp_rate <- obs_rates$obs_resp_rate[obs_rates$subjects == subj & obs_rates$SSD == ssd]
    obs_rt <- mean_rts$obs_mean_rt[mean_rts$subjects == subj & mean_rts$SSD == ssd]

    # Posterior predicted values
    pred_resp_vals <- pred_rates$pred_resp_rate[pred_rates$subjects == subj & pred_rates$SSD == ssd]
    pred_rt_vals <- pred_rts$pred_mean_rt[pred_rts$subjects == subj & pred_rts$SSD == ssd]

    # Misfit calculations
    p_gt_resp <- mean(obs_resp_rate > pred_resp_vals)
    p_lt_resp <- mean(obs_resp_rate < pred_resp_vals)
    misfit_resp <- 2 * min(p_gt_resp, p_lt_resp)

    p_gt_rt <- mean(obs_rt > pred_rt_vals)
    p_lt_rt <- mean(obs_rt < pred_rt_vals)
    misfit_rt <- 2 * min(p_gt_rt, p_lt_rt)

    data.frame(
      subject = subj,
      SSD = ssd,
      obs_resp_rate = obs_resp_rate,
      obs_mean_rt = obs_rt,
      prob_obs_gt_pred_resp = p_gt_resp,
      prob_obs_lt_pred_resp = p_lt_resp,
      misfit_resp = misfit_resp,
      prob_obs_gt_pred_rt = p_gt_rt,
      prob_obs_lt_pred_rt = p_lt_rt,
      misfit_rt = misfit_rt
    )
  }))

  # remove only NA rows

  return(results)
}
