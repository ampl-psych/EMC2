#' Plot Inhibition Functions
#'
#' Plots panels of the inhibition functions (probability of responding; P(R)) for each
#' level of specified factors of stop-signal data as a function of user-defined
#' SSD bins/categories. Optionally, posterior and/or prior predictive
#' inhibition functions can be overlaid.
#'
#' Per default, the SSD-categories are defined in terms of the percentiles of the
#' SSD distribution for each participant, and then averaged over participants (see `use_global_quantiles`).
#'
#' If no data posterior predictives are supplied, the data is plotted with error bars
#' (plus/minus the standard error per SSD bin/category)
#'
#' @param input Either an emc object or a stop-signal data frame, or a *list* of such objects. SSD column in data required.
#' @param probs Numeric vector of probabilities with values in [0,1] that defines SSD bins/categories.
#' @param post_predict Optional posterior predictive data (matching columns) or *list* thereof.
#' @param prior_predict Optional prior predictive data (matching columns) or list thereof.
#' @param factors Character vector of factor names to aggregate over; defaults to plotting full data set ungrouped by factors if NULL.
#' @param within_plot Character indicating factor for which inhibition functions are plotted in the same panel
#' @param use_global_quantiles If set to `TRUE`, SSDs are pooled over participants before calculating percentiles, so
#' the same absolute SSD range is used to get P(R) for each participant,
#' and then these probabilities are averaged over participants.
#' @param subject Subset the data to a single subject (by index or name).
#' @param quants Numeric vector of credible interval bounds (e.g. c(0.025, 0.975)).
#' @param functions A function (or list of functions) that create new columns in the datasets or predictives
#' @param n_cores Number of CPU cores to use if generating predictives from an emc object.
#' @param n_post Number of posterior draws to simulate if needed for predictives.
#' @param layout Numeric vector used in par(mfrow=...); use NA for auto-layout.
#' @param to_plot Character vector: any of "data", "posterior", "prior".
#' @param use_lim Character vector controlling which source(s) define axis limits
#' @param legendpos Character vector controlling the positions of the legends
#' @param posterior_args Optional list of graphical parameters for posterior lines/ribbons.
#' @param prior_args Optional list of graphical parameters for prior lines/ribbons.
#' @param ... Other graphical parameters for the real data lines.
#'
#' @return Returns NULL invisibly
#'
plot_pr <- function(input,
                    probs = seq(0,1, length.out=5),
                    post_predict = NULL,
                    prior_predict = NULL,
                    factors = NULL,
                    within_plot = NULL,
                    use_global_quantiles = FALSE,
                    subject = NULL,
                    quants = c(0.025, 0.975), functions = NULL,
                    n_cores = 1,
                    n_post = 50,
                    layout = NA,
                    to_plot = c('data','posterior','prior')[1:2],
                    use_lim = c('data','posterior','prior')[1:2],
                    legendpos = c('top', 'topright'),
                    posterior_args = list(),
                    prior_args = list(),
                    ...) {

  # 1) prep_data_plot
  check <- prep_data_plot(input, post_predict, prior_predict, to_plot, use_lim,
                          factors, within_plot, subject, n_cores, n_post, remove_na = FALSE, functions)
  data_sources <- check$datasets
  sources <- check$sources

  # Basic definitions
  dots <- add_defaults(list(...), col = c("black",  "#A9A9A9", "#666666"))
  posterior_args <- add_defaults(posterior_args, col = c("darkgreen",  "#0000FF", "#008B8B"))
  prior_args <- add_defaults(prior_args, col = c("red", "#800080", "#CC00FF"))

  unique_group_keys <- levels(factor(data_sources[[1]]$group_key))

  if (is.null(within_plot)) {
    within_levels <- "All"  # fallback
  } else {
    within_levels <- levels(factor(data_sources[[1]][[within_plot]]))

  }
  line_types <- seq_along(within_levels)

  # Single-p_resp results or multi-postn quantile results are stored in these:
  p_resp_list        <- list() # single
  p_resp_quants_list <- list() # multi

  # keep track of a global maximum in the vertical dimension
  # so that we can set a consistent y-lim across all panels
  y_max <- 1
  y_min <- 0

  x_min <- 0
  x_max <- 1

  # -------------------------------------------------------------------
  # 2) FIRST BIG LOOP: compute p_resp or p_resp-quantiles for each dataset
  # -------------------------------------------------------------------
  for (k in seq_along(data_sources)) {
    df   <- data_sources[[k]]
    styp <- sources[k]       # "data","posterior","prior"
    sname <- names(sources)[k]  # the name of this dataset in the list

    if (is.null(df) || !nrow(df)) {
      # skip empty
      next
    }

    # group by group_key
    splitted <- split(df, df$group_key)

    # Check if stop-signal dataset
    if(is.null(df$SSD)){
      stop("No SSD column in data.")
    }

    # set type of SSD binning method
    if(use_global_quantiles){
      quantile_fun <- get_response_probability_by_global_ssd_quantile
      global_ssd_pool <- df$SSD[is.finite(df$SSD)]
      global_ssd_breaks <- quantile(global_ssd_pool, probs = probs, na.rm = TRUE)
      dots$global_ssd_breaks <- global_ssd_breaks
    } else {
      quantile_fun <- get_response_probability_by_individual_ssd_quantile
    }

    # If there's a "postn" column => multiple draws => compute quantiles
    if ("postn" %in% names(df)) {
      # p_resp_list[[sname]] => list of (group_key => list of postn => single-p_resp)
      p_resp_list[[sname]] <- lapply(splitted, function(sub_grp) {
        #  further split by postn
        postn_splits <- split(sub_grp, sub_grp$postn)
        lapply(postn_splits, quantile_fun,
               group_factor = within_plot, probs = probs, dots = dots)
      })

      # derive p_resp_quants_list from p_resp_list
      p_resp_quants_list[[sname]] <- lapply(p_resp_list[[sname]], function(postn_list) {
        # postn_list => e.g. 100 draws => each draw is a named list of factor-level => cbind(x,y)
        out <- list()
        for (lev in within_levels) {
          # gather all x and y columns across draws
          x_mat <- do.call(cbind, lapply(postn_list, function(lst) lst[[lev]][,"x_plot"]))
          y_mat <- do.call(cbind, lapply(postn_list, function(lst) lst[[lev]][,"p_response"]))

          # Combine them in a matrix with 4 rows => y_lower, y_upper, y_median, x_plot
          qy <- apply(y_mat, 1, quantile, probs = sort(c(quants, 0.5)), na.rm = TRUE)
          xm <- apply(x_mat, 1, unique, na.rm=TRUE)
          out[[lev]] <- rbind(qy, xm)
        }
        out
      })

      # If this dataset is used to define y-limit
      if (styp %in% use_lim) {
        # maximum of x_plot row
        possible_vals <- unlist(lapply(p_resp_quants_list[[sname]], function(group_val) {
          # group_val => list( factor_level => matrix(4 x length-of-probs) )
          sapply(group_val, function(mat4) {
            if (is.null(mat4)) return(0)
            max(mat4[nrow(mat4), ], na.rm=TRUE)
          })
        }))
        x_max <- max(x_max, possible_vals, na.rm=TRUE)

        # y-limits (min/max across y_lower and y_upper)
        y_vals <- unlist(lapply(p_resp_quants_list[[sname]], function(group_val) {
          lapply(group_val, function(mat4) {
            if (is.null(mat4)) return(NULL)
            range(c(mat4[1, ], mat4[3, ]), na.rm = TRUE)  # y_lower and y_upper
          })
        }))
        y_min <- max(y_min, y_vals, na.rm = TRUE)
        y_max <- max(y_max, y_vals, na.rm = TRUE)
      }

    } else {

      # single dataset => p_resp_list[[sname]] => group_key => get_def_cdf => named list by factor level
      p_resp_list[[sname]] <- lapply(splitted, quantile_fun,
                                     group_factor = within_plot, probs = probs, dots = dots)

      # If this dataset is used for y-limit, find max
      if (styp %in% use_lim) {
        # find max y across all group_keys & factor-levels
        all_x_vals <- c()
        all_y_vals <- c()

        for (grp_name in names(p_resp_list[[sname]])) {
          p_resp_grp <- p_resp_list[[sname]][[grp_name]]

          if (!is.null(p_resp_grp)) {
            for (draw in p_resp_grp) {
              # extract x_plot and p_response if present
              if (!is.null(draw[["x_plot"]]) && !is.null(draw[["p_response"]])) {
                all_x_vals <- c(all_x_vals, draw$x_plot)
                all_y_vals <- c(all_y_vals, draw$p_response)
              }
            }
          }
        }

        # Now calculate global min/max
        x_min <- min(all_x_vals, na.rm = TRUE)
        x_max <- max(all_x_vals, na.rm = TRUE)
        y_max <- max(all_y_vals, na.rm = TRUE)
      }
    }
  }

  # -------------------------------------------------------------------
  # 3) SECOND BIG LOOP: Plot one panel per group_key
  # -------------------------------------------------------------------
  if(!is.null(layout)){
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }
  # layout
  if (any(is.na(layout))) {
    par(mfrow = coda_setmfrow(Nchains=1, Nparms=length(unique_group_keys), nplots=1))
  } else {
    par(mfrow = layout)
  }

  # define a global y-limit (with a bit of headroom)
  if (!is.finite(y_max) || y_max <= 0) y_max <- 1
  ylim <- c(0, y_max*1.5)

  if (!is.finite(x_min) || !is.finite(x_max) || x_min >= x_max) {
    xlim <- NULL  # let R handle it if bad limits
  } else {
    x_buffer <- 0.05 * (x_max - x_min)
    xlim <- c(x_min - x_buffer, x_max + x_buffer)
  }

  for (group_key in unique_group_keys) {
    tmp_dots <- dots
    tmp_posterior_args <- posterior_args
    tmp_prior_args <- prior_args

    # blank plot
    plot_args <- add_defaults(dots, ylim=ylim, xlim=xlim,
                              main=group_key, xlab="", ylab="")
    plot_args <- fix_dots_plot(plot_args)
    plot_args$axes <- FALSE  # prevent auto-drawing axes
    do.call(plot, c(list(NA), plot_args))
    # Draw y-axis manually (left side)
    axis(2)

    mtext("P(R)", side=2, line=3)

    # draw lines for each dataset
    legend_map <- character(0)  # to store source name -> color
    lwd_map <- numeric()
    for (k in seq_along(data_sources)) {
      styp  <- sources[k]       # "data","posterior","prior"
      sname <- names(sources)[k]
      if(styp == "data") {
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
      # if no postn => single p_resp => p_resp_list[[sname]][[group_key]] => factor-level => matrix(x,y)
      # if postn => p_resp_quants_list[[sname]][[group_key]] => factor-level => matrix(4 x length-probs)

      # There might be no p_resp if that group_key wasn't present in the data
      # so check if p_resp_list has an entry
      if (!is.null(p_resp_list[[sname]])) {
        # check for group group_key
        if (!("postn" %in% names(data_sources[[k]]))) {
          # single
          p_resp_df <- p_resp_list[[sname]][[group_key]]
          if (!is.null(p_resp_df)) {
            ilev <- 1
            for (df in p_resp_df) {
              lines_args <- add_defaults(src_args, lty = line_types[ilev])
              lines_args <- fix_dots_plot(lines_args)

              do.call(lines, c(list(x = df$x_plot, y = df$p_response), lines_args))
              do.call(points, c(list(x = df$x_plot, y = df$p_response), lines_args))

              if(is.null(post_predict)){
                # Add confidence intervals
                arrows(x0 = df$x_plot, y0 = df$p_response - df$se,
                       x1 = df$x_plot, y1 = df$p_response + df$se,
                       angle = 90, code = 3, length = 0.05)
              }


              ilev <- ilev + 1
            }
          }

          # Add x-axis label and tick labels (after drawing data points)
          if (styp == "data" && !is.null(p_resp_list[[sname]][[group_key]])) {
            first_level <- names(p_resp_list[[sname]][[group_key]])[1]
            tick_data <- p_resp_list[[sname]][[group_key]][[first_level]]

            if (!is.null(tick_data) && "x_plot" %in% names(tick_data)) {
              if (use_global_quantiles) {
                # Global bins → label with SSD values
                axis(1,
                     at = tick_data$x_plot,
                     labels = tick_data$ssd,
                     cex.axis = 0.8)
                title(xlab = "SSD Bin (global SSD values)", line = 2.5)
              } else {
                # Subject-specific quantiles → label with percentages
                quantile_labels <- paste0(round(tick_data$x_plot * 100), "%")
                axis(1,
                     at = tick_data$x_plot,
                     labels = quantile_labels,
                     cex.axis = 0.8)
                title(xlab = "SSD Quantile (subject-specific)", line = 2.5)
              }
            }
          }

        } else {
          # multi draws => p_resp_quants_list
          if (!is.null(p_resp_quants_list[[sname]])) {
            p_resp_quants_for_group <- p_resp_quants_list[[sname]][[group_key]]
            if (!is.null(p_resp_quants_for_group)) {
              ilev <- 1
              for (lev in within_levels) {
                mat4 <- p_resp_quants_for_group[[lev]]
                if (!is.null(mat4)) {
                  # mat4 => e.g. 4 rows x (length(probs)) columns
                  y_lower <- mat4[1,]
                  y_med   <- mat4[2,]
                  y_upper <- mat4[3,]
                  x_plot<- mat4[4,]

                  lines_args <- add_defaults(src_args, lty=line_types[ilev])
                  lines_args <- fix_dots_plot(lines_args)
                  do.call(lines, c(list(x=x_plot, y=y_med), lines_args))

                  # polygon for the ribbon
                  adj_color <- do.call(adjustcolor, fix_dots(add_defaults(src_args, alpha.f=0.2), adjustcolor))
                  poly_args <- src_args
                  poly_args$col <- adj_color
                  poly_args <- fix_dots_plot(poly_args)
                  do.call(polygon, c(list(
                    y = c(y_lower, rev(y_upper)),
                    x = c(x_plot, rev(x_plot)),
                    border = NA
                  ), poly_args))
                }
                ilev <- ilev+1
              }
            }
          }
        }
      }
    }

    # Factor-level legend, only if more than one defective level
    if (!is.na(legendpos[1]) && !(length(within_levels) == 1 && within_levels == "All")) {
      legend(legendpos[1], legend=within_levels, lty=line_types, col="black",
             title=within_plot, bty="n")
    }


    # If multiple data sources, show source legend
    if (length(data_sources) > 1) {
      if(!is.na(legendpos[2])){
        legend(legendpos[2], legend=names(legend_map), lty=1, col=legend_map,
               title="Source", bty="n", lwd = lwd_map)
      }
    }
  } # end for each group_key

  invisible(NULL)
}

get_response_probability_by_individual_ssd_quantile <- function(x, group_factor, probs, dots) {
  # Filter: only rows with finite SSD and subject info
  x <- x[is.finite(x[["SSD"]]) & !is.na(x[["subjects"]]), ]
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
      df <- data_subset[data_subset$subject == s, ]

      ssd_quants <- quantile(df$SSD, probs = probs, na.rm = TRUE, type = 1)
      if (anyDuplicated(ssd_quants)) {
        stop("Duplicate quantile values detected. Please use fewer bins.")
      }
      ssd_quants <- unique(ssd_quants)

      df$ssd_bin <- cut(df$SSD, breaks = ssd_quants, include.lowest = TRUE, right = TRUE)

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
      bin_stats$x_plot <- probs[-1]

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
    } else {
      summary_stats <- subj_stats[[1]]
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
  bin_labels <- levels(cut(x$SSD, breaks = global_ssd_breaks, include.lowest = TRUE, right = TRUE))

  # Function to compute stats for a subset of data
  compute_group_stats <- function(data_subset) {
    subj_stats <- lapply(unique(data_subset$subjects), function(s) {
      df <- data_subset[data_subset$subjects == s, ]

      df$ssd_bin <- cut(df$SSD, breaks = global_ssd_breaks, include.lowest = TRUE, right = TRUE)

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
                                   se_p <- sqrt(mean_p * (1 - mean_p) / n)
                                   c(mean = mean_p, se = se_p, n = n)
                                 },
                                 na.action = na.pass)
      summary_stats <- do.call(data.frame, summary_stats)

      # Add x_plot by matching back from all_stats
      x_plot_lookup <- aggregate(x_plot ~ ssd, data = all_stats, FUN = function(x) unique(x)[1])
      summary_stats <- merge(summary_stats, x_plot_lookup, by = "ssd")

      summary_stats <- summary_stats[order(summary_stats$x_plot), ]
      names(summary_stats) <- c("ssd", "p_response", "se", "n", "x_plot")
    } else {
      summary_stats <- subj_stats[[1]]
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

#' Plot Mean SRRT
#'
#' Plots panels of the mean signal-respond response time (SRRT) as a function of user-defined
#' SSD bins/categories for each level of specified factors of stop-signal data.
#' Optionally, posterior and/or prior predictive inhibition functions can be overlaid.
#'
#' Per default, SSDs are pooled over participants before calculating percentiles, so
#' the same absolute SSD range is used to get mean SRRT for each participant,
#' and then these probabilities are averaged over participants (see `use_global_quantiles`).
#'
#' If no data posterior predictives are supplied, the data is plotted with error bars
#' (plus/minus the standard error per SSD bin/category)
#'
#' @param input Either an emc object or a stop-signal data frame, or a list of such objects. SSD column in data required.
#' @param probs Numeric vector of probabilities with values in 0,1 that defines SSD bins/categories.
#' @param factors Character vector of factor names to aggregate over; defaults to plotting full data set ungrouped by factors if NULL.
#' @param within_plot Character indicating factor for which inhibition functions are plotted in the same panel
#' @param use_global_quantiles If set to FALSE, the SSD-categories are defined in terms of the percentiles of the
#' SSD distribution for each participant, and then averaged over participants.
#' @param post_predict Optional posterior predictive data (matching columns) or list thereof.
#' @param prior_predict Optional prior predictive data (matching columns) or list thereof.
#' @param subject Subset the data to a single subject (by index or name).
#' @param quants Numeric vector of credible interval bounds (e.g. c(0.025, 0.975)).
#' @param functions A function (or list of functions) that create new columns in the datasets or predictives
#' @param n_cores Number of CPU cores to use if generating predictives from an emc object.
#' @param n_post Number of posterior draws to simulate if needed for predictives.
#' @param layout Numeric vector used in par(mfrow=...); use NA for auto-layout.
#' @param to_plot Character vector: any of "data", "posterior", "prior".
#' @param use_lim Character vector controlling which source(s) define axis limits
#' @param legendpos Character vector controlling the positions of the legends
#' @param posterior_args Optional list of graphical parameters for posterior lines/ribbons.
#' @param prior_args Optional list of graphical parameters for prior lines/ribbons.
#' @param ... Other graphical parameters for the real data lines.
#'
#' @return Returns NULL invisibly
#'
plot_srrt <- function(input,
                      probs = seq(0,1, length.out=5),
                      factors = NULL,
                      within_plot = NULL,
                      use_global_quantiles = FALSE,
                      post_predict = NULL,
                      prior_predict = NULL,
                      subject = NULL,
                      quants = c(0.025, 0.975), functions = NULL,
                      n_cores = 1,
                      n_post = 50,
                      layout = NA,
                      to_plot = c('data','posterior','prior')[1:2],
                      use_lim = c('data','posterior','prior')[1:2],
                      legendpos = c('top', 'topright'),
                      posterior_args = list(),
                      prior_args = list(),
                      ...) {

  # 1) prep_data_plot
  check <- prep_data_plot(input, post_predict, prior_predict, to_plot, use_lim,
                          factors, within_plot, subject, n_cores, n_post, remove_na = FALSE, functions)
  data_sources <- check$datasets
  sources <- check$sources

  # Basic definitions
  dots <- add_defaults(list(...), col = c("black",  "#A9A9A9", "#666666"))
  posterior_args <- add_defaults(posterior_args, col = c("darkgreen",  "#0000FF", "#008B8B"))
  prior_args <- add_defaults(prior_args, col = c("red", "#800080", "#CC00FF"))

  unique_group_keys <- levels(factor(data_sources[[1]]$group_key))

  if (is.null(within_plot)) {
    within_levels <- "All"  # fallback
  } else {
    within_levels <- levels(factor(data_sources[[1]][[within_plot]]))

  }
  line_types <- seq_along(within_levels)

  # store single-p_resp results or multi-postn quantile results in these:
  p_resp_list        <- list() # single
  p_resp_quants_list <- list() # multi

  # keep track of a global maximum in the vertical dimension
  # so a consistent y-lim across all panels can be set
  y_max <- 1
  y_min <- 0

  x_min <- 0
  x_max <- 1

  # -------------------------------------------------------------------
  # 2) FIRST BIG LOOP: compute p_resp or p_resp-quantiles for each dataset
  # -------------------------------------------------------------------
  for (k in seq_along(data_sources)) {
    df   <- data_sources[[k]]
    styp <- sources[k]       # "data","posterior","prior"
    sname <- names(sources)[k]  # the name of this dataset in the list

    if (is.null(df) || !nrow(df)) {
      # skip empty
      next
    }

    # group by group_key
    splitted <- split(df, df$group_key)

    # Check if it is a stop-signal dataset
    if(is.null(df$SSD)){
      stop("No SSD column in data.")
    }

    # type of SSD binning method
    if(use_global_quantiles){
      quantile_fun <- get_srrt_by_global_ssd_quantile
      global_ssd_pool <- df$SSD[is.finite(df$SSD)]
      global_ssd_breaks <- quantile(global_ssd_pool, probs = probs, na.rm = TRUE)
      dots$global_ssd_breaks <- global_ssd_breaks
    } else {
      quantile_fun <- get_srrt_by_individual_ssd_quantile
    }

    # If there's a "postn" column => multiple draws => compute quantiles
    if ("postn" %in% names(df)) {
      # p_resp_list[[sname]] => list of (group_key => list of postn => single-p_resp)
      p_resp_list[[sname]] <- lapply(splitted, function(sub_grp) {
        # sub_grp is all rows for a single group_key
        # further split by postn
        postn_splits <- split(sub_grp, sub_grp$postn)
        lapply(postn_splits, quantile_fun,
               group_factor = within_plot, probs = probs, dots = dots)
      })

      # derive p_resp_quants_list from p_resp_list
      p_resp_quants_list[[sname]] <- lapply(p_resp_list[[sname]], function(postn_list) {
        # postn_list => e.g. 100 draws => each draw is a named list of factor-level => cbind(x,y)
        out <- list()
        for (lev in within_levels) {
          # gather all x and y columns across draws
          x_mat <- do.call(cbind, lapply(postn_list, function(lst) lst[[lev]][,"x_plot"]))
          y_mat <- do.call(cbind, lapply(postn_list, function(lst) lst[[lev]][,"srrt"]))

          qy <- apply(y_mat, 1, quantile, probs = sort(c(quants, 0.5)), na.rm = TRUE)
          xm <- apply(x_mat, 1, unique, na.rm=TRUE)
          out[[lev]] <- rbind(qy, xm)
        }
        out
      })

      # If this dataset is used to define y-limit
      if (styp %in% use_lim) {
        possible_vals <- unlist(lapply(p_resp_quants_list[[sname]], function(group_val) {
          # group_val => list( factor_level => matrix(4 x length-of-probs) )
          sapply(group_val, function(mat4) {
            if (is.null(mat4)) return(0)
            max(mat4[nrow(mat4), ], na.rm=TRUE)
          })
        }))
        x_max <- max(x_max, possible_vals, na.rm=TRUE)

        # y-limits (min/max across y_lower and y_upper)
        y_vals <- unlist(lapply(p_resp_quants_list[[sname]], function(group_val) {
          lapply(group_val, function(mat4) {
            if (is.null(mat4)) return(NULL)
            range(c(mat4[1, ], mat4[3, ]), na.rm = TRUE)  # y_lower and y_upper
          })
        }))
        y_min <- min(y_min, y_vals, na.rm = TRUE)
        y_max <- max(y_max, y_vals, na.rm = TRUE)
      }

    } else {

      # single dataset => p_resp_list[[sname]] => group_key => get_def_cdf => named list by factor level
      p_resp_list[[sname]] <- lapply(splitted, quantile_fun,
                                     group_factor = within_plot, probs = probs, dots = dots)

      # If this dataset for y-limit, find max
      if (styp %in% use_lim) {
        # find max y across all group_keys & factor-levels
        all_x_vals <- c()
        all_y_vals <- c()

        for (grp_name in names(p_resp_list[[sname]])) {
          p_resp_grp <- p_resp_list[[sname]][[grp_name]]

          if (!is.null(p_resp_grp)) {
            for (draw in p_resp_grp) {
              # extract x_plot and srrt if present
              if (!is.null(draw[["x_plot"]]) && !is.null(draw[["srrt"]])) {
                all_x_vals <- c(all_x_vals, draw$x_plot)
                all_y_vals <- c(all_y_vals, draw$srrt)
              }
            }
          }
        }

        # calculate global min/max
        x_min <- min(all_x_vals, na.rm = TRUE)
        x_max <- max(all_x_vals, na.rm = TRUE)
        y_max <- max(all_y_vals, na.rm = TRUE)
      }
    }
  }

  # -------------------------------------------------------------------
  # 3) SECOND BIG LOOP: Plot one panel per group_key
  # -------------------------------------------------------------------
  if(!is.null(layout)){
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }
  # layout
  if (any(is.na(layout))) {
    par(mfrow = coda_setmfrow(Nchains=1, Nparms=length(unique_group_keys), nplots=1))
  } else {
    par(mfrow = layout)
  }

  # define a global y-limit (with a bit of headroom)
  if (!is.finite(y_max) || y_max <= 0) y_max <- 1
  ylim <- c(0, y_max*1.1)

  if (!is.finite(x_min) || !is.finite(x_max) || x_min >= x_max) {
    xlim <- NULL  # let R handle it if bad limits
  } else {
    x_buffer <- 0.05 * (x_max - x_min)
    xlim <- c(x_min - x_buffer, x_max + x_buffer)
  }

  for (group_key in unique_group_keys) {
    tmp_dots <- dots
    tmp_posterior_args <- posterior_args
    tmp_prior_args <- prior_args

    # blank plot
    plot_args <- add_defaults(dots, ylim=ylim, xlim=xlim,
                              main=group_key, xlab="", ylab="")
    plot_args <- fix_dots_plot(plot_args)
    plot_args$axes <- FALSE  # prevent auto-drawing axes
    do.call(plot, c(list(NA), plot_args))
    # Draw y-axis manually (left side)
    axis(2)
    mtext("Mean SRRT", side=2, line=3)

    # draw lines for each dataset
    legend_map <- character(0)  # to store source name -> color
    lwd_map <- numeric()
    for (k in seq_along(data_sources)) {
      styp  <- sources[k]       # "data","posterior","prior"
      sname <- names(sources)[k]
      if(styp == "data") {
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
      # if no postn => single p_resp => p_resp_list[[sname]][[group_key]] => factor-level => matrix(x,y)
      # if postn => p_resp_quants_list[[sname]][[group_key]] => factor-level => matrix(4 x length-probs)

      # THere might not be p_resp if that group_key wasn't present in the data
      # so check if p_resp_list has an entry
      if (!is.null(p_resp_list[[sname]])) {
        # check for group group_key
        if (!("postn" %in% names(data_sources[[k]]))) {
          # single
          p_resp_df <- p_resp_list[[sname]][[group_key]]
          if (!is.null(p_resp_df)) {
            ilev <- 1
            for (df in p_resp_df) {
              lines_args <- add_defaults(src_args, lty = line_types[ilev])
              lines_args <- fix_dots_plot(lines_args)

              do.call(lines, c(list(x = df$x_plot, y = df$srrt), lines_args))
              do.call(points, c(list(x = df$x_plot, y = df$srrt), lines_args))

              if(is.null(post_predict)){
                # Add confidence intervals (error bars)
                arrows(x0 = df$x_plot, y0 = df$srrt - df$se,
                       x1 = df$x_plot, y1 = df$srrt + df$se,
                       angle = 90, code = 3, length = 0.05)
              }


              ilev <- ilev + 1
            }
          }

          # Add x-axis label and tick labels (after drawing data points)
          if (styp == "data" && !is.null(p_resp_list[[sname]][[group_key]])) {
            first_level <- names(p_resp_list[[sname]][[group_key]])[1]
            tick_data <- p_resp_list[[sname]][[group_key]][[first_level]]

            if (!is.null(tick_data) && "x_plot" %in% names(tick_data)) {
              if (use_global_quantiles) {
                # Global bins --> label with SSD values
                axis(1,
                     at = tick_data$x_plot,
                     labels = tick_data$ssd,
                     cex.axis = 0.8)
                title(xlab = "SSD Bin (global SSD values)", line = 2.5)
              } else {
                # Subject-specific quantiles --> label with percentages
                quantile_labels <- paste0(round(tick_data$x_plot * 100), "%")
                axis(1,
                     at = tick_data$x_plot,
                     labels = quantile_labels,
                     cex.axis = 0.8)
                title(xlab = "SSD Quantile (subject-specific)", line = 2.5)
              }
            }
          }

        } else {
          # multi draws => p_resp_quants_list
          if (!is.null(p_resp_quants_list[[sname]])) {
            p_resp_quants_for_group <- p_resp_quants_list[[sname]][[group_key]]
            if (!is.null(p_resp_quants_for_group)) {
              ilev <- 1
              for (lev in within_levels) {
                mat4 <- p_resp_quants_for_group[[lev]]
                if (!is.null(mat4)) {
                  y_lower <- mat4[1,]
                  y_med   <- mat4[2,]
                  y_upper <- mat4[3,]
                  x_plot<- mat4[4,]

                  lines_args <- add_defaults(src_args, lty=line_types[ilev])
                  lines_args <- fix_dots_plot(lines_args)
                  do.call(lines, c(list(x=x_plot, y=y_med), lines_args))

                  # polygon for the ribbon
                  adj_color <- do.call(adjustcolor, fix_dots(add_defaults(src_args, alpha.f=0.2), adjustcolor))
                  poly_args <- src_args
                  poly_args$col <- adj_color
                  poly_args <- fix_dots_plot(poly_args)
                  do.call(polygon, c(list(
                    y = c(y_lower, rev(y_upper)),
                    x = c(x_plot, rev(x_plot)),
                    border = NA
                  ), poly_args))
                }
                ilev <- ilev+1
              }
            }
          }
        }
      }
    }

    # Factor-level legend, only if more than one defective level
    if (!is.na(legendpos[1]) && !(length(within_levels) == 1 && within_levels == "All")) {
      legend(legendpos[1], legend=within_levels, lty=line_types, col="black",
             title=within_plot, bty="n")
    }


    # If multiple data sources, show source legend
    if (length(data_sources) > 1) {
      if(!is.na(legendpos[2])){
        legend(legendpos[2], legend=names(legend_map), lty=1, col=legend_map,
               title="Source", bty="n", lwd = lwd_map)
      }
    }
  } # end for each group_key

  invisible(NULL)
}

get_srrt_by_individual_ssd_quantile <- function(x, group_factor, probs, dots) {
  # Filter: only rows with finite SSD and subject info
  x <- x[is.finite(x[["SSD"]]) & !is.na(x[["subjects"]]), ]
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
      df <- data_subset[data_subset$subject == s, ]

      ssd_quants <- quantile(df$SSD, probs = probs, na.rm = TRUE, type = 1)
      if (anyDuplicated(ssd_quants)) {
        stop("Duplicate quantile values detected. Please use fewer bins.")
      }
      ssd_quants <- unique(ssd_quants)

      df$ssd_bin <- cut(df$SSD, breaks = ssd_quants, include.lowest = TRUE, right = TRUE)

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
      bin_stats$x_plot <- probs[-1]

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
    } else {
      summary_stats <- subj_stats[[1]][, c("x_plot", "srrt", "ssd", "se")]
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
  bin_labels <- levels(cut(x$SSD, breaks = global_ssd_breaks, include.lowest = TRUE, right = TRUE))

  # Function to compute stats for a subset of data
  compute_group_stats <- function(data_subset) {
    subj_stats <- lapply(unique(data_subset$subjects), function(s) {
      df <- data_subset[data_subset$subjects == s, ]

      df$ssd_bin <- cut(df$SSD, breaks = global_ssd_breaks, include.lowest = TRUE, right = TRUE)


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
    } else {
      summary_stats <- subj_stats[[1]][, c("x_plot", "srrt", "ssd", "se")]
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
