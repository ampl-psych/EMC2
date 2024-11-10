create_group_key <- function(df, factors) {
  if (length(factors) == 0) return(rep("All Data", nrow(df)))
  apply(df[, factors, drop = FALSE], 1, function(x) paste(paste(factors, x, sep = "="), collapse = " "))
}

#' Plot defective cumulative distribution functions
#'
#' Plots panels that contain a set of cumulative distribution functions (CDFs) for each level of the specified defective factor in the data.
#' These CDFs are defective; their maximum values are relative to the respective
#' proportions of the defective factor levels. Across all levels, the maximum value sums to 1.
#' Optionally, posterior predictive CDFs can be overlaid.
#'
#' @param input. Either an emc object or a data frame. The data frame must be in EMC2 format with at least subjects, rt, and the defective factor columns.
#' @param subject An integer or character string selecting a subject from the data.
#' If specified, only that subject is plotted. By default (NULL), all subjects are plotted.
#' @param factors A character vector of the factor names in the design to aggregate across.
#' Defaults to all factors (i.e., NULL).
#' @param defective_factor A character string specifying the factor to compute defective CDFs for. Defaults to "R".
#' @param layout A vector indicating which layout to use as in par(mfrow = layout). If NA, will automatically generate an appropriate layout.
#' @param post_predict A data frame containing posterior predictive datasets. Should have the same columns as data, with an additional postn column indicating the index of the predicted dataset.
#' @param quants A numeric vector specifying the quantiles for the credible intervals (defaults to c(0.025, 0.975)).
#' @param n_cores An integer specifying the number of cores with which to calculate predictives if specified.
#' @param n_post An optional integer specifying how many posterior predictives to generate if the input is an emc object.
#' @param predicted_plot_args A list of plotting arguments to be used for the predicted data (e.g., col, lty).
#' @param ... Optional arguments that can be passed to plotting functions for the real data.
#' @return Invisibly returns a data frame containing RT quantiles and mean proportions per factor combination and defective factor level.
#' If post_predict is provided or the input is an emc object, predicted values are included.
#' @examples
#' # Plot defective cdf for each subject and the factor combination in the design:
#' plot_cdf(forstmann)
#' # or for one subject:
#' plot_cdf(forstmann, subject = 1)
#' # Now collapsing across subjects and using a different defective factor:
#' plot_cdf(forstmann, factors = "S", defective_factor = "E")
#' # Or plot posterior predictives with an emc object
#' plot_cdf(samples_LNR, n_post = 10)
#' @export
plot_cdf <- function(input, subject = NULL, factors = NULL,
                     defective_factor = "R",
                     layout = NA, post_predict = NULL, quants = c(0.025, 0.975),
                     n_cores = 1, n_post = 100, predict_args = list(),
                     ...) {
  dots <- add_defaults(list(...), col = "black")
  check <- prep_data_plot(input, post_predict, factors, defective_factor, subject, n_cores = n_cores, n_post = n_post)
  data <- check$data
  post_predict <- check$post_predict
  # Get defective_levels
  defective_levels <- levels(factor(data[[defective_factor]]))
  line_types <- seq_along(defective_levels)

  # Initialize variables
  y_max <- 0
  quantiles <- sort(c(quants, 0.5))
  # Compute xlim based on quantiles
  x_lim_probs <- c(0.001, 0.995)
  dat_quants <- aggregate(rt ~ group_key, data, quantile, probs = x_lim_probs)
  xlim_data <- range(unlist(dat_quants$rt))
  if (!is.null(post_predict)) {
    pred_quants <- aggregate(rt ~ group_key, post_predict, quantile, probs = x_lim_probs)
    xlim_pred <- range(unlist(pred_quants$rt))
    xlim <- range(xlim_data, xlim_pred)
  } else {
    xlim <- xlim_data
  }
  predict_args <- add_defaults(predict_args, col = "red")
  # Compute densities and y_max
  dens_data <- lapply(split(data, data$group_key), get_def_cdf, defective_factor, dots)
  y_max <- max(sapply(unlist(dens_data, recursive = FALSE), function(x) return(max(x[,'y']))))
  # Process post_predict densities
  if (!is.null(post_predict)) {
    dens_pred <- lapply(split(post_predict, post_predict$group_key), FUN = function(x){
      lapply(split(x, x$postn), get_def_cdf, defective_factor, predict_args)
    })
    quants <- lapply(dens_pred, function(x){
      out <- list()
      for(deflev in defective_levels){
        tmp <- apply(do.call(cbind, lapply(x, function(q) return(q[[deflev]][,'x']))), 1, quantile, probs = quantiles)
        out[[deflev]] <- rbind(tmp, apply(do.call(cbind, lapply(x, function(q) return(q[[deflev]][,'y']))), 1, median))
      }
      return(out)
    })
    y_max <- max(y_max, sapply(unlist(quants, FALSE), function(x) return(max(x[4,]))))
  }

  # Set global y-limits
  ylim_global <- c(0, y_max*1.1)

  # Set up plotting layout
  unique_group_keys <- unique(data$group_key)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  if (any(is.na(layout))) {
    par(mfrow = coda_setmfrow(Nchains = 1, Nparms = length(unique_group_keys), nplots = 1))
  } else {
    par(mfrow = layout)
  }

  # Plotting loop
  for (group_key in unique_group_keys) {
    # Set plot arguments with consistent y-limits
    plot_args <- add_defaults(dots, xlim = xlim, ylim = ylim_global, main = group_key, xlab = "RT", ylab = "Defective Density")
    plot_args <- fix_dots_plot(plot_args)
    do.call(plot, c(list(NA), plot_args))
    cur_data <- dens_data[[group_key]]
    # Plot observed densities
    for (i in seq_along(defective_levels)) {
      level <- defective_levels[i]
      if (!is.null(cur_data[[level]])) {
        lines_args <- add_defaults(dots, lty = line_types[i])
        lines_args <- fix_dots_plot(lines_args)
        do.call(lines, c(list(x = cur_data[[level]][,'x'], y = cur_data[[level]][,'y']), lines_args))
      }
    }

    # Plot predicted densities
    if (!is.null(post_predict)) {
      cur_quants <- quants[[group_key]]
      for (i in seq_along(defective_levels)) {
        level <- defective_levels[i]
        if (!is.null(cur_quants[[level]])) {
          lines_args <- add_defaults(predict_args, lty = line_types[i])
          lines_args <- fix_dots_plot(lines_args)
          do.call(lines, c(list(x = cur_quants[[level]][2,], y = cur_quants[[level]][4,]), lines_args))
          # Add ribbon
          adj_color <- do.call(adjustcolor, fix_dots(add_defaults(predict_args, alpha.f = .2), adjustcolor))
          polygon_args <- predict_args
          polygon_args$col <- adj_color
          do.call(polygon, c(list(y = c(cur_quants[[level]][4,], rev(cur_quants[[level]][4,])), border = NA,
                                  x = c(cur_quants[[level]][1,], rev(cur_quants[[level]][3,]))), fix_dots_plot(polygon_args)))
        }
      }
    }

    # Add legends
    legend(x = "topright", legend = defective_levels, lty = line_types, col = "black", title = defective_factor, bty = "n")
    if (!is.null(post_predict)) {
      legend(x = "top", legend = c("Real", "Pred"), lty = 1, col = c(dots$col, predict_args$col), title = "Source", bty = "n")
    }
  }
}


#' Plot statistics on data
#'
#' Plots panels that contain a set of densities for each level of the specified `factor`
#' The densities represent the predicted data across the posterior, the vertical lines represent the real data.
#'
#' @inheritParams plot_cdf
#' @param stat_fun A function that can be applied to the data and returns a single value or a vector of values.
#' @param stat_name The name of the calculated quantity
#' @param ...
#'
#' @return an invisible data frame with the stat applied to the real data, posterior predictives and/or prior predictives
#' @export
#'
#' @examples
#' # For example plot the observed and predicted response accuracy
#' correct_fun <- function(data) mean(data$S == data$R)
#' plot_stat(samples_LNR, correct_fun, n_post = 10)
#' # Can also apply more sophisticated statistics
#' drt <- function(data) diff(tapply(data$rt,data[,c("E")],mean))
#' plot_stat(samples_LNR, drt, n_post = 10, stat_name = "RT diff Speed - A/N")
#'
plot_stat <- function(input, stat_fun, stat_name = NULL, factors = NULL, post_predict = NULL,
                      quants = c(0.025, 0.5, 0.975), do_plot = TRUE, layout = NA,
                      n_post = 100, n_cores = 1, predict_args = list(), ...) {
  # Prepare data
  if (inherits(input, "emc")) {
    data <- get_data(input)
    if (is.null(post_predict)) {
      post_predict <- predict(input, n_post = n_post, n_cores = n_cores)
    }
  } else {
    data <- input
  }

  if (!is.null(factors) && !all(factors %in% names(data))) {
    stop("Factors must be columns in data")
  }

  # Add group keys to data
  data$group_key <- create_group_key(data, factors)
  split_data <- split(data, data$group_key)

  # Apply stat_fun to each group
  stat_results <- lapply(split_data, function(data_group) {
    stats <- stat_fun(data_group)
    factors_values <- data_group[1, factors, drop = FALSE]
    result <- cbind(factors_values, as.data.frame(t(stats)))
    return(result)
  })

  # Combine results into summary data frame
  summary_df <- do.call(rbind, stat_results)
  rownames(summary_df) <- NULL
  stat_names <- names(summary_df)[!(names(summary_df) %in% factors)]

  # Process post_predict if provided
  if (!is.null(post_predict)) {
    if (!all(c(factors, 'postn') %in% names(post_predict))) {
      stop("post_predict must have columns: ", paste(c(factors, 'postn'), collapse = ', '))
    }
    post_predict$group_key <- create_group_key(post_predict, factors)
    split_post <- split(post_predict, list(post_predict$group_key, post_predict$postn))

    # Apply stat_fun to each group in post_predict
    stat_post_results <- lapply(split_post, function(data_group) {
      stats <- stat_fun(data_group)
      factors_values <- data_group[1, factors, drop = FALSE]
      postn <- data_group$postn[1]
      group_key <- data_group$group_key[1]
      result <- cbind(factors_values, group_key = group_key, postn = postn, as.data.frame(t(stats)))
      return(result)
    })

    # Combine and compute quantiles
    post_stats_df <- do.call(rbind, stat_post_results)
    quantile_results_list <- lapply(stat_names, function(sname) {
      agg_result <- aggregate(as.formula(paste0("`", sname, "` ~ ", paste(c('group_key', factors), collapse = '+'))),
                              data = post_stats_df,
                              FUN = function(x) quantile(x, probs = quants))
      quantile_matrix <- agg_result[[sname]]
      colnames(quantile_matrix) <- paste0(sname, "_pred_", quants)
      result <- cbind(agg_result[, c('group_key', factors), drop = FALSE], quantile_matrix)
      return(result)
    })

    # Merge quantile results into summary_df
    quantile_results_df <- Reduce(function(x, y) merge(x, y, by = c('group_key', factors)),
                                  quantile_results_list)
    summary_df <- merge(summary_df, quantile_results_df, by = factors)
  }

  # Plotting if requested
  if (do_plot) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    if (any(is.na(layout))) {
      par(mfrow = coda_setmfrow(Nchains = 1, Nparms = length(summary_df$group_key), nplots = 1))
    } else {
      par(mfrow = layout)
    }
    lty <- 1:length(stat_names)
    max_y <- 1
    xlim <- range(summary_df[,stat_names])
    if(!is.null(post_predict)){
      predict_args <- add_defaults(predict_args, lty = lty, col = "red")
      dens <- lapply(split(post_stats_df, post_stats_df$group_key), function(x){
        lapply(stat_names, function(y){
          do.call(density, c(list(x[,y]), fix_dots(predict_args, density.default, consider_dots = F)))
        })
      })
      max_y <- max(max_y, max(sapply(unlist(dens, recursive = F), function(x) return(max(x$y)))))
      xlim <- range(xlim, post_stats_df[,stat_names])
    }
    ylim <- c(0, max_y)
    dots <- add_defaults(list(...), lty = lty, ylim = ylim, xlim = xlim, col = "black")

    for (group_key in summary_df$group_key) {
      do.call(plot, c(list(NA), fix_dots_plot(add_defaults(dots, main = paste(group_key, "-", stat_name), ylab = "Density", xlab = " "))))
      do.call(abline, c(list(v = summary_df[summary_df$group_key == group_key, stat_names]), fix_dots_plot(dots)))

      if (!is.null(post_predict)) {
        lines_args <- fix_dots_plot(predict_args)
        mapply(dens[[group_key]], lines_args$lty, FUN = function(x, y) do.call(lines, c(list(x), lty = y,
                                                                                        lines_args[names(lines_args) != "lty"])))
        legend_labels <- c("Predicted", "Observed")
        legend_cols <- c(predict_args$col, dots$col)
        legend("topright", legend = legend_labels, col = legend_cols, lty = 1, bty = "n")
      }
      if(length(stat_names) > 1){
        legend(x = "topleft", legend = stat_names, lty = lty, col = "black", bty = "n")
      }
    }
  }
  return(invisible(summary_df[,colnames(summary_df) != "group_key"]))
}

get_def_dens <- function(x, defective_factor, dots){
  p_defective <- prop.table(table(x[[defective_factor]]))
  out <- mapply(split(x, x[[defective_factor]]), p_defective, FUN = function(y, z){
    dens <- do.call(density, c(list(y$rt), fix_dots(dots, density.default, consider_dots = FALSE)))
    tmp <- dens$y * z
    return(list(tmp))})
  return(out)
}

get_def_cdf <- function(x, defective_factor, dots){
  probs <- seq(0.01, 0.99, by = 0.01)
  p_defective <- prop.table(table(x[[defective_factor]]))
  out <- mapply(split(x, x[[defective_factor]]), p_defective, FUN = function(inp, z){
    rtquants <- quantile(inp$rt, probs = probs, type = 1, na.rm = TRUE)
    yvals <- probs * z
    return(list(cbind(x = rtquants, y = yvals)))})
  return(out)
}

prep_data_plot <- function(input = NULL, post_predict = NULL, plot_prior = FALSE, prior_predict = NULL, factors = NULL, defective_factor = NULL, subject = NULL,
                           n_cores, n_post){

  check_predictives <- function(predictives, defective_factor, subject, factors){
    # Check required columns
    required_cols_post <- c("rt", "subjects", "postn", defective_factor)
    missing_cols_post <- setdiff(required_cols_post, names(predictives))
    if (length(missing_cols_post) > 0) {
      stop("Ensure predictives has columns: ", paste(missing_cols_post, collapse = ", "))
    }

    # Handle subject argument
    if (!is.null(subject)) {
      predictives <- predictives[predictives$subjects %in% subject, ]
      predictives$subjects <- droplevels(predictives$subjects)
    }

    # Remove missing or infinite rt
    predictives <- predictives[is.finite(predictives$rt), ]

    # Create group_keys for predictives
    predictives$group_key <- create_group_key(predictives, factors)
    return(predictives)
  }
  if(is.null(input) & is.null(post_predict) & is.null(prior_predict)){
    stop("At least input, post_predict or prior_predict needs to be defined")
  }
  # Prepare data
  if (inherits(input, "emc")) {
    data <- get_data(input)
    if (is.null(post_predict)) {
      post_predict <- predict(input, n_post = n_post, n_cores = n_cores)
    }
    if(is.null(prior_predict) & plot_prior){
      prior_predict <- predict(get_prior(input), n_post, n_cores)
    }
  } else if(inherits(input, "emc.prior")){
    data <- NULL
    prior_predict <- predict(input, n_post, n_cores)
  } else{
    data <- input
  }

  # Check required columns
  required_cols <- c("rt", "subjects", defective_factor)
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Ensure data has columns: ", paste(missing_cols, collapse = ", "))
  }

  # Handle subject argument
  if (!is.null(subject)) {
    snams <- unique(data$subjects)
    if (is.numeric(subject)) subject <- snams[subject]
    if (!all(subject %in% snams)) stop("Subject not present\n")
    data <- data[data$subjects %in% subject, ]
    data$subjects <- droplevels(data$subjects)
  }

  # Remove missing or infinite rt
  data <- data[is.finite(data$rt), ]

  # Handle factors argument
  if (!is.null(factors) && !all(factors %in% names(data))) {
    stop("factors must name factors in data")
  }

  # Create group_keys
  data$group_key <- create_group_key(data, factors)

  # Prepare post_predict
  if (!is.null(post_predict)) {
    post_predict <- check_predictives(post_predict, defective_factor, subject, factors)
  }
  if(!is.null(prior_predict)){
    prior_predict <- check_predictives(prior_predict, defective_factor, subject, factors)
  }
  return(list(data = data, post_predict = post_predict))
}

#' Plot defective densities
#'
#' Plots panels that contain a set of densities for each level of the specified defective factor in the data.
#' These densities are defective; their areas are relative to the respective
#' proportions of the defective factor levels. Across all levels, the area sums to 1.
#' Optionally, posterior predictive densities can be overlaid.
#'
#' @inheritParams plot_cdf
#' @examples
#' # Plot defective densities for each subject and the factor combination in the design:
#' plot_density(forstmann)
#' # or for one subject:
#' plot_density(forstmann, subject = 1)
#' # Now collapsing across subjects and using a different defective factor:
#' plot_density(forstmann, factors = "S", defective_factor = "E")
#' # Or plot posterior predictives
#' plot_density(samples_LNR, n_post = 10)
#' @export
plot_density <- function(input, subject = NULL, factors = NULL,
                             defective_factor = "R",
                             layout = NA, post_predict = NULL, quants = c(0.025, 0.975),
                             predict_args = list(), n_cores = 1, n_post = 100,
                             ...) {
  dots <- add_defaults(list(...), col = "black")
  check <- prep_data_plot(input, post_predict, factors, defective_factor, subject, n_cores = n_cores, n_post = n_post)
  data <- check$data
  post_predict <- check$post_predict
  # Get defective_levels
  defective_levels <- levels(factor(data[[defective_factor]]))
  line_types <- seq_along(defective_levels)

  # Initialize variables
  y_max <- 0
  quantiles <- sort(c(quants, 0.5))
  # Compute xlim based on quantiles
  x_lim_probs <- c(0.001, 0.995)
  dat_quants <- aggregate(rt ~ group_key, data, quantile, probs = x_lim_probs)
  xlim_data <- range(unlist(dat_quants$rt))
  if (!is.null(post_predict)) {
    pred_quants <- aggregate(rt ~ group_key, post_predict, quantile, probs = x_lim_probs)
    xlim_pred <- range(unlist(pred_quants$rt))
    xlim <- range(xlim_data, xlim_pred)
  } else {
    xlim <- xlim_data
  }
  x_d_range <- range(c(post_predict$rt, data$rt))
  x_grid <- seq(x_d_range[1], x_d_range[2], length.out = 512) # 512 is the default number of density evaluations
  dots <- add_defaults(dots, from = x_d_range[1], to = x_d_range[2])
  predict_args <- add_defaults(predict_args, from = x_d_range[1], to = x_d_range[2], col = "red")
  # Compute densities and y_max
  dens_data <- lapply(split(data, data$group_key), get_def_dens, defective_factor, dots)
  y_max <- max(sapply(unlist(dens_data, recursive = FALSE), max))
  # Process post_predict densities
  if (!is.null(post_predict)) {
    dens_pred <- lapply(split(post_predict, post_predict$group_key), FUN = function(x){
      lapply(split(x, x$postn), get_def_dens, defective_factor, predict_args)
    })
    quants <- lapply(dens_pred, function(x){
      out <- list()
      for(deflev in defective_levels){
        out[[deflev]] <- apply(do.call(cbind, lapply(x, function(y) return(y[[deflev]]))), 1, quantile, probs = quantiles)
      }
      return(out)
    })
    y_max <- max(sapply(unlist(quants, FALSE), max), y_max)
  }

  # Set global y-limits
  ylim_global <- c(0, y_max)

  # Set up plotting layout
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  unique_group_keys <- unique(data$group_key)
  if (any(is.na(layout))) {
    par(mfrow = coda_setmfrow(Nchains = 1, Nparms = length(unique_group_keys), nplots = 1))
  } else {
    par(mfrow = layout)
  }

  # Plotting loop
  for (group_key in unique_group_keys) {
    # Set plot arguments with consistent y-limits
    plot_args <- add_defaults(dots, xlim = xlim, ylim = ylim_global, main = group_key, xlab = "RT", ylab = "Defective Density")
    plot_args <- fix_dots_plot(plot_args)
    do.call(plot, c(list(NA), plot_args))
    cur_data <- dens_data[[group_key]]
    # Plot observed densities
    for (i in seq_along(defective_levels)) {
      level <- defective_levels[i]
      if (!is.null(cur_data[[level]])) {
        lines_args <- add_defaults(dots, lty = line_types[i])
        lines_args <- fix_dots_plot(lines_args)
        do.call(lines, c(list(x = x_grid, y = cur_data[[level]]), lines_args))
      }
    }

    # Plot predicted densities
    if (!is.null(post_predict)) {
      cur_quants <- quants[[group_key]]
      for (i in seq_along(defective_levels)) {
        level <- defective_levels[i]
        if (!is.null(cur_quants[[level]])) {
          lines_args <- add_defaults(predict_args, lty = line_types[i])
          lines_args <- fix_dots_plot(lines_args)
          do.call(lines, c(list(x = x_grid, y = cur_quants[[level]][2,]), lines_args))
          # Add ribbon
          adj_color <- do.call(adjustcolor, fix_dots(add_defaults(predict_args, alpha.f = .2), adjustcolor))
          polygon_args <- predict_args
          polygon_args$col <- adj_color
          do.call(polygon, c(list(x = c(x_grid, rev(x_grid)), border = NA,
                                  y = c(cur_quants[[level]][1,], rev(cur_quants[[level]][3,]))), fix_dots_plot(polygon_args)))
        }
      }
    }

    # Add legends
    legend(x = "topright", legend = defective_levels, lty = line_types, col = "black", title = defective_factor, bty = "n")
    if (!is.null(post_predict)) {
      legend(x = "top", legend = c("Real", "Pred"), lty = 1, col = c(dots$col, predict_args$col), title = "Source", bty = "n")
    }
  }
}
