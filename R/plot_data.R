create_group_key <- function(df, factors) {
  if (length(factors) == 0) return(rep("All Data", nrow(df)))
  apply(df[, factors, drop = FALSE], 1, function(x) paste(paste(factors, x, sep = "="), collapse = " "))
}


prep_data_plot <- function(input = NULL, post_predict = NULL,
                           prior_predict = NULL, to_plot, limits,
                           factors = NULL, defective_factor = NULL, subject = NULL,
                           n_cores, n_post, make_xlim = TRUE){
  check_data_plot <- function(data, defective_factor, subject, factors){
    # Check required columns
    required_cols_post <- c("rt", "subjects", defective_factor)
    missing_cols_post <- setdiff(required_cols_post, names(data))
    if (length(missing_cols_post) > 0) {
      stop("Ensure data has columns: ", paste(missing_cols_post, collapse = ", "))
    }
    # Handle factors argument
    if (!is.null(factors) && !all(factors %in% names(data))) {
      stop("factors must name factors in data")
    }
    # Handle subject argument
    if (!is.null(subject)) {
      if(is.numeric(subject)){
        data <- data[data$subjects %in% unique(data$subjects)[subject], ]
      } else{
        data <- data[data$subjects %in% subject, ]
      }
      data$subjects <- droplevels(data$subjects)
    }

    # Remove missing or infinite rt
    data <- data[is.finite(data$rt), ]

    # Create group_keys for data
    data$group_key <- create_group_key(data, factors)
    return(data)
  }
  if(is.null(input) & is.null(post_predict) & is.null(prior_predict)){
    stop("At least input, post_predict or prior_predict needs to be defined")
  }
  # Prepare data
  if (inherits(input, "emc")) {
    # samples object, first get data
    data <- get_data(input)
    # Generate posterior predictives
    if (is.null(post_predict) & ('posterior' %in% to_plot)) {
      post_predict <- predict(input, n_post = n_post, n_cores = n_cores)
    }
    # See if we also want to generate prior predictives
    if(is.null(prior_predict) & ('prior' %in% to_plot)){
      prior_predict <- predict(get_prior(input), data = data, n_post = n_post, n_cores = n_cores)
    }
  } else{
    if('postn' %in% colnames(input)){
      if(all(to_plot == 'posterior')){
        post_predict <- input
      } else if(all(to_plot == 'prior')){
        stop("for prior predictives we still need the data to replicate the design factors
             please provide the data as the input argument and the prior predictives as the prior_predict argument")
      }
    } else{
      data <- input
    }
  }
  # If prior or posterior predictives were provided probably makes sense to plot them ;)
  to_plot <- c('real', 'posterior', 'prior')[!c(is.null(data), is.null(post_predict), is.null(prior_predict))]
  to_plot <- unique(c(to_plot))
  xlim <- NULL
  # Compute xlim based on quantiles and perform checks
  x_lim_probs <- c(0.0001, 0.995)
  if(!is.null(data)){
    data <- check_data_plot(data, defective_factor, subject, factors)
    if('real' %in% limits & make_xlim){
      dat_quants <- aggregate(rt ~ group_key, data, quantile, probs = x_lim_probs)
      xlim <- range(xlim, unlist(dat_quants$rt))
    }
  }
  # Prepare post_predict
  if (!is.null(post_predict)) {
    post_predict <- check_data_plot(post_predict, defective_factor, subject, factors)
    if('posterior' %in% limits & make_xlim){
      pred_quants <- aggregate(rt ~ group_key, post_predict, quantile, probs = x_lim_probs)
      xlim <- range(xlim, range(unlist(pred_quants$rt)))
    }
  }
  if(!is.null(prior_predict)){
    prior_predict <- check_data_plot(prior_predict, defective_factor, subject, factors)
    if('prior' %in% limits & make_xlim){
      xlim_probs_prior <- c(0.001, 0.95) # Slightly tighter plotting range for prior predictives
      prior_quants <- aggregate(rt ~ group_key, prior_predict, quantile, probs = xlim_probs_prior)
      xlim <- range(xlim, range(unlist(prior_quants$rt)))
    }
  }
  return(list(data = data, post_predict = post_predict, prior_predict = prior_predict, xlim = xlim, to_plot = to_plot))
}


#' Plot defective cumulative distribution functions
#'
#' Plots panels that contain a set of cumulative distribution functions (CDFs) for each level of the specified defective factor in the data.
#' These CDFs are defective; their maximum values are relative to the respective
#' proportions of the defective factor levels. Across all levels, the maximum value sums to 1.
#' Optionally, posterior predictive and/or prior predictive CDFs can be overlaid.
#'
#' @param input Either an emc object or a data frame. The data frame must be in EMC2 format with at least subjects, rt, and the defective factor columns.
#' @param post_predict An optional data frame containing posterior predictive datasets.
#' Should have the same columns as data, with an additional postn column indicating the index of the predicted dataset.
#' @param prior_predict An optional data frame containing prior predictive datasets.
#' Should have the same columns as data, with an additional postn column indicating the index of the predicted dataset.
#' @param subject An integer or character string selecting a subject from the data.
#' If specified, only that subject is plotted. By default (NULL), all subjects are plotted.
#' @param quants A numeric vector specifying the quantiles for the credible intervals.
#' @param factors A character vector of the factor names in the design to aggregate across.
#' Defaults to all factors (i.e., NULL).
#' @param defective_factor A character string specifying the factor to compute defective CDFs for. Defaults to "R".
#' @param n_cores An integer specifying the number of cores with which to calculate predictives if specified.
#' @param n_post An optional integer specifying how many posterior predictives to generate if the input is an emc object.
#' @param layout A vector indicating which layout to use as in par(mfrow = layout). If NA, will automatically generate an appropriate layout.
#' @param to_plot A character vector indicating whether to plot the real data, the posterior predictives, prior_predictives or a combination
#' @param use_lim A character vector indicating whether to base the plot limits on the real data, the posterior predictives, prior_predictives or a combination
#' @param posterior_args Optional arguments that can be passed to plotting functions for the posterior predictives.
#' @param prior_args  Optional arguments that can be passed to plotting functions for the prior predictives.
#' @param ... Optional arguments that can be passed to plotting functions for the real data.
#'
#' @return Invisibly returns a data frame containing RT quantiles and mean proportions per factor combination and defective factor level.
#' If post_predict is provided or the input is an emc object, predicted values are included.
#' @examples
#' # Plot defective cdf for each subject and the factor combination in the design:
#' plot_cdf(forstmann)
#' # or for one subject:
#' plot_cdf(forstmann, subject = 1)
#' # Now collapsing across subjects and using a different defective factor:
#' plot_cdf(forstmann, factors = "S", defective_factor = "E")
#' # Or plot posterior/prior predictives with an emc object
#' # potentially also add 'prior' to to_plot
#' plot_cdf(samples_LNR, to_plot = c('real', 'posterior'), n_post = 10)
#' @export
plot_cdf <- function(input, post_predict = NULL, prior_predict = NULL, subject = NULL, quants = c(0.025, 0.975),
                         factors = NULL, defective_factor = "R", n_cores = 1, n_post = 100,
                         layout = NA, to_plot = c('real', 'posterior', 'prior')[1:2], use_lim = c('real', 'posterior', 'prior')[1:2],
                         posterior_args = list(), prior_args = list(), ...) {
  dots <- add_defaults(list(...), col = "black")
  check <- prep_data_plot(input, post_predict, prior_predict, to_plot, use_lim,
                          factors, defective_factor, subject, n_cores, n_post)
  data <- check$data
  post_predict <- check$post_predict
  prior_predict <- check$prior_predict
  to_plot <- check$to_plot
  xlim <- check$xlim
  # Get defective_levels
  defective_levels <- levels(factor(data[[defective_factor]]))
  line_types <- seq_along(defective_levels)

  # Initialize variables
  y_max <- 0
  quantiles <- sort(c(quants, 0.5))
  # Compute xlim based on quantiles
  posterior_args <- add_defaults(posterior_args, col = "red")
  prior_args <- add_defaults(prior_args, col = "darkgreen")
  # Compute densities and y_max
  if(!is.null(data) & c('real' %in% to_plot)){
    dens_data <- lapply(split(data, data$group_key), get_def_cdf, defective_factor, dots)
    if('real' %in% use_lim){
      y_max <- max(sapply(unlist(dens_data, recursive = FALSE), function(x) return(max(x[,'y']))))
    }
  }
  # Process post_predict densities
  if (!is.null(post_predict) & c('posterior' %in% to_plot)) {
    dens_pred <- lapply(split(post_predict, post_predict$group_key), FUN = function(x){
      lapply(split(x, x$postn), get_def_cdf, defective_factor, posterior_args)
    })
    quants <- lapply(dens_pred, function(x){
      out <- list()
      for(deflev in defective_levels){
        tmp <- apply(do.call(cbind, lapply(x, function(q) return(q[[deflev]][,'x']))), 1, quantile, probs = quantiles)
        out[[deflev]] <- rbind(tmp, apply(do.call(cbind, lapply(x, function(q) return(q[[deflev]][,'y']))), 1, median))
      }
      return(out)
    })
    if('posterior' %in% use_lim){
      y_max <- max(y_max, sapply(unlist(quants, FALSE), function(x) return(max(x[4,]))))
    }
  }
  if (!is.null(prior_predict) & c('prior' %in% to_plot)) {
    dens_prior <- lapply(split(prior_predict, prior_predict$group_key), FUN = function(x){
      lapply(split(x, x$postn), get_def_cdf, defective_factor, prior_args)
    })
    quants_prior <- lapply(dens_prior, function(x){
      out <- list()
      for(deflev in defective_levels){
        tmp <- apply(do.call(cbind, lapply(x, function(q) return(q[[deflev]][,'x']))), 1, quantile, probs = quantiles)
        out[[deflev]] <- rbind(tmp, apply(do.call(cbind, lapply(x, function(q) return(q[[deflev]][,'y']))), 1, median))
      }
      return(out)
    })
    if('prior' %in% use_lim){
      y_max <- max(y_max, sapply(unlist(quants_prior, FALSE), function(x) return(max(x[4,]))))
    }
  }

  # Set global y-limits, add a little bit of leeway
  ylim_global <- c(0, y_max*1.1)

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
    if (!is.null(data) & c('real' %in% to_plot)) {
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
    }

    # Plot predicted densities
    if (!is.null(post_predict) & c('posterior' %in% to_plot)) {
      cur_quants <- quants[[group_key]]
      for (i in seq_along(defective_levels)) {
        level <- defective_levels[i]
        if (!is.null(cur_quants[[level]])) {
          lines_args <- add_defaults(posterior_args, lty = line_types[i])
          lines_args <- fix_dots_plot(lines_args)
          do.call(lines, c(list(x = cur_quants[[level]][2,], y = cur_quants[[level]][4,]), lines_args))
          # Add ribbon
          adj_color <- do.call(adjustcolor, fix_dots(add_defaults(posterior_args, alpha.f = .2), adjustcolor))
          polygon_args <- posterior_args
          polygon_args$col <- adj_color
          do.call(polygon, c(list(y = c(cur_quants[[level]][4,], rev(cur_quants[[level]][4,])), border = NA,
                                  x = c(cur_quants[[level]][1,], rev(cur_quants[[level]][3,]))), fix_dots_plot(polygon_args)))
        }
      }
    }
    if (!is.null(prior_predict) & c('prior' %in% to_plot)) {
      cur_quants <- quants_prior[[group_key]]
      for (i in seq_along(defective_levels)) {
        level <- defective_levels[i]
        if (!is.null(cur_quants[[level]])) {
          lines_args <- add_defaults(prior_args, lty = line_types[i])
          lines_args <- fix_dots_plot(lines_args)
          do.call(lines, c(list(x = cur_quants[[level]][2,], y = cur_quants[[level]][4,]), lines_args))
          # Add ribbon
          adj_color <- do.call(adjustcolor, fix_dots(add_defaults(prior_args, alpha.f = .2), adjustcolor))
          polygon_args <- prior_args
          polygon_args$col <- adj_color
          do.call(polygon, c(list(y = c(cur_quants[[level]][4,], rev(cur_quants[[level]][4,])), border = NA,
                                  x = c(cur_quants[[level]][1,], rev(cur_quants[[level]][3,]))), fix_dots_plot(polygon_args)))
        }
      }
    }
    # Add legends
    legend(x = "topright", legend = defective_levels, lty = line_types, col = "black", title = defective_factor, bty = "n")
    if (length(to_plot) > 1) {
      cols <- c(dots$col, posterior_args$col, prior_args$col)[c('real', 'posterior', 'prior') %in% to_plot]
      legend(x = "top", legend = to_plot, lty = 1, col = cols, title = "Source", bty = "n")
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
#' @return an invisible data frame with the stat applied to the real data, posterior predictives and/or prior predictives
#' @export
#'
#' @examples
#' # For example plot the observed and predicted response accuracy
# correct_fun <- function(data) mean(data$S == data$R)
# plot_stat(samples_LNR, stat_fun = correct_fun, n_post = 10)
#' # Can also apply more sophisticated statistics
#' drt <- function(data) diff(tapply(data$rt,data[,c("E")],mean))
#' plot_stat(samples_LNR, stat_fun = drt, n_post = 10, stat_name = "RT diff Speed - A/N")
#'
plot_stat <- function(input, post_predict = NULL, prior_predict = NULL, stat_fun, stat_name = NULL,
                      subject = NULL, factors = NULL, n_cores = 1, n_post = 100,
                      quants = c(0.025, 0.5, 0.975),
                      layout = NA, to_plot = c('real', 'posterior', 'prior')[1:2], use_lim = c('real', 'posterior', 'prior')[1:2],
                      posterior_args = list(), prior_args = list(), ...) {
  check <- prep_data_plot(input, post_predict, prior_predict, to_plot, use_lim,
                          factors, defective_factor = NULL, subject, n_cores, n_post,
                          make_xlim = FALSE)
  data <- check$data
  post_predict <- check$post_predict
  prior_predict <- check$prior_predict
  to_plot <- check$to_plot
  summary_df <- NULL
  if(!is.null(data)){
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
  }

  # Process post_predict if provided
  if (!is.null(post_predict)) {
    split_post <- split(post_predict, list(post_predict$group_key, post_predict$postn))

    # Apply stat_fun to each group in post_predict
    stat_post <- lapply(split_post, function(data_group) {
      stats <- stat_fun(data_group)
      factors_values <- data_group[1, factors, drop = FALSE]
      postn <- data_group$postn[1]
      group_key <- data_group$group_key[1]
      result <- cbind(factors_values, group_key = group_key, postn = postn, as.data.frame(t(stats)))
      return(result)
    })

    # Combine and compute quantiles
    post_stats_df <- do.call(rbind, stat_post)
    quantile_post <- lapply(stat_names, function(sname) {
      agg_result <- aggregate(as.formula(paste0("`", sname, "` ~ ", paste(c('group_key', factors), collapse = '+'))),
                              data = post_stats_df,
                              FUN = function(x) quantile(x, probs = quants))
      quantile_matrix <- agg_result[[sname]]
      colnames(quantile_matrix) <- paste0(sname, "_posterior_", quants)
      result <- cbind(agg_result[, c('group_key', factors), drop = FALSE], quantile_matrix)
      return(result)
    })

    # Merge quantile results into summary_df
    quantile_post_df <- Reduce(function(x, y) merge(x, y, by = c('group_key', factors)),
                               quantile_post)
    summary_df <- merge(summary_df, quantile_post_df)
  }

  if (!is.null(prior_predict)) {
    split_prior <- split(prior_predict, list(prior_predict$group_key, prior_predict$postn))

    # Apply stat_fun to each group in prior_predict
    stat_prior <- lapply(split_prior, function(data_group) {
      stats <- stat_fun(data_group)
      factors_values <- data_group[1, factors, drop = FALSE]
      postn <- data_group$postn[1]
      group_key <- data_group$group_key[1]
      result <- cbind(factors_values, group_key = group_key, postn = postn, as.data.frame(t(stats)))
      return(result)
    })

    # Combine and compute quantiles
    prior_stats_df <- do.call(rbind, stat_prior)
    quantile_prior <- lapply(stat_names, function(sname) {
      agg_result <- aggregate(as.formula(paste0("`", sname, "` ~ ", paste(c('group_key', factors), collapse = '+'))),
                              data = prior_stats_df,
                              FUN = function(x) quantile(x, probs = quants))
      quantile_matrix <- agg_result[[sname]]
      colnames(quantile_matrix) <- paste0(sname, "_prior_", quants)
      result <- cbind(agg_result[, c('group_key', factors), drop = FALSE], quantile_matrix)
      return(result)
    })

    # Merge quantile results into summary_df
    quantile_prior_df <- Reduce(function(x, y) merge(x, y, by = c('group_key', factors)),
                                  quantile_prior)
    summary_df <- merge(summary_df, quantile_prior_df)
  }


  # Plotting loop
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
  # Check limits
  if(!is.null(post_predict) & c('posterior' %in% to_plot)){
    posterior_args <- add_defaults(posterior_args, lty = lty, col = "red")
    dens <- lapply(split(post_stats_df, post_stats_df$group_key), function(x){
      lapply(stat_names, function(y){
        do.call(density, c(list(x[,y]), fix_dots(posterior_args, density.default, consider_dots = F)))
      })
    })
    if('posterior' %in% use_lim){
      max_y <- max(max_y, max(sapply(unlist(dens, recursive = F), function(x) return(max(x$y)))))
      xlim <- range(xlim, post_stats_df[,stat_names])
    }

  }

  if(!is.null(prior_predict) & c('prior' %in% to_plot)){
    prior_args <- add_defaults(prior_args, lty = lty, col = "darkgreen")
    dens_prior <- lapply(split(prior_stats_df, prior_stats_df$group_key), function(x){
      lapply(stat_names, function(y){
        do.call(density, c(list(x[,y]), fix_dots(prior_args, density.default, consider_dots = F)))
      })
    })
    if('prior' %in% use_lim){
      max_y <- max(max_y, max(sapply(unlist(dens_prior, recursive = F), function(x) return(max(x$y)))))
      xlim <- range(xlim, quantile(prior_stats_df[,stat_names], probs = c(0.025, 0.975)))
    }

  }

  ylim <- c(0, max_y)
  dots <- add_defaults(list(...), lty = lty, ylim = ylim, xlim = xlim, col = "black")

  for (group_key in summary_df$group_key) {
    do.call(plot, c(list(NA), fix_dots_plot(add_defaults(dots, main = paste(group_key, "-", stat_name), ylab = "Density", xlab = " "))))
    if('real' %in% to_plot){
      do.call(abline, c(list(v = summary_df[summary_df$group_key == group_key, stat_names]), fix_dots_plot(dots)))
    }
    if (!is.null(post_predict)) {
      lines_args <- fix_dots_plot(posterior_args)
      mapply(dens[[group_key]], lines_args$lty, FUN = function(x, y) do.call(lines, c(list(x), lty = y,
                                                                                      lines_args[names(lines_args) != "lty"])))
    }
    if (!is.null(prior_predict)) {
      lines_args <- fix_dots_plot(prior_args)
      mapply(dens_prior[[group_key]], lines_args$lty, FUN = function(x, y) do.call(lines, c(list(x), lty = y,
                                                                                      lines_args[names(lines_args) != "lty"])))
    }
    if(length(stat_names) > 1){
      legend(x = "topleft", legend = stat_names, lty = lty, col = "black", bty = "n")
    }
    if (length(to_plot) > 1) {
      cols <- c(dots$col, posterior_args$col, prior_args$col)[c('real', 'posterior', 'prior') %in% to_plot]
      legend(x = "topright", legend = to_plot, lty = 1, col = cols, title = "Source", bty = "n")
    }
  }
  return(invisible(summary_df[,!grepl("group_key", colnames(summary_df))]))
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
#' Plot defective densities
#'
#' Plots panels that contain a set of densities for each level of the specified defective factor in the data.
#' These densities are defective; their areas are relative to the respective
#' proportions of the defective factor levels. Across all levels, the area sums to 1.
#' Optionally, posterior/prior predictive densities can be overlaid.
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
plot_density <- function(input, post_predict = NULL, prior_predict = NULL, subject = NULL, quants = c(0.025, 0.975),
                         factors = NULL, defective_factor = "R", n_cores = 1, n_post = 100,
                         layout = NA, to_plot = c('real', 'posterior', 'prior')[1:2], use_lim = c('real', 'posterior', 'prior')[1:2],
                         posterior_args = list(), prior_args = list(), ...) {
  check <- prep_data_plot(input, post_predict, prior_predict, to_plot, use_lim,
                          factors, defective_factor, subject, n_cores, n_post)
  data <- check$data
  post_predict <- check$post_predict
  prior_predict <- check$prior_predict
  to_plot <- check$to_plot
  xlim <- check$xlim
  # Get defective_levels
  defective_levels <- levels(factor(data[[defective_factor]]))
  line_types <- seq_along(defective_levels)

  # Initialize variables
  y_max <- 0
  quantiles <- sort(c(quants, 0.5))
  dots <- add_defaults(list(...), col = "black")
  posterior_args <- add_defaults(posterior_args, col = "red")
  prior_args <- add_defaults(prior_args, col = "darkgreen")

  # Compute densities and y_max
  if(!is.null(data) & c('real' %in% to_plot)){
    x_d_range_real <- range(data$rt)
    x_grid_real <- seq(x_d_range_real[1], x_d_range_real[2], length.out = 512) # 512 is the default number of density evaluations
    dots <- add_defaults(dots, from = x_d_range_real[1], to = x_d_range_real[2])
    dens_data <- lapply(split(data, data$group_key), get_def_dens, defective_factor, dots)
    if('real' %in% use_lim){
      y_max <- max(sapply(unlist(dens_data, recursive = FALSE), max))
    }
  }

  # Process post_predict densities
  if (!is.null(post_predict) & c('posterior' %in% to_plot)) {
    x_d_range_predict <- range(post_predict$rt)
    x_grid_predict <- seq(x_d_range_predict[1], x_d_range_predict[2], length.out = 512)
    posterior_args <- add_defaults(posterior_args, from = x_d_range_predict[1], to = x_d_range_predict[2])
    dens_pred <- lapply(split(post_predict, post_predict$group_key), FUN = function(x){
      lapply(split(x, x$postn), get_def_dens, defective_factor, posterior_args)
    })
    quants <- lapply(dens_pred, function(x){
      out <- list()
      for(deflev in defective_levels){
        out[[deflev]] <- apply(do.call(cbind, lapply(x, function(y) return(y[[deflev]]))), 1, quantile, probs = quantiles)
      }
      return(out)
    })
    if('posterior' %in% use_lim){
      y_max <- max(sapply(unlist(quants, FALSE), max), y_max)
    }
  }
  if (!is.null(prior_predict) & c('prior' %in% to_plot)) {
    x_d_range_prior <- quantile(prior_predict$rt, probs = c(0.001, 0.975))
    x_grid_prior <- seq(x_d_range_prior[1], x_d_range_prior[2], length.out = 512)
    prior_args <- add_defaults(prior_args, from = x_d_range_prior[1], to = x_d_range_prior[2])
    dens_prior <- lapply(split(prior_predict, prior_predict$group_key), FUN = function(x){
      lapply(split(x, x$postn), get_def_dens, defective_factor, posterior_args)
    })
    quants_prior <- lapply(dens_prior, function(x){
      out <- list()
      for(deflev in defective_levels){
        out[[deflev]] <- apply(do.call(cbind, lapply(x, function(y) return(y[[deflev]]))), 1, quantile, probs = quantiles)
      }
      return(out)
    })
    if('prior' %in% use_lim){
      y_max <- max(sapply(unlist(quants_prior, FALSE), max), y_max)
    }
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
    # Plot observed densities
    if (!is.null(data) & c('real' %in% to_plot)) {
      cur_data <- dens_data[[group_key]]
      for (i in seq_along(defective_levels)) {
        level <- defective_levels[i]
        if (!is.null(cur_data[[level]])) {
          lines_args <- add_defaults(dots, lty = line_types[i])
          lines_args <- fix_dots_plot(lines_args)
          do.call(lines, c(list(x = x_grid_real, y = cur_data[[level]]), lines_args))
        }
      }
    }
    # Plot predicted densities
    if (!is.null(post_predict) & c('posterior' %in% to_plot)) {
      cur_quants <- quants[[group_key]]
      for (i in seq_along(defective_levels)) {
        level <- defective_levels[i]
        if (!is.null(cur_quants[[level]])) {
          lines_args <- add_defaults(posterior_args, lty = line_types[i])
          lines_args <- fix_dots_plot(lines_args)
          do.call(lines, c(list(x = x_grid_predict, y = cur_quants[[level]][2,]), lines_args))
          # Add ribbon
          adj_color <- do.call(adjustcolor, fix_dots(add_defaults(posterior_args, alpha.f = .2), adjustcolor))
          polygon_args <- posterior_args
          polygon_args$col <- adj_color
          do.call(polygon, c(list(x = c(x_grid_predict, rev(x_grid_predict)), border = NA,
                                  y = c(cur_quants[[level]][1,], rev(cur_quants[[level]][3,]))), fix_dots_plot(polygon_args)))
        }
      }
    }
    # Plot prior densities
    if (!is.null(prior_predict) & c('prior' %in% to_plot)) {
      cur_quants <- quants_prior[[group_key]]
      for (i in seq_along(defective_levels)) {
        level <- defective_levels[i]
        if (!is.null(cur_quants[[level]])) {
          lines_args <- add_defaults(prior_args, lty = line_types[i])
          lines_args <- fix_dots_plot(lines_args)
          do.call(lines, c(list(x = x_grid_prior, y = cur_quants[[level]][2,]), lines_args))
          # Add ribbon
          adj_color <- do.call(adjustcolor, fix_dots(add_defaults(prior_args, alpha.f = .2), adjustcolor))
          polygon_args <- prior_args
          polygon_args$col <- adj_color
          do.call(polygon, c(list(x = c(x_grid_prior, rev(x_grid_prior)), border = NA,
                                  y = c(cur_quants[[level]][1,], rev(cur_quants[[level]][3,]))), fix_dots_plot(polygon_args)))
        }
      }
    }

    # Add legends
    legend(x = "topright", legend = defective_levels, lty = line_types, col = "black", title = defective_factor, bty = "n")
    if (length(to_plot) > 1) {
      cols <- c(dots$col, posterior_args$col, prior_args$col)[c('real', 'posterior', 'prior') %in% to_plot]
      legend(x = "top", legend = to_plot, lty = 1, col = cols, title = "Source", bty = "n")
    }
  }
}

