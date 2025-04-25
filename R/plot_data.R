create_group_key <- function(df, factors) {
  if (length(factors) == 0) return(rep("All Data", nrow(df)))
  apply(df[, factors, drop = FALSE], 1, function(x) paste(paste(factors, x, sep = "="), collapse = " "))
}


check_data_plot <- function(data, defective_factor, subject, factors) {

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
  n_bins <- 4
  for(fact in factors){
    if(is.numeric(data[,fact])){
      if(length(unique(data[,fact])) > 6){
        quartile_breaks <- quantile(data[,fact], probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
        # Bin the data into quartiles using these breakpoints
        data[,fact] <- cut(data[,fact], breaks = quartile_breaks, include.lowest = TRUE, labels = paste0("Q", 1:n_bins))
      }
    }
  }
  # Handle subject argument
  if (!is.null(subject)) {
    if(is.numeric(subject)) {
      data <- data[data$subjects %in% unique(data$subjects)[subject], ]
    } else {
      data <- data[data$subjects %in% subject, ]
    }
    data$subjects <- factor(data$subjects)
  }

  # Remove missing or infinite rt
  data <- data[is.finite(data$rt), ]

  # --- Faster group_key creation when postn is present ---
  grp_cols <- unique(c("subjects", defective_factor, factors))
  grp_cols <- intersect(grp_cols, names(data))  # Just to be safe

  # 2. Extract the rows that are unique with respect to these columns (excluding 'postn')
  tmp_unique <- data[!duplicated(data[, grp_cols, drop = FALSE]), grp_cols, drop = FALSE]

  # 3. Create one group_key per unique row
  tmp_unique$group_key <- create_group_key(tmp_unique, factors)

  # 4. Merge back group_key into the full data (by the grouping columns)
  data <- merge(
    x = data,
    y = tmp_unique,
    by = grp_cols,
    all.x = TRUE,
    sort = FALSE
  )
  return(data)
}

get_emc_functions <- function(emc){
  out <- list()
  design <- get_design(emc)[[1]]
  out$lM <- design$matchfun
  out[names(design$Ffunctions)] <- design$Ffunctions
  return(out)
}

calc_functions <- function(functions, input){
  if(is.null(functions) || length(functions) == 0) return(input)
  fnams <- names(functions)
  input[['lR']] <- input$R
  for(i in 1:length(functions)){
    input[[fnams[i]]] <- functions[[i]](input)
  }
  return(input)
}

prep_data_plot <- function(input, post_predict, prior_predict, to_plot, limits,
                           factors = NULL, defective_factor = NULL, subject = NULL,
                           n_cores, n_post, functions){
  if(!is.data.frame(input) && !inherits(input, "emc") && !is.null(post_predict) && length(input) != length(post_predict)){
    stop("If input is a list, post_predict must be a list of the same length")
  }
  if(!is.data.frame(input) && !inherits(input, "emc") && !is.null(prior_predict) && length(input) != length(prior_predict)){
    stop("If input is a list, prior_predict must be a list of the same length")
  }
  datasets <- list()
  sources <- c()
  # Check for regular input
  if(!is.data.frame(input) && !inherits(input, "emc")){
    if(is.null(names(input))) stop("If input is a list, it must have names")
    is_emc <- lapply(input, function(x) inherits(x, "emc"))
    if(!all(is_emc) || all(is_emc)){
      stop("Input must be either all emc objects or no emc objects")
    }
  } else{
    input <- list(data = input)
  }
  # Check for post_predict and prior_predict
  if(!is.data.frame(post_predict) && is.list(post_predict)){
    if(is.null(names(post_predict))) stop("If post_predict is a list, it must have names")
    datasets[names(post_predict)] <- post_predict
    sources[names(post_predict)] <- 'posterior'
  } else if(!is.null(post_predict)){
    datasets[['posterior']] <- post_predict
    sources['posterior'] <- 'posterior'
  } else{
    post_predict <- vector("list", length(input))
  }
  if(!is.data.frame(prior_predict)  && is.list(prior_predict)){
    if(is.null(names(prior_predict))) stop("If prior_predict is a list, it must have names")
    datasets[names(prior_predict)] <- prior_predict
    sources[names(prior_predict)] <- 'prior'
  } else if(!is.null(prior_predict)){
    datasets[['prior']] <- prior_predict
    sources['prior'] <- 'prior'
  }  else{
    prior_predict <- vector("list", length(input))
  }

  all_data <- list()
  # Check if regular input is data or emc
  for(k in 1:length(input)){
    # Prepare data
    if (inherits(input[[k]], "emc")) {
      all_data[[names(input)[k]]] <- get_data(input[[k]])
      functions <- c(get_emc_functions(input[[k]]), functions)
    } else{
      all_data[names(input)[k]] <- input[k]
    }
  }
  if(length(unique(all_data)) == 1){
    all_data <- all_data[1]
    datasets['data'] <- all_data
    sources['data'] <- 'data'
  } else{
    datasets[names(input)] <- all_data
    sources[names(input)] <- 'data'
  }
  # Check if posterior or prior predictives need to be generated
  for(k in 1:length(input)){
    # Generate posterior predictives
    if (is.null(post_predict[[k]]) & ('posterior' %in% to_plot) & inherits(input[[k]], "emc")) {
      # In this case we provided multiple emc datasets but they have the same data
      # Just give them the name of the emc data if there is no prior to be plotted
      # as well.
      if(length(all_data) > 1 && !('prior' %in% to_plot)){
        datasets[[paste0(names(input)[k], ' - posterior')]] <- predict(input[[k]], n_post = n_post, n_cores = n_cores)
        sources[paste0(names(input)[k], ' - posterior')] <- 'posterior'
      } else{
        datasets[['posterior']] <- predict(input[[k]], n_post = n_post, n_cores = n_cores)
        sources['posterior'] <- 'posterior'
      }
    }
    # See if we also want to generate prior predictives
    if (is.null(prior_predict[[k]]) & ('prior' %in% to_plot) & inherits(input[[k]], "emc")) {
      # Same logic as for posterior
      if(length(all_data) > 1 && !('posterior' %in% to_plot)){
        datasets[[paste0(names(input)[k], ' - prior')]] <- predict(input[[k]], n_post = n_post, n_cores = n_cores)
        sources[paste0(names(input)[k], ' - prior')] <- 'prior'
      } else{
        datasets[['prior']] <- predict(get_prior(input[[k]]), data = get_data(input[[k]]), n_post = n_post, n_cores = n_cores)
        sources['prior'] <- 'prior'
      }
    }
  }
  xlim <- NULL
  # Compute xlim based on quantiles and perform checks
  for(j in 1:length(datasets)){
    datasets[[j]] <- calc_functions(functions, datasets[[j]])
    datasets[[j]] <- check_data_plot(datasets[[j]], defective_factor, subject, factors)
    if(sources[j] %in% limits){
      if(sources[j] == "prior"){
        x_lim_probs <- c(0, 0.95)
      } else{
        x_lim_probs <- c(0, 0.99)
      }
      quants <- aggregate(rt ~ group_key, datasets[[j]], quantile, probs = x_lim_probs)
      xlim <- range(xlim, unlist(quants$rt))
    }
  }
  datasets <- lapply(datasets, function(x){
    x <- x[x$rt > xlim[1] & x$rt < xlim[2],]
    return(x)
  })
  return(list(datasets = datasets, sources = sources, xlim = xlim))
}

#' Plot Statistics on Data
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
                      subject = NULL, factors = NULL, n_cores = 1, n_post = 50,
                      quants = c(0.025, 0.5, 0.975), functions = NULL,
                      layout = NA, to_plot = c('data', 'posterior', 'prior')[1:2], use_lim = c('data', 'posterior', 'prior')[1:2],
                      legendpos = c('topleft', 'top'), posterior_args = list(), prior_args = list(), ...) {
  check <- prep_data_plot(input, post_predict, prior_predict, to_plot, use_lim,
                          factors, defective_factor = NULL, subject, n_cores, n_post, functions)
  data_sources <- check$datasets
  sources <- check$sources

  # Basic definitions
  dots <- add_defaults(list(...), col = c("black",  "#A9A9A9", "#666666"))
  posterior_args <- add_defaults(posterior_args, col = c("darkgreen",  "#0000FF", "#008B8B"))
  prior_args <- add_defaults(prior_args, col = c("red", "#800080", "#CC00FF"))

  line_sources <- list()
  dens_sources <- list()
  summary_sources <- list()
  max_y <- .5
  xlim <- NULL

  # First big loop - calculating stats
  for (k in 1:length(data_sources)) {
    src_type <- sources[k]
    src_data <- data_sources[[k]]
    src_name <- names(sources)[k]
    if('postn' %in% colnames(src_data)){
      split_src <- split(src_data, list(src_data$group_key, src_data$postn))
      src_args <- ifelse(src_type == "posterior", posterior_args, prior_args)
      # Apply stat_fun to each group in src_data
      stats <- lapply(split_src, function(data_group) {
        stats <- stat_fun(data_group)
        factors_values <- data_group[1, factors, drop = FALSE]
        postn <- data_group$postn[1]
        group_key <- data_group$group_key[1]
        result <- cbind(factors_values, group_key = group_key, postn = postn, as.data.frame(t(stats)))
        return(result)
      })

      # Combine and compute quantiles
      src_stats_df <- do.call(rbind, stats)
      stat_names <- names(src_stats_df)[!(names(src_stats_df) %in% c('postn', 'group_key', factors))]

      quantile_stats <- lapply(stat_names, function(sname) {
        agg_result <- aggregate(as.formula(paste0("`", sname, "` ~ ", paste(c('group_key', factors), collapse = '+'))),
                                data = src_stats_df,
                                FUN = function(x) quantile(x, probs = quants))
        quantile_matrix <- agg_result[[sname]]
        colnames(quantile_matrix) <- paste0(sname, "_", src_name, "_", quants)
        result <- cbind(agg_result[, c('group_key', factors), drop = FALSE], quantile_matrix)
        return(result)
      })

      summary_sources[[src_name]] <- Reduce(function(x, y) merge(x, y, by = c('group_key', factors)),
                                 quantile_stats)
      dens <- lapply(split(src_stats_df, src_stats_df$group_key), function(x){
        lapply(stat_names, function(y){
          do.call(density, c(list(x[,y]), fix_dots(dots, density.default, consider_dots = F)))
        })
      })
      dens_sources[[src_name]] <- dens
      if(src_type %in% use_lim){
        max_y <- max(max_y, max(sapply(unlist(dens, recursive = F), function(x) return(max(x$y)))))
        xlim <- range(xlim, range(sapply(unlist(dens, recursive = F), function(x) return(quantile(x$x, probs = c(0.01, 0.99))))))
      }
    } else{
      split_src <- split(src_data, src_data$group_key)
      # Apply stat_fun to each group
      stat_results <- lapply(split_src, function(data_group) {
        stats <- stat_fun(data_group)
        factors_values <- data_group[1, factors, drop = FALSE]
        result <- cbind(factors_values, as.data.frame(t(stats)))
        return(result)
      })

      # Combine results into summary data frame
      summary_df <- do.call(rbind, stat_results)
      summary_df <- cbind(group_key = rownames(summary_df), summary_df)
      rownames(summary_df) <- NULL
      stat_names <- names(summary_df)[!(names(summary_df) %in% c('group_key', factors))]

      if('data' %in% use_lim){
        xlim <- range(xlim, summary_df[,stat_names])
      }
      line_sources[[src_name]] <- summary_sources[[src_name]] <- summary_df
    }
  }

  # Second big loop - Plotting loop
  first_data <- data_sources[[1]]
  if (is.null(first_data)) return(invisible(NULL))
  unique_group_keys <- unique(first_data$group_key)
  if(!is.null(layout)){
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }
  if (any(is.na(layout))) {
    par(mfrow = coda_setmfrow(Nchains = 1, Nparms = length(unique_group_keys), nplots = 1))
  } else {
    par(mfrow = layout)
  }
  lty <- 1:length(stat_names)

  ylim <- c(0, max_y)
  dots <- add_defaults(dots, lty = lty, ylim = ylim, xlim = xlim)
  for (group_key in unique_group_keys) {
    tmp_dots <- dots
    tmp_posterior_args <- posterior_args
    tmp_prior_args <- prior_args
    legend_map <- c()
    lwd_map <- numeric()
    main_in <- ifelse(is.null(stat_name), group_key, paste(group_key, "-", stat_name))
    do.call(plot, c(list(NA), fix_dots_plot(add_defaults(dots, main = main_in, ylab = "Density", xlab = " "))))
    for (k in 1:length(data_sources)) {
      src_type <- sources[k]
      src_data <- data_sources[[k]]
      src_name <- names(sources)[k]
      if(src_type == "data") {
        src_args <- tmp_dots
        tmp_dots$col <- tmp_dots$col[-1]
      } else if (src_type == "posterior") {
        src_args <- tmp_posterior_args
        tmp_posterior_args$col <- tmp_posterior_args$col[-1]
      } else if (src_type == "prior") {
        src_args <- tmp_prior_args
        tmp_prior_args$col <- tmp_prior_args$col[-1]
      }
      legend_map[src_name] <- src_args$col[1]
      lwd_map[src_name] <- ifelse(is.null(src_args$lwd), 1, src_args$lwd)
      lines_args <- fix_dots_plot(src_args)
      lines_args$col <- lines_args$col[1]
      if(src_type == 'data'){
        summary_df <- line_sources[[src_name]]
        do.call(abline, c(list(v = summary_df[summary_df$group_key == group_key, stat_names]), lines_args))
      } else{
        cur_dens <- dens_sources[[src_name]]
        mapply(cur_dens[[group_key]], lty, FUN = function(x, y) do.call(lines, c(list(x), lty = y,
                                                                                        lines_args[names(lines_args) != "lty"])))
      }
    }
    if(length(stat_names) > 1){
      legend(x = legendpos[1], legend = stat_names, lty = lty, col = "black", bty = "n")
    }
    if (length(data_sources) > 1) {
      legend(legendpos[2], legend = names(legend_map), lty = 1, col = legend_map, title = "Source", bty = "n", lwd = lwd_map)
    }
  }
  return(invisible(summary_df[,!grepl("group_key", colnames(summary_df))]))
}




# A small function to compute the defective densities across factor levels
compute_def_dens <- function(dat, defective_factor, dargs) {
  p_defective <- prop.table(table(dat[[defective_factor]]))
  # We'll call density() on each subset of rt, then multiply by proportion
  # so that the sum across factor levels is 1
  # We'll use the from/to in dargs
  by_deflev <- split(dat, dat[[defective_factor]])
  out <- list()
  for (lev in names(by_deflev)) {
    subdat <- by_deflev[[lev]]
    if (nrow(subdat) < 2) {
      # avoid error
      out[[lev]] <- rep(0, 512)
    } else {
      dd <- do.call(density, c(list(x = subdat$rt), fix_dots(dargs, density.default, consider_dots = FALSE)))
      out[[lev]] <- dd$y * p_defective[lev]
    }
  }
  out
}


#' Plot Defective Densities
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
plot_density <- function(input, post_predict = NULL, prior_predict = NULL,
                         subject = NULL, quants = c(0.025, 0.975), functions = NULL,
                         factors = NULL, defective_factor = "R", n_cores = 1, n_post = 50,
                         layout = NA,
                         to_plot = c('data', 'posterior', 'prior')[1:2],
                         use_lim = c('data', 'posterior', 'prior')[1:2],
                         legendpos = c("topright", "top"),
                         posterior_args = list(), prior_args = list(), ...) {
  # 1) prep_data_plot
  check <- prep_data_plot(input, post_predict, prior_predict, to_plot, use_lim,
                          factors, defective_factor, subject, n_cores, n_post, functions)
  data_sources <- check$datasets
  sources <- check$sources
  xlim <- check$xlim

  # Basic definitions
  dots <- add_defaults(list(...), col = c("black",  "#A9A9A9", "#666666"))
  posterior_args <- add_defaults(posterior_args, col = c("darkgreen",  "#0000FF", "#008B8B"))
  prior_args <- add_defaults(prior_args, col = c("red", "#800080", "#CC00FF"))

  data <- data_sources[[which(sources == "data")[1]]]
  defective_levels <- if (!is.null(data)) levels(factor(data[[defective_factor]])) else character(0)
  line_types <- seq_along(defective_levels)
  quantiles <- sort(c(quants, 0.5))

  # We'll keep track of the defective density results:
  # For single dataset => a single vector of y for each level
  # For postn => we keep them all, then compute quantiles
  dens_list <- list()        # For main
  dens_quants_list <- list() # For quantiles
  y_max <- 0

  # 2) First big loop: compute main density or quantiles
  for (k in 1:length(data_sources)) {
    src_data <- data_sources[[k]]
    src_type <- sources[k]
    src_name <- names(sources)[k]
    if (is.null(src_data)) next

    # If 'postn' in colnames => multiple sets => need quantiles
    if ("postn" %in% names(src_data)) {
      # We'll compute from/to from data range
      rng <-  quantile(src_data$rt, .99)
      dargs <- switch(
        src_type,
        "posterior" = add_defaults(posterior_args, to = rng),
        "prior"     = {
          rng2 <- quantile(src_data$rt, c(0.975))
          add_defaults(prior_args, to = rng2)
        },
        dots
      )

      splitted <- split(src_data, src_data$group_key)
      # splitted[group_key][[postn]] => we also need to sub-split by postn
      dens_list[[src_name]] <- lapply(splitted, function(dg) {
        postn_splits <- split(dg, dg$postn)
        lapply(postn_splits, function(dsub) {
          compute_def_dens(dsub, defective_factor, dargs)
        })
      })

      # now compute quantiles over postn
      dens_quants_list[[src_name]] <- lapply(dens_list[[src_name]], function(d_list_of_lists) {
        # d_list_of_lists => list of postn => each is a named list by factor level
        out <- list()
        for (lev in defective_levels) {
          # gather each postn's vector of y for this level => cbind
          mat_y <- do.call(cbind, lapply(d_list_of_lists, function(ll) ll[[lev]]))
          # compute row-wise quantile
          # mat_y is 512 x #postn
          if (!is.null(mat_y)) {
            # dimension checks
            q_y <- apply(mat_y, 1, quantile, probs = quantiles, na.rm = TRUE)
            out[[lev]] <- q_y
          } else {
            out[[lev]] <- NULL
          }
        }
        out
      })

      # if we use this source for y-lim
      if (src_type %in% use_lim) {
        # find max of upper quantile
        all_vals <- unlist(lapply(dens_quants_list[[src_name]], function(lv) {
          sapply(lv, function(m) if (!is.null(m)) max(m, na.rm = TRUE) else 0)
        }))
        y_max <- max(y_max, all_vals)
      }

    } else {
      # Single dataset
      rng <- quantile(src_data$rt, .99)
      dargs <- add_defaults(dots, to = rng)
      splitted <- split(src_data, src_data$group_key)
      dens_list[[src_name]] <- lapply(splitted, compute_def_dens, defective_factor, dargs)

      if (src_type %in% use_lim) {
        # find max
        out_max <- max(unlist(dens_list[[src_name]]), na.rm = TRUE)
        y_max <- max(y_max, out_max)
      }
    }
  }

  # 3) Second big loop: plotting
  if(!is.null(layout)){
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
  }
  # figure out group keys from whichever data source is not null
  # pick the first in data_sources
  first_data <- data_sources[[1]]
  if (is.null(first_data)) {
    # If there's absolutely no data, just return
    return(invisible(NULL))
  }
  unique_group_keys <- unique(first_data$group_key)

  if (any(is.na(layout))) {
    par(mfrow = coda_setmfrow(Nchains = 1, Nparms = length(unique_group_keys), nplots = 1))
  } else {
    par(mfrow = layout)
  }
  ylim_global <- c(0, y_max)

  # We also need an x-grid for each data source
  # For single dataset (no postn), we can store the computed from/to in dens_list.
  # But we haven't stored x-coordinates. We can reconstruct from the 'from','to' used in density.
  # The simplest approach is to do the same method as the original code: create a grid from the min to max.
  # We'll rely on the from/to used above.

  # function to retrieve the from/to arguments
  get_range_args <- function(dargs) {
    c(dargs$from, dargs$to)
  }

  # We'll define a small function to figure out the x-grid used for each data source
  # We do that at plotting time to keep the same 512 points
  build_x_grid <- function(src_type, src_name, group_key = NULL) {
    # for postn or single data we found a from/to in a dictionary
    if (src_type %in% c("data", 'posterior')) {
      # if data was used, dens_list[['data']][[group_key]] => each level => length=512
      # but we didn't store the x's. We'll reconstruct from the dots or from
      rng <- quantile(data_sources[[src_name]]$rt, c(0, 0.99))
      return(seq(rng[1], rng[2], length.out = 512))
    } else if (src_type == "prior") {
      rng2 <- quantile(data_sources[[src_name]]$rt, c(0.001, 0.975))
      return(seq(rng2[1], rng2[2], length.out = 512))
    }
  }


  for (group_key in unique_group_keys) {
    tmp_dots <- dots
    tmp_posterior_args <- posterior_args
    tmp_prior_args <- prior_args
    # empty plot
    plot_args <- add_defaults(dots, xlim = xlim, ylim = ylim_global,
                              main = group_key, xlab = "RT", ylab = "Defective Density")
    plot_args <- fix_dots_plot(plot_args)
    do.call(plot, c(list(NA), plot_args))
    legend_map <- c()
    lwd_map <- numeric()
    for (k in 1:length(data_sources)) {
      src_type <- sources[k]
      src_data <- data_sources[[k]]
      src_name <- names(sources)[k]
      if(src_type == "data") {
        src_args <- tmp_dots
        tmp_dots$col <- tmp_dots$col[-1]
      } else if (src_type == "posterior") {
        src_args <- tmp_posterior_args
        tmp_posterior_args$col <- tmp_posterior_args$col[-1]
      } else if (src_type == "prior") {
        src_args <- tmp_prior_args
        tmp_prior_args$col <- tmp_prior_args$col[-1]
      }
      legend_map[src_name] <- src_args$col[1]
      lwd_map[src_name] <- ifelse(is.null(src_args$lwd), 1, src_args$lwd)
      # if no postn => single dataset
      if (!("postn" %in% colnames(src_data))) {
        # dens_list[[src_name]][[group_key]] => list of defective_levels => vector of length=512
        cur_dens <- dens_list[[src_name]][[group_key]]
        if (!is.null(cur_dens)) {
          x_grid <- build_x_grid(src_type, src_name, group_key)
          for (i in seq_along(defective_levels)) {
            lev <- defective_levels[i]
            yvals <- cur_dens[[lev]]
            if (!is.null(yvals)) {
              lines_args <- add_defaults(src_args, lty = line_types[i])
              lines_args <- fix_dots_plot(lines_args)
              do.call(lines, c(list(x = x_grid, y = yvals), lines_args))
            }
          }
        }
      } else {
        # we have postn => look up dens_quants_list
        cur_quants <- dens_quants_list[[src_name]][[group_key]]
        if (!is.null(cur_quants)) {
          x_grid <- build_x_grid(src_type, src_name, group_key)
          for (i in seq_along(defective_levels)) {
            lev <- defective_levels[i]
            m <- cur_quants[[lev]]
            if (!is.null(m)) {
              # m => matrix of size length(quantiles) x 512
              # typically row 1 => lower, row 2 => median, row 3 => upper
              # Let's define row indexes:
              y_lower <- m[1,]
              y_med   <- m[2,]
              y_upper <- m[3,]

              lines_args <- add_defaults(src_args, lty = line_types[i])
              lines_args <- fix_dots_plot(lines_args)
              do.call(lines, c(list(x = x_grid, y = y_med), lines_args))

              # polygon
              adj_color <- do.call(adjustcolor, fix_dots(add_defaults(src_args, alpha.f = .2), adjustcolor))
              polygon_args <- src_args
              polygon_args$col <- adj_color
              polygon_args <- fix_dots_plot(polygon_args)
              do.call(polygon, c(list(
                x = c(x_grid, rev(x_grid)),
                y = c(y_lower, rev(y_upper)),
                border = NA
              ), polygon_args))
            }
          }
        }
      }
    }
    # Add legends
    legend(legendpos[1], legend = defective_levels, lty = line_types, col = "black",
           title = defective_factor, bty = "n")
    if (length(data_sources) > 1) {
      legend(legendpos[2], legend = names(legend_map), lty = 1, col = legend_map, title = "Source", bty = "n",
             lwd = lwd_map)
    }
  }
}

###############################################################################
## Helper: get_def_cdf
###############################################################################
get_def_cdf <- function(x, defective_factor, dots) {
  # Computes a single defective CDF for each level of 'defective_factor'
  # across the RT distribution (0.01 to 0.99).
  probs <- seq(0.01, 0.99, by = 0.01)
  p_defective <- prop.table(table(x[[defective_factor]]))

  # For each level, compute the empirical CDF and scale by that level's proportion
  out <- mapply(
    split(x, x[[defective_factor]]),
    p_defective,
    FUN = function(inp, prop_share) {
      # quantile of RTs
      rtquants <- quantile(inp$rt, probs = probs, type = 1, na.rm = TRUE)
      # defective cdf = empirical cdf (probs) times the proportion
      yvals <- probs * prop_share
      cbind(x = rtquants, y = yvals)
    },
    SIMPLIFY = FALSE
  )
  # 'out' is now a list whose names are the factor levels; each entry is a 2-col matrix (x,y).
  return(out)
}

###############################################################################
## Plot Defective CDFs
###############################################################################
#' Plot Defective Cumulative Distribution Functions
#'
#' Plots panels of cumulative distribution functions (CDFs) for each level of the specified
#' defective factor in the data. The CDFs are *defective*; each factor level's CDF
#' scales only up to that level's proportion. Summed across levels, the maximum is 1.
#' Optionally, posterior and/or prior predictive CDFs can be overlaid.
#'
#' @param input Either an `emc` object or a data frame, or a *list* of such objects.
#' @param post_predict Optional posterior predictive data (matching columns) or *list* thereof.
#' @param prior_predict Optional prior predictive data (matching columns) or *list* thereof.
#' @param subject Subset the data to a single subject (by index or name).
#' @param quants Numeric vector of credible interval bounds (e.g. `c(0.025, 0.975)`).
#' @param functions A function (or list of functions) that create new columns in the datasets or predictives
#' @param factors Character vector of factor names to aggregate over;
#' defaults to plotting full data set ungrouped by factors if `NULL`.
#' @param defective_factor Name of the factor used for the defective CDF (default "R").
#' @param n_cores Number of CPU cores to use if generating predictives from an `emc` object.
#' @param n_post Number of posterior draws to simulate if needed for predictives.
#' @param layout Numeric vector used in `par(mfrow=...)`; use `NA` for auto-layout.
#' @param to_plot Character vector: any of `"data"`, `"posterior"`, `"prior"`.
#' @param use_lim Character vector controlling which source(s) define `xlim`.
#' @param legendpos Character vector controlling the positions of the legends
#' @param posterior_args Optional list of graphical parameters for posterior lines/ribbons.
#' @param prior_args Optional list of graphical parameters for prior lines/ribbons.
#' @param ... Other graphical parameters for the real data lines.
#'
#' @return Returns `NULL` invisibly.
#' @examples
#' # Plot defective CDF for data only
#' # plot_cdf(forstmann, to_plot = "data")
#' #
#' # Plot with posterior predictions
#' # plot_cdf(samples_LNR, to_plot = c("data","posterior"), n_post=10)
#' #
#' # Or a list of multiple emc objects ...
#' @export
plot_cdf <- function(input,
                     post_predict = NULL,
                     prior_predict = NULL,
                     subject = NULL,
                     quants = c(0.025, 0.975), functions = NULL,
                     factors = NULL,
                     defective_factor = "R",
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
                          factors, defective_factor, subject, n_cores, n_post, functions)
  data_sources <- check$datasets
  sources <- check$sources
  xlim <- check$xlim

  # Basic definitions
  dots <- add_defaults(list(...), col = c("black",  "#A9A9A9", "#666666"))
  posterior_args <- add_defaults(posterior_args, col = c("darkgreen",  "#0000FF", "#008B8B"))
  prior_args <- add_defaults(prior_args, col = c("red", "#800080", "#CC00FF"))

  defective_levels <- levels(factor(data_sources[[1]][[defective_factor]]))
  unique_group_keys <- levels(factor(data_sources[[1]]$group_key))

  if (is.null(defective_levels) || length(defective_levels) == 0) {
    defective_levels <- "Level1"  # fallback
  }
  line_types <- seq_along(defective_levels)

  # We'll store single-CDF results or multi-postn quantile results in these:
  cdf_list        <- list() # single
  cdf_quants_list <- list() # multi

  # We'll keep track of a global maximum in the vertical dimension
  # so that we can set a consistent y-lim across all panels
  y_max <- 0

  # -------------------------------------------------------------------
  # 2) FIRST BIG LOOP: compute CDF or CDF-quantiles for each dataset
  # -------------------------------------------------------------------
  for (k in seq_along(data_sources)) {
    df   <- data_sources[[k]]
    styp <- sources[k]       # "data","posterior","prior"
    sname <- names(sources)[k]  # the name of this dataset in the list

    if (is.null(df) || !nrow(df)) {
      # skip empty
      next
    }

    # We'll group by group_key
    splitted <- split(df, df$group_key)

    # If there's a "postn" column => multiple draws => compute quantiles
    if ("postn" %in% names(df)) {
      # cdf_list[[sname]] => list of (group_key => list of postn => single-cdf)
      cdf_list[[sname]] <- lapply(splitted, function(sub_grp) {
        # sub_grp is all rows for a single group_key
        # we further split by postn
        postn_splits <- split(sub_grp, sub_grp$postn)
        lapply(postn_splits, get_def_cdf, defective_factor, dots)
      })

      # Now we derive cdf_quants_list from cdf_list
      cdf_quants_list[[sname]] <- lapply(cdf_list[[sname]], function(postn_list) {
        # postn_list => e.g. 100 draws => each draw is a named list of factor-level => cbind(x,y)
        out <- list()
        for (lev in defective_levels) {
          # gather all x and y columns across draws
          x_mat <- do.call(cbind, lapply(postn_list, function(lst) lst[[lev]][,"x"]))
          y_mat <- do.call(cbind, lapply(postn_list, function(lst) lst[[lev]][,"y"]))

          # row-wise quantiles for x, plus median for y
          # Just as in your older code: we do quantiles on x across draws, median on y
          # Or you might do quantiles on both x and y.
          # Typically, to replicate your older approach:
          # quantile of x at each index, and median of y at each index
          # We'll include 50% in quants to do median for x as well.
          # Then we combine them in a matrix with 4 rows => x_lower, x_median, x_upper, y_median
          qx <- apply(x_mat, 1, quantile, probs = sort(c(quants, 0.5)), na.rm = TRUE)
          ym <- apply(y_mat, 1, median, na.rm=TRUE)

          # rbind them
          # row 1 => x_lower
          # row 2 => x_mid (0.5)
          # row 3 => x_upper
          # row 4 => y_median
          # (If you used 3 quantiles, that's 3 rows for x, plus 1 for y.)
          out[[lev]] <- rbind(qx, ym)
        }
        out
      })

      # If we use this dataset to define y-limit ...
      if (styp %in% use_lim) {
        # maximum of y_median row
        possible_vals <- unlist(lapply(cdf_quants_list[[sname]], function(group_val) {
          # group_val => list( factor_level => matrix(4 x length-of-probs) )
          sapply(group_val, function(mat4) {
            if (is.null(mat4)) return(0)
            max(mat4[nrow(mat4), ], na.rm=TRUE)  # last row => y_median
          })
        }))
        y_max <- max(y_max, possible_vals, na.rm=TRUE)
      }

    } else {
      # single dataset => cdf_list[[sname]] => group_key => get_def_cdf => named list by factor level
      cdf_list[[sname]] <- lapply(splitted, get_def_cdf, defective_factor, dots)

      # If we use this dataset for y-limit, find max
      if (styp %in% use_lim) {
        # find max y across all group_keys & factor-levels
        max_val <- 0
        for (grp_name in names(cdf_list[[sname]])) {
          cdf_grp <- cdf_list[[sname]][[grp_name]]
          if (!is.null(cdf_grp)) {
            max_val <- max(max_val, unlist(lapply(cdf_grp, function(m) max(m[,"y"]))), na.rm=TRUE)
          }
        }
        y_max <- max(y_max, max_val)
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

  for (group_key in unique_group_keys) {
    tmp_dots <- dots
    tmp_posterior_args <- posterior_args
    tmp_prior_args <- prior_args

    # blank plot
    plot_args <- add_defaults(dots, xlim=xlim, ylim=ylim,
                              main=group_key, xlab="RT", ylab="Defective CDF")
    plot_args <- fix_dots_plot(plot_args)
    do.call(plot, c(list(NA), plot_args))

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
      # if no postn => single cdf => cdf_list[[sname]][[group_key]] => factor-level => matrix(x,y)
      # if postn => cdf_quants_list[[sname]][[group_key]] => factor-level => matrix(4 x length-probs)

      # We might not have a cdf if that group_key wasn't present in the data
      # so check if cdf_list has an entry
      if (!is.null(cdf_list[[sname]])) {
        # check for group group_key
        if (!("postn" %in% names(data_sources[[k]]))) {
          # single
          cdf_for_group <- cdf_list[[sname]][[group_key]]
          if (!is.null(cdf_for_group)) {
            # cdf_for_group => e.g. list( factor-level => matrix of x,y )
            ilev <- 1
            for (lev in defective_levels) {
              cmat <- cdf_for_group[[lev]]
              if (!is.null(cmat)) {
                # lines
                lines_args <- add_defaults(src_args, lty=line_types[ilev])
                lines_args <- fix_dots_plot(lines_args)
                do.call(lines, c(list(x=cmat[,"x"], y=cmat[,"y"]), lines_args))
              }
              ilev <- ilev+1
            }
          }
        } else {
          # multi draws => cdf_quants_list
          if (!is.null(cdf_quants_list[[sname]])) {
            cdf_quants_for_group <- cdf_quants_list[[sname]][[group_key]]
            if (!is.null(cdf_quants_for_group)) {
              ilev <- 1
              for (lev in defective_levels) {
                mat4 <- cdf_quants_for_group[[lev]]
                if (!is.null(mat4)) {
                  # mat4 => e.g. 4 rows x (length(probs)) columns: row1 => x_lower, row2 => x_mid, row3 => x_upper, row4 => y_median
                  x_lower <- mat4[1,]
                  x_med   <- mat4[2,]
                  x_upper <- mat4[3,]
                  y_median<- mat4[4,]

                  lines_args <- add_defaults(src_args, lty=line_types[ilev])
                  lines_args <- fix_dots_plot(lines_args)
                  do.call(lines, c(list(x=x_med, y=y_median), lines_args))

                  # polygon for the ribbon
                  adj_color <- do.call(adjustcolor, fix_dots(add_defaults(src_args, alpha.f=0.2), adjustcolor))
                  poly_args <- src_args
                  poly_args$col <- adj_color
                  poly_args <- fix_dots_plot(poly_args)
                  # We connect (x_lower, y_median) -> (x_upper, rev(y_median))
                  # (the original code uses horizontal ribbons).
                  do.call(polygon, c(list(
                    x = c(x_lower, rev(x_upper)),
                    y = c(y_median, rev(y_median)),
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

    # Factor-level legend
    legend(legendpos[1], legend=defective_levels, lty=line_types, col="black",
           title=defective_factor, bty="n")

    # If multiple data sources, show source legend
    if (length(data_sources) > 1) {
      legend(legendpos[2], legend=names(legend_map), lty=1, col=legend_map,
             title="Source", bty="n", lwd = lwd_map)
    }
  } # end for each group_key

  invisible(NULL)
}

