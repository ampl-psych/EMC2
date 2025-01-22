# ----------------------------
# 1) Utility Functions
# ----------------------------

# 1.1 Compute Parameter Means by Factor
get_param_means <- function(pars, factor_name, param_name) {
  if (is.null(factor_name)) {
    # Single global mean
    setNames(mean(pars[[param_name]], na.rm = TRUE), "global")
  } else {
    # Means by factor level
    means <- tapply(pars[[param_name]], pars[[factor_name]], mean, na.rm = TRUE)
    # Handle cases where tapply might return NA due to all NA values
    means[is.na(means)] <- mean(pars[[param_name]], na.rm = TRUE)
    means
  }
}

# 1.2 Assign Colors to Factor Levels
assign_colors <- function(levels, single_col = "#036ffc") {
  if (length(levels) <= 1) {
    return(setNames(single_col, names(levels)))
  }

  base_cols <- c("#8C2E89", "#009414", "#D68E07", "#200BA5",
                 "#02C6B2", "#C60202", "#788714", "#D33A7C")
  n <- length(levels)
  if (n > length(base_cols)) {
    base_cols <- rep(base_cols, length.out = n)
  }
  setNames(base_cols[1:n], names(levels))
}

# 1.3 Assign Line Types to Factor Levels
assign_ltys <- function(levels) {
  if (length(levels) <= 1) {
    return(setNames(1, names(levels)))
  }

  ltys <- c(1, 3:(length(levels) + 1))  # Start with solid, then dash
  setNames(ltys, names(levels))
}

# 1.4 Check Presence of Factor Level in Subset
level_in_subset <- function(dat, factor_name, level_name) {
  if (is.null(factor_name)) return(TRUE)
  if (!factor_name %in% colnames(dat)) return(TRUE)
  any(dat[[factor_name]] == level_name, na.rm = TRUE)
}

# 1.5 Determine Parameter Value for a Cell
get_cell_param <- function(cell_idx, factor_name, param_vals, data_sub) {
  if (is.null(factor_name)) {
    return(param_vals["global"])
  }
  if (!factor_name %in% colnames(data_sub)) {
    return(mean(param_vals, na.rm = TRUE))
  }
  factor_levels <- unique(data_sub[[factor_name]][cell_idx])
  factor_level <- if (length(factor_levels) > 1) factor_levels[1] else factor_levels
  if (!factor_level %in% names(param_vals)) {
    return(mean(param_vals, na.rm = TRUE))
  } else {
    return(param_vals[factor_level])
  }
}

# 1.6 Determine Plot Ranges Across Subsets
# 1.6 Determine Plot Ranges Across Subsets
get_subset_ranges <- function(
    data,
    b_factor, v_factor, t0_factor,
    b_vals,
    squash_dens = 1,
    split_factor = NULL
) {
  # If empty => trivial
  if (nrow(data) == 0) {
    return(list(
      max_rt_sub     = 1,
      min_y_sub      = 0,
      max_y_sub      = 0,
      defect_scalars = numeric(0)
    ))
  }

  # 1) Time dimension => 99th percentile
  max_rt_sub <- quantile(data$rt, 0.99, na.rm = TRUE)
  if (!is.finite(max_rt_sub)) {
    max_rt_sub <- 1
  }

  # 2) Identify factor columns (including split_factor!)
  all_factors_sub <- unique(c(v_factor, b_factor, t0_factor, split_factor))
  all_factors_sub <- all_factors_sub[!is.na(all_factors_sub)]
  data_factors_in_data <- intersect(all_factors_sub, colnames(data))

  # 3) Build "cells" (unmirrored) => includes split_factor if present
  if (length(data_factors_in_data) > 0) {
    cells <- apply(data[, data_factors_in_data, drop=FALSE], 1, function(r) {
      paste(mapply(function(f, v) paste0(f, "=", v),
                   data_factors_in_data, r),
            collapse=" ")
    })
  } else {
    cells <- rep("global", nrow(data))
  }
  unique_cells <- unique(cells)

  n_data         <- nrow(data)
  defect_scalars <- numeric(length(unique_cells))
  names(defect_scalars) <- unique_cells

  min_y_sub <-  Inf
  max_y_sub <- -Inf

  for (uc in unique_cells) {
    cell_idx <- which(cells == uc)
    d_rt     <- data$rt[cell_idx]

    # density
    if (length(d_rt) < 2) {
      d <- list(x=c(0, d_rt), y=c(0,0))
    } else {
      d <- density(d_rt, na.rm=TRUE)
    }

    # b_val for the cell (unmirrored)
    if (is.null(b_factor)) {
      b_val <- mean(b_vals, na.rm = TRUE)
    } else {
      factor_levels <- unique(data[[b_factor]][cell_idx])
      factor_level  <- if (length(factor_levels)>1) factor_levels[1] else factor_levels
      if (!is.null(factor_level) && factor_level %in% names(b_vals)) {
        b_val <- b_vals[factor_level]
      } else {
        b_val <- mean(b_vals, na.rm=TRUE)
      }
    }

    # proportion => defective scaling
    cell_prop           <- length(cell_idx) / n_data
    defect_scalars[uc]  <- cell_prop

    # scale => no sign multiplier
    # d$y => (density * squash_dens * cell_prop) + b_val
    d$y <- d$y * squash_dens * cell_prop + b_val

    local_min <- min(d$y, na.rm=TRUE)
    local_max <- max(d$y, na.rm=TRUE)

    if (local_min < min_y_sub) min_y_sub <- local_min
    if (local_max > max_y_sub) max_y_sub <- local_max
  }

  # if still infinite => set them to 0
  if (!is.finite(min_y_sub)) min_y_sub <- 0
  if (!is.finite(max_y_sub)) max_y_sub <- 0

  list(
    max_rt_sub     = max_rt_sub,
    min_y_sub      = min_y_sub,
    max_y_sub      = max_y_sub,
    defect_scalars = defect_scalars
  )
}




# 1.7 Determine t0 y-position
get_t0_ypos <- function(t0_val, t0_vals, data_sub, t0_factor,
                        t0_levels_sorted, vertical_step) {
  tn <- names(t0_vals)[ which(t0_vals==t0_val) ]
  if (!level_in_subset(data_sub, t0_factor, tn)) {
    return(0)
  }
  idx <- match(tn, names(t0_levels_sorted))
  if (is.na(idx)) return(0)
  t0_in_sub_sorted <- names(t0_levels_sorted)[
    sapply(names(t0_levels_sorted), function(x) level_in_subset(data_sub, t0_factor, x))
  ]
  idx2 <- match(tn, t0_in_sub_sorted) -1
  if (is.na(idx2) || idx2<0) return(0)
  idx2*vertical_step
}

# ----------------------------
# 2) Single DDM Plot Function
# ----------------------------
single_DDM_plot <- function(
    data_sub,
    main_sub,
    plot_legend_sub,
    x_lim,
    y_lim,
    v_factor,
    b_factor,
    t0_factor,
    v_levels,
    b_vals,
    t0_vals,
    v_colors,
    t0_ltys,
    split,
    draw_axis_labels = TRUE,
    within_noise     = TRUE,
    squash_dens,
    defect_scalars = NULL
) {
  ###########################################################
  ## (A) Split data by 'split' => mirror densities
  ###########################################################
  if (!is.null(split) && split %in% names(data_sub)) {
    split_levels <- levels(data_sub[[split]])
    if (is.null(split_levels)) {
      split_levels <- sort(unique(data_sub[[split]]))
    }
  } else {
    split_levels <- "ALL_SPLIT"
  }
  bottom_level <- split_levels[1]
  top_level    <- if (length(split_levels) > 1) split_levels[2] else bottom_level

  data_bottom <- data_sub[data_sub[[split]] == bottom_level, ]
  data_top    <- data_sub[data_sub[[split]] == top_level,    ]

  ###########################################################
  ## (B) Build mirrored densities
  ###########################################################
  build_cell_info_DDM <- function(dat, sign_mult) {
    # we incorporate the 'split' factor if it exists
    all_factors_sub <- c(v_factor, b_factor, t0_factor, split)
    all_factors_sub <- unique(all_factors_sub[!is.na(all_factors_sub)])
    data_factors_in_data <- intersect(all_factors_sub, colnames(dat))

    if (length(data_factors_in_data) > 0) {
      row_cells <- apply(dat[, data_factors_in_data, drop=FALSE], 1, function(r) {
        paste(mapply(function(f, v) paste0(f, "=", v),
                     data_factors_in_data, r),
              collapse=" ")
      })
    } else {
      row_cells <- rep("global", nrow(dat))
    }
    unique_cells <- unique(row_cells)
    cell_list   <- lapply(unique_cells, function(uc) {
      idx <- which(row_cells == uc)
      tmp_data <- dat[idx, ]

      v_cell  <- get_cell_param(idx, v_factor,  v_levels, dat)
      b_cell  <- get_cell_param(idx, b_factor,  b_vals,   dat)
      t0_cell <- get_cell_param(idx, t0_factor, t0_vals,  dat)
      med_rt  <- median(tmp_data$rt, na.rm=TRUE)

      # Build density
      if (nrow(tmp_data) < 2) {
        d <- list(x = c(0, tmp_data$rt), y = c(0, 0))
      } else {
        d <- density(tmp_data$rt, na.rm=TRUE)
      }

      # Mirror => multiply y by sign_mult, shift by sign_mult*b_cell
      cell_prop <- defect_scalars[[uc]]
      d$y <- sign_mult * (d$y * squash_dens * cell_prop) + sign_mult*b_cell
      list(
        cell    = uc,
        data    = tmp_data,
        v       = v_cell,
        b       = b_cell,
        t0      = t0_cell,
        density = d,
        med_rt  = med_rt,
        sign    = sign_mult
      )
    })
    cell_list
  }

  cell_info_bottom <- build_cell_info_DDM(data_bottom, -1)
  cell_info_top    <- build_cell_info_DDM(data_top,    +1)
  all_cell_info    <- c(cell_info_bottom, cell_info_top)

  ###########################################################
  ## (C) Setup Plot
  ###########################################################
  plot(0, 0, type="n",
       xlim=x_lim, ylim=y_lim,
       main=main_sub,
       axes=FALSE, xlab="", ylab="")

  # x-axis near bottom
  base_line_y <- y_lim[1] - 0.05
  arrows(x0=0, x1=x_lim[2],
         y0=base_line_y, y1=base_line_y,
         lwd=3, length=0.1)

  ###########################################################
  ## (D) ±b lines
  ###########################################################
  # For each b factor level, draw +b, -b
  use_b <- setNames(rep(T, length(b_vals)), names(b_vals))
  for (bn in names(b_vals)) {
    if (!level_in_subset(data_sub, b_factor, bn)){
      use_b[bn] <- F
      next
    }
    b_here <- b_vals[bn]
    segments(0, +b_here, x_lim[2], +b_here,
             lwd=1.5, col="black")
    segments(0, -b_here, x_lim[2], -b_here,
             lwd=1.5, col="black")
  }
  b_vals <- b_vals[use_b]
  ###########################################################
  ## (E) Sort by median RT => layering
  ###########################################################
  med_rts_sub  <- sapply(all_cell_info, `[[`, "med_rt")
  rt_order     <- order(med_rts_sub)
  all_cell_info<- all_cell_info[rt_order]

  min_diff_sub <- 0.05*x_lim[2]
  for (i in seq_along(all_cell_info)[-1]) {
    if ((all_cell_info[[i]]$med_rt - all_cell_info[[i-1]]$med_rt) < min_diff_sub) {
      all_cell_info[[i]]$med_rt <- all_cell_info[[i-1]]$med_rt + min_diff_sub
    }
  }

  ###########################################################
  ## (F) t0 lines
  ###########################################################
  t0_levels_sorted <- sort(t0_vals)
  vertical_step    <- 0.015*(y_lim[2] - y_lim[1])
  i <- 0
  t0_y_vals <- numeric(length(t0_levels_sorted))
  for (tn in names(t0_levels_sorted)) {
    if (!level_in_subset(data_sub, t0_factor, tn)) next
    y_i  <- i*vertical_step
    t0_i <- t0_vals[tn]
    segments(0, y_i, t0_i, y_i,
             lwd=3, lty=t0_ltys[tn], col="black")
    i <- i+1
    t0_y_vals[i] <- y_i
  }
  if (i>0) {
    text(mean(t0_vals)/3, (i-1)*vertical_step + vertical_step/2,
         "t0", cex=1.5, col="black", adj=c(0,0))
  }

  ###########################################################
  ## (G) Densities
  ###########################################################
  alpha<-0.2
  for (ci in all_cell_info) {
    # line type from t0
    t0_level_name <- names(t0_vals)[ which(t0_vals==ci$t0) ]
    c_lty <- if (length(t0_level_name)>0) t0_ltys[t0_level_name] else 1

    # color from v
    v_level_name <- names(v_levels)[ which(v_levels==ci$v) ]
    if (is.null(v_factor)) {
      c_col   <- "black"
      fill_col<- "#036ffc"
    } else {
      if (length(v_level_name)==0) {
        c_col   <- "black"
        fill_col<- "#036ffc"
      } else {
        c_col   <- v_colors[v_level_name]
        fill_col<- v_colors[v_level_name]
      }
    }

    lines(ci$density, lwd=3, col=c_col, lty=c_lty)
    polygon(ci$density, col=adjustcolor(fill_col, alpha.f=alpha), lty=c_lty)
  }

  #################################################################
  ## (H) Drift lines => "device-based slope"
  ##
  ## Start:
  ##   x0 = max(t0_vals),
  ##   y0 = mean(t0_vals)
  ##
  ## End:
  ##   if v>0 => y1 = +min(b_vals)
  ##   if v<0 => y1 = -min(b_vals)
  ##
  ## slope_device = drift_scale * v
  ## slope_data   = slope_device * (w_plot/h_plot)
  ## dx           = (y1 - y0)/ slope_data
  ## x1           = x0 + dx
  #################################################################
  w_plot <- par("pin")[1] / diff(par("usr")[1:2])
  h_plot <- par("pin")[2] / diff(par("usr")[3:4])
  drift_scale <- 1
  boundary_val <- if (length(b_vals)>0) min(b_vals) else 0

  for (ci in all_cell_info) {
    v_here <- ci$v
    if (is.na(v_here) || abs(v_here) < 1e-12) next

    # color
    v_level_name <- names(v_levels)[ which(v_levels==v_here) ]
    c_col <- if (length(v_level_name)==0) "black" else v_colors[v_level_name]

    x0 <- max(t0_vals)
    y0 <- mean(t0_y_vals)
    # End y
    # if v>0 => y1= +boundary_val
    # if v<0 => y1= -boundary_val
    y1 <- if (v_here>0) boundary_val else -boundary_val

    # slope_device
    slope_device <- drift_scale * v_here
    slope_data   <- slope_device * (w_plot / h_plot)

    dy <- (y1 - y0)
    if (abs(slope_data) < 1e-12) next
    dx <- dy / slope_data
    x1 <- x0 + dx
    arrows(x0, y0 = y0, x1 = x1, y1 = y1,
           lwd=3, length=0.1, col=c_col)
    if (within_noise) {
      draw_noise_paths_DDM(
        x0 = x0, x1 = x1,
        y0 = y0, y1 = y1,
        col=c_col
      )
    }
  }

  #################################################################
  ## (I) Legend
  #################################################################
  if (plot_legend_sub) {
    x_legend<- x_lim[2]*0.3
    y_legend<- min(b_vals) * .8
    if (length(v_levels)>1) {
      legend(x=x_legend, y=y_legend,
             legend=names(v_levels),
             title="v",
             col=v_colors[names(v_levels)],
             lwd=3, lty=1, bty="n", cex = .9)
    }
    if (length(t0_vals)>1) {
      x_legend<- x_lim[2]*0.6
      legend(x=x_legend, y=y_legend,
             legend=names(t0_vals),
             title="t0",
             col="black", lwd=3,
             lty=t0_ltys[names(t0_vals)],
             bty="n", cex = .9)
    }
  }

  #################################################################
  ## (J) Axes
  #################################################################
  arrows(x0=-0.025, x1=-0.025,
         y0=y_lim[1], y1=y_lim[2],
         lwd=3, length=0.1)
  mtext("evidence", side=2, line=0.5, cex=1.2)

  # if multiple b-levels, label them
  all_b_levels_sub <- unique(sapply(all_cell_info, "[[", "b"))
  if (length(all_b_levels_sub)>1 && draw_axis_labels) {
    for (bn in names(b_vals)) {
      val <- b_vals[bn]
      text(par("usr")[1], +val, labels=bn,
           adj=c(1.1,0.5), xpd=NA)
      text(par("usr")[1], -val, labels=bn,
           adj=c(1.1,0.5), xpd=NA)
    }
  }

  # x-axis
  title(xlab="RT", line=0, cex.lab=1.5)
}






# ----------------------------
# 2) Single race Plot Function
# ----------------------------
############################################################
## 1) Stand-alone function for drawing noise paths
############################################################
draw_noise_paths_race <- function(x0, x1, y0, y1, s = 5, n_paths = 1, col = "black") {
  # Degenerate case: vertical line in data space (x1 == x0)
  if (x1 == x0) {
    return()
  }
  # 1) Get the slope
  slope_data <- (y1 - y0) / (x1 - x0)

  # 3) Decide if threshold is below or above y0
  # (b) Determine min, max for vertical clamping
  y_min <- min(y0, y1)
  y_max <- max(y0, y1)

  # (3) Generate noise paths
  n_steps <- 2000
  dx      <- 0.01

  for (i in seq_len(n_paths)) {
    x_vals <- numeric(n_steps + 1)
    y_vals <- numeric(n_steps + 1)

    x_vals[1] <- x0
    y_vals[1] <- y0

    for (j in seq_len(n_steps)) {
      x_vals[j + 1] <- x_vals[j] + dx

      # Brownian step in y
      new_y <- y_vals[j] + slope_data * dx + s * sqrt(dx) * rnorm(1)

      # If it goes below the lower threshold => clamp
      if (new_y < y_min) {
        new_y <- y_min
      }

      # If it goes above the upper threshold => stop early
      if (new_y > y_max) {
        # Truncate the path at this point
        x_vals <- x_vals[1:j]
        y_vals <- y_vals[1:j]
        break
      }

      # Otherwise, keep the new y
      y_vals[j + 1] <- new_y
    }

    # Draw the path
    lines(x_vals, y_vals,
          col = adjustcolor(col, alpha.f = 0.5),
          lwd = 0.3)
  }
}

draw_noise_paths_DDM <- function(x0, x1, y0, y1,
                                 s       = 2,
                                 n_paths = 1,
                                 col     = "black") {
  # Degenerate case: vertical line in data space (x1 == x0)
  if (x1 == x0) {
    return()
  }

  # 1) Compute data-space slope
  slope_data <- (y1 - y0) / (x1 - x0)

  # 2) Determine boundaries
  y_min <- -abs(y1)
  y_max <- abs(y1)

  # 3) Generate noise paths
  n_steps <- 2000
  dx      <- 0.01

  for (i in seq_len(n_paths)) {
    x_vals <- numeric(n_steps + 1)
    y_vals <- numeric(n_steps + 1)

    x_vals[1] <- x0
    y_vals[1] <- y0

    for (j in seq_len(n_steps)) {
      x_vals[j + 1] <- x_vals[j] + dx

      # Brownian step in y:
      new_y <- y_vals[j] + slope_data * dx + s * sqrt(dx) * rnorm(1)


      # Store new_y
      y_vals[j + 1] <- new_y

      # If below lower boundary => stop
      if (new_y < y_min) {
        x_vals <- x_vals[1:(j+1)]
        y_vals <- y_vals[1:(j+1)]
        break
      }
      # If above upper boundary => stop
      if (new_y > y_max) {
        x_vals <- x_vals[1:(j+1)]
        y_vals <- y_vals[1:(j+1)]
        break
      }
    }
    y_vals <- pmax(y_min, pmin(y_vals, y_max))
    # Draw the path
    lines(x_vals, y_vals,
          col = adjustcolor(col, alpha.f = 0.5),
          lwd = 0.5)
  }
}


############################################################
## 2) Modified single_race_plot function calling draw_noise_paths_race()
############################################################
single_race_plot <- function(data_sub, main_sub, plot_legend_sub,
                            x_lim, y_lim,
                            v_factor, b_factor, t0_factor,
                            v_levels, b_vals, t0_vals,
                            v_colors, t0_ltys,
                            draw_axis_labels = TRUE,
                            within_noise = TRUE,
                            squash_dens,
                            defect_scalars = NULL) {

  # 1) Identify unique cells
  all_factors_sub <- unique(c(v_factor, b_factor, t0_factor))
  all_factors_sub <- all_factors_sub[!is.na(all_factors_sub)]
  data_factors_in_data <- intersect(all_factors_sub, colnames(data_sub))

  if (length(data_factors_in_data) > 0) {
    cells <- apply(data_sub[, data_factors_in_data, drop=FALSE], 1, function(r) {
      paste(mapply(function(f,v) paste0(f,"=",v),
                   data_factors_in_data, r),
            collapse=" ")
    })
  } else {
    cells <- rep("global", nrow(data_sub))
  }
  unique_cells <- unique(cells)

  # 2) Build Cell Info (unchanged: just collecting b,v,t0, plus data density)
  alpha       <- 0.2
  cell_info <- lapply(unique_cells, function(uc) {
    cell_idx <- (cells == uc)
    tmp_data <- data_sub[cell_idx,]

    v_cell <- get_cell_param(cell_idx, v_factor,  v_levels, data_sub)
    b_cell <- get_cell_param(cell_idx, b_factor,  b_vals,   data_sub)
    t0_cell<- get_cell_param(cell_idx, t0_factor, t0_vals,  data_sub)
    med_rt <- median(tmp_data$rt, na.rm=TRUE)

    if (length(tmp_data$rt) < 2) {
      d <- list(x = c(0, tmp_data$rt), y = c(0, 0))
    } else {
      d <- density(tmp_data$rt, na.rm=TRUE)
    }
    d$y <- d$y * squash_dens * defect_scalars[[uc]] + b_cell

    list(
      cell    = uc,
      data    = tmp_data,
      v       = v_cell,
      b       = b_cell,
      t0      = t0_cell,
      density = d,
      med_rt  = med_rt
    )
  })

  # 3) Setup Plot
  base_y <- -0.1
  plot(0, 0, type="n", xlim=x_lim, ylim=y_lim,
       xlab="", ylab="", axes=FALSE, main=main_sub)
  arrows(x0=0, x1=x_lim[2], y0=base_y, y1=base_y,
         lwd=3, length=0.1)

  # 4) Sort by median RT for layering only
  med_rts_sub <- sapply(cell_info, "[[", "med_rt")
  rt_order    <- order(med_rts_sub)
  cell_info   <- cell_info[rt_order]

  min_diff_sub<- 0.05*x_lim[2]
  for (i in seq_along(cell_info)[-1]) {
    if ((cell_info[[i]]$med_rt - cell_info[[i-1]]$med_rt) < min_diff_sub) {
      cell_info[[i]]$med_rt <- cell_info[[i-1]]$med_rt + min_diff_sub
    }
  }

  # 5) Draw b-lines
  all_b_levels_sub <- unique(sapply(cell_info, "[[", "b"))

  use_b <- setNames(rep(T, length(b_vals)), names(b_vals))
  for (bn in names(b_vals)) {
    if (!level_in_subset(data_sub, b_factor, bn)){
      use_b[bn] <- F
      next
    }
    y_b <- b_vals[bn]
    segments(0, y_b, x_lim[2], y_b, lwd=1.5, lty=1, col="black")
  }
  # b_vals <- b_vals[use_b]

  # if (length(all_b_levels_sub) > 0) {
  #   text(0, min(all_b_levels_sub), "b", cex=1.5, col="black", adj=c(0,2))
  # }

  # 6) t0 lines
  t0_levels_sorted <- sort(t0_vals)
  vertical_step    <- 0.015 * (y_lim[2] - y_lim[1])
  i <- 0
  t0_y_vals <- numeric(length(t0_levels_sorted))
  for (tn in names(t0_levels_sorted)) {
    if (!level_in_subset(data_sub, t0_factor, tn)) next
    y_i  <- i*vertical_step
    t0_i <- t0_vals[tn]
    segments(0, y_i, t0_i, y_i,
             lwd=3, lty=t0_ltys[tn], col="black")
    i <- i+1
    t0_y_vals[i] <- y_i
  }
  if (i > 0) {
    text(mean(t0_vals)/3, (i - 1)*vertical_step + vertical_step/2,
         "t0", cex=1.5, col="black", adj=c(0,0))
  }

  # 7) Densities
  for (ci in cell_info) {
    t0_level_name <- names(t0_vals)[which(t0_vals == ci$t0)]
    c_lty <- if (length(t0_level_name) > 0) t0_ltys[t0_level_name] else 1

    v_level_name <- names(v_levels)[which(v_levels == ci$v)]
    if (is.null(v_factor)) {
      c_col   <- "black"
      fill_col<- "#036ffc"
    } else {
      if (length(v_level_name) == 0) {
        c_col   <- "black"
        fill_col<- "#036ffc"
      } else {
        c_col   <- v_colors[v_level_name]
        fill_col<- v_colors[v_level_name]
      }
    }

    lines(ci$density, lwd=3, col=c_col, lty=c_lty)
    polygon(ci$density, col=adjustcolor(fill_col, alpha.f=0.2), lty=c_lty)
  }

  #################################################################
  ## (h) v lines + optional noise
  ##     Using device-based slope => bigger v => steeper line visually
  #################################################################
  # The arrow starts at x0= max(t0_vals), y0= chosen_t0_offset
  # and ends at y1= (lowest threshold - 0.02)
  # We'll compute x1 so that the slope in *device* space is alpha * v
  #################################################################
  w_plot <- par("pin")[1] / diff(par("usr")[1:2])  # inches / data unit (x)
  h_plot <- par("pin")[2] / diff(par("usr")[3:4])  # inches / data unit (y)

  # scaling factor for how strongly v influences slope:
  drift_scale   <- 1  # tweak as needed to make slopes bigger or smaller

  min_b_sub <- if (length(all_b_levels_sub)) min(all_b_levels_sub) else 0
  t0_vals   <- sort(t0_vals)

  for (ci in cell_info) {
    v_here <- ci$v

    # pick color
    v_level_name <- names(v_levels)[which(v_levels == v_here)]
    if (length(v_level_name)==0) {
      c_col <- "black"
    } else {
      c_col <- v_colors[v_level_name]
    }

    # define start (x0,y0)
    x0 <- max(t0_vals)
    y0 <- mean(t0_y_vals)
    # define end y
    y1 <- min_b_sub - 0.02

    # device slope = alpha * v_here
    # slope_data = slope_device * (w_plot / h_plot)
    slope_device <- drift_scale * v_here
    slope_data   <- slope_device * (w_plot / h_plot)

    dy      <- (y1 - y0)  # negative if y1 < y0
    if (abs(slope_data) < 1e-12) {
      # avoid division by zero if v=0
      next
    }

    dx      <- dy / slope_data
    x1      <- x0 + dx

    # draw arrow
    arrows(x0 = x0, x1 = x1,
           y0 = y0, y1 = y1,
           lwd=3, length=0.1, col=c_col, lty=1)

    if (within_noise) {
      draw_noise_paths_race(
        x0 = x0,
        x1 = x1,
        y0 = y0,
        y1 = y1,
        col = c_col
      )
    }
  }

  # (i) Legend
  if (plot_legend_sub) {
    x_legend<- x_lim[2]*0.5
    y_legend<- min(b_vals)*0.95
    if (length(v_levels)>1) {
      legend(x=x_legend, y=y_legend,
             legend=names(v_levels),
             title="v",
             col=v_colors[names(v_levels)],
             lwd=3, lty=1, bty="n")
    }
    if (length(t0_vals)>1) {
      y_legend<- min(b_vals)*0.55
      legend(x=x_legend, y=y_legend,
             legend=names(t0_vals),
             title="t0",
             col="black", lwd=3,
             lty=t0_ltys[names(t0_vals)],
             bty="n")
    }
  }

  # (j) Axis
  mtext("evidence", side=2, line=.25, cex=1.2, adj = .3)
  arrows(x0=-0.025, x1=-0.025, y0=-.1, y1=max(b_vals)*1.2, lwd=3, length=0.1)

  if (length(all_b_levels_sub)>1) {
    b_factor_levels <- names(b_vals)
    for (bn in b_factor_levels) {
      val <- b_vals[bn]
      if (draw_axis_labels) {
        text(x=par("usr")[1], y=val,
             labels=bn, srt=0,
             adj=c(1.1,0.4), xpd=NA)
      }
    }
  }

  # (k) X-axis
  title(xlab="RT", line=0, cex.lab=1.5)
}




# ----------------------------
# 3) Main Plot Function
# ----------------------------
make_design_plot <- function(
    data, pars, factors, main,
    type, split = NULL,
    plot_legend = TRUE,
    plot_factor = NULL,
    within_noise = TRUE
) {
  # 1) Identify factor columns
  v_factor  <- if ("v"  %in% names(factors)) factors$v  else NULL
  if ("B"  %in% names(factors)){
    b_factor <- factors$B
  } else if("a" %in% names(factors)){
    b_factor <- factors$a
  } else {
    b_factor <- NULL
  }
  t0_factor <- if ("t0" %in% names(factors)) factors$t0 else NULL

  # rename 'a' => 'b' in 'pars'
  colnames(pars)[colnames(pars) == "a"] <- "b"

  # 2) Param means (global)
  v_vals <- get_param_means(pars, v_factor, "v")
  b_vals <- get_param_means(pars, b_factor, "b")
  t0_vals<- get_param_means(pars, t0_factor,"t0")

  # 3) If effectively global => set factor to NULL
  if (!is.null(v_factor) && length(unique(v_vals))==1) {
    v_factor <- NULL
    v_vals   <- get_param_means(pars, NULL, "v")
  }
  if (!is.null(b_factor) && length(unique(b_vals))==1) {
    b_factor <- NULL
    b_vals   <- get_param_means(pars, NULL, "b")
  }
  if (!is.null(t0_factor) && length(unique(t0_vals))==1) {
    t0_factor <- NULL
    t0_vals   <- get_param_means(pars, NULL, "t0")
  }

  # 4) Aesthetics
  v_colors <- assign_colors(v_vals)
  t0_ltys  <- assign_ltys(t0_vals)

  # 5) We do a single global call => ignoring plot_factor
  #    to get min_y_sub, max_y_sub, defect_scalars, etc.
  #    but we DO pass 'split' as 'split_factor' so it’s included in cells
  squash_dens <- length(v_vals)*.8
  if(type != "DDM"){
    squash_dens <- squash_dens * 2
  }
  global_rng <- get_subset_ranges(
    data         = data,
    b_factor     = b_factor,
    v_factor     = v_factor,
    t0_factor    = t0_factor,
    b_vals       = b_vals,
    squash_dens  = squash_dens,
    split_factor = split  # <--- important
  )

  # from global_rng:
  # global_rng$max_rt_sub
  # global_rng$min_y_sub
  # global_rng$max_y_sub
  # global_rng$defect_scalars

  # We'll define x_lims, y_lims. For DDM => symmetrical about 0
  if (type=="race") {
    x_lims <- c(0, global_rng$max_rt_sub)
    y_lims <- c(-0.15, global_rng$max_y_sub + 0.1)
  } else {
    range_val <- max(abs(global_rng$min_y_sub), abs(global_rng$max_y_sub))
    x_lims <- c(0, global_rng$max_rt_sub)
    y_lims <- c(-range_val - 0.05, range_val + 0.1)
  }

  # 6) handle plot_factor => how many subplots
  if (!is.null(plot_factor) && plot_factor %in% colnames(data)) {
    levels_plot <- unique(data[[plot_factor]])
  } else {
    levels_plot <- "ALL"
  }
  n_plots <- length(levels_plot)

  # 7) Setup multi-panel
  old_par <- par(no.readonly=TRUE)
  on.exit(par(old_par), add=TRUE)
  par(mfrow=c(1, n_plots), mar=c(5,6,3,0))

  # 8) loop over factor levels
  for (i in seq_along(levels_plot)) {
    lvl <- levels_plot[i]
    if (lvl=="ALL") {
      data_subset <- data
      main_sub    <- main
    } else {
      data_subset <- data[data[[plot_factor]]==lvl,]
      main_sub    <- paste0(main, "\n", plot_factor,"=",lvl)
    }

    # subset 'pars' similarly if needed
    if (!is.null(plot_factor) && plot_factor %in% colnames(pars) && lvl!="ALL") {
      pars_subset <- pars[pars[[plot_factor]]==lvl,]
    } else {
      pars_subset <- pars
    }

    plot_legend_sub      <- (i==n_plots)
    draw_axis_labels_sub <- (i==1)

    if (type=="race") {
      single_race_plot(
        data_sub        = data_subset,
        main_sub        = main_sub,
        plot_legend_sub = plot_legend_sub && plot_legend,
        x_lim           = x_lims,
        y_lim           = y_lims,
        v_factor        = v_factor,
        b_factor        = b_factor,
        t0_factor       = t0_factor,
        v_levels        = v_vals,
        b_vals          = b_vals,
        t0_vals         = t0_vals,
        v_colors        = v_colors,
        t0_ltys         = t0_ltys,
        draw_axis_labels= draw_axis_labels_sub,
        within_noise    = within_noise,
        squash_dens     = squash_dens,
        defect_scalars  = global_rng$defect_scalars  # pass global
      )
    } else {
      single_DDM_plot(
        data_sub        = data_subset,
        main_sub        = main_sub,
        plot_legend_sub = plot_legend_sub && plot_legend,
        x_lim           = x_lims,
        y_lim           = y_lims,
        v_factor        = v_factor,
        b_factor        = b_factor,
        t0_factor       = t0_factor,
        v_levels        = v_vals,
        b_vals          = b_vals,
        t0_vals         = t0_vals,
        v_colors        = v_colors,
        t0_ltys         = t0_ltys,
        split           = split,     # so we still do top/bottom in single_DDM_plot
        draw_axis_labels= draw_axis_labels_sub,
        within_noise = within_noise,
        squash_dens     = squash_dens,
        defect_scalars  = global_rng$defect_scalars
      )
    }
  }
}


