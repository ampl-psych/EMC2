###############################################################################
### (1) combine_factors()
###############################################################################
combine_factors <- function(df, factor_names, new_col_name) {
  factor_names <- factor_names[factor_names %in% colnames(df)]
  if (length(factor_names) == 0) {
    df[[new_col_name]] <- "global"
    return(df)
  }
  df[[new_col_name]] <- apply(df[, factor_names, drop=FALSE], 1, function(rowvals) {
    paste(rowvals, collapse=":")
  })
  return(df)
}

###############################################################################
### (2) get_param_means()
###############################################################################
get_param_means <- function(pars, factor_name, param_name) {
  if (is.null(factor_name)) {
    return(setNames(mean(pars[[param_name]], na.rm = TRUE), "global"))
  }
  if (!factor_name %in% colnames(pars)) {
    return(setNames(mean(pars[[param_name]], na.rm = TRUE), "global"))
  }
  means <- tapply(pars[[param_name]], pars[[factor_name]], mean, na.rm = TRUE)
  means[is.na(means)] <- mean(pars[[param_name]], na.rm = TRUE)
  return(means)
}

###############################################################################
### (3) assign_colors() & assign_ltys()
###############################################################################
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

assign_ltys <- function(levels) {
  if (length(levels) <= 1) {
    return(setNames(1, names(levels)))
  }
  ltys <- c(1, 3:(length(levels) + 1))
  setNames(ltys, names(levels))
}

###############################################################################
### (4) get_cell_param()
###############################################################################
get_cell_param <- function(cell_idx, factor_name, param_vals, data_sub) {
  if (is.null(factor_name)) {
    return(param_vals["global"])
  }
  if (!factor_name %in% colnames(data_sub)) {
    return(mean(param_vals, na.rm = TRUE))
  }
  factor_levels <- unique(data_sub[[factor_name]][cell_idx])
  if (length(factor_levels) == 0) {
    return(mean(param_vals, na.rm = TRUE))
  }
  factor_level <- if (length(factor_levels) > 1) factor_levels[1] else factor_levels
  if (!factor_level %in% names(param_vals)) {
    return(mean(param_vals, na.rm = TRUE))
  } else {
    return(param_vals[factor_level])
  }
}
###############################################################################
### (5) get_defect_scalars() => single global pass
###############################################################################
###############################################################################
### (5) get_defect_scalars() => single global pass
###############################################################################
get_defect_scalars <- function(
    data,
    plot_cols_union,
    b_factor, b_vals,
    squash_dens = 1
) {
  plot_cols_union <- intersect(plot_cols_union, colnames(data))
  if (length(plot_cols_union) == 0) {
    cells <- rep("global", nrow(data))
  } else {
    cells <- apply(data[, plot_cols_union, drop=FALSE], 1, function(r) {
      paste(r, collapse=" ")
    })
  }
  unique_cells <- unique(cells)

  n_data         <- nrow(data)
  defect_scalars <- numeric(length(unique_cells))
  names(defect_scalars) <- unique_cells

  min_y_sub <- Inf
  max_y_sub <- -Inf
  for (uc in unique_cells) {
    idx <- which(cells == uc)
    d_rt<- data$rt[idx]
    if (length(d_rt)<2) {
      d <- list(x=c(0,d_rt), y=c(0,0))
    } else {
      d <- density(d_rt, na.rm=TRUE)
    }
    b_val <- get_cell_param(idx, b_factor, b_vals, data)
    cell_prop <- length(idx)/ n_data
    defect_scalars[uc] <- cell_prop

    d$y <- d$y*squash_dens*cell_prop + b_val
    local_min <- min(d$y,na.rm=TRUE)
    local_max <- max(d$y,na.rm=TRUE)
    if(local_min<min_y_sub) min_y_sub<- local_min
    if(local_max>max_y_sub) max_y_sub<- local_max
  }
  if(!is.finite(min_y_sub)) min_y_sub<-0
  if(!is.finite(max_y_sub)) max_y_sub<-0

  max_rt_sub <- quantile(data$rt, 0.95, na.rm=TRUE)
  if(!is.finite(max_rt_sub)) max_rt_sub<-3

  list(
    defect_scalars= defect_scalars,
    min_y_sub     = min_y_sub,
    max_y_sub     = max_y_sub,
    max_rt_sub    = max_rt_sub
  )
}


###############################################################################
### 7) Noise-drawing helpers
###############################################################################
draw_noise_paths_DDM <- function(x0, x1, y0, y1,
                                 s       = 3,
                                 n_paths = 5,
                                 col     = "black") {
  if (x1 == x0) {
    return()
  }
  slope_data <- (y1 - y0) / (x1 - x0)
  y_min <- -abs(y1)
  y_max <- abs(y1)
  n_steps <- 2000
  dx      <- 0.005

  for (i in seq_len(n_paths)) {
    x_vals <- numeric(n_steps + 1)
    y_vals <- numeric(n_steps + 1)

    x_vals[1] <- x0
    y_vals[1] <- y0

    for (j in seq_len(n_steps)) {
      x_vals[j + 1] <- x_vals[j] + dx
      new_y <- y_vals[j] + slope_data * dx + s * sqrt(dx) * rnorm(1)
      y_vals[j + 1] <- new_y
      if (new_y < y_min || new_y > y_max) {
        x_vals <- x_vals[1:(j+1)]
        y_vals <- y_vals[1:(j+1)]
        break
      }
    }
    y_vals <- pmax(y_min, pmin(y_vals, y_max))
    lines(x_vals, y_vals, col=adjustcolor(col, alpha.f=0.5), lwd=0.6)
  }
}

draw_noise_paths_race <- function(x0, x1, y0, y1, s = 1.5, n_paths = 5, col = "black") {
  if (x1 == x0) {
    return()
  }
  slope_data <- (y1 - y0) / (x1 - x0)
  y_min <- min(y0, y1)
  y_max <- max(y0, y1)
  n_steps <- 2000
  dx      <- 0.005

  for (i in seq_len(n_paths)) {
    x_vals <- numeric(n_steps + 1)
    y_vals <- numeric(n_steps + 1)
    x_vals[1] <- x0
    y_vals[1] <- y0

    for (j in seq_len(n_steps)) {
      x_vals[j + 1] <- x_vals[j] + dx
      new_y <- y_vals[j] + slope_data * dx + s * sqrt(dx) * rnorm(1)
      y_vals[j + 1] <- new_y
      if (new_y > y_max) {
        x_vals <- x_vals[1:(j+1)]
        y_vals <- y_vals[1:(j+1)]
        break
      }
    }
    y_vals <- pmin(y_vals, y_max)
    y_vals[y_vals < -0.1] <- -.1
    lines(x_vals, y_vals, col=adjustcolor(col, alpha.f=0.5), lwd=0.6)
  }
}



###############################################################################
### (8) single_DDM_plot()
###############################################################################
###############################################################################
### single_DDM_plot()
###############################################################################
single_DDM_plot <- function(
    data_sub,
    main_sub,
    plot_legend_sub,
    x_lim,
    y_lim,

    # factor columns for numeric param
    v_factor, b_factor, t0_factor,
    v_vals,   b_vals,   t0_vals,

    # factor columns for color
    v_factor_color = v_factor,
    t0_factor_color= t0_factor,
    v_vals_color   = v_vals,
    t0_vals_color  = t0_vals,

    # union of columns => used to build cell names
    plot_cols_union,

    v_colors,
    t0_ltys,

    split                = NULL,
    defect_scalars       = NULL,
    squash_dens          = 1,
    within_noise         = TRUE,
    draw_axis_labels     = TRUE
) {
  #### (A) Split top/bottom
  if (!is.null(split) && split %in% names(data_sub)) {
    split_levels <- sort(unique(data_sub[[split]]))
  } else {
    split_levels <- "ALL_SPLIT"
  }
  bottom_level <- split_levels[1]
  top_level    <- if (length(split_levels) > 1) split_levels[2] else bottom_level

  data_bottom <- data_sub[data_sub[[split]]==bottom_level,]
  data_top    <- data_sub[data_sub[[split]]==top_level,]

  #### (B) Build mirrored densities => same grouping as get_defect_scalars
  build_cell_info <- function(dat, sign_mult) {
    # unify => union of factor columns (including color,split, etc.)
    cols_in_data <- intersect(plot_cols_union, colnames(dat))
    if (length(cols_in_data)==0) {
      row_cells <- rep("global", nrow(dat))
    } else {
      row_cells <- apply(dat[,cols_in_data,drop=FALSE], 1, function(r) {
        paste(r, collapse=" ")
      })
    }
    unique_cells <- unique(row_cells)

    lapply(unique_cells, function(uc) {
      idx <- which(row_cells==uc)
      tmp_data <- dat[idx,]

      v_cell  <- get_cell_param(idx, v_factor,  v_vals, dat)
      b_cell  <- get_cell_param(idx, b_factor,  b_vals, dat)
      t0_cell <- get_cell_param(idx, t0_factor, t0_vals,dat)

      if (nrow(tmp_data)<2) {
        d <- list(x=c(0,tmp_data$rt), y=c(0,0))
      } else {
        d <- density(tmp_data$rt, from = x_lim[1], to = x_lim[2]*1.2)
      }
      # find proportion => defect_scalars
      cell_prop <- if (!is.null(defect_scalars[[uc]])) defect_scalars[[uc]] else 1
      d$y <- sign_mult*(d$y*squash_dens*cell_prop)+sign_mult*b_cell

      list(
        cell    = uc,
        data    = tmp_data,
        v       = v_cell,
        b       = b_cell,
        t0      = t0_cell,
        density = d,
        med_rt  = median(tmp_data$rt,na.rm=TRUE),
        sign    = sign_mult
      )
    })
  }

  info_bottom <- build_cell_info(data_bottom, -1)
  info_top    <- build_cell_info(data_top,    +1)
  all_info    <- c(info_bottom, info_top)

  #### (C) Setup plot
  plot(0,0,type="n", xlim=x_lim, ylim=y_lim, main=main_sub,
       axes=FALSE, xlab="", ylab="")
  base_line_y <- y_lim[1]-0.05
  arrows(x0=0, x1=x_lim[2], y0=base_line_y, y1=base_line_y, lwd=3, length=0.1)

  #### (D) ±b lines
  for (bn in names(b_vals)) {
    b_here <- b_vals[bn]
    segments(0,+b_here, x_lim[2], +b_here, lwd=1.5, col="black")
    segments(0,-b_here, x_lim[2], -b_here, lwd=1.5, col="black")
  }

  #### (E) Sort by median RT => layering
  med_rts <- sapply(all_info, `[[`, "med_rt")
  o       <- order(med_rts)
  all_info<- all_info[o]
  min_diff_sub <- 0.05*x_lim[2]
  for (i in seq_along(all_info)[-1]) {
    if ((all_info[[i]]$med_rt - all_info[[i-1]]$med_rt)<min_diff_sub) {
      all_info[[i]]$med_rt <- all_info[[i-1]]$med_rt+min_diff_sub
    }
  }

  #### (F) t0 lines
  t0_sorted <- sort(t0_vals)
  vertical_step <- 0.015*(y_lim[2]-y_lim[1])
  i <-0
  t0_y_vals <- numeric(length(t0_sorted))
  for (tn in names(t0_sorted)) {
    y_i <- i*vertical_step
    t0_i<- t0_vals[tn]
    lty_i<- if(tn %in% names(t0_ltys)) t0_ltys[tn] else 1
    segments(0,y_i,t0_i,y_i, lwd=3, col="black", lty=lty_i)
    t0_y_vals[i+1] <- y_i
    i<-i+1
  }

  #### (G) Densities
  alpha <-0.2
  for (ci in all_info) {
    # line type => from t0_factor_color => label from data
    idx_ci <- seq_len(nrow(ci$data))
    if (!is.null(t0_factor_color) && t0_factor_color %in% colnames(ci$data)) {
      t0_lab <- unique(ci$data[[t0_factor_color]][idx_ci])
      if (length(t0_lab)==1 && t0_lab %in% names(t0_ltys)) {
        c_lty <- t0_ltys[t0_lab]
      } else {
        c_lty <-1
      }
    } else c_lty<-1

    # color => from v_factor_color
    if (!is.null(v_factor_color) && v_factor_color %in% colnames(ci$data)) {
      v_lab <- unique(ci$data[[v_factor_color]][idx_ci])
      if (length(v_lab)==1 && v_lab %in% names(v_colors)) {
        c_col   <- v_colors[v_lab]
        fill_col<- v_colors[v_lab]
      } else {
        c_col   <- "black"
        fill_col<- "#036ffc"
      }
    } else {
      c_col   <- "black"
      fill_col<- "#036ffc"
    }

    lines(ci$density, lwd=3, col=c_col, lty=c_lty)
    polygon(ci$density, col=adjustcolor(fill_col, alpha.f=alpha), lty=c_lty)
  }

  #### (H) Drift lines => from ci$v
  df_v <- data.frame(
    idx   = seq_along(all_info),
    v_orig= sapply(all_info, function(ci) ci$v),
    stringsAsFactors = FALSE
  )

  # We'll also store the color label from data
  # so that if two cells differ in color factor, we keep them as separate lines.
  color_labels <- character(length(all_info))
  for (i in seq_along(all_info)) {
    ci     <- all_info[[i]]
    idx_ci <- seq_len(nrow(ci$data))
    # default => "global"
    color_here <- "global"
    if (!is.null(v_factor_color) && v_factor_color %in% names(ci$data)) {
      # if there's a unique color factor label in ci$data, pick it
      lab_vec <- unique(ci$data[[v_factor_color]][idx_ci])
      if (length(lab_vec)==1 && is.character(lab_vec)) {
        color_here <- lab_vec
      }
    }
    color_labels[i] <- color_here
  }
  df_v[["color_lbl"]] <- color_labels

  # absolute & sign
  df_v[["v_abs"]]  <- abs(df_v$v_orig)
  df_v[["v_sign"]] <- sign(df_v$v_orig)

  # 2) Sort by v_abs descending, then remove duplicates only if
  #    (v_abs, v_sign, color_lbl) are the same
  df_v <- df_v[order(df_v$v_abs, decreasing=TRUE), ]
  # build a "key" for uniqueness
  df_v[["uniq_key"]] <- paste0(
    round(df_v$v_abs, 8), "_sign", df_v$v_sign, "_col", df_v$color_lbl
  )
  # skip duplicates on that key
  df_v <- df_v[!duplicated(df_v$uniq_key), ]

  # 3) define the time range => min, max => quantile 0.95 maybe
  min_rt_sub <- suppressWarnings(quantile(data_sub$rt, probs=0.01, na.rm=TRUE))
  max_rt_sub <- suppressWarnings(quantile(data_sub$rt, probs=0.95, na.rm=TRUE))
  if (!is.finite(min_rt_sub)) min_rt_sub<- 0
  if (!is.finite(max_rt_sub) || max_rt_sub<=min_rt_sub) max_rt_sub<- min_rt_sub+1

  total_xrange <- (max_rt_sub - min_rt_sub)

  # We'll define x0= max(t0_vals)
  x0 <- max(t0_vals)

  # 4) ratio logic => descending order
  # first => x1= x0+ dist_02
  # next => x1[i]= x1[i-1] + dist_02*(1 - ratio)
  # ratio= v_abs[i]/ v_abs[i-1], clamp to >=0.85 unless ratio=1 => keep ratio=1
  xPositions <- numeric(nrow(df_v))
  if (nrow(df_v)>0) {
    xPositions[1] <- min_rt_sub + 0.2 * total_xrange
    if(nrow(df_v)>1) {
      for (i in seq(2, nrow(df_v))) {
        vA_cur <- df_v$v_abs[i]
        vA_prv <- df_v$v_abs[i-1]
        ratio_i<- vA_cur / vA_prv
        # If ratio_i <0.85 => ratio_i=0.85, unless ratio_i==1 => keep it
        if (abs(ratio_i -1)< 1e-9) {
          # ratio=1 => do nothing
        } else if (ratio_i > 0.95) {
          ratio_i <- 0.95
        }
        xPositions[i] <- xPositions[i-1] +  1/2*total_xrange*(1 - ratio_i)
      }
    }
  }
  df_v[["x1"]] <- xPositions

  # 5) final => draw lines in that order
  for (i in seq_len(nrow(df_v))) {
    v_sign  <- df_v$v_sign[i]
    x1      <- df_v$x1[i]
    if (!is.finite(x1)) next

    # vertical => if sign>=0 => + boundary, else => - boundary
    y0 <- mean(t0_y_vals)
    y1 <- if (v_sign>=0) min(b_vals) else -min(b_vals)

    # color => we can use df_v$color_lbl[i], or just "black"
    c_col <- v_colors[df_v$color_lbl[i]]
    # If you want each color distinct => do:
    # c_col <- if (df_v$color_lbl[i]=="global") "black" else "red"

    # arrow
    arrows(x0, y0, x1, y1, lwd=3, length=0.1, col=c_col)

    # noise function => same as before
    if (within_noise) {
      draw_noise_paths_DDM(x0, x1, y0, y1, col=c_col)
    }
  }

  #### (I) Legend
  if(plot_legend_sub){
    x_legend<- x_lim[2]*0.5
    if(length(b_vals)) {
      y_legend<- min(b_vals)*0.5
    } else y_legend<- y_lim[1]+0.2*(y_lim[2]-y_lim[1])
    # v-color legend
    if(length(v_colors)>1) {
      legend(x=x_legend, y=y_legend,
             legend=names(v_colors),
             title="v",
             col=v_colors,
             lwd=3,lty=1,bty="n")
    }
    # t0-lty legend
    if(length(t0_ltys)>1) {
      x_legend<- x_lim[2]*0.75
      legend(x=x_legend,y=y_legend,
             legend=names(t0_ltys),
             title="t0",
             col="black",lwd=3,
             lty=t0_ltys,
             bty="n")
    }
  }

  #### (J) Axes
  arrows(x0=-0.015,x1=-0.015, y0=y_lim[1],y1=y_lim[2], lwd=3, length=0.1)
  arrows(x0=-0.015,x1=-0.015, y0=y_lim[2],y1=y_lim[1], lwd=3, length=0.1)
  mtext("evidence",side=2,line=0.5,cex=1.2)

  # label b-levels if multiple
  all_b_nums <- unique(sapply(all_info, `[[`,"b"))
  if(length(all_b_nums)>1 && draw_axis_labels) {
    for(bn in names(b_vals)) {
      val <- b_vals[bn]
      text(par("usr")[1], +val, labels=bn, adj=c(1.1,0.5), xpd=NA)
      text(par("usr")[1], -val, labels=bn, adj=c(1.1,0.5), xpd=NA)
    }
  }
  title(xlab="RT",line=0, cex.lab=1.5)
}



###############################################################################
### (9) single_race_plot()
###############################################################################
single_race_plot <- function(
    data_sub,
    main_sub,
    plot_legend_sub,
    x_lim,
    y_lim,

    # "full" factor => numeric param
    v_factor, b_factor, t0_factor,
    v_vals, b_vals, t0_vals,

    # color factor => for legend
    v_factor_color = v_factor,
    t0_factor_color= t0_factor,
    v_vals_color   = v_vals,
    t0_vals_color  = t0_vals,

    # unified columns
    plot_cols_union,

    v_colors,
    t0_ltys,

    defect_scalars = NULL,
    squash_dens    = 1,
    within_noise   = TRUE,
    draw_axis_labels=TRUE
) {
  # (A) Identify cells => same columns as get_defect_scalars
  cols_in_data <- intersect(plot_cols_union, colnames(data_sub))
  if(length(cols_in_data)==0) {
    row_cells <- rep("global",nrow(data_sub))
  } else {
    row_cells <- apply(data_sub[,cols_in_data,drop=FALSE],1,function(r){
      paste(r,collapse=" ")
    })
  }
  unique_cells <- unique(row_cells)

  # (B) Build cell info
  cell_info <- lapply(unique_cells, function(uc){
    idx <- which(row_cells==uc)
    tmp_data <- data_sub[idx,]

    v_cell <- get_cell_param(idx, v_factor,v_vals,  data_sub)
    b_cell <- get_cell_param(idx, b_factor,b_vals,  data_sub)
    t0_cell<- get_cell_param(idx, t0_factor,t0_vals,data_sub)
    med_rt <- median(tmp_data$rt, na.rm=TRUE)

    if(nrow(tmp_data)<2) {
      d <- list(x=c(0,tmp_data$rt), y=c(0,0))
    } else {
      d <- density(tmp_data$rt, from = x_lim[1], to = x_lim[2]*1.2)
    }
    cell_prop <- if(!is.null(defect_scalars[[uc]])) defect_scalars[[uc]] else 1
    d$y <- d$y*squash_dens*cell_prop + b_cell

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

  # (C) Setup Plot
  base_y<- -0.1
  plot(0,0, type="n", xlim=x_lim, ylim=y_lim,
       xlab="", ylab="", axes=FALSE, main=main_sub)
  arrows(x0=0, x1=x_lim[2], y0=base_y,y1=base_y, lwd=3, length=0.1)

  # (D) Sort by median RT => layering
  med_rts_sub <- sapply(cell_info, `[[`,"med_rt")
  o           <- order(med_rts_sub)
  cell_info   <- cell_info[o]
  min_diff_sub<- 0.05*x_lim[2]
  for(i in seq_along(cell_info)[-1]) {
    if((cell_info[[i]]$med_rt - cell_info[[i-1]]$med_rt)<min_diff_sub){
      cell_info[[i]]$med_rt <- cell_info[[i-1]]$med_rt+min_diff_sub
    }
  }

  # (E) Draw b-lines
  for(bn in names(b_vals)) {
    y_b <- b_vals[bn]
    segments(0,y_b, x_lim[2],y_b, lwd=1.5, col="black")
  }

  # (F) t0 lines
  t0_sorted <- sort(t0_vals)
  vertical_step<- 0.015*(y_lim[2]-y_lim[1])
  i<-0
  t0_y_vals<- numeric(length(t0_sorted))
  for(tn in names(t0_sorted)){
    y_i<- i*vertical_step
    t0_i<- t0_vals[tn]
    lty_i<- if(tn%in%names(t0_ltys)) t0_ltys[tn] else 1
    segments(0,y_i, t0_i,y_i, lwd=3,lty=lty_i, col="black")
    t0_y_vals[i+1]<- y_i
    i<- i+1
  }

  # (G) Densities
  alpha<- 0.2
  for(ci in cell_info){
    idx_ci <- seq_len(nrow(ci$data))

    # line type => from t0_factor_color => label from data
    if(!is.null(t0_factor_color) && t0_factor_color%in%colnames(ci$data)){
      t0_lab <- unique(ci$data[[t0_factor_color]][idx_ci])
      if(length(t0_lab)==1 && t0_lab %in% names(t0_ltys)){
        c_lty<- t0_ltys[t0_lab]
      } else c_lty<-1
    } else {
      c_lty<-1
    }

    # color => from v_factor_color => label from data
    if(!is.null(v_factor_color) && v_factor_color%in%colnames(ci$data)){
      v_lab <- unique(ci$data[[v_factor_color]][idx_ci])
      if(length(v_lab)==1 && v_lab%in%names(v_colors)){
        c_col   <- v_colors[v_lab]
        fill_col<- v_colors[v_lab]
      } else {
        c_col   <- "black"
        fill_col<- "#036ffc"
      }
    } else {
      c_col   <- "black"
      fill_col<- "#036ffc"
    }

    lines(ci$density, lwd=3, col=c_col, lty=c_lty)
    polygon(ci$density, col=adjustcolor(fill_col,alpha.f=alpha), lty=c_lty)
  }

  # (H) v lines + optional noise
  #### (H) Race drift lines => skip duplicates if (abs,value,sign,color) match, ratio≥0.85, etc.
  # 1) build df of (idx, v_orig, color_lbl, v_abs, v_sign)
  min_b_sub<- if(length(b_vals)) min(b_vals) else 0

  df_v <- data.frame(
    idx   = seq_along(cell_info),
    v_orig= sapply(cell_info, function(ci) ci$v),
    stringsAsFactors=FALSE
  )
  # color
  color_vec <- character(length(cell_info))
  for (i in seq_along(cell_info)) {
    ci     <- cell_info[[i]]
    idx_ci <- seq_len(nrow(ci$data))
    lbl    <- "global"
    if (!is.null(v_factor_color) && v_factor_color %in% names(ci$data)) {
      labs <- unique(ci$data[[v_factor_color]][idx_ci])
      if (length(labs)==1 && is.character(labs)) lbl<- labs
    }
    color_vec[i] <- lbl
  }
  df_v[["color_lbl"]] <- color_vec

  df_v[["v_abs"]]  <- abs(df_v$v_orig)
  df_v[["v_sign"]] <- sign(df_v$v_orig)

  # 2) sort descending by v_abs => remove duplicates if same (v_abs, v_sign, color_lbl)
  df_v <- df_v[order(df_v$v_abs, decreasing=TRUE), ]
  df_v[["uniq_key"]] <- paste0(
    round(df_v$v_abs,8), "_sign", df_v$v_sign, "_col", df_v$color_lbl
  )
  df_v <- df_v[!duplicated(df_v$uniq_key), ]

  # 3) define time range => 95th percentile for max
  min_rt_sub <- suppressWarnings(quantile(data_sub$rt, probs = 0.01, na.rm=TRUE))
  max_rt_sub <- suppressWarnings(quantile(data_sub$rt, probs=0.95, na.rm=TRUE))
  if (!is.finite(min_rt_sub)) min_rt_sub<-0
  if (!is.finite(max_rt_sub)|| max_rt_sub<=min_rt_sub) max_rt_sub<- min_rt_sub+1

  total_xrange <- (max_rt_sub - min_rt_sub)

  x0<- max(t0_vals)
  y0<- mean(t0_y_vals)
  xPositions <- numeric(nrow(df_v))
  if (nrow(df_v)>0) {
    xPositions[1] <- min_rt_sub + .2*total_xrange
    if(nrow(df_v) >1){
      for(i in seq(2,nrow(df_v))) {
        ratio_i <- df_v$v_abs[i]/ df_v$v_abs[i-1]
        if (abs(ratio_i -1)<1e-9) {
          # ratio=1 => no clamp
        } else if (ratio_i>0.95) {
          ratio_i<- 0.95
        }
        xPositions[i] <- xPositions[i-1] + 1/2*total_xrange*(1 - ratio_i)
      }
    }
  }
  df_v[["x1"]] <- xPositions

  # 4) final => draw lines
  for (i in seq_len(nrow(df_v))) {
    x1    <- df_v$x1[i]
    # color
    c_col <- v_colors[df_v$color_lbl[i]]

    # vertical => e.g. y0= mean(t0_y_vals), y1= y0 - 0.2*(i)
    y0 <- mean(t0_y_vals)
    y1<- min_b_sub - 0.02
    arrows(x0,y0, x1,y1, lwd=3, length=0.1, col=c_col)

    if(within_noise){
      draw_noise_paths_race(x0, x1, y0, y1, col=c_col)
    }
  }

  # (I) Legend
  if(plot_legend_sub){
    x_legend<- x_lim[2]*0.6
    if(length(b_vals)) {
      y_legend<- min(b_vals)*0.95
    } else {
      y_legend<- y_lim[1]+0.2*(y_lim[2]-y_lim[1])
    }
    if(length(v_colors)>1){
      legend(x=x_legend, y=y_legend,
             legend=names(v_colors),
             title="v",
             col=v_colors,
             lwd=3,lty=1,bty="n")
    }
    if(length(t0_ltys)>1){
      y_legend<- y_legend*0.55
      legend(x=x_legend, y=y_legend,
             legend=names(t0_ltys),
             title="t0",
             col="black", lwd=3,
             lty=t0_ltys, bty="n")
    }
  }

  # (J) Axis
  mtext("evidence", side=2,line=0.25, cex=1.2, adj=.3)
  arrows(x0=-0.015,x1=-0.015, y0=-.1, y1=max(b_vals,0)*1.2, lwd=3, length=0.1)

  if(length(b_vals)>1 && draw_axis_labels){
    for(bn in names(b_vals)){
      val<- b_vals[bn]
      text(par("usr")[1], val, labels=bn, adj=c(1.1,0.4), xpd=NA)
    }
  }
  title(xlab="RT", line=0, cex.lab=1.5)
}


###############################################################################
### (10) single_LNR_plot()
###############################################################################
single_LNR_plot <- function(
    data_sub,
    main_sub,
    plot_legend_sub,
    x_lim,
    y_lim,

    # "full" factor => numeric param
    v_factor,
    t0_factor,
    v_vals,
    t0_vals,

    # color factor => for legend
    v_factor_color  = v_factor,
    t0_factor_color = t0_factor,
    v_vals_color    = v_vals,
    t0_vals_color   = t0_vals,

    # union of columns => used to build cell names
    plot_cols_union,

    v_colors,
    t0_ltys,

    defect_scalars   = NULL,
    squash_dens      = 1,
    within_noise     = TRUE,     # Not actually used for LNR
    draw_axis_labels = TRUE
) {
  ##################################################################
  ### (A) Figure out anchoring & arrow positions
  ##################################################################
  # 1) Horizontal arrow near bottom:
  arrow_y <- 0  # We will draw an arrow at y=0
  # 2) We want to shift densities up by 0.1.
  anchor_y <- 0.1
  # 3) The user asked to anchor x at mean(t0_vals):
  mean_t0  <- mean(t0_vals, na.rm=TRUE)
  # We will shift each density so its minimum x lands at mean_t0.

  ##################################################################
  ### (B) Identify unique "cells" by combining the relevant columns
  ##################################################################
  cols_in_data <- intersect(plot_cols_union, colnames(data_sub))
  if (length(cols_in_data) == 0) {
    row_cells <- rep("global", nrow(data_sub))
  } else {
    row_cells <- apply(data_sub[, cols_in_data, drop=FALSE], 1, function(r) {
      paste(r, collapse=" ")
    })
  }
  unique_cells <- unique(row_cells)

  ##################################################################
  ### (C) Build per-cell info (density, param values, median RT, etc.)
  ##################################################################
  cell_info <- lapply(unique_cells, function(uc) {
    idx <- which(row_cells == uc)
    tmp_data <- data_sub[idx, , drop=FALSE]

    v_cell  <- get_cell_param(idx, v_factor,  v_vals,  data_sub)
    t0_cell <- get_cell_param(idx, t0_factor, t0_vals, data_sub)
    med_rt  <- median(tmp_data$rt, na.rm=TRUE)

    # Density
    if (nrow(tmp_data) < 2) {
      d <- list(x = c(0, tmp_data$rt), y = c(0, 0))
    } else {
      d <- density(tmp_data$rt, from = x_lim[1], to = x_lim[2]*1.2)
    }

    # Scale the density by defect_scalars
    cell_prop <- if (!is.null(defect_scalars[[uc]])) defect_scalars[[uc]] else 1
    d$y <- d$y * squash_dens * cell_prop

    # We do NOT shift here; we store that for the next step
    list(
      cell    = uc,
      data    = tmp_data,
      v       = v_cell,
      t0      = t0_cell,
      density = d,
      med_rt  = med_rt
    )
  })

  ##################################################################
  ### (D) Set up the plot
  ##################################################################
  plot(0, 0, type="n", xlim=x_lim, ylim=y_lim,
       xlab="", ylab="", axes=FALSE, main=main_sub)

  # Draw the horizontal arrow at arrow_y
  arrows(x0 = 0, x1 = x_lim[2], y0 = arrow_y, y1 = arrow_y, lwd=3, length=0.1)
  title(xlab="RT",line=0, cex.lab=1.5)

  ##################################################################
  ### (E) Sort cells by median RT => layering of density polygons
  ##################################################################
  med_rts_sub <- sapply(cell_info, `[[`, "med_rt")
  o           <- order(med_rts_sub)
  cell_info   <- cell_info[o]

  # Optional: ensure minimal spacing in layering
  min_diff_sub <- 0.05 * x_lim[2]
  for (i in seq_along(cell_info)[-1]) {
    if ((cell_info[[i]]$med_rt - cell_info[[i-1]]$med_rt) < min_diff_sub) {
      cell_info[[i]]$med_rt <- cell_info[[i-1]]$med_rt + min_diff_sub
    }
  }

  ##################################################################
  ### (F) t0 lines
  ### We'll plot them offset in x by mean_t0, and in y by anchor_y
  ##################################################################
  if (length(t0_vals) > 0) {
    t0_sorted <- sort(t0_vals)
    vertical_step <- 0.1
    i <- 0
    for (tn in names(t0_sorted)) {
      y_i  <- anchor_y + i*vertical_step
      # We used to do segments(0, y_i, t0_i, y_i)
      # Now let's anchor them at x= mean_t0:
      t0_i <- t0_sorted[tn]
      lty_i <- if (tn %in% names(t0_ltys)) t0_ltys[tn] else 1

      # If you'd rather "start" the line at mean_t0 and extend by t0_i
      # that might be: segments(mean_t0, y_i, mean_t0 + t0_i, y_i)
      # but that can push them far out. Another approach is to do:
      # segments(mean_t0 - t0_i, y_i, mean_t0, y_i), etc.
      # For simplicity, let's do from (mean_t0) to (mean_t0 + t0_i):
      segments(0, y_i, 0 + t0_i, y_i, lwd=3, lty=lty_i, col="black")

      i <- i + 1
    }
  }

  ##################################################################
  ### (G) Draw the densities themselves
  ### Shift them so their minimum x is at mean_t0, and shift y by anchor_y.
  ##################################################################
  alpha <- 0.2
  for (ci in cell_info) {
    d <- ci$density
    # SHIFT in x so min(d$x) => mean_t0
    d$x <- d$x + mean_t0 + .001

    # SHIFT in y by anchor_y
    d$y <- d$y + anchor_y + (length(t0_vals) -1)*0.1*.5

    # Determine color from v_factor_color
    idx_ci <- seq_len(nrow(ci$data))
    c_col   <- "black"
    fill_col<- "#036ffc"
    if (!is.null(v_factor_color) && v_factor_color %in% names(ci$data)) {
      v_lab <- unique(ci$data[[v_factor_color]][idx_ci])
      if (length(v_lab) == 1 && v_lab %in% names(v_colors)) {
        c_col   <- v_colors[v_lab]
        fill_col<- v_colors[v_lab]
      }
    }

    # Determine lty from t0_factor_color
    c_lty <- 1
    if (!is.null(t0_factor_color) && t0_factor_color %in% names(ci$data)) {
      t0_lab <- unique(ci$data[[t0_factor_color]][idx_ci])
      if (length(t0_lab) == 1 && t0_lab %in% names(t0_ltys)) {
        c_lty <- t0_ltys[t0_lab]
      }
    }

    # Draw the line + polygon
    lines(d, lwd=3, col=c_col, lty=c_lty)
    polygon(d, col=adjustcolor(fill_col, alpha.f=alpha), lty=c_lty)
  }

  ##################################################################
  ### (H) Legend
  ##################################################################
  if (plot_legend_sub) {
    x_legend <- x_lim[2] * 0.6
    y_legend <- y_lim[2] * .8

    # v-color legend
    if (length(v_colors) > 1) {
      legend(
        x      = x_legend,
        y      = y_legend,
        legend = names(v_colors),
        title  = "v",
        col    = v_colors,
        lwd    = 3,
        lty    = 1,
        bty    = "n"
      )
    }

    # t0-lty legend
    if (length(t0_ltys) > 1) {
      y_legend <- y_lim[2] * 0.6
      legend(
        x      = x_legend,
        y      = y_legend,
        legend = names(t0_ltys),
        title  = "t0",
        col    = "black",
        lwd    = 3,
        lty    = t0_ltys,
        bty    = "n"
      )
    }
  }
}



###############################################################################
### (11) make_design_plot()
###############################################################################
make_design_plot <- function(
    data,
    pars,
    factors,
    main,
    type,
    layout = NA,
    split        = NULL,
    plot_legend  = TRUE,
    plot_factor  = NULL,
    within_noise = TRUE
) {
  #### (1) Rename "a" => "b"
  if ("a" %in% colnames(pars)) {
    colnames(pars)[colnames(pars) == "a"] <- "b"
  }
  if ("m" %in% colnames(pars)) {
    colnames(pars)[colnames(pars) == "m"] <- "v"
  }
  if(type == "LNR"){
    pars$b <- 0
  }

  #### (2) Identify factor sets
  if ("v" %in% names(factors)) {
    v_factors <- factors$v
  } else if ("m" %in% names(factors)) {
    v_factors <- factors$m
  } else {
    v_factors <- NULL
  }
  if ("B" %in% names(factors)) {
    b_factors <- factors$B
  } else if ("a" %in% names(factors)) {
    b_factors <- factors$a
  } else {
    b_factors <- NULL
  }
  t0_factors <- if ("t0" %in% names(factors)) factors$t0 else NULL

  #### A helper to build param columns => "x_factor_full","x_factor_color"
  build_param_columns <- function(df, param_factors, param_name, ignore_factor) {
    param_full_col  <- paste0(param_name, "_factor_full")
    param_color_col <- paste0(param_name, "_factor_color")

    df <- combine_factors(df, param_factors, param_full_col)

    # color => remove plot_factor
    if (!is.null(ignore_factor)) {
      param_factors_for_color <- setdiff(param_factors, ignore_factor)
    } else {
      param_factors_for_color <- param_factors
    }
    if (length(param_factors_for_color) == 0) {
      df[[param_color_col]] <- "global"
    } else {
      df <- combine_factors(df, param_factors_for_color, param_color_col)
    }
    df
  }

  #### (3) Build combined columns in both pars & data
  pars <- build_param_columns(pars, v_factors,  "v",  plot_factor)
  pars <- build_param_columns(pars, b_factors,  "b",  plot_factor)
  pars <- build_param_columns(pars, t0_factors, "t0", plot_factor)

  data <- build_param_columns(data, v_factors,  "v",  plot_factor)
  data <- build_param_columns(data, b_factors,  "b",  plot_factor)
  data <- build_param_columns(data, t0_factors, "t0", plot_factor)

  # full columns
  v_factor_full  <- if (!is.null(v_factors))  "v_factor_full"  else NULL
  b_factor_full  <- if (!is.null(b_factors))  "b_factor_full"  else NULL
  t0_factor_full <- if (!is.null(t0_factors)) "t0_factor_full" else NULL

  # color columns
  v_factor_color  <- if (!is.null(v_factors))  "v_factor_color"  else NULL
  b_factor_color  <- if (!is.null(b_factors))  "b_factor_color"  else NULL
  t0_factor_color <- if (!is.null(t0_factors)) "t0_factor_color" else NULL

  #### (4) A union of columns => for consistent cell naming
  plot_cols_union <- unique(c(
    v_factor_full,   v_factor_color,
    b_factor_full,   b_factor_color,
    t0_factor_full,  t0_factor_color,
    split
  ))
  if(!is.null(plot_factor) && plot_factor %in% colnames(data)){
    plot_cols_union <- unique(c(plot_cols_union, plot_factor))
  }
  plot_cols_union <- plot_cols_union[!is.na(plot_cols_union)]

  #### (5) For global axes => do a single get_defect_scalars
  # We need b param means for y-scaling
  b_vals_full <- if(!is.null(b_factor_full)){
    get_param_means(pars,b_factor_full,"b")
  } else {
    setNames(mean(pars[["b"]],na.rm=TRUE),"global")
  }
  # compute
  global_out <- get_defect_scalars(
    data            = data,
    plot_cols_union = plot_cols_union,
    b_factor        = b_factor_full,
    b_vals          = b_vals_full,
    squash_dens     = length(b_vals_full)*0.8
  )
  # from global_out => defect_scalars, min_y_sub, max_y_sub, max_rt_sub
  x_lims <- c(0, global_out$max_rt_sub)
  # Determine an ideal ratio of evidence plot width to density width
  # Twice as much plot_width
  if(type != "LNR"){
    target_ratio <- .4
    current_ratio <- (global_out$max_y_sub - max(b_vals_full))/max(b_vals_full)
    # Multiply defect scalars by ratio
    global_out$defect_scalars <- global_out$defect_scalars*(target_ratio/current_ratio)
    global_out$max_y_sub <- max(b_vals_full) + (global_out$max_y_sub - max(b_vals_full))*(target_ratio/current_ratio)
  } else{
    multiplier <- 1.5/max(global_out$defect_scalars)
    global_out$defect_scalars <- global_out$defect_scalars*multiplier
    global_out$max_y_sub <- global_out$max_y_sub*multiplier
  }

  if(type=="race"){
    y_lims <- c(-0.15, global_out$max_y_sub+0.1)
  } else if(type == "LNR"){
    y_lims <- c(-0.15, global_out$max_y_sub*1.2)
  } else {
    range_val<- max(abs(global_out$min_y_sub), abs(global_out$max_y_sub))
    y_lims   <- c(-range_val-0.05, range_val+0.1)
  }

  #### (6) Possibly multiple subplots => levels_plot
  if(!is.null(plot_factor) && plot_factor%in% colnames(data)){
    levels_plot <- unique(data[[plot_factor]])
  } else {
    levels_plot <- "ALL"
  }
  n_plots <- length(levels_plot)

  if (any(is.na(layout))) {
    par(mfrow = coda_setmfrow(Nchains = 1, Nparms = n_plots, nplots = 1))
  } else if(!is.null(layout)){
    par(mfrow = layout)
  }
  par(mar = c(4, 4, 2, 1))

  #### (7) for each level => subset data, but do param wise subsetting for v,b,t0
  for (i in seq_along(levels_plot)) {
    lvl <- levels_plot[i]
    if(lvl=="ALL"){
      data_subset <- data
      main_sub    <- main
    } else {
      data_subset <- data[data[[plot_factor]]==lvl, ]
      main_sub    <- paste0(main,"\n",plot_factor,"=",lvl)
    }

    # v => if (plot_factor %in% v_factors) subset => else keep all
    if(!is.null(v_factors) && plot_factor%in%v_factors && lvl!="ALL"){
      v_pars_sub <- pars[pars[[plot_factor]]==lvl,]
    } else {
      v_pars_sub <- pars
    }
    # b => if (plot_factor %in% b_factors) subset => else keep all
    if(!is.null(b_factors) && plot_factor%in%b_factors && lvl!="ALL"){
      b_pars_sub <- pars[pars[[plot_factor]]==lvl,]
    } else {
      b_pars_sub <- pars
    }
    # t0 => likewise
    if(!is.null(t0_factors) && plot_factor%in%t0_factors && lvl!="ALL"){
      t0_pars_sub <- pars[pars[[plot_factor]]==lvl,]
    } else {
      t0_pars_sub <- pars
    }

    # Now compute param means for each param subset
    v_vals_sub <- if(!is.null(v_factor_full)){
      get_param_means(v_pars_sub, v_factor_full, "v")
    } else setNames(mean(v_pars_sub[["v"]],na.rm=TRUE),"global")

    b_vals_sub <- if(!is.null(b_factor_full)){
      get_param_means(b_pars_sub, b_factor_full, "b")
    } else setNames(mean(b_pars_sub[["b"]],na.rm=TRUE),"global")

    t0_vals_sub<- if(!is.null(t0_factor_full)){
      get_param_means(t0_pars_sub,t0_factor_full,"t0")
    } else setNames(mean(t0_pars_sub[["t0"]],na.rm=TRUE),"global")

    # color => do the same logic if you want numeric
    v_vals_color_sub <- if(!is.null(v_factor_color)){
      if(plot_factor%in%v_factors && lvl!="ALL"){
        # v depends on plot_factor => use v_pars_sub
        get_param_means(v_pars_sub,v_factor_color,"v")
      } else {
        # no subsetting
        get_param_means(pars, v_factor_color,"v")
      }
    } else c(global=mean(pars[["v"]]))

    t0_vals_color_sub<- if(!is.null(t0_factor_color)){
      get_param_means(t0_pars_sub,t0_factor_color,"t0")
    } else c(global=mean(pars[["t0"]]))

    # build color sets
    v_colors_map <- assign_colors(v_vals_color_sub)
    t0_ltys_map  <- assign_ltys(t0_vals_color_sub)

    plot_legend_sub      <- (i==n_plots)
    draw_axis_labels_sub <- (i==1)

    if(type=="race"){
      single_race_plot(
        data_sub        = data_subset,
        main_sub        = main_sub,
        plot_legend_sub = plot_legend_sub && plot_legend,
        x_lim           = x_lims,
        y_lim           = y_lims,

        v_factor  = v_factor_full,
        b_factor  = b_factor_full,
        t0_factor = t0_factor_full,

        v_vals   = v_vals_sub,
        b_vals   = b_vals_sub,
        t0_vals  = t0_vals_sub,

        v_factor_color = v_factor_color,
        t0_factor_color= t0_factor_color,
        v_vals_color   = v_vals_color_sub,
        t0_vals_color  = t0_vals_color_sub,

        plot_cols_union = plot_cols_union,

        v_colors  = v_colors_map,
        t0_ltys   = t0_ltys_map,

        defect_scalars  = global_out$defect_scalars,
        squash_dens     = length(b_vals_full)*0.8,
        within_noise    = within_noise,
        draw_axis_labels= draw_axis_labels_sub
      )
    } else if(type == "DDM"){
      single_DDM_plot(
        data_sub        = data_subset,
        main_sub        = main_sub,
        plot_legend_sub = plot_legend_sub && plot_legend,
        x_lim           = x_lims,
        y_lim           = y_lims,

        v_factor  = v_factor_full,
        b_factor  = b_factor_full,
        t0_factor = t0_factor_full,

        v_vals   = v_vals_sub,
        b_vals   = b_vals_sub,
        t0_vals  = t0_vals_sub,

        v_factor_color = v_factor_color,
        t0_factor_color= t0_factor_color,
        v_vals_color   = v_vals_color_sub,
        t0_vals_color  = t0_vals_color_sub,

        plot_cols_union = plot_cols_union,

        v_colors  = v_colors_map,
        t0_ltys   = t0_ltys_map,
        split     = split,

        defect_scalars  = global_out$defect_scalars,
        squash_dens     = length(b_vals_full)*0.8,
        within_noise    = within_noise,
        draw_axis_labels= draw_axis_labels_sub
      )
    } else{
      # Then call single_LNR_plot(...) with no b.
      single_LNR_plot(
        data_sub         = data,
        main_sub         = main,
        plot_legend_sub  = plot_legend,
        x_lim            = x_lims,
        y_lim            = y_lims,
        v_factor         = v_factor_full,
        t0_factor        = t0_factor_full,
        v_vals           = v_vals_sub,
        t0_vals          = t0_vals_sub,
        v_factor_color   = v_factor_color,
        t0_factor_color  = t0_factor_color,
        v_vals_color     = v_vals_color_sub,
        t0_vals_color    = t0_vals_color_sub,
        plot_cols_union  = plot_cols_union,
        v_colors         = v_colors_map,
        t0_ltys          = t0_ltys_map,
        defect_scalars   = global_out$defect_scalars,
        squash_dens      = length(b_vals_sub)*0.8,  # or a fixed factor
        within_noise     = within_noise,
        draw_axis_labels = draw_axis_labels_sub
      )
    }
  }
}
