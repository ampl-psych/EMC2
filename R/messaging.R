# Nicely format a difftime (seconds / minutes / hours)
format_duration <- function(dt) {
  secs <- as.numeric(dt, units = "secs")
  if (secs < 60) {
    sprintf("%.1f s", secs)
  } else if (secs < 3600) {
    sprintf("%.1f min", secs / 60)
  } else {
    sprintf("%.1f h", secs / 3600)
  }
}

# Remaining time across all stages, based on the duration of the *last* try
# Stage assumptions:
#   preburn: 1-1
#   burn   : 1-max_tries
#   adapt  : 1-2
#   sample : 10-max_tries
estimate_remaining_total_time <- function(stage, tries_done, elapsed_dt, max_tries = 20L,
                                          current_iters = NULL, target_iters = NULL,
                                          step_size = NULL) {
  stage_order <- c("preburn", "burn", "adapt", "sample")

  min_total <- c(preburn = 1L, burn = 1L, adapt = 1L, sample = 10L)
  max_total <- c(preburn = 1L, burn = max_tries, adapt = 2L, sample = max_tries)

  if (!stage %in% stage_order) {
    min_rem_tries <- max(0L, 1L - tries_done)
    max_rem_tries <- max(0L, max_tries - tries_done)
    return(list(
      min_time = min_rem_tries * elapsed_dt,
      max_time = max_rem_tries * elapsed_dt
    ))
  }

  idx <- match(stage, stage_order)

  # --- Remaining in later stages (always try-count based) ---
  if (idx < length(stage_order)) {
    later_stages  <- stage_order[(idx + 1L):length(stage_order)]
    min_later_rem <- sum(min_total[later_stages])
    max_later_rem <- sum(max_total[later_stages])
  } else {
    min_later_rem <- 0L
    max_later_rem <- 0L
  }

  # --- Remaining in current stage ---
  max_current_rem <- max(0L, max_total[stage] - tries_done)

  # For sample stage: use iteration-based minimum if we have the info
  if (stage == "sample" &&
      !is.null(current_iters) && !is.null(target_iters) && !is.null(step_size)) {
    iters_remaining   <- max(0L, target_iters - current_iters)
    time_per_iter     <- elapsed_dt / step_size
    min_current_time  <- iters_remaining * time_per_iter
    min_time <- min_current_time + min_later_rem * elapsed_dt
  } else {
    min_current_rem <- max(0L, min_total[stage] - tries_done)
    min_time <- (min_current_rem + min_later_rem) * elapsed_dt
  }

  max_time <- (max_current_rem + max_later_rem) * elapsed_dt

  list(min_time = min_time, max_time = max_time)
}


accept_progress_bar <- function(min = 0, max = 1) {
  .val <- 0
  .killed <- FALSE
  .nb <- 0L
  .pc <- -1L # This ensures the initial value is displayed
  .ex <- 0
  component <- list(
    pchar = "=",
    prog_start = " |",
    prog_end = "| ",
    percent = "%3d%%",
    acc_sep = " | ",
    acc_msg = "New(%3d%%)"
  )
  width <- c(1,2,2,4,3,9) # previous code was giving warnings
  width <- split(unname(width), names(component))
  width$extras <- sum(unlist(width)) - width$pchar
  width$term <- getOption("width")
  width$progress <- trunc((width$term - width$extras) / width$pchar)

  if (max <= min) stop("must have 'max' > 'min'")

  # Handles an update to the progress bar
  up <- function(value, extra = 0) {
    if (!is.finite(value) || value < min || value > max) {
      return()
    }
    .val <<- value
    nb <- round(width$progress * (value - min) / (max - min))
    pc <- round(100 * (value - min) / (max - min))
    extra <- round(100 * extra)
    if (nb == .nb && pc == .pc && .ex == extra) {
      return()
    }
    # Clear the current progress bar
    cat(paste0("\r", strrep(" ", width$term)))
    # Write the updated progress bar
    cat(paste0(
      "\r",
      component$prog_start,
      strrep(component$pchar, nb),
      strrep(" ", width$pchar * (width$progress - nb)),
      component$prog_end,
      sprintf(component$percent, pc),
      component$acc_sep,
      sprintf(component$acc_msg, extra)))
    utils::flush.console()
    .nb <<- nb
    .pc <<- pc
    .ex <<- extra
  }

  get_value <- function() .val
  kill <- function() {
    if (!.killed) {
      cat("\n")
      utils::flush.console()
      .killed <<- TRUE
    }
  }
  up(0) # will check if in range

  structure(list(getVal = get_value, up = up, kill = kill),
            class = c("accept_progress_bar", "txtProgressBar"))
}

update_progress_bar <- function(pb, value, extra = 0) {
  if (!inherits(pb, "txtProgressBar")) {
    stop(gettextf(
      "'pb' is not from class %s",
      dQuote("txtProgressBar")
    ),
    domain = NA
    )
  }
  oldval <- pb$getVal()
  pb$up(value, extra)
  invisible(oldval)
}

accept_rate <- function(pmwgs, window_size = 200) {
  n_samples <- pmwgs$samples$idx
  if (is.null(n_samples) || n_samples < 3) {
    return(array(0, dim(pmwgs$samples$alpha)[2]))
  }
  if (n_samples <= window_size) {
    start <- 1
    end <- n_samples
  } else {
    start <- n_samples - window_size + 1
    end <- n_samples
  }
  vals <- pmwgs$samples$alpha[1, , start:end]
  if (is.null(dim(vals))) return(mean(diff(vals)!=0))
  apply(
    apply(vals, 1, diff) != 0, # If diff != 0
    2,
    mean
  )
}
