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
  width <- nchar(lapply(component, gettextf, 100), "w")
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
  if (is.null(dim(vals))) return(mean(diff(vals)==0))
  apply(
    apply(vals, 1, diff) != 0, # If diff != 0
    2,
    mean
  )
}
