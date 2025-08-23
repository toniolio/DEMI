### Code for plotting and visualizing TraceLab figures ###


### Import required packages ###

library(tidyr)
library(ggplot2)



### Load per-project settings ###

source("./_Scripts/_settings.R")



### Utility Functions ###

# Extracts the data for a given trial from a dataframe

get_trial <- function(dat, p_num, s_num, b_num, t_num) {
  subset(dat, id == p_num & session == s_num & block == b_num & trial == t_num)
}



### Functions for plotting single trials ###

# Plots a stimulus or tracing as a series of points, w/ optional colour coding

plot_figure_pts <- function(px, py, col = FALSE, path = FALSE) {
  hide_legend <- length(col) == 1
  dat <- data.frame(x = px, y = py, color = col)
  plt <- ggplot(dat, aes(x = x, y = y))
  if (path & nrow(dat) > 1) {
    plt <- plt + geom_path(color = "grey80")
  }
  plt +
    geom_point(alpha = 0.5, aes(color = color)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, screen_res[1])) +
    scale_y_reverse(expand = c(0, 0), limits = c(screen_res[2], 0)) +
    coord_fixed() +
    theme(legend.position = ifelse(hide_legend, "none", "right"))
}


# Plots a stimulus or tracing as a continuous path of line segments

plot_figure_path <- function(px, py) {
  dat <- data.frame(x = px, y = py)
  ggplot(dat, aes(x = x, y = y)) +
    geom_path() +
    scale_x_continuous(expand = c(0, 0), limits = c(0, screen_res[1])) +
    scale_y_reverse(expand = c(0, 0), limits = c(screen_res[2], 0)) +
    coord_fixed()
}


# Plots the stimulus and tracing points from a given trial together in
# different colours for comparison of overall shape

plot_tracing_pts <- function(px, py, tx, ty, path = FALSE) {

  max_len <- max(length(px), length(tx))
  length(px) <- length(py) <- length(tx) <- length(ty) <- max_len
  dat <- data.frame(
    n = seq_len(max_len),
    fig.x = px, fig.y = py,
    trace.x = tx, trace.y = ty
  )
  dat <- as_tibble(dat) %>%
    gather(2:5, key = "name", value = "val") %>%
    separate(name, c("type", "coord"), sep = "\\.") %>%
    spread(coord, val)

  plt <- ggplot(dat, aes(x = x, y = y, color = type))
  if (path) {
    plt <- plt + geom_path(alpha = 0.5)
  }
  plt +
    geom_point(alpha = 0.5) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, screen_res[1])) +
    scale_y_reverse(expand = c(0, 0), limits = c(screen_res[2], 0)) +
    coord_fixed()
}


# Visualizes the error in a tracing using translucent lines to connect tracing
# points to their corresponding stimulus frames. Tracing and stimulus must have
# an equal number of frames for this to work properly.

plot_tracing_err <- function(px, py, tx, ty) {
  dat <- data.frame(
    n = seq_along(px),
    fig.x = px, fig.y = py,
    trace.x = tx, trace.y = ty
  )
  dat2 <- as_tibble(dat) %>%
    gather(2:5, key = "name", value = "val") %>%
    separate(name, c("type", "coord"), sep = "\\.") %>%
    spread(coord, val)
  ggplot() +
    geom_segment(
      data = dat, alpha = 0.5,
      aes(x = fig.x, y = fig.y, xend = trace.x, yend = trace.y)
    ) +
    geom_point(data = dat2, aes(x = x, y = y, color = type), alpha = 0.5) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, screen_res[1])) +
    scale_y_reverse(expand = c(0, 0), limits = c(screen_res[2], 0)) +
    coord_fixed()
}



### Functions for plotting multiple trials ###

# Renders PDFs of all trials in a summary data frame to a given output path,
# colour coding the sample points using the 'color_code' argument

plot_trials <- function(trials, samples, color_code = NULL, outdir = ".") {

  # Create output dir, deleting first if it already exists
  if (outdir != ".") {
    if (dir.exists(outdir)) {
      unlink(outdir, recursive = TRUE)
    }
    dir.create(outdir, recursive = TRUE)
  }

  # Render and save a PDF plot for each trial in trials
  for (i in 1:nrow(trials)) {
    vals <- trials[i, ]
    fname <- paste0(paste(
      paste0("p", vals$id), paste0("s", vals$session),
      paste0("b", vals$block), paste0("t", vals$trial),
      sep = "_"
    ), ".pdf")
    tmp <- get_trial(samples, vals$id, vals$session, vals$block, vals$trial)
    if (is.null(color_code)) {
      plt <- plot_figure_pts(tmp$x, tmp$y, path = TRUE)
    } else {
      tmp$col <- with(tmp, eval(parse(text = color_code)))
      plt <- plot_figure_pts(tmp$x, tmp$y, tmp$col, path = TRUE)
    }
    ggsave(fname, plt, path = outdir, width = 11, height = 6)
  }
}


# Renders PDFs of all trials in a summary data frame to a given output path,
# overlaying the trial's tracing points/path on top of its figure animation
# frames/path.

plot_trial_paths <- function(trials, frames, samples, outdir = ".") {

  # Create output dir, deleting first if it already exists
  if (outdir != ".") {
    if (dir.exists(outdir)) {
      unlink(outdir, recursive = TRUE)
    }
    dir.create(outdir, recursive = TRUE)
  }

  # Render and save a PDF plot for each trial in trials
  for (i in 1:nrow(trials)) {
    vals <- trials[i, ]
    fname <- paste0(paste(
      paste0("p", vals$id), paste0("s", vals$session),
      paste0("b", vals$block), paste0("t", vals$trial),
      sep = "_"
    ), ".pdf")
    tf <- get_trial(frames, vals$id, vals$session, vals$block, vals$trial)
    tt <- get_trial(samples, vals$id, vals$session, vals$block, vals$trial)
    suppressWarnings({
      plt <- plot_tracing_pts(tf$x, tf$y, tt$x, tt$y, path = TRUE)
      ggsave(fname, plt, path = outdir, width = 11, height = 6)
    })
  }
}
