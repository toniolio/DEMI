### Physical tracing preprocessing and analysis functions for TraceLab ###


### Import required packages ###

library(tibble)
library(vegan)
library(dtw)



### Utility functions ###

# Calculates standard deviation of paired differences

paired.sd <- function(x) {
  SS <- sum(x ** 2)
  df <- length(x) - 1
  sqrt(SS / df)
}


# Reinterpolates a figure or tracing, either by time or by distance

reinterpolate <- function(x, y, time, equidistant = FALSE, fps = 60, n = NA) {

  # If equidistant, interpolate constant-velocity points along path. Otherwise,
  # use timestamps to generate points at same variable speed as input tracing.
  if (equidistant) {
    dists <- sqrt((x - lag(x))^2 + (y - lag(y))^2)
    dists[1] <- 0
    steps <- cumsum(dists)
  } else {
    steps <- time
  }

  # Figure out number of frames to generate based on tracing time and framerate,
  # unless number of frames is explicitly provided
  if (!is.na(n)) {
    out_num <- n
  } else {
    out_num <- round(max(time) / (1 / fps)) + 1 # add 1 because first time is 0
  }

  # Return new interpolated points
  interpolated <- tibble(
    x = approx(steps, x, n = out_num, method = "linear")$y,
    y = approx(steps, y, n = out_num, method = "linear")$y
  )

  interpolated
}


# Matches length of stimulus to response via reinterpolation

match_lengths <- function(x, y, tx, ty, equidist = FALSE) {

  # Get path lengths of stim and tracing
  stim_n <- max(which(!is.na(x)))
  trace_n <- max(which(!is.na(tx)))

  # Reinterpolate stim to make it equal to tracing in length
  stime <- seq(1, stim_n) * (1 / 60)
  newstim <- reinterpolate(x[!is.na(x)], y[!is.na(y)], stime, n = trace_n)
  matched_df <- tibble(
    x = newstim$x, y = newstim$y,
    trace.x = tx[!is.na(tx)], trace.y = ty[!is.na(ty)]
  )

  # If enabled, reinterpolate tracing so all points are equidistant
  if (equidist) {
    ttime <- seq(1, trace_n) * (1 / 60)
    newtrace <- reinterpolate(
      tx[!is.na(tx)], ty[!is.na(ty)],
      ttime, equidistant = TRUE, n = trace_n
    )
    matched_df$trace.x <- newtrace$x
    matched_df$trace.y <- newtrace$y
  }

  matched_df
}



### Transformation functions ###

procrustes2df <- function(x, y, tx, ty) {

  stim <- cbind(x, y)
  resp <- cbind(tx, ty)

  # Do procrustes transformation and get translated dataframe
  proc <- procrustes(stim, resp, scale = TRUE, symmetric = FALSE)
  proc_mat <- predict(proc, resp)
  proc_df <- as_tibble(proc_mat, .name_repair = ~ c("proc.x", "proc.y"))

  # Bind stim columns to procrustes data
  proc_df <- add_column(proc_df, x = x, .before = 1)
  proc_df <- add_column(proc_df, y = y, .before = 2)

  # Get procrustes metrics and add to dataframe
  proc_dx <- proc$xmean[1] - mean(tx)
  proc_dy <- proc$xmean[2] - mean(ty)
  proc_df$translation <- sqrt(proc_dx ** 2 + proc_dy ** 2)
  proc_df$scale <- proc$scale
  proc_df$rotation <- acos(proc$rotation[1, 1])

  proc_df
}


dtw2df <- function(x, y, tx, ty, step = asymmetric) {

  # Do dynamic time warping and get warped dataframe
  warped <- dtw(
    x = cbind(tx[!is.na(tx)], ty[!is.na(ty)]),
    y = cbind(x[!is.na(x)], y[!is.na(y)]),
    dist.method = "Euclidean",
    step.pattern = step,
    window.type = "none",
    keep.internals = TRUE,
    distance.only = FALSE,
    open.end = FALSE,
    open.begin = FALSE
  )

  # Return dataframe with warped stimulus and response points
  dtw_df <- tibble(
    x_w = x[warped$index2],
    y_w = y[warped$index2],
    trace.x_w = tx[warped$index1],
    trace.y_w = ty[warped$index1]
  )

  dtw_df
}
