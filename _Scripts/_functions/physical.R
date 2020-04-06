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


# Matches lengths of stimuli and responses via downsampling

match_lengths <- function(x1, y1, x2, y2) {

  # Get longest and shortest path lengths for stim and tracing
  pathlens <- c(max(which(!is.na(x1))), max(which(!is.na(x2))))
  longest <- max(pathlens)
  shortest <- min(pathlens)

  # Get indices of frames to keep from shorter of two paths
  keep <- round(seq(from = 1, to = longest, by = longest / shortest))

  # Drop frames from longer path to make stim and tracing equal in frame count
  if (pathlens[1] == longest) {
    matched_df <- tibble(
      x = x1[keep], y = y1[keep],
      trace.x = x2[1:shortest], trace.y = y2[1:shortest]
    )
  } else {
    matched_df <- tibble(
      x = x1[1:shortest], y = y1[1:shortest],
      trace.x = x2[keep], trace.y = y2[keep]
      )
  }

  matched_df
}



### Transformation functions ###

# NOTE: rename to "do_procrustes" and "do_dtw" for better readability?

procrustes2df <- function(x1, y1, x2, y2) {

  stim <- cbind(x1, y1)
  resp <- cbind(x2, y2)

  # Do procrustes transformation and get translated dataframe
  proc <- procrustes(stim, resp, scale = TRUE, symmetric = FALSE)
  proc_mat <- predict(proc, resp)
  proc_df <- as_tibble(proc_mat, .name_repair = ~ c("proc.x", "proc.y"))

  # Bind stim columns to procrustes data
  proc_df <- add_column(proc_df, x = x1, .before = 1)
  proc_df <- add_column(proc_df, y = y1, .before = 2)

  # Get procrustes metrics and add to dataframe
  proc_dx <- proc$xmean[1] - mean(x2)
  proc_dy <- proc$xmean[2] - mean(y2)
  proc_df$translation <- sqrt(proc_dx ** 2 + proc_dy ** 2)
  proc_df$scale <- proc$scale
  proc_df$rotation <- acos(proc$rotation[1, 1])

  proc_df
}


dtw2df <- function(x1, y1, x2, y2) {

  # Do dynamic time warping and get warped dataframe
  warped <- dtw(
    x = cbind(x2, y2), y = cbind(x1, y1),
    dist.method = "Euclidean",
    step.pattern = symmetric1,
    window.type = "none",
    keep.internals = TRUE,
    distance.only = FALSE,
    open.end = FALSE,
    open.begin = FALSE
  )

  # Return dataframe with warped stimulus and response points
  dtw_df <- tibble(
    x_w = x1[warped$index2],
    y_w = y1[warped$index2],
    trace.x_w = x2[warped$index1],
    trace.y_w = y2[warped$index1]
  )

  dtw_df
}
