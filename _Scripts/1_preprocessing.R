##########################################
### Tracelab data preprocessing script ###
##########################################


### Import required packages and functions ###

library(dplyr)
library(TSEntropies)

source("./_Scripts/_functions/complexity.R")
source("./_Scripts/_functions/physical.R")



### Import figure data ###

source("./_Scripts/0_import.R")



### Analyze figure stimuli ###

# First, get bezier path length by summing distances between frames

pathlen_summary <- frames %>%
  group_by(id, session, block, trial) %>%
  mutate(seglen = line_length(lag(x), lag(y), x, y)) %>%
  summarize(PLstim = sum(seglen, na.rm = TRUE))


# Calculate real bezier path length, sinuosity, and total absolute curvature

# NOTE: get bezier bounds for height/width and proximity to screen edges?

figsummary <- segments %>%
  mutate(
    curve_len = bezier_length(start.x, start.y, end.x, end.y, ctrl.x, ctrl.y),
    line_len = line_length(start.x, start.y, end.x, end.y)
  ) %>%
  group_by(id, session, block, trial) %>%
  summarize(
    real_length = sum(curve_len),
    total_point_dist = sum(line_len),
    totabscurv = total_abs_curvature(
      start.x, start.y,
      end.x, end.y,
      ctrl.x, ctrl.y
    )
  ) %>%
  mutate(sinuosity = real_length / total_point_dist) %>%
  select(c(1:trial, real_length, sinuosity, totabscurv))


# Calcuate turning angle for entropy metrics

points_per_segment <- 60
tvals <- seq(0, points_per_segment - 1) / points_per_segment

turnangledat <- segments %>%
  group_by(id, session, block, trial) %>%
  group_modify(
    ~ get_fig_points(tvals,
      .x$start.x, .x$start.y,
      .x$end.x, .x$end.y,
      .x$ctrl.x, .x$ctrl.y
    )
  ) %>%
  mutate(
    theta = get_angle_diffs(x - lag(x), y - lag(y))
  )


# Calculate entropy metrics from turning angle data

entropy_summary <- turnangledat %>%
  filter(!is.na(theta)) %>%
  summarize(
    turnangle_sum = sum(abs(theta)),
    turnangle_mean = mean(abs(theta)),
    turnangle_sd = sd(abs(theta)), # should absolute value be taken?
    approx_en = ApEn(theta),
    sample_en = SampEn(theta)
  )


# Join path length and entropy metrics to rest of figure data

figsummary <- pathlen_summary %>%
  left_join(figsummary, by = c("id", "session", "block", "trial")) %>%
  left_join(entropy_summary, by = c("id", "session", "block", "trial"))



### Preprocess tracing data prior to accuracy analysis ###

# Perform initial data cleaning

responsedat <- tracings %>%
  group_by(id, session, block, trial) %>%
  mutate(
    time = time - time[1], # ensure first time is 0.0
    dist = line_length(lag(x), lag(y), x, y)
  ) %>%
  filter(is.na(dist) | dist > 0) %>%  # drop repeated points
  mutate(mt_clip = max(time))


# Get tracing time and path length summary for all figures

tracesummary <- responsedat %>%
  summarize(
    mt_clip = mt_clip[1],
    PLresp = sum(dist, na.rm = TRUE)
  ) %>%
  mutate(vresp = PLresp / mt_clip)



### Perform accuracy analyses for tracing responses ###

# Join figure and tracing data into a single data frame

stim_temp <- frames %>%
  group_by(id, session, block, trial) %>%
  mutate(frame = 1:n(), time = time - time[1])

resp_temp <- responsedat %>%
  group_by(id, session, block, trial) %>%
  filter(n() >= 20) %>% # drop tracings w/ less than 20 unique frames
  mutate(frame = 1:n()) %>% # NOTE: doesn't account for dropped/filtered frames
  rename(trace.x = x, trace.y = y, trace.time = time)

figtrace <- stim_temp %>%
  full_join(resp_temp, by = c("id", "session", "block", "trial", "frame")) %>%
  filter(!is.na(trace.time[1])) %>%
  select(1:trial, frame, x, y, time, trace.x, trace.y, trace.time)


# Downsample paths to be equal lengths (shorten whichever is longest)

figtrace <- figtrace %>%
  group_by(id, session, block, trial) %>%
  group_modify(~ match_lengths(.$x, .$y, .$trace.x, .$trace.y))


# Get raw error measures

err_raw <- figtrace %>%
  mutate(err = line_length(x, y, trace.x, trace.y)) %>%
  summarize(
    raw_err_tot = sum(err),
    raw_err_mean = mean(err),
    raw_err_sd = sd(err),
    raw_err_paired_sd = paired.sd(err)
  )


# Get Procrustes error measures

err_proc <- figtrace %>%
  group_modify(~ procrustes2df(.$x, .$y, .$trace.x, .$trace.y)) %>%
  mutate(err = line_length(x, y, proc.x, proc.y)) %>%
  summarize(
    translation = translation[1],
    scale = scale[1],
    rotation = rotation[1],
    shape_err_tot = sum(err),
    shape_err_mean = mean(err),
    shape_err_sd = sd(err),
    shape_err_paired_sd = paired.sd(err)
  )


# Get dynamic time warping error measures

figtrace_dtw <- figtrace %>%
  group_modify(~ dtw2df(.$x, .$y, .$trace.x, .$trace.y))

err_dtw <- figtrace_dtw %>%
  mutate(err = line_length(x_w, y_w, trace.x_w, trace.y_w)) %>%
  summarize(
    dtw_err_tot = sum(err),
    dtw_err_mean = mean(err),
    dtw_err_sd = sd(err),
    dtw_err_paired_sd = paired.sd(err)
  )


# Get dynamic time warping + Procrustes error measures

err_dtw_proc <- figtrace_dtw %>%
  group_modify(~ procrustes2df(.$x_w, .$y_w, .$trace.x_w, .$trace.y_w)) %>%
  mutate(err = line_length(x, y, proc.x, proc.y)) %>%
  summarize(
    translation_dtw = translation[1],
    scale_dtw = scale[1],
    rotation_dtw = rotation[1],
    shape_dtw_err_tot = sum(err),
    shape_dtw_err_mean = mean(err),
    shape_dtw_err_sd = sd(err),
    shape_dtw_err_paired_sd = paired.sd(err)
  )


# Join all tracing summary data together

tracesummary <- tracesummary %>%
  left_join(err_raw, by = c("id", "session", "block", "trial")) %>%
  left_join(err_proc, by = c("id", "session", "block", "trial")) %>%
  left_join(err_dtw, by = c("id", "session", "block", "trial")) %>%
  left_join(err_dtw_proc, by = c("id", "session", "block", "trial"))



### Merge summarized figure data with task data ###

# Generate proper id key for joining figure data to task data

taskdat <- taskdat %>%
  mutate(fig_id = as.numeric(gsub("^p(\\d+)_.*", "\\1", figure_file)))

merge_key <- c(
  "fig_id" = "id",
  "session_num" = "session",
  "block_num" = "block",
  "trial_num" = "trial"
)


# Actually join data

fulldat <- taskdat %>%
  left_join(figsummary, by = merge_key) %>%
  left_join(tracesummary, by = merge_key)

fulldat$fig_id <- NULL
