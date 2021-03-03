##########################################
### Tracelab data preprocessing script ###
##########################################


### Import required packages and functions ###

library(dplyr)
library(TSEntropies)

source("./_Scripts/_settings.R")
source("./_Scripts/_functions/complexity.R")
source("./_Scripts/_functions/filters.R")
source("./_Scripts/_functions/physical.R")
source("./_Scripts/_functions/visualization.R")


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



### Filter tracing data prior to accuracy analysis ###

# Get origin point for each trial

origins <- points %>%
  group_by(id, session, block, trial) %>%
  summarize(origin.x = x[1], origin.y = y[1])


# Calculate useful metrics for data filtering and drop all repeated points

responsedat <- tracings %>%
  left_join(origins, by = c("id", "session", "block", "trial")) %>%
  group_by(id, session, block, trial) %>%
  mutate(
    time = time - time[1], # ensure first time is 0.0
    dist = line_length(lag(x), lag(y), x, y),
    origin.dist = line_length(origin.x, origin.y, x, y)
  ) %>%
  filter(is.na(dist) | dist > 0)


# Flag and drop any trials with no formed shape

fig_info <- frames %>%
  group_by(id, session, block, trial) %>%
  summarize(fig_max_size = max(c(max(x) - min(x), max(y) - min(y))))

no_shape_trials <- responsedat %>%
  summarize(
    duration = max(time),
    max_size = max(c(max(x) - min(x), max(y) - min(y)))
  ) %>%
  left_join(fig_info, by = c("id", "session", "block", "trial")) %>%
  mutate(
    size_ratio = fig_max_size / max_size,
    no_shape = size_ratio > no_shape_params$min_size_ratio
  ) %>%
  filter(no_shape)

if (plot_filters) {
  plot_trials(no_shape_trials, responsedat, outdir = "./filters/no_shape")
}
trial_key <- c("id", "session", "block", "trial")
responsedat <- anti_join(responsedat, no_shape_trials, by = trial_key)


# Trim extra points following failed trial end

responsedat <- responsedat %>%
  mutate(timediff = ifelse(is.na(lag(time)), 0, time - lag(time))) %>%
  mutate(done = trial_done(origin.dist, timediff, done_filter_params))

failed_end_trials <- responsedat %>%
  summarize(
    samples = n(),
    num_done = sum(done),
    mt_diff = max(time[!done]) - max(time),
    prop_remaining = 1 - (sum(done) / n())
  ) %>%
  filter(num_done > 0)

if (plot_filters) {
  plot_trials(failed_end_trials, responsedat, "done", "./filters/done")
}
responsedat <- subset(responsedat, !done)


# Flag and remove glitch points during tracings

responsedat <- responsedat %>%
  mutate(angle_diff = (get_angle_diffs(x - lag(x), y - lag(y)) / pi) * 180) %>%
  mutate(
    glitch = is_glitch(x, y, angle_diff, origin.dist, glitch_filter_params)
  )

glitch_trials <- responsedat %>%
  summarize(
    samples = n(),
    glitches = sum(glitch)
  ) %>%
  filter(glitches > 0)

if (plot_filters) {
  plot_trials(glitch_trials, responsedat, "glitch", "./filters/glitch")
}
responsedat <- subset(responsedat, !glitch)


# Flag and remove false start samples

responsedat <- responsedat %>%
  mutate(timediff = ifelse(is.na(lag(time)), 0, time - lag(time))) %>%
  mutate(false_start = false_start(origin.dist, timediff, false_start_params))

false_starts <- responsedat %>%
  summarize(
    samples = n(),
    flagged = sum(false_start)
  ) %>%
  filter(flagged > 0)

if (plot_filters) {
  plot_trials(false_starts, responsedat, "false_start", "./filters/false_start")
}
responsedat <- subset(responsedat, !false_start)


# Try to flag and remove hand noise samples

responsedat <- responsedat %>%
  mutate(
    timediff = ifelse(is.na(lag(time)), 0, time - lag(time)),
    angle_diff = (get_angle_diffs(x - lag(x), y - lag(y)) / pi) * 180,
    dist = line_length(lag(x), lag(y), x, y)
  ) %>%
  mutate(
    hnoise = hand_noise(
      x, y, timediff, angle_diff, origin.dist, hand_noise_params
    )
  )

hand_noise_trials <- responsedat %>%
  summarize(
    samples = n(),
    flagged = sum(hnoise)
  ) %>%
  filter(flagged > 0)

if (plot_filters) {
  plot_trials(hand_noise_trials, responsedat, "hnoise", "./filters/hand_noise")
}
responsedat <- subset(responsedat, !hnoise)


# Flag and drop trials with incomplete tracings

origin_x <- screen_res[1] / 2

fig_info <- frames %>%
  group_by(id, session, block, trial) %>%
  summarize(
    fig_max_size = max(c(max(x) - min(x), max(y) - min(y))),
    fig_lat_shift = ((min(x) - origin_x) / (max(x) - min(x)) + 0.5) * 2
  ) %>%
  left_join(pathlen_summary, by = c("id", "session", "block", "trial"))

incomplete_trials <- responsedat %>%
  mutate(seglen = line_length(lag(x), lag(y), x, y)) %>%
  summarize(
    end_gap = abs(origin.dist[n()] - origin.dist[1]),
    max_size = max(c(max(x) - min(x), max(y) - min(y))),
    lat_shift = ((min(x) - origin_x) / (max(x) - min(x)) + 0.5) * 2,
    trace_len = sum(seglen, na.rm = TRUE)
  ) %>%
  left_join(fig_info, by = c("id", "session", "block", "trial")) %>%
  mutate(
    size_ratio = fig_max_size / max_size,
    shift_diff = abs(lat_shift) - abs(fig_lat_shift),
    len_ratio = PLstim / trace_len
  ) %>%
  mutate(
    incomplete = is_incomplete(
      end_gap, size_ratio, len_ratio, shift_diff, incomplete_params
    )
  ) %>%
  filter(incomplete)

if (plot_filters) {
  plot_trial_paths(
    incomplete_trials, frames, responsedat, "./filters/incomplete"
  )
}
trial_key <- c("id", "session", "block", "trial")
responsedat <- anti_join(responsedat, incomplete_trials, by = trial_key)


# Flag and drop trials with excessive time or distance gaps

responsedat <- responsedat %>%
  mutate(
    timediff = ifelse(is.na(lag(time)), 0, time - lag(time)),
    angle_diff = (get_angle_diffs(x - lag(x), y - lag(y)) / pi) * 180,
    dist = line_length(lag(x), lag(y), x, y)
  ) %>%
  mutate(turn_sum = abs(angle_diff) + abs(lead(angle_diff))) %>%
  mutate(
    gap = is_gap(dist, timediff, turn_sum, gap_filter_params)
  )

gap_trials <- responsedat %>%
  summarize(
    samples = n(),
    longest_pause = max(timediff),
    longest_jump = max(dist, na.rm = TRUE),
    any_gaps = any(gap)
  ) %>%
  filter(any_gaps)

if (plot_filters) {
  plot_trials(gap_trials, responsedat, "gap", "./filters/large_gap")
}
responsedat <- anti_join(responsedat, gap_trials, by = trial_key)


# Flag (but don't remove) remaining trials where tracing hits edge of screen

edge_trials <- responsedat %>%
  summarize(
    hit_top_edge = any(y == 0),
    hit_bottom_edge = any(y == (screen_res[2] - 1)),
    hit_left_edge = any(x == 0),
    hit_right_edge = any(x == (screen_res[1] - 1))
  ) %>%
  filter(hit_top_edge | hit_bottom_edge | hit_left_edge | hit_right_edge)



### Prepare tracing & figure data for accuracy analysis ###

# Get duration and path length for all tracings

tracesummary <- responsedat %>%
  mutate(mt_clip = max(time)) %>%
  summarize(
    mt_clip = mt_clip[1],
    PLresp = sum(dist, na.rm = TRUE)
  ) %>%
  mutate(vresp = PLresp / mt_clip)


# Reinterpolate dropped/filtered points in figures by time

resp_interpolated <- responsedat %>%
  group_modify(~ reinterpolate(.$x, .$y, .$time))


# Join figure and tracing data into a single data frame

stim_temp <- frames %>%
  group_by(id, session, block, trial) %>%
  mutate(frame = 1:n(), time = time - time[1])

resp_temp <- resp_interpolated %>%
  mutate(frame = 1:n()) %>%
  rename(trace.x = x, trace.y = y)

figtrace <- stim_temp %>%
  full_join(resp_temp, by = c("id", "session", "block", "trial", "frame")) %>%
  filter(!is.na(trace.x[1])) %>%
  select(1:trial, frame, x, y, trace.x, trace.y)



### Perform accuracy analyses for raw tracing responses ###

# Resample stimulus frames to match trace lengths

figtrace_eq <- figtrace %>%
  group_by(id, session, block, trial) %>%
  group_modify(~ match_lengths(.$x, .$y, .$trace.x, .$trace.y))


# Get raw error measures

err_raw <- figtrace_eq %>%
  mutate(err = line_length(x, y, trace.x, trace.y)) %>%
  summarize(
    raw_err_tot = sum(err),
    raw_err_mean = mean(err),
    raw_err_sd = sd(err),
    raw_err_paired_sd = paired.sd(err)
  )


# Get Procrustes error measures

err_proc <- figtrace_eq %>%
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



### Perform accuracy analyses for equidistant tracing responses ###

# Resample stimulus frames to match trace lengths & reinterpolate tracing points
# to be evenly spaced along path

figtrace_equidist <- figtrace %>%
  group_by(id, session, block, trial) %>%
  group_modify(~ match_lengths(.$x, .$y, .$trace.x, .$trace.y, equidist = TRUE))


# Get raw error measures

err_eqd <- figtrace_equidist %>%
  mutate(err = line_length(x, y, trace.x, trace.y)) %>%
  summarize(
    eqd_err_tot = sum(err),
    eqd_err_mean = mean(err),
    eqd_err_sd = sd(err),
    eqd_err_paired_sd = paired.sd(err)
  )


# Get Procrustes error measures

err_eqd_proc <- figtrace_equidist %>%
  group_modify(~ procrustes2df(.$x, .$y, .$trace.x, .$trace.y)) %>%
  mutate(err = line_length(x, y, proc.x, proc.y)) %>%
  summarize(
    translation_eqd = translation[1],
    scale_eqd = scale[1],
    rotation_eqd = rotation[1],
    shape_eqd_err_tot = sum(err),
    shape_eqd_err_mean = mean(err),
    shape_eqd_err_sd = sd(err),
    shape_eqd_err_paired_sd = paired.sd(err)
  )



### Perform accuracy analyses for dtw-processed tracing responses ###

# Resample stimulus frames to match trace lengths using dtw

figtrace_dtw <- figtrace %>%
  group_by(id, session, block, trial) %>%
  group_modify(~ dtw2df(.$x, .$y, .$trace.x, .$trace.y))


# Get dynamic time warping error measures

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



### Merge summarized figure data with task data ###

# Join all tracing summary data together

tracesummary <- tracesummary %>%
  left_join(err_raw, by = c("id", "session", "block", "trial")) %>%
  left_join(err_proc, by = c("id", "session", "block", "trial")) %>%
  left_join(err_eqd, by = c("id", "session", "block", "trial")) %>%
  left_join(err_eqd_proc, by = c("id", "session", "block", "trial")) %>%
  left_join(err_dtw, by = c("id", "session", "block", "trial")) %>%
  left_join(err_dtw_proc, by = c("id", "session", "block", "trial"))


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
