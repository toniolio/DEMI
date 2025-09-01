### Custom EEG processing functions for DEMI TraceLab analysis ###

### Import required libraries ###

library(dplyr)
library(tidyr)
library(tibble)
library(eeguana)
library(biwavelet)



### Utility Functions ###

# Prints percentage progress for slow operations in loops

print_progress_pct <- function(n, v, every = 10) {

  index <- match(n, v)
  pct <- 100 * index / length(v)
  last_pct <- 100 * (index - 1) / length(v)
  step_change <- floor(pct / every) != floor(last_pct / every)

  if (index == 1) {
    cat("Progress: ")
  } else if (index == length(v)) {
    cat("100%\n")
  } else if (index > 1 & step_change) {
    cat(paste0(floor(pct / every) * every, "%... "))
  }
}


### Wavelet Transformation Functions ###

# Converts frequencies in Hz into the step format required by biwavelet::wt()

hz_steps_to_period <- function(time, freqs) {

  # Get fourier factor for data (needed for Hz -> period calculations)
  n <- length(time)
  dt <- biwavelet::check.datum(cbind(time, time))$dt
  k <- seq_len(floor(0.5 * n)) * 2 * pi / (n * dt)
  k <- c(0, k, -k[floor(0.5 * (n - 1)):1])
  ffac <- biwavelet::wt.bases("morlet", k, dt * 2, -1)$fourier.factor

  # Convert range of frequency values into step sizes for wt()
  freqs <- sort(freqs, decreasing = TRUE)
  first_step <- 1 / (ffac * freqs[1])
  numerator <- log(1 / (ffac * freqs * first_step))
  period_steps <- numerator / (log(2) * 0:(length(freqs) - 1))
  period_steps[1] <- first_step

  period_steps
}


# Performs wavelet decomposition on a given signal, optionally trimming padding

wavelet_transform <- function(time, signal, freqs, trim = c(1, 1)) {

  # Transform requested frequencies (in Hz) to wavelet step sizes (in period)
  steps <- hz_steps_to_period(time, freqs)

  # Actually do wavelet transform using biwavelet
  wt_results <- biwavelet::wt(
    cbind(time, signal),
    dj = steps,
    s0 = steps[1],
    J1 = length(freqs) - 1,
    do.sig = FALSE
  )

  # Convert period to Hz
  freq <- round(1 / wt_results$period, 8)

  # Construct data table from biwavelet results
  dt_wt <- as.data.table(t(wt_results$power))
  setnames(dt_wt, as.character(freq))
  dt_wt$time <- time

  # Trim padding from data after transformation
  new_start <- min(dt_wt$time) + trim[1]
  new_end <- max(dt_wt$time) - trim[2]
  dt_wt <- dt_wt[time >= new_start & time <= new_end, ]
  if (new_start >= new_end) {
    stop("'trim' parameter values remove all data (check units?)")
  }

  dt_wt
}


# Performs wavelet decomposition for all channels and trials of a specific
# epoch time for a single participant, optionally trimming padding and
# downsampling data by a factor of 10. Can be very slow without downsampling.

wavelet_transform_id <- function(
  eeg_signal, freqs, trim, baseline = NULL, downsample = TRUE) {

  trials <- unique(eeg_signal$trial)
  ch_names <- setdiff(names(eeg_signal), c("trial", "epoch", "time"))

  # For each trial and channel, perform wavelet decomposition
  wt_by_trial <- lapply(trials, function(tnum) {

    print_progress_pct(tnum, trials)
    in_trial <- eeg_signal$trial == tnum

    wt_by_ch <- lapply(ch_names, function(ch) {

      # Perform wavelet transform on channel & trim padding
      dt_wt <- wavelet_transform(
        time = eeg_signal$time[in_trial],
        signal = eeg_signal[[ch]][in_trial],
        freqs = freqs,
        trim = trim
      )

      # If a baseline window was provided, get & append the mean baseline power
      if (length(baseline) == 2) {
        t1 <- baseline[1]
        t2 <- baseline[2]
        baseline_pwr <- dt_wt[time >= t1 & time <= t2, lapply(.SD, mean)]
        baseline_pwr$time <- NA  # NA here means baseline average
        dt_wt <- data.table::rbindlist(list(dt_wt, baseline_pwr))
      }

      # Downsample wavelet transformed data by 10x
      if (downsample) {
        keep <- dt_wt[, .I %% 10 == 1] | is.na(dt_wt$time)
        dt_wt <- dt_wt[keep]
      }

      dt_wt
    })
    names(wt_by_ch) <- ch_names

    # Bind the transformed data for all channels into a single frame
    data.table::rbindlist(wt_by_ch, idcol = "chan")

  })
  names(wt_by_trial) <- as.character(trials)

  # Bind the transformed data from each trial together into a single frame
  cat("Merging transformed data and reshaping to long format...\n")
  out <- data.table::rbindlist(wt_by_trial, idcol = "trial")

  # Transform merged data from wide to long format
  out[, trial := as.numeric(trial)]
  out <- data.table::melt(
    out,
    id.vars = c("trial", "chan", "time"),
    variable.factor = FALSE,
    variable.name = "freq",
    value.name = "power"
  )
  out[, freq := as.numeric(freq)]

  # If baseline provided, calculate dB-normalized power per trial/chan/frequency
  if (length(baseline) == 2) {
    freq_key <- c("trial", "chan", "freq")
    baseline_pwr <- out[is.na(time)]
    setnames(baseline_pwr, "power", "avg_power")
    out <- out[!is.na(time)][baseline_pwr, on = freq_key, baseline := avg_power]
    out[, powerdb := 10 * (log10(power) - log10(baseline))]
    out[, baseline := NULL]
  }

  out
}



### EEG Trigger Processing Functions ###

# Processes and updates EEG triggers, dropping triggers from practice trials and
# trials already excluded in the behavioural data and adding new
# 'real_trace_start' and 'real_trace_end' triggers based on the actual tracing
# starts and ends (as determined by the filtered tracing data)

update_events <- function(eeg_events, epoch_offsets) {

  # Gather all triggers from EEG data and compute trial numbers & figure speeds

  triggers <- as_tibble(eeg_events) %>%
    mutate(
      bad_start = .description == "red_on" & lag(.description) != "stim_on",
      duration = ifelse(bad_start, NA, .initial - lag(.initial)),
      trial = cumsum(.description == "stim_on" | bad_start)
    ) %>%
    group_by(trial) %>%
    mutate(
      fig_duration = duration[.description == "red_on"],
      practice = !is.na(fig_duration) & fig_duration > 3000,
      bad_start = any(bad_start)
    )


  # Drop triggers from practice trials and trials missing 'stim_on' trigger

  triggers <- subset(triggers, !practice) %>%
    ungroup() %>%
    mutate(
      trial = cumsum(is.na(lag(trial)) | trial > lag(trial))
    ) %>%
    filter(!bad_start) %>%
    select(-c(practice, bad_start))


  # Calculate actual tracing starts/end epochs using behavioural data

  trace_epoch_info <- triggers %>%
    filter(.description %in% c("trace_start", "trace_end")) %>%
    select(c(.id, .type, .channel, trial, .description, .initial)) %>%
    spread(key = .description, value = .initial) %>%
    rename(start = trace_start, end = trace_end) %>%
    right_join(epoch_offsets, by = c("trial" = "trial_count")) %>%
    filter(!is.na(end))

  new_starts <- trace_epoch_info %>%
    mutate(
      .initial = ifelse(physical, start + real_start * 1000, start),
      .final = .initial,
      .description = "real_trace_start"
    ) %>%
    select(c(.id, .type, .description, .initial, .final, .channel, trial))

  new_ends <- trace_epoch_info %>%
    mutate(
      .initial = ifelse(physical, start + real_end * 1000, end),
      .final = .initial,
      .description = "real_trace_end"
    ) %>%
    select(c(.id, .type, .description, .initial, .final, .channel, trial))


  # Sanity-check that expected epoch lengths line up with actual ones

  bad_epoch_estimates <- trace_epoch_info %>%
    mutate(
      end_diff = abs(end - (start + end_trigger * 1000)),
      bad_estimate = end_diff > 150
    ) %>%
    filter(bad_estimate)

  if (nrow(bad_epoch_estimates) > 0) {
    n_bad <- nrow(bad_epoch_estimates)
    cat(paste0(
      "\nNote: Encountered ", n_bad, " bad tracing epoch estimates.\n\n"
    ))
  }


  # Make sure figure durations match up between task and EEG data

  fig_durations_compare <- triggers %>%
    filter(.description == "trace_start") %>%
    left_join(epoch_offsets, by = c("trial" = "trial_count")) %>%
    mutate(
      stimulus_mt = round(stimulus_mt * 1000),
      mt_diff = fig_duration - stimulus_mt,
      large_diff = !is.na(mt_diff) & abs(mt_diff) > 250
    ) %>%
    filter(large_diff)

  if (nrow(fig_durations_compare) > 0) {
    cat("Large differences in stimulus duration between EDF and task data.\n")
    cat("The files may not be aligned properly.\n")
    for (i in seq_len(nrow(fig_durations_compare))) {
      cat(paste0(
        "  - Trial ", fig_durations_compare$trial[i], ": ",
        "Expected = ", fig_durations_compare$stimulus_mt[i], " ms, ",
        "EEG = ", fig_durations_compare$fig_duration[i], " ms\n"
      ))
    }
  }


  # Bind new tracing start/end triggers to full EEG events table

  new_events <- triggers %>%
    select(c(.id, .type, .description, .initial, .final, .channel, trial)) %>%
    mutate(.initial = as.integer(.initial), .final = as.integer(.final)) %>%
    bind_rows(new_starts) %>%
    bind_rows(new_ends) %>%
    arrange(.initial)


  # Filter out triggers from any trials already dropped from behavioural data

  new_events <- new_events %>%
    group_by(trial) %>%
    mutate(bad_trial = !any(.description == "real_trace_end")) %>%
    ungroup() %>%
    filter(!bad_trial) %>%
    select(-c(bad_trial)) %>%
    arrange(.initial)


  # Coerce updated events data back into eeguana's custom 'events_tbl' class

  eeg_hz <- attributes(eeg_events$.initial)$sampling_rate
  new_events <- new_events %>%
    mutate(
      .initial = as_sample_int(.initial, sampling_rate = eeg_hz, .unit = "ms"),
      .final = as_sample_int(.final, sampling_rate = eeg_hz, .unit = "ms")
    )

  new_events <- data.table(new_events)
  data.table::setattr(new_events, "class", c("events_tbl", class(new_events)))

  new_events
}
