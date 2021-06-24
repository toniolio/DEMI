### Custom EMG processing functions for DEMI TraceLab analysis ###

### Import required libraries ###

library(dplyr)
library(tidyr)
library(tibble)
library(eeguana)


### EMG preprocessing functions ###

emg_bandpass <- function(x, low_freq, high_freq, srate) {

  # NOTE: Filter design and parameters are based on recommendations in
  # "Systematic evaluation of EMG data" (Schepens, 2018)

  # Create separate low/highpass butterworth filters for doing initial bandpass
  highpass_w <- low_freq / (srate / 2)
  lowpass_w <- high_freq / (srate / 2)
  bfilt_hi <- signal::butter(9, W = highpass_w, type = "high")
  bfilt_lo <- signal::butter(3, W = lowpass_w, type = "low")

  # Apply both filters to the EMG signal
  signal::filtfilt(bfilt_lo, signal::filtfilt(bfilt_hi, x))
}


emg_rectify <- function(x) {
  abs(x - median(x))
}


emg_smoothing <- function(x, high_freq, srate) {

  # Create a lowpass filter to apply to the rectified EMG signal
  lowpass_w <- high_freq / (srate / 2)
  bfilt_lo <- signal::butter(6, W = lowpass_w, type = "low")

  # Apply the filter to the EMG signal
  signal::filtfilt(bfilt_lo, x)
}



### eeguana wrangling functions ###

append_epochs <- function(signal, events) {

  # Get event names and onsets per trial and join to signal
  join_key <- c(".id" = ".id", ".sample" = ".initial")
  events <- events %>% select(.id, .initial, .description)
  signal <- left_join(signal, events, by = join_key)

  # Create a numeric column for epoch
  signal <- signal %>%
    group_by(.id) %>%
    filter(.sample < max(.sample)) %>%  # drop last sample for each trial
    mutate(epoch_num = cumsum(!is.na(.description)))

  # Get name/number map for epochs
  name_map <- signal %>%
    filter(!is.na(.description)) %>%
    mutate(epoch = .description) %>%
    select(c(.id, epoch_num, epoch))

  # Convert the numeric epoch column into a named factor column
  signal <- ungroup(signal) %>%
    left_join(name_map, by = c(".id", "epoch_num")) %>%
    mutate(epoch = as.factor(epoch)) %>%
    select(-c(epoch_num, .description))

  signal
}
