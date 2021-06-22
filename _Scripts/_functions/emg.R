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
