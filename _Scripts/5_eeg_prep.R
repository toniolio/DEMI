############################################
### DEMI EEG processing & merging script ###
############################################

### Import required packages & functions ###

library(readr)
library(dplyr)
library(data.table)

source("./_Scripts/_settings.R")
source("./_Scripts/_functions/eeg.R")



### Load cached behavioural and EEG data ###

bdat <- readRDS("./_Scripts/_rds/bdat2.rds")
eeg_data <- readRDS("./_Scripts/_rds/eeg_imported.Rds")

# NOTE: should check and see how this differs from the MNE coords we use
eegcoords <- read_csv("./_Data/eeg/BESA-81.csv", col_types = cols())



### Perform preprocessing on EEG files to prepare for GAMs ###

# Create output path if it doesn't already exist
participants_rds_path <- "./_Scripts/_rds/participants/"
if (!dir.exists(participants_rds_path)) {
  dir.create(participants_rds_path)
}

subject_ids <- names(eeg_data)

for (id in subject_ids) {

  # Load in epoched and downsampled EEG data from list
  id_num <- as.numeric(gsub("\\D", "", id))
  cat("\n### Processing EEG signal for participant", id_num, "###\n")
  epoched <- eeg_data[[id]]$eeg

  # Convert "time" to seconds
  epoched$time <- as.double(epoched$time) / 1000

  # Perform wavelet decomposition on baseline data
  cat("\n# Performing wavelet decomposition on baseline epochs\n")
  is_baseline <- epoched$epoch == "baseline"
  baseline_wt <- wavelet_transform_id(
    eeg_signal = epoched[is_baseline, ],
    freqs = wt_frequencies,
    trim = c(1, 1),
    downsample = TRUE
  )

  # Perform wavelet decomposition on tracing data
  cat("\n# Performing wavelet decomposition on tracing epochs\n")
  is_tracing <- epoched$epoch == "tracing"
  tracing_wt <- wavelet_transform_id(
    eeg_signal = epoched[is_tracing, ],
    freqs = wt_frequencies,
    trim = c(1, 1),
    downsample = TRUE
  )

  # Perform wavelet decomposition on post-tracing data
  cat("\n# Performing wavelet decomposition on post-tracing epochs\n")
  is_post_trace <- epoched$epoch == "post_trace"
  post_trace_wt <- wavelet_transform_id(
    eeg_signal = epoched[is_post_trace, ],
    freqs = wt_frequencies,
    trim = c(1, 1),
    downsample = TRUE
  )

  # Merge wavelet-decomposed data from the tracing and post-tracing epochs
  epoched_wt <- data.table::rbindlist(list(
    tracing = tracing_wt,
    post_trace = post_trace_wt
  ), idcol = "epoch")

  # Remove intermediate objects to free up memory
  rm(epoched, tracing_wt, post_trace_wt)

  # Decibel-normalize the power data using baseline epoch
  # Note: If using CSD, resulting units will be (mV/m^2)^2/Hz (dB)
  cat("\n# Decibel-normalizing power using mean baseline epoch power...\n")
  freq_key <- c("trial", "chan", "freq")
  baseline_pwr <- baseline_wt[, .(avg_pwr = mean(power)), by = freq_key]
  epoched_wt <- epoched_wt[baseline_pwr, on = freq_key, baseline_pwr := avg_pwr]
  epoched_wt[, powerdb := 10 * (log10(power) - log10(baseline_pwr))]
  epoched_wt[, baseline_pwr := NULL]

  # Save Rds of data for future modelling and clear data objects from memory
  cat("\n# Saving data to .Rds...\n")
  outfile <- paste0(participants_rds_path, id, "_eeg_processed.rds")
  setcolorder(epoched_wt, c("trial", "epoch", "chan", "freq"))
  saveRDS(epoched_wt, file = outfile)
  rm(epoched_wt, baseline_wt, baseline_pwr)
  cat("\n### Participant", id_num, "successfully processed! ###\n\n")
}



### Combine full dataset ###

# Free up any memory taken up by wavelet transformation

rm(eeg_data)


# Get list of all wavelet-transformed EEG .Rds files

eeg_wt_files <- list.files(
  participants_rds_path, pattern = "*.rds",
  full.names = TRUE
)


# Merge wavelet-transformed EEG data from all participants

all_dat <- NULL

cat("\n### Combining wavelet-transformed EEG for all participants ###\n")
for (f in eeg_wt_files) {

  id_num <- as.integer(gsub("\\D", "", basename(f)))
  cat("\n# Merging EEG data from participant", id_num, "...\n")

  if (f == eeg_wt_files[1]) {
    all_dat <- readRDS(f)[, id := id_num]

  } else {
    all_dat <- data.table::rbindlist(list(
      current = all_dat,
      new = readRDS(f)[, id := id_num]
    ))

  }
}


# Make some variable integers to reduce memory demands

all_dat[, trial := as.integer(trial)]
all_dat[, freq := as.integer(freq)]
all_dat[, epoch := ifelse(epoch == "tracing", 1L, 2L)]


# Merge behavioural data with EEG data

cols_to_merge <- c(
  "id", "group", "trial", "condition", "rep", "complexity",
  "avg_velocity", "error", "vresp", "accuracy_rating"
)

bdat_merge <- bdat %>%
  mutate(
    id = as.integer(as.character(participant)),
    trial = as.integer(trial_num + (block_num - 1) * 20),
    rep = as.integer(rep),
    accuracy_rating = as.integer(accuracy_rating)
  ) %>%
  select(all_of(cols_to_merge)) %>%
  arrange(id) %>%
  as.data.table()

all_dat <- all_dat[bdat_merge,
  on = c("id", "trial"),
  `:=`(
    group = group,
    condition = condition,
    rep = rep,
    complexity = complexity,
    avg_velocity = avg_velocity,
    error = error,
    vresp = vresp,
    accuracy_rating = accuracy_rating
  )
]

setcolorder(all_dat, c("id", "group", "trial", "condition", "rep"))


# Add spatial electrode coordinate info for each channel

# NOTE: Shouldn't do this here: different scripts use coords in different
# formats, so would be better to join coords in the correct format in
# scripts before they're needed

eegcoords <- eegcoords %>%
  mutate(
    lat = 90 - (asin(z) * (180 / pi)),
    long = atan2(y, x) * (180 / pi)
  )

all_dat <- all_dat[eegcoords, on = c("chan"), `:=`(lat = lat, long = long)]


# Write out giant merged data frame to Rds for modelling and plotting

saveRDS(all_dat, file = "./_Scripts/_rds/all_dat.rds")
