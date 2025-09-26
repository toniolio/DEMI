######################################################
### DEMI EEG processing — 05 EEG prep for analysis ###
######################################################

# performs emg and eeg analyses
# emg analysis: detects whether there are imagery trials with excessive movement
# eeg analysis: wavelet analyses
# writes per participant data for full time-frequency analyses (plotted in 06)
# also writes participant-aggregated band data (drops time) for scripts 07+

suppressPackageStartupMessages({
	library(data.table)
	library(dplyr)
	library(readr)
	library(purrr)
	library(tidyr)
	library(tibble)
})

# Unconditional settings & helpers (hard fail if missing)
source("./_Scripts/_settings.R")
source("./_Scripts/_functions/eeg.R")
source("./_Scripts/_functions/emg.R")

options(dplyr.summarise.inform = FALSE)

# ---- Inputs (hard-fail if missing) ----
if (!file.exists("./_Scripts/_rds/eeg_imported.Rds")) stop("Missing: ./_Scripts/_rds/eeg_imported.Rds")
if (!file.exists("./_Scripts/_rds/bdat2.rds"))        stop("Missing: ./_Scripts/_rds/bdat2.rds")
if (!file.exists("./_Data/eeg/BESA-81.csv"))          stop("Missing: ./_Data/eeg/BESA-81.csv")

eeg_data <- readRDS("./_Scripts/_rds/eeg_imported.Rds")
bdat     <- readRDS("./_Scripts/_rds/bdat2.rds")
coords_path <- "./_Data/eeg/BESA-81.csv"

# ---- Outputs ----
participants_rds_path <- "./_Scripts/_rds/participants/"
agg_dir               <- "./_Scripts/_rds/participants_agg/"
out_bands_path        <- "./_Scripts/_rds/dat_bands_notime.rds"

if (!dir.exists(participants_rds_path)) dir.create(participants_rds_path, recursive = TRUE)
if (!dir.exists(agg_dir))               dir.create(agg_dir, recursive = TRUE)

# Optional (OFF): build monolithic all_dat.rds (discouraged)
also_build_all_dat <- FALSE

# -------- 0) EMG analysis & trial exclusion --------

cat("\n### Preprocessing EMG data for imagery participants ###\n")
emg_activity <- map_df(names(eeg_data), function(id) {
	emg <- eeg_data[[id]]$emg
	if (all(is.na(emg))) return(NULL)  # non‑imagery participants
	id_num <- as.numeric(gsub("\\D", "", id))
	cat(" - Processing EMG signal for participant", id_num, "...\n")

	# Clean up epoch labels and append to signal
	cleaned_events <- emg$.events %>%
		filter(.description != "trace_start") %>%
		mutate(.description = gsub("real_", "", .description)) %>%
		filter(.description != "trace_end")
	emg$.signal <- append_epochs(emg$.signal, cleaned_events)

	# Wrangle, rectify, and smooth EMG channels; support underscore/dot names
	has_us <- all(c("EMG_L","EMG_A") %in% names(emg$.signal))
	has_dot <- all(c("EMG.L","EMG.A") %in% names(emg$.signal))
	if (!has_us && !has_dot) stop("Expected EMG channels EMG_L/EMG_A or EMG.L/EMG.A not found")

	emg_cleaned <- emg$.signal %>%
		{
			if (has_us) select(., trial, .sample, epoch, EMG_L, EMG_A) %>%
				rename(Lateral = EMG_L, Anterior = EMG_A) else
					select(., trial, .sample, epoch, EMG.L, EMG.A) %>%
				rename(Lateral = EMG.L, Anterior = EMG.A)
		} %>%
		gather("deltoid", "signal", Lateral, Anterior) %>%
		group_by(trial, deltoid) %>%
		mutate(cleaned = emg_rectify(signal)) %>%
		mutate(cleaned = channel_dbl(emg_smoothing(cleaned, 30, srate = 1000)))

	# Median absolute deviation per trial × epoch; label trial type
	emg_summary <- emg_cleaned %>%
		group_by(trial, deltoid, epoch) %>%
		summarize(amplitude = mad(cleaned)) %>%
		mutate(trial_type = ifelse(trial >= 100, "physical", "imagery")) %>%
		add_column(id = id_num, .before = "trial")

	emg_summary
})

# Difference scores between epochs (trace_start vs stim_on), per id × deltoid
emg_diffs <- emg_activity %>%
	spread("epoch", "amplitude") %>%
	group_by(id, deltoid) %>%
	mutate(diff_score = trace_start - stim_on)

# Flag imagery trials with high activity relative to physical trials (Q1 cutoff)
emg_bads <- emg_diffs %>%
	group_by(id, deltoid) %>%
	mutate(bad_threshold = quantile(diff_score[trial_type == "physical"], prob = 0.25)) %>%
	group_by(id, trial) %>%
	summarize(bad = all(trial_type == "imagery" & diff_score >= bad_threshold)) %>%
	filter(bad == TRUE)

# -------- 1) Per‑participant time×freq shards --------

# Baseline window driven by _settings.R
if (use_pretrial_baseline) {
	baseline_window <- NULL
} else {
	baseline_window <- c(-0.5, -0.2)
}

subject_ids <- names(eeg_data)
if (!length(subject_ids)) stop("eeg_imported.Rds has no participants.")

for (id in subject_ids) {
	id_num <- as.integer(gsub("\\D", "", id))
	cat("\n### Processing participant", id_num, "###\n")

	epoched <- eeg_data[[id]]$eeg
	if (is.null(epoched) || !is.data.frame(epoched))
		stop("Unexpected structure for eeg_data[['", id, "']]$eeg")

	# EMG‑based exclusion of imagery trials
	if (exclude_bad_by_emg) {
		bad_img_trials <- emg_bads %>% filter(id == id_num) %>% pull(trial)
		bad_img_rows <- epoched$trial %in% bad_img_trials
		epoched <- epoched[!bad_img_rows, ]
	}

	# seconds
	epoched$time <- as.double(epoched$time) / 1000

	# Optional pre‑trial baseline decomposition (no dB here)
	if (use_pretrial_baseline) {
		cat("\n# Wavelet: baseline epochs\n")
		is_baseline <- epoched$epoch == "baseline"
		baseline_wt <- wavelet_transform_id(
			eeg_signal = epoched[is_baseline, ],
			freqs      = wt_frequencies,
			trim       = c(1, 1),
			baseline   = NULL,
			downsample = TRUE
		)
	}

	# Tracing & post‑tracing with within‑epoch baseline if not using pre‑trial
	cat("\n# Wavelet: tracing epochs\n")
	is_tracing <- epoched$epoch == "tracing"
	tracing_wt <- wavelet_transform_id(
		eeg_signal = epoched[is_tracing, ],
		freqs      = wt_frequencies,
		trim       = c(1, 1),
		baseline   = baseline_window,
		downsample = TRUE
	)

	cat("\n# Wavelet: post‑tracing epochs\n")
	is_post <- epoched$epoch == "post_trace"
	post_trace_wt <- wavelet_transform_id(
		eeg_signal = epoched[is_post, ],
		freqs      = wt_frequencies,
		trim       = c(1, 1),
		baseline   = baseline_window,
		downsample = TRUE
	)

	# Merge task epochs
	epoched_wt <- data.table::rbindlist(list(
		tracing = tracing_wt,
		post_trace = post_trace_wt
	), idcol = "epoch")

	# dB relative to pre‑trial baseline if used
	if (use_pretrial_baseline) {
		cat("\n# dB normalization using pre‑trial baseline\n")
		freq_key <- c("trial", "chan", "freq")
		baseline_pwr <- baseline_wt[, .(avg_pwr = mean(power)), by = freq_key]
		epoched_wt <- epoched_wt[baseline_pwr, on = freq_key, baseline_pwr := avg_pwr]
		epoched_wt[, powerdb := 10 * (log10(power) - log10(baseline_pwr))]
		epoched_wt[, baseline_pwr := NULL]
		rm(baseline_wt)
	}

	# Canonicalize epoch names in shard just before saving
	# (expect only tracing/post_trace above; map to during/after for all downstream use)
	epoched_wt[, epoch := as.character(epoch)]
	u_epoch <- sort(unique(epoched_wt$epoch))
	if (!all(u_epoch %in% c("tracing","post_trace")))
		stop("Unexpected epoch labels in shard: ", paste(u_epoch, collapse=", "))
	epoched_wt[epoch == "tracing",    epoch := "during"]
	epoched_wt[epoch == "post_trace", epoch := "after"]

	# Persist shard, GC (match original front column order; others retained)
	outfile <- file.path(participants_rds_path, paste0(id, "_eeg_processed.rds"))
	setcolorder(epoched_wt, c("trial","epoch","chan","freq"))
	saveRDS(epoched_wt, file = outfile)
	rm(epoched_wt, tracing_wt, post_trace_wt, epoched); gc()
	cat("✓ Saved ", outfile, "\n", sep="")
}

# (Optional) Build all_dat.rds (not recommended)
if (also_build_all_dat) {
	pp_files <- list.files(participants_rds_path, pattern = "_eeg_processed\\.rds$", full.names = TRUE)
	if (!length(pp_files)) stop("No per-participant shards found in ", participants_rds_path)
	all_dat <- rbindlist(lapply(pp_files, readRDS), use.names = TRUE, fill = TRUE)
	saveRDS(all_dat, "./_Scripts/_rds/all_dat.rds")
}

# -------- 2) Aggregate shards -> dat_bands_notime.rds --------

# Load channel coords (expect chan,x,y,z) and compute lat/long
coords <- data.table::fread(coords_path)
req_c <- c("chan","x","y","z")
miss_c <- setdiff(req_c, names(coords))
if (length(miss_c)) stop("Coord file missing: ", paste(miss_c, collapse=", "))
coords[, lat  := 90 - (asin(z) * (180 / pi))]
coords[, long := atan2(y, x) * (180 / pi)]
coords <- coords[, .(chan, lat, long)]

# Band helper
assign_band_inplace <- function(DT) {
	if (!("freq" %in% names(DT))) { DT[, band := NA_character_]; return(invisible(NULL)) }
	DT[, band := fifelse(freq >=  4 & freq <=  8, "theta",
						 fifelse(freq >=  9 & freq <= 12, "alpha",
						 		fifelse(freq >= 13 & freq <= 30, "beta", NA_character_))) ]
}

# Robust rep normalization
norm_rep <- function(x) {
	rep_chr <- tolower(as.character(x))
	rep_num <- suppressWarnings(as.numeric(as.character(x)))
	dplyr::case_when(
		!is.na(rep_chr) & rep_chr %in% c("random","repeated") ~ rep_chr,
		!is.na(rep_num) & rep_num == 1                        ~ "random",
		!is.na(rep_num) & rep_num %in% c(0,2)                 ~ "repeated",
		TRUE                                                  ~ NA_character_
	)
}

# Aggregate each participant shard -> *_agg.rds
pp_files <- list.files(participants_rds_path, pattern = "_eeg_processed\\.rds$", full.names = TRUE)
if (!length(pp_files)) stop("No per-participant EEG files found in ", participants_rds_path)

process_one_agg <- function(fpath) {
	id_chr <- gsub("_eeg_processed\\.rds$", "", basename(fpath))
	id_num <- as.integer(gsub("\\D", "", id_chr))

	dt <- readRDS(fpath); setDT(dt)

	need <- c("trial","epoch","chan","freq","powerdb")
	if (!all(need %in% names(dt))) stop("Missing required columns in ", fpath)

	# Attach participant from filename and merge behaviour
	dt[, participant := id_num]
	bdat_merge <- bdat %>%
		mutate(participant = as.integer(as.character(participant)),
			   trial = as.integer(trial_num + (block_num - 1L) * 20L)) %>%
		select(participant, trial, block_num, group, condition, rep,
			   complexity, avg_velocity, error, vresp, accuracy_rating) %>%
		as.data.table()
	dt <- bdat_merge[dt, on = .(participant, trial)]

	# Channel coords
	dt <- coords[dt, on = .(chan)]

	# Map freq → band and drop out-of-band
	assign_band_inplace(dt)
	dt <- dt[!is.na(band)]

	# Aggregate mean powerdb over time within epoch
	keep_cov <- c("group","condition","rep","complexity","avg_velocity","error",
				  "vresp","accuracy_rating","lat","long")
	keep_cov <- intersect(keep_cov, names(dt))

	agg <- dt[, .(powerdb = mean(powerdb, na.rm = TRUE)),
			  by = c("participant","trial","epoch","chan","band", keep_cov)]

	agg[, band := factor(band, levels = c("theta","alpha","beta"))]
	agg[, rep  := norm_rep(rep)]

	out_shard <- file.path(agg_dir, sprintf("%03d_agg.rds", id_num))
	saveRDS(agg, out_shard)
	invisible(out_shard)
}

invisible(lapply(pp_files, function(fp) { process_one_agg(fp); gc() }))

# Stack -> final compact
agg_files <- list.files(agg_dir, pattern = "_agg\\.rds$", full.names = TRUE)
if (!length(agg_files)) stop("No aggregated shard files found in ", agg_dir)

dat_bands <- rbindlist(lapply(agg_files, readRDS), use.names = TRUE, fill = TRUE)
dat_bands[, rep := norm_rep(rep)]
if (anyNA(dat_bands$rep)) dat_bands <- dat_bands[!is.na(rep)]

setcolorder(dat_bands, intersect(
	c("participant","group","trial","condition","rep","epoch","chan","band",
	  "powerdb","complexity","avg_velocity","error","vresp","accuracy_rating",
	  "lat","long"),
	names(dat_bands)
))

saveRDS(dat_bands, out_bands_path)
message("Wrote compact file: ", out_bands_path, " (rows = ", nrow(dat_bands), ")")
