######################
### 06_eeg_plots.R ###
######################
# plot for every participant (as well as a grand mean) the complete
# time-frequency plot at every sensor. It's a gigantic plot, and not
# meant for publication — but for inspection of the data by the user.

suppressPackageStartupMessages({
	library(tidyverse)
	library(data.table)
	library(ggplot2)
	library(grid)
	library(readr)
	library(purrr)
})

options(clustermq.scheduler = "multicore")

# --------
# Helpers
# --------
scale_to_0range <- function(x, range=1){
	x <- x - min(x, na.rm = TRUE)
	x <- x / max(x, na.rm = TRUE)
	x <- x - .5
	x * range
}

# Build one channel raster; returns list with grobs/limits (throws if nothing drawable)
get_raster_for_channel <- function(this_channel_data, all_channels_powerdb_range){
	if (nrow(this_channel_data) == 0L || all(is.na(this_channel_data$powerdb))) {
		stop("Channel panel has no drawable data (all NA or empty).")
	}

	temp <- this_channel_data %>%
		ggplot() +
		geom_raster(aes(x = to_plot_x, y = to_plot_y, fill = powerdb)) +
		scale_fill_viridis_c(option = 'C',
							 limits = all_channels_powerdb_range,
							 na.value = 'transparent') +
		theme(legend.position = "top", legend.title = element_blank())

	gb   <- ggplot_build(temp)
	gt   <- ggplot_gtable(gb)
	lg_ix <- which(sapply(gt$grobs, function(x) x$name) == "guide-box")
	legend_grob <- if (length(lg_ix)) gt$grobs[[lg_ix[1]]] else NULL

	layers <- layer_grob(temp)
	topo <- tryCatch({
		if (!is.null(layers$`1`$raster)) layers$`1`$raster else layers[[1]]$raster
	}, error = function(e) NULL)
	if (is.null(topo)) stop("Failed to extract raster layer for a channel panel.")

	list(
		legend = legend_grob,
		topo_raster = topo,
		xmin = min(this_channel_data$to_plot_x, na.rm = TRUE),
		xmax = max(this_channel_data$to_plot_x, na.rm = TRUE),
		ymin = min(this_channel_data$to_plot_y, na.rm = TRUE),
		ymax = max(this_channel_data$to_plot_y, na.rm = TRUE)
	)
}

# Per-participant page plot (returns ggplot object; caller prints it to the active device)
plot_participant <- function(this_participant_data, all_participants_powerdb_range){
	if (this_participant_data$participant[1] == 0) {
		this_participant_powerdb_range <- range(this_participant_data$powerdb, na.rm = TRUE)
	} else {
		this_participant_powerdb_range <- all_participants_powerdb_range
	}

	channel_rasters <- this_participant_data %>%
		group_by(x, y) %>% group_split() %>%
		purrr::map(
			.f = get_raster_for_channel,
			all_channels_powerdb_range = this_participant_powerdb_range
		)

	if (!length(channel_rasters)) {
		stop("No drawable rasters for participant ",
			 this_participant_data$participant[1],
			 " (group: ", paste(unique(this_participant_data$group), collapse = ","), ")")
	}

	# Canvas
	p <- ggplot() +
		scale_x_continuous(limits = range(this_participant_data$to_plot_x, na.rm = TRUE)) +
		scale_y_continuous(limits = range(this_participant_data$to_plot_y, na.rm = TRUE)) +
		coord_equal() +
		labs(title = this_participant_data$participant[1]) +
		theme(
			axis.title = element_blank(),
			axis.ticks = element_blank(),
			axis.text  = element_blank(),
			panel.grid = element_blank(),
			panel.background = element_rect(fill='transparent', colour='grey90')
		)

	# Add rasters
	for (cr in channel_rasters){
		p <- p + annotation_raster(
			raster = cr$topo_raster,
			xmin = cr$xmin, xmax = cr$xmax,
			ymin = cr$ymin, ymax = cr$ymax
		)
	}

	# Legend
	last_leg <- NULL
	for (cr in channel_rasters) if (!is.null(cr$legend)) last_leg <- cr$legend
	if (!is.null(last_leg)) {
		p <- p +
			annotation_custom(
				last_leg,
				xmin = max(this_participant_data$to_plot_x, na.rm = TRUE) - .2 * diff(range(this_participant_data$to_plot_x, na.rm = TRUE)),
				xmax = max(this_participant_data$to_plot_x, na.rm = TRUE),
				ymin = min(this_participant_data$to_plot_y, na.rm = TRUE) - .1 * diff(range(this_participant_data$to_plot_y, na.rm = TRUE)),
				ymax = min(this_participant_data$to_plot_y, na.rm = TRUE) + .2 * diff(range(this_participant_data$to_plot_y, na.rm = TRUE))
			)
	}

	# Axis guides
	x_tick_locs <- this_participant_data %>%
		group_by(time_plot) %>% filter(all(is.na(powerdb))) %>% pull(time_scaled) %>% unique() %>% sort()

	y_tick_locs <- this_participant_data %>%
		group_by(freq) %>% filter(all(is.na(powerdb))) %>% pull(freq_scaled) %>% unique() %>% sort()
	y_tick_locs <- c(-.5, y_tick_locs, .5)
	y_label_locs <- y_tick_locs[1:(length(y_tick_locs)-1)] + diff(y_tick_locs)/2

	y_axis_y_offset <- -4.5
	y_axis_x_offset <- -5
	y_axis_dat_ticks <- tibble(
		y_scaled = y_tick_locs,
		to_plot_y = y_scaled + y_axis_y_offset,
		to_plot_x = rep(0, length(y_scaled)) + y_axis_x_offset,
		label = 1:length(y_scaled)
	)
	y_axis_dat_labels <- tibble(
		label = c('delta','theta','alpha','beta','gamma'),
		y_scaled = y_label_locs,
		to_plot_y = y_scaled + y_axis_y_offset,
		to_plot_x = rep(0, length(y_scaled)) + y_axis_x_offset
	)

	x_axis_y_offset <- y_axis_y_offset - .5
	x_axis_x_offset <- y_axis_x_offset + .5

	if (length(x_tick_locs) < 3) {
		stop("Expected three x tick locations (epoch boundaries) were not found in data.")
	}

	x_axis_dat_small_ticks <- tibble(
		label = c('0','0'),
		x_scaled = x_tick_locs[c(1,3)],
		to_plot_x = x_scaled + x_axis_x_offset,
		to_plot_y = rep(0, length(label)) + y_axis_y_offset
	)
	x_axis_dat_big_ticks <- tibble(
		x_scaled = c(-.5, x_tick_locs[2], .5),
		to_plot_x = x_scaled + x_axis_x_offset,
		to_plot_y = rep(0, length(x_scaled)) + y_axis_y_offset,
		label = 1:length(x_scaled)
	)
	x_axis_dat_big_labels <- tibble(
		label = c('During','After'),
		x_scaled = c(-.25, .25),
		to_plot_x = x_scaled + x_axis_x_offset,
		to_plot_y = rep(0, length(label)) + y_axis_y_offset - .05
	)

	axis_title_dat <- tibble(
		label = c('Band','Time'),
		x = c(y_axis_x_offset - .25, x_axis_x_offset),
		y = c(y_axis_y_offset,      x_axis_y_offset - .4),
		angle = c(90, 0),
		hjust = c('center','center'),
		vjust = c('bottom','top')
	)

	p <- p +
		# y-axis line
		geom_line(data = y_axis_dat_ticks, aes(x = to_plot_x, y = to_plot_y)) +
		# y-axis ticks
		geom_line(
			data = (
				y_axis_dat_ticks %>%
					mutate(xmin = to_plot_x - .05, xmax = to_plot_x) %>%
					select(-to_plot_x) %>%
					pivot_longer(cols = c(xmin, xmax), values_to = 'to_plot_x')
			),
			aes(x = to_plot_x, y = to_plot_y, group = label),
			size = .25
		) +
		# y-axis labels
		geom_text(
			data = y_axis_dat_labels,
			aes(x = to_plot_x - .1, y = to_plot_y, label = label),
			hjust = 'right', size = 1, parse = TRUE
		) +
		# x-axis line
		geom_line(data = x_axis_dat_big_ticks, aes(x = to_plot_x, y = to_plot_y)) +
		# x-axis big ticks
		geom_line(
			data = (
				x_axis_dat_big_ticks %>%
					mutate(ymin = to_plot_y - .1, ymax = to_plot_y) %>%
					select(-to_plot_y) %>%
					pivot_longer(cols = c(ymin, ymax), values_to = 'to_plot_y')
			),
			aes(x = to_plot_x, y = to_plot_y, group = interaction(label, x_scaled)),
			size = .25
		) +
		# x-axis small ticks
		geom_line(
			data = (
				x_axis_dat_small_ticks %>%
					mutate(ymin = to_plot_y - .05, ymax = to_plot_y) %>%
					select(-to_plot_y) %>%
					pivot_longer(cols = c(ymin, ymax), values_to = 'to_plot_y')
			),
			aes(x = to_plot_x, y = to_plot_y, group = interaction(label, x_scaled)),
			size = .25
		) +
		# x-axis small labels
		geom_text(
			data = x_axis_dat_small_ticks,
			aes(x = to_plot_x, y = to_plot_y - .08, label = label),
			vjust = 'top', size = 1
		) +
		# extra x-axis big labels
		geom_text(
			data = x_axis_dat_big_labels,
			aes(x = to_plot_x, y = to_plot_y - .15, label = label),
			vjust = 'top', size = 2, parse = TRUE
		) +
		# axis titles
		geom_text(
			data = axis_title_dat,
			aes(x = x, y = y, label = label, angle = angle, hjust = hjust, vjust = vjust),
			size = 3
		) +
		# example panel
		annotation_raster(
			raster = channel_rasters[[1]]$topo_raster,
			xmin = y_axis_x_offset, xmax = y_axis_x_offset + 1,
			ymin = y_axis_y_offset - .5, ymax = y_axis_y_offset + .5
		)

	return(p)
}

# ---------------
# Paths & inputs
# ---------------
path_participants <- "_Scripts/_rds/participants/"
path_bdat         <- "_Scripts/_rds/bdat2.rds"
path_coords_csv   <- "_Data/eeg/BESA-81.csv"
path_mean_shards  <- "_Scripts/_rds/participants_mean/"

# everything saves under _Plots/
plot_root     <- "_Plots"
path_fig_root <- file.path(plot_root, "raw_power")

if (!dir.exists(path_participants)) stop("Missing folder: ", path_participants)
if (!file.exists(path_bdat))       stop("Missing: ", path_bdat)
if (!file.exists(path_coords_csv)) stop("Missing: ", path_coords_csv)
if (!dir.exists(path_mean_shards)) dir.create(path_mean_shards, recursive = TRUE)
if (!dir.exists(path_fig_root))    dir.create(path_fig_root,    recursive = TRUE)

# Behaviour, groups, handedness
bdat <- readRDS(path_bdat) %>% mutate(participant = as.character(participant))
handedness_map <- bdat %>%
	group_by(participant) %>%
	summarize(handedness = handedness[1], .groups = "drop")

bdat_trials <- bdat %>%
	transmute(
		participant = as.character(participant),
		trial = as.integer(trial_num + (block_num - 1L) * 20L),
		group, condition, rep
	)

# Groups we expect to plot (from bdat)
groups_all <- sort(unique(bdat$group))
if (!length(groups_all)) stop("No groups found in bdat2.rds")

# Preload electrode coords
coords_df <- readr::read_csv(path_coords_csv, col_types = cols()) %>%
	mutate(
		lat  = 90 - (asin(z) * (180 / pi)),
		long = atan2(y, x) * (180 / pi)
	) %>%
	select(chan, lat, long)

# ---------------------------------------
# PASS 1: stream each participant shard -> per-participant trial-mean TF (preserve epoch)
#         compute per-group color ranges; cache which groups each file has
# ---------------------------------------
pp_files <- list.files(path_participants, pattern = "_eeg_processed\\.rds$", full.names = TRUE)
if (!length(pp_files)) stop("No participant shards found in ", path_participants)

group_ranges <- list()          # per-group min/max
mean_file_groups <- list()      # key: mean file path -> character vector of groups

for (f in pp_files) {
	id_num <- as.integer(gsub("\\D", "", basename(f)))
	message("PASS1: ", basename(f))
	dt <- readRDS(f) %>% as_tibble()

	need <- c("trial","epoch","chan","freq","time","powerdb")
	if (!all(need %in% names(dt))) stop("Missing columns in ", f, ": ", paste(setdiff(need, names(dt)), collapse=", "))

	# Merge behaviour
	dt <- dt %>%
		mutate(participant = as.character(id_num)) %>%
		left_join(bdat_trials,    by = c("participant","trial"))

	# Filter out channels & conditions first
	dt <- dt %>%
		filter(!(group == 'imagery' & condition == 'physical'),
			   !(chan %in% c('M1','M2')))

	# Join coords and enforce strict presence AFTER filtering
	dt <- dt %>% left_join(coords_df, by = "chan")
	miss_ch <- dt %>% filter(is.na(lat) | is.na(long)) %>% distinct(chan) %>% pull(chan)
	if (length(miss_ch)) {
		stop(
			"Missing coordinates for channels in ", f, ": ",
			paste(unique(miss_ch), collapse=", "),
			". Fix BESA-81.csv or channel naming before plotting."
		)
	}

	# Handedness & cartesian
	dt <- dt %>%
		left_join(handedness_map, by = "participant") %>%
		mutate(
			x = lat * cos(long * (pi/180)),
			y = lat * sin(long * (pi/180)),
			x = case_when(handedness == 'l' ~ -x, TRUE ~ x)
		)

	# Epoch -> numeric (1=tracing, 2=post_trace); preserve epoch in summary
	dt <- dt %>%
		mutate(epoch_num = case_when(
			epoch %in% c("tracing", 1L, "1") ~ 1L,
			epoch %in% c("post_trace", 2L, "2") ~ 2L,
			TRUE ~ NA_integer_
		))
	if (any(is.na(dt$epoch_num))) {
		bad_epochs <- dt %>% filter(is.na(epoch_num)) %>% distinct(epoch) %>% pull(epoch)
		stop("Unrecognized epoch values in ", f, ": ", paste(unique(bad_epochs), collapse=", "))
	}

	# per-participant mean across trials *by epoch* (key to keep both halves)
	dt_p <- dt %>%
		group_by(group, participant, chan, x, y, time, freq, epoch_num) %>%
		summarize(powerdb = mean(powerdb, na.rm = TRUE), .groups = 'drop')

	# Determine groups present for this participant
	g_here <- unique(dt_p$group)
	if (!length(g_here)) {
		stop("Participant ", id_num, " contributes to no groups after filtering. Check 'group'/'condition' in bdat2.rds.")
	}

	# Update per-group ranges using per-epoch means
	for (g in g_here) {
		rng <- range(dt_p$powerdb[dt_p$group == g], na.rm = TRUE)
		if (!all(is.finite(rng))) next
		if (is.null(group_ranges[[g]])) group_ranges[[g]] <- rng
		else {
			group_ranges[[g]][1] <- min(group_ranges[[g]][1], rng[1], na.rm = TRUE)
			group_ranges[[g]][2] <- max(group_ranges[[g]][2], rng[2], na.rm = TRUE)
		}
	}

	# Save compact shard and remember which groups it has
	out_f <- file.path(path_mean_shards, sprintf("%03d_mean.rds", id_num))
	saveRDS(dt_p, out_f)
	mean_file_groups[[out_f]] <- g_here

	rm(dt, dt_p); gc()
}

# Assert that every expected group has at least one participant shard
for (g in groups_all) {
	has_any <- any(vapply(mean_file_groups, function(gs) any(gs == g), logical(1)))
	if (!has_any) {
		stop("No participants found for expected group '", g, "'. Check earlier steps (bdat2.rds, filtering rules).")
	}
}

# ---------------------------------------
# PASS 2: for each group, open a single PDF, draw all participant pages
#         (epochs concatenated along time), then draw the group-mean page.
#         Also save mean as PNG. Hard-fail on any device error.
# ---------------------------------------
for (g in groups_all) {
	group_dir <- file.path(path_fig_root, g)
	if (!dir.exists(group_dir)) dir.create(group_dir, recursive = TRUE)

	candidate_files <- names(Filter(function(gs) any(gs == g), mean_file_groups))
	if (!length(candidate_files)) {
		stop("Internal error: group '", g, "' had no candidate files after PASS 1 checks.")
	}

	# Single PDF per group
	pdf_file <- file.path(group_dir, "group.pdf")
	pdf(file = pdf_file, height = 8, width = 8, onefile = TRUE, useDingbats = FALSE)

	# Participant pages
	for (mf in candidate_files) {
		dt_p <- readRDS(mf) %>% as_tibble() %>% filter(group == g)
		if (!nrow(dt_p)) stop("Empty participant mean shard for group '", g, "': ", mf)

		# Concatenate epochs side-by-side along time
		dat <- dt_p %>%
			mutate(
				time_plot  = time + 2 * (epoch_num - 1L),  # epoch 1 at [..], epoch 2 shifted by +2s
				x_scaled    = scale_to_0range(x, 8),
				y_scaled    = scale_to_0range(y, 8.5),
				time_scaled = scale_to_0range(time_plot, 1),
				freq_scaled = scale_to_0range(freq, 1),
				to_plot_x   = x_scaled + time_scaled,
				to_plot_y   = y_scaled + freq_scaled
			)

		# NA stripes at band/time boundaries
		unique_freqs    <- unique(dat$freq)
		band_boundaries <- c(
			unique_freqs[which.min(abs(unique_freqs -  3.5))],
			unique_freqs[which.min(abs(unique_freqs -  8.5))],
			unique_freqs[which.min(abs(unique_freqs - 12.5))],
			unique_freqs[which.min(abs(unique_freqs - 30.5))]
		)
		time_boundaries <- c(0, 1.5, 2)  # 0 & 1.5 within first epoch; 2 splits epochs

		dat <- dat %>%
			mutate(
				powerdb = if_else(!(freq %in% band_boundaries), powerdb, NA_real_),
				powerdb = if_else(!(time_plot %in% time_boundaries), powerdb, NA_real_)
			)

		p_page <- plot_participant(dat, all_participants_powerdb_range = group_ranges[[g]])
		tryCatch({
			print(p_page)
		}, error = function(e) {
			dev.off()
			stop("Failed to write participant page for group '", g, "' (", basename(mf), "): ", conditionMessage(e))
		})
		rm(dt_p, dat, p_page); gc()
	}

	# Group-mean page (participant 0)
	dt_g <- rbindlist(lapply(candidate_files, function(mf) {
		readRDS(mf) %>% as_tibble() %>% filter(group == g)
	}), use.names = TRUE, fill = TRUE)

	if (!nrow(dt_g)) stop("No rows to compute group mean for '", g, "'.")

	dt_g <- dt_g %>%
		group_by(group, chan, x, y, time, freq, epoch_num) %>%
		summarize(powerdb = mean(powerdb, na.rm = TRUE), .groups = 'drop') %>%
		mutate(participant = 0)

	dat_mean <- dt_g %>%
		mutate(
			time_plot  = time + 2 * (epoch_num - 1L),
			x_scaled    = scale_to_0range(x, 8),
			y_scaled    = scale_to_0range(y, 8.5),
			time_scaled = scale_to_0range(time_plot, 1),
			freq_scaled = scale_to_0range(freq, 1),
			to_plot_x   = x_scaled + time_scaled,
			to_plot_y   = y_scaled + freq_scaled
		)

	unique_freqs    <- unique(dat_mean$freq)
	band_boundaries <- c(
		unique_freqs[which.min(abs(unique_freqs -  3.5))],
		unique_freqs[which.min(abs(unique_freqs -  8.5))],
		unique_freqs[which.min(abs(unique_freqs - 12.5))],
		unique_freqs[which.min(abs(unique_freqs - 30.5))]
	)
	time_boundaries <- c(0, 1.5, 2)

	dat_mean <- dat_mean %>%
		mutate(
			powerdb = if_else(!(freq %in% band_boundaries), powerdb, NA_real_),
			powerdb = if_else(!(time_plot %in% time_boundaries), powerdb, NA_real_)
		)

	p_mean <- plot_participant(dat_mean, all_participants_powerdb_range = group_ranges[[g]])
	tryCatch({
		print(p_mean)  # last page in the PDF
	}, error = function(e) {
		dev.off()
		stop("Failed to write group-mean page for group '", g, "': ", conditionMessage(e))
	})

	# Close the PDF before writing PNG
	dev.off()

	# Also write a PNG of the mean page (separate device)
	tryCatch({
		ggsave(filename = file.path(group_dir, "group_mean.png"),
			   plot = p_mean, width = 8, height = 8, dpi = 300)
	}, error = function(e) {
		stop("Failed to write PNG for group '", g, "': ", conditionMessage(e))
	})

	rm(dt_g, dat_mean, p_mean); gc()
}
