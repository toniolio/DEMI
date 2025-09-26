######################
### 06_eeg_plots.R ###
######################
# Plots (per participant + group mean) the full time-frequency map at every sensor.
# It's a big, busy plot for visual inspection of all the EEG data.

suppressPackageStartupMessages({
	library(tidyverse)
	library(data.table)
	library(ggplot2)
	library(grid)
	library(readr)
	library(purrr)
})

options(clustermq.scheduler = "multicore")

# -------------------------
# User options
# -------------------------
use_cached_means <- TRUE  # reuse _Scripts/_rds/participants_mean/*.rds if present

# Channel-grid spread (moves sensor centers farther apart/closer together)
#   ↑ Increase to spread centers out; decrease to bring them closer.
scalp_scale_x <- 10.5
scalp_scale_y <- 14.5

# Panel footprint (how big each time×freq tile is around its center)
#   ↓ Smaller values = smaller tiles → more whitespace between panels.
panel_scale_x <- 1.00
panel_scale_y <- 1.00

# Channel label styling/offset (relative to panel size)
label_size_mm <- 3.0       # bigger labels
label_dx_frac <- -0.08     # left offset (fraction of panel width; negative = left)
label_dy_frac <-  0.08     # up   offset (fraction of panel height)

# Verbose progress logging
verbose <- TRUE
vmsg <- function(...) if (isTRUE(verbose)) { message(sprintf(...)); flush.console() }

# -------------
# Helpers
# -------------
scale_to_0range <- function(x, range = 1){
	x <- x - min(x, na.rm = TRUE)
	x <- x / max(x, na.rm = TRUE)
	(x - .5) * range
}

scale_to_bounds <- function(x, lo, hi, range = 1){
	if (!is.finite(lo) || !is.finite(hi) || lo >= hi) stop("Bad bounds for scaling")
	xnorm <- (x - lo) / (hi - lo)   # 0..1
	(xnorm - .5) * range
}

# lightweight colorbar grob with correct limits/breaks/palette
make_colorbar_grob <- function(lims){
	stopifnot(length(lims) == 2, is.finite(lims[1]), is.finite(lims[2]))
	df <- tibble(x = seq(0, 1, length.out = 256),
				 y = 0,
				 z = seq(lims[1], lims[2], length.out = 256))
	p_legend <- ggplot(df, aes(x, y, fill = z)) +
		geom_raster() +
		scale_fill_viridis_c(
			option = 'C',
			limits = lims,
			breaks = scales::pretty_breaks(5)
		) +
		guides(fill = guide_colorbar(
			title = NULL,
			direction = "horizontal",
			barheight = grid::unit(3, "mm"),
			barwidth  = grid::unit(35, "mm"),
			label.position = "bottom"
		)) +
		theme_void() +
		theme(legend.position = "bottom")

	gt <- ggplotGrob(p_legend)
	gt$grobs[[which(sapply(gt$grobs, function(g) g$name) == "guide-box")[1]]]
}

# Returns per-channel raster grob and its bbox (legend disabled for PDF stability)
get_raster_for_channel <- function(this_channel_data, all_channels_powerdb_range){
	if (nrow(this_channel_data) == 0L || all(is.na(this_channel_data$powerdb))) {
		stop("Channel panel has no drawable data (all NA or empty).")
	}
	chan_name <- as.character(this_channel_data$chan[1])

	temp <- this_channel_data %>%
		ggplot() +
		geom_raster(aes(x = to_plot_x, y = to_plot_y, fill = powerdb)) +
		scale_fill_viridis_c(option = 'C',
							 limits = all_channels_powerdb_range,
							 na.value = 'transparent') +
		theme(legend.position = "none")

	layers <- layer_grob(temp)
	topo <- tryCatch({
		if (!is.null(layers$`1`$raster)) layers$`1`$raster else layers[[1]]$raster
	}, error = function(e) NULL)
	if (is.null(topo)) stop("Failed to extract raster layer for a channel panel.")

	xmin <- min(this_channel_data$to_plot_x, na.rm = TRUE)
	xmax <- max(this_channel_data$to_plot_x, na.rm = TRUE)
	ymin <- min(this_channel_data$to_plot_y, na.rm = TRUE)
	ymax <- max(this_channel_data$to_plot_y, na.rm = TRUE)

	list(
		chan        = chan_name,
		topo_raster = topo,
		xmin = xmin, xmax = xmax,
		ymin = ymin, ymax = ymax
	)
}

# Builds one participant page (including the “key” tile, labels, and colorbar)
plot_participant <- function(this_participant_data, all_participants_powerdb_range,
							 X_ABS, Y_MIN, Y_MAX){

	# palette range: group range for individuals; within-page range for mean (0)
	palette_lims <-
		if (this_participant_data$participant[1] == 0) {
			range(this_participant_data$powerdb, na.rm = TRUE)
		} else {
			all_participants_powerdb_range
		}

	# panel layout coords (centers fixed by canonical head coords; tile footprint via panel_scale_*)
	dat <- this_participant_data %>%
		mutate(
			x_scaled    = scale_to_bounds(x, -X_ABS,  X_ABS,  range = scalp_scale_x),
			y_scaled    = scale_to_bounds(y,  Y_MIN,  Y_MAX,  range = scalp_scale_y),
			time_scaled = panel_scale_x * scale_to_0range(time_plot, range = 1),
			freq_scaled = panel_scale_y * scale_to_0range(freq,      range = 1),
			to_plot_x   = x_scaled + time_scaled,
			to_plot_y   = y_scaled + freq_scaled
		)

	# NA stripes (band/time boundaries)
	unique_freqs    <- unique(dat$freq)
	band_boundaries <- c(
		unique_freqs[which.min(abs(unique_freqs -  3.5))],
		unique_freqs[which.min(abs(unique_freqs -  8.5))],
		unique_freqs[which.min(abs(unique_freqs - 12.5))],
		unique_freqs[which.min(abs(unique_freqs - 30.5))]
	)
	time_boundaries <- c(0, 1.5, 2)

	dat <- dat %>%
		mutate(
			powerdb = if_else(!(freq %in% band_boundaries), powerdb, NA_real_),
			powerdb = if_else(!(time_plot %in% time_boundaries), powerdb, NA_real_)
		)

	# per-channel rasters (full page, both epochs concatenated in time_plot)
	channel_rasters <- dat %>%
		group_by(chan) %>% group_split() %>%
		purrr::map(
			.f = get_raster_for_channel,
			all_channels_powerdb_range = palette_lims
		)
	if (!length(channel_rasters)) {
		stop("No drawable rasters for participant ",
			 dat$participant[1],
			 " (group: ", paste(unique(dat$group), collapse = ","), ")")
	}

	# median panel bbox (in plot coords) — use as our unit for the key & axis geometry
	w_med <- median(vapply(channel_rasters, function(cr) cr$xmax - cr$xmin, numeric(1)), na.rm = TRUE)
	h_med <- median(vapply(channel_rasters, function(cr) cr$ymax - cr$ymin, numeric(1)), na.rm = TRUE)

	# canonical canvas limits
	x_scaled_master <- scale_to_bounds(c(-X_ABS, X_ABS), -X_ABS, X_ABS, range = scalp_scale_x)
	y_scaled_master <- scale_to_bounds(c(Y_MIN, Y_MAX),  Y_MIN,  Y_MAX, range = scalp_scale_y)
	x_lim <- range(x_scaled_master) + c(-.5, +.5)
	y_lim <- range(y_scaled_master) + c(-.5, +.5)

	# --- Anchor the key/axes near the lower-left corner using panel-sized padding ---
	x_pad_units <- 1.2 * w_med
	y_pad_units <- 1.5 * h_med

	# y-axis anchor is the LEFT edge of the key; y-axis Y value is the *center* of the key
	y_axis_x_offset <- x_lim[1] + x_pad_units
	y_axis_y_offset <- y_lim[1] + y_pad_units + 0.5*h_med

	# x-axis baseline is the BOTTOM edge of the key
	x_axis_x_offset <- y_axis_x_offset + 0.5*w_med
	x_axis_y_offset <- y_axis_y_offset - 0.5*h_med

	p <- ggplot() +
		scale_x_continuous(limits = x_lim) +
		scale_y_continuous(limits = y_lim) +
		coord_equal() +
		labs(title = dat$participant[1]) +
		theme(
			axis.title = element_blank(),
			axis.ticks = element_blank(),
			axis.text  = element_blank(),
			panel.grid = element_blank(),
			panel.background = element_rect(fill='transparent', colour='grey90'),
			plot.margin = margin(10, 30, 10, 30)
		)

	# draw rasters
	for (cr in channel_rasters){
		p <- p + annotation_raster(
			raster = cr$topo_raster,
			xmin = cr$xmin, xmax = cr$xmax,
			ymin = cr$ymin, ymax = cr$ymax
		)
	}

	# channel labels (anchored to each panel bbox)
	labels_df <- tibble(
		chan = vapply(channel_rasters, function(cr) cr$chan, character(1)),
		xmin = vapply(channel_rasters, function(cr) cr$xmin, numeric(1)),
		xmax = vapply(channel_rasters, function(cr) cr$xmax, numeric(1)),
		ymin = vapply(channel_rasters, function(cr) cr$ymin, numeric(1)),
		ymax = vapply(channel_rasters, function(cr) cr$ymax, numeric(1))
	) %>%
		mutate(
			w = xmax - xmin, h = ymax - ymin,
			label_x = xmin + label_dx_frac * w,
			label_y = ymax + label_dy_frac * h
		)

	p <- p +
		geom_text(
			data = labels_df,
			aes(x = label_x, y = label_y, label = chan),
			color = "black", size = label_size_mm,
			hjust = 1, vjust = 1
		)

	# ===== Axis guides (all sizes proportional to panel w/h so the key axes align) =====

	# x “tick” locations come from NA stripes
	x_tick_locs <- dat %>%
		group_by(time_plot) %>% filter(all(is.na(powerdb))) %>%
		pull(time_scaled) %>% unique() %>% sort()
	if (length(x_tick_locs) < 3) stop("Expected three x tick locations (epoch boundaries).")

	# y “tick” locations come from NA stripes
	y_tick_locs <- dat %>%
		group_by(freq) %>% filter(all(is.na(powerdb))) %>%
		pull(freq_scaled) %>% unique() %>% sort()
	y_tick_locs <- c(-.5, y_tick_locs, .5)
	y_label_locs <- y_tick_locs[1:(length(y_tick_locs)-1)] + diff(y_tick_locs)/2

	# vertical axis line at left (across full key height)
	y_axis_dat_ticks <- tibble(
		to_plot_x = rep(y_axis_x_offset, length(y_tick_locs)),
		to_plot_y = y_tick_locs * h_med + y_axis_y_offset  # scale panel [-.5..+.5] into key height
	)

	# y-axis tick segments (short horizontals)
	y_ticks_seg <- tibble(
		x    = y_axis_x_offset - 0.05 * w_med,
		xend = y_axis_x_offset,
		y    = y_tick_locs * h_med + y_axis_y_offset,
		yend = y_tick_locs * h_med + y_axis_y_offset
	)

	# y-axis labels
	y_axis_dat_labels <- tibble(
		label = c('delta','theta','alpha','beta','gamma'),
		to_plot_x = y_axis_x_offset - 0.10 * w_med,
		to_plot_y = y_label_locs * h_med + y_axis_y_offset
	)

	# x-axis line (bottom of key)
	x_axis_dat_big_ticks <- tibble(
		to_plot_x = c(-.5, x_tick_locs[2], .5) * w_med + x_axis_x_offset,
		to_plot_y = rep(x_axis_y_offset, 3)
	)

	# x-axis tick segments (verticals)
	x_big_ticks_seg <- tibble(
		x    = c(-.5, x_tick_locs[2], .5) * w_med + x_axis_x_offset,
		xend = c(-.5, x_tick_locs[2], .5) * w_med + x_axis_x_offset,
		y    = x_axis_y_offset - 0.10 * h_med,
		yend = x_axis_y_offset
	)

	x_small_ticks_seg <- tibble(
		x    = x_tick_locs[c(1,3)] * w_med + x_axis_x_offset,
		xend = x_tick_locs[c(1,3)] * w_med + x_axis_x_offset,
		y    = x_axis_y_offset - 0.05 * h_med,
		yend = x_axis_y_offset
	)

	# x-axis small labels (“0”)
	x_axis_dat_small_ticks <- tibble(
		label = c('0','0'),
		to_plot_x = x_tick_locs[c(1,3)] * w_med + x_axis_x_offset,
		to_plot_y = rep(x_axis_y_offset, 2)
	)

	# x-axis big labels (“During/After”)
	x_axis_dat_big_labels <- tibble(
		label = c('During','After'),
		to_plot_x = c(-.25, .25) * w_med + x_axis_x_offset,
		to_plot_y = rep(x_axis_y_offset - 0.05 * h_med, 2)
	)

	# axis titles (drop “Time” slightly more)
	axis_title_dat <- tibble(
		label = c('Band','Time'),
		x = c(y_axis_x_offset - 0.25 * w_med, x_axis_x_offset),
		y = c(y_axis_y_offset,               x_axis_y_offset - 0.40 * h_med),
		angle = c(90, 0),
		hjust = c('center','center'),
		vjust = c('bottom','top')
	)

	# draw guides
	p <- p +
		# y-axis line
		geom_line(data = y_axis_dat_ticks, aes(x = to_plot_x, y = to_plot_y), linewidth = .25) +
		# y-axis ticks (segments)
		geom_segment(data = y_ticks_seg,
					 aes(x = x, y = y, xend = xend, yend = yend),
					 linewidth = .25) +
		# y-axis labels
		geom_text(data = y_axis_dat_labels,
				  aes(x = to_plot_x, y = to_plot_y, label = label),
				  hjust = 'right', size = 1, parse = TRUE) +
		# x-axis baseline
		geom_line(data = x_axis_dat_big_ticks, aes(x = to_plot_x, y = to_plot_y), linewidth = .25) +
		# x-axis tick segments
		geom_segment(data = x_big_ticks_seg,
					 aes(x = x, y = y, xend = xend, yend = yend),
					 linewidth = .25) +
		geom_segment(data = x_small_ticks_seg,
					 aes(x = x, y = y, xend = xend, yend = yend),
					 linewidth = .25) +
		# x-axis labels
		geom_text(data = x_axis_dat_small_ticks,
				  aes(x = to_plot_x, y = to_plot_y - 0.08 * h_med, label = label),
				  vjust = 'top', size = 1) +
		geom_text(data = x_axis_dat_big_labels,
				  aes(x = to_plot_x, y = to_plot_y - 0.15 * h_med, label = label),
				  vjust = 'top', size = 2, parse = TRUE) +
		# axis titles
		geom_text(data = axis_title_dat,
				  aes(x = x, y = y, label = label, angle = angle, hjust = hjust, vjust = vjust),
				  size = 3)

	# ===== Example “key” tile — BOTH epochs; size matches a true channel panel =====
	key_chan <- if ("C3" %in% dat$chan) "C3" else unique(dat$chan)[1]
	dat_key  <- dat %>% filter(chan == key_chan)
	if (nrow(dat_key)) {
		key_panel <- get_raster_for_channel(dat_key, all_channels_powerdb_range = palette_lims)

		p <- p + annotation_raster(
			raster = key_panel$topo_raster,
			xmin = y_axis_x_offset,
			xmax = y_axis_x_offset + w_med,
			ymin = y_axis_y_offset - 0.5 * h_med,
			ymax = y_axis_y_offset + 0.5 * h_med
		)
	}

	# ===== Bottom-right colorbar that matches this page’s palette limits (always on) =====
	leg_grob <- make_colorbar_grob(palette_lims)
	leg_w <- 0.20 * diff(x_lim)
	leg_h <- 0.22 * diff(y_lim)
	p <- p + annotation_custom(
		leg_grob,
		xmin = x_lim[2] - leg_w, xmax = x_lim[2],
		ymin = y_lim[1],         ymax = y_lim[1] + leg_h
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

# Behavior, groups, handedness
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

groups_all <- sort(unique(bdat$group))
if (!length(groups_all)) stop("No groups found in bdat2.rds")

# Canonical electrode coords → fixed layout bounds
coords_df <- readr::read_csv(path_coords_csv, col_types = cols()) %>%
	mutate(
		lat  = 90 - (asin(z) * (180 / pi)),
		long = atan2(y, x) * (180 / pi),
		x_raw = lat * cos(long * (pi/180)),
		y_raw = lat * sin(long * (pi/180))
	) %>%
	select(chan, lat, long, x_raw, y_raw)

X_ABS <- max(abs(coords_df$x_raw), na.rm = TRUE)
Y_MIN <- min(coords_df$y_raw,      na.rm = TRUE)
Y_MAX <- max(coords_df$y_raw,      na.rm = TRUE)

# ---------------------------------------
# PASS 1: build or reuse per-participant mean TF (by epoch)
# ---------------------------------------
pp_files <- list.files(path_participants, pattern = "_eeg_processed\\.rds$", full.names = TRUE)
if (!length(pp_files)) stop("No participant shards found in ", path_participants)

group_ranges     <- list()
mean_file_groups <- list()

if (use_cached_means) {
	mean_files <- list.files(path_mean_shards, pattern = "_mean\\.rds$", full.names = TRUE)
	if (!length(mean_files)) {
		stop("use_cached_means=TRUE but no cached means found in ", path_mean_shards,
			 ". Set use_cached_means <- FALSE to rebuild.")
	}
	for (mf in mean_files) {
		dt_p <- readRDS(mf)
		if (!all(c("group","chan","time","freq","epoch_num","powerdb") %in% names(dt_p))) {
			stop("Cached mean shard missing required columns: ", basename(mf),
				 ". Delete cached means and rebuild with use_cached_means <- FALSE.")
		}
		g_here <- unique(dt_p$group)
		mean_file_groups[[mf]] <- g_here
		for (g in g_here) {
			rng <- range(dt_p$powerdb[dt_p$group == g], na.rm = TRUE)
			if (!all(is.finite(rng))) next
			if (is.null(group_ranges[[g]])) group_ranges[[g]] <- rng
			else {
				group_ranges[[g]][1] <- min(group_ranges[[g]][1], rng[1], na.rm = TRUE)
				group_ranges[[g]][2] <- max(group_ranges[[g]][2], rng[2], na.rm = TRUE)
			}
		}
		rm(dt_p); gc()
	}
} else {
	for (f in pp_files) {
		id_num <- as.integer(gsub("\\D", "", basename(f)))
		vmsg("PASS1 (rebuild): %s", basename(f))
		dt <- readRDS(f) %>% as_tibble()

		need <- c("trial","epoch","chan","freq","time","powerdb")
		if (!all(need %in% names(dt))) stop("Missing columns in ", f, ": ", paste(setdiff(need, names(dt)), collapse=", "))

		dt <- dt %>%
			mutate(participant = as.character(id_num)) %>%
			left_join(bdat_trials, by = c("participant","trial")) %>%
			filter(!(group == 'imagery' & condition == 'physical'),
				   !(chan %in% c('M1','M2')))

		if (!nrow(dt)) next

		dt <- dt %>% left_join(coords_df, by = "chan")
		miss_ch <- dt %>% filter(is.na(lat) | is.na(long)) %>% distinct(chan) %>% pull(chan)
		if (length(miss_ch)) {
			stop("Missing coordinates for channels in ", f, ": ",
				 paste(unique(miss_ch), collapse=", "),
				 ". Fix BESA-81.csv or channel naming before plotting.")
		}

		dt <- dt %>%
			left_join(handedness_map, by = "participant") %>%
			mutate(
				x = if_else(handedness == 'l', -x_raw, x_raw),
				y = y_raw
			) %>%
			mutate(epoch_num = if_else(epoch == "during", 1L,
									   if_else(epoch == "after",  2L, NA_integer_))) %>%
			{ if (any(is.na(.$epoch_num))) stop("Non-canonical epoch present in 06 input."); . }

		dt_p <- dt %>%
			group_by(group, participant, chan, x, y, time, freq, epoch_num) %>%
			summarize(powerdb = mean(powerdb, na.rm = TRUE), .groups = 'drop')

		g_here <- unique(dt_p$group)
		if (!length(g_here)) {
			stop("Participant ", id_num, " contributes to no groups after filtering. Check 'group'/'condition' in bdat2.rds.")
		}
		for (g in g_here) {
			rng <- range(dt_p$powerdb[dt_p$group == g], na.rm = TRUE)
			if (!all(is.finite(rng))) next
			if (is.null(group_ranges[[g]])) group_ranges[[g]] <- rng
			else {
				group_ranges[[g]][1] <- min(group_ranges[[g]][1], rng[1], na.rm = TRUE)
				group_ranges[[g]][2] <- max(group_ranges[[g]][2], rng[2], na.rm = TRUE)
			}
		}

		out_f <- file.path(path_mean_shards, sprintf("%03d_mean.rds", id_num))
		saveRDS(dt_p, out_f)
		mean_file_groups[[out_f]] <- g_here

		rm(dt, dt_p); gc()
	}
	mean_files <- list.files(path_mean_shards, pattern = "_mean\\.rds$", full.names = TRUE)
}

for (g in groups_all) {
	has_any <- if (use_cached_means && length(mean_file_groups) == 0) {
		mean_files <- list.files(path_mean_shards, pattern = "_mean\\.rds$", full.names = TRUE)
		any(vapply(mean_files, function(mf){ any(readRDS(mf)$group == g) }, logical(1)))
	} else {
		any(vapply(mean_file_groups, function(gs) any(gs == g), logical(1)))
	}
	if (!has_any) stop("No participants found for expected group '", g, "'.")
}

# ---------------------------------------
# PASS 2: per group → one PDF of all participants + mean page; and a mean PNG
# ---------------------------------------
if (!exists("mean_files") || !length(mean_files)) {
	mean_files <- list.files(path_mean_shards, pattern = "_mean\\.rds$", full.names = TRUE)
}
files_by_group <- lapply(groups_all, function(g) {
	Filter(function(mf) {
		dt_p <- readRDS(mf)
		any(dt_p$group == g)
	}, mean_files)
})
names(files_by_group) <- groups_all

for (g in groups_all) {
	group_dir <- file.path(path_fig_root, g)
	if (!dir.exists(group_dir)) dir.create(group_dir, recursive = TRUE)

	candidate_files <- files_by_group[[g]]
	if (!length(candidate_files)) stop("Internal error: no files for group '", g, "'.")

	pdf_file <- file.path(group_dir, "group.pdf")
	vmsg("Opening PDF device: %s", pdf_file)
	grDevices::pdf(file = pdf_file, width = 8, height = 8, onefile = TRUE, useDingbats = FALSE)

	tryCatch({
		# Participant pages
		n <- length(candidate_files)
		for (i in seq_along(candidate_files)) {
			mf <- candidate_files[[i]]
			id_str <- gsub("_mean\\.rds$", "", basename(mf))
			dt_p <- readRDS(mf) %>% as_tibble() %>% filter(group == g)
			if (!nrow(dt_p)) stop("Empty participant mean shard for group '", g, "': ", mf)

			dat <- dt_p %>% mutate(time_plot = time + 2 * (epoch_num - 1L))
			vmsg("Group %s — participant %s (%d/%d): plotting page...", g, id_str, i, n)
			p_page <- plot_participant(dat,
									   all_participants_powerdb_range = group_ranges[[g]],
									   X_ABS = X_ABS, Y_MIN = Y_MIN, Y_MAX = Y_MAX)
			print(p_page)
			rm(dt_p, dat, p_page); gc()
		}

		# Group mean (participant 0) using canonical coords
		vmsg("Group %s — computing & plotting group mean...", g)
		dt_g <- rbindlist(lapply(candidate_files, function(mf) {
			readRDS(mf) %>% as_tibble() %>% filter(group == g)
		}), use.names = TRUE, fill = TRUE)
		if (!nrow(dt_g)) stop("No rows to compute group mean for '", g, "'.")

		dt_g <- dt_g %>%
			group_by(group, chan, time, freq, epoch_num) %>%
			summarize(powerdb = mean(powerdb, na.rm = TRUE), .groups = 'drop') %>%
			left_join(coords_df, by = "chan") %>%
			mutate(
				x = x_raw,
				y = y_raw,
				participant = 0
			)

		dat_mean <- dt_g %>% mutate(time_plot = time + 2 * (epoch_num - 1L))
		p_mean <- plot_participant(dat_mean,
								   all_participants_powerdb_range = group_ranges[[g]],
								   X_ABS = X_ABS, Y_MIN = Y_MIN, Y_MAX = Y_MAX)
		print(p_mean)

		# Also write PNG of the mean page (separate device)
		vmsg("Group %s — saving group mean PNG...", g)
		ggsave(filename = file.path(group_dir, "group_mean.png"),
			   plot = p_mean, width = 8, height = 8, dpi = 300)

		rm(dt_g, dat_mean, p_mean); gc()
	},
	finally = {
		dev.off()
		ok <- file.exists(pdf_file)
		sz <- if (ok) file.info(pdf_file)$size else 0L
		if (!ok || is.na(sz) || sz <= 0L) {
			stop("Failed to write non-empty PDF: ", pdf_file)
		} else {
			message("Wrote: ", pdf_file, " (", format(sz, big.mark=","), " bytes)")
		}
	})
}
