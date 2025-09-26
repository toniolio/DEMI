#####################
### 07_gam_prep.R ###
#####################
# Authors: Tony Ingram & Mike Lawrence
# Purpose: Prepare tidy trial-level data for GAMs from the
#          pre-aggregated file dat_bands_notime.rds:
#          - already band- & epoch-averaged power (powerdb)
#          - left-hander mirroring of coordinates
#          - scaled accuracy_rating
#          - factor levels & labels
#          - plot a PMBR visual for sanity check
# Output : _Scripts/_rds/dat_gam.rds

suppressPackageStartupMessages({
	library(tidyverse)
})

# -------------------- User parameters --------------------
path_rds  <- "_Scripts/_rds/"
path_misc <- "_Scripts/_misc/"

# Channels to drop (e.g., references)
drop_chans <- c("M1","M2")

# Round coordinates to nearest integer grid (as in previous code)
round_coords <- TRUE
# ---------------------------------------------------------

# Helper: safe readRDS
read_rds_safe <- function(p) {
	if (!file.exists(p)) stop("Missing file: ", p)
	readRDS(p)
}

# 1) Load the pre-aggregated dataset (trial × epoch × band × sensor)
# Expected columns include:
# participant, trial, epoch, chan, band, powerdb,
# group, condition, rep, complexity, avg_velocity, error, vresp,
# accuracy_rating, lat, long
bands_path <- file.path(path_rds, "dat_bands_notime.rds")
dat <- read_rds_safe(bands_path) %>%
	as_tibble() %>%
	dplyr::mutate(epoch = factor(as.character(epoch), levels = c("during","after")))

# Normalize 'rep' defensively (in case upstream files predate the fix)
rep_chr <- tolower(as.character(dat$rep))
rep_num <- suppressWarnings(as.numeric(as.character(dat$rep)))
dat <- dat %>%
	mutate(
		rep = case_when(
			!is.na(rep_chr) & rep_chr %in% c("random","repeated") ~ rep_chr,
			!is.na(rep_num) & rep_num == 1                        ~ "random",
			!is.na(rep_num) & rep_num %in% c(0,2)                 ~ "repeated",
			TRUE                                                  ~ NA_character_
		)
	)
# We expect no NA after the 5 fix; warn if found
if (any(is.na(dat$rep))) {
	warning("Found NA in 'rep' after normalization in 07_gam_prep. Consider re-running 5_eeg_prep.R.")
}

# 2) Basic sanity checks (adapted to aggregated file)
req_cols <- c("participant","group","trial","condition","rep","epoch",
			  "chan","lat","long","band","powerdb","accuracy_rating")
missing_cols <- setdiff(req_cols, names(dat))
if (length(missing_cols)) {
	stop("Missing required columns: ", paste(missing_cols, collapse=", "))
}

# 3) Map epoch to labeled factor if coded numerically (1=during, 2=after)
if (is.numeric(dat$epoch)) {
	dat <- dat %>%
		mutate(epoch = factor(if_else(epoch == 1, "during", "after"),
							  levels = c("during","after")))
} else {
	dat <- dat %>%
		mutate(epoch = factor(as.character(epoch), levels = c("during","after")))
}

# 4) Drop bad/reference channels if present
if (!is.null(drop_chans) && length(drop_chans)) {
	dat <- dat %>% filter(!(chan %in% drop_chans))
}

# 5) Scale accuracy rating (z-score)
dat <- dat %>% mutate(accuracy = as.numeric(scale(accuracy_rating)))

# 6) Left-hander mirroring (about the midline) using lat/long → (x,y) → mirror → lat/long
#    Requires a handedness lookup file with columns: participant, handedness ('l','r','a'…)
handedness_path <- file.path(path_rds, "bdat2.rds")
if (file.exists(handedness_path)) {
	handed <- read_rds_safe(handedness_path) %>%
		mutate(participant = as.character(participant)) %>%
		group_by(participant) %>%
		summarize(handedness = first(handedness), .groups = "drop")

	dat <- dat %>%
		mutate(participant = as.character(participant)) %>%
		left_join(handed, by = "participant") %>%
		# convert to planar for mirroring (x,y) on circle of radius=lat
		mutate(
			x = lat * cos(long * (pi/180)),
			y = lat * sin(long * (pi/180)),
			x = if_else(handedness == 'l', -x, x),
			long = atan2(y, x) * (180/pi)
		) %>%
		select(-x, -y, -handedness)
}

# 7) Coordinate housekeeping
if (round_coords) {
	dat <- dat %>%
		mutate(
			lat  = round(lat),
			long = round(long),
			long = if_else(lat == 0, 0, long) # keep poles tidy
		)
}

# 8) Factor levels & labels
halfsum_contrasts <- function (...) contr.sum(...) * 0.5

dat <- dat %>%
	mutate(
		group       = factor(group, levels = c("physical","imagery")),
		# rep may be coded as 0/1 or strings; normalize to levels
		rep         = case_when(
			rep %in% c(1,"1","random")   ~ "random",
			rep %in% c(0,"0","repeated") ~ "repeated",
			TRUE                          ~ as.character(rep)
		),
		rep         = factor(rep, levels = c("random","repeated")),
		band        = factor(band, levels = c("theta","alpha","beta")),
		epoch       = factor(epoch, levels = c("during","after")),
		participant = factor(participant),
		chan        = factor(chan)
	)

# Optional: derive 1..6 block id from trial index if not present
if (!"block" %in% names(dat)) {
	dat <- dat %>%
		mutate(block = case_when(
			trial >=   1 & trial <=  20 ~ 1L,
			trial >=  21 & trial <=  40 ~ 2L,
			trial >=  41 & trial <=  60 ~ 3L,
			trial >=  61 & trial <=  80 ~ 4L,
			trial >=  81 & trial <= 100 ~ 5L,
			trial >= 101 & trial <= 120 ~ 6L,
			TRUE ~ NA_integer_
		))
}

# 9) Save tidy GAM-ready data
out_path <- file.path(path_rds, "dat_gam.rds")
saveRDS(dat, out_path)
message("Wrote: ", out_path)

# 10) Minimal checks
message("N participants: ", n_distinct(dat$participant))
# robust unique trials across participants:
message("N trials (unique pid×trial): ", n_distinct(paste(dat$participant, dat$trial)))
message("N rows: ", nrow(dat))


# ===================== PMBR QC PLOT ======================

# Check if you see canonical PMBR, otherwise your time averaging in script 05 is bad!

# Assumes `dat` is already loaded above and `epoch` is factor(during, after)

# Output location (same root as script 06)
plot_root <- "_Plots"
if (!dir.exists(plot_root)) dir.create(plot_root, recursive = TRUE)

# Tunables to match the original 09 defaults
chan_of_interest <- "C3"
band_of_interest <- "beta"
split_by_group   <- TRUE    # color by group
split_by_rep     <- TRUE    # and by rep (interaction)

# Filter to channel & band
d0 <- dat %>%
	dplyr::filter(.data$chan == chan_of_interest,
				  .data$band == band_of_interest)

if (nrow(d0) == 0L) {
	stop("QC plot: no rows for chan=", chan_of_interest, " & band=", band_of_interest,
		 ". Check channel/band labels or drop list.")
}

# Participant means, then t-based 95% CI across participants
grp <- c("epoch")
if (split_by_group) grp <- c(grp, "group")
if (split_by_rep)   grp <- c(grp, "rep")

d_part <- d0 %>%
	dplyr::group_by(dplyr::across(dplyr::all_of(c("participant", grp)))) %>%
	dplyr::summarise(powerdb = mean(.data$powerdb, na.rm = TRUE), .groups = "drop")

t_mult <- function(n, level = 0.95) stats::qt(1 - (1 - level)/2, df = pmax(n - 1, 1))

d_sum <- d_part %>%
	dplyr::group_by(dplyr::across(dplyr::all_of(grp))) %>%
	dplyr::summarise(
		n    = dplyr::n_distinct(.data$participant),
		mean = mean(.data$powerdb, na.rm = TRUE),
		sd   = stats::sd(.data$powerdb, na.rm = TRUE),
		se   = .data$sd / sqrt(.data$n),
		tcrit= t_mult(.data$n),
		ci_lo= .data$mean - .data$tcrit * .data$se,
		ci_hi= .data$mean + .data$tcrit * .data$se,
		.groups = "drop"
	)

pos <- ggplot2::position_dodge(width = 0.5)

# Build plot (replicates original 09 aesthetics and labeling)
if (split_by_group && split_by_rep) {
	d_sum$col <- interaction(d_sum$group, d_sum$rep, drop = TRUE, sep = " × ")
	p <- ggplot2::ggplot(d_sum, ggplot2::aes(x = .data$epoch, y = .data$mean, color = .data$col)) +
		ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$ci_lo, ymax = .data$ci_hi), width = 0.15, position = pos) +
		ggplot2::geom_point(size = 3, position = pos) +
		ggplot2::labs(
			title    = paste0("PMBR QC at ", chan_of_interest, " (", band_of_interest, ")"),
			subtitle = "Participant means; error bars = 95% CI across participants",
			x = "Epoch", y = "Power (dB)", color = "Group × Rep"
		) +
		ggplot2::theme_minimal(base_size = 13)
} else if (split_by_group) {
	p <- ggplot2::ggplot(d_sum, ggplot2::aes(x = .data$epoch, y = .data$mean, color = .data$group)) +
		ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$ci_lo, ymax = .data$ci_hi), width = 0.15, position = pos) +
		ggplot2::geom_point(size = 3, position = pos) +
		ggplot2::labs(
			title    = paste0("PMBR QC at ", chan_of_interest, " (", band_of_interest, ")"),
			subtitle = "Participant means; error bars = 95% CI across participants",
			x = "Epoch", y = "Power (dB)", color = "Group"
		) +
		ggplot2::theme_minimal(base_size = 13)
} else if (split_by_rep) {
	p <- ggplot2::ggplot(d_sum, ggplot2::aes(x = .data$epoch, y = .data$mean, color = .data$rep)) +
		ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$ci_lo, ymax = .data$ci_hi), width = 0.15, position = pos) +
		ggplot2::geom_point(size = 3, position = pos) +
		ggplot2::labs(
			title    = paste0("PMBR QC at ", chan_of_interest, " (", band_of_interest, ")"),
			subtitle = "Participant means; error bars = 95% CI across participants",
			x = "Epoch", y = "Power (dB)", color = "Rep"
		) +
		ggplot2::theme_minimal(base_size = 13)
} else {
	p <- ggplot2::ggplot(d_sum, ggplot2::aes(x = .data$epoch, y = .data$mean)) +
		ggplot2::geom_errorbar(ggplot2::aes(ymin = .data$ci_lo, ymax = .data$ci_hi), width = 0.15) +
		ggplot2::geom_point(size = 3) +
		ggplot2::labs(
			title    = paste0("PMBR QC at ", chan_of_interest, " (", band_of_interest, ")"),
			subtitle = "Participant means; error bars = 95% CI across participants",
			x = "Epoch", y = "Power (dB)"
		) +
		ggplot2::theme_minimal(base_size = 13)
}

print(p)
outfile <- file.path(plot_root, sprintf("pmbr_%s_%s.png", chan_of_interest, band_of_interest))
ggplot2::ggsave(outfile, p, width = 7.5, height = 4.8, dpi = 300)
message("Wrote QC figure: ", outfile)
