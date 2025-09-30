###############################################
### 10_posthoc_group_maps.R (GP + acc_bin) ###
###############################################
suppressPackageStartupMessages({
	library(tidyverse)
	library(mgcv)
	library(readr)
	library(jsonlite)
})

message("\n=== 10_posthoc_group_maps.R (descriptive per-group maps) ===")

# -------- Paths & options --------
path_rds   <- "_Scripts/_rds/"
path_misc  <- "_Scripts/_misc/"
out_plots  <- "_Plots/"
out_tables <- "_Tables/"
dir.create(out_plots,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_tables, recursive = TRUE, showWarnings = FALSE)

# Plot options
t_cutoff        <- 2.0      # per-sensor |t| threshold for a black ring (descriptive only)
add_chan_labels <- TRUE    # set TRUE to draw channel codes next to points
common_scale    <- TRUE    # TRUE => one symmetric color scale across all 6 maps

# -------- Load model & data --------
gam_full <- readRDS(file.path(path_rds, "gam_full.rds"))
dat_gam  <- readRDS(file.path(path_rds, "dat_gam.rds")) %>% as_tibble()

# Factor housekeeping
dat_gam <- dat_gam %>%
	mutate(
		group = factor(group, levels = c("physical","imagery")),
		band  = factor(band,  levels = c("theta","alpha","beta")),
		epoch = factor(epoch, levels = c("during","after")),
		rep   = factor(rep,   levels = c("random","repeated")),
		chan  = factor(chan),
		participant = factor(participant)
	) %>% droplevels()

# Channel map (lat/long) used for plotting
map_file <- file.path(path_misc, "sensor_latlong_chan_map.csv")
if (file.exists(map_file)) {
	latlong_map <- read_csv(map_file, show_col_types = FALSE) %>% distinct(chan, lat, long)
} else {
	latlong_map <- dat_gam %>%
		group_by(chan) %>% summarise(lat = first(lat), long = first(long), .groups = "drop")
}

# ----- Build x,y used by the GP smooths (same transform as in 08_run_gams.R GP) -----
xy_all <- dat_gam %>%
	transmute(x_raw = lat * cos(long * (pi/180)),
			  y_raw = lat * sin(long * (pi/180)))
rmax <- max(sqrt(xy_all$x_raw^2 + xy_all$y_raw^2), na.rm = TRUE)

chan_grid <- latlong_map %>%
	mutate(
		x = (lat * cos(long * (pi/180))) / rmax,
		y = (lat * sin(long * (pi/180))) / rmax
	)

# ----- Levels from the fitted model (to guarantee matching) -----
lev_group <- levels(gam_full$model$group)
lev_band  <- levels(gam_full$model$band)
lev_epoch <- levels(gam_full$model$epoch)
lev_rep   <- levels(gam_full$model$rep)
lev_part  <- levels(gam_full$model$participant)

# Detect cell structure (with or without rep)
if (!"cell" %in% names(gam_full$model))
	stop("Model frame has no 'cell'. This script expects GP + acc_bin factor-smooths with a 'cell' factor.")
lev_cell <- levels(gam_full$model$cell)
cell_has_rep <- any(grepl("random|repeated", lev_cell))

# Helper to build a valid 'cell' label (tries with rep first, then without)
make_cell_level <- function(g, b, e, accb, r = NULL) {
	# with rep
	if (!is.null(r)) {
		c1 <- as.character(interaction(
			factor(g, levels = lev_group),
			factor(b, levels = lev_band),
			factor(e, levels = lev_epoch),
			factor(accb, levels = c("Low","High")),
			factor(r, levels = lev_rep),
			drop = TRUE
		))
		if (c1 %in% lev_cell) return(c1)
	}
	# without rep
	c2 <- as.character(interaction(
		factor(g, levels = lev_group),
		factor(b, levels = lev_band),
		factor(e, levels = lev_epoch),
		factor(accb, levels = c("Low","High")),
		drop = TRUE
	))
	if (c2 %in% lev_cell) return(c2)
	stop("Could not match a 'cell' level for: ",
		 paste(g, b, e, accb, r, collapse = " / "))
}

# Build newdata rows at one channel
mk_row <- function(base_xy, g, b, e, r, accb, participant = lev_part[1]) {
	cell_char <- make_cell_level(g, b, e, accb, if (cell_has_rep) r else NULL)
	base_xy %>%
		mutate(
			group       = factor(g, levels = lev_group),
			band        = factor(b, levels = lev_band),
			epoch       = factor(e, levels = lev_epoch),
			rep         = factor(r, levels = lev_rep),
			acc_bin     = factor(accb, levels = c("Low","High")),
			cell        = factor(cell_char, levels = lev_cell),
			participant = factor(participant, levels = lev_part)
		)
}

# Design row helper (exclude subject RE)
X_of <- function(nd) {
	predict.gam(gam_full, newdata = nd, type = "lpmatrix",
				exclude = c("s(participant)"), discrete = FALSE)
}

# ----- Contrasts for the 6 maps -----
contrast_LminusH_group <- function(chan_row, g, b, e) {
	base <- chan_row %>% transmute(x, y)
	rows <- list(
		lo_rand = mk_row(base, g, b, e, "random",   "Low"),
		lo_rep  = mk_row(base, g, b, e, "repeated", "Low"),
		hi_rand = mk_row(base, g, b, e, "random",   "High"),
		hi_rep  = mk_row(base, g, b, e, "repeated", "High")
	)
	Xs <- lapply(rows, X_of)
	as.numeric( (Xs$lo_rand + Xs$lo_rep)/2 - (Xs$hi_rand + Xs$hi_rep)/2 )
}

contrast_PMBR_group <- function(chan_row, g) {
	base <- chan_row %>% transmute(x, y)
	b <- "beta"
	# during
	d_lo_rand <- mk_row(base, g, b, "during", "random",   "Low")
	d_lo_rep  <- mk_row(base, g, b, "during", "repeated", "Low")
	d_hi_rand <- mk_row(base, g, b, "during", "random",   "High")
	d_hi_rep  <- mk_row(base, g, b, "during", "repeated", "High")
	# after
	a_lo_rand <- mk_row(base, g, b, "after", "random",    "Low")
	a_lo_rep  <- mk_row(base, g, b, "after", "repeated",  "Low")
	a_hi_rand <- mk_row(base, g, b, "after", "random",    "High")
	a_hi_rep  <- mk_row(base, g, b, "after", "repeated",  "High")

	X <- lapply(list(d_lo_rand,d_lo_rep,d_hi_rand,d_hi_rep,
					 a_lo_rand,a_lo_rep,a_hi_rand,a_hi_rep), X_of)
	names(X) <- c("d_lo_r","d_lo_p","d_hi_r","d_hi_p","a_lo_r","a_lo_p","a_hi_r","a_hi_p")

	D <- 0.5*((X$d_lo_r + X$d_lo_p)/2 + (X$d_hi_r + X$d_hi_p)/2)
	A <- 0.5*((X$a_lo_r + X$a_lo_p)/2 + (X$a_hi_r + X$a_hi_p)/2)
	as.numeric(A - D)
}

# ---- Manifest for 6 maps ----
manifest <- tibble::tribble(
	~id, ~title,                                                   ~fn,                    ~args,
	"OM_theta_during_LH", "OM: (Low-High) @ Theta, During",        contrast_LminusH_group, list(g="physical", b="theta", e="during"),
	"MI_theta_during_LH", "MI: (Low-High) @ Theta, During",        contrast_LminusH_group, list(g="imagery",  b="theta", e="during"),
	"OM_alpha_after_LH",  "OM: (Low-High) @ Alpha, After",         contrast_LminusH_group, list(g="physical", b="alpha", e="after"),
	"MI_alpha_after_LH",  "MI: (Low-High) @ Alpha, After",         contrast_LminusH_group, list(g="imagery",  b="alpha", e="after"),
	"OM_beta_PMBR",       "OM: (After-During) @ Beta (mid acc)",   contrast_PMBR_group,    list(g="physical"),
	"MI_beta_PMBR",       "MI: (After-During) @ Beta (mid acc)",   contrast_PMBR_group,    list(g="imagery")
)

# ---- Compute ĉ, SE, t per map ----
beta_hat <- coef(gam_full)
V        <- vcov(gam_full, unconditional = TRUE)

one_map <- function(man_row) {
	# Extract the function and argument list cleanly from a 1-row tibble:
	fn   <- man_row$fn[[1]]
	args <- man_row$args[[1]]

	L_rows <- do.call(rbind, lapply(split(chan_grid, chan_grid$chan), function(one){
		do.call(fn, c(list(one), args))
	}))
	rownames(L_rows) <- chan_grid$chan

	c_hat <- as.numeric(L_rows %*% beta_hat)
	Sigma <- L_rows %*% V %*% t(L_rows)
	se    <- sqrt(pmax(diag(Sigma), 0))
	t_obs <- ifelse(se > 0, c_hat / se, 0)

	tibble(chan = chan_grid$chan,
		   c_hat = c_hat,
		   se    = se,
		   t_obs = t_obs)
}

# Build all tables
maps <- manifest
maps$tab <- purrr::map(seq_len(nrow(manifest)), ~ one_map(manifest[.x, , drop = FALSE]))

# Optionally compute a common color scale across all six
if (common_scale) {
	lim_all <- max(abs(unlist(lapply(maps$tab, function(tt) tt$c_hat))), na.rm = TRUE)
} else {
	lim_all <- NA_real_
}

# ---- Plot helpers ----
scale_to_0range <- function(x, range = 8) {
	x <- x - min(x); x <- x / max(x); (x - .5) * range
}

plot_grid <- chan_grid %>%
	mutate(
		x_plot = latlong_map$lat[match(chan, latlong_map$chan)] * cos(latlong_map$long[match(chan, latlong_map$chan)] * (pi/180)),
		y_plot = latlong_map$lat[match(chan, latlong_map$chan)] * sin(latlong_map$long[match(chan, latlong_map$chan)] * (pi/180))
	) %>%
	mutate(
		x_scaled = scale_to_0range(x_plot, 8),
		y_scaled = scale_to_0range(y_plot, 8.5)
	) %>% select(chan, x_scaled, y_scaled)

make_plot <- function(id, title, tab, lim_all) {
	dfp <- tab %>% left_join(plot_grid, by = "chan") %>%
		mutate(sig = abs(t_obs) >= t_cutoff)

	lim <- if (is.finite(lim_all)) lim_all else max(abs(dfp$c_hat), na.rm = TRUE)
	if (!is.finite(lim) || lim <= 0) lim <- 1

	p <- ggplot(dfp, aes(x = x_scaled, y = y_scaled)) +
		geom_point(aes(fill = c_hat), shape = 21, size = 6, stroke = 0.5, color = "grey30") +
		geom_point(data = dfp %>% filter(sig),
				   aes(x = x_scaled, y = y_scaled),
				   shape = 21, size = 8, stroke = 1.1, color = "black", fill = NA) +
		{ if (add_chan_labels)
			geom_text(data = dfp, aes(label = chan), nudge_y = 0.35, size = 2.8, color = "grey20")
			else NULL } +
		scale_fill_gradient2(name = "c_hat", limits = c(-lim, lim), midpoint = 0) +
		coord_equal() +
		theme_minimal(base_size = 12) +
		theme(axis.title = element_blank(), axis.text = element_blank(),
			  axis.ticks = element_blank(), panel.grid = element_blank(),
			  legend.position = "right",
			  plot.title = element_text(face = "bold")) +
		labs(title = title,
			 subtitle = sprintf("Descriptive rings: |t| >= %.2f (no multiplicity adj.)", t_cutoff))

	out_csv  <- file.path(out_tables, paste0("posthoc_", id, ".csv"))
	out_plot <- file.path(out_plots,  paste0("posthoc_", id, ".pdf"))
	write.csv(dfp %>% select(chan, c_hat, se, t_obs, sig), out_csv, row.names = FALSE)
	ggsave(out_plot, p, width = 7, height = 6)

	message("Saved: ", out_csv, " ; ", out_plot)
	invisible(dfp)
}

# ---- Render all six maps ----
for (i in seq_len(nrow(manifest))) {
	make_plot(manifest$id[i], manifest$title[i], maps$tab[[i]], lim_all)
}

message("\nDone: 6 descriptive per-group maps (tables + PDFs).")
