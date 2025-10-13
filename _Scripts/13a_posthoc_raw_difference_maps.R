#########################################
### 13a_posthoc_raw_difference_maps.R ###
#########################################
suppressPackageStartupMessages({
	library(tidyverse)
	library(mgcv)
	library(readr)
	library(grid)      # unit()
	library(patchwork) # combining panels
})

message("\n=== 13a_posthoc_raw_maps.R ===")

# -------------------- Paths --------------------
path_rds   <- "_Scripts/_rds/"
path_misc  <- "_Scripts/_misc/"
out_plots  <- "_Plots/"
out_tables <- "_Tables/"
dir.create(out_plots,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_tables, recursive = TRUE, showWarnings = FALSE)

# -------------------- Controls --------------------
t_thresh <- 2.0        # uncorrected “rings” on sensors (descriptive only)

# Smoother for continuous field built from SUBJECT-MEAN sensor values
smooth_bs <- "tp"       # try "gp" if you prefer; "tp" is robust
k_xy      <- 60         # upper limit; *adaptive per facet* below

# Plot knobs
sensor_size     <- 6.0
outline_size    <- 1.0
add_chan_labels <- FALSE
label_size      <- 2.5
label_nudge_y   <- 0.35

x_span_units <- 8.0
y_span_units <- 8.5

panel_spacing_x_lines <- 1.3
panel_spacing_y_lines <- 0.8

share_color_scale <- FALSE   # FALSE = per-facet auto-scale; TRUE = shared
clip_to_pct       <- NA      # e.g., 95 to clip; NA = no clipping

# Output
save_png_also <- TRUE
png_dpi       <- 300

# -------------------- Helpers --------------------
scale_to_0range <- function(x, span = 8) {
	xr <- range(x, na.rm = TRUE)
	if (!is.finite(xr[1]) || !is.finite(xr[2]) || xr[1] == xr[2]) return(rep(0, length(x)))
	x0 <- (x - xr[1]) / (xr[2] - xr[1])
	(x0 - 0.5) * span
}

# Use the same Low/High cutpoints as in 08_run_gams
cutfile <- file.path(path_rds, "accuracy_bin_cutpoints.rds")
if (!file.exists(cutfile)) stop("Missing accuracy cutpoints: ", cutfile)
acc_info <- readRDS(cutfile)
q <- as.numeric(acc_info$q)  # 20th & 80th percentiles

mk_acc_bin <- function(acc) case_when(
	acc <= q[1] ~ "Low",
	acc >= q[2] ~ "High",
	TRUE ~ NA_character_
)

# -------------------- Load data --------------------
dat <- readRDS(file.path(path_rds, "dat_gam.rds")) %>% as_tibble() %>%
	mutate(
		group = factor(group, levels = c("physical","imagery")),
		band  = factor(band,  levels = c("theta","alpha","beta")),
		epoch = factor(epoch, levels = c("during","after")),
		rep   = factor(rep,   levels = c("random","repeated")),
		chan  = factor(chan),
		participant = factor(participant),
		acc_bin = mk_acc_bin(accuracy)
	) %>% droplevels()

# Channel coords
map_file <- file.path(path_misc, "sensor_latlong_chan_map.csv")
if (!file.exists(map_file)) stop("Missing sensor_latlong_chan_map.csv at: ", map_file)
latlong_map <- read_csv(map_file, show_col_types = FALSE) %>% distinct(chan, lat, long)

# Plotting coordinates
coords <- latlong_map %>%
	mutate(
		x = lat * cos(long * (pi/180)),
		y = lat * sin(long * (pi/180)),
		x_scaled = scale_to_0range(x, x_span_units),
		y_scaled = scale_to_0range(y, y_span_units)
	)

# -------------------- Subject-level effects (per group) --------------------
# H1: Low-High @ theta,during
summ_h1 <- dat %>%
	filter(band=="theta", epoch=="during", !is.na(acc_bin)) %>%
	group_by(participant, group, chan, acc_bin) %>%
	summarise(m = mean(powerdb, na.rm = TRUE), .groups = "drop_last") %>%
	summarise(effect = m[acc_bin=="Low"] - m[acc_bin=="High"], .groups = "drop") %>%
	group_by(group, chan) %>%
	summarise(n = n(),
			  mean_eff = mean(effect, na.rm = TRUE),
			  sd_eff   = sd(effect,   na.rm = TRUE),
			  se_eff   = sd_eff / sqrt(n),
			  t_eff    = ifelse(se_eff>0, mean_eff/se_eff, 0),
			  .groups = "drop") %>%
	mutate(h = "H1: (Low-High) @ Theta, During")

# H2: Low-High @ alpha,after
summ_h2 <- dat %>%
	filter(band=="alpha", epoch=="after", !is.na(acc_bin)) %>%
	group_by(participant, group, chan, acc_bin) %>%
	summarise(m = mean(powerdb, na.rm = TRUE), .groups = "drop_last") %>%
	summarise(effect = m[acc_bin=="Low"] - m[acc_bin=="High"], .groups = "drop") %>%
	group_by(group, chan) %>%
	summarise(n = n(),
			  mean_eff = mean(effect, na.rm = TRUE),
			  sd_eff   = sd(effect,   na.rm = TRUE),
			  se_eff   = sd_eff / sqrt(n),
			  t_eff    = ifelse(se_eff>0, mean_eff/se_eff, 0),
			  .groups = "drop") %>%
	mutate(h = "H2: (Low-High) @ Alpha, After")

# H3: PMBR ≈ (After - During) at mid ≈ 0.5·Low + 0.5·High
summ_h3 <- dat %>%
	filter(band=="beta", !is.na(acc_bin)) %>%
	group_by(participant, group, chan, epoch, acc_bin) %>%
	summarise(m = mean(powerdb, na.rm = TRUE), .groups = "drop_last") %>%
	summarise(mid = 0.5*m[acc_bin=="Low"] + 0.5*m[acc_bin=="High"], .groups = "drop_last") %>%
	summarise(effect = mid[epoch=="after"] - mid[epoch=="during"], .groups = "drop") %>%
	group_by(group, chan) %>%
	summarise(n = n(),
			  mean_eff = mean(effect, na.rm = TRUE),
			  sd_eff   = sd(effect,   na.rm = TRUE),
			  se_eff   = sd_eff / sqrt(n),
			  t_eff    = ifelse(se_eff>0, mean_eff/se_eff, 0),
			  .groups = "drop") %>%
	mutate(h = "H3: PMBR (After - During) @ Beta (mid ≈ 0.5·Low + 0.5·High)")

summ_all <- bind_rows(summ_h1, summ_h2, summ_h3) %>%
	left_join(coords, by = "chan") %>%
	mutate(
		facet = paste0(h, " — ", ifelse(group=="physical", "OM", "MI")),
		facet = factor(facet, levels = c(
			"H1: (Low-High) @ Theta, During — OM",
			"H1: (Low-High) @ Theta, During — MI",
			"H2: (Low-High) @ Alpha, After — OM",
			"H2: (Low-High) @ Alpha, After — MI",
			"H3: PMBR (After - During) @ Beta (mid ≈ 0.5·Low + 0.5·High) — OM",
			"H3: PMBR (After - During) @ Beta (mid ≈ 0.5·Low + 0.5·High) — MI"
		))
	)

# Save per-sensor summaries for record
write_csv(summ_all %>% select(group, h, chan, n, mean_eff, se_eff, t_eff),
		  file.path(out_tables, "posthoc_raw_sensor_summaries.csv"))
message("Saved: ", file.path(out_tables, "posthoc_raw_sensor_summaries.csv"))

# -------------------- Build prediction grid --------------------
xg <- seq(min(coords$x_scaled), max(coords$x_scaled), length.out = 121)
yg <- seq(min(coords$y_scaled), max(coords$y_scaled), length.out = 121)
grid <- expand.grid(x_scaled = xg, y_scaled = yg)
r <- sqrt((grid$x_scaled / (x_span_units/2))^2 + (grid$y_scaled / (y_span_units/2))^2)
grid <- grid %>% mutate(in_head = r <= 1.02)

# Adaptive k per facet + prediction
fit_predict_field <- function(df, facet_name) {
	df_fit <- df %>% select(x_scaled, y_scaled, mean_eff) %>% distinct()
	n_unique <- nrow(df_fit)
	# choose k adaptively: <= n_unique - 3, with a floor of 10
	k_use <- max(10L, min(k_xy, n_unique - 3L))
	message(sprintf("Facet '%s': n_unique=%d, using k=%d (bs=%s)",
					facet_name, n_unique, k_use, smooth_bs))
	# fit smoother with adaptive k; fallback with smaller k if needed
	m <- try(
		mgcv::gam(mean_eff ~ s(x_scaled, y_scaled, bs = smooth_bs, k = k_use), data = df_fit),
		silent = TRUE
	)
	if (inherits(m, "try-error")) {
		k_use2 <- max(10L, floor(n_unique/2))
		message(sprintf("  retry with k=%d", k_use2))
		m <- mgcv::gam(mean_eff ~ s(x_scaled, y_scaled, bs = smooth_bs, k = k_use2), data = df_fit)
	}
	grid$z <- predict(m, newdata = grid, type = "response")
	grid
}

# Shared color limits (optional)
if (share_color_scale) {
	vals <- summ_all$mean_eff
	if (is.na(clip_to_pct)) {
		lim <- max(abs(vals), na.rm = TRUE)
	} else {
		lo  <- quantile(vals, (100 - clip_to_pct)/200, na.rm = TRUE)
		hi  <- quantile(vals, 1 - (100 - clip_to_pct)/200, na.rm = TRUE)
		lim <- max(abs(c(lo, hi)))
	}
	if (!is.finite(lim) || lim <= 0) lim <- 1
	limits_shared <- c(-lim, lim)
} else {
	limits_shared <- NULL
}

make_plot <- function(df_one, title_tag) {
	GG <- fit_predict_field(df_one, title_tag)
	# per-facet limits or shared
	if (is.null(limits_shared)) {
		lim <- max(abs(df_one$mean_eff), na.rm = TRUE)
		if (!is.finite(lim) || lim <= 0) lim <- 1
		limits <- c(-lim, lim)
	} else {
		limits <- limits_shared
	}
	ggplot() +
		geom_raster(
			data = GG %>% filter(in_head),
			aes(x = x_scaled, y = y_scaled, fill = z),
			interpolate = TRUE
		) +
		scale_fill_gradient2(name = "Raw mean effect", midpoint = 0, limits = limits) +
		geom_point(
			data = df_one,
			aes(x = x_scaled, y = y_scaled),
			shape = 21, size = sensor_size, stroke = 0.5, color = "grey30", fill = "white"
		) +
		geom_point(
			data = df_one %>% filter(abs(t_eff) >= t_thresh),
			aes(x = x_scaled, y = y_scaled),
			shape = 21, size = sensor_size + 1.2, stroke = outline_size,
			color = "black", fill = NA
		) +
		{ if (add_chan_labels)
			geom_text(data = df_one, aes(x = x_scaled, y = y_scaled, label = chan),
					  nudge_y = label_nudge_y, size = label_size, color = "grey20")
			else NULL } +
		coord_equal() +
		theme_minimal(base_size = 12) +
		theme(
			axis.title = element_blank(),
			axis.text  = element_blank(),
			axis.ticks = element_blank(),
			panel.grid = element_blank(),
			panel.background = element_rect(fill = "white", colour = NA),
			plot.background  = element_rect(fill = "white", colour = NA),
			legend.background= element_rect(fill = "white", colour = NA),
			plot.title = element_text(face = "bold")
		) +
		labs(title = title_tag)
}

# Build six panels
df_plot <- summ_all %>% arrange(facet)
panels <- df_plot %>%
	group_split(facet, .keep = TRUE) %>%
	purrr::map(~ make_plot(.x, unique(.x$facet)))

# 2x3 layout
figure <- (panels[[1]] | panels[[2]] | panels[[3]]) /
	(panels[[4]] | panels[[5]] | panels[[6]])

figure <- figure + plot_layout(guides = "collect") &
	theme(legend.position = "right",
		  plot.margin = margin(8, 12, 8, 12),
		  panel.spacing.x = unit(panel_spacing_x_lines, "lines"),
		  panel.spacing.y = unit(panel_spacing_y_lines, "lines"))

outfile_pdf <- file.path(out_plots, "figure_posthoc_raw_maps.pdf")
try(suppressWarnings(file.remove(outfile_pdf)), silent = TRUE)
ggsave(outfile_pdf, figure, width = 12, height = 7.8, device = grDevices::pdf, bg = "white")
message("Saved: ", outfile_pdf)

if (save_png_also) {
	outfile_png <- file.path(out_plots, "figure_posthoc_raw_maps.png")
	ggsave(outfile_png, figure, width = 12, height = 7.8, dpi = png_dpi, bg = "white")
	message("Saved: ", outfile_png)
}
