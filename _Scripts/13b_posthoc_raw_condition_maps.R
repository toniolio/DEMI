########################################
### 13b_posthoc_raw_condition_maps.R ###
########################################
suppressPackageStartupMessages({
	library(tidyverse)
	library(mgcv)
	library(readr)
	library(grid)      # unit()
	library(patchwork) # combine panels + collect legend
	# Optional: uncomment if you want non-overlapping labels
	# library(ggrepel)
})

message("\n=== 13b_posthoc_raw_condition_maps.R ===")

# -------------------- Paths --------------------
path_rds   <- "_Scripts/_rds/"
path_misc  <- "_Scripts/_misc/"
out_plots  <- "_Plots/"
out_tables <- "_Tables/"
dir.create(out_plots,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_tables, recursive = TRUE, showWarnings = FALSE)

# -------------------- Controls --------------------
# What to show for each hypothesis QC:
acc_show_H1 <- "Low"   # theta,during: show Low-accuracy raw map (OM & MI)
acc_show_H2 <- "Low"   # alpha,after : show Low-accuracy raw map (OM & MI)

# PMBR "after" accuracy strategy:
#   "mid" = 0.5*Low + 0.5*High,  "low" = Low only,  "high" = High only,  "all" = all accuracies pooled
acc_strategy_H3 <- "mid"

# Uncorrected t-ring threshold (purely descriptive)
t_thresh <- 2.0

# Smoother for continuous field (built from subject-mean sensor values)
smooth_bs <- "tp"     # "tp" is robust; you can try "gp"
k_xy      <- 60       # upper bound; we adapt per panel: k_use = min(k_xy, n_unique - 3), floor 10

# Plot knobs
sensor_size      <- 6.0
outline_size     <- 1.0
add_chan_labels  <- FALSE   # TRUE to draw channel codes
label_size       <- 2.5
label_nudge_y    <- 0.35

# Layout spacing
x_span_units     <- 8.0
y_span_units     <- 8.5
panel_spacing_x_lines <- 1.5
panel_spacing_y_lines <- 1.0

# Color scaling
share_color_scale <- FALSE   # one symmetric color scale across all six panels
clip_to_pct       <- NA     # e.g., 95 to clip to central 95%; NA = no clipping

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

# Channel coords (for plotting & smoother input)
map_file <- file.path(path_misc, "sensor_latlong_chan_map.csv")
if (!file.exists(map_file)) stop("Missing sensor_latlong_chan_map.csv at: ", map_file)
coords <- read_csv(map_file, show_col_types = FALSE) %>%
	distinct(chan, lat, long) %>%
	mutate(
		x = lat * cos(long * (pi/180)),
		y = lat * sin(long * (pi/180)),
		x_scaled = scale_to_0range(x, x_span_units),
		y_scaled = scale_to_0range(y, y_span_units)
	)

# -------------------- Build subject-level RAW means for selected conditions ----
# H1: theta, during, acc_show_H1
h1_raw <- dat %>%
	filter(band=="theta", epoch=="during", acc_bin %in% acc_show_H1) %>%
	group_by(participant, group, chan) %>%
	summarise(raw = mean(powerdb, na.rm = TRUE), .groups = "drop") %>%
	group_by(group, chan) %>%
	summarise(n = n(),
			  mean_raw = mean(raw, na.rm = TRUE),
			  sd_raw   = sd(raw,   na.rm = TRUE),
			  se_raw   = sd_raw / sqrt(n),
			  t_raw    = ifelse(se_raw>0, mean_raw/se_raw, 0),
			  .groups = "drop") %>%
	mutate(h = sprintf("H1 QC: Theta, During, %s", acc_show_H1))

# H2: alpha, after, acc_show_H2
h2_raw <- dat %>%
	filter(band=="alpha", epoch=="after", acc_bin %in% acc_show_H2) %>%
	group_by(participant, group, chan) %>%
	summarise(raw = mean(powerdb, na.rm = TRUE), .groups = "drop") %>%
	group_by(group, chan) %>%
	summarise(n = n(),
			  mean_raw = mean(raw, na.rm = TRUE),
			  sd_raw   = sd(raw,   na.rm = TRUE),
			  se_raw   = sd_raw / sqrt(n),
			  t_raw    = ifelse(se_raw>0, mean_raw/se_raw, 0),
			  .groups = "drop") %>%
	mutate(h = sprintf("H2 QC: Alpha, After, %s", acc_show_H2))

# H3: beta, AFTER only, accuracy strategy
h3_raw <- dat %>%
	filter(band=="beta", epoch=="after", !is.na(acc_bin)) %>%
	group_by(participant, group, chan, acc_bin) %>%
	summarise(m = mean(powerdb, na.rm = TRUE), .groups = "drop_last") %>%
	group_by(participant, group, chan) %>%
	summarise(
		raw = case_when(
			acc_strategy_H3 == "low" ~ m[acc_bin=="Low"],
			acc_strategy_H3 == "high" ~ m[acc_bin=="High"],
			acc_strategy_H3 == "mid" ~ 0.5*m[acc_bin=="Low"] + 0.5*m[acc_bin=="High"],
			acc_strategy_H3 == "all" ~ mean(m, na.rm = TRUE),
			TRUE ~ NA_real_
		),
		.groups = "drop"
	) %>%
	group_by(group, chan) %>%
	summarise(n = n(),
			  mean_raw = mean(raw, na.rm = TRUE),
			  sd_raw   = sd(raw,   na.rm = TRUE),
			  se_raw   = sd_raw / sqrt(n),
			  t_raw    = ifelse(se_raw>0, mean_raw/se_raw, 0),
			  .groups = "drop") %>%
	mutate(h = sprintf("H3 QC: Beta, After (%s accuracy)", acc_strategy_H3))

summ_all <- bind_rows(h1_raw, h2_raw, h3_raw) %>%
	left_join(coords, by = "chan") %>%
	mutate(
		facet = paste0(h, " — ", ifelse(group=="physical", "OM", "MI")),
		facet = factor(facet, levels = c(
			sprintf("H1 QC: Theta, During, %s — OM", acc_show_H1),
			sprintf("H1 QC: Theta, During, %s — MI", acc_show_H1),
			sprintf("H2 QC: Alpha, After, %s — OM", acc_show_H2),
			sprintf("H2 QC: Alpha, After, %s — MI", acc_show_H2),
			sprintf("H3 QC: Beta, After (%s accuracy) — OM", acc_strategy_H3),
			sprintf("H3 QC: Beta, After (%s accuracy) — MI", acc_strategy_H3)
		))
	)

# Persist per-sensor summaries (for manuscript supplement if needed)
write_csv(summ_all %>% select(group, h, chan, n, mean_raw, se_raw, t_raw),
		  file.path(out_tables, "posthoc_raw_condition_sensor_summaries.csv"))
message("Saved: ", file.path(out_tables, "posthoc_raw_condition_sensor_summaries.csv"))

# -------------------- Grid & adaptive smoother per panel --------------------
xg <- seq(min(coords$x_scaled), max(coords$x_scaled), length.out = 121)
yg <- seq(min(coords$y_scaled), max(coords$y_scaled), length.out = 121)
grid <- expand.grid(x_scaled = xg, y_scaled = yg)
r <- sqrt((grid$x_scaled / (x_span_units/2))^2 + (grid$y_scaled / (y_span_units/2))^2)
grid <- grid %>% mutate(in_head = r <= 1.02)

# Optional shared color scale across all panels
if (share_color_scale) {
	vals <- summ_all$mean_raw
	if (is.na(clip_to_pct)) {
		lim <- max(abs(vals), na.rm = TRUE)
	} else {
		lo <- quantile(vals, (100 - clip_to_pct)/200, na.rm = TRUE)
		hi <- quantile(vals, 1 - (100 - clip_to_pct)/200, na.rm = TRUE)
		lim <- max(abs(c(lo, hi)))
	}
	if (!is.finite(lim) || lim <= 0) lim <- 1
	limits_shared <- c(-lim, lim)
} else {
	limits_shared <- NULL
}

fit_predict_field <- function(df, facet_name) {
	df_fit <- df %>% select(x_scaled, y_scaled, mean_raw) %>% distinct()
	n_unique <- nrow(df_fit)
	k_use <- max(10L, min(k_xy, n_unique - 3L))
	message(sprintf("Facet '%s': n_unique=%d, using k=%d (bs=%s)",
					facet_name, n_unique, k_use, smooth_bs))
	m <- try(
		mgcv::gam(mean_raw ~ s(x_scaled, y_scaled, bs = smooth_bs, k = k_use), data = df_fit),
		silent = TRUE
	)
	if (inherits(m, "try-error")) {
		k_use2 <- max(10L, floor(n_unique/2))
		message(sprintf("  retry with k=%d", k_use2))
		m <- mgcv::gam(mean_raw ~ s(x_scaled, y_scaled, bs = smooth_bs, k = k_use2), data = df_fit)
	}
	grid$z <- predict(m, newdata = grid, type = "response")
	grid
}

make_panel <- function(df_one, title_tag) {
	GG <- fit_predict_field(df_one, title_tag)
	limits <- if (is.null(limits_shared)) {
		lim <- max(abs(df_one$mean_raw), na.rm = TRUE); if (!is.finite(lim) || lim <= 0) lim <- 1
		c(-lim, lim)
	} else limits_shared

	base <- ggplot() +
		geom_raster(
			data = GG %>% filter(in_head),
			aes(x = x_scaled, y = y_scaled, fill = z),
			interpolate = TRUE
		) +
		scale_fill_gradient2(name = "Raw mean (dB)", midpoint = 0, limits = limits) +
		geom_point(
			data = df_one,
			aes(x = x_scaled, y = y_scaled),
			shape = 21, size = sensor_size, stroke = 0.5, color = "grey30", fill = "white"
		) +
		geom_point(
			data = df_one %>% filter(abs(t_raw) >= t_thresh),
			aes(x = x_scaled, y = y_scaled),
			shape = 21, size = sensor_size + 1.2, stroke = outline_size,
			color = "black", fill = NA
		) +
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

	if (add_chan_labels) {
		base <- base +
			# For non-overlapping labels, swap geom_text for geom_text_repel() from ggrepel
			geom_text(aes(x = x_scaled, y = y_scaled, label = chan),
					  data = df_one, nudge_y = label_nudge_y, size = label_size, color = "grey20")
		# + ggrepel::geom_text_repel(..., max.overlaps = Inf)  # optional
	}
	base
}

# Build six panels (OM, MI for each QC)
df_plot <- summ_all %>% arrange(facet)
panels <- df_plot %>% group_split(facet, .keep = TRUE) %>%
	purrr::map(~ make_panel(.x, unique(.x$facet)))

# Layout: 3 columns × 2 rows (horizontal groupings)
figure <- (panels[[1]] | panels[[2]] | panels[[3]]) /
	(panels[[4]] | panels[[5]] | panels[[6]])

figure <- figure + plot_layout(guides = "collect") &
	theme(legend.position = "right",
		  plot.margin = margin(8, 12, 8, 12),
		  panel.spacing.x = unit(panel_spacing_x_lines, "lines"),
		  panel.spacing.y = unit(panel_spacing_y_lines, "lines"))

outfile_pdf <- file.path(out_plots, "figure_posthoc_raw_conditions.pdf")
try(suppressWarnings(file.remove(outfile_pdf)), silent = TRUE)
ggsave(outfile_pdf, figure, width = 12, height = 7.8, device = grDevices::pdf, bg = "white")
message("Saved: ", outfile_pdf)

if (save_png_also) {
	outfile_png <- file.path(out_plots, "figure_posthoc_raw_conditions.png")
	ggsave(outfile_png, figure, width = 12, height = 7.8, dpi = png_dpi, bg = "white")
	message("Saved: ", outfile_png)
}
