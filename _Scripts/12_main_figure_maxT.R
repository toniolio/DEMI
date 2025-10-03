#############################
### 12_main_figure_maxT.R ###
#############################
# plots figure for the main result

suppressPackageStartupMessages({
	library(tidyverse)
	library(readr)
	library(ggplot2)
	library(grid)   # for unit()
})

message("\n=== 12_main_figure_maxT.R ===")

# -------------------- Paths --------------------
path_tables <- "_Tables/"
path_misc   <- "_Scripts/_misc/"
out_plots   <- "_Plots/"
dir.create(out_plots, recursive = TRUE, showWarnings = FALSE)

# -------------------- Plot controls --------------------
sensor_size      <- 10     # point size for sensors
outline_size     <- 1.1     # stroke for the black ring (significant sensors)
add_chan_labels  <- FALSE    # TRUE -> print channel codes
label_size       <- 2.7
label_nudge_y    <- 0.5

# Layout spacing (positive): spreads points within each scalp
x_span_units     <- 6.0
y_span_units     <- 6.8

# Facet spacing (distance between the three panels)
panel_spacing_x_lines <- 2   # horizontal gap between panels (in "lines")
panel_spacing_y_lines <- 0.6   # vertical gap (not important in 1 row, but set anyway)

share_color_scale <- TRUE   # one symmetric color scale across H1–H3
clip_to_pct       <- NA     # e.g., 95 to clip to central 95%; NA = no clipping

# Output
save_png_also     <- TRUE
png_width         <- 12.0
png_height        <- 5
dpi_png           <- 300

# Always use the base PDF device (no Cairo)
pdf_width         <- 12.0
pdf_height        <- 5

# -------------------- Helpers --------------------
scale_to_0range <- function(x, span = 8) {
	xr <- range(x, na.rm = TRUE)
	if (!is.finite(xr[1]) || !is.finite(xr[2]) || xr[1] == xr[2]) return(rep(0, length(x)))
	x0 <- (x - xr[1]) / (xr[2] - xr[1])
	(x0 - 0.5) * span
}

read_one <- function(csv_path, facet_title) {
	if (!file.exists(csv_path)) stop("Missing: ", csv_path)
	read_csv(csv_path, show_col_types = FALSE) %>%
		transmute(
			chan = chan,
			c_hat = effect_c_hat,
			t_obs = t_obs,
			sig   = as.logical(sig_FWER05_global),
			facet = facet_title
		)
}

safe_save_pdf_base <- function(filename, plot, width, height) {
	# If the PDF is open in a viewer, overwrite may fail silently. Try to remove first.
	try(suppressWarnings(file.remove(filename)), silent = TRUE)
	ggplot2::ggsave(
		filename, plot,
		width = width, height = height,
		device = grDevices::pdf, # base device, no Cairo
		bg = "white"             # ensure white background
	)
	if (file.exists(filename) && file.info(filename)$size > 0) {
		message("Saved: ", filename)
	} else {
		warning("Failed to save PDF at: ", filename,
				" (is it open in a viewer or write-protected?)")
	}
}

# -------------------- Load result tables --------------------
paths <- c(
	h1 = file.path(path_tables, "maxT_results_h1.csv"),
	h2 = file.path(path_tables, "maxT_results_h2.csv"),
	h3 = file.path(path_tables, "maxT_results_h3.csv")
)

df_h <- bind_rows(
	read_one(paths["h1"], "H1: (Low-High)_OM - (Low-High)_MI @ Theta, During"),
	read_one(paths["h2"], "H2: (Low-High)_OM - (Low-High)_MI @ Alpha, After"),
	read_one(paths["h3"], "H3: (After-During)_OM - (After-During)_MI @ Beta (mid ~= 0.5*Low + 0.5*High)")
)

# -------------------- Channel coordinates --------------------
map_file <- file.path(path_misc, "sensor_latlong_chan_map.csv")
if (!file.exists(map_file)) stop("Missing sensor_latlong_chan_map.csv at: ", map_file)

latlong_map <- read_csv(map_file, show_col_types = FALSE) %>%
	distinct(chan, lat, long)

missing_ch <- df_h$chan[!(df_h$chan %in% latlong_map$chan)]
if (length(missing_ch)) {
	warning("Channels in results but missing from coord map: ",
			paste(sort(unique(missing_ch)), collapse = ", "))
}

# Build plotting coordinates
dfp <- df_h %>%
	left_join(latlong_map, by = "chan") %>%
	mutate(
		x = lat * cos(long * (pi/180)),
		y = lat * sin(long * (pi/180)),
		x_scaled = scale_to_0range(x, x_span_units),
		y_scaled = scale_to_0range(y, y_span_units),
		facet = factor(
			facet,
			levels = c(
				"H1: (Low-High)_OM - (Low-High)_MI @ Theta, During",
				"H2: (Low-High)_OM - (Low-High)_MI @ Alpha, After",
				"H3: (After-During)_OM - (After-During)_MI @ Beta (mid ~= 0.5*Low + 0.5*High)"
			)
		)
	)

# -------------------- Color limits --------------------
if (share_color_scale) {
	if (is.na(clip_to_pct)) {
		lim <- max(abs(dfp$c_hat), na.rm = TRUE)
	} else {
		lo  <- quantile(dfp$c_hat, (100 - clip_to_pct) / 200, na.rm = TRUE)
		hi  <- quantile(dfp$c_hat, 1 - (100 - clip_to_pct) / 200, na.rm = TRUE)
		lim <- max(abs(c(lo, hi)))
	}
	if (!is.finite(lim) || lim <= 0) lim <- 1
	fill_limits <- c(-lim, lim)
} else {
	fill_limits <- NULL
}

# -------------------- Plot (HORIZONTAL: 1 row, 3 columns) --------------------
p <- ggplot(dfp, aes(x = x_scaled, y = y_scaled)) +
	geom_point(aes(fill = c_hat),
			   shape = 21, size = sensor_size, stroke = 0.5, color = "grey30") +
	geom_point(data = dfp %>% filter(sig),
			   shape = 21, size = sensor_size + 1.5, stroke = outline_size,
			   color = "black", fill = NA) +
	{ if (add_chan_labels)
		geom_text(data = dfp, aes(label = chan),
				  nudge_y = label_nudge_y, size = label_size, color = "grey20")
		else NULL } +
	scale_fill_gradient2(
		name     = "Effect (c_hat)",
		midpoint = 0,
		limits   = fill_limits
	) +
	coord_equal() +
	facet_wrap(~ facet, ncol = 3, scales = "fixed") +
	theme_minimal(base_size = 12) +
	theme(
		axis.title = element_blank(),
		axis.text  = element_blank(),
		axis.ticks = element_blank(),
		panel.grid = element_blank(),
		legend.position = "right",
		strip.text = element_text(face = "bold"),
		# White backgrounds for dark-mode-friendly viewing
		panel.background = element_rect(fill = "white", colour = NA),
		plot.background  = element_rect(fill = "white", colour = NA),
		legend.background = element_rect(fill = "white", colour = NA),
		# Increase spacing between the 3 panels
		panel.spacing.x = unit(panel_spacing_x_lines, "lines"),
		panel.spacing.y = unit(panel_spacing_y_lines, "lines"),
		# Extra outer margins if you need a bit more breathing room
		plot.margin = margin(10, 14, 10, 14)
	)

# -------------------- Save --------------------
outfile_pdf <- file.path(out_plots, "figure_main_maxT.pdf")
safe_save_pdf_base(outfile_pdf, p, width = pdf_width, height = pdf_height)

if (save_png_also) {
	outfile_png <- file.path(out_plots, "figure_main_maxT.png")
	ggplot2::ggsave(outfile_png, p, width = png_width, height = png_height, dpi = dpi_png, bg = "white")
	message("Saved: ", outfile_png)
}
