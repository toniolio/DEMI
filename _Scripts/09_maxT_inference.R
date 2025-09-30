#############################################
### 09_maxT_inference_global.R (GP, robust)
#############################################
suppressPackageStartupMessages({
	library(tidyverse)
	library(mgcv)
	library(mvtnorm)
	library(readr)
})

message("\n=== 11_maxT_inference_global.R (GP, auto cell/dev) ===")

path_rds   <- "_Scripts/_rds/"
out_plots  <- "_Plots/"
out_tables <- "_Tables/"
dir.create(out_plots,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_tables, recursive = TRUE, showWarnings = FALSE)

alpha_fwer <- 0.05
B          <- 20000L
use_t_null <- TRUE
df_t       <- NULL

add_chan_labels <- FALSE
scale_by_ci95   <- FALSE

# ---- Load model & base data ----
gam_full <- readRDS(file.path(path_rds, "gam_full.rds"))
dat_gam  <- readRDS(file.path(path_rds, "dat_gam.rds")) %>% as_tibble()

dat_gam <- dat_gam %>%
	mutate(
		group = factor(group, levels = c("physical","imagery")),
		band  = factor(band,  levels = c("theta","alpha","beta")),
		epoch = factor(epoch, levels = c("during","after")),
		rep   = factor(rep,   levels = c("random","repeated")),
		chan  = factor(chan),
		participant = factor(participant)
	) %>% droplevels()

# channel grid (from dat_gam)
chan_grid <- dat_gam %>%
	group_by(chan) %>%
	summarise(lat = first(lat), long = first(long), .groups = "drop") %>%
	arrange(chan)

# ---- Match the GP x,y transform used at fit time ----
xy_all <- dat_gam %>%
	transmute(x_raw = lat * cos(long * (pi/180)),
			  y_raw = lat * sin(long * (pi/180)))
rmax <- max(sqrt(xy_all$x_raw^2 + xy_all$y_raw^2), na.rm = TRUE)

chan_grid <- chan_grid %>%
	mutate(
		x = (lat * cos(long * (pi/180))) / rmax,
		y = (lat * sin(long * (pi/180))) / rmax
	)

# ---- Levels from the fitted model frame ----
lev_group <- levels(gam_full$model$group)
lev_band  <- levels(gam_full$model$band)
lev_epoch <- levels(gam_full$model$epoch)
lev_rep   <- levels(gam_full$model$rep)
lev_part  <- levels(gam_full$model$participant)

has_cell <- "cell" %in% names(gam_full$model)
has_dev4 <- all(c("dev_theta_during_L","dev_theta_during_H",
				  "dev_alpha_after_L","dev_alpha_after_H") %in% names(gam_full$model))
has_dev8 <- all(c("dev_theta_during_L_phys","dev_theta_during_H_phys",
				  "dev_alpha_after_L_phys","dev_alpha_after_H_phys",
				  "dev_theta_during_L_imag","dev_theta_during_H_imag",
				  "dev_alpha_after_L_imag","dev_alpha_after_H_imag") %in% names(gam_full$model))

if (!has_cell && !has_dev4 && !has_dev8) {
	stop("Model has neither 'cell' nor expected dev_* indicators. Please refit with GP_v3 script.")
}
lev_cell <- if (has_cell) levels(gam_full$model$cell) else NULL

# ---- Helper: build a one-row newdata with the right indicators ----
mk_row <- function(base_xy, g, b, e, r, accbin, participant = lev_part[1]) {
	nd <- base_xy %>%
		mutate(
			group       = factor(g, levels = lev_group),
			band        = factor(b, levels = lev_band),
			epoch       = factor(e, levels = lev_epoch),
			rep         = factor(r, levels = lev_rep),
			participant = factor(participant, levels = lev_part)
		)

	if (has_cell) {
		cell_val <- interaction(nd$group, nd$band, nd$epoch,
								factor(accbin, levels = c("Low","High")),
								drop = TRUE)
		nd <- nd %>%
			mutate(
				acc_bin = factor(accbin, levels = c("Low","High")),
				cell    = factor(as.character(cell_val), levels = lev_cell)
			)
	} else if (has_dev4) {
		# zero all dev_*, then set the one matching (b,e,accbin) to 1
		nd$acc_bin <- factor(accbin, levels = c("Low","High"))
		nd$dev_theta_during_L <- as.numeric(b=="theta" & e=="during" & accbin=="Low")
		nd$dev_theta_during_H <- as.numeric(b=="theta" & e=="during" & accbin=="High")
		nd$dev_alpha_after_L  <- as.numeric(b=="alpha"  & e=="after"  & accbin=="Low")
		nd$dev_alpha_after_H  <- as.numeric(b=="alpha"  & e=="after"  & accbin=="High")
	} else if (has_dev8) {
		nd$acc_bin <- factor(accbin, levels = c("Low","High"))
		# zero all 8
		for (nm in c("dev_theta_during_L_phys","dev_theta_during_H_phys",
					 "dev_alpha_after_L_phys","dev_alpha_after_H_phys",
					 "dev_theta_during_L_imag","dev_theta_during_H_imag",
					 "dev_alpha_after_L_imag","dev_alpha_after_H_imag")) {
			nd[[nm]] <- 0
		}
		# turn on the matching one
		if (g == "physical") {
			nd$dev_theta_during_L_phys <- as.numeric(b=="theta" & e=="during" & accbin=="Low")
			nd$dev_theta_during_H_phys <- as.numeric(b=="theta" & e=="during" & accbin=="High")
			nd$dev_alpha_after_L_phys  <- as.numeric(b=="alpha"  & e=="after"  & accbin=="Low")
			nd$dev_alpha_after_H_phys  <- as.numeric(b=="alpha"  & e=="after"  & accbin=="High")
		} else {
			nd$dev_theta_during_L_imag <- as.numeric(b=="theta" & e=="during" & accbin=="Low")
			nd$dev_theta_during_H_imag <- as.numeric(b=="theta" & e=="during" & accbin=="High")
			nd$dev_alpha_after_L_imag  <- as.numeric(b=="alpha"  & e=="after"  & accbin=="Low")
			nd$dev_alpha_after_H_imag  <- as.numeric(b=="alpha"  & e=="after"  & accbin=="High")
		}
	}
	nd
}

X_of <- function(nd) {
	predict.gam(gam_full, newdata = nd, type = "lpmatrix",
				exclude = c("s(participant)"), discrete = FALSE)
}

# ---- Contrast per channel ----
contrast_row <- function(chan_row, spec) {
	base <- chan_row %>% transmute(x = x, y = y)

	if (spec$type == "acc_gxe") {
		b <- spec$band; e <- spec$epoch
		rows <- list(
			om_lo_rand = mk_row(base,"physical", b, e,"random","Low"),
			om_lo_rep  = mk_row(base,"physical", b, e,"repeated","Low"),
			om_hi_rand = mk_row(base,"physical", b, e,"random","High"),
			om_hi_rep  = mk_row(base,"physical", b, e,"repeated","High"),

			mi_lo_rand = mk_row(base,"imagery",  b, e,"random","Low"),
			mi_lo_rep  = mk_row(base,"imagery",  b, e,"repeated","Low"),
			mi_hi_rand = mk_row(base,"imagery",  b, e,"random","High"),
			mi_hi_rep  = mk_row(base,"imagery",  b, e,"repeated","High")
		)
		Xs <- lapply(rows, X_of)
		l <- (Xs$om_lo_rand + Xs$om_lo_rep)/2 - (Xs$om_hi_rand + Xs$om_hi_rep)/2 -
			((Xs$mi_lo_rand + Xs$mi_lo_rep)/2 - (Xs$mi_hi_rand + Xs$mi_hi_rep)/2)
		as.numeric(l)

	} else if (spec$type == "pmb_gxe") {
		b <- "beta"
		rows <- list(
			# OM during
			om_dur_lo_rand = mk_row(base, "physical", b, "during", "random",   "Low"),
			om_dur_lo_rep  = mk_row(base, "physical", b, "during", "repeated", "Low"),
			om_dur_hi_rand = mk_row(base, "physical", b, "during", "random",   "High"),
			om_dur_hi_rep  = mk_row(base, "physical", b, "during", "repeated", "High"),
			# OM after
			om_aft_lo_rand = mk_row(base, "physical", b, "after",  "random",   "Low"),
			om_aft_lo_rep  = mk_row(base, "physical", b, "after",  "repeated", "Low"),
			om_aft_hi_rand = mk_row(base, "physical", b, "after",  "random",   "High"),
			om_aft_hi_rep  = mk_row(base, "physical", b, "after",  "repeated", "High"),
			# MI during
			mi_dur_lo_rand = mk_row(base, "imagery",  b, "during", "random",   "Low"),
			mi_dur_lo_rep  = mk_row(base, "imagery",  b, "during", "repeated", "Low"),
			mi_dur_hi_rand = mk_row(base, "imagery",  b, "during", "random",   "High"),
			mi_dur_hi_rep  = mk_row(base, "imagery",  b, "during", "repeated", "High"),
			# MI after
			mi_aft_lo_rand = mk_row(base, "imagery",  b, "after",  "random",   "Low"),
			mi_aft_lo_rep  = mk_row(base, "imagery",  b, "after",  "repeated", "Low"),
			mi_aft_hi_rand = mk_row(base, "imagery",  b, "after",  "random",   "High"),
			mi_aft_hi_rep  = mk_row(base, "imagery",  b, "after",  "repeated", "High")
		)
		Xs <- lapply(rows, X_of)

		# Mid ≈ 0.5*Low + 0.5*High per (g,e), then (After−During)_OM − (After−During)_MI
		OM_dur <- 0.5*((Xs$om_dur_lo_rand + Xs$om_dur_lo_rep)/2 + (Xs$om_dur_hi_rand + Xs$om_dur_hi_rep)/2)
		OM_aft <- 0.5*((Xs$om_aft_lo_rand + Xs$om_aft_lo_rep)/2 + (Xs$om_aft_hi_rand + Xs$om_aft_hi_rep)/2)
		MI_dur <- 0.5*((Xs$mi_dur_lo_rand + Xs$mi_dur_lo_rep)/2 + (Xs$mi_dur_hi_rand + Xs$mi_dur_hi_rep)/2)
		MI_aft <- 0.5*((Xs$mi_aft_lo_rand + Xs$mi_aft_lo_rep)/2 + (Xs$mi_aft_hi_rand + Xs$mi_aft_hi_rep)/2)

		l <- (OM_aft - OM_dur) - (MI_aft - MI_dur)
		as.numeric(l)

	} else stop("Unknown spec$type")
}

# ---- Hypothesis manifest ----
manifest <- tibble::tribble(
	~hid, ~type,     ~band,   ~epoch,
	"h1", "acc_gxe", "theta", "during",
	"h2", "acc_gxe", "alpha", "after",
	"h3", "pmb_gxe", NA_character_, NA_character_
) %>%
	mutate(
		title = case_when(
			hid == "h1" ~ "H1: (Low-High)_OM - (Low-High)_MI @ Theta, During (avg rep)",
			hid == "h2" ~ "H2: (Low-High)_OM - (Low-High)_MI @ Alpha, After (avg rep)",
			hid == "h3" ~ "H3: (After-During)_OM - (After-During)_MI @ Beta (mid≈0.5·Low+0.5·High, avg rep)"
		)
	)

# ---- Coefs, covariance, L ----
beta_hat <- coef(gam_full)
V        <- vcov(gam_full, unconditional = TRUE)

L_list <- list()
for (i in seq_len(nrow(manifest))) {
	h <- manifest[i, ]
	spec <- list(type = h$type, band = h$band, epoch = h$epoch)
	L_h <- do.call(rbind, lapply(split(chan_grid, chan_grid$chan),
								 function(one) contrast_row(one, spec)))
	rownames(L_h) <- paste(h$hid, chan_grid$chan, sep = "::")
	L_list[[h$hid]] <- L_h
}
L_stack <- do.call(rbind, L_list)

# Observed stats
obs <- list()
for (hid in names(L_list)) {
	L <- L_list[[hid]]
	c_hat <- as.numeric(L %*% beta_hat)
	Sigma <- L %*% V %*% t(L)
	se    <- sqrt(pmax(diag(Sigma), 0))
	t_obs <- ifelse(se > 0, c_hat / se, 0)
	obs[[hid]] <- list(c_hat = c_hat, se = se, t = t_obs, Sigma = Sigma)
}

# Global null across H×sensors
Sigma_stack <- L_stack %*% V %*% t(L_stack)
if (use_t_null) {
	if (is.null(df_t)) df_t <- max(1, floor(gam_full$df.residual))
	message("Using multivariate-t null with df = ", df_t)
	Z <- mvtnorm::rmvt(B, sigma = Sigma_stack, df = df_t, type = "shifted")
} else {
	message("Using multivariate normal null.")
	Z <- mvtnorm::rmvnorm(B, sigma = Sigma_stack)
}
se_stack <- sqrt(pmax(diag(Sigma_stack), 0))
Tstar    <- sweep(Z, 2, se_stack, "/")
maxabs   <- apply(abs(Tstar), 1, max, na.rm = TRUE)
thr_glob <- as.numeric(quantile(maxabs, probs = 1 - alpha_fwer, na.rm = TRUE))
message(sprintf("Global max-T threshold across H×sensors = %.3f", thr_glob))

# ---- Plotting helpers ----
scale_to_0range <- function(x, range = 8) {
	x <- x - min(x); x <- x / max(x); (x - .5) * range
}
plot_grid <- chan_grid %>%
	mutate(
		x_plot = lat * cos(long * (pi/180)),
		y_plot = lat * sin(long * (pi/180)),
		x_scaled = scale_to_0range(x_plot, 8),
		y_scaled = scale_to_0range(y_plot, 8.5)
	)

# ---- Save per-H tables/plots ----
for (hid in names(obs)) {
	h <- manifest %>% filter(hid == !!hid) %>% slice(1)
	c_hat <- obs[[hid]]$c_hat
	se    <- obs[[hid]]$se
	t_obs <- obs[[hid]]$t
	nm    <- rownames(L_list[[hid]])
	chans <- sub("^.*::", "", nm)

	pmax_glob <- sapply(abs(t_obs), function(tt) mean(maxabs >= tt))
	sig_glob  <- abs(t_obs) >= thr_glob

	tab <- tibble(
		chan = chans,
		effect_c_hat = c_hat,
		se           = se,
		t_obs        = t_obs,
		pmax_global  = pmax_glob,
		sig_FWER05_global = sig_glob
	) %>% arrange(desc(abs(t_obs)))
	write.csv(tab, file.path(out_tables, paste0("maxT_results_", hid, ".csv")), row.names = FALSE)

	dfp <- tibble(
		chan = chans, c_hat = c_hat, se = se, t_obs = t_obs, sig = sig_glob
	) %>% left_join(plot_grid, by = "chan")

	lim <- if (scale_by_ci95) {
		lo <- quantile(dfp$c_hat, 0.025, na.rm = TRUE)
		hi <- quantile(dfp$c_hat, 0.975, na.rm = TRUE)
		max(abs(c(lo, hi)))
	} else max(abs(dfp$c_hat), na.rm = TRUE)
	if (!is.finite(lim) || lim <= 0) lim <- 1

	p <- ggplot(dfp, aes(x = x_scaled, y = y_scaled)) +
		geom_point(aes(fill = c_hat), shape = 21, size = 6, stroke = 0.5, color = "grey30") +
		geom_point(data = dfp %>% filter(sig),
				   aes(x = x_scaled, y = y_scaled),
				   shape = 21, size = 8, stroke = 1.1, color = "black", fill = NA) +
		{ if (add_chan_labels)
			geom_text(data = dfp, aes(label = chan), nudge_y = 0.35, size = 2.6, color = "grey20")
			else NULL } +
		scale_fill_gradient2(name = "c_hat", limits = c(-lim, lim), midpoint = 0) +
		coord_equal() +
		theme_minimal(base_size = 12) +
		theme(axis.title = element_blank(), axis.text = element_blank(),
			  axis.ticks = element_blank(), panel.grid = element_blank(),
			  legend.position = "right",
			  plot.title = element_text(face = "bold")) +
		labs(title = h$title,
			 subtitle = sprintf("Global FWER max-T threshold = %.3f", thr_glob))

	ggsave(file.path(out_plots, paste0("maxT_topo_", hid, ".pdf")), p, width = 7, height = 6)
	message("Saved: ", file.path(out_tables, paste0("maxT_results_", hid, ".csv")),
			" ; ", file.path(out_plots, paste0("maxT_topo_", hid, ".pdf")))
}

message("\nDone.")
