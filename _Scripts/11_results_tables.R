#######################################
### 11_results_tables.R (summaries) ###
#######################################
# create a main table for manuscript

suppressPackageStartupMessages({
	library(tidyverse)
	library(readr)
})

message("\n=== 11_results_tables.R ===")

path_rds   <- "_Scripts/_rds/"
out_tables <- "_Tables/"

# ---------- Helpers ----------
pairs_lr <- tribble(
	~L, ~R,
	"Fp1","Fp2","F7","F8","F3","F4","FT7","FT8","FC3","FC4","T7","T8",
	"C3","C4","TP7","TP8","CP3","CP4","P3","P4","P7","P8","O1","O2"
)

laterality_stats <- function(tab_effects) {
	# tab_effects: data.frame with cols chan, effect_c_hat
	x <- tab_effects %>% select(chan, effect_c_hat)
	lr <- pairs_lr %>%
		rowwise() %>%
		mutate(
			cL = x$effect_c_hat[match(L, x$chan)],
			cR = x$effect_c_hat[match(R, x$chan)],
			diff = cL - cR
		) %>% ungroup()
	tibble(
		mean_LminusR = mean(lr$diff, na.rm = TRUE),
		sd_LminusR   = sd(lr$diff,   na.rm = TRUE)
	)
}

# Tiny channel metadata for hemisphere/midline counts (32-quickcap style)
hemi_map <- tibble(
	chan = c("Fp1","F7","F3","FT7","FC3","T7","C3","TP7","CP3","P3","P7","O1",
			 "Fp2","F8","F4","FT8","FC4","T8","C4","TP8","CP4","P4","P8","O2",
			 "Fz","FCz","Cz","CPz","Pz","Oz"),
	hemi = c(rep("L",12), rep("R",12), rep("M",6))
)

# ---------- Table 1: Main hypotheses (FWER) ----------
manifest_main <- tribble(
	~hid, ~label,
	"h1", "H1: (Low-High)_OM - (Low-High)_MI @ Theta, During",
	"h2", "H2: (Low-High)_OM - (Low-High)_MI @ Alpha, After",
	"h3", "H3: (After-During)_OM - (After-During)_MI @ Beta (mid acc)"
)

# read threshold from RDS if present
thr_glob <- NA_real_
rds_file <- file.path(path_rds, "maxT_results_global.rds")
if (file.exists(rds_file)) {
	allr <- readRDS(rds_file)
	# allr is a list by hid; thresholds should match across hids, take first
	if (is.list(allr) && length(allr) > 0 && !is.null(allr[[1]]$thr_global))
		thr_glob <- allr[[1]]$thr_global
}

summ_main <- purrr::map_dfr(manifest_main$hid, function(hid) {
	f <- file.path(out_tables, paste0("maxT_results_", hid, ".csv"))
	if (!file.exists(f)) stop("Missing: ", f)
	tab <- read_csv(f, show_col_types = FALSE)

	# significant set and counts by hemisphere
	sig <- tab %>% filter(sig_FWER05_global)
	hemi_counts <- sig %>%
		left_join(hemi_map, by = "chan") %>%
		count(hemi, name = "n") %>%
		tidyr::complete(hemi = c("L","R","M"), fill = list(n = 0)) %>%
		pivot_wider(names_from = hemi, values_from = n, names_prefix = "n_")

	# peak by absolute effect
	peak <- tab %>% slice_max(abs(effect_c_hat), n = 1, with_ties = FALSE)

	# laterality quick stat on full map
	lat <- laterality_stats(tab)

	tibble(
		hid = hid,
		label = manifest_main$label[manifest_main$hid == hid],
		thr_global = thr_glob,
		n_sig = nrow(sig),
		n_total = nrow(tab),
		prop_sig = nrow(sig)/nrow(tab),
		peak_sensor = peak$chan,
		peak_c_hat  = peak$effect_c_hat,
		peak_t      = peak$t_obs,
		mean_LminusR = lat$mean_LminusR,
		sd_LminusR   = lat$sd_LminusR,
		n_left  = hemi_counts$n_L,
		n_mid   = hemi_counts$n_M,
		n_right = hemi_counts$n_R
	)
})

# make it pretty order
summ_main <- summ_main %>%
	select(hid, label, thr_global, n_sig, n_total, prop_sig,
		   peak_sensor, peak_c_hat, peak_t,
		   mean_LminusR, sd_LminusR, n_left, n_mid, n_right)

write.csv(summ_main, file.path(out_tables, "main_hypotheses_summary.csv"), row.names = FALSE)
message("Wrote: ", file.path(out_tables, "main_hypotheses_summary.csv"))
print(summ_main, n = nrow(summ_main))

# ---------- Table 2: Post-hoc maps (descriptive) ----------
# match ids used in 12_posthoc_group_maps.R
manifest_posthoc <- tribble(
	~id, ~label,
	"OM_theta_during_LH", "OM: (Low-High) @ Theta, During",
	"MI_theta_during_LH", "MI: (Low-High) @ Theta, During",
	"OM_alpha_after_LH",  "OM: (Low-High) @ Alpha, After",
	"MI_alpha_after_LH",  "MI: (Low-High) @ Alpha, After",
	"OM_beta_PMBR",       "OM: (After-During) @ Beta (mid acc)",
	"MI_beta_PMBR",       "MI: (After-During) @ Beta (mid acc)"
)

posthoc_cutoff <- 2.0  # keep in sync with 12_posthoc_group_maps.R

one_posthoc <- function(id) {
	f <- file.path(out_tables, paste0("posthoc_", id, ".csv"))
	if (!file.exists(f)) stop("Missing: ", f)
	tab <- read_csv(f, show_col_types = FALSE)  # columns: chan, c_hat, se, t_obs, sig
	n_ring <- sum(abs(tab$t_obs) >= posthoc_cutoff, na.rm = TRUE)

	top3 <- tab %>%
		arrange(desc(abs(t_obs))) %>%
		slice_head(n = 3) %>%
		transmute(top = paste0(chan, " (ĉ=", sprintf("%.3f", c_hat),
							   ", t=", sprintf("%.2f", t_obs), ")")) %>%
		pull(top)

	tibble(
		id = id,
		n_ge_tcut = n_ring,
		top1 = top3[1] %||% NA_character_,
		top2 = top3[2] %||% NA_character_,
		top3 = top3[3] %||% NA_character_
	)
}

summ_posthoc <- purrr::map_dfr(manifest_posthoc$id, one_posthoc) %>%
	left_join(manifest_posthoc, by = "id") %>%
	select(id, label, n_ge_tcut, top1, top2, top3)

write.csv(summ_posthoc, file.path(out_tables, "posthoc_maps_summary.csv"), row.names = FALSE)
message("Wrote: ", file.path(out_tables, "posthoc_maps_summary.csv"))
print(summ_posthoc, n = nrow(summ_posthoc))

message("\nDone.")
