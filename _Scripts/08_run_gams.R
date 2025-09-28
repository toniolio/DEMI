###################
### 08_run_gams ###
###################

# Run the GAMs! See methods section in paper.

suppressPackageStartupMessages({
	library(tidyverse)
	library(mgcv)
})

message("\n=== 08_run_gams (GP baseline + factor-by cell deviations) ===")

# -------------------- User choices --------------------
path_rds <- "_Scripts/_rds/"
drop_final_block <- TRUE
drop_ref_chans   <- TRUE
ref_chans        <- c("M1","M2")

# Runtime/detail knobs
n_threads     <- 1L
discrete_fit  <- TRUE     # set FALSE if you truly want longer (exact) fits
k_xy          <- 80L      # upper limit for GP basis per surface; bump to 100 if you want
include_rep_in_cell <- FALSE  # TRUE -> deviations also differ by rep (doubles the number of surfaces)

# Accuracy bins
acc_q <- c(0.20, 0.80)

# -------------------- Helpers --------------------
read_rds_safe <- function(p) { if (!file.exists(p)) stop("Missing file: ", p); readRDS(p) }
safe_set_halfsum <- function (f) { n <- nlevels(f); if (n > 1) contr.sum(n)*0.5 else NULL }
size_summary <- function(g) tibble(
	n_coef   = length(coef(g)),
	n_smooth = length(g$smooth),
	smooth   = sapply(g$smooth, `[[`, "label"),
	bs_dim   = sapply(g$smooth, `[[`, "bs.dim")
)

# 1) Load GAM-ready data -------------------------------------------------------
dat <- read_rds_safe(file.path(path_rds, "dat_gam.rds")) %>% as_tibble()

req <- c("powerdb","participant","group","band","epoch","rep","trial","chan","lat","long","accuracy")
miss <- setdiff(req, names(dat)); if (length(miss)) stop("dat_gam.rds missing: ", paste(miss, collapse=", "))

if (drop_final_block && "block" %in% names(dat)) {
	dat <- dat %>% filter(block < max(block, na.rm = TRUE))
}
if (drop_ref_chans) {
	dat <- dat %>% filter(!(chan %in% ref_chans))
}

# Factors / contrasts
dat <- dat %>%
	mutate(
		participant = factor(participant),
		group  = factor(group, levels = c("physical","imagery")),
		band   = factor(band,  levels = c("theta","alpha","beta")),
		epoch  = factor(epoch, levels = c("during","after")),
		rep    = factor(rep,   levels = c("random","repeated")),
		chan   = factor(chan)
	) %>% droplevels()

for (nm in c("group","band","epoch","rep")) {
	contrasts(dat[[nm]]) <- safe_set_halfsum(dat[[nm]])
}

core_na <- sapply(dat[, c("powerdb","accuracy","group","band","epoch","rep","lat","long")], \(x) sum(is.na(x)))
if (any(core_na > 0)) { print(core_na); stop("NAs present in core fields; fix upstream.") }

message("\nParticipants (by group):")
print(dat %>% distinct(group, participant) %>% count(group, name = "n_participants"))
message("\nData size after filters:")
print(dat %>% summarise(rows = n(), sensors = n_distinct(chan), participants = n_distinct(participant)))

# 2) Planar coords for GP (same transform used later in inference) ------------
xy <- transmute(dat,
				x = lat * cos(long * (pi/180)),
				y = lat * sin(long * (pi/180))
)
rmax <- max(sqrt(xy$x^2 + xy$y^2))
dat  <- bind_cols(dat, mutate(xy, x = x / rmax, y = y / rmax))

eps <- 1e-8
r   <- sqrt(dat$x^2 + dat$y^2)
if (any(r > 1 + 1e-6)) stop("Some (x,y) lie outside unit disk; check lat/long scaling.")
if (any(r >= 1)) {
	dat <- dat %>% mutate(
		x = if_else(r >= 1, x * (1 - eps), x),
		y = if_else(r >= 1, y * (1 - eps), y)
	)
}

# 3) Accuracy bins -------------------------------------------------------------
q <- quantile(dat$accuracy, acc_q)
dat <- dat %>%
	mutate(acc_bin = case_when(
		accuracy <= q[1] ~ "Low",
		accuracy >= q[2] ~ "High",
		TRUE ~ NA_character_
	)) %>%
	filter(!is.na(acc_bin)) %>%
	mutate(acc_bin = factor(acc_bin, levels = c("Low","High")))
saveRDS(list(q = q, acc_q = acc_q), file.path(path_rds, "accuracy_bin_cutpoints.rds"))

# 4) Deviation structure as a FACTOR 'cell' -----------------------------------
#    (group × band × epoch × acc_bin); optional: include rep for more nuance.
if (include_rep_in_cell) {
	dat <- dat %>%
		mutate(cell = interaction(group, band, epoch, rep, acc_bin, drop = TRUE))
} else {
	dat <- dat %>%
		mutate(cell = interaction(group, band, epoch,       acc_bin, drop = TRUE))
}
message("Number of cell levels: ", nlevels(dat$cell))

# 5) Formula: baseline GP + factor-by GP deviations + participant RE ----------
n_participants <- nlevels(dat$participant)
k_participant  <- max(1L, n_participants - 1L)

form_full <- as.formula(paste(
	"powerdb ~",
	"group * band * epoch * rep +",
	# Baseline spatial field (GP)
	"s(x, y, bs='gp', k=", k_xy, ") +",
	# One deviation GP per cell (FACTOR-BY; each gets its own λ because no 'id=')
	"s(x, y, bs='gp', k=", k_xy, ", by=cell) +",
	# Participant random intercept
	"s(participant, bs='re', k=", k_participant, ")"
))

message("\n--- Fitting FULL (GP baseline + by=cell deviations) ---")
t_full <- system.time(
	gam_full <- mgcv::bam(
		formula   = form_full,
		data      = dat,
		method    = "fREML",
		discrete  = discrete_fit,   # set FALSE for slower / exact fits
		nthreads  = n_threads,
		select    = TRUE,           # redundant surfaces shrink away
		gc.level  = 0
	)
)
print(t_full)

message("\n--- Summarizing ---")
t_sum <- system.time(gam_full_summary <- summary(gam_full))
print(t_sum)

saveRDS(gam_full,           file.path(path_rds, "gam_full.rds"))
saveRDS(gam_full_summary,   file.path(path_rds, "gam_full_summary.rds"))
saveRDS(gam_full$sp,        file.path(path_rds, "gam_full_sp.rds"))

message("\nSize summary:")
print(size_summary(gam_full))

# Optional deeper checks:
# mgcv::gam.check(gam_full)      # k-index etc. (can be slow)
# print(mgcv::concurvity(gam_full, full = FALSE))

message("\n--- Done. Saved to: ", path_rds, " ---")
