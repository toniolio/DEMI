######################
### DEMI gam tests ###
######################

library(tidyverse)
library(mgcv)

path <- "_Scripts/_rds/"

dat_PP_train = readRDS(paste0(path, "dat_PP_train.rds"))
dat_MI_train = readRDS(paste0(path, "dat_MI_train.rds"))
dat_PP_test = readRDS(paste0(path, "dat_PP_test.rds"))
dat_MI_test = readRDS(paste0(path, "dat_MI_test.rds"))
# gam_PP = readRDS(paste0(path, "gam_PP.rds"))
# gam_MI = readRDS(paste0(path, "gam_MI.rds"))
# gam_PP_summary = readRDS(paste0(path, "gam_PP_summary.rds"))
# gam_MI_summary = readRDS(paste0(path, "gam_MI_summary.rds"))

gam = readRDS(paste0(path, "gam.rds"))
# gam_summary = readRDS(paste0(path, "gam_summary.rds"))

# sapply(gam_PP$smooth, "[[",  "label")
# sapply(gam_MI$smooth, "[[",  "label")

dat <- rbind(dat_PP_train
			 , dat_PP_test
			 , dat_MI_train
			 , dat_MI_test
			 )
rm(dat_PP_train
   , dat_PP_test
   , dat_MI_train
   , dat_MI_test
   )
(
	dat
	%>% dplyr::filter(
		block != max(dat$block) # drop 6th block
	)
) -> dat

halfsum_contrasts = function (...) contr.sum(...) * 0.5
contrasts(dat$group) = halfsum_contrasts
contrasts(dat$band) = halfsum_contrasts
contrasts(dat$epoch) = halfsum_contrasts
contrasts(dat$rep) = halfsum_contrasts
dat$block = dat$block-3

#### separate groups; mgcv version; no subject 're' ####

# (
# 	dat # has both training and testing datasets, indicated by column `set`
# 	%>% dplyr::ungroup()
# 	%>% dplyr::mutate(
# 		PPpreds =
# 			(
# 				mgcv::predict.bam(
# 					gam_PP
# 					, newdata = dat
# 					, discrete = F
# 					# , n.threads = parallel::detectCores()
# 				)
# 			)
# 		, MIpreds =
# 			(
# 				mgcv::predict.bam(
# 					gam_MI
# 					, newdata = dat
# 					, discrete = F
# 					# , n.threads = parallel::detectCores()
# 				)
# 			)
# 	)
# 	%>% tidyr::pivot_longer(
# 		cols = c(PPpreds,MIpreds)
# 		, names_to = 'model'
# 		, values_to = 'preds'
# 	)
# 	%>% dplyr::group_by(
# 		model
# 		, group
# 		, set
# 	)
# 	%>% dplyr::mutate(
# 		model_squared_error = (powerdb-preds)^2
# 		, null_squared_error = (powerdb-mean(powerdb))^2
# 	)
# 	%>% dplyr::summarise(
# 		model_sse = sum(model_squared_error)
# 		, null_sse = sum(null_squared_error)
# 		, model_null_ratio = model_sse/null_sse
# 		, prop_accounted = 1-model_null_ratio
# 	)
# )

# NOTE: as a sanity check, the prop_accounted number for models using it's own
# training data should show the same number as summary() "Deviance explained".

#### single model; mgcv version; no subject 're' ####

(
	dat # has both training and testing datasets, indicated by column `set`
	%>% dplyr::ungroup()
	%>% dplyr::mutate(
		preds =
			(
				mgcv::predict.bam(
					gam
					, newdata = dat
					, discrete = F
					, n.threads = parallel::detectCores()
					# , exclude = "s(participant)"
				)
			)
	)
	%>% dplyr::group_by(
		set
	)
	%>% dplyr::mutate(
		model_squared_error = (powerdb-preds)^2
		, null_squared_error = (powerdb-mean(powerdb))^2
	)
	%>% dplyr::summarise(
		model_sse = sum(model_squared_error)
		, null_sse = sum(null_squared_error)
		, model_null_ratio = model_sse/null_sse
		, prop_accounted = 1-model_null_ratio
	)
)
