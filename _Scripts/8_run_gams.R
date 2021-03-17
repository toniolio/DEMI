#####################
### DEMI Run GAMs ###
#####################

## Detecting Errors in Motor Imagery (DEMI) ##
## RUN GAMS ##

library(tidyverse)
library(mgcv)
library(brms)

path <- "_Scripts/_rds/"

# make a column "accuracy" ("high" vs. "low") from the upper and lower quantities of accuracy_rating

#### train PP model ####

if(!file.exists(paste(path, "gam_PP_train_summary.rds", sep=""))){
	dat_PP_train = readRDS(paste0(path, "dat_PP_train.rds"))
	dat_PP_train$band <- as.factor(dat_PP_train$band)
	dat_PP_train$epoch <- as.factor(dat_PP_train$epoch)
	dat_PP_train$rep <- as.factor(dat_PP_train$rep)
	(
		dat_PP_train
		%>% dplyr::mutate(
			block = dplyr::case_when(
				(trial>=1) & (trial<=20) ~ 1
				, (trial>=21) & (trial<=40) ~ 2
				, (trial>=41) & (trial<=60) ~ 3
				, (trial>=61) & (trial<=80) ~ 4
				, (trial>=81) & (trial<=100) ~ 5
				, (trial>=101) & (trial<=120) ~ 6
			)
		)
	) -> dat_PP_train
	system.time(
		gam_PP_train <-
			mgcv::bam(
				formula = powerdb ~ (
					band * epoch * rep + t2(
						lat , long , accuracy , block
						, by=interaction( band , epoch , rep )
									  , d = c(2,1,1)
									  , bs = c("sos","cr","cr")
									  , k = c(32,10,5)
									  )
					)
				, data = dat_PP_train
				, method = "fREML"
				, discrete = TRUE
				, nthreads = floor(parallel::detectCores()/3)
				, gc.level = 0
			)
	)
	# system.time(
	# 	gam_PP_train_brms <- brms::brm(
	# 		formula = brms::bf(
	# 			powerdb ~ (
	# 				band * epoch * rep + t2(
	# 					lat , long , accuracy , block
	# 					, by=interaction( band , epoch , rep )
	# 					, d = c(2,1,1)
	# 					, bs = c("sos","cr","cr")
	# 					, k = c(31,9,5)
	# 					)
	# 				)
	# 			)
	# 			, data = dat_PP_train
	# 			, silent = F
	# 			, refresh = 20
	# 			, iter = 2000
	# 			, chains = floor(parallel::detectCores()/3)
	# 			, cores = floor(parallel::detectCores()/3)
	# 		)
	# 	)
	saveRDS(gam_PP_train,file=paste0(path, "gam_PP_train.rds"))
	system.time(
		gam_PP_train_summary <- summary(gam_PP_train)
	)
	saveRDS(gam_PP_train_summary,file=paste0(path, "gam_PP_train_summary.rds"))
	rm(list=ls())
	gc()
}

#### train MI model ####

if(!file.exists(paste(path, "gam_MI_train_summary.rds", sep=""))){
	dat_MI_train = readRDS(paste0(path, "dat_MI_train.rds"))
	dat_MI_train$band <- as.factor(dat_MI_train$band)
	dat_MI_train$epoch <- as.factor(dat_MI_train$epoch)
	dat_MI_train$rep <- as.factor(dat_MI_train$rep)
	(
		dat_MI_train
		%>% dplyr::mutate(
			block = dplyr::case_when(
				(trial>=1) & (trial<=20) ~ 1
				, (trial>=21) & (trial<=40) ~ 2
				, (trial>=41) & (trial<=60) ~ 3
				, (trial>=61) & (trial<=80) ~ 4
				, (trial>=81) & (trial<=100) ~ 5
				, (trial>=101) & (trial<=120) ~ 6
			)
		)
	) -> dat_MI_train
	system.time(
		gam_MI_train <-
			mgcv::bam(
				formula = powerdb ~ (
					band * epoch * rep + t2(
						lat , long , accuracy , block
						, by=interaction( band , epoch , rep )
						, d = c(2,1,1)
						, bs = c("sos","cr","cr")
						, k = c(32,10,5)
					)
				)
				, data = dat_MI_train
				, method = "fREML"
				, discrete = TRUE
				, nthreads = floor(parallel::detectCores()/3)
				, gc.level = 0
			)
	)
	# system.time(
	# 	gam_MI_train_brms <- brms::brm(
	# 		formula = brms::bf(
	# 			powerdb ~ (
	# 				band * epoch * rep + t2(
	# 					lat , long , accuracy , block
	# 					, by=interaction( band , epoch , rep )
	# 					, d = c(2,1,1)
	# 					, bs = c("sos","cr","cr")
	# 					, k = c(31,9,5)
	# 					)
	# 				)
	# 			)
	# 			, data = dat_MI_train
	# 			, silent = F
	# 			, refresh = 20
	# 			, iter = 2000
	# 			, chains = floor(parallel::detectCores()/3)
	# 			, cores = floor(parallel::detectCores()/3)
	# 		)
	# 	)
	saveRDS(gam_MI_train,file=paste0(path, "gam_MI_train.rds"))
	system.time(
		gam_MI_train_summary <- summary(gam_MI_train)
	)
	saveRDS(gam_MI_train_summary,file=paste0(path, "gam_MI_train_summary.rds"))
	rm(list=ls())
	gc()
}
