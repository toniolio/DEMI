#####################
### DEMI Run GAMs ###
#####################

## Detecting Errors in Motor Imagery (DEMI) ##
## RUN GAMS ##

library(tidyverse)
library(mgcv)
# library(brms)

# set path for each model below.

#### train PP model ####

if(!file.exists(paste(path, "gam_PP_summary.rds", sep=""))){

	# load training data
	path <- "_Scripts/_rds/"
	dat_PP_train = readRDS(paste0(path, "dat_PP_train.rds"))

	#-------------------------------#
	# mgcv version; no subject 're' #
	#-------------------------------#

	# system.time(
	# 	gam_PP <-
	# 		mgcv::bam(
	# 			formula = powerdb ~ (
	# 				band * epoch * rep + te(
	# 					lat , long , accuracy , block
	# 					, by=interaction( band , epoch , rep )
	# 					, d = c(2,1,1)
	# 					, bs = c("sos","cr","cr")
	# 					, k = c(31 # k = 32 locations - 1
	# 							,9 # k = 10 possible responses - 1
	# 							,5) # k = imagery only has 5 blocks - 1
	# 					)
	# 				)
	# 			, data = dat_PP_train
	# 			, method = "fREML"
	# 			, discrete = TRUE
	# 			, nthreads = floor(parallel::detectCores())
	# 			, gc.level = 0
	# 		)
	# )
	# saveRDS(gam_PP,file=paste0(path, "gam_PP.rds"))
	# system.time(
	# 	gam_PP_summary <- summary(gam_PP)
	# )
	# saveRDS(gam_PP_summary,file=paste0(path, "gam_PP_summary.rds"))
	# gam_PP_summary
	# rm(list=ls())
	# gc()

	#---------------------------------#
	# mgcv version; with subject 're' #
	#---------------------------------#

	# system.time(
	# 	gam_PPre <-
	# 		mgcv::bam(
	# 			formula = powerdb ~ (
	# 				band * epoch * rep + te(
	# 					lat , long , accuracy , block
	# 					, by=interaction( band , epoch , rep )
	# 					, d = c(2,1,1)
	# 					, bs = c("sos","cr","cr")
	# 					, k = c(31 # k = 32 locations - 1
	# 							,9 # k = 10 possible responses - 1
	# 							,4) # k = imagery only has 5 blocks - 1
	# 					) + t2(
	# 						lat , long , accuracy , block, participant
	# 						, by=interaction( band , epoch , rep )
	# 						, d = c(2,1,1,1)
	# 						, bs = c("sos","cr","cr","re")
	# 						, k = c(31 # k = 32 locations - 1
	# 								,9 # k = 10 possible responses - 1
	# 								,4 # k = imagery only has 5 blocks - 1
	# 								,38) # k = N of smaller group (imagery) - 1
	# 						, full = TRUE
	# 					)
	# 				)
	# 			, data = dat_PP_train
	# 			, method = "fREML"
	# 			, discrete = TRUE
	# 			, nthreads = floor(parallel::detectCores())
	# 			, gc.level = 0
	# 		)
	# )
	# saveRDS(gam_PPre,file=paste0(path, "gam_PPre.rds"))
	# system.time(
	# 	gam_PPre_summary <- summary(gam_PPre)
	# )
	# gam_PPre_summary
	# saveRDS(gam_PPre_summary,file=paste0(path, "gam_PPre_summary.rds"))
	# rm(list=ls())
	# gc()

	#-------------------------------#
	# brms version; no subject 're' #
	#-------------------------------#

	# system.time(
	# 	gam_PP_brms <- brms::brm(
	# 		formula = brms::bf(
	# 			powerdb ~ (
	# 				band * epoch * rep + t2(
	# 					lat , long , accuracy , block
	# 					, by=interaction( band , epoch , rep )
	# 					, d = c(2,1,1)
	# 					, bs = c("sos","cr","cr")
	# 					, k = c(31 # k = 32 locations - 1
	# 							,9 # k = 10 possible responses - 1
	# 							,5) # k = imagery only has 5 blocks - 1
	# 					)
	# 				)
	# 			)
	# 			, data = dat_PP_train
	# 			, silent = F
	# 			, refresh = 20
	# 			, iter = 2000
	# 			, chains = floor(parallel::detectCores())
	# 			, cores = floor(parallel::detectCores())
	# 		)
	# 	)
	# saveRDS(gam_PP_brms,file=paste0(path, "gam_PP_brms.rds"))
	# system.time(
	# 	gam_PP_brms_summary <- summary(gam_PP_brms)
	# )
	# saveRDS(gam_PP_brms_summary,file=paste0(path, "gam_PP_brms_summary.rds"))
	# rm(list=ls())
	# gc()
}

#### train MI model ####

if(!file.exists(paste(path, "gam_MI_summary.rds", sep=""))){

	# load training data
	path <- "_Scripts/_rds/"
	dat_MI_train = readRDS(paste0(path, "dat_MI_train.rds"))

	#-------------------------------#
	# mgcv version; no subject 're' #
	#-------------------------------#

	# system.time(
	# 	gam_MI <-
	# 		mgcv::bam(
	# 			formula = powerdb ~ (
	# 				band * epoch * rep + te(
	# 					lat , long , accuracy , block
	# 					, by=interaction( band , epoch , rep )
	# 					, d = c(2,1,1)
	# 					, bs = c("sos","cr","cr")
	# 					, k = c(31 # k = 32 locations - 1
	# 							,9 # k = 10 possible responses - 1
	# 							,5) # k = imagery only has 5 blocks - 1
	# 				)
	# 			)
	# 			, data = dat_MI_train
	# 			, method = "fREML"
	# 			, discrete = TRUE
	# 			, nthreads = floor(parallel::detectCores())
	# 			, gc.level = 0
	# 		)
	# )
	# saveRDS(gam_MI,file=paste0(path, "gam_MI.rds"))
	# system.time(
	# 	gam_MI_summary <- summary(gam_MI)
	# )
	# saveRDS(gam_MI_summary,file=paste0(path, "gam_MI_summary.rds"))
	# gam_MI_summary
	# rm(list=ls())
	# gc()

	#---------------------------------#
	# mgcv version; with subject 're' #
	#---------------------------------#

	# system.time(
	# 	gam_MIre <-
	# 		mgcv::bam(
	# 			formula = powerdb ~ (
	# 				band * epoch * rep + te(
	# 					lat , long , accuracy , block
	# 					, by=interaction( band , epoch , rep )
	# 					, d = c(2,1,1)
	# 					, bs = c("sos","cr","cr")
	# 					, k = c(31 # k = 32 locations - 1
	# 							,9 # k = 10 possible responses - 1
	# 							,4) # k = imagery only has 5 blocks - 1
	# 				) + t2(
	# 					lat , long , accuracy , block, participant
	# 					, by=interaction( band , epoch , rep )
	# 					, d = c(2,1,1,1)
	# 					, bs = c("sos","cr","cr","re")
	# 					, k = c(31 # k = 32 locations - 1
	# 							,9 # k = 10 possible responses - 1
	# 							,4 # k = imagery only has 5 blocks - 1
	# 							,38) # k = N of smaller group (imagery) - 1
	# 					, full = TRUE
	# 				)
	# 			)
	# 			, data = dat_MI_train
	# 			, method = "fREML"
	# 			, discrete = TRUE
	# 			, nthreads = floor(parallel::detectCores())
	# 			, gc.level = 0
	# 		)
	# )
	# saveRDS(gam_MIre,file=paste0(path, "gam_MIre.rds"))
	# system.time(
	# 	gam_MIre_summary <- summary(gam_MIre)
	# )
	# gam_MIre_summary
	# saveRDS(gam_MIre_summary,file=paste0(path, "gam_MIre_summary.rds"))
	# rm(list=ls())
	# gc()

	#-------------------------------#
	# brms version; no subject 're' #
	#-------------------------------#

	# system.time(
	# 	gam_MI_brms <- brms::brm(
	# 		formula = brms::bf(
	# 			powerdb ~ (
	# 				band * epoch * rep + t2(
	# 					lat , long , accuracy , block
	# 					, by=interaction( band , epoch , rep )
	# 					, d = c(2,1,1)
	# 					, bs = c("sos","cr","cr")
	# 					, k = c(31 # k = 32 locations - 1
	# 							,9 # k = 10 possible responses - 1
	# 							,5) # k = imagery only has 5 blocks - 1
	# 					)
	# 				)
	# 			)
	# 			, data = dat_MI_train
	# 			, silent = F
	# 			, refresh = 20
	# 			, iter = 2000
	# 			, chains = floor(parallel::detectCores())
	# 			, cores = floor(parallel::detectCores())
	# 		)
	# 	)
	# saveRDS(gam_MI_brms,file=paste0(path, "gam_MI_brms.rds"))
	# system.time(
	# 	gam_MI_brms_summary <- summary(gam_MI_brms)
	# )
	# saveRDS(gam_MI_brms_summary,file=paste0(path, "gam_MI_brms_summary.rds"))
	# rm(list=ls())
	# gc()
}

#### run both group model ####

# load training data
path <- "_Scripts/_rds/"
dat_PP_train = readRDS(paste0(path, "dat_PP_train.rds"))
dat_MI_train = readRDS(paste0(path, "dat_MI_train.rds"))

dat_train <- rbind(dat_PP_train
			 , dat_MI_train
)
rm(dat_PP_train,dat_MI_train)

# consider moving this to previous script

# dat_train$rep <- as.factor(ifelse(dat_train$rep==0, "random", "repeated"))
# dat_train$epoch <- as.factor(ifelse(dat_train$epoch==1, "during", "after"))
halfsum_contrasts = function (...) contr.sum(...) * 0.5
contrasts(dat_train$group) = halfsum_contrasts
contrasts(dat_train$band) = halfsum_contrasts
contrasts(dat_train$epoch) = halfsum_contrasts
contrasts(dat_train$rep) = halfsum_contrasts
dat_train$block = dat_train$block-3

#-------------------------------#
# mgcv version; no subject 're' #
#-------------------------------#

system.time(
	gam <- mgcv::bam(
		formula = powerdb ~ (
			group * band * epoch * rep + te(
				lat , long , accuracy , block
				, by=interaction( group , band , epoch , rep )
				, d = c(2,1,1)
				, bs = c("sos","cr","cr")
				, k = c(31 # k = 32 locations - 1
						,9 # k = 10 possible responses - 1
						,5) # k = imagery only has 5 blocks - 1
				)
			)
		, data = dat_train
		, method = "fREML"
		, discrete = TRUE
		, nthreads = floor(parallel::detectCores())
		, gc.level = 0
		)
	)
saveRDS(gam,file=paste0(path, "gam.rds"))
system.time(
	gam_summary <- summary(gam)
)
gam_summary
saveRDS(gam_summary,file=paste0(path, "gam_summary.rds"))
rm(list=ls())
gc()

#-------------------------------#
# brms version; no subject 're' #
#-------------------------------#

# system.time(
# 	gam_brms <- brms::brm(
# 		formula = powerdb ~ (
# 			group * band * epoch * rep + t2(
# 				lat , long , accuracy , block
# 				, by=interaction( group , band , epoch , rep )
# 				, d = c(2,1,1)
# 				, bs = c("sos","cr","cr")
# 				, k = c(31 # k = 32 locations - 1
# 						,9 # k = 10 possible responses - 1
# 						,5) # k = imagery only has 5 blocks - 1
# 				)
# 			)
# 		, data = dat_train
# 		, silent = F
# 		, refresh = 20
# 		, iter = 2000
# 		, chains = floor(parallel::detectCores())
# 		, cores = floor(parallel::detectCores())
# 		)
# 	)
# saveRDS(gam_brms,file=paste0(path, "gam_brms.rds"))
# system.time(
# 	gam_brms_summary <- summary(gam_brms)
# )
# saveRDS(gam_brms_summary,file=paste0(path, "gam_brms_summary.rds"))
# rm(list=ls())
# gc()
