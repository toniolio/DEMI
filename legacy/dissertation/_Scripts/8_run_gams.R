#####################
### DEMI Run GAMs ###
#####################

## Detecting Errors in Motor Imagery (DEMI) ##
## RUN GAMS ##

library(tidyverse)
library(mgcv)
# library(brms)

# set path for each model below.

# #### train PP only model ####
#
# if(!file.exists(paste(path, "gam_PP_summary.rds", sep=""))){
#
# 	# load training data
# 	path <- "_Scripts/_rds/"
# 	dat_PP_train = readRDS(paste0(path, "dat_PP_train.rds"))
#
# 	#-------------------------------#
# 	# mgcv version; no subject 're' #
# 	#-------------------------------#
#
# 	# system.time(
# 	# 	gam_PP <-
# 	# 		mgcv::bam(
# 	# 			formula = powerdb ~ (
# 	# 				band * epoch * rep + te(
# 	# 					lat , long , accuracy , block
# 	# 					, by=interaction( band , epoch , rep )
# 	# 					, d = c(2,1,1)
# 	# 					, bs = c("sos","cr","cr")
# 	# 					, k = c(31 # k = 32 locations - 1
# 	# 							,9 # k = 10 possible responses - 1
# 	# 							,5) # k = imagery only has 5 blocks - 1
# 	# 					)
# 	# 				)
# 	# 			, data = dat_PP_train
# 	# 			, method = "fREML"
# 	# 			, discrete = TRUE
# 	# 			, nthreads = floor(parallel::detectCores())
# 	# 			, gc.level = 0
# 	# 		)
# 	# )
# 	# saveRDS(gam_PP,file=paste0(path, "gam_PP.rds"))
# 	# system.time(
# 	# 	gam_PP_summary <- summary(gam_PP)
# 	# )
# 	# saveRDS(gam_PP_summary,file=paste0(path, "gam_PP_summary.rds"))
# 	# gam_PP_summary
# 	# rm(list=ls())
# 	# gc()
#
# 	#---------------------------------#
# 	# mgcv version; with subject 're' #
# 	#---------------------------------#
#
# 	# system.time(
# 	# 	gam_PPre <-
# 	# 		mgcv::bam(
# 	# 			formula = powerdb ~ (
# 	# 				band * epoch * rep + te(
# 	# 					lat , long , accuracy , block
# 	# 					, by=interaction( band , epoch , rep )
# 	# 					, d = c(2,1,1)
# 	# 					, bs = c("sos","cr","cr")
# 	# 					, k = c(31 # k = 32 locations - 1
# 	# 							,9 # k = 10 possible responses - 1
# 	# 							,4) # k = imagery only has 5 blocks - 1
# 	# 					) + t2(
# 	# 						lat , long , accuracy , block, participant
# 	# 						, by=interaction( band , epoch , rep )
# 	# 						, d = c(2,1,1,1)
# 	# 						, bs = c("sos","cr","cr","re")
# 	# 						, k = c(31 # k = 32 locations - 1
# 	# 								,9 # k = 10 possible responses - 1
# 	# 								,4 # k = imagery only has 5 blocks - 1
# 	# 								,38) # k = N of smaller group (imagery) - 1
# 	# 						, full = TRUE
# 	# 					)
# 	# 				)
# 	# 			, data = dat_PP_train
# 	# 			, method = "fREML"
# 	# 			, discrete = TRUE
# 	# 			, nthreads = floor(parallel::detectCores())
# 	# 			, gc.level = 0
# 	# 		)
# 	# )
# 	# saveRDS(gam_PPre,file=paste0(path, "gam_PPre.rds"))
# 	# system.time(
# 	# 	gam_PPre_summary <- summary(gam_PPre)
# 	# )
# 	# gam_PPre_summary
# 	# saveRDS(gam_PPre_summary,file=paste0(path, "gam_PPre_summary.rds"))
# 	# rm(list=ls())
# 	# gc()
#
# 	#-------------------------------#
# 	# brms version; no subject 're' #
# 	#-------------------------------#
#
# 	# system.time(
# 	# 	gam_PP_brms <- brms::brm(
# 	# 		formula = brms::bf(
# 	# 			powerdb ~ (
# 	# 				band * epoch * rep + t2(
# 	# 					lat , long , accuracy , block
# 	# 					, by=interaction( band , epoch , rep )
# 	# 					, d = c(2,1,1)
# 	# 					, bs = c("sos","cr","cr")
# 	# 					, k = c(31 # k = 32 locations - 1
# 	# 							,9 # k = 10 possible responses - 1
# 	# 							,5) # k = imagery only has 5 blocks - 1
# 	# 					)
# 	# 				)
# 	# 			)
# 	# 			, data = dat_PP_train
# 	# 			, silent = F
# 	# 			, refresh = 20
# 	# 			, iter = 2000
# 	# 			, chains = floor(parallel::detectCores())
# 	# 			, cores = floor(parallel::detectCores())
# 	# 		)
# 	# 	)
# 	# saveRDS(gam_PP_brms,file=paste0(path, "gam_PP_brms.rds"))
# 	# system.time(
# 	# 	gam_PP_brms_summary <- summary(gam_PP_brms)
# 	# )
# 	# saveRDS(gam_PP_brms_summary,file=paste0(path, "gam_PP_brms_summary.rds"))
# 	# rm(list=ls())
# 	# gc()
# }
#
# #### train MI only model ####
#
# if(!file.exists(paste(path, "gam_MI_summary.rds", sep=""))){
#
# 	# load training data
# 	path <- "_Scripts/_rds/"
# 	dat_MI_train = readRDS(paste0(path, "dat_MI_train.rds"))
#
# 	#-------------------------------#
# 	# mgcv version; no subject 're' #
# 	#-------------------------------#
#
# 	# system.time(
# 	# 	gam_MI <-
# 	# 		mgcv::bam(
# 	# 			formula = powerdb ~ (
# 	# 				band * epoch * rep + te(
# 	# 					lat , long , accuracy , block
# 	# 					, by=interaction( band , epoch , rep )
# 	# 					, d = c(2,1,1)
# 	# 					, bs = c("sos","cr","cr")
# 	# 					, k = c(31 # k = 32 locations - 1
# 	# 							,9 # k = 10 possible responses - 1
# 	# 							,5) # k = imagery only has 5 blocks - 1
# 	# 				)
# 	# 			)
# 	# 			, data = dat_MI_train
# 	# 			, method = "fREML"
# 	# 			, discrete = TRUE
# 	# 			, nthreads = floor(parallel::detectCores())
# 	# 			, gc.level = 0
# 	# 		)
# 	# )
# 	# saveRDS(gam_MI,file=paste0(path, "gam_MI.rds"))
# 	# system.time(
# 	# 	gam_MI_summary <- summary(gam_MI)
# 	# )
# 	# saveRDS(gam_MI_summary,file=paste0(path, "gam_MI_summary.rds"))
# 	# gam_MI_summary
# 	# rm(list=ls())
# 	# gc()
#
# 	#---------------------------------#
# 	# mgcv version; with subject 're' #
# 	#---------------------------------#
#
# 	# system.time(
# 	# 	gam_MIre <-
# 	# 		mgcv::bam(
# 	# 			formula = powerdb ~ (
# 	# 				band * epoch * rep + te(
# 	# 					lat , long , accuracy , block
# 	# 					, by=interaction( band , epoch , rep )
# 	# 					, d = c(2,1,1)
# 	# 					, bs = c("sos","cr","cr")
# 	# 					, k = c(31 # k = 32 locations - 1
# 	# 							,9 # k = 10 possible responses - 1
# 	# 							,4) # k = imagery only has 5 blocks - 1
# 	# 				) + t2(
# 	# 					lat , long , accuracy , block, participant
# 	# 					, by=interaction( band , epoch , rep )
# 	# 					, d = c(2,1,1,1)
# 	# 					, bs = c("sos","cr","cr","re")
# 	# 					, k = c(31 # k = 32 locations - 1
# 	# 							,9 # k = 10 possible responses - 1
# 	# 							,4 # k = imagery only has 5 blocks - 1
# 	# 							,38) # k = N of smaller group (imagery) - 1
# 	# 					, full = TRUE
# 	# 				)
# 	# 			)
# 	# 			, data = dat_MI_train
# 	# 			, method = "fREML"
# 	# 			, discrete = TRUE
# 	# 			, nthreads = floor(parallel::detectCores())
# 	# 			, gc.level = 0
# 	# 		)
# 	# )
# 	# saveRDS(gam_MIre,file=paste0(path, "gam_MIre.rds"))
# 	# system.time(
# 	# 	gam_MIre_summary <- summary(gam_MIre)
# 	# )
# 	# gam_MIre_summary
# 	# saveRDS(gam_MIre_summary,file=paste0(path, "gam_MIre_summary.rds"))
# 	# rm(list=ls())
# 	# gc()
#
# 	#-------------------------------#
# 	# brms version; no subject 're' #
# 	#-------------------------------#
#
# 	# system.time(
# 	# 	gam_MI_brms <- brms::brm(
# 	# 		formula = brms::bf(
# 	# 			powerdb ~ (
# 	# 				band * epoch * rep + t2(
# 	# 					lat , long , accuracy , block
# 	# 					, by=interaction( band , epoch , rep )
# 	# 					, d = c(2,1,1)
# 	# 					, bs = c("sos","cr","cr")
# 	# 					, k = c(31 # k = 32 locations - 1
# 	# 							,9 # k = 10 possible responses - 1
# 	# 							,5) # k = imagery only has 5 blocks - 1
# 	# 					)
# 	# 				)
# 	# 			)
# 	# 			, data = dat_MI_train
# 	# 			, silent = F
# 	# 			, refresh = 20
# 	# 			, iter = 2000
# 	# 			, chains = floor(parallel::detectCores())
# 	# 			, cores = floor(parallel::detectCores())
# 	# 		)
# 	# 	)
# 	# saveRDS(gam_MI_brms,file=paste0(path, "gam_MI_brms.rds"))
# 	# system.time(
# 	# 	gam_MI_brms_summary <- summary(gam_MI_brms)
# 	# )
# 	# saveRDS(gam_MI_brms_summary,file=paste0(path, "gam_MI_brms_summary.rds"))
# 	# rm(list=ls())
# 	# gc()
# }

#### run both group model ####

# # load training data only
# path <- "_Scripts/_rds/"
# dat_PP_train = readRDS(paste0(path, "dat_PP_train.rds"))
# dat_MI_train = readRDS(paste0(path, "dat_MI_train.rds"))
#
# dat_train <- rbind(dat_PP_train
# 			 , dat_MI_train
# )
# rm(dat_PP_train,dat_MI_train)

# load all data
path <- "_Scripts/_rds/"
dat_train = readRDS(paste0(path, "dat_gam.rds"))

# consider moving all of below to previous script

# important: ensure "geographic" coordinates
dat_train$lat = 90 - dat_train$lat

# dat_train$rep <- as.factor(ifelse(dat_train$rep==0, "random", "repeated"))
# dat_train$epoch <- as.factor(ifelse(dat_train$epoch==1, "during", "after"))
halfsum_contrasts = function (...) contr.sum(...) * 0.5
contrasts(dat_train$group) = halfsum_contrasts
contrasts(dat_train$condition) = halfsum_contrasts
contrasts(dat_train$band) = halfsum_contrasts
contrasts(dat_train$epoch) = halfsum_contrasts
contrasts(dat_train$rep) = halfsum_contrasts
dat_train$block = dat_train$block-3

# remove final block:
dat_train <- (
	dat_train
	%>% dplyr::filter(
		block < max(block)
	)
)
# consider removing M1 and M2 as these are typically just reference sensors
dat_train <- (
	dat_train
	%>% dplyr::filter(
		chan != 'M1'
		, chan != 'M2'
	)
)

#### Participant Characterization ####

# NOTE from "3_behav_analysis.R":
bdat <- readRDS("_Scripts/_rds/bdat2.rds")
unique(sort(bdat$participant))
# lost due to experimenter error (wrong experiment): 9, 10, 12, 14, 20,
# didn't complete experiment (tech issue): 89, 96, 100
# lost due to very bad EEG: 24
# id skipped: 13, 26, 36, 78

# then for further eeg analysis:
unique(sort(dat_train$participant))
# experimenter error (eeg not recorded): 6,7,11
# technical error (triggers missing): 54, 56
# technical error (software crash): 65
# technical error (VEOL did not record): 86

# In total we lose 7 more due to EEG related technical errors. For the
# behavioral analysis, we had 96 recruited, and 9 dropped due to technical
# issues with the experimental setup.

# Therefore in the paper we describe that we recruited 96 and lost 16 due to
# technical errors. Total remianing = 80.

# remove these from bdat:
bdat <- (
	bdat
	%>% dplyr::filter(
		!(participant) %in% c(6,7,11,54,56,65,86)
	)
)

unique(sort(subset(bdat, group =='physical')$participant))
unique(sort(subset(bdat, group =='imagery')$participant))

# age:

print("mean (sd) age for all participants:")
mean(bdat[!duplicated(bdat$participant),]$age)
sd(bdat[!duplicated(bdat$participant),]$age)

print("mean (sd) age for each group:")
bdat[!duplicated(bdat$participant),] %>%
	group_by(group) %>%
	summarise( count = length(unique(participant))
			   , agemean = mean(age)
			   , agesd = sd(age)
	)

# bio sex

bdat[!duplicated(bdat$participant),] %>%
	group_by(group) %>%
	count(sex)

# handedness

bdat[!duplicated(bdat$participant),] %>%
	group_by(group) %>%
	count(handedness)



#### run model ####

#---------------------------------#
# full version; with subject 're' #
#---------------------------------#

system.time(
	gam_re_1 <- mgcv::bam(
		formula = powerdb ~ (
			group * band * epoch * rep + te(
				lat , long , accuracy , block
				, by=interaction( group , band , epoch , rep )
				, d = c(2,1,1)
				, bs = c("sos","cr","cr")
				, k = c(29 # k = 30 locations - 1
						,9 # k = 10 possible responses - 1
						,4) # k = imagery only has 5 blocks - 1
			) + s(
				participant
				, bs = "re"
				, k = 79 # k = 80 participants - 1
			)
		)
		, data = dat_train
		, method = "fREML"
		, discrete = TRUE
		, nthreads = floor(parallel::detectCores())
		, gc.level = 0
	)
)
saveRDS(gam_re_1,file=paste0(path, "gam_re_1.rds"))
system.time(
	gam_re_1_summary <- summary(gam_re_1)
)
gam_re_1_summary
saveRDS(gam_re_1_summary,file=paste0(path, "gam_re_1_summary.rds"))
rm(gam_re_1, gam_re_1_summary)
gc()

#----------------------------------#
# without block; with subject 're' #
#----------------------------------#

system.time(
	gam_re_2 <- mgcv::bam(
		formula = powerdb ~ (
			group * band * epoch * rep + te(
				lat , long , accuracy
				, by=interaction( group , band , epoch , rep )
				, d = c(2,1)
				, bs = c("sos","cr")
				, k = c(29 # k = 30 locations - 1
						,9) # k = 10 possible responses - 1
			) + s(
				participant
				, bs = "re"
				, k = 79 # k = 80 participants - 1
			)
		)
		, data = dat_train
		, method = "fREML"
		, discrete = TRUE
		, nthreads = floor(parallel::detectCores())
		, gc.level = 0
	)
)
saveRDS(gam_re_2,file=paste0(path, "gam_re_2.rds"))
system.time(
	gam_re_2_summary <- summary(gam_re_2)
)
gam_re_2_summary
saveRDS(gam_re_2_summary,file=paste0(path, "gam_re_2_summary.rds"))
rm(gam_re_2, gam_re_2_summary)
gc()

#-----------------------------------------#
# without block or rep; with subject 're' #
#-----------------------------------------#

system.time(
	gam_re_3 <- mgcv::bam(
		formula = powerdb ~ (
			group * band * epoch + te(
				lat , long , accuracy
				, by=interaction( group , band , epoch )
				, d = c(2,1)
				, bs = c("sos","cr")
				, k = c(29 # k = 30 locations - 1
						,9) # k = 10 possible responses - 1
			) + s(
				participant
				, bs = "re"
				, k = 79 # k = 80 participants - 1
			)
		)
		, data = dat_train
		, method = "fREML"
		, discrete = TRUE
		, nthreads = floor(parallel::detectCores())
		, gc.level = 0
	)
)
saveRDS(gam_re_3,file=paste0(path, "gam_re_3.rds"))
system.time(
	gam_re_3_summary <- summary(gam_re_3)
)
gam_re_3_summary
saveRDS(gam_re_3_summary,file=paste0(path, "gam_re_3_summary.rds"))
rm(gam_re_3, gam_re_3_summary)
gc()

#---------------------------------------------------------#
# without block or rep; with subject 're' and subject 'ti #
#---------------------------------------------------------#

system.time(
	gam_re_4 <- mgcv::bam(
		formula = powerdb ~ (
			group * band * epoch + te(
				lat , long , accuracy
				, by=interaction( group , band , epoch )
				, d = c(2,1)
				, bs = c("sos","cr")
				, k = c(29 # 30 locations - 1
						,9 # 10 possible responses - 1
						)
			) + s(
				participant
				, bs = "re"
				, k = 79 # 80 participants - 1
			) + ti(
				lat , long , accuracy, participant
				# , by=interaction( group , band , epoch )
				, d = c(2,1,1)
				, bs = c("sos","cr","re")
				, k = c(29 # 30 locations - 1
						,9 # 10 possible responses - 1
						,79 # 80 participants - 1
						)
			)
		)
		, data = dat_train
		, method = "fREML"
		, discrete = TRUE
		, nthreads = floor(parallel::detectCores())
		, gc.level = 0
	)
)
saveRDS(gam_re_4,file=paste0(path, "gam_re_4.rds"))
system.time(
	gam_re_4_summary <- summary(gam_re_4)
)
gam_re_4_summary
saveRDS(gam_re_4_summary,file=paste0(path, "gam_re_4_summary.rds"))
rm(gam_re_4, gam_re_4_summary)
gc()

#--------------------------------------------------#
# without block; with subject 're' and subject 'ti #
#--------------------------------------------------#

system.time(
	gam_re_5 <- mgcv::bam(
		formula = powerdb ~ (
			group * band * epoch * rep + te(
				lat , long , accuracy
				, by=interaction( group , band , epoch, rep )
				, d = c(2,1)
				, bs = c("sos","cr")
				, k = c(29 # 30 locations - 1
						,9 # 10 possible responses - 1
				)
			) + s(
				participant
				, bs = "re"
				, k = 79 # 80 participants - 1
			) + ti(
				lat , long , accuracy, participant
				# , by=interaction( group , band , epoch , rep )
				, d = c(2,1,1)
				, bs = c("sos","cr","re")
				, k = c(29 # 30 locations - 1
						,9 # 10 possible responses - 1
						,79 # 80 participants - 1
				)
			)
		)
		, data = dat_train
		, method = "fREML"
		, discrete = TRUE
		, nthreads = floor(parallel::detectCores())
		, gc.level = 0
	)
)
saveRDS(gam_re_5,file=paste0(path, "gam_re_5.rds"))
system.time(
	gam_re_5_summary <- summary(gam_re_5)
)
gam_re_5_summary
saveRDS(gam_re_5_summary,file=paste0(path, "gam_re_5_summary.rds"))
rm(gam_re_5, gam_re_5_summary)
gc()
