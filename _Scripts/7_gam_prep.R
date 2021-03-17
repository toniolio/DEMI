#############################
### Prepare data for GAMs ###
#############################

# written by Tony Ingram & Mike Lawrence
# correspondence: tony.ingram@dal.ca

library(tidyverse)

path <- "_Scripts/_rds/"

if(!file.exists(paste(path, "dat_MI_test.rds", sep=""))){

	#### prepare data for GAMs ####

	if(!file.exists(paste(path, "dat_bands_notime.rds", sep=""))){

		# load dat:
		dat <- as_tibble(readRDS(paste0(path, "all_dat.rds")))
		# dat <- as_tibble(readRDS(paste0(path, "all_dat_first_participant.rds")))

		# reminders: group = group assignment (0 = physical; 1 = imagery),
		# condition = whether a trial is physical vs. imagery. Note: recall
		# that the physical group performed only the physical task for all 6
		# blocks of 20 trials; while the imagery group performed the imagery
		# task for 5 blocks, then a final block of the physical task.

		## remove bad channels? ##

		# # drop bad channels
		# bads <- eeg_info$interpolated[subject] # as determined using PREP algorithm
		# bads <- as.list(strsplit(bads, " ")[[1]])
		# epoched.long <- epoched.long[!(epoched.long$chan %in% bads),]
		# rm(bads)

		## average power across bands ##

		# average both frequencies and time to produce 1 power value per trial
		# for each frequency band of interest (theta, alpha, beta)
		(
			dat
			%>% dplyr::select(-power)
			%>% dplyr::mutate(
				band = dplyr::case_when(
					#theta (4-8), alpha (9-12), and beta (13-30) bands
					(freq>=4) & (freq<=8) ~ 'theta'
					, (freq>=9) & (freq<=12) ~ 'alpha'
					, (freq>=13) & (freq<=30) ~ 'beta'
					#all others will be NA
				)
			)
			%>% dplyr::select(-freq) # don't need freq anymore
			%>% dplyr::filter(
				!is.na(band)
				, time >= 0
				, time <= 1
			)
			%>% dplyr::group_by(
				band # computing a mean power across time and frequency in each band
				# also need to group by all the variables we want to KEEP
				, participant
				, group
				, trial
				, condition
				, rep
				, epoch
				, chan
				, complexity
				, avg_velocity
				, error
				, vresp
				, accuracy_rating
				, lat
				, long
			)
			%>% dplyr::summarize(
				powerdb = mean(powerdb)
				, .groups = 'drop'
			)
		) -> dat
		saveRDS(dat,file=paste(path, "dat_bands_notime.rds", sep=""))
	}else{
		dat = readRDS(paste(path, "dat_bands_notime.rds", sep=""))
	}

	#### scale accuracy rating ####

	dat <- (
		dat
		%>% dplyr::mutate(
			accuracy = scale(accuracy_rating)[,1]
		)
	)

	#### set up training and test sets ####

	# separate PP vs MI data sets
	dat_PP <-
		dat %>% dplyr::filter(
			group == 'physical'
		)
	dat_MI <-
		dat %>% dplyr::filter(
			group == 'imagery'
			, condition != 'physical' # do not train the imagery model on physical trials
		)

	# set up training and testing data
	test_proportion = .2
	set.seed(1) #for reproducibility
	dat_PP =
		(
			dat_PP
			#first collapse to a list of trials per participant
			%>% dplyr::group_by(
				participant
				, trial
			)
			%>% dplyr::summarise(.groups='keep')
			# modify trial number to be guaranteed sequential 1:num_trials_tot
			%>% dplyr::ungroup(trial)
			%>% dplyr::arrange(participant,trial)
			%>% dplyr::mutate(
				num_trials_tot = n()
				, num_test = floor(num_trials_tot*test_proportion)
				, test_trial_list = list(sample(unique(trial),num_test))
				, set = case_when(
					trial %in% unlist(test_trial_list) ~ 'test'
					, TRUE ~ 'train'
				)
			)
			# %>% View())
			# check the %s as expected (see console print)
			%>% (function(x){
				print(mean(x$set=='test'))
				return(x)
			})
			# join back with original data
			%>% dplyr::right_join(
				dat_PP
				, by = c('participant','trial')
			)
			%>% (function(x){
				print(mean(x$set=='test'))
				return(x)
			})
		)
	dat_MI =
		(
			dat_MI
			# first collapse to a list of trials per participant
			%>% dplyr::group_by(
				participant
				, trial
			)
			%>% dplyr::summarise(.groups='keep')
			# modify trial number to be guaranteed sequential 1:num_trials_tot
			%>% dplyr::ungroup(trial)
			%>% dplyr::arrange(participant,trial)
			%>% dplyr::mutate(
				num_trials_tot = n()
				, num_test = floor(num_trials_tot*test_proportion)
				, test_trial_list = list(sample(unique(trial),num_test))
				, set = case_when(
					trial %in% unlist(test_trial_list) ~ 'test'
					, TRUE ~ 'train'
				)
			)
			# %>% View())
			# check the %s as expected (see console print)
			%>% (function(x){
				print(mean(x$set=='test'))
				return(x)
			})
			# join back with original data
			%>% dplyr::right_join(
				dat_MI
				, by = c('participant','trial')
			)
			%>% (function(x){
				print(mean(x$set=='test'))
				return(x)
			})
		)
	# separate training and testing sets for each group:
	dat_PP_train <- (
		dat_PP
		%>% dplyr::filter(
			set == 'train'
		)
	)
	saveRDS(dat_PP_train, paste(path, "dat_PP_train.rds", sep=""))

	dat_PP_test <- (
		dat_PP
		%>% dplyr::filter(
			set == 'test'
		)
	)
	saveRDS(dat_PP_test, paste(path, "dat_PP_test.rds", sep=""))

	dat_MI_train <- (
		dat_MI
		%>% dplyr::filter(
			set == 'train'
		)
	)
	saveRDS(dat_MI_train, paste(path, "dat_MI_train.rds", sep=""))

	dat_MI_test <- (
		dat_MI
		%>% dplyr::filter(
			set == 'test'
		)
	)
	saveRDS(dat_MI_test, paste(path, "dat_MI_test.rds", sep=""))
}
