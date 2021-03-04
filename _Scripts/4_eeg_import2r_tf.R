	## Detecting Errors in Motor Imagery (DEMI) ##
			## import processed EEG data ##
			## perform wavelet transform ##
		## and combine with behavioural data ##

# written by Tony Ingram & Mike Lawrence
# correspondence: tony.ingram@dal.ca

# # useful code:
# rm(list=setdiff(ls(), c())) # clear environment
# graphics.off() # clear figures
# cat("\014") # clear console

# code for timing this script:
ptm <- proc.time()

library(tidyverse)
library(eeguana)
library(biwavelet)
# library(multidplyr)

#### Load behavioural data and set up ####

# load cleaned behavioural data (bad trials removed):
bdat <- readRDS("_Scripts/_rds/bdat2.rds")

# important reminder: group = true group assignment,
# condition = the task condition (physical vs. imagery)
# where the imagery group will have trials of each,
# because they perform physically in the final block.

max_mt <- max(bdat$mt) # in seconds

#### load in EEG data, downsampling, get frequency data ####

# # read in info from PREP algorithm:
# eeg_info <- read.csv("prep_inf_inter_csd.csv")

# get eeg file names
path <- "_Data/eeg/"
filenames <- dir(path, pattern=".edf")

# what if you already ran this and some participants are already done?
if(file.exists("_Scripts/_rds/alldat.rds")){
	alldat <- readRDS("_Scripts/_rds/alldat.rds")
	ps <- unique(as.numeric(levels(alldat$participant))[alldat$participant])
	donefiles <- sprintf("sub-%03d_eeg_prepped.edf",ps)
	filenames <- setdiff(filenames,donefiles)
}

# for testing:
# filenames <- c(filenames[2], filenames[5], filenames[length(filenames)])

# define function for importing and preprocessing EDFs
edf2epochs <- function(f, event, end = NA, window=c(-1000, 1000)) {
	# Actually read in EDF data
	eeg <- read_edf(f, event_ch = "EDF Annotations")
	# Find practice trials using figure animation times
	annot_by_trial <- as_tibble(eeg$.events) %>%
		mutate(
			bad_start = .description == "red_on" & lag(.description) != "stim_on",
			trial = cumsum(.description == "stim_on" | bad_start)
		)
	fig_durations <- annot_by_trial %>%
		select(c(trial, .initial, .description)) %>%
		spread(key = .description, value = .initial) %>%
		mutate(fig_duration = red_on - stim_on)
	practice_trials <- fig_durations %>%
		filter(fig_duration > 3000 & !is.na(fig_duration))
	# Drop practice trial annotations from events table
	practice_annots <- annot_by_trial %>%
		filter(trial %in% practice_trials$trial)
	eeg$.events <- subset(
		eeg$.events, !(.initial %in% practice_annots$.initial)
	)
	# Drop EMG and EOG channels from data
	emg <- eeg %>% select(EMG.L, EMG.A) # make a backup
	eeg <- eeg %>% select(-HEO, -VEO, -EMG.L, -EMG.A)
	# Epoch eeg data
	if (is.na(end)) {
		epoched <- eeg %>%
			eeg_segment(.description == event, lim = window, unit = "ms")
	} else {
		epoched <- eeg %>%
			eeg_segment(.description == event, end = .description == end)
	}
	return(epoched)
}

# define a function to process a given level
get_wt = function(x){
	(
		x
		%>% dplyr::select(time,eeg)
		%>% as.matrix()
		%>% biwavelet::wt()
		%>% (function(x){
			(
				x$power
				%>% t()
				%>% tibble::as_tibble()
				%>% dplyr::mutate(
					time = x$t
				)
				%>% tidyr::pivot_longer(
					-time
					, names_to = 'period_index'
					, names_prefix = 'V'
				)
				%>% dplyr::mutate(
					period = x$period[as.numeric(period_index)]
				)
				%>% dplyr::select(-period_index)
			)
		})
	)
}

for (i in 1:length(filenames)) {
	# identify subject and get path to data
	subject <- -parse_number(filenames[i])
	subjectstr <- sprintf("%03d", subject)
	subjectpath <- paste(path, "sub-", subjectstr,"_eeg_prepped.edf", sep="")

	print(paste("Participant:", subject))

	# load eeg data
	print(paste("... loading and epoching EEG"))
	# Notes: .sample = 1 sample / millisecond, and .id = trial; longest trials
	# is about 8 seconds, so take 1.5 seconds extra, and baseline period
	# is -2000 to -1000. Then add 1000ms on each end for edge effects of wt():
	epoched <- edf2epochs(subjectpath, "trace_start", window = c(-3000, 10500))
	epoched <- epoched$.signal # units: mV/m^2 (if CSD was performed)

	# # OLD but potentially useful:
	# epoch.rest <- edf2epochs(subjectpath, "stim_on", end = "trace_start")
	# epoch.rest <- epoch.rest$.signal
	# epoch.task <- edf2epochs(subjectpath, "trace_start", end = "trace_end")
	# epoch.task <- epoch.task$.signal
	#
	# # create column to easily identify task vs. rest, and combine
	# epoch.rest$epoch <- rep(0, nrow(epoch.rest))
	# epoch.task$epoch <- rep(1, nrow(epoch.task))
	# epoched <- rbind(epoch.rest, epoch.task)
	# # fix "time" column (here named .sample)
	# epoched <- (
	# 	epoched
	# 	%>% dplyr::group_by(.id) # .id = trial
	# 	%>% dplyr::mutate(
	# 		.sample = 1:length(.sample)
	# 		)
	# 	)
	# rm(epoch.rest, epoch.task)
	# gc()

	# downsample time round 1 (before wavelet decomposition)
	# note, assuming Nyquist of 2.56, and the highest frequency band we are
	# interested in for this experiment is Beta (upper range 31Hz), we will
	# downsample to a sampling rate of now lower than ~ 80Hz.
	downsample_by_factor = 10 # gives 100Hz, because we sampled at 1000Hz originally
	print(paste("... downsampling round 1"))
	epoched =
		(
			epoched
			#keep only Nth sample
			# `==1` ensures that sample_nm==1 is kept, which is what you want to do
			# when the sample count starts at 1. If there is a sample_num==0, you'd
			# want to have `==0`
			%>% dplyr::filter(
				(.sample %% downsample_by_factor) == 0
			)
			#re-express sample number
			%>% dplyr::mutate(
				#again, this assumes that sample numbering starts at 1
				.sample = (.sample/downsample_by_factor)
			)
		)
	# make sample number represent time (in seconds)
	epoched$.sample <- (epoched$.sample*downsample_by_factor)/1000

	# # if no downsampling, convert ms to seconds:
	# epoched$.sample <- epoched$.sample/1000

	# check new sample rate (also look at actual time values for gut check)
	# length(unique(epoched$.sample))/(max(unique(epoched$.sample))-min(unique(epoched$.sample)))

	# make eeg long
	epoched.long <- gather(epoched
						   , chan
						   , eeg
						   , c(-.id,-.sample)
						   , factor_key=TRUE)
	names(epoched.long)[names(epoched.long) == '.id'] <- 'trial' # more intuitive
	names(epoched.long)[names(epoched.long) == '.sample'] <- 'time' # more intuitive

	# TEMPORARY:
	epoched.long$eeg <- epoched.long$eeg/1000

	# # drop bad channels
	# bads <- eeg_info$interpolated[subject] # as determined using PREP algorithm
	# bads <- as.list(strsplit(bads, " ")[[1]])
	# epoched.long <- epoched.long[!(epoched.long$chan %in% bads),]
	# rm(bads)

	rm(epoched) # free up memory
	gc()

	# # depending on how you epoched, you may have to shift time so that active
	# # epoch starts at 0 seconds (rest is negative time):
	# epoched.long <- (
	# 	epoched.long
	# 	%>% dplyr::group_by(trial, chan)
	# 	%>% dplyr::arrange(time)
	# 	%>% dplyr::mutate(
	# 		time = time - min(time[epoch==1])
	# 	)
	# 	%>% dplyr::ungroup()
	# )
	#
	# # drop epoch column (don't need anymore)
	# epoched.long <- (
	# 	epoched.long
	# 	%>% dplyr::select(-epoch)
	# )

	# wavelet decomposition to get time-freq info
	# new units will be: (mV/m^2)^2/Hz
	print(paste("... performing wavelet decomposition"))
	dat_wt =
		(
			epoched.long
			%>% dplyr::group_by(trial, chan)
			%>% filter(n() >= 10) # reject very short trials
			%>% dplyr::summarise(
				get_wt(cur_data())
				, .groups = 'drop'
			)
		)

	dat_wt =
		(
			dat_wt
			%>% dplyr::full_join(
				(
					epoched.long
					#group by same groups as when computing dat_wt
					%>% dplyr::group_by(trial,chan)
					#grabs the first row in each group
					%>% dplyr::slice(1)
					#toss the time column
					%>% dplyr::select(-time)
					#ungroup (just in case)
					%>% dplyr::ungroup()
				)
			)
		)
	# note: eeg column appears not to have been added back properly

	rm(epoched.long) # free up memory
	gc()

	# rename output of wt() to "power"
	names(dat_wt)[names(dat_wt) == 'value'] <- 'power' # (mV/m^2)^2/Hz

	# change period to Hz to align with EEG literature (note: be sure time is in seconds)
	dat_wt <- (
		dat_wt
		%>% dplyr::mutate(
			freq = 1/period
		)
		%>% dplyr::select(-period)
	)

	# # check some things:
	# range(unique(dat_wt$freq))
	# length(unique(dat_wt$freq))
	# length(unique(dat_wt$time))/(max(unique(dat_wt$time))-min(dat_wt$time))
	# range(dat_wt$time)

	# # plot a channel and a trial — how does it look?
	# (
	# 	dat_wt
	# 	%>% dplyr::filter(
	# 		chan == "C3"
	# 		, trial == 2
	# 	)
	# 	%>% ggplot()
	# 	+ geom_tile(
	# 		aes(
	# 			x = time
	# 			, y = as.factor(freq) # yuk
	# 			, fill = power
	# 		)
	# 	)
	# )

	# decibel normalize the power data
	# resulting units will be: (mV/m^2)^2/Hz (dB)
	# using baseline period before active part of trial
	print(paste("... computing decibel normalization"))
	dat_wt <- (
		dat_wt
		%>% dplyr::group_by(trial, chan, freq)
		%>% dplyr::mutate(
			powerdb = 10*(
				log10(power) - log10( mean(power[ (time >= -2) & (time < -1) ]) ) # in seconds
			)
		)
		%>% dplyr::ungroup()
	)

	# # plot decibel normalized power — look OK?
	# (
	# 	dat_wt
	# 	%>% dplyr::filter(
	# 		chan == "C3"
	# 		, trial == 2
	# 	)
	# 	%>% ggplot()
	# 	+ geom_tile(
	# 		aes(
	# 			x = time
	# 			, y = as.factor(freq) # yuk
	# 			, fill = powerdb
	# 		)
	# 	)
	# )

	# remove extra second on each end since wavelet transform is now complete
	dat_wt <- (
		dat_wt
		%>% dplyr::filter(
			time >= -2
			, time <= 9.5
		)
	)

	# we can downsample again — both time and frequency bins
	print(paste("... downsampling round 2"))

	# downsample time round 2
	downsample_by_factor = 10
	dat_wt <- (
		dat_wt
		%>% dplyr::group_by(trial, chan, freq)
		# ensure sorted by time
		%>% dplyr::arrange(time)
		# add the sample num
		%>% dplyr::mutate(sample_num = 1:n())
		# keep only Nth sample
		# ==1 ensures that sample_nm==1 is kept, which is what you want to do
		# when the sample count starts at 1. If there is a sample_num==0, you'd
		# want to have ==0
		%>% dplyr::filter(
			(sample_num %% downsample_by_factor) == 1
		)
	)

	# downsampling frequency bins
	downsample_by_factor = 2
	dat_wt <- (
		dat_wt
		%>% dplyr::ungroup()
		%>% dplyr::group_by(trial, chan, time)
		# ensure sorted by period
		%>% dplyr::arrange(desc(freq))
		# add the sample num
		%>% dplyr::mutate(sample_num = 1:n())
		# keep only Nth sample
		# ==1 ensures that sample_nm==1 is kept, which is what you want to do
		# when the sample count starts at 1. If there is a sample_num==0, you'd
		# want to have ==0
		%>% dplyr::filter(
			(sample_num %% downsample_by_factor) == 1
		)
	)

	# # check some things:
	# length(unique(dat_wt$time))/(max(unique(dat_wt$time))-min(dat_wt$time)) # should be 10 Hz
	# range(dat_wt$time) # it seems downsampling is dropping edge times but that's OK
	# length(unique(dat_wt$time)) # this is the number of time samples (not sampling rate)
	# length(unique(dat_wt$freq)) # this is the number of frequency bins

	# get behavioural data for current eeg participant
	bdata <- bdat[bdat$participant == subject,]

	# create trial column (1-120)
	bdata <- mutate(bdata, trial = bdata$trial_num + (bdata$block_num-1)*20)

	# take only columns needed for model to reduce size
	bdata <- subset(bdata, select = c(participant, group, condition, trial, rep, complexity, avg_velocity, error, vresp, accuracy_rating)) # pick only what columns we really want

	# combine behavioural and eeg data by trial
	print(paste("... merging dat_wt and bdata"))
	dat <- merge(dat_wt, bdata, by.x = "trial", by.y = "trial", all.x = TRUE)

	rm(bdata, dat_wt)
	gc()

	# read eeg sensor coordinates (x,y,z)
	eegcoords <- read_csv("BESA-81.csv")

	# assume sphere and derive lat and long
	eegcoords$lat = 90 - (asin(eegcoords$z) * (180/pi))
	eegcoords$long = atan2(eegcoords$y,eegcoords$x) * (180/pi)
	eegcoords <- eegcoords[ eegcoords$chan %in% unique(dat$chan), ]

	# add coords to data frame
	print(paste("... merging dat and eegcoords"))
	dat <- merge(dat, eegcoords, by.x = "chan", by.y = "chan", all.x = TRUE)

	# clean up / reduce size:
	dat <- dplyr::filter(
		.data = dat
		# , is.na(eeg) == FALSE # there may be errors in eeg data
		, is.na(power) == FALSE
		, is.na(accuracy_rating) == FALSE
		, is.na(complexity) == FALSE # shouldn't happen, but just in case
		, is.na(time) == FALSE # shouldn't happen, but just in case
		, is.na(freq) == FALSE # shouldn't happen, but just in case
		, rep == 0 # random shape trials only = HUGE reduction in file size
	)

	# save data so far
	print(paste("... saving..."))
	saveRDS(dat,file=paste0('_Scipts/_rds/participants/',subject,'_participant_data.rds'))

	# add this participants data to all data:
	# invisible(ifelse(exists("alldat"), alldat <- rbind(alldat, dat), alldat <- dat))

	# free up memory
	rm(eegcoords, dat)
	gc()

	# # save data so far
	# print(paste("... saving..."))
	# fil <- "alldat.rds"
	# saveRDS(alldat, fil)

	print(paste("... participant", subject, "DONE"))
}

print(paste("All participants done!"))
Rtime <- proc.time() - ptm
print(Rtime)

#### combine all participants ####

print(paste("Combining all participants..."))
ptm <- proc.time()

library(tidyverse)

options(warn=1)
(
	list.files(
		path = 'rds_interpolated'
		, full.names = T
		, pattern = 'participant_data.rds'
	)
	%>% purrr::map_dfr(
		.f = function(x){
			dat = NULL
			try(dat<-readRDS(x))
			if(is.null(dat)){
				warning('Failed file: ',x)
			}
			return(dat)
		}
	)
) -> alldat

print(paste("... saving all participants in one object..."))
saveRDS(alldat, "alldat.rds")

print(paste("... ALL DONE!"))

Rtime <- proc.time() - ptm
print(Rtime)
