##########################
### Plot raw EEG data  ###
##########################


library(tidyverse)
options(clustermq.scheduler = "multicore")

scale_to_0range = function(x,range=1){
	x = x - min(x)
	x = x/max(x)
	x = x-.5
	x = x * range
	return(x)
}

(
	'_Scripts/_rds/all_dat.rds'
	# '_Scripts/_rds/all_dat_first_participant.rds'
	%>% readRDS()
	%>% as_tibble()
	%>% mutate(
		time = time + (2 * (epoch - 1))
	)
) -> dat

#add handedness & flip lefties
(
	readRDS('_Scripts/_rds/bdat2.rds')
	%>% mutate(participant=as.character(participant))
	%>% group_by(
		participant
	)
	%>% summarize(
		handedness = handedness[1]
		, .groups = 'drop'
	)
	%>% right_join(
		dat	%>% mutate(participant=as.character(participant))
		, by = 'participant'
	)
	%>% mutate(

		#polar to cartesian:
		x = lat*cos(long*(pi/180))
		, y = lat*sin(long*(pi/180))

		#flip lefties
		, x = case_when(
			handedness=='l' ~ -x
			, T ~ x
		)
	)
) -> dat


#for each participant, compute mean across trials
(
	dat
	%>% group_by(group,participant,x,y,time,freq)
	%>% summarize(
		powerdb = mean(powerdb)
		, .groups = 'drop'
	)
) -> dat

#for each group, compute mean across participants and append as 0th participant
(
	dat
	%>% group_by(group,x,y,time,freq)
	%>% summarize(
		powerdb = mean(powerdb)
		, .groups = 'drop'
	)
	%>% mutate(participant='0')
	%>% dplyr::bind_rows(
		(
			dat
			%>% select(group,participant,x,y,time,freq,powerdb)
		)
	)
	%>% mutate(participant=as.numeric(participant))
) -> dat



#determine freqs to NA for band boundaries
unique_freqs = unique(dat$freq)
band_boundaries = c(
	unique_freqs[which.min(abs(unique_freqs-4))]
	, unique_freqs[which.min(abs(unique_freqs-8))]
	, unique_freqs[which.min(abs(unique_freqs-16))]
	, unique_freqs[which.min(abs(unique_freqs-32))]
)
unique_times = unique(dat$time)
# zero_time = unique_times[which.min(abs(unique_times))]
time_boundaries = c(
	0, 1.5, 2
)

#add stuff
(
	dat
	%>% ungroup()
	%>% mutate(

		#NA powerdb on band boundaries
		powerdb = case_when(
			!(freq %in% band_boundaries) ~ powerdb
		)

		#NA powerdb at time==0
		, powerdb = case_when(
			# time!=zero_time ~ powerdb
			!(time %in% time_boundaries) ~ powerdb
		)

		#rescaling
		, x_scaled = scale_to_0range(x,10)
		, y_scaled = scale_to_0range(y,10)
		, time_scaled = scale_to_0range(time,1)
		, freq_scaled = scale_to_0range(freq,1)

		#sum to get final positions
		, to_plot_x = x_scaled + time_scaled
		, to_plot_y = y_scaled + freq_scaled

	)
) -> dat

get_raster_for_channel = function(this_channel_data,all_channels_powerdb_range){
	# (
	# 	dat
	# 	%>% filter(
	# 		participant==participant[1]
	# 		, chan==chan[1]
	# 		, trial==0
	# 	)
	# )->this_channel_data
	# all_channels_powerdb_range = range(this_channel_data$powerdb)
	(
		this_channel_data
		%>% ggplot()
		+ geom_raster(
			mapping = aes(
				x = to_plot_x
				, y = to_plot_y
				, fill = powerdb
			)
		)
		+ scale_fill_viridis_c(
			option='C'
			, limits = all_channels_powerdb_range
			, na.value = 'transparent'
		)
		+ theme(
			legend.position = "top"
			, legend.title = element_blank()
		)
	) -> temp
	tmp = ggplot_gtable(ggplot_build(temp))
	return(list(
		legend = tmp$grobs[[which(sapply(tmp$grobs, function(x) x$name) == "guide-box") ]]
		, topo_raster = layer_grob(temp)$`1`$raster
		, xmin = min(this_channel_data$to_plot_x)
		, xmax = max(this_channel_data$to_plot_x)
		, ymin = min(this_channel_data$to_plot_y)
		, ymax = max(this_channel_data$to_plot_y)
	))
}

plot_participant = function(this_participant_data,all_participants_powerdb_range){
	# (
	# 	dat
	# 	%>% filter(
	# 		participant==participant[1]
	# 		, participant==0
	# 		# , chan==chan[1]
	# 	)
	# ) -> this_participant_data

	#determine which range to use
	if(this_participant_data$participant[1]==0){
		this_participant_powerdb_range = range(this_participant_data$powerdb,na.rm=T)
	}else{
		this_participant_powerdb_range = all_participants_powerdb_range
	}
	#get the rasters per channel
	(
		this_participant_data
		%>% group_by(x,y)
		%>% group_split()
		%>% purrr::map(
			.f = get_raster_for_channel
			, all_channels_powerdb_range = this_participant_powerdb_range
		)
	) -> channel_rasters

	#construct the canvas with limits & theme stuff
	(
		ggplot()
		+ scale_x_continuous(limits = range(this_participant_data$to_plot_x))
		+ scale_y_continuous(limits = range(this_participant_data$to_plot_y))
		+ coord_equal()
		+ labs(
			title = this_participant_data$participant[1]
		)
		+ theme(
			axis.title = element_blank()
			, axis.ticks = element_blank()
			, axis.text = element_blank()
			, panel.grid = element_blank()
			, panel.background = element_rect(fill='transparent',colour='grey90')
		)
	) -> p

	#add each channel's raster
	for(chan in channel_rasters){
		(
			p
			+ annotation_raster(
				raster = chan$topo_raster
				, xmin = chan$xmin
				, xmax = chan$xmax
				, ymin = chan$ymin
				, ymax = chan$ymax
			)
		)->p
	}

	#add legend
	(
		p
		+ annotation_custom(
			chan$legend
			, xmin = max(this_participant_data$to_plot_x)-.2*diff(range(this_participant_data$to_plot_x))
			, xmax = max(this_participant_data$to_plot_x)
			, ymin = min(this_participant_data$to_plot_y)-.1*diff(range(this_participant_data$to_plot_y))
			, ymax = min(this_participant_data$to_plot_y)+.2*diff(range(this_participant_data$to_plot_y))
		)
	) -> p

	#plot
	print(p)
	return(NULL)
}

plot_group_participants = function(this_group_data,path){
	# (
	# 	dat
	# 	%>% filter(
	# 		group==group[1]
	# 		# , trial==0
	# 		# , chan==chan[1]
	# 	)
	# ) -> this_group_data

	#prepare the directory structure
	path_list = c('plots','raw_power',this_group_data$group[1])
	for(i in 1:length(path_list)){
		path = paste(path_list[1:i],collapse='/')
		if(!dir.exists(path)){
			dir.create(path)
		}
	}

	pdf(
		file = paste0(path,'/group.pdf')
		, height = 8
		, width = 8
	)
	(
		this_group_data
		%>% group_by(participant)
		%>% group_split()
		%>% purrr::map(
			.f = plot_participant
			, all_participants_powerdb_range = range(this_group_data$powerdb,na.rm=T)
		)
	)
	dev.off()
	return(NULL)
}


(
	dat
	%>% group_by(group)
	%>% group_split()
	%>% purrr::map(
		.f = plot_group_participants
	)
)

