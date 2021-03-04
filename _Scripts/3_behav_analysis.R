	## Detecting Errors in Motor Imagery (DEMI) ##
	 ## descriptives and behavioural analyses ##

library(tidyverse)

# load data
dat <- readRDS('_Scripts/_rds/bdat.rds')
bad_imagery_trials <- readRDS('_Scripts/_rds/diagnostics/bad_imagery_trials.rds')
# note that bad physical trials have error metrics and mt_clip as NA

#### SETUP DATA ####

# remove pilot participants and bad trials
dat <- (
	dat
	%>% dplyr::filter(
		participant < 500 # remove pilot participants
		, !(condition == 'physical' & is.na(mt_clip)) # remove bad physical trials
		, !(figure_file %in% as.list(bad_imagery_trials$figure_file)) # remove bad imagery trials
	)
)
rm(bad_imagery_trials)

# make participant a factor
dat$participant = factor(dat$participant)

# clarify group vs condition
dat$group <- factor(dat$exp_condition)
dat$condition <- factor(dat$condition)

# create potentially useful dummy variables
# note: recall difference between 'group' and 'condition'
dat$grp <- ifelse(dat$group=="imagery",1,0) # physical group = 0 ; imagery group = 1
dat$MI <- ifelse(dat$condition=="imagery", 1, 0) # is condition imagery = 1
dat$PP <- ifelse(dat$condition=="physical", 1, 0) # is condition physical = 0
dat$rep <- ifelse(dat$figure_type=="repeated", 1, 0) # random = 0 ; repeated = 1

# choose measure of complexity
dat$complexity <- dat$sinuosity
# choose measure of physical error
# dat$error <- dat$raw_err_tot

# update mt with mt_clip where applicable
dat$mt <- ifelse(!is.na(dat$mt_clip),dat$mt_clip,dat$mt)

# save for later analyses
saveRDS(dat, '_Scripts/_rds/bdat2.rds')

#### Participant Characterization ####

# age:

print("mean (sd) age for all participants:")
mean(dat[!duplicated(dat$participant),]$age)
sd(dat[!duplicated(dat$participant),]$age)

print("mean (sd) age for each group:")
dat[!duplicated(dat$participant),] %>%
	group_by(group) %>%
	summarise( count = length(unique(participant))
			   , agemean = mean(age)
			   , agesd = sd(age)
			   )

# bio sex

# TO DO

# handedness

# TO DO

#### Mental Chronometry ####

# visualize matching of MT for each condition:

ggplot(data = dat
	   , mapping = aes(
	   	x = stimulus_mt
	   	, y = mt
	   	, color = factor(block_num)
	   )) + geom_point(na.rm = TRUE, alpha = .5) +
	facet_grid(. ~ condition) +
	geom_smooth(na.rm = TRUE) +
	# geom_density2d(na.rm = TRUE) +
	theme_minimal() #+
	# lims(x = c(0, 3), y = c(0, 7)) # note outliers

# to do

#### Error and Accuracy Analyses ####

# visualize error

ggplot(data = subset(dat, (group == 'physical'))
	   , mapping = aes(
	   	x = vresp #
	   	, y = (shape_dtw_err_mean*0.2715) # convert to mm
	   )) + geom_point(na.rm = TRUE, alpha = .5) +
	# facet_grid(. ~ figure_type) +
	facet_wrap(~participant) +
	geom_smooth(method = lm, na.rm = TRUE) +
	# geom_density2d(na.rm = TRUE) +
	theme_minimal() +
	theme(aspect.ratio=4/4)


# to do


