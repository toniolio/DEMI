	## Detecting Errors in Motor Imagery (DEMI) ##
	 ## descriptives and behavioural analyses ##

library(tidyverse)
# library(lme4)
# library(lmerTest)
library(brms)

# TO DO: ifexists bdat2.rds stuff

# load data
dat <- readRDS('_Scripts/_rds/bdat.rds')
bad_imagery_trials <- readRDS('_Scripts/_rds/bad_imagery_trials.rds')
# note that bad physical trials have error metrics and mt_clip as NA

r2.corr.mer <- function(m) {
	lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
	summary(lmfit)$r.squared
} # a very crude R2 for linear mixed models
# from: http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html

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
dat$rep <- ifelse(dat$figure_type=="random", 0, 1) # random = 0 ; repeated = 1

# choose measure of complexity
dat$complexity <- dat$sinuosity # note publication used totabscurv
# choose measure of physical error
dat$error <- dat$dtw_err_mean

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

dat[!duplicated(dat$participant),] %>%
	group_by(group) %>%
	count(sex)

# handedness

dat[!duplicated(dat$participant),] %>%
	group_by(group) %>%
	count(handedness)

#### Mental Chronometry ####

# visualize matching of MT for each group and condition:

ggplot(data = dat
	   , mapping = aes(
	   	x = stimulus_mt
	   	, y = mt
	   	, color = factor(block_num)
	   )) + geom_point(na.rm = TRUE, alpha = .5) +
	facet_grid(group ~ condition) +
	geom_smooth(na.rm = TRUE) +
	# geom_density2d(na.rm = TRUE) +
	theme_minimal() #+
	# lims(x = c(0, 3), y = c(0, 7)) # note outliers

# looks great

	#------------------------#
	## statistical analyses ##
	#------------------------#

# NOTE: if we do bayes later, do bayes here.

#----------------#
# physical group #
#----------------#

# MC.lmm.1 <-(
# 	dat
# 	%>% dplyr::filter(
# 		group == 'physical'
# 	)
# 	%>% lmer(
# 		mt ~ stimulus_mt + (1 + stimulus_mt | participant)
# 		, data = .
# 	)
# )
# summary(MC.lmm.1)
# r2.corr.mer(MC.lmm.1)

MC.blmm.1 <- (
	dat
	%>% dplyr::filter(
		group == 'physical'
	)
	%>% brms::brm(
		formula = (
			mt ~ stimulus_mt + (1 + stimulus_mt | participant)
			)
		, data = .
		, silent = F
		, refresh = 20
		, chains = 4
		, iter = 2000
		, cores = 4
	)
)
summary(MC.blmm.1)

# yes, they scale their response to speed


#-----------------------------------#
# imagery group — imagery condition #
#-----------------------------------#

# MC.lmm.2 <-(
# 	dat
# 	%>% dplyr::filter(
# 		group == 'imagery'
# 		, condition == 'imagery'
# 	)
# 	%>% lmer(
# 		mt ~ stimulus_mt + (1 + stimulus_mt | participant)
# 		, data = .
# 	)
# )
# summary(MC.lmm.2)
# r2.corr.mer(MC.lmm.2)

MC.blmm.2 <-(
	dat
	%>% dplyr::filter(
		group == 'imagery'
		, condition == 'imagery'
	)
	%>% brms::brm(
		formula = (
			mt ~ stimulus_mt + (1 + stimulus_mt | participant)
		)
		, data = .
		, silent = F
		, refresh = 20
		, chains = 4
		, iter = 2000
		, cores = 4
	)
)
summary(MC.blmm.2)

# yes, the imagery group scales their speed


#------------------------------------#
# imagery group — physical condition #
#------------------------------------#

# MC.lmm.3 <-(
# 	dat
# 	%>% dplyr::filter(
# 		group == 'imagery'
# 		, condition == 'physical'
# 	)
# 	%>% lmer(
# 		mt ~ stimulus_mt + (1 + stimulus_mt | participant)
# 		, data = .
# 	)
# )
# summary(MC.lmm.3)
# r2.corr.mer(MC.lmm.3)

MC.blmm.3 <-(
	dat
	%>% dplyr::filter(
		group == 'imagery'
		, condition == 'physical'
	)
	%>% brms::brm(
		formula = (
			mt ~ stimulus_mt + (1 + stimulus_mt | participant)
		)
		, data = .
		, silent = F
		, refresh = 20
		, chains = 4
		, iter = 2000
		, cores = 4
	)
)
summary(MC.blmm.3)

#-------------------------------------#
# imagery group — physical vs imagery #
#-------------------------------------#

# MC.lmm.4 <-(
# 	dat
# 	%>% dplyr::filter(
# 		group == 'imagery'
# 	)
# 	%>% lmer(
# 		mt ~ stimulus_mt * condition + (1 + stimulus_mt * condition | participant)
# 		, data = .
# 		, REML = TRUE
# 	)
# )
# summary(MC.lmm.4)
# r2.corr.mer(MC.lmm.4)

MC.blmm.4 <-(
	dat
	%>% dplyr::filter(
		group == 'imagery'
	)
	%>% brms::brm(
		formula = (
			mt ~ stimulus_mt * condition + (1 + stimulus_mt * condition | participant)
		)
		, data = .
		, silent = F
		, refresh = 20
		, chains = 4
		, iter = 2000
		, cores = 4
	)
)
summary(MC.blmm.4)

# no meaningful difference between conditions, and while
# stimulus * condition interaction is "significant" the estimate is
# so small that it's not meaningful either.

# it seems imagery participants perform their physical trials at
# roughly the same speeds.

#--------------------------------#
# compare groups — physical only #
#--------------------------------#

# MC.lmm.5 <-(
# 	dat
# 	%>% dplyr::filter(
# 		condition == 'physical'
# 	)
# 	%>% lmer(
# 		mt ~ stimulus_mt * group + (1 + stimulus_mt * group | participant)
# 		, data = .
# 		, REML = TRUE
# 	)
# )
# summary(MC.lmm.5)
# r2.corr.mer(MC.lmm.5)

MC.blmm.5 <-(
	dat
	%>% dplyr::filter(
		condition == 'physical'
	)
	%>% brms::brm(
		formula = (
			mt ~ stimulus_mt * group + (1 + stimulus_mt * group | participant)
		)
		, data = .
		, silent = F
		, refresh = 20
		, chains = 4
		, iter = 2000
		, cores = 4
	)
)
summary(MC.blmm.5)

# when performing physical trials, physical group appears to be faster
# (lower intercept) than imagery group, and a less steep slope but that
# "significant" effect is tiny and not very meaningful

# it may be that the imagery task causes participants to approach
# physical trials slower, as some sort of after-effect.

#-------------------------------------------------------#
# compare groups — by conditions (imagery vs. physical) #
#-------------------------------------------------------#

# that is, compare groups but remove physical trials from
# the imagery group.

# MC.lmm.6 <-(
# 	dat
# 	%>% dplyr::filter(
# 		!(group == 'imagery' & condition == 'physical')
# 	)
# 	%>% lmer(
# 		mt ~ stimulus_mt * group + (1 + stimulus_mt * group | participant)
# 		, data = .
# 		, REML = TRUE
# 	)
# )
# summary(MC.lmm.6)
# r2.corr.mer(MC.lmm.6)

MC.blmm.6 <-(
	dat
	%>% dplyr::filter(
		!(group == 'imagery' & condition == 'physical')
	)
	%>% brms::brm(
		formula = (
			mt ~ stimulus_mt * group + (1 + stimulus_mt * group | participant)
		)
		, data = .
		, silent = F
		, refresh = 20
		, chains = 4
		, iter = 2000
		, cores = 4
	)
)
summary(MC.blmm.6)

# the physical group (and condition) seems to move faster than the imagery
# group (and condition), but it's not a huge effect

# overall, the expectation holds that imagery participants match the movement
# times of the stimulus well.


#### Error and Accuracy Analyses ####

#----------------------------#
# error ~ speed * complexity #
#----------------------------#

# visualize error as it is related to speed:

ggplot(data = subset(dat, (group == 'physical'))
	   , mapping = aes(
	   	x = vresp
	   	, y = (dtw_err_mean*0.2715) # convert to mm
	   )) + geom_point(na.rm = TRUE, alpha = .5) +
	# facet_grid(. ~ figure_type) +
	facet_wrap(~participant) +
	geom_smooth(method = lm, na.rm = TRUE) +
	# geom_density2d(na.rm = TRUE) +
	theme_minimal()
	# theme(aspect.ratio=4/4)

# visualize error as it is related to complexity:

ggplot(data = subset(dat, (group == 'physical' & rep == 0))
	   , mapping = aes(
	   	x = sinuosity
	   	, y = (dtw_err_mean*0.2715) # convert to mm
	   )) + geom_point(na.rm = TRUE, alpha = .5) +
	# facet_grid(. ~ figure_type) +
	facet_wrap(~participant) +
	geom_smooth(method = lm, na.rm = TRUE) +
	# geom_density2d(na.rm = TRUE) +
	theme_minimal()
	# theme(aspect.ratio=4/4)

# model relationship: error ~ speed * complexity
# to confirm that error measure is sensible.

# AA.lmm.1 <- (
# 	dat
# 	%>% dplyr::filter(
# 		condition == 'physical'
# 		, rep == 0
# 	)
# 	%>% lmer(
# 		scale(dtw_err_mean) ~ scale(sinuosity) * scale(vresp) + (1 + scale(sinuosity) * scale(vresp) | participant)
# 		, data = .
# 		# , REML = TRUE
# 		# , control=glmerControl(optimizer="NM2")
# 	)
# )
# summary(AA.lmm.1)
# r2.corr.mer(AA.lmm.1)

AA.blmm.1 <- (
	dat
	%>% dplyr::filter(
		condition == 'physical'
		, rep == 0
	)
	%>% brms::brm(
		formula = (
			scale(dtw_err_mean) ~ scale(sinuosity) * scale(vresp) + (1 + scale(sinuosity) * scale(vresp) | participant)
		)
		, data = .
		, silent = F
		, refresh = 20
		, chains = 4
		, iter = 2000
		, cores = 4
	)
)
summary(AA.blmm.1)

#-----------------------------------------#
# perceived accuracy ~ speed * complexity #
#-----------------------------------------#

# for both groups, but leaving out repeated trials (as they don't vary complexity much)

# visualize perceived accuracy predicted by speed:

ggplot(data = subset(dat, (rep == 0))
	   , mapping = aes(
	   	x = vresp # convert to mm
	   	, y = scale(accuracy_rating) # NOTE: consider scaling!
	   )) + geom_point(na.rm = TRUE, alpha = .5) +
	# facet_grid(. ~ figure_type) +
	facet_wrap(~participant) +
	geom_smooth(method = lm, na.rm = TRUE) +
	# geom_density2d(na.rm = TRUE) +
	theme_minimal() #+
	# lims(y = c(0, 11)) # note outliers
	# theme(aspect.ratio=4/4)

# visualize perceived accuracy predicted by complexity:

ggplot(data = subset(dat, (rep == 0))
	   , mapping = aes(
	   	x = sinuosity # convert to mm
	   	, y = scale(accuracy_rating) # NOTE: consider scaling!
	   )) + geom_point(na.rm = TRUE, alpha = .5) +
	# facet_grid(. ~ figure_type) +
	facet_wrap(~participant) +
	geom_smooth(method = lm, na.rm = TRUE) +
	# geom_density2d(na.rm = TRUE) +
	theme_minimal() #+
	# lims(y = c(0, 11)) # note outliers
	# theme(aspect.ratio=4/4)

# model perceived accuracy ~ speed * complexity

# AA.lmm.2 <- (
# 	dat
# 	%>% dplyr::filter(
# 		rep == 0
# 	)
# 	%>% lmer(
# 		scale(accuracy_rating) ~ scale(sinuosity) * scale(vresp) + (1 + scale(sinuosity) * scale(vresp) | participant)
# 		, data = .
# 		# , REML = TRUE
# 		# , control=glmerControl(optimizer="NM2")
# 	)
# )
# summary(AA.lmm.2)
# r2.corr.mer(AA.lmm.2)

AA.blmm.2 <- (
	dat
	%>% dplyr::filter(
		rep == 0
	)
	%>% brms::brm(
		formula = (
			scale(accuracy_rating) ~ scale(sinuosity) * scale(vresp) + (1 + scale(sinuosity) * scale(vresp) | participant)
		)
		, data = .
		, silent = F
		, refresh = 20
		, chains = 4
		, iter = 2000
		, cores = 4
	)
)
summary(AA.blmm.2)

# like actual error, perceived accuracy is predicted by both complexity and
# speed but there is no interaction between complexity and speed. This is good
# evidence that perceived accuracy is capturing a similar phenomenon as error,
# one that is dependent on known drivers of task performance.

# according to the R2 estimate... accuracy_rating is actually better predicted
# by speed and complexity than actual error? this estimate of R2 must leave out
# random effects / correlations because participants obviously vary more here
# than observed in the actual error data.

#-------------------------------------------------#
# perceived accuracy ~ error (physical condition) #
#-------------------------------------------------#

# visualize perceived accuracy predicted by actual error:

ggplot(data = subset(dat, (condition == 'physical'))
	   , mapping = aes(
	   	x = (dtw_err_mean*0.2715) # convert to mm
	   	, y = scale(accuracy_rating) # NOTE: consider scaling!
	   )) + geom_point(na.rm = TRUE, alpha = .5) +
	# facet_grid(. ~ figure_type) +
	facet_wrap(~participant) +
	geom_smooth(method = lm, na.rm = TRUE) +
	# geom_density2d(na.rm = TRUE) +
	theme_minimal() #+
	# lims(y = c(0, 11)) # note outliers
	# theme(aspect.ratio=4/4)

# model perceived accuracy ~ error

# AA.lmm.3 <- (
# 	dat
# 	%>% dplyr::filter(
# 		condition == 'physical'
# 	)
# 	%>% lmer(
# 		scale(accuracy_rating) ~ scale(dtw_err_mean) + (1 + scale(dtw_err_mean) | participant)
# 		, data = .
# 		# , REML = TRUE
# 		# , control=glmerControl(optimizer="NM2")
# 	)
# )
# summary(AA.lmm.3)
# r2.corr.mer(AA.lmm.3)

AA.blmm.3 <- (
	dat
	%>% dplyr::filter(
		condition == 'physical'
	)
	%>% brms::brm(
		formula = (
			scale(accuracy_rating) ~ scale(dtw_err_mean) + (1 + scale(dtw_err_mean) | participant)
		)
		, data = .
		, silent = F
		, refresh = 20
		, chains = 4
		, iter = 2000
		, cores = 4
	)
)
summary(AA.blmm.3)

# R2 of .55 is quite good validation of this likert scale, I think.


#-------------------------------------------------------------#
# perceived accuracy ~ speed * complexity (imagery condition) #
#-------------------------------------------------------------#

# note that must use stimulus velocity, and not participant actual velocity
# try again with stimulus_mt?

AA.blmm.4 <- (
	dat
	%>% dplyr::filter(
		condition == 'imagery'
		# , rep == 0 # this doesn't change results much, but note complexity in rep == 1 doesn't vary much
	)
	%>% brms::brm(
		formula = (
			scale(accuracy_rating) ~ scale(sinuosity) * scale(avg_velocity) + (1 + scale(sinuosity) * scale(avg_velocity) | participant)
		)
		, data = .
		, silent = F
		, refresh = 20
		, chains = 4
		, iter = 2000
		, cores = 4
	)
)
summary(AA.blmm.4)

# Results are as expected — not different from above models (perceived accuracy ~ speed * complexity)
# showing that, just like actual error, perceived accuracy is predicted by
# both speed and complexity, but speed and complexity do not interact.
