#############################################
### Descriptives and behavioural analyses ###
#############################################

library(tidyverse)
# library(lme4)
# library(lmerTest)
library(brms)
library(bayestestR)

# not used but potentially useful:
# r2.corr.mer <- function(m) {
# 	lmfit <-  lm(model.response(model.frame(m)) ~ fitted(m))
# 	summary(lmfit)$r.squared
# } # a very crude R2 for linear mixed models
# # from: http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html

path <- "_Scripts/_rds/" # for saving models

#### SETUP DATA ####

if(!file.exists(paste0("_Scripts/_rds/bdat2.rds"))){

	# load data
	dat <- readRDS('_Scripts/_rds/bdat.rds')
	bad_imagery_trials <- readRDS('_Scripts/_rds/bad_imagery_trials.rds')
	# note that bad physical trials have error metrics and mt_clip as NA

	# remove pilot participants
	dat <- (
		dat
		%>% dplyr::filter(
			participant < 500
			)
	)

	# make participant a factor
	dat$participant = factor(dat$participant)

	# clarify group vs condition
	dat$group <- factor(dat$exp_condition)
	dat$condition <- factor(dat$condition)

	# back fill the vividness ratings for each block
	dat <- (
		dat
		%>% group_by(participant, block_num)
		%>% mutate(
			vivid = ifelse(condition == 'imagery',vividness_rating[!is.na(vividness_rating)],NA)
		)
		%>% ungroup()
	)

	# remove bad trials
	dat <- (
		dat
		%>% dplyr::filter(
			!(condition == 'physical' & is.na(mt_clip)) # remove bad physical trials
			, !(figure_file %in% as.list(bad_imagery_trials$figure_file)) # remove bad imagery trials
		)
	)
	rm(bad_imagery_trials)

	# make rep a factor:
	dat$rep <- ifelse(dat$figure_type=="random", 0, 1) # random = 0 ; repeated = 1
	dat$rep <- as.factor(ifelse(dat$rep==0, "random", "repeated")) # better label

	# choose measure of complexity
	dat$complexity <- dat$sinuosity # note prev publication used totabscurv

	# choose measure of physical error
	dat$error <- dat$dtw_err_mean

	# update mt with mt_clip where applicable
	dat$mt <- ifelse(!is.na(dat$mt_clip),dat$mt_clip,dat$mt)

	# save for later analyses
	saveRDS(dat, '_Scripts/_rds/bdat2.rds')

} else {
	dat <- readRDS("_Scripts/_rds/bdat2.rds")
}

# for these stats specifically
halfsum_contrasts = function (...) contr.sum(...) * 0.5
contrasts(dat$group) = halfsum_contrasts
contrasts(dat$condition) = halfsum_contrasts
contrasts(dat$rep) = halfsum_contrasts
dat$trial = (dat$trial_num + (dat$block_num - 1) * 20)
dat$blockz = dat$block_num-3


#### Participant Characterization ####

# note this is for behavioural analysis only
# for the participants in the EEG paper, see later scrips: e.g., "run_gams.R"

unique(sort(dat$participant))
# lost due to experimenter error: 9, 10, 12, 14, 20,
# didn't complete experiment (tech issue): 89, 96, 100
# lost due to very bad EEG: 24
# id skipped: 13, 26, 36, 78

# therefore, while highest participant number is 99, in fact we recruited p100
# but skipped 4 ids = 96 were recruited, and 5 dropped because of experimenter
# issue (used wrong experimental protocol), and then 3 (tech) + 1 (bad) = 4
# were dropped due to a technical issue with our equipment. The "bad EEG"
# participant completed the behavioral experiment but the EEG was so bad that
# it we didn't trust that the experiment was completed properly so decided to
# take a conservative approach and drop them completely.

# in the paper, we will describe 96 recruited, and 9 dropped due to technical
# issues with the experimental setup.

# which participants in each group (check)

unique(sort(subset(dat, group =='physical')$participant))
unique(sort(subset(dat, group =='imagery')$participant))

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

#### Hypothesis 1 ####

# If motor imagery requires an existing representation to predict error, then
# perceived accuracy during a novel motor task will have a poor relationship
# with expected error (given the speed and complexity), compared to the PP
# groups actual relationship between their perceived accuracy and actual error.

# build a model using all physical trials (including from the imagery group)
# that explains actual error from speed, complexity, etc.

# define "skill"
dat <- (
	dat
	%>% dplyr::mutate(
		skill = vresp/error # flipped to make performance improvements positive
	) # scale everything manually
	%>% dplyr::mutate(
		skillz = scale(skill)[,1]
		, trialz = scale(trial)[,1]
		, complexityz = scale(complexity)[,1]
		, self_accz = scale(accuracy_rating)[,1]
		, vividz = scale(vivid)[,1]
	)
)

# check some assumptions:
# does vividz and self_accz vary similarly?

# plot distribution:
plot(density(dat$vividz, na.rm =  TRUE)
	 , main = "Vividness"
	 , xlab = "Response (z-score)"
	 )
plot(density(dat$self_accz, na.rm =  TRUE)
	 , main = "Self-report Accuracy"
	 , xlab = "Response (z-score)"
	 )
# calculate coefficient of variation:
(sd(dat$vividness_rating, na.rm =  TRUE)/abs(mean(dat$vividness_rating, na.rm =  TRUE)))*100 # as percentage
(sd(dat$accuracy_rating, na.rm =  TRUE)/abs(mean(dat$accuracy_rating, na.rm =  TRUE)))*100 # as percentage

H1 <- (
	dat
	%>% dplyr::filter(
		condition == 'physical' # include both groups (only 6th block for MI)
	) # note it's important to include imagery group physical trials for prediction step below
	%>% brms::brm(
		formula = (
			skillz ~ trialz * rep + complexityz
			+ (1 + trialz * rep + complexityz | participant)
		)
		, data = .
		, prior = c(prior(normal(0,10), class = b))
		, silent = F
		, refresh = 20
		, chains = floor(parallel::detectCores())
		, iter = 10000
		, cores = floor(parallel::detectCores())
	)
)
saveRDS(H1, paste0(path, "H1.rds"))
# H1 <- readRDS(paste0(path, "H1.rds"))
summary(H1)
bayes_R2(H1)
rope(H1
	 , range = c(-0.1,0.1) # assume small effect size negligible
	 , ci = 1 # "full ROPE"
	 , ci_method = "HDI"
	 )

plot1brm <- plot(conditional_effects(H1
								  , effects = "trialz:rep"
								  , re_formula = NA
								  , prob = 0.95
								  )
			  )[[1]]


# create nice plot for publication
plot1df <- as.data.frame(plot1brm$data)
levels(plot1df$effect2__) <- c("Random","Repeated")
(
	ggplot(
		plot1df
		, aes(
			x = (effect1__ * attr(scale(dat$trial), "scaled:scale"))
			+ attr(scale(dat$trial), "scaled:center")
		   	, y = estimate__
		   	, color = effect2__
		   )
	)
	+ geom_line(
		aes(
			color = effect2__
			, linetype = effect2__
		)
		, size = .5
	)
	+ geom_ribbon(
		aes(
			ymin=lower__
			, ymax=upper__
			, fill = effect2__
			, alpha = effect2__
			)
		, color = NA
		)
	+ scale_linetype_manual(
		values = c(
			"Random" = "dashed"
			, "Repeated" = "solid"
		)
	)
	+ scale_alpha_manual(
		values = c(
			"Random" = .2
			, "Repeated" = .4
		)
	)
	+ labs(
		x = "Trial"
		, y = "Actual Performance (z-score)"
		, color  = "Trial Type"
		, linetype = "Trial Type"
		, fill = "Trial Type"
		, alpha = "Trial Type"
		)
	+ scale_x_continuous(
		breaks=c(1,20,40,60,80,100,120)
		)
	+ theme_minimal()
	+ theme(
		text=element_text(
			family="Times"
			, size = 7
		)
		# , legend.position="none"
	)
	+ theme(axis.title.x = element_text(margin = margin(t = 5)))
	+ theme(legend.title=element_blank())
)
ggsave(
	file = paste0('_Scripts/_plots/Fig2_90mm.jpeg')
	, plot = last_plot()
	, units = "mm"
	, width = 90
	, height = 60
)

# side bar (requested by reviewer): compare actual performance in final block
# between groups (imagery vs. overt)
HR1 <- (
	dat
	%>% dplyr::filter(
		condition == 'physical' # include both groups (only 6th block for MI)
		, block_num == 6 # also restrict to 6th block of physical group
	)
	%>% brms::brm(
		formula = (
			skillz ~ group * rep
			+ (1 + group * rep | participant)
		)
		, data = .
		, prior = c(prior(normal(0,10), class = b))
		, silent = F
		, refresh = 20
		, chains = floor(parallel::detectCores())
		, iter = 10000
		, cores = floor(parallel::detectCores())
	)
)
saveRDS(HR1, paste0(path, "HR1.rds"))
# HR1 <- readRDS(paste0(path, "HR1.rds"))
summary(HR1)
bayes_R2(HR1)

# use H1 model to derive “expected skill”

# note that for imagery trials, this will use "avg_velocity" instead of "vresp"
# which means rather than using the speed of actual execution (vresp), it uses
# the speed of the stimulus (note results from mental chronometry)
expected_skillz <- predict(H1
						  , newdata = dat # use all data
						  , re_formula = NULL # include participant effects
						  , summary = TRUE
						  , cores = floor(parallel::detectCores())
)

dat$expected_skillz <- expected_skillz[,1]
rm(expected_skillz)
gc()

# takes a while to run, so save a temporary version of data with expected_skillz
saveRDS(dat, paste0(path, "bdat_temp.rds"))
# dat <- readRDS(paste0(path, "bdat_temp.rds"))

# Finally, test how well perceived accuracy ratings are correlated with:
# a. Physical condition actual skill (this is more of a sanity check)
# b. Physical condition expected skill
# c. Imagery condition expected skill

H1a <- (
	dat
	%>% dplyr::filter(
		condition == 'physical' # include both groups (an extra 6th block)
	)
	%>% brms::brm(
		formula = (
			self_accz ~
				skillz # actual skill
			+ (1 + skillz | participant)
		)
		, data = .
		, prior = c(prior(normal(0,10), class = b))
		, silent = F
		, refresh = 20
		, chains = floor(parallel::detectCores())
		, iter = 5000
		, cores = floor(parallel::detectCores())
	)
)
saveRDS(H1a, paste0(path, "H1a.rds"))
# H1a <- readRDS(paste0(path, "H1a.rds"))
summary(H1a)
bayes_R2(H1a)
plot(conditional_effects(H1a
						 , re_formula = NA
						 )
	 )

H1b <- (
	dat
	%>% dplyr::filter(
		condition == 'physical' # include both groups (an extra 6th block)
	)
	%>% brms::brm(
		formula = (
			self_accz ~
				expected_skillz # expected skill
			+ (1 + expected_skillz | participant)
		)
		, data = .
		, prior = c(prior(normal(0,10), class = b))
		, silent = F
		, refresh = 20
		, chains = floor(parallel::detectCores())
		, iter = 5000
		, cores = floor(parallel::detectCores())
	)
)
saveRDS(H1b, paste0(path, "H1b.rds"))
# H1b <- readRDS(paste0(path, "H1b.rds"))
summary(H1b)
bayes_R2(H1b)
plot(conditional_effects(H1b
						 , re_formula = NA
						 )
	 )

H1c <- (
	dat
	%>% dplyr::filter(
		condition == 'imagery' # just the imagery group
	)
	%>% brms::brm(
		formula = (
			self_accz ~
				expected_skillz # expected skill
			+ (1 + expected_skillz | participant)
		)
		, data = .
		, prior = c(prior(normal(0,10), class = b))
		, silent = F
		, refresh = 20
		, chains = floor(parallel::detectCores())
		, iter = 5000
		, cores = floor(parallel::detectCores())
	)
)
saveRDS(H1c, paste0(path, "H1c.rds"))
# H1c <- readRDS(paste0(path, "H1c.rds"))
summary(H1c)
bayes_R2(H1c)
plot(conditional_effects(H1c
						 , re_formula = NA
						 )
	 )

H1c.v <- (
	dat
	%>% dplyr::filter(
		condition == 'imagery' # just the imagery group
	)
	%>% brms::brm(
		formula = (
			self_accz ~
				vividz * expected_skillz # expected skill
			+ (1 + vividz * expected_skillz | participant)
		)
		, data = .
		, prior = c(prior(normal(0,10), class = b))
		, silent = F
		, refresh = 20
		, chains = floor(parallel::detectCores())
		, iter = 5000
		, cores = floor(parallel::detectCores())
	)
)
saveRDS(H1c.v, paste0(path, "H1cv.rds"))
# H1c.v <- readRDS(paste0(path, "H1cv.rds"))
summary(H1c.v)
bayes_R2(H1c.v)
plot(conditional_effects(H1c.v
						 , re_formula = NA
						 )
	 )

H1d <- (
	dat
	%>% brms::brm(
		formula = (
			self_accz ~
				condition * expected_skillz # expected skill
			+ (1 + condition * expected_skillz | participant)
		)
		, data = .
		, prior = c(prior(normal(0,10), class = b))
		, silent = F
		, refresh = 20
		, chains = floor(parallel::detectCores())
		, iter = 5000
		, cores = floor(parallel::detectCores())
	)
)
saveRDS(H1d, paste0(path, "H1d.rds"))
# H1d <- readRDS(paste0(path, "H1d.rds"))
summary(H1d)
bayes_R2(H1d)
plot(conditional_effects(H1d
						 , effects = 'expected_skillz:condition'
						 , re_formula = NA
						 )
	 )


# all those things again for comparison:

# a. physical condition: accuracy rating ~ actual skill
summary(H1a)
# b. physical condition: accuracy rating ~ expected skill
summary(H1b)
# c. imagery condition: accuracy rating ~ expected skill
summary(H1c)
# c.v, imagery condition: accuracy rating ~ expected skill * vividness
summary(H1c.v)

bayes_R2(H1a)
bayes_R2(H1b)
bayes_R2(H1c)
bayes_R2(H1c.v)

# a and b have nearly the same R2 value (not significantly different) but
# c has a significantly lower R2 — but still a positive correlation and not
# that different from the physical condition

plot(conditional_effects(H1a
						 , re_formula = NA
						 )
	 )
plot(conditional_effects(H1b
						 , re_formula = NA
						 )
	 )
plot(conditional_effects(H1c
						 , re_formula = NA
						 )
	 )
plot(conditional_effects(H1c.v
						 , effects = 'expected_skillz:vividz'
						 , re_formula = NA
						 )
	 )[[1]] + theme_minimal()

# look at between condition difference:

# d. imagery vs. physical: accuracy rating ~ expected skill
summary(H1d)
bayes_R2(H1d)
plot2brm <- plot(conditional_effects(H1d
						 , effects = 'expected_skillz:condition'
						 , re_formula = NA
						 )
				 )[[1]]

# create nice plot for publication
plot2df <- as.data.frame(plot2brm$data)
levels(plot2df$effect2__) <- c("Imagery","Overt")
(
	ggplot(
		plot2df
		, aes(
			x = effect1__
			, y = estimate__
			, color = effect2__
		)
	)
	+ geom_line(
		aes(
			color = effect2__
			, linetype = effect2__
		)
		, size = .5
	)
	+ geom_ribbon(
		aes(
			ymin=lower__
			, ymax=upper__
			, fill = effect2__
			, alpha = effect2__
		)
		, color = NA
	)
	+ scale_linetype_manual(
		values = c(
			"Imagery" = "dashed"
			, "Overt" = "solid"
		)
	)
	+ scale_alpha_manual(
		values = c(
			"Imagery" = .2
			, "Overt" = .4
		)
	)
	+ labs(
		x = "Expected Performance (z-score)"
		, y = "Self-rated Accuracy (z-score)"
		, color  = "Group"
		, linetype = "Group"
		, fill = "Group"
		, alpha = "Group"
	)
	# + scale_x_continuous(
	# 	breaks=c(1,20,40,60,80,100,120)
	# )
	+ theme_minimal()
	+ theme(
		text=element_text(
			family="Times"
			, size = 7
		)
		# , legend.position="none"
	)
	+ theme(axis.title.x = element_text(margin = margin(t = 5)))
	+ theme(legend.title=element_blank())
)
ggsave(
	file = paste0('_Scripts/_plots/Fig3_90mm.jpeg')
	, plot = last_plot()
	, units = "mm"
	, width = 90
	, height = 60
)


#### Hypothesis 2 ####

# If motor representations cannot be formed internally, then performing a novel
# task via imagery repeatedly should not result in changing accuracy ratings
# over time compared to random tasks (the random condition) — but during
# physical trials, it will. After all, no representation can be formed and
# therefore updated from some kind of error signal.

H2 <- (
	dat
	%>% brms::brm(
		formula = (
			self_accz ~ trialz * rep * condition
			+ (1 + trialz * rep * condition | participant)
		)
		, data = .
		, prior = c(prior(normal(0,10), class = b))
		, silent = F
		, refresh = 20
		, chains = floor(parallel::detectCores())
		, iter = 10000
		, cores = floor(parallel::detectCores())
		, control = list(adapt_delta = .95)
	)
)
saveRDS(H2, paste0(path, "H2.rds"))
summary(H2)
bayes_R2(H2)
# H2 <- readRDS(paste0(path, "H2.rds"))

conditions <- data.frame(condition = c('imagery','physical'))
plot3brm <- conditional_effects(H2
								, effects = "trialz:rep"
								, conditions = conditions
								, re_formula = NA
								)

# create nice plot for publication
plot3df <- as.data.frame(plot3brm[[1]])
levels(plot3df$effect2__) <- c("Random","Repeated")
cond_names <- c(
	'imagery' = "Imagery"
	, 'physical' = "Overt"
)
(
	ggplot(
		plot3df
		, aes(
			x = (effect1__ * attr(scale(dat$trial), "scaled:scale"))
			+ attr(scale(dat$trial), "scaled:center")
			, y = estimate__
			, color = effect2__
		)
	)
	+ facet_grid(
		cols = vars(condition)
		, labeller = as_labeller(cond_names)
	)
	+ geom_line(
		aes(
			color = effect2__
			, linetype = effect2__
		)
		, size = .5
	)
	+ geom_ribbon(
		aes(
			ymin=lower__
			, ymax=upper__
			, fill = effect2__
			, alpha = effect2__
		)
		, color = NA
	)
	+ scale_linetype_manual(
		values = c(
			"Random" = "dashed"
			, "Repeated" = "solid"
		)
	)
	+ scale_alpha_manual(
		values = c(
			"Random" = .2
			, "Repeated" = .4
		)
	)
	+ labs(
		x = "Trial"
		, y = "Self-rated Accuracy (z-score)"
		, color  = "Trial Type"
		, linetype = "Trial Type"
		, fill = "Trial Type"
		, alpha = "Trial Type"
	)
	+ scale_x_continuous(
		breaks=c(1,20,40,60,80,100,120)
	)
	+ theme_minimal()
	+ theme(
		text=element_text(
			family="Times"
			, size = 7
		)
		# , legend.position="none"
	)
	+ theme(axis.title.x = element_text(margin = margin(t = 5)))
	+ theme(legend.title=element_blank())
)
ggsave(
	file = paste0('_Scripts/_plots/Fig4_90mm.jpeg')
	, plot = last_plot()
	, units = "mm"
	, width = 90
	, height = 60
)

#### Hypothesis 3 ####

# scale everything manually
dat <- (
	dat
	%>% dplyr::mutate(
		mt_z = scale(mt)[,1]
		, stim_mt_z = scale(stimulus_mt)[,1]
		, complexityz = scale(complexity)[,1]
	)
)

H3 <-(
	dat # using all data
	%>% brms::brm(
		formula = (
			mt_z ~ condition * stim_mt_z * complexityz
			+ (1 + condition * stim_mt_z * complexityz | participant)
		)
		, data = .
		, prior = c(prior(normal(0,10), class = b))
		, silent = F
		, refresh = 20
		, chains = floor(parallel::detectCores())
		, iter = 10000
		, cores = floor(parallel::detectCores())
	)
)
saveRDS(H3, paste0(path, "H3.rds"))
# H3 <- readRDS(paste0(path, "H3.rds"))
summary(H3)
bayes_R2(H3)

plot(conditional_effects(H3
						 , re_formula = NA
						 )
	 , ask = FALSE
	 )

conditions <- data.frame(condition = c('physical','imagery'))
plot4brm <- plot(conditional_effects(H3
									 , effects = 'stim_mt_z:complexityz'
									 , conditions = conditions
									 , re_formula = NA
									 )
				 , ask = FALSE
				 )[[1]]

# create nice plot for publication
plot4df <- as.data.frame(plot4brm$data)
levels(plot4df$effect2__) <- c("1 SD","Mean","-1 SD")
cond_names <- c(
	'imagery' = "Imagery"
	, 'physical' = "Overt"
)
(
	ggplot(
		plot4df
		, aes(
			# x = effect1__
			# , y = estimate__
			x = effect1__ * attr(scale(dat$stimulus_mt), "scaled:scale") + attr(scale(dat$stimulus_mt), "scaled:center")
			, y = estimate__ * attr(scale(dat$mt), "scaled:scale") + attr(scale(dat$mt), "scaled:center")
			, color = effect2__
		)
	)
	+ facet_grid(
		cols = vars(condition)
		, labeller = as_labeller(cond_names)
	)
	+ geom_line(
		aes(
			color = effect2__
			, linetype = effect2__
		)
		, size = .5
	)
	+ geom_ribbon(
		aes(
			ymin=lower__ * attr(scale(dat$mt), "scaled:scale") + attr(scale(dat$mt), "scaled:center")
			, ymax=upper__ * attr(scale(dat$mt), "scaled:scale") + attr(scale(dat$mt), "scaled:center")
			, fill = effect2__
			, alpha = effect2__
		)
		, color = NA
	)
	+ scale_linetype_manual(
		values = c(
			"1 SD" = "dashed"
			, "Mean" = "solid"
			, "-1 SD" = "dotted"
		)
	)
	+ scale_alpha_manual(
		values = c(
			"1 SD" = .4
			, "Mean" = .3
			, "-1 SD" = .2
		)
	)
	+ labs(
		x = "Stimulus Animation Time (sec)"
		, y = "Participant Movement Time (sec)"
		, color  = "Complexity"
		, linetype = "Complexity"
		, fill = "Complexity"
		, alpha = "Complexity"
	)
	# + scale_x_continuous(
	# 	breaks=c(1,20,40,60,80,100,120)
	# )
	+ theme_minimal()
	# + theme(
	# 	legend.position = c(.9,.15)
	# )
	+ theme(
		text=element_text(
			family="Times"
			, size = 7
			)
		# , legend.position="none"
	)
	+ theme(axis.title.x = element_text(margin = margin(t = 5)))
	# + theme(legend.title=element_blank())
)
# saves last plot
ggsave(
	file = paste0('_Scripts/_plots/Fig5_90mm.jpeg')
	, plot = last_plot()
	, units = "mm"
	, width = 90
	, height = 60
)
