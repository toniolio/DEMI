
library(tidyverse)
library(mgcv) #gam stuff
library(modelr) #for data_grid
library(Rfast) #for mat.mult & rmvnorm


#generating predictions (and samples) ----

#get the model
path <- "_Scripts/_rds/"
gam = readRDS(paste0(path, "gam_re.rds"))

#extract the data and get the unique combinations of predictors
(
	gam$model
	%>% modelr::data_grid(
		group
		, participant
		, band
		, epoch
		, rep
		, block
		# get just the first and last quintile points for accuracy
		, accuracy = quantile(accuracy,probs=c(.2,.8))
		# nesting gets just those combos in the data
		, nesting(
			#round both to eliminate numeric drift from ealier conversions
			round(lat)
			, round(long)
		)
	)
	#renaming
	%>% rename(
		lat = `round(lat)`
		, long = `round(long)`
	)
	%>% filter(
		!((lat==0) & (long==180))
		, participant==participant[1]
	)
) -> preds_dat

#double-checking the locations
(
	preds_dat
	%>% group_keys(lat,long)
	%>% arrange(lat,long)
	%>% print(n=40)
	%>% invisible()
)

# get the model matrix implied by the preds & model structure
mm = mgcv::predict.bam(
	gam
	, newdata = preds_dat
	, type = 'lpmatrix'
	, discrete = FALSE
	, exclude = 's(participant)'
)

#get the coefficients and covariance matrix
f = coef(gam)
v = vcov(gam)

# get model predictions for the condition means
# (weird name will make sense later)
preds_dat$`...0` = Rfast::mat.mult(mm , matrix(f,nrow=length(f)))[,1]

# sample the model to make uncertainty intervals easier
num_samples = 1e3
system.time(sample_coefs <- Rfast::rmvnorm(num_samples, f, v))
system.time(sample_vals <- Rfast::mat.mult(mm,t(sample_coefs) ) )

# bind together with preds_dat
(
	preds_dat
	%>% bind_cols(
		as_tibble(sample_vals,.name_repair='unique')
		#names on this new tibble are `...X` where X is 10:(num_samples-10)
	)
	%>% pivot_longer(
		cols = contains('...')
		, names_to = 'sample'
		, names_prefix = '...'
		, names_transform = list(sample=as.integer)
	)
) -> preds_dat

#ensuring geographic lat/long (important for plots)
preds_dat$lat = 90-preds_dat$lat

#save to file
saveRDS(
	preds_dat
	, paste0(path,'preds_dat_re.rds')
)
