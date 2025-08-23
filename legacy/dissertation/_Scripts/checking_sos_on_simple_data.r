library(tidyverse)
library(mgcv)

dat = readRDS('_rds/dat_PP_train.rds')

(
	dat
	# %>% filter(
	# 	band=='beta'
	# )
	%>% group_by(
		epoch
		, band
		, lat
		, long
		, participant
	)
	%>% summarise(
		value = mean(powerdb)
	)
) -> dat


#viz
(
	dat
	%>% group_by(
		epoch
		, band
		, lat
		, long
	)
	%>% summarise(
		value = mean(value)
	)
	%>% group_by(band,lat,long)
	%>% summarise(
		value = value[epoch=='during'] - value[epoch=='after']
	)
	%>% mutate(
		x = lat*cos(long*(pi/180))
		, y = lat*sin(long*(pi/180))
	)
	%>% ggplot()
	+ facet_wrap(~band)
	+ geom_point(
		aes(
			x = x
			, y = y
			, colour = value
		)
		, size = 4
	)
)


fit = gam(
	formula = value ~ (
		epoch * band
		+ te(
			lat,long
			, by = interaction(epoch,band)
			, bs = 'sos'
			, d = 2
			, k = 32
		)
	)
	, data = dat
)

(
	dat
	%>% group_by(epoch,band,lat,long)
	%>% group_keys()
	%>% ungroup()
	%>% mutate(
		value = predict(fit,newdata=.)
	)
	%>% group_by(band,lat,long)
	%>% summarise(
		value = value[epoch=='during'] - value[epoch=='after']
	)
	%>% mutate(
		x = lat*cos(long*(pi/180))
		, y = lat*sin(long*(pi/180))
	)
	%>% ggplot()
	+ facet_wrap(~band)
	+ geom_point(
		aes(
			x = x
			, y = y
			, colour = value
		)
		, size = 4
	)

)
