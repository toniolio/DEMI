library(tidyverse)
library(mgcv) #gam stuff
library(modelr) #for data_grid
library(Rfast) #for mat.mult & rmvnorm
library(plotly)
library(circular)

scale_to_0range = function(x,range=1){
	x = x - min(x)
	x = x/max(x)
	x = x-.5
	x = x * range
	return(x)
}

#get the preds:
preds_dat = readRDS('_rds/preds_dat_re.rds')

#make fake data as before:
hot_lat = 50
hot_long = 152
(
	preds_dat
	%>% modelr::data_grid(
		nesting(
			lat
			, long
		)
	)
	%>% mutate(
		x = lat*cos(long*(pi/180))
		, y = lat*sin(long*(pi/180))
		, value = -sqrt((
			((x-x[(lat==hot_lat)&(long==hot_long)])^2)
			+ ((y-y[(lat==hot_lat)&(long==hot_long)])^2)
		))
		, value = value - min(value)
		, value = value/max(value)
	)
) -> topo_dat

#viz 2D
(
	topo_dat
	%>% ggplot()
	+ geom_jitter(
		aes(
			x = x
			, y = y
			, colour = value
		)
		, size = 4
		# , alpha = .5
		, width = 5
		, height = 5
	)
	+ geom_text(
		aes(
			x = x
			, y = y
			# , label = paste(round(x),round(y))
			, label = paste(lat,long)
		)
	)
)


#show the *proper* 3D:
(
	topo_dat
	%>% mutate(lat=90-lat) #this is the correction
	%>% mutate(
		x = cos(rad(lat)) * cos(rad(long))
		, y = cos(rad(lat)) * sin(rad(long))
		, z = sin(rad(lat))
	)
	%>% plot_ly(
		x = ~x
		, y = ~y
		, z = ~z
		, marker = list(color=~value)
	)
)

#show what mgcv saw with the *bad* 3D:
(
	topo_dat
	%>% mutate(
		x = cos(rad(lat)) * cos(rad(long))
		, y = cos(rad(lat)) * sin(rad(long))
		, z = sin(rad(lat))
	)
	%>% plot_ly(
		x = ~x
		, y = ~y
		, z = ~z
		, marker = list(color=~value)
	)
)
