hot_lat = 50
hot_long = 152

# hot_lat = 92
# hot_long = 162

# hot_lat = 50
# hot_long = 28

# hot_lat = 50
# hot_long = -28

# hot_lat = 92
# hot_long = -90

# hot_lat = 0
# hot_long = 0

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

#viz
(
	topo_dat
	%>% ggplot()
	+ geom_point(
		aes(
			x = x
			, y = y
			, colour = value
		)
		, size = 4
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

#add noisy replicates
(
	topo_dat
	%>% expand_grid(
		rep = 1:10
	)
	%>% mutate(
		value = value + rnorm(n(),0,.1)
	)

) -> topo_dat_reps

#viz
(
	topo_dat_reps
	%>% ggplot()
	+ geom_jitter(
		aes(
			x = x
			, y = y
			, colour = value
		)
		, size = 2
		, alpha = .5
		, width = 5
		, height = 5
	)
)

#fit
library(mgcv)
fit = gam(
	data = topo_dat_reps
	, formula = value ~ te(lat,long,bs='sos',d=2,k=30)
)

#viz
(
	topo_dat
	%>% mutate(
		preds = predict(fit,newdata=topo_dat)
	)
	%>% ggplot()
	+ geom_point(
		aes(
			x = x
			, y = y
			, colour = value
		)
		, size = 4
	)
)
