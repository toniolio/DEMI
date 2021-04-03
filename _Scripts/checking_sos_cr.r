# hot_lat = 50
# hot_long = 152

hot_lat = 92
hot_long = 162

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
	%>% expand_grid(
		accuracy = seq(-1,1,by=.1)
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
		, value = value*dnorm(accuracy,0,.3)
	)
) -> topo_dat

#viz
(
	topo_dat
	%>% mutate(

		####
		# Content of this section will should remain the same for all topo plots
		####

		#polar to cartesian then re-scaled versions
		x = lat*cos(long*(pi/180))
		, y = lat*sin(long*(pi/180))
		, x_scaled = scale_to_0range(x,9) #9 makes the plot 10x10 bc 1x1 panels will be centered on these
		, y_scaled = scale_to_0range(y,9)

		#here we find the min & max if we were to plot all the data in one panel
		, min_ = min(value)
		, range_ = max(value) - min_
		, value_scaled = (value-min_)/range_ - .5


		# now get the global y-position given the subpanel location and subpanel's scaled y-axis data
		, to_plot_y = y_scaled + value_scaled

		####
		# Content of this section will change depending on variables in the plot
		####

		#accuracy is going to be mapped to the x-axis of each sub-panel, so re-map it to have a range of 1
		#now make them equally spaced between -.5 and .5 (but not AT those values)
		, accuracy_scaled = scale_to_0range(accuracy)

		# just like we did above for the y-axis, get the global x-axis position given the subpanel location and scaled x-axis data
		, to_plot_x = x_scaled + accuracy_scaled

	)
	%>% ggplot()
	+ geom_line(
		aes(
			x = to_plot_x
			, y = to_plot_y
			, group = interaction(lat,long)
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
	%>% mutate(

		####
		# Content of this section will should remain the same for all topo plots
		####

		#polar to cartesian then re-scaled versions
		x = lat*cos(long*(pi/180))
		, y = lat*sin(long*(pi/180))
		, x_scaled = scale_to_0range(x,9) #9 makes the plot 10x10 bc 1x1 panels will be centered on these
		, y_scaled = scale_to_0range(y,9)

		#here we find the min & max if we were to plot all the data in one panel
		, min_ = min(value)
		, range_ = max(value) - min_
		, value_scaled = (value-min_)/range_ - .5


		# now get the global y-position given the subpanel location and subpanel's scaled y-axis data
		, to_plot_y = y_scaled + value_scaled

		####
		# Content of this section will change depending on variables in the plot
		####

		#accuracy is going to be mapped to the x-axis of each sub-panel, so re-map it to have a range of 1
		#now make them equally spaced between -.5 and .5 (but not AT those values)
		, accuracy_scaled = scale_to_0range(accuracy)

		# just like we did above for the y-axis, get the global x-axis position given the subpanel location and scaled x-axis data
		, to_plot_x = x_scaled + accuracy_scaled

	)
	%>% ggplot()
	+ geom_line(
		aes(
			x = to_plot_x
			, y = to_plot_y
			, group = interaction(lat,long,rep)
		)
	)
)

#fit
library(mgcv)
fit = gam(
	data = topo_dat_reps
	, formula = value ~ te(lat,long,accuracy,bs=c('sos','cr'),d=c(2,1),k=c(30,10))
)

#viz
(
	topo_dat
	%>% mutate(
		value = predict(fit,newdata=topo_dat)

		####
		# Content of this section will should remain the same for all topo plots
		####

		#polar to cartesian then re-scaled versions
		, x = lat*cos(long*(pi/180))
		, y = lat*sin(long*(pi/180))
		, x_scaled = scale_to_0range(x,9) #9 makes the plot 10x10 bc 1x1 panels will be centered on these
		, y_scaled = scale_to_0range(y,9)

		#here we find the min & max if we were to plot all the data in one panel
		, min_ = min(value)
		, range_ = max(value) - min_
		, value_scaled = (value-min_)/range_ - .5


		# now get the global y-position given the subpanel location and subpanel's scaled y-axis data
		, to_plot_y = y_scaled + value_scaled

		####
		# Content of this section will change depending on variables in the plot
		####

		#accuracy is going to be mapped to the x-axis of each sub-panel, so re-map it to have a range of 1
		#now make them equally spaced between -.5 and .5 (but not AT those values)
		, accuracy_scaled = scale_to_0range(accuracy)

		# just like we did above for the y-axis, get the global x-axis position given the subpanel location and scaled x-axis data
		, to_plot_x = x_scaled + accuracy_scaled

	)
	%>% ggplot()
	+ geom_line(
		aes(
			x = to_plot_x
			, y = to_plot_y
			, group = interaction(lat,long)
		)
	)
)
