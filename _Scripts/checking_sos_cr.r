library(tidyverse)
library(mgcv) #gam stuff
library(modelr) #for data_grid
library(Rfast) #for mat.mult & rmvnorm

scale_to_0range = function(x,range=1){
	x = x - min(x)
	x = x/max(x)
	x = x-.5
	x = x * range
	return(x)
}


preds_dat = readRDS('_rds/preds_dat.rds')

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




# get the model matrix implied by the preds & model structure
mm = mgcv::predict.bam(
	fit
	, newdata = topo_dat
	, type = 'lpmatrix'
	, discrete = FALSE
)

#get the coefficients and covariance matrix
f = coef(fit)
v = vcov(fit)

# get model predictions for the condition means
# (weird name will make sense later)

(
	topo_dat
	%>% mutate(
		value = Rfast::mat.mult(mm , matrix(f,nrow=length(f)))[,1]
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



topo_dat$`...0` = Rfast::mat.mult(mm , matrix(f,nrow=length(f)))[,1]


# sample the model to make uncertainty intervals easier
num_samples = 1e4
system.time(sample_coefs <- Rfast::rmvnorm(num_samples, f, v))
system.time(sample_vals <- Rfast::mat.mult(mm,t(sample_coefs) ) )

# bind together with preds_dat
(
	topo_dat
	%>% select(-value)
	%>% bind_cols(
		as_tibble(sample_vals,.name_repair='unique')
		, .name_repair = 'check_unique'
	)
	%>% pivot_longer(
		cols = contains('...')
		, names_to = 'sample'
		, names_prefix = '...'
		, names_transform = list(sample=as.integer)
	)
) -> preds_dat


#get the data to plot
(
	# start with the full preds
	preds_dat
	# group by the variables you want AND sample
	%>% group_by(
		lat
		, long
		, accuracy
		# , sample
	)
	# # collapse to a mean, dropping sample from the grouping thereafter
	# %>% summarise(
	# 	value = mean(value)
	# 	, .groups = 'drop_last'
	# )
	# compute uncertainty intervals and midpoint (using sample==0 for midpoint)
	%>% summarise(
		lo = quantile(value,.03/2) #97%ile lower-bound
		, hi = quantile(value,1-.03/2) #97%ile upper-bound
		, mid = value[sample==min(sample)]
		, .groups = 'drop'
	)
	#save output so we don't have to re-compute if we make mistakes or tweaks below
) -> to_plot

# Now we have some mutations to apply to arrange things visually
(
	to_plot
	#always ungroup first! (min/max stuff assumes no grouping)
	%>% ungroup()
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
		, min_lo = min(lo)
		, range_ = max(hi) - min_lo
		, lo_scaled = (lo-min_lo)/range_ - .5
		, hi_scaled = (hi-min_lo)/range_ - .5
		, mid_scaled = (mid-min_lo)/range_ - .5


		# now get the global y-position given the subpanel location and subpanel's scaled y-axis data
		, to_plot_lo = y_scaled + lo_scaled
		, to_plot_hi = y_scaled + hi_scaled
		, to_plot_mid = y_scaled + mid_scaled

		####
		# Content of this section will change depending on variables in the plot
		####

		#accuracy is going to be mapped to the x-axis of each sub-panel, so re-map it to have a range of 1
		, accuracy_scaled = scale_to_0range(accuracy,1)

		# just like we did above for the y-axis, get the global x-axis position given the subpanel location and scaled x-axis data
		, to_plot_x = x_scaled + accuracy_scaled

	)
	#save as new object so we don't have to re-run the above if we make mistakes or tweaks below
) -> ready_to_plot



#start plotting!

(
	ready_to_plot
	%>% ggplot()
	#create a rect around each subpanel
	+ geom_rect(
		#use the data from the pipe (.) and use group_keys to reduce to info on the subpanels
		data = . %>% group_keys(lat,long,x_scaled,y_scaled)
		, aes(
			xmin = x_scaled-.5
			, xmax = x_scaled+.5
			, ymin = y_scaled-.5
			, ymax = y_scaled+.5
			, group = interaction(lat,long)
		)
		, fill = 'grey90'
		, colour = 'transparent'
	)

	# render the uncertainty intervals
	# (could use geom_errorbar instead)
	+ geom_ribbon(
		mapping = aes(
			x = to_plot_x
			, ymin = to_plot_lo
			, ymax = to_plot_hi
			, group = interaction(lat,long)
		)
		, alpha = .5
	)
	# render the predictions for the mean
	# (could use geom_point instead or in addition)
	+ geom_line(
		mapping = aes(
			x = to_plot_x
			, y = to_plot_mid
			, group = interaction(lat,long)
		)
		, alpha = .5
	)
	+ coord_equal() #important to make subpanel locations accurate
	+ theme(
		legend.position = 'none'
		, legend.justification = c(0,0)
		, legend.title = element_blank()
		, axis.title = element_blank()
		, axis.ticks = element_blank()
		, axis.text = element_blank()
		, panel.grid = element_blank()
		, panel.background = element_rect(fill='transparent',colour='grey90')
	)
)
