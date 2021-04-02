######################
### DEMI plot GAMs ###
######################

library(tidyverse)

scale_to_0range = function(x,range=1){
	x = x - min(x)
	x = x/max(x)
	x = x-.5
	x = x * range
	return(x)
}

#in case the preds haven't been loaded
preds_dat = readRDS('_rds/preds_dat.rds')



# example main-effect-topo plot: block topo ----

#get the data to plot
(
	# start with the full preds
	preds_dat
	# group by the variables you want AND sample
	%>% group_by(
		lat
		, long
		, block
		, sample
	)
	# collapse to a mean, dropping sample from the grouping thereafter
	%>% summarise(
		value = mean(value)
		, .groups = 'drop_last'
	)
	# compute uncertainty intervals and midpoint (using sample==0 for midpoint)
	%>% summarise(
		lo = quantile(value,.03/2) #97%ile lower-bound
		, hi = quantile(value,1-.03/2) #97%ile upper-bound
		, mid = value[sample==0]
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
		, x_scaled = scale_to_0range(x,10)
		, y_scaled = scale_to_0range(y,10)

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

		#block is going to be mapped to the x-axis of each sub-panel, so re-map it to have a range of 1
		, block_scaled = scale_to_0range(block,1)

		# just like we did above for the y-axis, get the global x-axis position given the subpanel location and scaled x-axis data
		, to_plot_x = x_scaled + block_scaled

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
		data = . %>% group_keys(lat,long,x,y,x_scaled,y_scaled)
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
			, fill = factor(1)# can't remember why this is here
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
	#now save (weird to have a + instead of %>%, I know)
	+ ggsave(
		file = '_Scripts/_plots/examples_block_topo.pdf'
		, width = 10
		, height = 10
	)
)
