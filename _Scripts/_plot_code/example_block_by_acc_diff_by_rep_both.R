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


#get the data to plot
(
	# start with the full preds
	preds_dat
	# filter to factors of interest
	%>% filter(
		band == 'beta'
		, epoch == 2
		, group == 'physical'
	)
	# group by the variables you want AND sample
	%>% group_by(
		lat
		, long
		, block
		, rep # this factor gets each plotted
		, sample
		, accuracy # this factor gets differenced
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

		#block is going to be mapped to the x-axis of each sub-panel, so re-map it to have a range of 1
		, block_scaled = scale_to_0range(block,1)

		# just like we did above for the y-axis, get the global x-axis position given the subpanel location and scaled x-axis data
		, to_plot_x = x_scaled + block_scaled

	)
	#save as new object so we don't have to re-run the above if we make mistakes or tweaks below
) -> ready_to_plot


#compute some quantities for the axes panel
#whole plot is 10x10 units centered on zero
y_axis_y_offset = -4
y_axis_x_offset = -4
y_axis_dat = tibble(
	label = seq(ready_to_plot$min_lo[1],ready_to_plot$min_lo[1]+ready_to_plot$range_[1],length.out=3)
	, y_scaled = c(0,.5,1)
	, to_plot_y = y_scaled -.5 + y_axis_y_offset
	, to_plot_x = rep(0,3) + y_axis_x_offset
)
x_axis_y_offset = y_axis_y_offset-.5
x_axis_x_offset = y_axis_x_offset+.5
x_axis_dat = tibble(
	label = c('1','','','','5','')
	, x_scaled = seq(0,1,length.out=6)
	, to_plot_x = x_scaled -.5 + x_axis_x_offset
	, to_plot_y = rep(0,length(label)) + x_axis_y_offset
)
axis_title_dat = tibble(
	label = c('Relative power\n(log-dB)','Block')
	, x = c(y_axis_x_offset-.5,x_axis_x_offset)
	, y = c(y_axis_y_offset,x_axis_y_offset-.3)
	, angle = c(90,0)
	, hjust = c('center','center')
	, vjust = c('bottom','top')
)

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

	# y-axis line
	+ geom_line(
		data = y_axis_dat
		, aes(
			x = to_plot_x
			, y = to_plot_y
		)
	)
	# y-axis ticks
	+ geom_line(
		data = (
			y_axis_dat
			%>% mutate(
				xmin = to_plot_x-.1
				, xmax = to_plot_x
			)
			%>% select(-to_plot_x)
			%>% pivot_longer(
				cols = c(xmin,xmax)
				, values_to = 'to_plot_x'
			)
		)
		, aes(
			x = to_plot_x
			, y = to_plot_y
			, group = label
		)
	)
	# y-axis labels
	+ geom_text(
		data = y_axis_dat
		, aes(
			x = to_plot_x-.1
			, y = to_plot_y
			, label = paste(format(label,digits=2),' ')
		)
		, hjust = 'right'
		, size = 2
	)
	# x-axis line
	+ geom_line(
		data = x_axis_dat
		, aes(
			x = to_plot_x
			, y = to_plot_y
		)
	)
	# x-axis ticks
	+ geom_line(
		data = (
			x_axis_dat
			%>% mutate(
				ymin = to_plot_y-.1
				, ymax = to_plot_y
			)
			%>% select(-to_plot_y)
			%>% pivot_longer(
				cols = c(ymin,ymax)
				, values_to = 'to_plot_y'
			)
		)
		, aes(
			x = to_plot_x
			, y = to_plot_y
			, group = interaction(label,x_scaled)
		)
	)
	#x-axis labels
	+ geom_text(
		data = x_axis_dat
		, aes(
			x = to_plot_x
			, y = to_plot_y-.15
			, label = label
		)
		, vjust = 'top'
		, size = 2
	)
	# axis titles
	+ geom_text(
		data = axis_title_dat
		, aes(
			x = x
			, y = y
			, label = label
			, angle = angle
			, hjust = hjust
			, vjust = vjust
		)
		, size = 3
	)

	# render the uncertainty intervals
	# (could use geom_errorbar instead)
	+ geom_ribbon(
		mapping = aes(
			x = to_plot_x
			, ymin = to_plot_lo
			, ymax = to_plot_hi
			, group = interaction(lat,long,rep)
			, fill = rep
		)
		, alpha = .5
	)
	# render the predictions for the mean
	# (could use geom_point instead or in addition)
	+ geom_line(
		mapping = aes(
			x = to_plot_x
			, y = to_plot_mid
			, group = interaction(lat,long,rep)
			, colour = rep
		)
		, alpha = .5
	)
	+ coord_equal() #important to make subpanel locations accurate
	+ theme(
		# legend.position = 'none'
		# , legend.justification = c(0,0)
		# , legend.title = element_blank()
		axis.title = element_blank()
		, axis.ticks = element_blank()
		, axis.text = element_blank()
		, panel.grid = element_blank()
		, panel.background = element_rect(fill='transparent',colour='grey90')
	)
)

#now save
ggsave(
	file = '_plots/physical_beta2_by_rep_topo_both.pdf'
	, width = 10
	, height = 10
)
