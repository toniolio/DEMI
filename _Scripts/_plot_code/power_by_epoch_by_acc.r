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
preds_dat = readRDS('_rds/preds_dat_re_2.rds')

grp = 'physical' # physical, imagery
bnd = 'theta' # theta, alpha, beta

approach = 1 # approach to inference
# 1 means non-overlapping error bars
# 2 means the difference not being credibly zero

plotCI = .05 # determines error bars by dividing by 2
sigCI = .05 # set "alpha" for "significance" tests (also divided by 2)
# so, .1 = 95% CI; .05 = 97.5%, .02 = 99%, .01 = 99.5%, .002 = 99.9%

#get the data to plot
(
	# start with the full preds
	preds_dat
	# re-label accuracy
	%>% mutate(
		accuracy = case_when(
			accuracy==max(accuracy) ~ 'High'
			, accuracy==min(accuracy) ~ 'Low'
		)
	)
	# filter to to conditions of interest
	%>%filter(
		group == grp # physical, imagery
		, band == bnd # theta, alpha, beta
		# , rep == 'repeated' # random, repeated
		# , block < max(block) # remove final block — better group comparison
	)
	# group by the variables you want AND sample
	%>% group_by(
		lat
		, long
		, epoch
		, accuracy
		, sample
	)
	# collapse to a mean, dropping sample from the grouping thereafter
	%>% summarise(
		value = mean(value)
		, .groups = 'drop_last'
	)
	# compute uncertainty intervals and midpoint (using sample==0 for midpoint)
	%>% summarise(
		lo = quantile(value,plotCI/2) # CI lower-bound
		, hi = quantile(value,1-plotCI/2) #CI upper-bound
		, mid = value[sample==0]
		, samples = list(tibble(value,sample)) #hack for significance labels
		, .groups = 'drop'
	)
	#save output so we don't have to re-compute if we make mistakes or tweaks below
) -> to_plot

#get "significance"
(
	to_plot
	%>% select(-lo,-hi,-mid)
	%>% unnest(samples)
	# group by all but the one you want to collapse to a difference
	%>% group_by(
		lat
		, long
		, epoch
		, sample
	)
	%>% summarise(
		value = value[accuracy=='High'] - value[accuracy=='Low']
		, .groups = 'drop_last'
	)
	%>% summarise(
		lo = quantile(value,sigCI/2) #97.5%ile lower-bound
		, hi = quantile(value,1-sigCI/2) #97.5%ile upper-bound
		, .groups = 'drop'
	)
	%>% mutate(
		diff_sig = case_when(
			(lo>0)|(hi<0) ~ T
			, T ~ F
		)
	)
	%>% group_by(
		lat
		, long
	)
	%>% summarise(
		diff_sig = any(diff_sig)
		, .groups = 'drop'
	)
	%>% right_join(to_plot,by=c('lat','long'))
	%>% select(-samples)
	%>% group_by(
		lat
		, long
		, epoch
	)
	%>% mutate(
		indiv_sig = case_when(
			(lo[accuracy=='High']>hi[accuracy=='Low'])
			| (lo[accuracy=='Low']>hi[accuracy=='High']) ~ T
			, T ~ F
		)
	)
	%>% group_by(
		lat
		, long
	)
	%>% mutate(
		indiv_sig = any(indiv_sig)
	)
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
		, zero_scaled = (0-min_lo)/range_ - .5 #############For when zero is interesting!

		# now get the global y-position given the subpanel location and subpanel's scaled y-axis data
		, to_plot_lo = y_scaled + lo_scaled
		, to_plot_hi = y_scaled + hi_scaled
		, to_plot_mid = y_scaled + mid_scaled
		, to_plot_zero = y_scaled + zero_scaled

		####
		# Content of this section will change depending on variables in the plot
		####

		#epoch is going to be mapped to the x-axis of each sub-panel, so re-map it to have a range of 1
		#now make them equally spaced between -.5 and .5 (but not AT those values)
		, epoch_scaled = case_when(
			epoch=="during" ~ -.25
			, epoch=="after" ~ .25
		)

		#also want to shift the x-axis locations by rep a smidge
		, acc_scaled = case_when(
			accuracy=="Low" ~ -.075
			, accuracy=="High" ~ .075
		)

		# just like we did above for the y-axis, get the global x-axis position given the subpanel location and scaled x-axis data
		, to_plot_x = x_scaled + epoch_scaled + acc_scaled

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
	label = c('During','After')
	, x_scaled = c(-.25,.25)
	, to_plot_x = x_scaled + x_axis_x_offset
	, to_plot_y = rep(0,length(label)) + x_axis_y_offset
)
axis_title_dat = tibble(
	label = c('Relative power\n(log-dB)','Epoch')
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
		data = . %>% group_keys(lat,long,x_scaled,y_scaled,diff_sig,indiv_sig)
		, aes(
			xmin = x_scaled-.5
			, xmax = x_scaled+.5
			, ymin = y_scaled-.5
			, ymax = y_scaled+.5
			, group = interaction(lat,long)
			, fill = if (approach==1) indiv_sig else diff_sig
		)
		# , fill = 'grey90'
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
	# (could use geom_ribbon instead)
	+ geom_errorbar(
		mapping = aes(
			x = to_plot_x
			, ymin = to_plot_lo
			, ymax = to_plot_hi
			, group = interaction(lat,long,accuracy)
			# , fill = accuracy
		)
		, alpha = .5
		, width = .125 #for bar, not ribbon
	)
	# render the predictions for the mean
	# (could use geom_line instead or in addition)
	+ geom_point(
		mapping = aes(
			x = to_plot_x
			, y = to_plot_mid
			, group = interaction(lat,long,accuracy)
			, colour = accuracy
		)
		, alpha = .5
	)
	# line at zero
	+ geom_line(
		data = (
			ready_to_plot
			%>% group_keys(lat,long,x_scaled,to_plot_zero)
			%>% mutate(
				xmin = x_scaled-.5
				, xmax = x_scaled+.5
			)
			%>% select(-x_scaled)
			%>% pivot_longer(
				cols = c(xmin,xmax)
				, values_to = 'to_plot_x'
			)
		)
		, aes(
			x = to_plot_x
			, y = to_plot_zero
			, group = interaction(lat,long)
		)
		, colour = 'white'
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
	+ labs(title = bnd
		   # , tag = grp
		   )
	+ scale_fill_manual(
		values = c('grey50','grey90')
		, breaks = c(T,F)
		, guide = F
	)
)

#now save
ggsave(
	file = paste0('_plots/',grp,'_',bnd,'_by_epoch_by_acc.pdf')
	, width = 10
	, height = 10
)
