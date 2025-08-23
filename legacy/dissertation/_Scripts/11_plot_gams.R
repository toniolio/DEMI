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

#see examples in `_plot_code`

source('_plot_code/example_block_topo.r')

library(diffr)
diffr(
	'_plot_code/example_block_by_rep_topo_both.r'
	, '_plot_code/example_block_by_accuracy_topo_both.r'
)
