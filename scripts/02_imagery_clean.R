###############################################
### Imagery clean up due to TraceLab errors ###
###############################################


library(tidyverse)
library(cmdstanr)

(
	'_Scripts/_rds/bdat.rds'
	%>% readRDS()
	%>% dplyr::filter(
		condition == "imagery"
	)
	%>% dplyr::rename(
		obs = mt
		, stim_speed = stimulus_mt
	)
) -> dat


(
	dat
	%>% ggplot()
	+ geom_point(
		mapping = aes(
			x = stim_speed
			, y = obs
			# , colour = is_fast_true
		)
		, alpha = .5
	)
)

(
	dat
	%>% dplyr::mutate(
		obs = base::scale(obs,center=F)
	)
	%>% ggplot()
	+ geom_point(
		mapping = aes(
			x = stim_speed
			, y = obs
			# , colour = is_fast_true
		)
		, alpha = .5
	)
	+ geom_vline(xintercept=2.4)
)
(
	dat
	%>% dplyr::mutate(
		obs = base::scale(obs,center=F)
	)
	%>% dplyr::filter(
		stim_speed>2.4
		, obs>.5
	)
	%>% dplyr::summarise(
		m = mean(obs)
		, s = sd(obs)
	)
)

mod = cmdstanr::cmdstan_model('_Scripts/_stan/outlier_mixture.stan')
fit = mod$sample(
	data = tibble::lst(
		obs = dat$obs
		, stim_speed = dat$stim_speed
		, num_obs = length(dat$obs)
	)
	, chains = parallel::detectCores()/2
	, parallel_chains = parallel::detectCores()/2
	, iter_warmup = 1e3
	, iter_sampling = 1e3
	, refresh = 20
	# , metric = 'dense_e'
)
fit$cmdstan_diagnose()


#get is_fast_est
(
	fit$summary(variables='is_fast','mean')
	%>% dplyr::mutate(
		is_fast_est = dplyr::case_when(
			mean>.5 ~ T
			, mean<=.5 ~ F
		)
	)
	%>% dplyr::select(is_fast_est)
	%>% dplyr::bind_cols(dat)
) -> dat_with_est

#quick viz:
(
	dat_with_est
	%>% (function(x){
		print(mean(x$is_fast_est))
		return(x)
	})
	%>% ggplot()
	+ geom_point(
		mapping = aes(
			x = stim_speed
			, y = obs
			, colour = is_fast_est
		)
		, alpha = .5
		, size = 3
	)
)

#add good_cdf
(
	fit$summary(variables='good_cdf','median')
	%>% dplyr::rename(good_cdf=median)
	%>% dplyr::select(good_cdf)
	%>% dplyr::bind_cols(dat_with_est)
) -> dat_with_est

# add the final label
criterion = 5
(
	dat_with_est
	%>% dplyr::mutate(
		q_unreasonable = abs(qnorm(good_cdf))>criterion
	)
	%>% (function(x){
		(
			x
			%>% dplyr::filter(!is_fast_est)
			%>% dplyr::summarize(
				observed = mean(q_unreasonable)*100
				, expected = ((1-pnorm(criterion))+(pnorm(-criterion)))*100
			)
			%>% print()
		)
		return(x)
	})
	%>% dplyr::mutate(
		good = !(is_fast_est | q_unreasonable)
	)
) -> dat_with_est

#plot:
(
	dat_with_est
	%>% dplyr::mutate(
		label = dplyr::case_when(
			is_fast_est ~ 'fast'
			, q_unreasonable ~ 'unreasonable'
			, T ~ 'good'
		)
		, label2 = label
		, stim_rounded = round(stim_speed*2)
		, stim_rounded = dplyr::case_when(
			stim_rounded==6 ~ 5
			, T ~ stim_rounded
		)
		, stim_rounded = stim_rounded/2
	)
	# %>% ggplot()
	# + geom_point(
	# 	mapping = aes(
	# 		x = stim_speed
	# 		, y = obs
	# 		, colour = factor(stim_rounded)
	# 	)
	# ))
	%>% dplyr::group_by(
		stim_rounded
	)
	%>% dplyr::mutate(
		good_min = min(obs[label=='good'])
		, good_max = max(obs[label=='good'])
	)
	%>% dplyr::ungroup()
	%>% (function(x){
		dplyr::bind_rows(
			x
			, (
				x
				%>% dplyr::mutate(
					label = 'all_data'
				)
			)
		)
	})
	%>% dplyr::mutate(
		to_NA = label %in% c('all_data','good')
		, good_min = dplyr::case_when(
			to_NA ~ as.numeric(NA)
			, T ~ good_min
		)
		, good_max = dplyr::case_when(
			to_NA ~ as.numeric(NA)
			, T ~ good_max
		)
	)
	# %>% select(label,to_NA,good_min)
	# %>% summary())
	# %>% View())
	%>% (function(x){
		(
			x
			%>% ggplot()
			+ facet_wrap(~label,ncol=4)
			+ geom_errorbar(
				data = (
					x
					%>% dplyr::group_by(label,stim_rounded)
					%>% dplyr::summarise(
						good_min = good_min[1]
						, good_max = good_max[1]
						, .groups = 'drop'
					)
				)
				, mapping = aes(
					x = stim_rounded
					, ymin = good_min
					, ymax = good_max
				)
			)
			+ geom_point(
				mapping = aes(
					x = stim_speed
					, y = obs
					# , colour = label2
					, fill = label2
				)
				, alpha = .5
				, size = 3
				, shape = 21
				# , colour = 'black'
				, colour = 'transparent'
				, position='jitter'
			)
			+labs(
				color = 'label'
				, fill = 'label'
			)
		)->p
		print(p)
		return(x)
	})
) -> dat_with_est

bad_imagery_trials <- (
	dat_with_est
	%>% dplyr::filter(
		label != 'all_data'
		, (label == 'fast' | label == 'unreasonable')
	)
	%>% dplyr::select(
		figure_file
	)
)
saveRDS(bad_imagery_trials, '_Scripts/_rds/bad_imagery_trials.rds')
