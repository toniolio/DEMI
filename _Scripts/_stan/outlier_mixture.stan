data{
	int num_obs ;
	vector[num_obs] obs ;
	vector[num_obs] stim_speed ;
}
transformed data{
	real obs_sd = sd(obs) ;
	vector[num_obs] scaled_obs = obs/obs_sd ;
}
parameters{
	real good_scale ;
	// real good_shift_intercept ;
	// real good_shift_stim_speed_effect ;
	real<lower=0> good_shift_intercept ;
	real<lower=0> good_shift_stim_speed_effect ;
	real fast_scale_intercept ;
	real fast_scale_stim_speed_effect ;
	real<lower=0,upper=1> p_fast ;
}
transformed parameters{
	vector[num_obs] fast_lpdf;
	vector[num_obs] good_lpdf;
	vector[num_obs] obs_fast_z = (
		(
			scaled_obs
		) ./ (
			(inv_logit(good_scale)/2)
			* inv_logit(
				fast_scale_intercept + fast_scale_stim_speed_effect*stim_speed
			)
		)
	) ;
	vector[num_obs] obs_good_z = (
		(
			scaled_obs
			// - inv_logit( inv_logit(good_shift_intercept) + good_shift_stim_speed_effect*stim_speed )
			- ( (good_shift_intercept) + good_shift_stim_speed_effect*stim_speed )
		) ./ (inv_logit(good_scale)/2)
	) ;
	for (this_obs in 1:num_obs) {
		fast_lpdf[this_obs] = std_normal_lpdf(obs_fast_z[this_obs]) ;
		good_lpdf[this_obs] = std_normal_lpdf(obs_good_z[this_obs]) ;
	}
	fast_lpdf += log(p_fast) ;
	good_lpdf += log1m(p_fast) ;
}
model{
	//implicit uniform(1,3) prior on shape
	target += std_normal_lpdf(good_scale) ;
	// target += std_normal_lpdf(good_shift_intercept) ;
	// target += std_normal_lpdf(good_shift_stim_speed_effect) ;
	target += normal_lpdf(good_shift_intercept|.5,.1) ;
	target += normal_lpdf(good_shift_stim_speed_effect|.25,.05) ;
	target += std_normal_lpdf(fast_scale_intercept) ;
	target += std_normal_lpdf(fast_scale_stim_speed_effect) ;
	target += normal_lpdf(p_fast|0,.2) ; //peaked-at-0, falls to 0 by about .5
	//Likelihood:
	target += log(exp(fast_lpdf)+exp(good_lpdf)) ;

}
generated quantities{
	// real good_mean_0 = obs_sd * inv_logit( inv_logit(good_shift_intercept) ) ;
	// real good_mean_4 = obs_sd * inv_logit( inv_logit(good_shift_intercept) + good_shift_stim_speed_effect*4 ) ;
	real good_mean_0 = obs_sd * ( (good_shift_intercept) ) ;
	real good_mean_2 = obs_sd * ( (good_shift_intercept) + good_shift_stim_speed_effect*2 ) ;
	real good_sd = obs_sd * (inv_logit(good_scale)/2) ;
	real fast_sd_0 = obs_sd * (
		inv_logit(good_scale)
		* inv_logit(
			fast_scale_intercept
		)
	) ;
	real fast_sd_4 = obs_sd * (
		inv_logit(good_scale)
		* inv_logit(
			fast_scale_intercept + fast_scale_stim_speed_effect*4
		)
	) ;
	array[num_obs] int is_fast ;
	vector[num_obs] good_cdf ;
	{
		vector[num_obs] lppr = fast_lpdf - good_lpdf ; // log posterior probability ratio
		for(this_obs in 1:num_obs){
			if(obs_good_z[this_obs]>0){
				is_fast[this_obs] = 0 ;
			}else{
				is_fast[this_obs] = ( ( lppr[this_obs] > 0 ) ? 1 : 0 ) ;
			}
			good_cdf[this_obs] = std_normal_cdf( obs_good_z[this_obs]);
		}
	}
}
