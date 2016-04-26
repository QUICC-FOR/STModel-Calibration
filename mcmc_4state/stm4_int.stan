functions
{
	vector interval_log(vector x, int[] dataInterval)
	{
		// returns the log probability of the event across the data interval given the 
		// probability on the target interval
		return log(1.0 - ((1.0 - x) ^ (dataInterval/targetInterval)));
	}

	// this is dumb, at least for this model;
	// all of these will look the same
	// should instead adapt it into a single model that applies to all macro params	
	// will need to update the likelihood below
	vector epsilon(e0, interval)
	{ return interval_log(inv_logit(e0), interval); }
	
	beta_b
	beta_t
	theta
	theta_t
	alpha_t
	alpha_b
}
data
{
	int<lower=1> targetInterval;	// the parameters describe transition probabilities for an interval of xx years
//	int<lower=1> n_env_vars;		// number of environmental variables
	
	// data for anystate -> r
	int<lower=1> n_tor;
	int<lower=5,upper=15> interval_tor[n_tor];
//	matrix[n_tor,n_env_vars] vars_tor;

	// data for t -> m
	int<lower=1> n_tm;
	int<lower=5,upper=15> interval_tm[n_tm];
	vector<lower=0,upper=1>[n_tm] prevB_tm;
	vector<lower=0,upper=1>[n_tm] prevM_tm;

	// data for t -> t
	int<lower=1> n_tt;
	int<lower=5,upper=15> interval_tt[n_tt];
	vector<lower=0,upper=1>[n_tt] prevB_tt;
	vector<lower=0,upper=1>[n_tt] prevM_tt;

	// data for b -> m
	int<lower=1> n_bm;
	int<lower=5,upper=15> interval_bm[n_bm];
	vector<lower=0,upper=1>[n_bm] prevT_bm;
	vector<lower=0,upper=1>[n_bm] prevM_bm;

	// data for b -> b
	int<lower=1> n_bb;
	int<lower=5,upper=15> interval_bb[n_bb];
	vector<lower=0,upper=1>[n_bb] prevT_bb;
	vector<lower=0,upper=1>[n_bb] prevM_bb;

	// data for m -> b
	int<lower=1> n_mb;
	int<lower=5,upper=15> interval_mb[n_mb];

	// data for m -> t
	int<lower=1> n_mt;
	int<lower=5,upper=15> interval_mt[n_mt];

	// data for m -> m
	int<lower=1> n_mm;
	int<lower=5,upper=15> interval_mm[n_mm];

	// data for r -> t
	int<lower=1> n_rt;
	int<lower=5,upper=15> interval_rt[n_rt];
	vector<lower=0,upper=1>[n_rt] prevT_rt;
	vector<lower=0,upper=1>[n_rt] prevB_rt;
	vector<lower=0,upper=1>[n_rt] prevM_rt;

	// data for r -> b
	int<lower=1> n_rb;
	int<lower=5,upper=15> interval_rb[n_rb];
	vector<lower=0,upper=1>[n_rb] prevT_rb;
	vector<lower=0,upper=1>[n_rb] prevB_rb;
	vector<lower=0,upper=1>[n_rb] prevM_rb;

	// data for r -> m
	int<lower=1> n_rm;
	int<lower=5,upper=15> interval_rm[n_rm];
	vector<lower=0,upper=1>[n_rm] prevT_rm;
	vector<lower=0,upper=1>[n_rm] prevB_rm;
	vector<lower=0,upper=1>[n_rm] prevM_rm;

	// data for r -> r
	int<lower=1> n_rr;
	int<lower=5,upper=15> interval_rr[n_rr];
	vector<lower=0,upper=1>[n_rr] prevT_rr;
	vector<lower=0,upper=1>[n_rr] prevB_rr;
	vector<lower=0,upper=1>[n_rr] prevM_rr;
}
parameters
{
	real e0;
//	vector[n_env_vars] e;

	real ab0;
	real at0;
	real bb0;
	real bt0;
	real th0;
	real tt0;
}
model
{
	// minimally informative Cauchy priors on the parameters
	e0 ~ cauchy(0,10);
	ab0 ~ cauchy(0,10);
	at0 ~ cauchy(0,10);
	bb0 ~ cauchy(0,10);
	bt0 ~ cauchy(0,10);
	th0 ~ cauchy(0,10);
	tt0 ~ cauchy(0,10);
//	e ~ cauchy(0,2.5);
	
	// an alternative to inv_logit is to use inv_Phi, which is the probit link

	// anystate -> r
//	increment_log_prob(interval_log(inv_logit(e0 + vars_tor*e), interval_tor));
	increment_log_prob(epsilon(e0, interval_tor));

	// t -> m
	increment_log_prob((prevB_tm + prevM.tm) * beta_b(bb0, interval_tm) * 
			(1.0 - epsilon(e0, interval_tm)));
	
	// t -> t
	increment_log_prob(1.0 - ((prevB_tt + prevM.tm) * beta_b(bb0, interval_tt) * 
			(1.0 - epsilon(e0, interval_tt))) - 
			epsilon(e0, interval_tt));

	// b -> m
	increment_log_prob((prevT_bm + prevM.tm) * beta_t(bt0, interval_bm) * 
			(1.0 - epsilon(e0, interval_bm)));

	// b -> b
	increment_log_prob(1.0 - ((prevT_bb + prevM.tm) * beta_t(bt0, interval_bb) * 
			(1.0 - epsilon(e0, interval_bb))) - 
			epsilon(e0, interval_bb));

	// m -> t
	increment_log_prob(theta(th0, interval_mt) * 
			theta_t(tt0, interval_mt) * 
			(1.0 - epsilon(e0, interval_mt)));
			
	// m -> b
	increment_log_prob(theta(th0, interval_mb) * 
			(1.0 - theta_t(tt0, interval_mb)) * 
			(1.0 - epsilon(e0, interval_mb)));
			
	// m -> m
	increment_log_prob(1.0 - 
		(theta(th0, interval_mm) * theta_t(tt0, interval_mm) * (1.0 - epsilon(e0, interval_mm))) - 
		(theta(th0, interval_mm) * (1.0 - theta_t(tt0, interval_mm)) * (1.0 - epsilon(e0, interval_mm))) - 
		epsilon(e0, interval_mm));

	// r -> t
	increment_log_prob(alpha_t(at0, interval_rt) * (prevT_rt + prevM_rt) * (1.0 - alpha_b(ab0, interval_rt) * (prevB_rt + prevM_rt)));
	
	// r -> b
	increment_log_prob(alpha_b(ab0, interval_rb) * (prevB_rb + prevM_rb) * (1.0 - alpha_t(at0, interval_rb) * (prevT_rb + prevM_rb)));

	// r -> m
	increment_log_prob(alpha_b(ab0, interval_rm) * (prevB_rm + prevM_rm) * alpha_t(at0, interval_rm) * (prevT_rm + prevM_rm));

	// r -> r
	increment_log_prob(1.0 - 
		(alpha_t(at0, interval_rr) * (prevT_rr + prevM_rr) * (1.0 - alpha_b(ab0, interval_rr) * (prevB_rr + prevM_rr))) - 
		(alpha_b(ab0, interval_rr) * (prevB_rr + prevM_rr) * (1.0 - alpha_t(at0, interval_rr) * (prevT_rr + prevM_rr))) - 
		(alpha_b(ab0, interval_rr) * (prevB_rr + prevM_rr) * alpha_t(at0, interval_rr) * (prevT_rr + prevM_rr)));

}