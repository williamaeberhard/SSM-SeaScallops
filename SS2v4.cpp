// Alternative SSM for SPA4 sea scallops fishery, Yin et al. (2018) CJS
// SSModel with mods (SS2v4):
// Summary:
//   * 3 process eq: biomass, recruits and natural mortality
//   * 4 obs eq: clappers, survey commercial, survey recruits, catch

#include <TMB.hpp>

template <class Type> 
Type sqr(Type x){
	return x*x;
}

template<class Type>
Type objective_function<Type>::operator() () {
	
	//----------------------------------------------------------------------------
	// Data
	//----------------------------------------------------------------------------

	DATA_VECTOR(I); // response: survey index commercial size (dim NY)
	DATA_VECTOR(IR); // response: survey index recruit size (dim NY)
	DATA_VECTOR(L); // response: survey estimate clappers (dim NY)
	DATA_VECTOR(C); // response: catch index (dim NY)

	DATA_VECTOR(g); // covariate: growth rate commercial size (dim NY)
	DATA_VECTOR(gR); // covariate: growth rate recruit size (dim NY)
	DATA_VECTOR(N); // covariate: survey estimate of live animals (dim NY)

	//----------------------------------------------------------------------------
	// Parameters
	//----------------------------------------------------------------------------

	PARAMETER(log_sigma_tau); // proc sd biomass
	PARAMETER(log_sigma_phi); // proc sd recruits
	PARAMETER(log_sigma_m); // proc sd nat mort
	PARAMETER(log_sigma_epsilon); // obs sd survey indices commercial size
	PARAMETER(log_sigma_upsilon); // obs sd survey indices recruits
	PARAMETER(log_sigma_kappa); // obs sd clappers
	PARAMETER(log_sigma_C); // obs sd effort dynamic for catches
	PARAMETER(log_q_I); // catchability for commercial size
	PARAMETER(log_q_R); // catchability for recruits
	PARAMETER(log_S); // dissolution rate in clappers eq
	PARAMETER(log_a); // in effort dynamic catch eq
	PARAMETER(log_chi); // in effort dynamic catch eq

	PARAMETER_VECTOR(log_B); // randeff: biomass (dim NY)
	PARAMETER_VECTOR(log_R); // randeff: recruits (dim NY)
	PARAMETER_VECTOR(log_m); // randeff: nat mort (dim NY)

	//----------------------------------------------------------------------------
	// Setup
	//----------------------------------------------------------------------------

	int NY = I.size(); // number of years

	Type sigma_tau = exp(log_sigma_tau);
	Type sigma_phi = exp(log_sigma_phi);
	Type sigma_m = exp(log_sigma_m);
	Type sigma_epsilon = exp(log_sigma_epsilon);
	Type sigma_upsilon = exp(log_sigma_upsilon);
	Type sigma_kappa = exp(log_sigma_kappa);
	Type sigma_C = exp(log_sigma_C);

	Type q_I = exp(log_q_I);
	Type q_R = exp(log_q_R);
	Type S = exp(log_S);
	Type a = exp(log_a);
	Type chi = exp(log_chi);
	
	vector<Type> B = exp(log_B);
	vector<Type> R = exp(log_R);
	vector<Type> m = exp(log_m);

	vector<Type> nll_comp(7); // nll components, break down contributions
	nll_comp.setZero(); // initialize 

	//----------------------------------------------------------------------------
	// Proc eq
	//----------------------------------------------------------------------------

	// RW for recruits
	for (int t=1; t<NY; t++){ // R[0] free, yet predicted
	  	Type mean_proc_R = R[t-1];
	  	nll_comp[0] -= dnorm(log_R[t], log(mean_proc_R)-sqr(sigma_phi)/2.0, sigma_phi, true);
	}
	
	// RW for nat mort
	for (int t=1; t<NY; t++){ // m[0] free, yet predicted
	  	Type mean_proc_m = m[t-1];
	  	nll_comp[1] -= dnorm(log_m[t], log(mean_proc_m)-sqr(sigma_m)/2.0, sigma_m, true);
	}

	// biomass
	for (int t=1; t<NY; t++){ // B[0] free, yet predicted
	  	Type mean_proc_B = exp(-m[t])*g[t-1]*(B[t-1]-C[t-1])+exp(-m[t])*gR[t-1]*R[t-1];
	  	nll_comp[2] -= dnorm(log_B[t], log(mean_proc_B)-sqr(sigma_tau)/2.0, sigma_tau, true);
	}
	
	//----------------------------------------------------------------------------
	// Obs eq
	//----------------------------------------------------------------------------

	// survey indices for commercial size
	for (int t=0; t<NY; t++){
		Type mean_obs_I = q_I*B[t];
  	nll_comp[3] -= dnorm(log(I[t]), log(mean_obs_I)-sqr(sigma_epsilon)/2.0, sigma_epsilon, true);
	}

	// survey indices for recruits
	for (int t=0; t<NY; t++){
		Type mean_obs_R = q_R*R[t];
  	nll_comp[4] -= dnorm(log(IR[t]), log(mean_obs_R)-sqr(sigma_upsilon)/2.0, sigma_upsilon, true);
	}

	// clappers
	Type mean_obs_L = m[0]*S*N[0]; // initial one different
	nll_comp[5] -= dnorm(log(L[0]), log(mean_obs_L)-sqr(sigma_kappa)/2.0, sigma_kappa, true);
	for (int t=1; t<NY; t++){
		Type mean_obs_L = m[t]*S*(S*N[t-1] + (2.0-S)*N[t])/2.0;
		nll_comp[5] -= dnorm(log(L[t]), log(mean_obs_L)-sqr(sigma_kappa)/2.0, sigma_kappa, true);
	}
	
	// effort dynamic eq for catch, nothing for C[0]
	Type denom = a*B[0]/2.0;
	for (int t=1; t<NY; t++){
		Type mean_obs_C = C[t-1]/B[t-1]*B[t]*pow(B[t-1]/denom,chi);
		nll_comp[6] -= dnorm(log(C[t]), log(mean_obs_C)-sqr(sigma_C)/2.0, sigma_C, true);
	}


	//----------------------------------------------------------------------------
	// Outputs
	//----------------------------------------------------------------------------
	
	REPORT(nll_comp);

	ADREPORT(sigma_tau);
	ADREPORT(sigma_phi);
	ADREPORT(sigma_m);
	ADREPORT(sigma_epsilon);
	ADREPORT(sigma_upsilon);
	ADREPORT(sigma_kappa);
	ADREPORT(sigma_C);
	ADREPORT(q_I);
	ADREPORT(q_R);
	ADREPORT(S);
	ADREPORT(a);
	ADREPORT(chi);
	
	ADREPORT(B);
	ADREPORT(R);
	ADREPORT(m);
	
	Type nll = nll_comp.sum();
	return nll;
}
