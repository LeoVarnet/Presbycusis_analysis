// full hierarchical model

data {
  int<lower=0> N;
  int<lower=0> Ntrials;
  int<lower=0> Ncenter;
  int<lower=0> center[N];
  int<lower=0> NCs[N];
  int<lower=0> NCn[N];
  int<lower=0,upper=1> prior_only;
  int<lower=0,upper=1> gender[N];
  vector[N] PTAz;
  vector[N] agez;
}

parameters {
  real<lower=0> beta_0;
  real gammaz_0 [Ncenter];
  real gammaz_cond [Ncenter];
  real gammaz_PTA [Ncenter];
  real gammaz_age [Ncenter];
  real gammaz_agecond [Ncenter];
  real gammaz_agePTA [Ncenter];
  real gammaz_condPTA [Ncenter];
  real gammaz_agecondPTA [Ncenter];
  real gammaz_gender[Ncenter];
  real<lower=0,upper=1> sigma_0;
  real<lower=0,upper=1> sigma_cond;
  real<lower=0,upper=1> sigma_PTA;
  real<lower=0,upper=1> sigma_age;
  real<lower=0,upper=1> sigma_agecond;
  real<lower=0,upper=1> sigma_agePTA;
  real<lower=0,upper=1> sigma_condPTA;
  real<lower=0,upper=1> sigma_agecondPTA;
  real<lower=0,upper=1> sigma_gender;
  real<upper=0> beta_PTA;
  real<upper=0> beta_age;
  real<upper=0> beta_cond;
  real beta_gender;
  real beta_agePTA;
  real beta_condPTA;
  real beta_agecond;
  real beta_agecondPTA;
  real<lower=0,upper=0.5> plapse[Ncenter];
}

transformed parameters {
  vector[N] eta_s;
  vector[N] p_s;
  vector[N] eta_n;
  vector[N] p_n;  
  real gamma_0 [Ncenter];
  real gamma_cond [Ncenter];
  real gamma_PTA [Ncenter];
  real gamma_age [Ncenter];
  real gamma_agecond [Ncenter];
  real gamma_agePTA [Ncenter];
  real gamma_condPTA [Ncenter];
  real gamma_agecondPTA [Ncenter];
  real gamma_gender[Ncenter];

  for (i in 1:Ncenter){
    gamma_0[i] = beta_0 + gammaz_0[i]*sigma_0;
    gamma_cond[i] = beta_cond + gammaz_cond[i]*sigma_cond;
    gamma_PTA[i] = beta_PTA + gammaz_PTA[i]*sigma_PTA;
    gamma_age[i] = beta_age + gammaz_age[i]*sigma_age;
    gamma_agecond[i] = beta_agecond + gammaz_agecond[i]*sigma_agecond;
    gamma_agePTA[i] = beta_agePTA + gammaz_agePTA[i]*sigma_agePTA;
    gamma_condPTA[i] = beta_condPTA + gammaz_condPTA[i]*sigma_condPTA;
    gamma_agecondPTA[i] = beta_agecondPTA + gammaz_agecondPTA[i]*sigma_agecondPTA;
    gamma_gender[i] = beta_gender + gammaz_gender[i]*sigma_gender;
  }
  for (i in 1:N){
    eta_s[i] = gamma_0[center[i]]
          + gamma_age[center[i]] * agez[i]
          + gamma_PTA[center[i]] * PTAz[i]
          + gamma_gender[center[i]] * gender[i]
          + gamma_agePTA[center[i]] * PTAz[i] .* agez[i];
    eta_n[i] = eta_s[i]
          + gamma_cond[center[i]]       
          + gamma_agecond[center[i]] * agez[i]
          + gamma_condPTA[center[i]] * PTAz[i]
          + gamma_agecondPTA[center[i]] * PTAz[i] .* agez[i];
    p_s[i] = 1.0/16.0 + (1-1.0/16.0)*inv_logit(eta_s[i]);
    p_n[i] = 1.0/16.0 + (1-1.0/16.0-plapse[center[i]])*inv_logit(eta_n[i]);
  }
  }

model {
  beta_0 ~ normal(2,1);//normal(0,1);//
  gammaz_0 ~ normal(0,1);
  sigma_0 ~ normal(0,0.05);//normal(0,0.05);
  beta_age ~ normal(0,1);//normal(0,2);
  gammaz_age ~ normal(0,1);//normal(0,2);
  sigma_age ~ normal(0,0.05);//normal(0,2);
  beta_cond ~ normal(0,1);//normal(0,2);
  gammaz_cond ~ normal(0,1);//normal(0,2);
  sigma_cond ~ normal(0,0.05);//normal(0,2);
  beta_gender ~ normal(0,1);//normal(0,2);
  gammaz_gender ~ normal(0,1);//normal(0,2);
  sigma_gender ~ normal(0,0.05);//normal(0,2);
  beta_PTA ~ normal(0,1);//normal(0,2);
  gammaz_PTA ~ normal(0,1);//normal(0,2);
  sigma_PTA ~ normal(0,0.05);//normal(0,2);
  beta_agecond ~ normal(0,1);//normal(0,2);
  gammaz_agecond ~ normal(0,1);//normal(0,2);
  sigma_agecond ~ normal(0,0.05);//normal(0,2);
  beta_condPTA ~ normal(0,1);//normal(0,2);
  gammaz_condPTA ~ normal(0,1);//normal(0,2);
  sigma_condPTA ~ normal(0,0.05);//normal(0,2);
  beta_agecondPTA ~ normal(0,1);//normal(0,2);
  gammaz_agecondPTA ~ normal(0,1);//normal(0,2);
  sigma_agecondPTA ~ normal(0,0.05);//normal(0,2);
  
  plapse ~ beta(3600.0-2756.0,3600.0);
  if(!prior_only)
  for (i in 1:N){
    if (NCn[i]!=3)//(NCn[i]>10)//
      {NCn[i] ~ binomial(Ntrials, p_n[i]);
      NCs[i] ~ binomial(Ntrials, p_s[i]);}
     else
       {0 ~ binomial(Ntrials, p_n[i]);
       NCs[i] ~ binomial(Ntrials, p_s[i]);}
  }
}

generated quantities {
  //vector[N] Rsquared;
  // for (i in 1:N){
    //   if (NCn[i]!=3)//(NCn[i]>10)//
    //   {Rsquared[i] = (eta_s[i]-logit(to_vector(NCs[i])/to_vector(Ntrials))^2}
    
    // Rsquared = (eta_s-logit(to_vector(NCs)/Ntrials)^2.
    
    vector[2*N] log_lik;
    for (i in 1:N){
    //silence
    log_lik[i] = binomial_lpmf(NCs[i] | Ntrials, p_s[i]);
    //noise
    log_lik[N+i] = binomial_lpmf(NCn[i] | Ntrials, p_n[i]);}
    
}
