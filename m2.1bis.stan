// simple model with only main effects

data {
  int<lower=0> N;
  int<lower=0> Ntrials;
  //int<lower=0> Ncenter;
  //int<lower=0> center[N];
  int<lower=0> NCs[N];
  int<lower=0> NCn[N];
  int<lower=0,upper=1> prior_only;
  int<lower=0, upper=1> gender[N];
  vector[N] PTAz;
  vector[N] agez;
}

parameters {
  real<lower=0> beta_0;
  real<upper=0> beta_PTA;
  real<upper=0> beta_age;
  real<upper=0> beta_cond;
  real beta_agecond;
  real beta_condPTA;
  real beta_agePTA;
  real beta_agecondPTA;
  real beta_gender;
  real<lower=0,upper=0.5> plapse;
}

transformed parameters {
  vector[N] eta_s;
  vector[N] p_s;
  vector[N] eta_n;
  vector[N] p_n;
    
  for (i in 1:N){
    eta_s[i] = beta_0
          + beta_age * agez[i]
          + beta_PTA * PTAz[i]
          + beta_gender * gender[i];
    eta_n[i] = eta_s[i]
          + beta_cond          
          + beta_agecond * agez[i]
          + beta_condPTA * PTAz[i]
          + beta_agecondPTA * PTAz[i] .* agez[i];
    //masking[i] = gamma_cond[center[i]]+ gamma_agecond[center[i]] * agez[i] + gamma_agecondPTA[center[i]] * PTAz[i] .* agez[i];
    p_s[i] = 1.0/16.0 + (1-1.0/16.0)*inv_logit(eta_s[i]);
    p_n[i] = 1.0/16.0 + (1-1.0/16.0-plapse)*inv_logit(eta_n[i]);
  }

  }

model {
  beta_0 ~ normal(2,1);//normal(0,2);
  beta_age ~ normal(0,1);//normal(0,2);
  beta_cond ~ normal(0,1);//normal(0,2);
  beta_gender ~ normal(0,1);//normal(0,2);
  beta_PTA ~ normal(0,10);//normal(0,2);
  beta_agecond ~ normal(0,1);//normal(0,2);
  beta_condPTA ~ normal(0,1);//normal(0,2);
  beta_agePTA ~ normal(0,1);//normal(0,2);
  beta_agecondPTA ~ normal(0,1);//normal(0,2);
  plapse ~ beta(1,1);
  if(!prior_only)
  for (i in 1:N){
    if (NCn[i]!=100)//(NCn[i]>10)//
      {NCn[i] ~ binomial(Ntrials, p_n[i]);
      NCs[i] ~ binomial(Ntrials, p_s[i]);}
     else
       {//1 ~ binomial(Ntrials, p_n[i]);
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
