// simple varying-intercepts hierarchical model with only main effects 

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
  real<lower=0,upper=1> sigma_0[Ncenter];
  real<upper=0> beta_PTA;
  //real gammaz_PTA[Ncenter];//real<upper=0> gamma_PTA[Ncenter];
  //real<lower=0> sigma_PTA;
  real<upper=0> beta_age;//real beta_age;//
  //real gammaz_age[Ncenter];//real<upper=0> gamma_age[Ncenter];//
  //real<lower=0> sigma_age;
  real<upper=0> beta_cond;
  //real gammaz_cond[Ncenter];
  //real<lower=0> sigma_cond;
  real beta_gender;
  //real gammaz_gender[Ncenter];
  //real<lower=0> sigma_gender;
  real<lower=0,upper=0.5> plapse[Ncenter];
}

transformed parameters {
  vector[N] eta_s;
  vector[N] p_s;
  vector[N] eta_n;
  vector[N] p_n;  
  real gamma_0[Ncenter];
  // real gamma_PTA[Ncenter];
  // real gamma_age[Ncenter];
  // real gamma_cond[Ncenter];
  // real gamma_gender[Ncenter];
    
  for (i in 1:Ncenter){
    gamma_0[i] = beta_0 + gammaz_0[i]*sigma_0[i];
    // gamma_PTA[i] = beta_PTA + gammaz_PTA[i]*sigma_PTA;
    // gamma_age[i] = beta_age + gammaz_age[i]*sigma_age;
    // gamma_cond[i] = beta_cond + gammaz_cond[i]*sigma_cond;
    // gamma_gender[i] = beta_gender + gammaz_gender[i]*sigma_gender;
  }
  for (i in 1:N){
    eta_s[i] = gamma_0[center[i]]
          + beta_age * agez[i]
          + beta_PTA * PTAz[i]
          + beta_gender;
    eta_n[i] = gamma_0[center[i]] 
          + beta_cond
          + beta_age * agez[i]
          + beta_PTA * PTAz[i]
          + beta_gender;
    p_s[i] = 1.0/16.0 + (1-1.0/16.0)*inv_logit(eta_s[i]);
    p_n[i] = 1.0/16.0 + (1-1.0/16.0-plapse[center[i]])*inv_logit(eta_n[i]);
  }

  }

model {
  beta_0 ~ normal(2,1);//normal(0,2);
  gammaz_0 ~ normal(0,1);
  sigma_0 ~ normal(0,0.05);
  beta_age ~ normal(0,1);//normal(0,2);
  // gammaz_age ~ normal(0,1);
  // sigma_age ~ normal(0,0.05);
  beta_cond ~ normal(0,1);//normal(0,2);
  // gammaz_cond ~ normal(0,1);
  // sigma_cond ~ normal(0,0.05);
  beta_gender ~ normal(0,1);//normal(0,2);
  // gammaz_gender ~ normal(0,1);
  // sigma_gender ~ normal(0,0.05);
  beta_PTA ~ normal(0,1);//normal(0,2);
  // gammaz_PTA ~ normal(0,1);
  // sigma_PTA ~ normal(0,0.05);
  
  plapse ~ beta(3600.0-2756.0,3600.0);
  if(!prior_only)
  for (i in 1:N){
    if (NCn[i]!=3)//(NCn[i]>10)//
      {NCn[i] ~ binomial(Ntrials, p_n[i]);
      NCs[i] ~ binomial(Ntrials, p_s[i]);}
     else
       {1 ~ binomial(Ntrials, p_n[i]);
       NCs[i] ~ binomial(Ntrials, p_s[i]);}
  }
}
