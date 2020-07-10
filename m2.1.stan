// model with PTA and cond centered hierarchical

data {
  int<lower=0> N;
  int<lower=0> Ntrials;
  int<lower=0> Ncenter;
  int<lower=0> center[N];
  int<lower=0> NCs[N];
  int<lower=0> NCn[N];
  int<lower=0,upper=1> prior_only;
  vector[N] PTAz;
  vector[N] agez;
}

parameters {
  real<lower=0> beta_0;
  real gammaz_0 [Ncenter];
  real<lower=0,upper=1> sigma_0;
  real<upper=0> beta_PTA;
  real gammaz_PTA[Ncenter];//real<upper=0> gamma_PTA[Ncenter];
  real<lower=0> sigma_PTA;
  real<upper=0> beta_age;//real beta_age;//
  real gammaz_age[Ncenter];//real<upper=0> gamma_age[Ncenter];//
  real<lower=0> sigma_age;
  real<upper=0> beta_cond;
  real gammaz_cond[Ncenter];
  real<lower=0> sigma_cond;
  real beta_agePTA;
  real gammaz_agePTA[Ncenter];
  real<lower=0> sigma_agePTA;
  real beta_condPTA;
  real gammaz_condPTA[Ncenter];
  real<lower=0> sigma_condPTA;
  real beta_agecond;
  real gammaz_agecond[Ncenter];
  real<lower=0> sigma_agecond;
  real beta_agecondPTA;
  real gammaz_agecondPTA[Ncenter];
  real<lower=0> sigma_agecondPTA;
  real<lower=0,upper=0.5> plapse[Ncenter];
  // real U_plapse;
  // real V_plapse;
}

transformed parameters {
  vector[N] eta_s;
  vector[N] p_s;
  vector[N] eta_n;
  vector[N] p_n;
  real alpha_plapse;
  real beta_plapse;
  real gamma_0[Ncenter];
  real gamma_PTA[Ncenter];
  real gamma_age[Ncenter];
  real gamma_cond[Ncenter];
  real gamma_agePTA[Ncenter];
  real gamma_condPTA[Ncenter];
  real gamma_agecond[Ncenter];
  real gamma_agecondPTA[Ncenter];
    
  for (i in 1:Ncenter){
    gamma_0[i] = beta_0 + gammaz_0[i]*sigma_0;
    gamma_PTA[i] = beta_PTA + gammaz_PTA[i]*sigma_PTA;
    gamma_age[i] = beta_age + gammaz_age[i]*sigma_age;
    gamma_cond[i] = beta_cond + gammaz_cond[i]*sigma_cond;
    gamma_agePTA[i] = beta_agePTA + gammaz_agePTA[i]*sigma_agePTA;
    gamma_condPTA[i] = beta_condPTA + gammaz_condPTA[i]*sigma_condPTA;
    gamma_agecond[i] = beta_agecond + gammaz_agecond[i]*sigma_agecond;
    gamma_agecondPTA[i] = beta_agecondPTA + gammaz_agecondPTA[i]*sigma_agecondPTA;
  }
  
  //real beta_PTA = 0;
  //real beta_age = 0;
  //real beta_agePTA = 0;
  //real beta_condPTA = 0;
  
  for (i in 1:N){
    eta_s[i] = gamma_0[center[i]]
          + gamma_age[center[i]] * agez[i]
          + gamma_PTA[center[i]] * PTAz[i]
          + gamma_agePTA[center[i]] * PTAz[i] .* agez[i];
    eta_n[i] = gamma_0[center[i]] + gamma_cond[center[i]]
          + (gamma_age[center[i]] + gamma_agecond[center[i]]) * agez[i]
          + (gamma_PTA[center[i]] + gamma_condPTA[center[i]]) * PTAz[i]
          + (gamma_agePTA[center[i]] + gamma_agecondPTA[center[i]]) * PTAz[i] .* agez[i];
    //masking[i] = gamma_cond[center[i]]+ gamma_agecond[center[i]] * agez[i] + gamma_agecondPTA[center[i]] * PTAz[i] .* agez[i];
    p_s[i] = 1.0/16.0 + (1-1.0/16.0)*inv_logit(eta_s[i]);
    p_n[i] = 1.0/16.0 + (1-1.0/16.0-plapse[center[i]])*inv_logit(eta_n[i]);
  }
  
  beta_plapse = 3600.0;//V_plapse - U_plapse*V_plapse;
  alpha_plapse = beta_plapse-2756.0;//U_plapse*V_plapse;
  
  
  }

model {
  beta_0 ~ normal(2,1);//normal(0,2);
  gammaz_0 ~ normal(0,1);
  sigma_0 ~ normal(0,0.05);//exponential(0.00001);//lognormal(-8,1);////lognormal(-3,1);//
  beta_age ~ normal(0,1);//normal(0,2);
  gammaz_age ~ normal(0,1);
  sigma_age ~ normal(0,0.05);//exponential(1.5);////lognormal(-3,1);//
  beta_cond ~ normal(0,1);//normal(0,2);
  gammaz_cond ~ normal(0,1);
  sigma_cond ~ normal(0,0.05);//exponential(1.5);////lognormal(-3,1);//
  beta_PTA ~ normal(0,1);//normal(0,2);
  gammaz_PTA ~ normal(0,1);
  sigma_PTA ~ normal(0,0.05);//exponential(1.5);////lognormal(-3,1);//
  beta_agePTA ~ normal(0,1);//normal(0,2);
  gammaz_agePTA ~ normal(0,1);
  sigma_agePTA ~ normal(0,0.05);//exponential(1.5);////lognormal(-3,1);//
  beta_condPTA ~ normal(0,1);//normal(0,2);
  gammaz_condPTA ~ normal(0,1);
  sigma_condPTA ~ normal(0,0.05);//exponential(1.5);////lognormal(-3,1);//
  beta_agecond ~ normal(0,1);//normal(0,2);
  gammaz_agecond ~ normal(0,1);
  sigma_agecond ~ normal(0,0.05);//exponential(1.5);////lognormal(-3,1);//
  beta_agecondPTA ~ normal(0,1);//normal(0,2);
  gammaz_agecondPTA ~ normal(0,1);
  sigma_agecondPTA ~ normal(0,0.05);//exponential(1.5);////lognormal(-3,1);//
  //plapse ~ beta(2,5);//normal(0,1);//normal(0,2);
  plapse ~ beta(alpha_plapse,beta_plapse);//normal(0,1);//normal(0,2);
  // U_plapse ~ uniform(0,1);
  // V_plapse ~ gamma(1,20);
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
