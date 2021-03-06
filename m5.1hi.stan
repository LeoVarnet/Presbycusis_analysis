// full varying-intercepts hierarchical model

data {
  int<lower=0> N;
  int<lower=0> Ntrials;
  int<lower=0> Ncenter;
  int<lower=0> center[N];
  int<lower=0> NCs[N];
  int<lower=0> NCn[N];
  int<lower=0> Nmissing;
  int<lower=0,upper=1> prior_only;
  int<lower=0,upper=1> gender[N];
  int<lower=0,upper=1> group[N];
  vector[N] ESII_sz;
  vector[N] ESII_nz;
  vector[N] agez;
}

parameters {
  real<lower=0> beta_0;
  real gammaz_0 [Ncenter];
  real<lower=0,upper=1> sigma_0;
  real<lower=0> beta_ESII;
  real<upper=0> beta_age;
  real beta_cond;
  real beta_gender;
  real beta_group;
  real beta_groupcond;
  real beta_groupage;
  //real beta_groupESII;
  real beta_ageESII;
  real beta_condESII;
  real beta_agecond;
  //real beta_groupageESII;
  //real beta_groupcondESII;
  //real beta_groupagecond;
  real beta_agecondESII;
  //real beta_groupagecondESII;
  real<lower=0,upper=0.5> plapse [2];
}

transformed parameters {
  vector[N] eta_s;
  vector[N] p_s;
  vector[N] eta_n;
  vector[N] p_n;  
  real gamma_0[Ncenter];

  for (i in 1:Ncenter){
    gamma_0[i] = beta_0 + gammaz_0[i]*sigma_0;
  }
  for (i in 1:N){
    eta_s[i] = gamma_0[center[i]]
          + beta_age * agez[i]
          + beta_ESII * ESII_sz[i]
          + beta_gender * gender[i]
          + beta_ageESII * ESII_sz[i] .* agez[i]
          + beta_group * group[i]
          + beta_groupage * group[i] .* agez[i];
          //+ beta_groupESII * group[i] .* ESII_sz[i];
          // + beta_groupageESII * group[i] * ESII_sz[i] .* agez[i];
    eta_n[i] = gamma_0[center[i]]
          + beta_age * agez[i]
          + beta_ESII * ESII_nz[i]
          + beta_gender * gender[i]
          + beta_ageESII * ESII_nz[i] .* agez[i]
          + beta_group * group[i]
          + beta_groupage * group[i] .* agez[i]
          //+ beta_groupESII * group[i] .* ESII_nz[i]
          // + beta_groupageESII * group[i] * ESII_nz[i] .* agez[i]
          + beta_cond
          + beta_agecond * agez[i]
          + beta_condESII * ESII_nz[i]
          + beta_agecondESII * ESII_nz[i] .* agez[i]
          + beta_groupcond * group[i];
          //+ beta_groupagecond * group[i] .* agez[i];
          // + beta_groupcondESII * group[i] .* ESII_nz[i]
          // + beta_groupagecondESII * group[i] * ESII_nz[i] .* agez[i];
    p_s[i] = 1.0/16.0 + (1-1.0/16.0)*inv_logit(eta_s[i]);
    p_n[i] = 1.0/16.0 + (1-1.0/16.0-plapse[group[i]+1])*inv_logit(eta_n[i]);
  }
  }

model {
  beta_0 ~ normal(0,1);//normal(2,1);//normal(0,1);//
  gammaz_0 ~ normal(0,1);
  sigma_0 ~ normal(0,0.05);//normal(0,0.05);
  beta_age ~ normal(0,1);//normal(0,2);
  beta_cond ~ normal(0,1);//normal(0,2);
  beta_gender ~ normal(0,1);//normal(0,2);
  beta_ESII ~ normal(0,1);//normal(0,2);
  beta_group ~ normal(0,1);
  beta_groupcond ~ normal(0,1);
  beta_groupage ~ normal(0,1);
  //beta_groupESII ~ normal(0,1);
  beta_agecond ~ normal(0,1);//normal(0,2);
  beta_condESII ~ normal(0,1);//normal(0,2);
  beta_ageESII ~ normal(0,1);//normal(0,2);
  //beta_groupageESII ~ normal(0,1);
  //beta_groupcondESII ~ normal(0,1);
  //beta_groupagecond ~ normal(0,1);
  beta_agecondESII ~ normal(0,1);//normal(0,2);
  //beta_groupagecondESII ~ normal(0,1);
  
  plapse ~ beta(1,1);//beta(3600.0-2756.0,3600.0);
  if(!prior_only){
  for (i in 1:N){
    if (NCn[i]!=100)
      {NCn[i] ~ binomial(Ntrials, p_n[i]);
      NCs[i] ~ binomial(Ntrials, p_s[i]);}
     else
       {0 ~ binomial(Ntrials, p_n[i]);
       NCs[i] ~ binomial(Ntrials, p_s[i]);}
  }
  }
}

generated quantities {
    vector[2*N-Nmissing] log_lik;
    vector[N] err_s;
  vector[N] err_n;
  
    int j=1;
  //residuals
  err_s = logit(p_s)-logit(to_vector(NCs)/Ntrials);
  err_n = logit(p_n)-logit(to_vector(NCn)/Ntrials);
    for (i in 1:N){
    log_lik[i] = binomial_lpmf(NCs[i] | Ntrials, p_s[i]);    //silence
    if (NCn[i]!=100)
    {log_lik[N+j] = binomial_lpmf(NCn[i] | Ntrials, p_n[i]); //noise
    j=j+1;}
    }
}
