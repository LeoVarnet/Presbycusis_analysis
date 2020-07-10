data {
int<lower=0> N;
vector[N] age;
vector[N] PTA;
vector[N] group;
}
parameters {
real beta_age;
real beta_group;
real beta_age_group;
real beta0;
real m_age_NH;
real<lower=0> sd_age_NH;
real m_age_HI;
real<lower=0> sd_age_HI;
real<lower=0> sigma;
}
model {
PTA ~ normal(beta0 + beta_age * age + beta_group * group + beta_age_group * age .* group, sigma);
age ~ normal(m_age_HI * group + m_age_NH * (1-group), sd_age_HI * group + sd_age_NH * (1-group));
}