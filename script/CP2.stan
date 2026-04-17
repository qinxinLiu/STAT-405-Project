data {
  int<lower=1> N_obs;
  int<lower=1> N_subj;
  array[N_obs] int<lower=1, upper=N_subj> id;
  vector<lower=0>[N_obs] time;
  vector<lower=0>[N_obs] dv;
  vector<lower=0>[N_subj] amt;
}

parameters {
  real log_CL;
  real log_V;
  real log_delta;   // delta = ka - ke > 0
  real log_sigma;
}

transformed parameters {
  real<lower=0> CL_pop = exp(log_CL);
  real<lower=0> V_pop  = exp(log_V);
  real<lower=0> ke     = CL_pop / V_pop;
  real<lower=0> delta  = exp(log_delta);
  real<lower=0> ka_pop = ke + delta;
  real<lower=0> sigma  = exp(log_sigma);
}

model {
  // priors on log scale
  log_CL    ~ normal(log(0.2), 0.5);
  log_V     ~ normal(log(3.5), 0.5);
  log_delta ~ normal(log(0.5), 0.7);
  log_sigma ~ normal(log(0.3), 0.5);

  for (i in 1:N_obs) {
    real dose = amt[id[i]];
    real pred;

    pred = dose * ka_pop / (V_pop * delta) *
           (exp(-ke * time[i]) - exp(-ka_pop * time[i]));

    dv[i] ~ lognormal(log(pred + 1e-12), sigma);
  }
}

generated quantities {
  vector[N_obs] dv_sim;

  for (i in 1:N_obs) {
    real dose = amt[id[i]];
    real pred;

    pred = dose * ka_pop / (V_pop * delta) *
           (exp(-ke * time[i]) - exp(-ka_pop * time[i]));

    dv_sim[i] = lognormal_rng(log(pred + 1e-12), sigma);
  }
}
