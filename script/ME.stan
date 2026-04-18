data {
  int<lower=1> N_obs;
  int<lower=1> N_subj;
  array[N_obs] int<lower=1, upper=N_subj> id;
  vector<lower=0>[N_obs] time;
  vector<lower=0>[N_obs] dv;
  vector<lower=0>[N_subj] amt;

  array[N_subj] int<lower=0, upper=1> gender;
}

parameters {
  // population parameters on log scale
  real log_CL_pop;
  real log_V_pop;
  real log_ka_pop;

  real<lower=0> sigma;

  // gender fixed effects
  real beta_CL;
  real beta_V;

  // random effects
  real<lower=0> omega_CL;
  real<lower=0> omega_V;

  vector[N_subj] eta_CL;
  vector[N_subj] eta_V;
}

transformed parameters {
    // back-transform population parameters
  real<lower=0> CL_pop = exp(log_CL_pop);
  real<lower=0> V_pop  = exp(log_V_pop);
  real<lower=0> ka_pop = exp(log_ka_pop);

  vector<lower=0>[N_subj] CL;
  vector<lower=0>[N_subj] V;

  for (i in 1:N_subj) {
    CL[i] = exp(log_CL_pop + beta_CL * gender[i] + omega_CL * eta_CL[i]);
    V[i]  = exp(log_V_pop  + beta_V  * gender[i] + omega_V  * eta_V[i]);
  }
}


model {
  // priors
  log_CL_pop ~ normal(log(0.2), 0.5);
  log_V_pop  ~ normal(log(3.5), 0.5);
  log_ka_pop ~ normal(log(1), 0.5);

  sigma    ~ normal(0, 0.3);
  omega_CL ~ normal(0, 0.3);
  omega_V  ~ normal(0, 0.3);

  beta_CL ~ normal(0, 0.5);
  beta_V  ~ normal(0, 0.5);

  eta_CL ~ std_normal();
  eta_V  ~ std_normal();

  for (i in 1:N_obs) {
    real ke = CL[id[i]] / V[id[i]];
    real dose = amt[id[i]];
    real ka_subj = ka_pop;
    real V_subj = V[id[i]];
    real pred;

    if (abs(ka_subj - ke) < 1e-6) {
      pred = (dose / V_subj) * ka_subj * time[i] * exp(-ke * time[i]);
    } else {
      pred = (dose * ka_subj / (V_subj * (ka_subj - ke))) *
             (exp(-ke * time[i]) - exp(-ka_subj * time[i]));
    }

    dv[i] ~ lognormal(log(pred + 1e-6), sigma);
  }
}

generated quantities {
  vector[N_obs] dv_sim;
  vector[N_obs] log_lik;

  // Odds ratio of gender 
  real CL_ratio = exp(beta_CL);
  real V_ratio  = exp(beta_V);

  for (i in 1:N_obs) {
    real ke = CL[id[i]] / V[id[i]];
    real dose = amt[id[i]];
    real ka_subj = ka_pop;
    real V_subj = V[id[i]];
    real pred;

    if (abs(ka_subj - ke) < 1e-6) {
      pred = (dose / V_subj) * ka_subj * time[i] * exp(-ke * time[i]);
    } else {
      pred = (dose * ka_subj / (V_subj * (ka_subj - ke))) *
             (exp(-ke * time[i]) - exp(-ka_subj * time[i]));
    }

    dv_sim[i] = lognormal_rng(log(pred + 1e-6), sigma);
    log_lik[i] = lognormal_lpdf(dv[i] | log(pred + 1e-6), sigma);
  }
}
