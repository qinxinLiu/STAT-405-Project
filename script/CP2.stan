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
  real log_ka;
  real<lower=0> sigma;
}

transformed parameters {
  real<lower=0> CL_pop = exp(log_CL);
  real<lower=0> V_pop  = exp(log_V);
  real<lower=0> ka_pop = exp(log_ka);
  real<lower=0> ke     = CL_pop / V_pop;
}

model {
  log_CL ~ normal(log(0.2), 0.5);
  log_V  ~ normal(log(3.5), 0.5);
  log_ka ~ normal(log(1), 0.5);
  sigma  ~ normal(0, 0.3);

  for (i in 1:N_obs) {
    real dose = amt[id[i]];
    real pred;

    // Numerical stablization for k_a = k_e.
    // Replace the direct computation with the limiting form of the original formula.
    if (abs(ka_pop - ke) < 1e-6) {
      pred = (dose / V_pop) * ka_pop * time[i] * exp(-ke * time[i]);
    } else {
      pred = (dose * ka_pop / (V_pop * (ka_pop - ke))) *
             (exp(-ke * time[i]) - exp(-ka_pop * time[i]));
    }

    dv[i] ~ lognormal(log(pred + 1e-6), sigma);
  }
}

generated quantities {
  vector[N_obs] dv_sim;
  vector[N_obs] log_lik;

  for (i in 1:N_obs) {
    real dose = amt[id[i]];
    real pred;

    if (abs(ka_pop - ke) < 1e-6) {
      pred = (dose / V_pop) * ka_pop * time[i] * exp(-ke * time[i]);
    } else {
      pred = (dose * ka_pop / (V_pop * (ka_pop - ke))) *
             (exp(-ke * time[i]) - exp(-ka_pop * time[i]));
    }

    dv_sim[i] = lognormal_rng(log(pred + 1e-6), sigma);
    log_lik[i] = lognormal_lpdf(dv[i] | log(pred + 1e-6), sigma);
  }
}
