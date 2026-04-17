data {
  int<lower=1> N_obs;
  int<lower=1> N_subj;
  array[N_obs] int<lower=1, upper=N_subj> id;
  vector[N_obs] time;
  vector[N_obs] dv;
  vector[N_subj] amt;
}

parameters {
  real<lower=0> CL_pop;
  real<lower=0> V_pop;
  real<lower=0> ka_pop;
  real<lower=0> sigma;
}

model {
  CL_pop ~ lognormal(log(0.2), 1); 
  V_pop ~ lognormal(log(3.5), 1);
  ka_pop ~ lognormal(log(1.0), 1);
  sigma ~ normal(0, 0.5);
  
  for (i in 1:N_obs) {
    real ke = CL_pop / V_pop;
    real dose = amt[id[i]];
    
    real pred = (dose*ka_pop/(V_pop*(ka_pop - ke)))*(exp(-ke*time[i]) - exp(-ka_pop*time[i]));
    
    dv[i] ~ lognormal(log(pred + 1e-6), sigma); 
  }
}

generated quantities {
  vector[N_obs] dv_sim;

  for (i in 1:N_obs) {
    real ke = CL_pop / V_pop;
    real dose = amt[id[i]];
    
    real pred = (dose*ka_pop/(V_pop*(ka_pop - ke)))*(exp(-ke*time[i]) - exp(-ka_pop*time[i]));
    
    dv_sim[i] = lognormal_rng(log(pred + 1e-6), sigma);
  }
}