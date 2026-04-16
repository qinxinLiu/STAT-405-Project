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

  real<lower=0> omega_CL;
  real<lower=0> omega_V;
  real<lower=0> omega_ka;

  vector[N_subj] eta_CL;
  vector[N_subj] eta_V;
  vector[N_subj] eta_ka;
}

transformed parameters {
  vector[N_subj] CL;
  vector[N_subj] V;
  vector[N_subj] ka;

  for (i in 1:N_subj) {
    CL[i] = CL_pop * exp(omega_CL * eta_CL[i]);
    V[i]  = V_pop * exp(omega_V * eta_V[i]);
    ka[i] = ka_pop * exp(omega_ka * eta_ka[i]);
  }
}

model {
  CL_pop ~ lognormal(log(0.1), 1); 
  V_pop ~ lognormal(log(8.0), 1);
  ka_pop ~ lognormal(log(1.0), 1);
  sigma ~ normal(0, 0.5);

  omega_CL ~ normal(0, 0.5);
  omega_V ~ normal(0, 0.5);
  omega_ka ~ normal(0, 0.5);

  eta_CL ~ std_normal();
  eta_V ~ std_normal();
  eta_ka ~ std_normal();

  for (i in 1:N_obs) {
    real ke = CL[id[i]] / V[id[i]];
    real dose = amt[id[i]];
    real ka_subj = ka[id[i]];
    real V_subj = V[id[i]];
    
    real pred = (dose*ka_subj/(V_subj*(ka_subj - ke)))*(exp(-ke*time[i]) - exp(-ka_subj*time[i]));
    
    dv[i] ~ lognormal(log(pred + 1e-6), sigma); 
  }
}