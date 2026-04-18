# Load required libraries
library(cmdstanr)
library(dplyr)
library(ggplot2)
library(bayesplot)
library(gridExtra)
library(tidyr)
library(posterior)
library(loo)

setwd("/Users/jessieliu/Desktop/Msc Statistics/2025 Winter Term 2/STAT 405/STAT-405-Project")

# DATA SETUP ----

# Load the dataset
load(file = "/Users/jessieliu/Desktop/Msc Statistics/2025 Winter Term 2/STAT 405/STAT-405-Project/warfarin.rda")
df <- warfarin

#Extract single dose amount per subject
dosing_data <- df |>
  filter(evid == 1) |>
  group_by(id) |>
  summarize(amt = max(amt), .groups = 'drop')

#Filter for observation data
obs_data <- df |> filter(evid == 0, dv > 0) |> arrange(id, time)

unique_ids <- unique(obs_data$id)
ids <- data.frame(id = unique_ids, new_id = 1:length(unique_ids))

obs_data <- obs_data |> left_join(ids, by = "id") |> filter(dvid == "cp")
dosing_data <- dosing_data |> left_join(ids, by = "id") |> arrange(new_id)

stan_data <- list(
  N_obs = nrow(obs_data),
  N_subj = length(unique_ids),
  id = obs_data$new_id,
  time = obs_data$time,
  dv = obs_data$dv,
  amt = dosing_data$amt
)


# COMPLETE POOLING ----


## Define initial value
init <- function() {
  list(
    log_CL = rnorm(1, log(0.2), 0.1),
    log_V  = rnorm(1, log(3.5), 0.1),
    log_ka = rnorm(1, log(1.0), 0.1),
    sigma  = runif(1, 0.1, 0.5)
  )
}

## Model Fitting ----

# Create a new CmdStan model object
mod_cp <- cmdstan_model("CP2.stan")

# Fit model
n_chains <- 4
hmc_fit_cp <- mod$sample(data = stan_data, 
                         chains = n_chains, 
                         init = init, 
                     iter_warmup = 2000, 
                     iter_sampling = 8000, 
                     show_messages = FALSE)


## Filter the unconverged posterior parameters.
hmc_fit_cp$summary() |>
  filter(rhat > 1.05) |>
  print(n = Inf)

## Trace plot ----
hmc_pars_cp <- c("log_CL", "log_V", "log_ka", "sigma")

mcmc_trace(
  hmc_fit_cp$draws(variables = hmc_pars_cp),
  facet_args = list(ncol = 2)
)


## Posterior predictive check ----

## Original posterior predictive plot
y_obs <- stan_data$dv
yrep <- as_draws_matrix(hmc_fit_cp$draws("dv_sim"))
yrep <- as.matrix(yrep)
ppc_dens_overlay(y = y_obs, yrep = yrep[1:50, ])


## Posterior concentration vs. time plot
hmc_sim_mat_cp <- hmc_fit_cp$draws(variables = "dv_sim", format = "matrix")

hmc_ppc_df_cp <- tibble(
  obs = seq_len(ncol(hmc_sim_mat_cp)),
  id = stan_data$id,
  time = stan_data$time,
  dv_obs = stan_data$dv,
  pred_med = apply(hmc_sim_mat_cp, 2, median),
  pred_lwr = apply(hmc_sim_mat_cp, 2, quantile, probs = 0.05),
  pred_upr = apply(hmc_sim_mat_cp, 2, quantile, probs = 0.95)
) |> arrange(id, time)


ggplot(hmc_ppc_df_cp, aes(x = time)) +
  geom_ribbon(aes(ymin = pred_lwr, ymax = pred_upr), alpha = 0.25) +
  geom_line(aes(y = pred_med), linewidth = 0.7) +
  geom_point(aes(y = dv_obs), size = 1.2) +
  facet_wrap(~ id, scales = "free_y") +
  labs(
    x = "Time",
    y = "Concentration",
    title = "Posterior predictive check: concentration vs time"
  ) +
  theme_bw()

## LOO-CV check ----
hmc_log_lik_mat_cp <- as_draws_matrix(hmc_fit_cp$draws("log_lik"))
dim(hmc_log_lik_mat_cp)

## Larger elpd_loo and smaller looic are better
hmc_loo_res_cp <- loo(hmc_log_lik_mat_cp) 
print(hmc_loo_res_cp)

## A good model should have more k concentrated below 0.5
pareto_k_table(hmc_loo_res_cp)
plot(hmc_loo_res_cp)

## Identify the bad points
hmc_k_vals_cp <- pareto_k_values(hmc_loo_res_cp)
hmc_bad_idx_cp <- which(hmc_k_vals_cp > 0.7)

## The bad points
hmc_bad_points_cp <- data.frame(
  obs  = hmc_bad_idx_cp,
  id   = stan_data$id[hmc_bad_idx_cp],
  time = stan_data$time[hmc_bad_idx_cp],
  dv   = stan_data$dv[hmc_bad_idx_cp],
  k    = hmc_k_vals_h[hmc_bad_idx_cp]
)





# HIERARCHICAL MODEL ----

## Model Fitting ----
## Create a stan object for Hierarchical model
mod_h <- cmdstan_model("Hierarchical.stan")

## Define initialization function
init_fun <- function() {
  list(
    log_CL_pop = rnorm(1, log(0.2), 0.05),
    log_V_pop  = rnorm(1, log(3.5), 0.05),
    log_ka_pop = rnorm(1, log(1.0), 0.05),
    sigma    = runif(1, 0.08, 0.2),
    omega_CL = runif(1, 0.05, 0.15),
    omega_V  = runif(1, 0.05, 0.15),
    omega_ka = runif(1, 0.05, 0.15),
    eta_CL = rnorm(stan_data$N_subj, 0, 0.02),
    eta_V  = rnorm(stan_data$N_subj, 0, 0.02),
    eta_ka = rnorm(stan_data$N_subj, 0, 0.02)
  )
}


hmc_fit_h <- mod_h$sample(
  data = stan_data,
  init = init_fun,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 8000,
  adapt_delta = 0.99,
  max_treedepth = 12
)

## Filter the unconverged posterior parameters.
hmc_fit_h$summary() |>
  filter(rhat > 1.05) |>
  print(n = Inf)

## Trace plot ---- 
pop_pars <- c(
  "lp__", #log posterior density
  "log_CL_pop",
  "log_V_pop",
  "log_ka_pop",
  "omega_CL",
  "omega_V"
)

subj_pars <- c(
  "eta_CL[1]",
  "eta_CL[15]",
  "eta_CL[29]",
  "eta_V[1]",
  "eta_V[15]",
  "eta_V[29]"
)

mcmc_trace(
  hmc_fit_h$draws(variables = pop_pars),
  facet_args = list(ncol = 2)
)

mcmc_trace(
  hmc_fit_h$draws(variables = subj_pars),
  facet_args = list(ncol = 2)
)

## Posterior predictive check ----

## Original posterior predictive plot
y_obs <- stan_data$dv
yrep <- as_draws_matrix(hmc_fit_h$draws("dv_sim"))
yrep <- as.matrix(yrep)
ppc_dens_overlay(y = y_obs, yrep = yrep[1:50, ])


## Posterior concentration vs. time plot
hmc_sim_mat_h <- hmc_fit_h$draws(variables = "dv_sim", format = "matrix")

hmc_ppc_df_h <- tibble(
  obs = seq_len(ncol(hmc_sim_mat_h)),
  id = stan_data$id,
  time = stan_data$time,
  dv_obs = stan_data$dv,
  pred_med = apply(hmc_sim_mat_h, 2, median),
  pred_lwr = apply(hmc_sim_mat_h, 2, quantile, probs = 0.05),
  pred_upr = apply(hmc_sim_mat_h, 2, quantile, probs = 0.95)
) |> arrange(id, time)


ggplot(hmc_ppc_df_h, aes(x = time)) +
  geom_ribbon(aes(ymin = pred_lwr, ymax = pred_upr), alpha = 0.25) +
  geom_line(aes(y = pred_med), linewidth = 0.7) +
  geom_point(aes(y = dv_obs), size = 1.2) +
  facet_wrap(~ id, scales = "free_y") +
  labs(
    x = "Time",
    y = "Concentration",
    title = "Posterior predictive check: concentration vs time"
  ) +
  theme_bw()

## LOO-CV check ----
hmc_log_lik_mat_h <- as_draws_matrix(hmc_fit_h$draws("log_lik"))
dim(hmc_log_lik_mat_h)

## Larger elpd_loo and smaller looic are better
hmc_loo_res_h <- loo(hmc_log_lik_mat_h) 
print(hmc_loo_res_h)

## A good model should have more k concentrated below 0.5
pareto_k_table(hmc_loo_res_h)
plot(hmc_loo_res_h)

## Identify the bad points
hmc_k_vals_h <- pareto_k_values(hmc_loo_res_h)
hmc_bad_idx_h <- which(hmc_k_vals_h > 0.7)

## The bad points
hmc_bad_points_h <- data.frame(
  obs  = hmc_bad_idx_h,
  id   = stan_data$id[hmc_bad_idx_h],
  time = stan_data$time[hmc_bad_idx_h],
  dv   = stan_data$dv[hmc_bad_idx_h],
  k    = hmc_k_vals_h[hmc_bad_idx_h]
)



# MODEL EXTENSION

# HMC MODEL COMPARISON 
