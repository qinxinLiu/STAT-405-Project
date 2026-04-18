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
init_fun <- function() {
  list(
    log_CL = rnorm(1, log(0.2), 0.5),
    log_V  = rnorm(1, log(3.5), 0.5),
    log_ka = rnorm(1, log(1.0), 0.5),
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
                         init = init_fun, 
                     iter_warmup = 2000, 
                     iter_sampling = 8000, 
                     show_messages = FALSE)


## Filter the unconverged posterior parameters.
hmc_fit_cp$summary() |>
  filter(rhat > 1.05) |>
  print(n = Inf)

hmc_pars_cp <- c("CL_pop", "V_pop", "ka_pop", "sigma")
hmc_fit_cp$summary(variables = hmc_pars_cp)

## Trace plot ----

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
    y = "Concentration"
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

## Define initialization function
init_fun <- function() {
  list(
    log_CL_pop = rnorm(1, log(0.2), 0.05),
    log_V_pop  = rnorm(1, log(3.5), 0.05),
    log_ka_pop = rnorm(1, log(1.0), 0.05),
    sigma    = runif(1, 0.08, 0.3),
    omega_CL = runif(1, 0.05, 0.3),
    omega_V  = runif(1, 0.05, 0.3),
    eta_CL = rnorm(stan_data$N_subj, 0, 0.5),
    eta_V  = rnorm(stan_data$N_subj, 0, 0.5)
  )
}

# Trace plot of add random effect to all pk parameters

## Create a stan object for Full Random-Effect Hierarchical model
mod_h_notwork <- cmdstan_model("VI.stan")
hmc_fit_h_notwork <- mod_h_notwork$sample(
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
hmc_fit_h_notwork$summary() |>
  filter(rhat > 1.05) |>
  print(n = Inf)


## Create a stan object for Hierarchical model
mod_h <- cmdstan_model("Hierarchical.stan")
init_fun <- function() {
  list(
    log_CL_pop = rnorm(1, log(0.2), 0.05),
    log_V_pop  = rnorm(1, log(3.5), 0.05),
    log_ka_pop = rnorm(1, log(1.0), 0.05),
    sigma    = runif(1, 0.08, 0.3),
    omega_CL = runif(1, 0.05, 0.3),
    omega_V  = runif(1, 0.05, 0.3),
    eta_CL = rnorm(stan_data$N_subj, 0, 0.5),
    eta_V  = rnorm(stan_data$N_subj, 0, 0.5)
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

pop_pars <- c(
  "CL_pop",
  "V_pop",
  "ka_pop",
  "sigma",
  "omega_CL",
  "omega_V"
)

# Randomly drawn subject-specific parameters
subj_pars <- c(
  "eta_CL[1]", #male
  "eta_CL[14]", #female
  "eta_V[1]",
  "eta_V[14]" 
)


hmc_fit_h$summary(variables = pop_pars)
hmc_fit_h$summary(variables = subj_pars)

## Trace plot ---- 
mcmc_trace(
  hmc_fit_h$draws(variables = pop_pars),
  facet_args = list(ncol = 3)
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



# MODEL EXTENSION ----

## Gender ----
sex <- obs_data |>
  group_by(new_id) |>
  slice(1) |>
  ungroup() |>
  pull(sex)
sex <- factor(sex, levels = c("male", "female"))
gender<- as.integer(sex) - 1



## Data for model extension
stan_data_e <- list(
    N_obs = nrow(obs_data),
    N_subj = length(unique_ids),
    id = obs_data$new_id,
    time = obs_data$time,
    dv = obs_data$dv,
    amt = dosing_data$amt,
    gender = gender
  )

## Model Fitting ----
## Create a stan object for model extension
mod_e <- cmdstan_model("ME.stan")

## Define initialization function
init_fun<- function() {
  list(
    log_CL_pop = rnorm(1, log(0.2), 0.5),
    log_V_pop  = rnorm(1, log(3.5), 0.5),
    log_ka_pop = rnorm(1, log(1.0), 0.5),
    
    sigma    = runif(1, 0.08, 0.3),
    
    beta_CL  = rnorm(1, 0, 0.5),
    beta_V   = rnorm(1, 0, 0.5),
    
    omega_CL = runif(1, 0.05, 0.3),
    omega_V  = runif(1, 0.05, 0.3),
    
    eta_CL   = rnorm(stan_data$N_subj, 0, 0.5),
    eta_V    = rnorm(stan_data$N_subj, 0, 0.5)
  )
}


hmc_fit_e <- mod_e$sample(
  data = stan_data_e,
  init = init_fun,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 2000,
  iter_sampling = 8000,
  adapt_delta = 0.99,
  max_treedepth = 12
)

## Filter the unconverged posterior parameters.
hmc_fit_e$summary() |>
  filter(rhat > 1.05) |>
  print(n = Inf)


pop_pars_e <- c(
  "log_CL_pop",
  "log_V_pop",
  "log_ka_pop",
  "sigma",
  "beta_CL",
  "beta_V",
  "omega_CL",
  "omega_V"
)

subj_pars_e <- c(
  "eta_CL[1]",
  "eta_CL[14]",
  "eta_V[1]",
  "eta_V[14]"
)

hmc_fit_e$summary(variables = pop_pars_e)
hmc_fit_e$summary(variables = subj_pars)

## Trace plot ---- 
mcmc_trace(
  hmc_fit_e$draws(variables = pop_pars_e),
  facet_args = list(ncol = 2)
)

mcmc_trace(
  hmc_fit_e$draws(variables = subj_pars_e),
  facet_args = list(ncol = 2)
)

## Gender effect
gender_pars_e <- c("beta_CL", "beta_V")

mcmc_trace(
  hmc_fit_e$draws(variables = gender_pars_e),
  facet_args = list(ncol = 2)
)

## Posterior predictive check ----

## Original posterior predictive plot
y_obs <- stan_data_e$dv
yrep <- as_draws_matrix(hmc_fit_e$draws("dv_sim"))
yrep <- as.matrix(yrep)

ppc_dens_overlay(y = y_obs, yrep = yrep[1:50, ])


## Posterior concentration vs. time plot
hmc_sim_mat_e <- hmc_fit_e$draws(variables = "dv_sim", format = "matrix")

hmc_ppc_df_e <- tibble(
  obs = seq_len(ncol(hmc_sim_mat_e)),
  id = stan_data_e$id,
  time = stan_data_e$time,
  dv_obs = stan_data_e$dv,
  gender_num = stan_data_e$gender[stan_data_e$id],
  pred_med = apply(hmc_sim_mat_e, 2, median),
  pred_lwr = apply(hmc_sim_mat_e, 2, quantile, probs = 0.05),
  pred_upr = apply(hmc_sim_mat_e, 2, quantile, probs = 0.95)
) |>
  mutate(
    gender = factor(gender_num, levels = c(0, 1), labels = c("Male", "Female"))
  ) |>
  arrange(id, time)

## Orignal plot
ggplot(hmc_ppc_df_e, aes(x = time)) +
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


## Plot colored by gender
ggplot(hmc_ppc_df_e, aes(x = time, color = gender, fill = gender)) +
  geom_ribbon(aes(ymin = pred_lwr, ymax = pred_upr), alpha = 0.2, colour = NA) +
  geom_line(aes(y = pred_med), linewidth = 0.7) +
  geom_point(aes(y = dv_obs), size = 1.1) +
  facet_wrap(~ id, scales = "free_y") +
  labs(
    x = "Time",
    y = "Concentration",
    title = "Posterior predictive check by subject and gender"
  ) +
  theme_bw()


hmc_ppc_df_e2 <- hmc_ppc_df_e |>
  group_by(gender, id) |>
  mutate(id_index = cur_group_id()) |>
  ungroup()

ggplot(hmc_ppc_df_e2, aes(x = time, group = id, color = id_index)) +
  geom_line(aes(y = pred_med), linewidth = 0.8, alpha = 0.9) +
  geom_point(aes(y = dv_obs), size = 1, alpha = 0.7) +
  facet_wrap(~ gender, ncol = 1, scales = "free_y") +
  scale_color_viridis_c() +
  labs(
    x = "Time",
    y = "Concentration",
    color = "Subject index",
    title = "Posterior predictive concentration-time profiles by gender"
  ) +
  theme_bw()




## gender-level summary band
gender_band <- hmc_ppc_df_e |>
  group_by(gender, time) |>
  summarise(
    band_med = median(pred_med),
    band_lwr = quantile(pred_med, 0.05),
    band_upr = quantile(pred_med, 0.95),
    .groups = "drop"
  )

## Single plot with gender
ggplot() +
  geom_ribbon(
    data = gender_band,
    aes(x = time, ymin = band_lwr, ymax = band_upr, fill = gender),
    alpha = 0.2
  ) +
  geom_line(
    data = hmc_ppc_df_e,
    aes(x = time, y = pred_med, group = id, color = gender),
    alpha = 0.25,
    linewidth = 0.5
  ) +
  geom_line(
    data = gender_band,
    aes(x = time, y = band_med, color = gender),
    linewidth = 1.2
  ) +
  labs(
    x = "Time",
    y = "Concentration",
    color = "Gender",
    fill = "Gender") +
  theme_bw()


## LOO-CV check ----
hmc_log_lik_mat_e <- as_draws_matrix(hmc_fit_e$draws("log_lik"))
dim(hmc_log_lik_mat_e)

## Larger elpd_loo and smaller looic are better
hmc_loo_res_e <- loo(hmc_log_lik_mat_e) 
print(hmc_loo_res_e)

## A good model should have more k concentrated below 0.5
pareto_k_table(hmc_loo_res_e)
plot(hmc_loo_res_e)

## Identify the bad points
hmc_k_vals_e <- pareto_k_values(hmc_loo_res_e)
hmc_bad_idx_e <- which(hmc_k_vals_e > 0.7)

## The bad points
hmc_bad_points_e <- data.frame(
  obs  = hmc_bad_idx_e,
  id   = stan_data_e$id[hmc_bad_idx_e],
  time = stan_data_e$time[hmc_bad_idx_e],
  dv   = stan_data_e$dv[hmc_bad_idx_e],
  k    = hmc_k_vals_e[hmc_bad_idx_e]
)


# HMC MODEL COMPARISON ----
