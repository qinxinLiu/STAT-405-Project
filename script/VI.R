# Load required libraries
library(cmdstanr)
library(dplyr)
library(ggplot2)
library(bayesplot)
library(gridExtra)
library(tidyr)


#Set rstan options for performance
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

load(file = "data/warfarin.rda")
df <- warfarin
#colnames(warfarin)
#Preprocessing

#Extract single dose amount per subject
dosing_data <- df |>
  filter(evid == 1) |>
  group_by(id) |>
  summarize(amt = max(amt), .groups = 'drop')

#Filter for observation data
obs_data <- df |> filter(evid == 0, dv > 0) |> arrange(id, time)

unique_ids <- unique(obs_data$id)
ids <- data.frame(id = unique_ids, new_id = 1:length(unique_ids))

obs_data <- obs_data |> left_join(ids, by = "id") |> filter(dvid == "cp") # remove pca values
dosing_data <- dosing_data |> left_join(ids, by = "id") |> arrange(new_id)

stan_data <- list(
  N_obs = nrow(obs_data),
  N_subj = length(unique_ids),
  id = obs_data$new_id,
  time = obs_data$time,
  dv = obs_data$dv,
  amt = dosing_data$amt
)

# HIERARCHICAL MODEL ----
model <- cmdstan_model(stan_file = "script/Hierarchical.stan")

#Variational Inference (ADVI)
set.seed(405)
vi_fit <- model$variational(
  seed = 1,
  refresh = 500,
  algorithm = "meanfield",
  output_samples = 2000,
  iter = 10000,
  data = stan_data,
  eval_elbo = 50,
  save_latent_dynamics = TRUE
)

# Extract the posterior samples
draws <- vi_fit$draws()
posterior <- posterior::as_draws_matrix(draws)

#Create individual plots for each parameter
p1 <- mcmc_areas(posterior, pars = "CL_pop", prob = 0.8) + ggtitle("Population Clearance (L/h)")

p2 <- mcmc_areas(posterior, pars = "V_pop", prob =0.8) + ggtitle("Population Volume (L)")

p3 <- mcmc_areas(posterior, pars="ka_pop", prob = 0.8) + ggtitle("Population Absorption Rate (1/h)")

grid.arrange(p1, p2, p3, ncol = 1,top = "Posterior Distributions from VI")


## Convergence: ELBO ----
# median elbo converges at 350 iterations according to model fit
vi_fit$latent_dynamics_files() # to get values
elbo <- as.data.frame(cbind(iter = c(seq(from = 50, to = 1100, by = 50)), 
                            ELBO = c(-977.78498
                                     -610.55066,
                                     -584.79822,
                                     -532.7852,
                                     -540.29097,
                                     -520.37681,
                                     -530.3298,
                                     -516.87348,
                                     -515.33438,
                                     -517.60476,
                                     -516.0683,
                                     -515.55344,
                                     -514.47138,
                                     -527.63497,
                                     -516.85296,
                                     -515.93813,
                                     -513.27726,
                                     -515.1719,
                                     -513.61062,
                                     -515.5289,
                                     -523.15042,
                                     -517.09038,
                                     -513.75776
                            )))

ggplot(data = elbo, aes(x = iter, y = ELBO)) +
  geom_line() +
  labs(x = "Number of iterations") +
  theme_light()


## Quality of fit: Posterior predictive checks ----
draws <- vi_fit$draws()
y_sim <- posterior::as_draws_matrix(draws[ , 7:253], variable = "^dv_sim")
y_obs <- stan_data$dv

bayesplot::ppc_dens_overlay(y_obs, y_sim[1:200, ])

# 90% predictive intervals
preds <- cbind(
  Estimate = colMeans(y_sim), 
  Q5 = apply(y_sim, 2, quantile, probs = 0.05),
  Q95 = apply(y_sim, 2, quantile, probs = 0.95)
)

ggplot(cbind(obs_data, preds), aes(x = time, y = Estimate)) +
  geom_smooth(aes(ymin = Q5, ymax = Q95), stat = "identity", linewidth = 0.5) +
  geom_point(aes(y = dv)) + 
  labs(
    y = "Concentration (mg/L)", 
    x = "Time (h)",
    title = "90% predictive intervals for Hierarchical"
  ) +
  facet_wrap(~ id, nrow = 5) +
  theme_light()



# COMPLETE POOLING MODEL ----
model_CP <- cmdstan_model(stan_file = "script/CP2.stan")

# Variational Inference (ADVI)
set.seed(405)
vi_fit_CP <- model_CP$variational(
  seed = 1,
  refresh = 500,
  algorithm = "meanfield",
  output_samples = 2000,
  iter = 10000,
  data = stan_data,
  eval_elbo = 50,
  save_latent_dynamics = TRUE
)

# Extract the posterior samples
draws_CP <- vi_fit_CP$draws()
posterior_CP <- posterior::as_draws_matrix(draws_CP)

# Create plots to compare VI posteriors for each parameter
posterior_all <- as.data.frame(rbind(posterior_CP[ , 3:5], posterior[ , 3:5]))
posterior_all$model <- c(rep("Complete Pooling", times = 2000), rep("Hierarchical", times = 2000))
posterior_long <- posterior_all %>%
  select(model, log_CL, log_V, log_ka) %>%
  pivot_longer(cols = -model, names_to = "parameter", values_to = "log_value")

ggplot(posterior_long, aes(x = exp(log_value), fill = model)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~parameter, scales = "free", ncol = 1) +
  labs(title = "Comparison of Approximate Posterior Distributions",
       x = "Parameter value",
       y = "Density",
       fill = "Model") +
  theme_light()


## Convergence: ELBO ----
# median elbo converges at 350 iterations according to model fit
vi_fit_CP$latent_dynamics_files() # to get values
elbo <- as.data.frame(cbind(iter = c(seq(from = 50, to = 550, by = 50)), 
             ELBO = c(-886.40897,
                      -739.84712,
                      -633.44807,
                      -559.17576,
                      -549.8991,
                      -547.9312,
                      -550.27883,
                      -549.76396,
                      -549.22639,
                      -547.50176,
                      -547.73766
             )))

ggplot(data = elbo, aes(x = iter, y = ELBO)) +
  geom_line() +
  labs(x = "Number of iterations") +
  theme_light()
  

## Quality of fit: Posterior predictive checks ----
draws_CP <- vi_fit_CP$draws()
y_sim <- posterior::as_draws_matrix(draws_CP[ , 7:253], variable = "^dv_sim")
y_obs <- stan_data$dv

bayesplot::ppc_dens_overlay(y_obs, y_sim[1:200, ])

# 90% predictive intervals
preds <- cbind(
  Estimate = colMeans(y_sim), 
  Q5 = apply(y_sim, 2, quantile, probs = 0.05),
  Q95 = apply(y_sim, 2, quantile, probs = 0.95)
)

ggplot(cbind(obs_data, preds), aes(x = time, y = Estimate)) +
  geom_smooth(aes(ymin = Q5, ymax = Q95), stat = "identity", linewidth = 0.5) +
  geom_point(aes(y = dv)) + 
  labs(
    y = "Concentration (mg/L)", 
    x = "Time (h)",
    title = "90% predictive intervals for CP"
  ) +
  facet_wrap(~ id, nrow = 5) +
  theme_light()


## Predictive performance: LOOCV ----


# MODEL EXTENSION ----
model_ext <- cmdstan_model(stan_file = "script/ME.stan")

# Gender
sex <- obs_data |>
  group_by(new_id) |>
  slice(1) |>
  ungroup() |>
  pull(sex)
sex <- factor(sex, levels = c("male", "female"))
gender<- as.integer(sex) - 1

# Data for model extension
stan_data_ext <- list(
  N_obs = nrow(obs_data),
  N_subj = length(unique_ids),
  id = obs_data$new_id,
  time = obs_data$time,
  dv = obs_data$dv,
  amt = dosing_data$amt,
  gender = gender
)

#Variational Inference (ADVI)
set.seed(405)
vi_fit_ext <- model_ext$variational(
  seed = 1,
  refresh = 500,
  algorithm = "meanfield",
  output_samples = 2000,
  iter = 10000,
  data = stan_data_ext,
  eval_elbo = 50,
  save_latent_dynamics = TRUE
)

# Extract the posterior samples
draws_ext <- vi_fit_ext$draws()
posterior_ext <- posterior::as_draws_matrix(draws_ext)


## Convergence: ELBO ----
# median elbo converges at 350 iterations according to model fit
vi_fit_ext$latent_dynamics_files() # to get values
elbo <- as.data.frame(cbind(iter = c(seq(from = 50, to = 1100, by = 50)), 
                            ELBO = c(-869.97425,
                                     -563.66845,
                                     -533.40084,
                                     -552.16894,
                                     -514.04264,
                                     -512.29236,
                                     -514.83782,
                                     -511.75054,
                                     -517.52108,
                                     -512.9051,
                                     -521.50707,
                                     -535.66194,
                                     -511.18953,
                                     -511.7233,
                                     -515.43701,
                                     -515.71281,
                                     -514.34929,
                                     -508.44704,
                                     -511.83623,
                                     -509.47202,
                                     -520.86266,
                                     -515.85726
                                     
                            )))

ggplot(data = elbo, aes(x = iter, y = ELBO)) +
  geom_line() +
  labs(x = "Number of iterations") +
  theme_light()


## Quality of fit: Posterior predictive checks
draws_ext <- vi_fit_ext$draws()
y_sim <- posterior::as_draws_matrix(draws_ext[ , 7:253], variable = "^dv_sim")
y_obs <- stan_data$dv

bayesplot::ppc_dens_overlay(y_obs, y_sim[1:200, ])

VI_ppc_df_e <- tibble(
  obs = seq_len(ncol(VI_gend)),
  id = stan_data_ext$id,
  time = stan_data_ext$time,
  dv_obs = stan_data_ext$dv,
  
  gender_num = stan_data_ext$gender[stan_data_ext$id],
  
  pred_med = apply(VI_gend, 2, median),
  pred_lwr = apply(VI_gend, 2, quantile, probs = 0.05),
  pred_upr = apply(VI_gend, 2, quantile, probs = 0.95)
) |>
  mutate(
    gender = factor(gender_num, levels = c(0, 1), labels = c("Male", "Female"))
  ) |>
  arrange(id, time)

gender_band <- VI_ppc_df_e |>
  group_by(time, gender) |>
  summarize(
    band_med = median(pred_med),
    band_lwr = median(pred_lwr), 
    band_upr = median(pred_upr),
    .groups = "drop"
  )

ggplot() +
  geom_ribbon(
    data = gender_band,
    aes(x = time, ymin = band_lwr, ymax = band_upr, fill = gender),
    alpha = 0.2
  ) +
  geom_line(
    data = VI_ppc_df_e,
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
    x = "Time (h)",
    y = "Concentration (mg/L)",
    color = "Gender",
    fill = "Gender",
    title = "VI Posterior Predictive Concentration Over Time by Gender"
  ) +
  theme_bw() +
  scale_fill_manual(values = c("Male" = "#F8766D", "Female" = "#00BFC4")) +
  scale_color_manual(values = c("Male" = "#F8766D", "Female" = "#00BFC4"))

# 90% predictive intervals
preds <- cbind(
  Estimate = colMeans(y_sim), 
  Q5 = apply(y_sim, 2, quantile, probs = 0.05),
  Q95 = apply(y_sim, 2, quantile, probs = 0.95)
)

ggplot(cbind(obs_data, preds), aes(x = time, y = Estimate)) +
  geom_smooth(aes(ymin = Q5, ymax = Q95), stat = "identity", linewidth = 0.5) +
  geom_point(aes(y = dv)) + 
  labs(
    y = "Concentration (mg/L)", 
    x = "Time (h)",
    title = "90% predictive intervals for CP"
  ) +
  facet_wrap(~ id, nrow = 5) +
  theme_light()


# Extract summary for the Complete Pooling model
cp_summary <- vi_fit_CP$summary(
  variables = c("sigma")
)
print(cp_summary[, c("variable", "mean", "sd", "q5", "q95")])

# Extract summary for the Hierarchical model
hier_summary <- vi_fit$summary(
  variables = c("sigma")
)
print(hier_summary[, c("variable", "mean", "sd", "q5", "q95")])


# Extract summary for the Extended model
ext_summary <- vi_fit_ext$summary(
  variables = c("beta_CL","beta_V")
)

print(ext_summary[, c("variable", "mean", "sd", "q5", "q95")])

time_cp <- vi_fit_CP$time()$total
time_hier <- vi_fit$time()$total
time_ext <- vi_fit_ext$time()$total
runtimes <- data.frame(
  Model = c("Complete Pooling", "Hierarchical", "Extended"),
  Runtime_Seconds = c(time_cp, time_hier, time_ext)
)
print(runtimes)


library(loo)
#Complete Pooling Model LOO
log_lik_cp <- vi_fit_CP$draws("log_lik", format = "matrix")
# Calculate LOO-CV
loo_cp <- loo(log_lik_cp)
pareto_k_table(loo_cp)
plot(loo_cp)
# Hierarchical Model LOO ---
log_lik_hier <- vi_fit$draws("log_lik", format = "matrix")
loo_hier <- loo(log_lik_hier)
pareto_k_table(loo_hier)
plot(loo_hier)
# Extended Model LOO ---
log_lik_ext <- vi_fit_ext$draws("log_lik", format = "matrix")
loo_ext <- loo(log_lik_ext)
pareto_k_table(loo_ext)
plot(loo_ext)

