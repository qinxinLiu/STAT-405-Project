# Load required libraries
library(rstan) # to remove probably
library(cmdstanr)
library(dplyr)
library(ggplot2)
library(bayesplot)
library(gridExtra)
library(tidyr)


#Set rstan options for performance
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

load(file = "/Users/griffinmcdonald/Documents/GitHub/STAT-405-Project/data/warfarin.rda")
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

obs_data <- obs_data |> left_join(ids, by = "id")
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
model <- stan_model(file = "/Users/griffinmcdonald/Documents/GitHub/STAT-405-Project/script/VI.stan")

#Variational Inference (ADVI)
set.seed(405)
vi_fit <- vb(model, data = stan_data, 
  algorithm = "meanfield",iter = 10000,output_samples = 2000,tol_rel_obj = 0.001)

#Check Output and Params
print(vi_fit, pars = c("CL_pop", "V_pop", "ka_pop", "omega_CL", "omega_V", "omega_ka", "sigma"))

#Extract and Plot ELBO Convergence
vi_summary <- vi_fit@sim$diagnostics

#Extract the posterior samples
posterior <- as.matrix(vi_fit)

#Create individual plots for each parameter
p1 <- mcmc_areas(posterior, pars = "CL_pop", prob = 0.8) + ggtitle("Population Clearance (L/h)")

p2 <- mcmc_areas(posterior, pars = "V_pop", prob =0.8) + ggtitle("Population Volume (L)")

p3 <- mcmc_areas(posterior, pars="ka_pop", prob = 0.8) + ggtitle("Population Absorption Rate (1/h)")

grid.arrange(p1, p2, p3, ncol = 1,top = "Posterior Distributions from VI")



# COMPLETE POOLING MODEL ----
model_CP <- cmdstan_model(stan_file = "script/VI_CP.stan")

# Variational Inference (ADVI)
set.seed(405)
vi_fit_CP = model_CP$variational(
  seed = 1,
  refresh = 500,
  algorithm = "meanfield",
  output_samples = 2000,
  iter = 10000,
  data = stan_data
)

# Extract and Plot ELBO Convergence
# rstan:::get_elbo(vi_fit_CP)
# vi_CP_diag <- vi_fit_CP@sim$diagnostics
# plot(vi_CP_diag$elbo, type = "l",
#      main = "ELBO convergence",
#      xlab = "Iteration",
#      ylab = "ELBO")

# Extract the posterior samples
posterior_CP <- as.matrix(vi_fit_CP)

# Create plots to compare VI posterior for each parameter
posterior_all <- as.data.frame(rbind(posterior_CP[ , 1:3], posterior[ , 1:3]))
posterior_all$model <- c(rep("Complete Pooling", times = 2000), rep("Hierarchical", times = 2000))
posterior_long <- posterior_all %>%
  select(model, CL_pop, V_pop, ka_pop) %>%
  pivot_longer(cols = -model, names_to = "parameter", values_to = "value")

ggplot(posterior_long, aes(x = value, fill = model)) +
  geom_density(alpha = 0.4) +
  facet_wrap(~parameter, scales = "free", ncol = 1) +
  labs(title = "Comparison of Approximate Posterior Distributions",
       x = "Parameter value",
       y = "Density",
       fill = "Model") +
  theme_light()


## Convergence: ELBO ----
# TO DO ####

## Quality of fit: Posterior predictive checks ----
draws_CP <- vi_fit_CP$draws()
y_sim <- posterior::as_draws_matrix(draws_CP[ , 7:485], variable = "^dv_sim")
y_obs <- stan_data$dv

bayesplot::ppc_dens_overlay(y_obs, y_sim[1:200, ])


## Predictive performance: LOOCV ----
