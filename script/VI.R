# Load required libraries
library(rstan)
library(dplyr)
library(ggplot2)
library(bayesplot)
library(gridExtra)


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

