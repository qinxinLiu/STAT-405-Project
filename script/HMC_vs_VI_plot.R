library(bayesplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(posterior)

## ---------- helper functions ----------

extract_common_pop_draws <- function(fit, cp_vi = FALSE) {
  dr <- posterior::as_draws_df(fit$draws())
  
  if (cp_vi) {
    out <- dr %>%
      transmute(
        CL_pop = exp(log_CL),
        V_pop  = exp(log_V),
        ka_pop = exp(log_ka),
        sigma  = sigma
      )
  } else {
    out <- dr %>%
      transmute(
        CL_pop = CL_pop,
        V_pop  = V_pop,
        ka_pop = ka_pop,
        sigma  = sigma
      )
  }
  
  posterior::as_draws_df(out)
}

build_interval_df <- function(fit, method, model, cp_vi = FALSE) {
  extract_common_pop_draws(fit, cp_vi = cp_vi) %>%
    bayesplot::mcmc_intervals_data() %>%
    mutate(
      method = method,
      model  = model
    )
}

## ---------- build combined summary data ----------

int_hmc_cp <- build_interval_df(hmc_fit_cp, "HMC", "Complete Pooling")
int_hmc_h  <- build_interval_df(hmc_fit_h,  "HMC", "Hierarchical")
int_hmc_e  <- build_interval_df(hmc_fit_e,  "HMC", "Extended")

int_vi_cp  <- build_interval_df(vi_fit_CP,  "VI",  "Complete Pooling", cp_vi = TRUE)
int_vi_h   <- build_interval_df(vi_fit,     "VI",  "Hierarchical")
int_vi_e   <- build_interval_df(vi_fit_ext, "VI",  "Extended")

interval_df <- bind_rows(
  int_hmc_cp, int_hmc_h, int_hmc_e,
  int_vi_cp,  int_vi_h,  int_vi_e
) %>%
  mutate(
    parameter = factor(
      parameter,
      levels = c("CL_pop", "V_pop", "ka_pop", "sigma")
    ),
    model = factor(
      model,
      levels = c("Complete Pooling", "Hierarchical", "Extended")
    ),
    method = factor(
      method,
      levels = c("HMC", "VI")
    )
  )

pd <- position_dodge(width = 0.5)

## ---------- 2x2 combined plot ----------

ggplot(
  interval_df,
  aes(
    y = model,
    x = m,
    color = model,
    shape = method,
    linetype = method
  )
) +
  geom_errorbarh(
    aes(xmin = ll, xmax = hh),
    position = pd,
    height = 0.12,
    linewidth = 0.6,
    alpha = 0.4
  ) +
  geom_errorbarh(
    aes(xmin = l, xmax = h),
    position = pd,
    height = 0.12,
    linewidth = 1.0
  ) +
  geom_point(
    position = pd,
    size = 2.3
  ) +
  facet_wrap(~ parameter, ncol = 2, scales = "free_x") +
  labs(
    x = "Posterior estimate",
    y = NULL,
    color = "Model",
    shape = "Method",
    linetype = "Method"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.background = element_blank()
  )