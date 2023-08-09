#' Plots of simulation results estimating mean methylation at two temperatures.
#' See `temp_diff_sims.R` for how they were generated

library("tidyverse")

# Import simulation results
sim_temp_diffs <- read_csv(
  "05_results/14_temperature_diff_sims/sim_temp_diffs.csv",
  col_types = 'ffnnnn'
  )

# Real mean methylation against realised values
# Facet grid shows increasing number of trials
# Means estimate true means well, with substantial spread for few samples
# Colours show error rates - there's no difference
pd <- position_dodge(0.5)
sim_temp_diffs %>%
  mutate(
    lambda = as.factor(lambda),
    real_theta_a = as.factor(real_theta_a),
    real_theta_b = as.factor(real_theta_b)
  ) %>%
  group_by(size, lambda, real_theta_a) %>%
  summarize(
    mean = mean(inf_theta_a),
    lower = quantile(inf_theta_a, 0.02),
    upper = quantile(inf_theta_a, 0.98)
  ) %>%
  ggplot(aes(x = real_theta_a, y = mean, colour = lambda)) +
  geom_point(position = pd) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = pd, width = 0) +
  labs(
    x = "True value",
    y = "Apparent mean"
  ) +
  theme_bw()+
  facet_grid(~size)

#' Mean and CIs of the difference in realised methylation against the true difference.
#' Rows of the facet grid show increasing number of trials
#' Colours show underlying mean at temp A.
sim_temp_diffs %>%
  filter(lambda == 1) %>%
  mutate(
    real_diff = round(real_theta_b - real_theta_a, 3),
    inf_diff  = inf_theta_b - inf_theta_a
  ) %>%
  group_by(size, lambda, real_theta_a, real_diff) %>%
  summarize(
    mean = mean(inf_diff),
    lower = quantile(inf_diff, 0.02),
    upper = quantile(inf_diff, 0.98)
  ) %>%
  ggplot(aes(x = as.factor(real_diff), y = mean, colour = as.factor(real_theta_a))) +
  geom_point(position = pd) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position = pd, width = 0) +
  labs(
    x = "True difference",
    y = "Apparent difference",
    colour = "True mean"
  ) +
  theme_bw()+
  facet_grid(~size)
