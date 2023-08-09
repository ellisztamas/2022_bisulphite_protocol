#' How many TEs do you need to sequence to get a reliable estimate of mean
#' methylation?
#'
#' Error rates on a 1000bp window are roughly beta distritbuted with a=17, b=220
#' with mean ~0.072
#' The average TE in TAIR10 is 798bp, which is roughly the same scale.
#' I don't expect to get a reliable estimate for any single TE, but if we average
#' over sufficient TEs we can.
#' This simulation aims to check how many TEs you would need to look at and at
#' what coverage to say something sensible about mean methylation, or a difference
#' in mean methylation between two temperatures.
#'
#' To do this:
#' - Simulate a ton of 1000bp TEs (=200 CHH positions)
#' - Simulate an error rate for each from rbeta(., 17,220). Also use 0 for comparison
#' - Draw methylated/methylated reads from a binomial with errors for individual
#'    TEs given its error rate
#' - Sum reads over two or more TEs
#' - Estimate mean methylation using the *mean* error rate
#'
#' Averaging over 100 loci should do the trick, even for pretty low coverage.


library('tidyverse')
library('ggpubr')

source("02_library/R/binomial_with_errors.R")
source("02_library/R/sim_mean_mC_over_loci.R")

# Using CMT2-targetted TEs as a guide, about 70% of cytosines are in the CHH
# context.
# Assuming cytosines are 1/4 of nucleotides, I should simulate 200 CHH sites
cmt2 <- read_csv(file = "05_results/13_identifying_te_methylation/output/mC_on_cmt2_TEs.csv")
cmt2 %>%
  filter(context != "total") %>%
  group_by(id) %>%
  reframe(
    context,
    ncytosines,
    pcytosines = ncytosines / sum(ncytosines)
  ) %>%
  ggplot(aes(x = ncytosines, y = pcytosines)) +
  geom_point() +
  facet_grid(~ context)

# Simulate data with no conversion errors
sim_loci_zero <- vector('list', length = length(coverage))
for(cov in 1:length(coverage)){
  sim_loci_zero[[cov]] <- sim_mean_mC_over_loci(
    shape1=0,
    shape2=220,
    theta = 0.15,
    ncytosines = 200,
    coverage = coverage[cov],
    ndraws = 200,
    mean_lambda1 = 0
  )
}
sim_loci_zero <- do.call('rbind', sim_loci_zero)
sim_loci_zero$data <- "No errors"

# Simulate mean lambda1=0.072
# Errors drawn from Beta(17,220)
coverage <- c(1, 2, 4, 8, 16, 32)
sim_loci_errors <- vector('list', length = length(coverage))
for(cov in 1:length(coverage)){
  sim_loci_errors[[cov]] <- sim_mean_mC_over_loci(
    shape1=17,
    shape2=220,
    theta = 0.15,
    ncytosines = 200,
    coverage = coverage[cov],
    ndraws = 200,
    mean_lambda1 = 0.072
  )
}
sim_loci_errors <- do.call('rbind', sim_loci_errors)
sim_loci_errors$data <- "With errors"
# Bind into one dataframe
sim_loci <- rbind(sim_loci_errors, sim_loci_zero)

write_csv(sim_loci, "05_results/14_temperature_diff_sims/average_over_TEs.csv")

#' Boxplot shows distribution of theta_hat for simulations with/without errors
#' for increasing number of loci
#' Facets show increasing coverage.
#' N. loci has a larger effect than coverage or errors
sim_loci %>%
  ggplot(aes(x = as.factor(loci_to_average), y = theta_hat, groups = data, fill=data)) +
  geom_boxplot() +
  labs(
    x = "Number of loci to average",
    y = "Estimated mean methylation"
  ) +
  theme_bw() +
  facet_grid(~coverage)

# The width of 96% confidence intervals
sim_loci %>%
  group_by(data, coverage, loci_to_average) %>%
  # summarise(
  #   diff = sd(theta_hat)
  # ) %>%
  summarise(
    lower = quantile(theta_hat, 0.02),
    upper = quantile(theta_hat, 0.98),
    diff  = upper - lower
  ) %>%
  ggplot(aes( x= coverage, y = diff, groups=data, colour=data))+
  geom_line() +
  labs(
    x = "Coverage",
    y = "Width of 96% confidence intervals"
  ) +
  # scale_y_continuous(trans = 'log10') +
  # scale_x_continuous(trans = 'log10') +
  facet_grid(~loci_to_average)

# Averaging over coveragem this corresponds to roughly 1.5-fold increase in CIs.
sim_loci %>%
  group_by(data, loci_to_average) %>%
  # summarise(
  #   diff = sd(theta_hat)
  # ) %>%
  summarise(
    lower = quantile(theta_hat, 0.02),
    upper = quantile(theta_hat, 0.98),
    diff  = upper - lower
  ) %>%
  select(loci_to_average, data, diff) %>%
  pivot_wider(names_from = data, values_from = diff) %>%
  mutate(
    diff = `With errors` / `No errors`
  )
