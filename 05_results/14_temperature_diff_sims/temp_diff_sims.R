#' Script to simulate differences in methylation at two temperatures
#'
#' Laura wants to characterise methylation levels of TEs in plants grown at cold
#' and warm temperatures. How well we can do this depends on:
#' 1. The absolute level
#' 2. The difference between the two temperatures
#' 3. Number of trials = Coverage per cytosine * number of cytosines
#' 4. Conversion errors
#'
#' This script simulates a series of mean methylation levels at one temperature,
#' a deviation from that mean for a second temperature, errors (zero or 5%), and
#' increasing number of trials.

source("02_library/R/binomial_with_errors.R")

# Number of trials is the product of how many cytosines, and coverage per C
# In TAIR10, the minimum TE size is 138bp, mean is 798, and max is 31000
n_cytosines <- c(1e2, 1e3, 1e4, 1e5)
coverage <- 1
# True mean methylation for temperature A
real_theta_a <- c(0.01, 0.05, 0.1, 0.15, 0.2)
# The deviation in mean methylation for temperture B from temp A
# i.e. methylation at temp B = methylation at temp A + these values.
diff <- c(0.02, 0.04, 0.08, 0.16)
# Number of replicate simulations for each set of combinations.
nreps <- 1000

# Empty list to store the results
nsims <- length(n_cytosines) * 2 * length(diff) * length(real_theta_a)
inferred_theta <- vector('list', length = nsims)

counter <- 0
for(n in n_cytosines){
  for(l in c(0,1)){
    for(d in diff){
      for(ra in real_theta_a){
        counter <- counter + 1
        cat(counter, " ")
        # difference in methylation for the two samples
        real_theta_b <- ra + d

        # Simulate methylated reads with conversion errors
        ya <- rbinom_with_errors(
          ra, lambda1 = l*0.072, lambda2 = 0, size = n * coverage, nreps = nreps
          # ra, lambda1 = l*rbeta(1,17,220), lambda2 = 0, size = n * coverage, nreps = nreps
        )
        yb <- rbinom_with_errors(
          # real_theta_b, lambda1 = l*rbeta(1,17,220), lambda2 = 0, size = n * coverage, nreps = nreps
          real_theta_b, lambda1 = l*0.072, lambda2 = 0, size = n * coverage, nreps = nreps
        )

        # Maximum-likelihood true means, accounting for errors
        inf_theta_a = sapply(1:nreps, function(i){
          mlbinom_with_errors(ya$mC[i], ya$n[i], lambda1 = l*0.072, lambda2 = 0)
        })
        inf_theta_b = sapply(1:nreps, function(i) {
          mlbinom_with_errors(yb$mC[i], yb$n[i], lambda1 = l*0.072, lambda2 = 0)
        })
        # Output the results for this set of parameters
        inferred_theta[[counter]] <- data.frame(
          size = n * coverage,
          lambda = l,
          real_theta_a = ra,
          real_theta_b = real_theta_b,
          inf_theta_a = inf_theta_a,
          inf_theta_b = inf_theta_b
        )
      }
    }
  }
}

inferred_theta <- do.call('rbind', inferred_theta)

write_csv(x = inferred_theta,
          "05_results/14_temperature_diff_sims/sim_temp_diffs.csv"
)
