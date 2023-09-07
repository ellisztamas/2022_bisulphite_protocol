#' Simulate multiple loci and estimate mean methylation across them
#'
#' `sim_mean_mC_over_loci` simulates reads for many (default=1e6) loci in
#' the presence of conversion errors drawn from a beta distribution. It then
#' takes a subset of these loci and calculates mean methylation over them
#' using a single *estimate* of the error distribution.
#'
#' @param theta True methylation level to be estimated between zero and one.
#' @param ncytosines Int number of cytosines.
#' @param coverage Float >0 giving mean coverage to cytosine.
#' @param shape1 `shape1` parameter for the Beta distribution.
#' @param shape2 `shape1` parameter for the Beta distribution.
#' @param mean_lambda1 Point estimate for the non-conversion error rate
#' @param ndraws Number of subsamples to draw.
#'
#' @return Dataframe giving input parameters, plus `theta_hat` the maximum-
#' likelihood estimate of mean methylation.
sim_mean_mC_over_loci <- function(theta, ncytosines, coverage, shape1, shape2, mean_lambda1, ndraws){
  # Simulate a very large number of loci
  total_loci <- 1e6
  # Simulate a vector of locus-specific error rates
  lambda1 <- rbeta(total_loci, shape1, shape2)
  # Simulate read counts at each locus.
  y <- sapply(1:total_loci, function(i){
    as.numeric(
      rbinom_with_errors(theta, size = ncytosines*coverage, lambda1 = lambda1[i], lambda2 = 0)
    )
  })
  y <- t(y)

  # Draw successively larger subsamples of loci
  sample_sizes <- c(2, 10^(1:5) )
  sim_loci <- vector('list', length = length(sample_sizes))
  for(i in 1:length(sample_sizes)){
    # Take a sample of all TEs
    these_draws <- matrix(NA, nrow = ndraws, ncol = 3)
    for(draw in 1:ndraws){
      ix <- sample(1:total_loci, sample_sizes[i], replace = FALSE)
      these_draws[draw,] <- colSums( y[ix,] )
    }
    # Estimate methylation given the error-rate point estimate
    theta_hat <- mlbinom_with_errors(these_draws[,1], these_draws[,3], lambda1 = mean_lambda1, lambda2 = 0)
    # Assemble and export parameters
    sim_loci[[i]] <- data.frame(
      theta_real = theta,
      ncytosines = ncytosines,
      coverage =coverage,
      shape1= shape1,
      shape2=shape2,
      mean_lambda1 = mean_lambda1,
      loci_to_average = sample_sizes[i],
      theta_hat
    )
  }
  sim_loci <- do.call('rbind', sim_loci)
}
