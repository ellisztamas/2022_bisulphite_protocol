#' Functions to work with binomial distributions when conversion errors are
#' present
#' Tom Ellis

# Sum over a vector of log values
logsumexp <- function(x){
  max(x) + log(sum(exp(x - max(x))))
}

#' (log) likelihood of a binomial model with errors
#'
#' @param y Int number of successes. Must be less than n
#' @param n Int number of trials.
#' @param p True mean of the underlying binomial process
#' @param lambda1 Float between 0 and 1 giving the probability that an
#' unmethylated read is appears methylated
#' @param lambda2 Float between 0 and 1 giving the probability that a methylated
#' read is appears unmethylated
#'
#' @return Log likelihood of the data given p
#' @author Tom Ellis
llbinom_with_errors <- function(y, n, p, lambda1, lambda2 = lambda1){
  loglik <- lchoose(n = n, k = y) +
    log( p*(1-lambda1) + (1-p) *    lambda2  ) * y  +
    log( p *  lambda1  + (1-p) * (1-lambda2) ) * (n-y)

  loglik[n==0] <- 0

  return(loglik)
}

#' Random number generator for the binomial distribution with errors
#'
#' @param p True mean of the underlying binomial process
#' @param lambda1 Float between 0 and 1 giving the probability that an
#' unmethylated read is appears methylated
#' @param lambda2 Float between 0 and 1 giving the probability that a methylated
#' read is appears unmethylated
#' @param size number of trials (zero or more)
#' @param nreps Number of draws (length of the output vector)
#'
#' @return Vector of integers
#' @author Tom Ellis
rbinom_with_errors <- function(p, lambda1, lambda2, size, nreps=1){
  sim_mean <- (p*(1-lambda1)) + ((1-p)*lambda2)
  mC <- rbinom(n=nreps, size=size, prob = sim_mean)
  uC <- size - mC

  list(
    mC = mC,
    uC = uC,
    n = mC + uC
  )
}

#' Maximum-likelihood value for a binomial process with errors
#'
#' Grid approximation
#'
#' @param y Int number of successes. Must be less than n
#' @param n Int number of trials.
#' @param lambda1 Float between 0 and 1 giving the probability that an
#' unmethylated read is appears methylated
#' @param lambda2 Float between 0 and 1 giving the probability that a methylated
#' read is appears unmethylated
#'
#' @return Float giving the maximum likelihood value.
mlbinom_with_errors <- function(y, n, lambda1, lambda2=0){
  p <- y/n
  theta <-  (lambda1 - p) / (lambda1 + lambda2 - 1)

  # If estimates fall outside [0,1], set them articifically.
  theta[theta<0] <- 0
  theta[theta>1] <- 1

  return(theta)
}

#' Maximum-likelihood value for a binomial process with uncertainty in errors
#'
#' Calculate the maximum-likelihood true mean for a binomial process in the
#' presence of false-positive successes. This integrates out uncertainty in the
#' true false positive rate by modelling it as a beta distribution with shape
#' parameters supplied by the used. At present this is only implemented for
#' false positive errors.
#'
#' `mlbinom_with_uncertainty` performs grid interpolation over possible values
#' of the error rate between zero and one, gets the maximum-likelihood value for
#' the true mean for each, and weights each estimate by the probability density
#' for the error-rate value. It returns the mean with the highest likelihood.
#'
#' @param y Int number of successes. Must be less than n
#' @param n Int number of trials.
#' @param alpha2,beta2 Floats giving the shape parameters for a beta
#' distribution describing uncertainty in false methylation.
#' @param precision Float passed to `by` argument of `seq` to state how fine a
#' grid to interpolate over.
#'
#' @return Float giving the maximum likelihood value.
#' @author Tom Ellis
mlbinom_with_uncertainty <- function(y, n, alpha2=NULL, beta2=NULL, precision=1e-3){
  p <- y/n
  lambda_vals <- seq(from = precision, to = 1-precision, by = precision)
  pr_lambda <- dbeta(lambda_vals, alpha2, beta2)
  theta_hat <- (lambda_vals-p) / (lambda_vals-1)
  lik <- (theta_hat * pr_lambda)

  theta_hat[ which(lik == max(lik)) ]
}
