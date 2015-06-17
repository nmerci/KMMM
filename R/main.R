get_mixture_sample <- function(data, mixture_probs)
{
  n <- nrow(data)
  
  result <- numeric()
  
  for(i in 1:n)
    result <- c(result, as.numeric(sample(data[i, ], 1, F, mixture_probs[i, ])))
  
  result
}

get_mixture_samples <- function(n, data, mixture_probs)
{
  m <- ncol(mixture_probs)
  
  mixture_samples <- list()
  
  accum <- 0
  for(i in 1:length(n))
  {
    mixture_samples[[i]] <- get_mixture_sample(data[(accum + 1):(accum + n[i]), ], 
                                               mixture_probs[rep(i, n[i]), ])
    accum <- accum + n[i]
  }
  
  mixture_samples
}

get_mixture_weights <- function(mixture_probs)
{
  solve(t(mixture_probs) %*% mixture_probs) %*% t(mixture_probs)  
}

get_weighted_mixture_weights <- function(mixture_probs, mixture_probs_weights)
{
  # Here weights corresponds to precision matrix
  
  # regularization
  mixture_probs_weights <- mixture_probs_weights / mean(mixture_probs_weights)
  
  solve(t(mixture_probs) %*% mixture_probs_weights %*% mixture_probs) %*% 
    t(mixture_probs) %*% mixture_probs_weights
}

get_censored_sample <- function(sample, censors)
{
  n <- length(sample)
  
  censored <- rep(1, n)
  censored[sample > censors] <- 0
  sample[sample > censors] <- censors[sample > censors]
  
  data.frame(time=sample, censored=censored)
}

get_censored_samples <- function(samples, censors)
{
  censored_samples <- list()
  
  accum <- 0
  for(i in 1:length(samples))
  {
    censored_samples[[i]] <- get_censored_sample(samples[[i]], 
                                                 censors[(accum + 1):(accum + n[i])])
    accum <- accum + n[i]
  }
  
  censored_samples
}

get_KM_GW_estimator <- function(censored_sample)
{
  n <- length(censored_sample$time)
  
  # sort sample
  permutation <- order(censored_sample$time)  
  censored_sample$time <- censored_sample$time[permutation]
  censored_sample$censored <- censored_sample$censored[permutation]
  
  # estimate survival function
  S <- 1 - censored_sample$censored[1] / n  
  for(i in 2:n)
    S[i] <- S[i - 1] * (1 - censored_sample$censored[i] / (n - i + 1))
  
  # estimate variance
  v <- censored_sample$censored[1] / (n * (n - censored_sample$censored[1]))
  for(i in 2:n)
    v[i] <- v[i - 1] + censored_sample$censored[i] / 
    ((n - i + 1) * (n - i + 1 - censored_sample$censored[i]))
  
  v[n] <- v[n - 1]
  
  # return
  data.frame(time=censored_sample$time[unique_times], F=1-S, Var=v * S^2)
}

get_KMM_estimator <- function(censored_sample, mixture_probs)
{
  n <- length(censored_sample$time)
  
  permutation <- order(censored_sample$time)
  
  censored_sample$time <- censored_sample$time[permutation]
  censored_sample$censored <- censored_sample$censored[permutation]
  mixture_weights <- get_mixture_weights(mixture_probs)[, permutation]
  
  mixture_weights_sums <- cbind(0, mixture_weights[, 1:(n-1)])
  for(i in 3:n)
    mixture_weights_sums[, i] <- mixture_weights_sums[, i] + mixture_weights_sums[, i - 1]
  
  H <- 1 - t(t(mixture_weights) * censored_sample$censored) / (1 - mixture_weights_sums)
  for(i in 2:n)
    H[, i] <- H[, i] * H[, i - 1]
  
  1 - H  
}

get_sup_norm <- function(y_estimate, y_true)
{
  n <- length(y_estimate)
  y_estimate <- c(0, y_estimate)
  
  max(c(abs(y_estimate[1:n] - y_true),
        abs(y_estimate[2:(n+1)] - y_true)))
}