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
  # Preconditioning
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
                                                 censors[(accum + 1):(accum + length(samples[[i]]))])
    accum <- accum + length(samples[[i]])
  }
  
  censored_samples
}

sort_censored_sample <- function(censored_sample, mixture_probs=numeric())
{
  permutation <- order(censored_sample$time) 

  censored_sample$time <- censored_sample$time[permutation]
  censored_sample$censored <- censored_sample$censored[permutation]
  
  if(length(mixture_probs))
  {
    mixture_weights <- get_mixture_weights(mixture_probs)[, permutation]
    censored_sample <- list(sample=censored_sample, weights=mixture_weights)
  }

  censored_sample
}

get_KM_GW_estimator <- function(censored_sample)
{
  n <- length(censored_sample$time)
  
  # sort sample
  censored_sample <- sort_censored_sample(censored_sample)
  
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
  data.frame(time=censored_sample$time, F=1-S, Var=v * S^2)
}

get_Ryzhov_estimator <- function(censored_samples, mixture_probs)
{
  KM_GW_estimators <- list()
  for(i in 1:length(censored_samples))
    KM_GW_estimators[[i]] <- get_KM_GW_estimator(censored_samples[[i]])
  
  time <- numeric()
  H <- numeric()

  for(i in 1:length(KM_GW_estimators))
  {
    for(j in 1:length(KM_GW_estimators[[i]]$time))
    {
      t0 <- KM_GW_estimators[[i]]$time[j]
      
      KM_F <- numeric()
      KM_Var <- numeric()
      
      for(k in 1:length(KM_GW_estimators))
      {
        t0_position <- match(F, t0 > KM_GW_estimators[[k]]$time)
        if(is.na(t0_position))
          t0_position <- length(KM_GW_estimators[[k]]$time)
        
        KM_F <- c(KM_F, KM_GW_estimators[[k]]$F[t0_position])
        KM_Var <- c(KM_Var, KM_GW_estimators[[k]]$Var[t0_position])
      }
      
      precision <- diag(1/KM_Var)
      if(max(precision) < Inf)
      {
        time <- c(time, t0)
        H <- cbind(H, get_weighted_mixture_weights(mixture_probs, precision) %*% KM_F)
      }
    }
  }
  
  permutation <- order(time)
  H <- H[, permutation]
 
  cbind(time[permutation], t(H))
}

get_KMMM_estimator <- function(censored_sample, mixture_probs)
{
  n <- length(censored_sample$time)
  
  # sort sample
  censored_sample <- sort_censored_sample(censored_sample, mixture_probs)
  mixture_weights <- censored_sample$weights
  censored_sample <- censored_sample$sample
  
  # estimate survival function
  mixture_weights_sums <- cbind(0, mixture_weights[, 1:(n-1)])
  for(i in 3:n)
    mixture_weights_sums[, i] <- mixture_weights_sums[, i] + mixture_weights_sums[, i - 1]
  
  H <- 1 - t(t(mixture_weights) * censored_sample$censored) / (1 - mixture_weights_sums)
  for(i in 2:n)
    H[, i] <- H[, i] * H[, i - 1]
  
  cbind(censored_sample$time, t(1 - H))
}

get_sup_norm <- function(y_estimate, y_true)
{
  n <- length(y_estimate)
  y_estimate <- c(0, y_estimate)
  
  max(c(abs(y_estimate[1:n] - y_true),
        abs(y_estimate[2:(n+1)] - y_true)))
}
