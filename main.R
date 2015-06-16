get_mixture_sample <- function(n, data, mixture_probs)
{
  result <- numeric()
  
  for(i in 1:n)
    result <- c(result, as.numeric(sample(data[i, ], 1, F, mixture_probs[i, ])))
  
  result
}

get_mixture_weights <- function(mixture_probs)
{
  solve(t(mixture_probs) %*% mixture_probs) %*% t(mixture_probs)  
}

get_censored_sample <- function(sample, censors)
{
  n <- length(sample)
  
  censored <- rep(1, n)
  censored[sample > censors] <- 0
  sample[sample > censors] <- censors[sample > censors]
  
  data.frame(time=sample, censored=censored)
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
  
  # remove duplicates
  unique_times <- !duplicated(censored_sample$time, fromLast=T)
  S <- S[unique_times]
  v <- v[unique_times]
  
  # return
  data.frame(time=censored_sample$time[unique_times], F=1-S, Var=v * S^2)
}

get_KMM_estimator <- function(censored_sample, mixture_probs, correction="n")
{
  n <- length(censored_sample$time)
  
  permutation <- order(censored_sample$time)
  
  censored_sample$time <- censored_sample$time[permutation]
  censored_sample$censored <- censored_sample$censored[permutation]
  mixture_weights <- get_mixture_weights(mixture_probs)[, permutation]
  
  # correction
  if(correction == "u")
  {
    redundancy <- 0
    for(i in 1:n)
    {
      mixture_weights[, i] <- mixture_weights[, i] + redundancy
      
      redundancy <- mixture_weights[, i]
      
      redundancy[redundancy > 0] <- 0
      mixture_weights[mixture_weights[, i] < 0, i] <- 0
    }
    
    redundancy <- rowSums(mixture_weights)
    for(i in n:1)
    {
      overlap <- redundancy > 1
      
      if(sum(overlap))
      {
        redundancy[overlap] <- redundancy[overlap] - mixture_weights[overlap, i]
        mixture_weights[overlap, i] <- 0
      }
      else
        break      
    }
  }
  
  if(correction == "l")
  {
    redundancy <- 0
    for(i in n:1)
    {
      mixture_weights[, i] <- mixture_weights[, i] + redundancy
      
      redundancy <- mixture_weights[, i]
      
      redundancy[redundancy > 0] <- 0
      mixture_weights[mixture_weights[, i] < 0, i] <- 0
    }
    
    redundancy <- rowSums(mixture_weights)
    for(i in 1:n)
    {
      overlap <- redundancy < 0
      
      if(sum(overlap))
      {
        redundancy[overlap] <- redundancy[overlap] - mixture_weights[overlap, i]
        mixture_weights[overlap, i] <- 0
      }
      else
        break      
    }
  }
  
  
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

test <- function(n, calc_norm=T, correction="n", plot_graphic=F)
{
  m <- 2
  
  data <- cbind(abs(rnorm(n, 0, 2)), rchisq(n, 3))
  
  mixture_probs <- matrix(c(1:n/n, 1 - 1:n/n), n, 2)
  
  mixture_sample <- get_mixture_sample(n, data, mixture_probs)
  
  censored_sample <- get_censored_sample(mixture_sample, rexp(n, 0.1))
  
  KMM <- get_KMM_estimator(censored_sample, mixture_probs, correction)
  
  sorted_sample <- sort(censored_sample$time)
  
  norm <- numeric()
  if(calc_norm == T)
  {
    #     # |N(0, 2)| quantile = 3.289709
    #     # Chi(3)    quantile = 6.251388
    #     KMM1 <- KMM[1, sorted_sample < 3.289709]
    #     sample1 <- sorted_sample[sorted_sample < 3.289709]
    #     
    #     KMM2 <- KMM[2, sorted_sample < 6.251388]
    #     sample2 <- sorted_sample[sorted_sample < 6.251388]
    #     
    #     norm <- c(norm, get_sup_norm(KMM1, pnorm(sample1, 0, 2) - pnorm(-sample1, 0, 2)))
    #     norm <- c(norm, get_sup_norm(KMM2, pchisq(sample2, 3)))    
    
    norm <- c(norm, get_sup_norm(KMM[1, ], pnorm(sorted_sample, 0, 2) - pnorm(-sorted_sample, 0, 2)))
    norm <- c(norm, get_sup_norm(KMM[2, ], pchisq(sorted_sample, 3))) 
    
    # point 1.85
    point <- match(T, sorted_sample > 1.85)
    norm <- c(norm, KMM[1, point] - (pnorm(1.85, 0, 2) - pnorm(-1.85, 0, 2)))
    norm <- c(norm, KMM[2, point] - pchisq(1.85, 3))
  }  
  
  
  if(plot_graphic==T) 
  {
    x <- 0:5000/100
    
    plot(x=sorted_sample, y=KMM[1, ], type="s", col="red")
    lines(x=x, y=pnorm(x, 0, 2) - pnorm(-x, 0, 2))
    
    plot(x=sorted_sample, y=KMM[2, ], type="s", col="red")
    lines(x=x, y=pchisq(x, 3))
  }
  
  norm  
}

run_test <- function(correction)
{
  medians <- numeric()
  iqrs <- numeric()
  
  n <- c(100, 250, 500, 1000, 2000)
  iter <- 1000
  
  for(i in 1:length(n))
  {  
    norm <- numeric()
    for(j in 1:iter)
      norm <- rbind(norm, test(n[i], T, correction))
    
    medians <- rbind(medians, colMedians(norm, T))
    iqrs <- rbind(iqrs, colIQRs(norm, T))
    
    #print("*")
  }
  
  print(correction)
  print(medians)
  print(iqrs)
}

library(matrixStats)

correction <- c("n", "l", "u")

for(i in 1:3)
  run_test(correction[i])


