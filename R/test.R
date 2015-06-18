library(matrixStats)
source("R/main.R")

test_KMMM <- function(n, calc_norm=T, plot_graphic=F)
{
  m <- 2
  
  data <- cbind(abs(rnorm(n, 0, 2)), rchisq(n, 3))
  
  mixture_probs <- matrix(c(1:n/n, 1 - 1:n/n), n, 2)
  
  mixture_sample <- get_mixture_sample(data, mixture_probs)
  
  censored_sample <- get_censored_sample(mixture_sample, rexp(n, 0.1))
  
  KMMM <- get_KMMM_estimator(censored_sample, mixture_probs)
  
  sorted_sample <- sort(censored_sample$time)
  
  norm <- numeric()
  if(calc_norm == T)
  {
    #     # |N(0, 2)| quantile = 3.289709
    #     # Chi(3)    quantile = 6.251388
    #     KMMM1 <- KMMM[1, sorted_sample < 3.289709]
    #     sample1 <- sorted_sample[sorted_sample < 3.289709]
    #     
    #     KMMM2 <- KMMM[2, sorted_sample < 6.251388]
    #     sample2 <- sorted_sample[sorted_sample < 6.251388]
    #     
    #     norm <- c(norm, get_sup_norm(KMMM1, pnorm(sample1, 0, 2) - pnorm(-sample1, 0, 2)))
    #     norm <- c(norm, get_sup_norm(KMMM2, pchisq(sample2, 3)))    
    
    norm <- c(norm, get_sup_norm(KMMM[1, ], pnorm(sorted_sample, 0, 2) - pnorm(-sorted_sample, 0, 2)))
    norm <- c(norm, get_sup_norm(KMMM[2, ], pchisq(sorted_sample, 3))) 
    
    # point 1.85
    point <- match(T, sorted_sample > 1.85)
    norm <- c(norm, KMMM[1, point] - (pnorm(1.85, 0, 2) - pnorm(-1.85, 0, 2)))
    norm <- c(norm, KMMM[2, point] - pchisq(1.85, 3))
  }  
  
  
  if(plot_graphic==T) 
  {
    x <- 0:5000/100
    
    plot(x=sorted_sample, y=KMMM[1, ], type="s", col="red")
    lines(x=x, y=pnorm(x, 0, 2) - pnorm(-x, 0, 2))
    
    plot(x=sorted_sample, y=KMMM[2, ], type="s", col="red")
    lines(x=x, y=pchisq(x, 3))
  }
  
  norm  
}

run_test <- function()
{
  medians <- numeric()
  iqrs <- numeric()
  
  n <- c(100, 250, 500, 1000, 2000)
  iter <- 1000
  
  for(i in 1:length(n))
  {  
    norm <- numeric()
    for(j in 1:iter)
      norm <- rbind(norm, test_KMMM(n[i]))
    
    medians <- rbind(medians, colMedians(norm, T))
    iqrs <- rbind(iqrs, colIQRs(norm, T))
    
    print("*")
  }
  
  print(medians)
  print(iqrs)
}

test_Ryzhov <- function()
{
  n <- c(30, 20, 50, 40, 30)
  m <- 2
  
  data <- cbind(runif(sum(n)), rchisq(sum(n), 1))
  
  mixture_probs <- matrix(runif(length(n) * m), length(n), m)
  mixture_probs <- mixture_probs / rowSums(mixture_probs)
  
  censors <- runif(sum(n))
  
  mixture_samples <- get_mixture_samples(n, data, mixture_probs)
  censored_samples <- get_censored_samples(mixture_samples, censors)
  
  ryzhov <- get_Ryzhov_estimator(censored_samples, mixture_probs)
}

generate_data <- function(n, data, censors)
{
  m <- ncol(data)
  
  mixture_probs <- matrix(runif(length(n) * m), length(n), m)
  mixture_probs <- mixture_probs / rowSums(mixture_probs)
  
  mixture_samples <- get_mixture_samples(n, data, mixture_probs)
  censored_samples <- get_censored_samples(mixture_samples, censors)
  
  list(censored_samples=get_censored_samples(mixture_samples, censors),
       mixture_probs=mixture_probs)
}

compare_KMMM_Ryzhov <- function(censored_samples, mixture_probs)
{
  n <- sapply(censored_samples, nrow)
  
  ryzhov <- get_Ryzhov_estimator(censored_samples, mixture_probs)
  
  KMMM <- get_KMMM_estimator(do.call("rbind", censored_samples),
                             mixture_probs[rep(1:length(n), times=n), ])
  
  list(KMMM, ryzhov)
}





