source("main.R")

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
    
    print("*")
  }
  
  print(correction)
  print(medians)
  print(iqrs)
}

library(matrixStats)

correction <- c("n", "l", "u")

for(i in 1:3)
  run_test(correction[i])
