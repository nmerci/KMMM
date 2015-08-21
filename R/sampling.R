source("R/extra.R")

#tested
sample_mixture_weights <- function(n, m)
{
  mixture_weights <- matrix(0)
  
  while(det(t(mixture_weights) %*% mixture_weights) < 0.2)
  {
    mixture_weights <- matrix(runif(n*m), n, m)
    mixture_weights <- mixture_weights / rowSums(mixture_weights)
  }
  
  mixture_weights
}

#tested
sample_mixture_distribution <- function(data, mixture_weights)
{
  n <- nrow(data)
  result <- numeric()
  for(i in 1:n)
    result <- c(result, sample(data[i, ], 1, prob=mixture_weights[i, ]))
  
  as.vector(result)
}

#tested
sample_censored_data <- function(data, censors)
{
  n <- length(data)
  
  is_censored <- rep(F, n)
  is_censored[data > censors] <- T
  data[is_censored] <- censors[is_censored]
  
  data.frame(value=data, is_censored=1-is_censored)
}

#tested
sample_data <- function(n, m)
{
  data <- numeric()
  for(i in 1:m)
    data <- cbind(data, rchisq(n, i))
  
  mixture_weights <- sample_mixture_weights(n, m)
  
  data <- sample_mixture_distribution(data, mixture_weights)
  censored_data <- sample_censored_data(data, runif(n, 0, 8))
  
  list(censored_data=censored_data, mixture_weights=mixture_weights)
}
