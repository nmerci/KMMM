source("R/extra.R")

#tested
sample_mixture_weights <- function(n, m)
{
  mixture_weights <- matrix(runif(n*m), n, m)
  mixture_weights <- mixture_weights / rowSums(mixture_weights)
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
sample_data <- function(n, m, clusters)
{
  #sample data from ChiSquared distribution
  #with different parameter
  #NOTE: another distribution could be used
  data <- numeric()
  for(i in 1:m)
    data <- cbind(data, rchisq(n, i))
  
  #sample mixture weights with respect to clusters
  k <- max(clusters)
  mixture_weights <- replicate_rows(sample_mixture_weights(k, m), get_frequencies(clusters))
  
  #sample data from mixture distribution
  #and censor it using Uniform distribution
  #NOTE: another distribution could be used
  data <- sample_mixture_distribution(data, mixture_weights)
  censored_data <- sample_censored_data(data, runif(n, 0, 10))
  
  list(censored_data=censored_data, mixture_weights=mixture_weights)
}