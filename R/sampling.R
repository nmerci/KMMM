
sample_mixture_weights <- function(n, m)
{
  #TODO: adopt for ryzhov estimator
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
censor_data <- function(data, censors)
{
  n <- length(data)
  
  is_censored <- rep(F, n)
  is_censored[data > censors] <- T
  data[is_censored] <- censors[is_censored]
  
  data.frame(value=data, is_censored=1-is_censored)
}
