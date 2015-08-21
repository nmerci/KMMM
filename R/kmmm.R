source("R/extra.R")

library(MASS)
library(stats)

correct_coefficients <- function(ls_coefficients)
{
  n <- ncol(ls_coefficients)
  ls_coefficients <- cbind(ls_coefficients, 0)
  
  for(i in 1:n)
  {
    ls_coefficients[, i] <- ls_coefficients[, i] + ls_coefficients[, n + 1]
    
    ls_coefficients[, n + 1] <- 0
    ls_coefficients <- swap_columns(ls_coefficients, i, n + 1, ls_coefficients[, i] < 0)
  }
  
  ls_coefficients[, 1:n]
}

#grouping works for nx2 weight matrices only
group_mixture_weights <- function(mixture_weights, n_groups)
{
  clusters <- rep(0, nrow(mixture_weights))
  
  for(i in 1:nrow(mixture_weights))
    clusters[i] <- match(T, mixture_weights[i, 1] < (1:n_groups) / n_groups)

  for(i in 1:n_groups)
  {
    mixture_weights[i, 1] <- i/n_groups - 0.5/n_groups
    mixture_weights[i, 2] <- 1 - mixture_weights[i, 1]
  }
  
  list(mixture_weights=mixture_weights[sort(unique(clusters)), ], clusters=clusters)
}

#tested
khizanov_maiboroda <- function(censored_data, mixture_weights)
{
  n <- nrow(censored_data)
  
  #sort sample and weights
  ord <- order(censored_data$value)
  censored_data <- censored_data[ord, ]
  mixture_weights <- mixture_weights[ord, ]
  
  #calculate least squares solution
  ls_coefficients <- ginv(mixture_weights)
  
  #correction of CDFs
  #ls_coefficients <- correct_coefficients(ls_coefficients)
  
  #calculate modified Kaplan-Meier estimator
  1 - apply(1 - t(t(ls_coefficients) * censored_data$is_censored) / 
              (1 - cbind(0, t(apply(ls_coefficients, 1, cumsum))[, -n])), 1, cumprod)
}

#tested
kaplan_meier <- function(censored_data)
{
  n <- nrow(censored_data)
  
  #sort sample
  censored_data <- censored_data[order(censored_data$value), ] 
  
  #Kaplan-Meier estimator
  s <- cumprod(1 - censored_data$is_censored / (n:1)) 
  
  #Greenwood estimator
  var <- s^2 * cumsum(censored_data$is_censored / ((n:1) * (n:1 - 1))) 
  
  data.frame(x=censored_data$value, cdf=1-s, var=var)
}

#tested
ryzhov <- function(censored_data, mixture_weights, clusters)
{
  n <- nrow(censored_data)
  k <- nrow(mixture_weights)
  
  #run Kaplan-Meier estimator for each cluster
  km <- list()
  for(i in 1:k)
    km[[i]] <- kaplan_meier(censored_data[clusters==i, ])

  #calculate generalized least squares solution for each point in data set
  censored_data$value <- sort(censored_data$value)
  result <- numeric()
  for(i in 1:n)
  {
    #construct response vector and covariance matrix
    y <- numeric()
    wt <- numeric()
    for(j in 1:k)
    {
      #determine position of the current point
      pos <- match(FALSE, censored_data$value[i] >= km[[j]]$x, nrow(km[[j]]))
      
      y <- c(y, km[[j]]$cdf[pos])
      wt <- c(wt, 1 / km[[j]]$var[pos])
    }
    
    #calculate GLS coefficients (if possible)
    if(!sum(is.na(wt)) && max(wt, T) < Inf) #this condition could be improved
      result <- rbind(result, lsfit(mixture_weights, y, wt, F)$coef)
    else
      result <- rbind(result, rep(NA, ncol(mixture_weights)))
  }
  
  colnames(result) <- NULL
  result
}
