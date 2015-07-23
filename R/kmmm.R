library(MASS) # ginv

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
khizanov_maiboroda <- function(censored_data, mixture_weights)
{
  n <- nrow(censored_data)
  
  #sort sample and weights
  ord <- order(censored_data$value)
  censored_data <- censored_data[ord, ]
  mixture_weights <- mixture_weights[ord, ]
  
  #calculate least squares solution
  ls_coefficients <- ginv(mixture_weights)
  
  #correction of CDFs can be made here
  # >>>
  
  #calculate modified Kaplan-Meier estimator
  1 - apply(1 - t(t(ls_coefficients) * censored_data$is_censored) / 
              (1 - cbind(0, t(apply(ls_coefficients, 1, cumsum))[, -n])), 1, cumprod)
}

#tested
ryzhov <- function(censored_data, mixture_weights, cluster)
{
  n <- nrow(mixture_weights)
  
  #clustering algorithm can be added here
  # >>>
  
  #for now clusters are assumed to be known
  k <- max(cluster)
  
  #run Kaplan-Meier estimator for each cluster
  #and construct feature matrix
  km <- list()
  x <- numeric()
  for(i in 1:k)
  {
    km[[i]] <- kaplan_meier(censored_data[cluster==i, ])
    x <- rbind(x, colMeans(mixture_weights[cluster==i, ]))
  }
  
  #calculate generalized least squares solution
  #for each point in data set
  censored_data$value <- sort(censored_data$value)
  
  result <- numeric()
  for(i in 1:n)
  {
    #construct response vector and covariance matrix
    y <- numeric()
    wt <- numeric()
    for(j in 1:k)
    {
      pos <- match(FALSE, censored_data$value[i] >= km[[j]]$x, nrow(km[[j]]))
      y <- c(y, km[[j]]$cdf[pos])
      wt <- c(wt, 1 / km[[j]]$var[pos])
    }
    
    #calculate GLS coefficients
    if(!sum(is.na(wt)) && max(wt, T) < Inf)
      result <- rbind(result, lsfit(x, y, wt, F)$coef)
    else
      result <- rbind(result, NA)
  }
  
  result
}
