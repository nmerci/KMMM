source("R/kmmm.R")
source("R/sampling.R")
source("R/extra.R")

#tested
plot_estimate <- function(x, y)
{
  plot(0, type="n", xlim=c(min(x), max(x)), ylim=c(0, 1))
  for(i in 1:ncol(y))
  {
    lines(x=x, y=y[, i], type="s", col="red")
    lines(x=x, y=pchisq(x, i), type="l", col="blue") #WARNING: this line depends on sampling
  }
}

#partially tested
compare_estimators <- function(n, m, clusters)
{
  data <- sample_data(n, m, clusters)
  
  km <- khizanov_maiboroda(data$censored_data, data$mixture_weights)
  r <- ryzhov(data$censored_data, data$mixture_weights, clusters)
  
  x <- sort(data$censored_data$value)
  #plotting
  #plot_estimate(x, km)
  #plot_estimate(x, r)
  
  #compare deviation
  true_value <- sapply(1:m, pchisq, q=x)
  true_value[is.na(r)] <- NA
  list(km=km - true_value, r=r - true_value)
}

compare_estimators_grouped <- function(n, m, n_groups)
{
  data <- sample_data(n, m)
  
  grouped <- group_mixture_weights(data$mixture_weights, n_groups)
  
  km <- khizanov_maiboroda(data$censored_data, data$mixture_weights)
  r <- ryzhov(data$censored_data, grouped$mixture_weights, grouped$clusters)
  
  x <- sort(data$censored_data$value)
#   #plotting
#   plot_estimate(x, km)
#   plot_estimate(x, r)
  
  #compare deviation
  true_value <- sapply(1:m, pchisq, q=x)
  true_value[is.na(r)] <- NA
  list(km=km - true_value, r=r - true_value)
}


