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
  
  #plotting
  plot_estimate(sort(data$censored_data$value), km)
  plot_estimate(sort(data$censored_data$value), r)
  
  #compare deviation
  #TODO: implement some kind of goodness of fit
}



