# TODO: implement testing framework
source("R/test.R")
source("R/extra.R")

n <- 1000
m <- 2
clusters <- replicate_numbers(2, 500)

compare_estimators(n, m, clusters)
