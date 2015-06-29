library(matrixStats)
source("R/main.R")

test_KMMM <- function(n, calc_norm=T, plot_graphic=F)
{
  m <- 2
  
  data <- cbind(abs(rnorm(n, 0, 2)), rchisq(n, 3))
  
  mixture_probs <- matrix(c(1:n/n, 1 - 1:n/n), n, 2)
  
  mixture_sample <- get_mixture_sample(data, mixture_probs)
  
  censored_sample <- get_censored_sample(mixture_sample, rexp(n, 0.05))
  
  KMMM <- get_KMMM_estimator(censored_sample, mixture_probs)
  
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
    
    norm <- c(norm, get_sup_norm(KMMM[, 2], pnorm(KMMM[, 1], 0, 2) - pnorm(-KMMM[, 1], 0, 2)))
    norm <- c(norm, get_sup_norm(KMMM[, 3], pchisq(KMMM[, 1], 3))) 
    
    # point 1.85
    point <- match(T, KMMM[, 1] > 1.85)
    norm <- c(norm, KMMM[point, 2] - (pnorm(1.85, 0, 2) - pnorm(-1.85, 0, 2)))
    norm <- c(norm, KMMM[point, 3] - pchisq(1.85, 3))
  }  
  
  
  if(plot_graphic==T) 
  {
    #x <- 0:5000/100
    
    plot(x=KMMM[, 1], y=KMMM[, 2] - (pnorm(KMMM[, 1], 0, 2) - pnorm(-KMMM[, 1], 0, 2)), type="s", col="red")
    #lines(x=x, y=pnorm(x, 0, 2) - pnorm(-x, 0, 2))
    
    plot(x=KMMM[, 1], y=KMMM[, 3] - pchisq(KMMM[, 1], 3), type="s", col="red")
    #lines(x=x, y=pchisq(x, 3))
  }
  
  norm  
}

run_test <- function()
{
  medians <- numeric()
  iqrs <- numeric()
  
  n <- c(100, 250, 500, 1000, 2000)
  iter <- 100
  
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

generate_ryzhov_data <- function()
{
  ryzhov_data <- read.csv("data/ryzhov_data.txt", T, "\t")
  
  censored_samples <- list()
  mixture_probs <- numeric()
  
  for(i in 1:27)
  {
    censored_samples[[i]] <- data.frame(time=ryzhov_data[ryzhov_data[, 4] == i, 1],
                                        censored=ryzhov_data[ryzhov_data[, 4] == i, 2])
    
    mixture_probs <- rbind(mixture_probs, c(sum(ryzhov_data[ryzhov_data[, 4] == i, 3] == 1),
                                            sum(ryzhov_data[ryzhov_data[, 4] == i, 3] == 2)))
  }
  
  mixture_probs <- mixture_probs / rowSums(mixture_probs)
  
  list(censored_samples=censored_samples, mixture_probs=mixture_probs)
}

compare_KMMM_Ryzhov <- function(censored_samples, mixture_probs, true_data=NA)
{
  n <- sapply(censored_samples, nrow)
  
  ryzhov <- get_Ryzhov_estimator(censored_samples, mixture_probs)
  
  KMMM <- get_KMMM_estimator(do.call("rbind", censored_samples),
                             mixture_probs[rep(1:length(n), times=n), ])
  KMMM <- KMMM[1:nrow(ryzhov), ]
  
  if(length(true_data) == 0)
    return(list(KMMM, ryzhov))
  else
  {
    true_data <- true_data[1:nrow(ryzhov), ]
    norm <- numeric()
    for(i in 1:ncol(mixture_probs))
      norm <- c(norm, c(max(abs(KMMM[, i + 1] - true_data[, i])), 
                        max(abs(ryzhov[, i + 1] - true_data[, i]))))
    
    return(list(KMMM, ryzhov, norm))
  }
}

test_KMMM_Ryzhov <- function(n, calc_norm=T, plot_graphic=F)
{
  data <- cbind(rchisq(sum(n), 1), rexp(sum(n), 1))
  censors <- runif(sum(n), 0, 5) # ~20% censored
  
  censored_data <- generate_data(n, data, censors)
  
  true_data <- sort(do.call("rbind", censored_data$censored_samples)[, 1])
  true_data <- cbind(pchisq(true_data, 1), pexp(true_data, 1))
  
  KMMM_Ryzhov <- compare_KMMM_Ryzhov(censored_data$censored_samples,
                                     censored_data$mixture_probs,
                                     true_data)
  
  if(plot_graphic == T)
  {
    plot (x=KMMM_Ryzhov[[1]][, 1],
          y=KMMM_Ryzhov[[2]][, 2] - pchisq(KMMM_Ryzhov[[1]][, 1], 1),
          type="s", col="blue")
    lines(x=KMMM_Ryzhov[[1]][, 1], 
          y=KMMM_Ryzhov[[1]][, 2] - pchisq(KMMM_Ryzhov[[1]][, 1], 1),
          type="s", col="red")
    
    plot (x=KMMM_Ryzhov[[1]][, 1],
          y=KMMM_Ryzhov[[2]][, 3] - pexp(KMMM_Ryzhov[[1]][, 1], 1),
          type="s", col="blue")
    lines(x=KMMM_Ryzhov[[1]][, 1], 
          y=KMMM_Ryzhov[[1]][, 3] - pexp(KMMM_Ryzhov[[1]][, 1], 1),
          type="s", col="red")
  }
  
  KMMM_Ryzhov[[3]]
}


iter <- 100
n <- rep(500, 2)
norm <- numeric()

for(i in 1:iter)
{
  norm <- rbind(norm, test_KMMM_Ryzhov(n))
}

colnames(norm) <- c("K1", "R1", "K2", "R2")

norm


par(oma=c(0,0,2,0)) 
par(mfrow=c(2,2))

main_titles <- c("Khizanov", "Ryzhov")
y_titles <- c("ChiSq(1)", "Exp(1)")

for(i in 1:2)
  for(j in 1:2)
    hist(x=norm[, (i - 1)*2 + j], breaks=10, 
         main=main_titles[j],
         xlab=paste("Med=", round(median(norm[, (i - 1)*2 + j], T), 5), 
                    "; IQR=", round(IQR(norm[, (i - 1)*2 + j], T), 5)),
         ylab=y_titles[i])



title(main="2 sampl. with 500 obs. (~20% censored)",outer=T)


