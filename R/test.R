library(matrixStats)
source("R/main.R")

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
    s_norm <- numeric()
    p_norm <- numeric()
    for(i in 1:ncol(mixture_probs))
    {
      s_norm <- c(s_norm, c(max(abs(KMMM[, i + 1] - true_data[, i])), 
                            max(abs(ryzhov[, i + 1] - true_data[, i]))))
      
      pos <- match(T, 1 < KMMM[, 1])
      p_norm <- c(p_norm, c(abs(KMMM[pos, i + 1] - true_data[pos, i]),
                            abs(ryzhov[pos, i + 1] - true_data[pos, i])))
    }
      
    
    return(list(KMMM, ryzhov, s_norm, p_norm))
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
          type="s", col="blue", main="ChiSq")
    lines(x=KMMM_Ryzhov[[1]][, 1], 
          y=KMMM_Ryzhov[[1]][, 2] - pchisq(KMMM_Ryzhov[[1]][, 1], 1),
          type="s", col="red")
    
    plot (x=KMMM_Ryzhov[[1]][, 1],
          y=KMMM_Ryzhov[[2]][, 3] - pexp(KMMM_Ryzhov[[1]][, 1], 1),
          type="s", col="blue", main="Exp")
    lines(x=KMMM_Ryzhov[[1]][, 1], 
          y=KMMM_Ryzhov[[1]][, 3] - pexp(KMMM_Ryzhov[[1]][, 1], 1),
          type="s", col="red")
  }
  
  c(KMMM_Ryzhov[[3]], KMMM_Ryzhov[[4]])
}


f <- function(size, groups)
{
  iter <- 1000
  n <- rep(size, groups)
  
  norm <- numeric()
  for(i in 1:iter)
    norm <- rbind(norm, test_KMMM_Ryzhov(n))
  
  colMedians(norm)
}

tab1 <- numeric()
tab2 <- numeric()
for(i in 1:5)
{
  tab1 <- rbind(tab1, f(50, i+1))
  tab2 <- rbind(tab2, f(50*i, 2))
}
  




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


