# TODO: implement testing framework
source("R/test.R")

library(MASS)
library(stats)
library(matrixStats)

result <- list()

for(s in 2:5)
{
  result[[s]] <- numeric()
  for(j in 0:9)
  {
    print(paste("j=", j, sep=""))
    
    n <- 100 + 100*j
    m <- 2
    
    sup_km <- numeric()
    sup_r <- numeric()
    for(i in 1:1000)
    {
      res <- sapply(compare_estimators_grouped(n, m, s), abs, simplify=F)
      #res <- sapply(compare_estimators(n, m, clusters), "^", 2, simplify = F)
      
      sup_km <- rbind(sup_km, apply(res$km, 2, max, na.rm=T))
      sup_r <- rbind(sup_r, apply(res$r, 2, max, na.rm=T))
    }
    
    sup_km[sup_km > max(sup_r)] <- NA
    
    result[[s]] <- rbind(result[[s]], cbind(t(colMedians(sup_km, na.rm=T)), t(colMedians(sup_r, na.rm=T))))
  }
}





#temp plot
n <- 100 * (1:10)
for(s in 2:5)
{
  plot(x=n, y=result[[s]][, 1], type="o", col="blue", main=paste("ChiSq(1), N =", s), xlab="n", ylab="E[sup|F_est-F|]", 
       ylim=c(min(result[[s]][, c(1, 3)]), max(result[[s]][, c(1, 3)])))
  lines(x=n, y=result[[s]][, 3], type="o", pch=22, lty=2, col="red")
  legend(x=750, y=0.18, c("P-L", "R"), col=c("blue", "red"), pch=c(21, 22), lty=c(1, 2))
  
  plot(x=n, y=result[[s]][, 2], type="o", col="blue", main=paste("ChiSq(2), N =", s), xlab="n", ylab="E[sup|F_est-F|]", 
       ylim=c(min(result[[s]][, c(2, 4)]), max(result[[s]][, c(2, 4)])))
  lines(x=n, y=result[[s]][, 4], type="o", pch=22, lty=2, col="red")
  legend(x=750, y=0.18, c("P-L", "R"), col=c("blue", "red"), pch=c(21, 22), lty=c(1, 2))
}



