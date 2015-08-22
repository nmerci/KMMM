# TODO: implement testing framework
source("R/test.R")

library(matrixStats)

result <- numeric()
for(j in 0:9)
{
  print(paste("j=", j, sep=""))
  
  n <- 100 + 100*j
  m <- 2
  
  sup_km <- numeric()
  sup_r <- numeric()
  for(i in 1:1000)
  {
    res <- sapply(compare_estimators_grouped(n, m, 5), abs, simplify=F)
    #res <- sapply(compare_estimators(n, m, clusters), "^", 2, simplify = F)
    
    sup_km <- rbind(sup_km, apply(res$km, 2, max, na.rm=T))
    sup_r <- rbind(sup_r, apply(res$r, 2, max, na.rm=T))
  }
  
  sup_km[sup_km > max(sup_r)] <- NA
  
  result <- rbind(result, cbind(t(colMedians(sup_km, na.rm=T)), t(colMedians(sup_r, na.rm=T))))
}

result


#temp plot
plot(x=n, y=result[, 2], type="o", col="blue", main="ChiSq(2)", xlab="n", ylab="E[sup|F_est-F|]", ylim=c(0.06, 0.2))
lines(x=n, y=result[, 4], type="o", pch=22, lty=2, col="red")
legend(x=750, y=0.18, c("P-L", "R"), col=c("blue", "red"), pch=c(21, 22), lty=c(1, 2))
