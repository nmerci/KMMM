# TODO: implement testing framework
source("R/test.R")

library(matrixStats)


result <- numeric()
for(j in 0:9)
{
  print(paste("j=", j, sep=""))
  
  n <- 100 + 100*j
  m <- 2
  clusters <- replicate_numbers(2 + j*2, 50)
  
  sup_km <- numeric()
  sup_r <- numeric()
  for(i in 1:1000)
  {
    res <- sapply(compare_estimators(n, m, clusters), abs, simplify=F)
    #res <- sapply(compare_estimators(n, m, clusters), "^", 2, simplify = F)
    
    sup_km <- rbind(sup_km, apply(res$km, 2, max, na.rm=T))
    sup_r <- rbind(sup_r, apply(res$r, 2, max, na.rm=T))
  }
  
  sup_km[sup_km > max(sup_r)] <- NA
  
  result <- rbind(result, cbind(t(colMeans(sup_km, na.rm=T)), t(colMeans(sup_r, na.rm=T))))
}

result


