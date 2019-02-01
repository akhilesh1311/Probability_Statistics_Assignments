

Bootstrap_gold <- function(vec0, statfunc, nboot = 10000, alpha = 0.1) {
  n0 <- length(vec0)
  mean0 <- statfunc(vec0)
  sd0 <- Jackknife(vec0, nboot, statfunc, alpha)
  bootvec <- NULL
  
  for(i in 1:nboot) {
    vecb <- sample(vec0, replace=T)
    meanb <- statfunc(vecb)
    sdb <- Jackknife(vecb, nboot, statfunc, alpha)
    bootvec <- c(bootvec, (meanb-mean0)/(sdb/sqrt(n0)))
  }
  lq <- quantile(bootvec, alpha/2)
  uq <- quantile(bootvec, 1-alpha/2)
  
  LB <- mean0-(sd0/sqrt(n0))*uq
  UB <- mean0-(sd0/sqrt(n0))*lq
  
  NLB <- mean0 - (sd0/sqrt(n0))*qnorm(1-alpha/2)
  NUB <- mean0 + (sd0/sqrt(n0))*qnorm(1-alpha/2)
  list(bootstrap.confidence.interval = c(LB, UB), 
       normal.confidence.interval = c(NLB, NUB))
}